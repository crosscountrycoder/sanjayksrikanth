import { MF_DATA, MF_SPECIES } from './mole-fractions.ts';
import { tempFromK, pressureFromPa } from './convert.ts';

// ── Section 1: Geodesy ────────────────────────────────────────────────────────

export const EARTH_RADIUS_M = 6356766; // m (USSA 1976)
const G0 = 9.80665; // standard gravity m/s²

export function geometricToGeopotential(z: number): number {
    return (EARTH_RADIUS_M * z) / (EARTH_RADIUS_M + z);
}

export function geopotentialToGeometric(H: number): number {
    return (EARTH_RADIUS_M * H) / (EARTH_RADIUS_M - H);
}

function _g(z: number): number {
    return G0 * (EARTH_RADIUS_M / (EARTH_RADIUS_M + z)) ** 2;
}

// ── Section 2: Temperature ────────────────────────────────────────────────────
// Below 20 km geopotential: sea-level temp T0 sets the baseline, cooling at the
// standard -6.5 K/km lapse rate until it hits 216.65 K, then holding there.
// 20 km to 86 km: piecewise linear lapse-rate layers in geopotential altitude
// (USSA 1976 Table 4), unaffected by T0.
// Above 86 km: USSA 1976 empirical formula defined in geometric altitude.

const KM_86 = geometricToGeopotential(86000);
const TEMP_LAYERS = [
    { H:     0, T: 288.15  },
    { H: 11000, T: 216.65  },
    { H: 20000, T: 216.65  },
    { H: 32000, T: 228.65  },
    { H: 47000, T: 270.65  },
    { H: 51000, T: 270.65  },
    { H: 71000, T: 214.65  },
    { H: KM_86, T: 186.8673 },
] as const;

const TROPOPAUSE_T = 216.65;   // K (-56.5 °C) — floor of the sea-level-driven layer
const LAPSE_RATE   = 6.5 / 1000; // K/m
const EPS = 1e-9; // floating-point slack for boundary comparisons
export const SEA_LEVEL_TEMP_MIN_K = 216.65 - EPS; // -56.5 °C
export const SEA_LEVEL_TEMP_MAX_K = 346.65 + EPS; // 73.5 °C

export function getTemperature(z: number, T0_K: number = 288.15): number {
    if (T0_K < SEA_LEVEL_TEMP_MIN_K || T0_K > SEA_LEVEL_TEMP_MAX_K) return NaN;
    if (z > 86000) {
        if (z < 91000)  return 186.8673;
        if (z < 110000) { const b = (z - 91000) / 19942.9; return 263.1905 - 76.3232 * Math.sqrt(1 - b * b); }
        if (z < 120000) return 240.0 + 0.012 * (z - 110000);
        const b = (EARTH_RADIUS_M + 120000) / (EARTH_RADIUS_M + z);
        return 1000 - 640 * Math.exp(-1.875e-5 * (z - 120000) * b);
    }
    const H = geometricToGeopotential(z);
    if (H <= 20000) {
        return Math.max(T0_K - LAPSE_RATE * H, TROPOPAUSE_T);
    }
    const L = TEMP_LAYERS;
    for (let i = 0; i < L.length - 1; i++) {
        if (H <= L[i + 1].H) {
            const t = (H - L[i].H) / (L[i + 1].H - L[i].H);
            return L[i].T + t * (L[i + 1].T - L[i].T);
        }
    }
    return L[L.length - 1].T;
}

// ── Section 3: Composition and molar mass ────────────────────────────────────
// Molar mass and mole fractions are read from the pre-computed table in
// mole-fractions.ts, which covers −6000 m to 1 000 000 m in 1 km steps.
// Values are linearly interpolated between table rows.

export const R = 8314.46261815324; // J/(kmol·K) — SI 2019; molar masses are in g/mol = kg/kmol
const k_B = 1.380649e-23;          // J/K        — SI 2019 exact

const MF_Z_MIN  = 0; // CSV starts at sea level; negative altitudes clamp to sea-level composition
const MF_Z_STEP = 1000;

// Interpolate a single column from MF_DATA at altitude z (m).
// col 1 = molar mass, col 2+ = mole fractions in MF_SPECIES order.
function _mfInterp(z: number, col: number): number {
    const z_lo = Math.floor(z / MF_Z_STEP) * MF_Z_STEP;
    const idx  = (z_lo - MF_Z_MIN) / MF_Z_STEP;
    // Clamp to table bounds
    if (idx < 0) return MF_DATA[0]![col];
    if (idx >= MF_DATA.length - 1) return MF_DATA[MF_DATA.length - 1]![col];
    const lo = MF_DATA[idx]!;
    if (z === z_lo) return lo[col];
    const hi = MF_DATA[idx + 1]!;
    return lo[col] + ((z - z_lo) / MF_Z_STEP) * (hi[col] - lo[col]);
}

// Mean molar mass (kg/mol) at geometric altitude z (m).
export function getMolarMass(z: number): number {
    return _mfInterp(z, 1);
}

// Mole fraction of a species at geometric altitude z (m).
export function getMoleFraction(species: string, z: number): number {
    const si = (MF_SPECIES as readonly string[]).indexOf(species);
    if (si < 0) return 0;
    return _mfInterp(z, si + 2);
}

// Partial pressure of a species (Pa) given its mole fraction and total pressure.
export function getPartialPressure(moleFraction: number, totalPressure: number): number {
    return moleFraction * totalPressure;
}

export function getDensity(P_Pa: number, T_K: number, z: number): number {
    return P_Pa * getMolarMass(z) / (R * T_K);
}

// ── Section 4: Pressure integration ──────────────────────────────────────────
// Hydrostatic equation: dP/dz = −P · M(z) · g(z) / (R · T(z))
// Unified for all altitudes from −6 km to 1000 km.
// M(z) is read from the pre-computed table via getMolarMass().

const Z_STEP = 50;   // integration step (m)
const Z_MIN  = -5000;
const Z_MAX  = 1000000;

// One RK4 step of the unified hydrostatic equation.
function _rk4P(z: number, P: number, h: number, T0_K: number): number {
    const f = (zi: number) => getMolarMass(zi) * _g(zi) / (R * getTemperature(zi, T0_K));
    const k1 = -P               * f(z);
    const k2 = -(P + 0.5*h*k1) * f(z + 0.5*h);
    const k3 = -(P + 0.5*h*k2) * f(z + 0.5*h);
    const k4 = -(P + h*k3)     * f(z + h);
    return P + (h / 6) * (k1 + 2*k2 + 2*k3 + k4);
}

// Integrate the hydrostatic equation from z0 to z1 starting at pressure P.
function _integrateP(z0: number, z1: number, P: number, T0_K: number): number {
    if (z0 === z1) return P;
    const dir   = z1 > z0 ? 1 : -1;
    const dist  = Math.abs(z1 - z0);
    const nFull = Math.floor(dist / Z_STEP);
    let Z = z0;
    for (let i = 0; i < nFull; i++) {
        P = _rk4P(Z, P, dir * Z_STEP, T0_K);
        Z += dir * Z_STEP;
    }
    const rem = z1 - Z;
    if (Math.abs(rem) > 1e-9) P = _rk4P(Z, P, rem, T0_K);
    return P;
}

export function getPressure(z: number, P0: number, T0_K: number): number {
    return _integrateP(0, z, P0, T0_K);
}

// Inverts the tiered temperature model: given an air temperature at z, finds
// the sea-level temperature that produces it. Only solvable at or below 20 km
// geopotential (T0 has no effect above that), and only for temperatures that
// are actually reachable (>= 216.65 K, the tropopause floor).

export function getSeaLevelTemperature(z: number, T_K: number): number {
    const H = geometricToGeopotential(z);
    if (H > 20000 || T_K < SEA_LEVEL_TEMP_MIN_K) return NaN;
    const T0_K = T_K + LAPSE_RATE * H;
    if (T0_K < SEA_LEVEL_TEMP_MIN_K || T0_K > SEA_LEVEL_TEMP_MAX_K) return NaN;
    return T0_K;
}

export function getSeaLevelPressure(z: number, P_Pa: number, T0_K: number): number {
    return P_Pa * 101325 / getPressure(z, 101325, T0_K);
}

export function getAltimeterSetting(z: number, P_Pa: number): number {
    return P_Pa * 101325 / getPressure(z, 101325, 288.15);
}

// ── Standard atmosphere lookup table ─────────────────────────────────────────
// Precomputed pressure and density at every 50 m from −5 km to 1000 km.
// Used by getPressureAltitude and getDensityAltitude.

const _STD_N   = ((Z_MAX - Z_MIN) / Z_STEP) + 1; // 10061 entries
const _stdP   = new Float64Array(_STD_N);
const _stdRho = new Float64Array(_STD_N);
(function _buildStdTable() {
    const i0 = (0 - Z_MIN) / Z_STEP; // index of z = 0
    _stdP[i0] = 101325;
    for (let i = i0; i < _STD_N - 1; i++) {
        _stdP[i + 1] = _rk4P(Z_MIN + i * Z_STEP, _stdP[i], Z_STEP, 288.15);
    }
    for (let i = i0; i > 0; i--) {
        _stdP[i - 1] = _rk4P(Z_MIN + i * Z_STEP, _stdP[i], -Z_STEP, 288.15);
    }
    for (let i = 0; i < _STD_N; i++) {
        const z = Z_MIN + i * Z_STEP;
        _stdRho[i] = _stdP[i] * getMolarMass(z) / (R * getTemperature(z));
    }
}());

function _bsearchDecreasing(table: Float64Array, value: number): number {
    let lo = 0, hi = _STD_N - 1;
    while (lo < hi - 1) {
        const mid = Math.floor((lo + hi) / 2);
        table[mid] >= value ? lo = mid : hi = mid;
    }
    return lo;
}

// Pressure altitude: altitude in the standard atmosphere with pressure P_Pa.
// Returns NaN if that altitude would fall outside -5 km to 1000 km.
export function getPressureAltitude(P_Pa: number): number {
    if (P_Pa > _stdP[0] || P_Pa < _stdP[_STD_N - 1]) return NaN;
    const i    = _bsearchDecreasing(_stdP, P_Pa);
    const z_lo = Z_MIN + i * Z_STEP;
    const P_lo = _stdP[i];
    let lo = 0, hi = Z_STEP;
    while (hi - lo > 1e-6) {
        const mid = (lo + hi) / 2;
        _rk4P(z_lo, P_lo, mid, 288.15) > P_Pa ? lo = mid : hi = mid;
    }
    return z_lo + lo;
}

// Density altitude: altitude in the standard atmosphere with density rho.
// Returns NaN if that altitude would fall outside -5 km to 1000 km.
export function getDensityAltitude(rho: number): number {
    if (rho > _stdRho[0] || rho < _stdRho[_STD_N - 1]) return NaN;
    const i    = _bsearchDecreasing(_stdRho, rho);
    const z_lo = Z_MIN + i * Z_STEP;
    const P_lo = _stdP[i];
    let lo = 0, hi = Z_STEP;
    while (hi - lo > 1e-6) {
        const mid = (lo + hi) / 2;
        const rho_mid = _rk4P(z_lo, P_lo, mid, 288.15) * getMolarMass(z_lo + mid) / (R * getTemperature(z_lo + mid));
        rho_mid > rho ? lo = mid : hi = mid;
    }
    return z_lo + lo;
}

// Standard-atmosphere pressure (Pa) at geometric altitude z, read directly from the
// precomputed table (O(1) interpolation, unlike getPressure which re-integrates from
// scratch). Returns NaN outside -5 km to 1000 km.
export function getStandardPressure(z: number): number {
    if (z < Z_MIN || z > Z_MAX) return NaN;
    const idx_f = (z - Z_MIN) / Z_STEP;
    const idx   = Math.min(Math.floor(idx_f), _STD_N - 2);
    const frac  = idx_f - idx;
    return _stdP[idx] + frac * (_stdP[idx + 1] - _stdP[idx]);
}

// Formats a temperature value for display: 6 significant figures counted against the
// *absolute* temperature (Kelvin for °C/K, Rankine for °F/°R), not against the display
// value directly — this is why an ordinary temperature like 15 °C or -17.4745 °C renders
// with 3 decimal places rather than 6 significant digits of its own magnitude. Shared by
// all three atmosphere calculator pages (calculator, table, graph) to avoid duplicating it.
// trailingZeros (default true) controls whether e.g. 15 renders as "15.000" or "15".
export function fmtTemp(v: number, unit: string, sigFigs = 6, trailingZeros = true): string {
    if (unit === 'K' || unit === 'R') {
        const s = v.toPrecision(sigFigs);
        return trailingZeros ? s : parseFloat(s).toString();
    }
    const T_abs = v + (unit === 'C' ? 273.15 : 459.67);
    const digits = Math.max(0, (sigFigs - 1) - Math.floor(Math.log10(Math.abs(T_abs))));
    const s = v.toFixed(digits);
    return trailingZeros ? s : parseFloat(s).toString();
}

// Builds the "Standard conditions" hint shown near the sea-level temperature/pressure
// inputs. When a field is in "Air temperature"/"Air pressure" mode, its standard-atmosphere
// reference value depends on altitude rather than always being 15 °C / 1 atm, so the note
// is phrased accordingly. Returns null if an altitude-dependent value is needed but z is
// invalid or outside the model's range.
export function getStandardConditionsNote(
    z: number,
    tempMode: 'sea-level' | 'air',
    pressMode: 'sea-level' | 'air' | 'altimeter',
    tempUnit: string,
    pressUnit: string,
): string | null {
    const tempNeedsAlt  = tempMode === 'air';
    const pressNeedsAlt = pressMode === 'air';
    if ((tempNeedsAlt || pressNeedsAlt) && isNaN(z)) return null;

    const stdT_K  = tempNeedsAlt  ? getTemperature(z) : 288.15;
    const stdP_Pa = pressNeedsAlt ? getStandardPressure(z) : 101325;
    if (isNaN(stdT_K) || isNaN(stdP_Pa)) return null;

    const tempUnitLabel = tempUnit === 'C' ? '°C' : tempUnit === 'F' ? '°F' : tempUnit === 'R' ? '°R' : 'K';
    const T = fmtTemp(tempFromK(stdT_K, tempUnit), tempUnit, 6, false);
    const P = parseFloat(pressureFromPa(stdP_Pa, pressUnit).toPrecision(6)).toString();

    if (!tempNeedsAlt && !pressNeedsAlt) {
        return `Standard sea-level conditions are ${T} ${tempUnitLabel} and ${P} ${pressUnit}.`;
    }
    if (tempNeedsAlt && pressNeedsAlt) {
        return `Standard conditions at this altitude are ${T} ${tempUnitLabel} and ${P} ${pressUnit}.`;
    }
    if (pressNeedsAlt) {
        return `Standard sea-level temperature is ${T} ${tempUnitLabel}; standard pressure at this altitude is ${P} ${pressUnit}.`;
    }
    const pressPhrase = pressMode === 'altimeter' ? 'altimeter setting' : 'sea level pressure';
    return `Standard air temperature at this altitude is ${T} ${tempUnitLabel}; standard ${pressPhrase} is ${P} ${pressUnit}.`;
}

export function getSpeedOfSound(T_K: number, M: number): number {
    return Math.sqrt(1.4 * R * T_K / M);
}

// Atmospheric scale height (m): RT/(Mg) — the altitude over which pressure/density would
// decay by a factor of e if temperature, molar mass, and gravity stayed constant at their
// values at z. Uses the local gravity at z, not standard sea-level gravity.
export function getScaleHeight(T_K: number, M: number, z: number): number {
    return R * T_K / (M * _g(z));
}

// Sutherland's formula for dynamic viscosity of air (Pa·s).
export function getViscosity(T_K: number): number {
    return 1.458e-6 * T_K ** 1.5 / (T_K + 110.4);
}

const SIGMA_MFP = 3.65e-10;
export function getMeanFreePath(P_Pa: number, T_K: number): number {
    const n = P_Pa / (k_B * T_K);
    return 1 / (Math.SQRT2 * Math.PI * SIGMA_MFP * SIGMA_MFP * n);
}

// ── Atmosphere profile builder ────────────────────────────────────────────────
// Integrates pressure continuously across a sorted altitude array (O(N) total).
// zPoints must be sorted ascending. Returns parallel P, T, M arrays.

export interface AtmosphereProfile {
    readonly P: Float64Array;
    readonly T: Float64Array;
    readonly M: Float64Array;
}

export function buildAtmosphereProfile(zPoints: number[], P0: number, T0: number): AtmosphereProfile {
    const n   = zPoints.length;
    const P_o = new Float64Array(n);
    const T_o = new Float64Array(n);
    const M_o = new Float64Array(n);
    if (n === 0) return { P: P_o, T: T_o, M: M_o };

    P_o[0] = _integrateP(0, zPoints[0], P0, T0);
    T_o[0] = getTemperature(zPoints[0], T0);
    M_o[0] = getMolarMass(zPoints[0]);
    for (let k = 1; k < n; k++) {
        P_o[k] = _integrateP(zPoints[k - 1], zPoints[k], P_o[k - 1], T0);
        T_o[k] = getTemperature(zPoints[k], T0);
        M_o[k] = getMolarMass(zPoints[k]);
    }
    return { P: P_o, T: T_o, M: M_o };
}

