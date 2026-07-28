import { MF_DATA, MF_SPECIES, type MFSpecies } from './mole-fractions.ts';

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
// Below 86 km: piecewise linear lapse-rate layers in geometric altitude.
// Above 86 km: USSA 1976 empirical formula with five segments.
// Extrapolates below 0 m using the tropospheric lapse rate.

const TEMP_LAYERS = [
    { z:     0, T: 288.15   },
    { z: 11000, T: 216.65   },
    { z: 20100, T: 216.65   },
    { z: 32200, T: 228.65   },
    { z: 47400, T: 270.65   },
    { z: 51400, T: 270.65   },
    { z: 71800, T: 214.65   },
    { z: 86000, T: 186.8673 },
] as const;

export function getTemperature(z: number): number {
    if (z > 86000) {
        if (z < 91000)  return 186.8673;
        if (z < 110000) { const b = (z - 91000) / 19942.9; return 263.1905 - 76.3232 * Math.sqrt(1 - b * b); }
        if (z < 120000) return 240.0 + 0.012 * (z - 110000);
        const b = (EARTH_RADIUS_M + 120000) / (EARTH_RADIUS_M + z);
        return 1000 - 640 * Math.exp(-1.875e-5 * (z - 120000) * b);
    }
    const L = TEMP_LAYERS;
    if (z <= L[0].z) {
        const dTdz = (L[1].T - L[0].T) / (L[1].z - L[0].z);
        return L[0].T + dTdz * (z - L[0].z);
    }
    for (let i = 0; i < L.length - 1; i++) {
        if (z <= L[i + 1].z) {
            const t = (z - L[i].z) / (L[i + 1].z - L[i].z);
            return L[i].T + t * (L[i + 1].T - L[i].T);
        }
    }
    return L[L.length - 1].T;
}

// ── Section 3: Composition and molar mass ────────────────────────────────────
// Molar mass and mole fractions are read from the pre-computed table in
// mole-fractions.ts, which covers −6000 m to 1 000 000 m in 1 km steps.
// Values are linearly interpolated between table rows.

export const R   = 8.31446261815324; // J/(mol·K) — SI 2019
const k_B = 1.380649e-23;            // J/K       — SI 2019 exact

// Sea-level dry-air composition (reference data). Molar masses in kg/mol.
// Sources: NOAA GML (CO2, CH4, N2O, 2024 annual means); USSA 1976 (all others).
export const LOWER_COMPOSITION = [
    { key: 'N2',  xi: 0.780859353, Mi: 0.02801340 },
    { key: 'O2',  xi: 0.209351263, Mi: 0.03199880 },
    { key: 'Ar',  xi: 0.009340000, Mi: 0.03994800 },
    { key: 'CO2', xi: 0.000422000, Mi: 0.04400950 },
    { key: 'Ne',  xi: 0.000018180, Mi: 0.02017970 },
    { key: 'He',  xi: 0.000005240, Mi: 0.00400260 },
    { key: 'CH4', xi: 0.000001900, Mi: 0.01604246 },
    { key: 'Kr',  xi: 0.000001140, Mi: 0.08379800 },
    { key: 'H2',  xi: 0.000000500, Mi: 0.00201588 },
    { key: 'N2O', xi: 0.000000337, Mi: 0.04401280 },
    { key: 'Xe',  xi: 0.000000087, Mi: 0.13129300 },
] as const;

export const M0 = LOWER_COMPOSITION.reduce((s, c) => s + c.xi * c.Mi, 0);

const MF_Z_MIN  = -5000;
const MF_Z_STEP = 1000;

// Interpolate a single column from MF_DATA at altitude z (m).
// col 1 = molar mass, col 2+ = mole fractions in MF_SPECIES order.
function _mfInterp(z: number, col: number): number {
    if (MF_DATA.length === 0) return col === 1 ? M0 : 0;
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
export function getMoleFraction(species: MFSpecies, z: number): number {
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

const Z_STEP = 100;   // integration step (m)
const Z_MIN  = -5000;
const Z_MAX  = 1000000;

// One RK4 step of the unified hydrostatic equation.
function _rk4P(z: number, P: number, h: number, dT: number): number {
    const f = (zi: number) => getMolarMass(zi) * _g(zi) / (R * (getTemperature(zi) + dT));
    const k1 = -P               * f(z);
    const k2 = -(P + 0.5*h*k1) * f(z + 0.5*h);
    const k3 = -(P + 0.5*h*k2) * f(z + 0.5*h);
    const k4 = -(P + h*k3)     * f(z + h);
    return P + (h / 6) * (k1 + 2*k2 + 2*k3 + k4);
}

// Integrate the hydrostatic equation from z0 to z1 starting at pressure P.
function _integrateP(z0: number, z1: number, P: number, dT: number): number {
    if (z0 === z1) return P;
    const dir   = z1 > z0 ? 1 : -1;
    const dist  = Math.abs(z1 - z0);
    const nFull = Math.floor(dist / Z_STEP);
    let Z = z0;
    for (let i = 0; i < nFull; i++) {
        P = _rk4P(Z, P, dir * Z_STEP, dT);
        Z += dir * Z_STEP;
    }
    const rem = z1 - Z;
    if (Math.abs(rem) > 1e-9) P = _rk4P(Z, P, rem, dT);
    return P;
}

export function getPressure(z: number, P0: number, T0: number): number {
    return _integrateP(0, z, P0, T0 - 288.15);
}

export function getSeaLevelTemperature(z: number, T_K: number): number {
    return 288.15 + T_K - getTemperature(z);
}

export function getSeaLevelPressure(z: number, P_Pa: number, T0_K: number): number {
    return P_Pa * 101325 / getPressure(z, 101325, T0_K);
}

export function getAltimeterSetting(z: number, P_Pa: number): number {
    return P_Pa * 101325 / getPressure(z, 101325, 288.15);
}

// ── Standard atmosphere lookup table ─────────────────────────────────────────
// Precomputed pressure and density at every 100 m from −6 km to 1000 km.
// Used by getPressureAltitude and getDensityAltitude.

const _STD_N   = ((Z_MAX - Z_MIN) / Z_STEP) + 1; // 10061 entries
const _stdP   = new Float64Array(_STD_N);
const _stdRho = new Float64Array(_STD_N);
(function _buildStdTable() {
    const i0 = (0 - Z_MIN) / Z_STEP; // index of z = 0
    _stdP[i0] = 101325;
    for (let i = i0; i < _STD_N - 1; i++) {
        _stdP[i + 1] = _rk4P(Z_MIN + i * Z_STEP, _stdP[i], Z_STEP, 0);
    }
    for (let i = i0; i > 0; i--) {
        _stdP[i - 1] = _rk4P(Z_MIN + i * Z_STEP, _stdP[i], -Z_STEP, 0);
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
export function getPressureAltitude(P_Pa: number, fast = false): number {
    const i    = _bsearchDecreasing(_stdP, P_Pa);
    const z_lo = Z_MIN + i * Z_STEP;
    const P_lo = _stdP[i];
    if (fast) {
        const t = Math.log(P_Pa / P_lo) / Math.log(_stdP[i + 1] / P_lo);
        return Math.max(Z_MIN, Math.min(z_lo + t * Z_STEP, Z_MAX));
    }
    let lo = 0, hi = Z_STEP;
    while (hi - lo > 1e-6) {
        const mid = (lo + hi) / 2;
        _rk4P(z_lo, P_lo, mid, 0) > P_Pa ? lo = mid : hi = mid;
    }
    return z_lo + lo;
}

// Density altitude: altitude in the standard atmosphere with density rho.
export function getDensityAltitude(rho: number, fast = false): number {
    const i    = _bsearchDecreasing(_stdRho, rho);
    const z_lo = Z_MIN + i * Z_STEP;
    const P_lo = _stdP[i];
    if (fast) {
        const t = Math.log(rho / _stdRho[i]) / Math.log(_stdRho[i + 1] / _stdRho[i]);
        return Math.max(Z_MIN, Math.min(z_lo + t * Z_STEP, Z_MAX));
    }
    let lo = 0, hi = Z_STEP;
    while (hi - lo > 1e-6) {
        const mid = (lo + hi) / 2;
        const rho_mid = _rk4P(z_lo, P_lo, mid, 0) * getMolarMass(z_lo + mid) / (R * getTemperature(z_lo + mid));
        rho_mid > rho ? lo = mid : hi = mid;
    }
    return z_lo + lo;
}

export function getSpeedOfSound(P_Pa: number, rho: number): number | null {
    if (P_Pa < 5e-4) return null; // anacoustic zone
    return Math.sqrt(1.4 * P_Pa / rho);
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

    const dT = T0 - 288.15;
    P_o[0] = _integrateP(0, zPoints[0], P0, dT);
    T_o[0] = getTemperature(zPoints[0]) + dT;
    M_o[0] = getMolarMass(zPoints[0]);
    for (let k = 1; k < n; k++) {
        P_o[k] = _integrateP(zPoints[k - 1], zPoints[k], P_o[k - 1], dT);
        T_o[k] = getTemperature(zPoints[k]) + dT;
        M_o[k] = getMolarMass(zPoints[k]);
    }
    return { P: P_o, T: T_o, M: M_o };
}

// ── Section 5: Water vapor and boiling point ──────────────────────────────────
// IAPWS-95 saturation pressure formula, accurate to 0.01°C.
// Below the triple point (611.657 Pa) liquid water cannot exist.

export function waterVaporPressure(T: number): number {
    const Tc = 647.096, Pc = 22064000, theta = 1 - T / Tc;
    return Pc * Math.exp((Tc / T) * (
        -7.85951783*theta + 1.84408259*theta**1.5 - 11.7866497*theta**3 +
        22.6807411*theta**3.5 - 15.9618719*theta**4 + 1.80122502*theta**7.5
    ));
}

export function getBoilingPoint(P_Pa: number): number {
    let lo = 273, hi = 647.096;
    for (let i = 0; i < 30; i++) {
        const mid = (lo + hi) / 2;
        waterVaporPressure(mid) < P_Pa ? lo = mid : hi = mid;
    }
    return hi;
}
