// USSA 1976 reference implementation using the standard analytical formulas.
// Constant molar mass, geopotential-based temperature layers, no integration.
// Sea-level temperature T0 only affects the troposphere/lower stratosphere
// (H < AIR_TEMP_MAX_GEOPOT); at or above that the profile always matches the fixed
// standard layers, matching the tiered model in src/lib/atmosphere.ts's getTemperature/getPressure.
// Range: −5000 m to 86000 m geometric altitude. Valid T0: -56.5 to 73.5 °C (216.65 to
// 346.65 K, SEA_LEVEL_TEMP_MIN_K/MAX_K in atmosphere.ts). T0 may be given in either unit —
// a value in [-56.5, 73.5] is read as °C (the two valid ranges never overlap numerically).
// Usage: node scripts/atm-test-simple.ts [z_m] [T0_K_or_C] [P0_Pa]

import * as atm from '../src/lib/atmosphere.ts';
import {getBoilingPoint} from '../src/lib/water-properties.ts';
import {convert} from '../src/lib/convert.ts';
import {roundSig} from '../src/lib/helpers.ts';

const G0 = 9.80665;               // standard gravity (m/s²)
const R  = 8314.46261815324;      // gas constant J/(kmol·K)
const M  = 28.9659;               // Sea-level molar mass (g/mol = kg/kmol), updated to 2026
const RE = 6356766;               // USSA 1976 Earth radius (m)

// Sea-level temp only shifts the troposphere/lower-stratosphere (H < AIR_TEMP_MAX_GEOPOT):
// cools at LAPSE from T0 until hitting TROPOPAUSE_T, then holds there up to 20 km. At or
// above 20 km the profile is the fixed standard-atmosphere layers, unaffected by T0.
const TROPOPAUSE_T = 216.65; // K (-56.5 °C)
const LAPSE        = 0.0065; // K/m

// Standard layers from 20 km up (geopotential altitude, base temperature). Lapse rate
// between consecutive breakpoints is derived where needed rather than stored.
// Unaffected by sea-level temperature — base temperatures here are fixed constants.
const KM_86 = atm.geometricToGeopotential(86000);
const LAYERS = [
    { H: 20000, T: 216.65  },
    { H: 32000, T: 228.65  },
    { H: 47000, T: 270.65  },
    { H: 51000, T: 270.65  },
    { H: 71000, T: 214.65  },
    { H: KM_86, T: 186.8673 },
] as const;

// Geometric altitude → geopotential altitude (m).
function _H(z: number): number {
    return (RE * z) / (RE + z);
}

// Geopotential altitude (m) where the troposphere lapse from T0_K reaches TROPOPAUSE_T.
function _Htrop(T0_K: number): number {
    return (T0_K - TROPOPAUSE_T) / LAPSE;
}

// Temperature at geopotential altitude H given sea-level temperature T0_K.
function _T(H: number, T0_K: number): number {
    if (H < atm.AIR_TEMP_MAX_GEOPOT) return Math.max(T0_K - LAPSE * H, TROPOPAUSE_T);
    for (let i = 0; i < LAYERS.length - 1; i++) {
        if (H <= LAYERS[i + 1].H) {
            const t = (H - LAYERS[i].H) / (LAYERS[i + 1].H - LAYERS[i].H);
            return LAYERS[i].T + t * (LAYERS[i + 1].T - LAYERS[i].T);
        }
    }
    return LAYERS[LAYERS.length - 1].T;
}

// Pressure at geopotential altitude H (Pa), analytical USSA 1976 formula.
// Isothermal layer: P = Pb · exp(−M·g₀·ΔH / (R·Tb))
// Gradient layer:   P = Pb · (Tb / T)^(M·g₀ / (R·L))
function _P(H: number, P0: number, T0_K: number): number {
    const Htrop = Math.min(_Htrop(T0_K), 20000);

    // Troposphere: gradient layer from H=0, base temp T0_K, lapse −LAPSE.
    const dH1 = Math.min(H, Htrop);
    let Pb = P0 * Math.pow(T0_K / (T0_K - LAPSE * dH1), -M * G0 / (R * LAPSE));
    if (H <= Htrop) return Pb;

    // Isothermal layer: Htrop to 20 km, at TROPOPAUSE_T.
    const dH2 = Math.min(H, 20000) - Htrop;
    Pb *= Math.exp(-M * G0 * dH2 / (R * TROPOPAUSE_T));
    if (H <= 20000) return Pb;

    // Above 20 km: standard layers, unaffected by T0. Lapse rate derived from consecutive
    // breakpoints' temperatures rather than stored, same as _T.
    for (let i = 0; i < LAYERS.length - 1; i++) {
        const Hb   = LAYERS[i].H;
        const Hnxt = LAYERS[i + 1].H;
        const Tb   = LAYERS[i].T;
        const L    = (LAYERS[i + 1].T - Tb) / (Hnxt - Hb);
        const dH   = Math.min(H, Hnxt) - Hb;
        if (Math.abs(L) < 1e-12) Pb *= Math.exp(-M * G0 * dH / (R * Tb));
        else                     Pb *= Math.pow(Tb / (Tb + L * dH), M * G0 / (R * L));
        if (H <= Hnxt) return Pb;
    }
    return Pb;
}

// Pressure altitude: geometric z where standard atmosphere pressure equals P_Pa.
function _pressureAltitude(P_Pa: number): number {
    let lo = -5000, hi = 86000;
    while (hi - lo > 1e-6) {
        const mid = (lo + hi) / 2;
        _P(_H(mid), 101325, 288.15) > P_Pa ? lo = mid : hi = mid;
    }
    return (lo + hi) / 2;
}

// Density altitude: geometric z where standard atmosphere density equals rho.
function _densityAltitude(rho: number): number {
    let lo = -5000, hi = 86000;
    while (hi - lo > 1e-6) {
        const mid  = (lo + hi) / 2;
        const Hm   = _H(mid);
        const rhoM = _P(Hm, 101325, 288.15) * M / (R * _T(Hm, 288.15));
        rhoM > rho ? lo = mid : hi = mid;
    }
    return (lo + hi) / 2;
}

// ── Main ──────────────────────────────────────────────────────────────────────

const args = process.argv.slice(2);
const z    = args[0] !== undefined ? parseFloat(args[0]) : 0;
let   T0   = args[1] !== undefined ? parseFloat(args[1]) : 288.15;
const P0   = args[2] !== undefined ? parseFloat(args[2]) : 101325;

if (isNaN(z) || isNaN(T0) || isNaN(P0)) {
    console.error('Usage: node scripts/atm-test-simple.ts [z_m] [T0_K or T0_C] [P0_Pa]');
    process.exit(1);
}

// T0 is unambiguous: the valid Celsius range (-56.5 to 73.5) and valid Kelvin range
// (216.65 to 346.65) never overlap, so a value in the Celsius range must be Celsius.
if (T0 >= -56.5 && T0 <= 73.5) T0 += 273.15;

if (T0 < atm.SEA_LEVEL_TEMP_MIN_K || T0 > atm.SEA_LEVEL_TEMP_MAX_K) {
    console.error(`Error: sea-level temperature ${T0} K is outside the valid range ` +
        `[216.65, 346.65] K (-56.5 to 73.5 °C).`);
    process.exit(1);
}

const H   = _H(z);
const T   = _T(H, T0);
const P   = _P(H, P0, T0);
const rho = P * M / (R * T);
const pa  = _pressureAltitude(P);
const da  = _densityAltitude(rho);
const alt = P * 101325 / _P(H, 101325, 288.15);  // altimeter setting
const mu  = atm.getViscosity(T);
const sos = atm.getSpeedOfSound(T, M);
const bp  = getBoilingPoint(P);
const mfp = atm.getMeanFreePath(P, T);
const Geo = atm.geometricToGeopotential(z);
const H_sc = atm.getScaleHeight(T, M, z);
const n    = atm.getNumberDensity(P, T);

function fmtC(K: number): string {
    const C      = K - 273.15;
    const digits = Math.max(0, 5 - Math.floor(Math.log10(Math.abs(K))));
    const s      = C.toFixed(digits);
    return C >= 0 ? `+${s}` : s;
}

function tempLine(label: string, K: number): void {
    console.log(`${label}: ${roundSig(K, 6, 1e-6, true)} K  (${fmtC(K)} °C)`);
}

console.log(`\nInputs`);
console.log(`  Geometric altitude : ${roundSig(z, 6, 1e-6, true)} m`);
tempLine(`  Sea-level temp     `, T0);
console.log(`  Sea-level pressure : ${roundSig(P0, 6, 1e-6, true)} Pa`);

console.log(`\nOutputs`);
console.log(`  Air pressure       : ${roundSig(P, 6, 1e-6, true)} Pa`);
console.log(`  Altimeter setting  : ${roundSig(alt, 6, 1e-6, true)} Pa`);
tempLine(`  Air temperature    `, T);
console.log(`  Air density        : ${roundSig(rho, 6, 1e-6, true)} kg/m³`);
console.log(`  Pressure altitude  : ${roundSig(pa, 6, 1e-6, true)} m`);
console.log(`  Density altitude   : ${roundSig(da, 6, 1e-6, true)} m`);
console.log(`  Geopotential alt   : ${roundSig(Geo, 6, 1e-6, true)} m`);
console.log(`  Molar mass         : ${roundSig(M, 6, 1e-6, true)} kg/kmol`);
console.log(`  Speed of sound     : ${roundSig(sos, 6, 1e-6, true)} m/s`);
console.log(`  Dynamic viscosity  : ${roundSig(mu, 6, 1e-6, true)} Pa·s`);
console.log(`  Mean free path     : ${roundSig(mfp, 6, 1e-6, true)} m`);
console.log(`  Scale height       : ${roundSig(H_sc, 6, 1e-6, true)} m`);
console.log(`  Number density     : ${roundSig(n, 6, 1e-6, true)} /m³  (${roundSig(convert(n, 'per_m3', 'kmol_m3'), 6, 1e-6, true)} kmol/m³)`);
if (bp !== null) tempLine(`  Boiling point      `, bp);
else console.log(`  Boiling point      : N/A (below triple point)`);
