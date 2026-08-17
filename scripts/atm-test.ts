import { MF_DATA, MF_SPECIES } from '../src/lib/mole-fractions.ts';
import {
    geometricToGeopotential,
    getPressure, getTemperature, getDensity,
    getAltimeterSetting,
    getPressureAltitude, getDensityAltitude,
    getMoleFraction, getMolarMass, getMeanFreePath,
    getSpeedOfSound,
} from '../src/lib/atmosphere.ts';
import { getBoilingPoint } from '../src/lib/water-properties.ts';

const args = process.argv.slice(2);

// ── verify-mole-fractions command ────────────────────────────────────────────

if (args[0] === 'verify-mole-fractions') {
    const EPS = 1e-13;
    let failures = 0;

    for (const row of MF_DATA) {
        const z   = row[0];
        const sum = MF_SPECIES.reduce((s, _, i) => s + row[i + 2], 0);
        const err = Math.abs(sum - 1);
        if (err > EPS) {
            console.log(`FAIL  z=${z} m  sum=${sum.toPrecision(15)}  |err|=${err.toExponential(3)}`);
            failures++;
        }
    }

    if (failures === 0) console.log(`OK  all ${MF_DATA.length} altitudes sum to 1 within ${EPS}`);
    else console.log(`${failures} failure(s)`);
    process.exit(failures > 0 ? 1 : 0);
}

// ── single-point calculator ───────────────────────────────────────────────────

const z    = args[0] !== undefined ? parseFloat(args[0]) : 0;
const T0   = args[1] !== undefined ? parseFloat(args[1]) : 288.15;
const P0   = args[2] !== undefined ? parseFloat(args[2]) : 101325;

const dT   = T0 - 288.15;
const P    = getPressure(z, P0, T0);
const T    = getTemperature(z) + dT;
const rho  = getDensity(P, T, z);
const H    = geometricToGeopotential(z);
const mu   = 1.458e-6 * T ** 1.5 / (T + 110.4);
const sos  = getSpeedOfSound(P, rho);
const bp   = getBoilingPoint(P);
const mfp  = getMeanFreePath(P, T);
const pa   = getPressureAltitude(P);
const da   = getDensityAltitude(rho);

function sf6(v: number): string {
    return parseFloat(v.toPrecision(6)).toString();
}

function fmtC(K: number): string {
    const C = K - 273.15;
    const digits = Math.max(0, 5 - Math.floor(Math.log10(Math.abs(K))));
    const s = C.toFixed(digits);
    return C >= 0 ? `+${s}` : s;
}

function tempLine(label: string, K: number): void {
    console.log(`${label}: ${parseFloat(K.toPrecision(6))} K  (${fmtC(K)} °C)`);
}

console.log(`\nInputs`);
console.log(`  Geometric altitude : ${sf6(z)} m`);
console.log(`  Sea-level temp     : ${parseFloat(T0.toPrecision(6))} K  (${fmtC(T0)} °C)`);
console.log(`  Sea-level pressure : ${sf6(P0)} Pa`);

console.log(`\nOutputs`);
console.log(`  Air pressure       : ${sf6(P)} Pa`);
console.log(`  Altimeter setting  : ${sf6(getAltimeterSetting(z, P))} Pa`);
tempLine(`  Air temperature    `, T);
console.log(`  Air density        : ${sf6(rho)} kg/m³`);
console.log(`  Pressure altitude  : ${sf6(pa)} m`);
console.log(`  Density altitude   : ${sf6(da)} m`);
console.log(`  Geopotential alt   : ${sf6(H)} m`);
console.log(`  Speed of sound     : ${sos !== null ? sf6(sos) + ' m/s' : 'N/A'}`);
console.log(`  Dynamic viscosity  : ${sf6(mu)} Pa·s`);
console.log(`  Mean free path     : ${sf6(mfp)} m`);
if (bp !== null) tempLine(`  Boiling point      `, bp);
else console.log(`  Boiling point      : N/A (below triple point)`);

console.log(`\nMole fractions`);
for (const species of MF_SPECIES) {
    const frac = getMoleFraction(species, z);
    if (frac > 1e-20) console.log(`  ${species.padEnd(4)}: ${sf6(frac)}`);
}
