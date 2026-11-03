import { MF_DATA, MF_SPECIES } from '../src/lib/mole-fractions.ts';
import * as atm from '../src/lib/atmosphere.ts';
import { getBoilingPoint } from '../src/lib/water-properties.ts';
import { convert } from '../src/lib/convert.ts';
import { formatTemp, formatAltitude, formatSig } from '../src/lib/helpers.ts';

const args = process.argv.slice(2);

// ── verify-mole-fractions command ────────────────────────────────────────────

if (args[0] === 'verify-mole-fractions') {
    const EPS = 1e-9;
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
let   T0   = args[1] !== undefined ? parseFloat(args[1]) : 288.15;
const P0   = args[2] !== undefined ? parseFloat(args[2]) : 101325;

// T0 is unambiguous: the valid Celsius range (-56.5 to 73.5) and valid Kelvin range
// (216.65 to 346.65) never overlap, so a value in the Celsius range must be Celsius.
if (T0 >= -56.5 && T0 <= 73.5) T0 += 273.15;

if (T0 < atm.SEA_LEVEL_TEMP_MIN_K || T0 > atm.SEA_LEVEL_TEMP_MAX_K) {
    console.error(`Error: sea-level temperature ${T0} K is outside the valid range ` +
        `[216.65, 346.65] K (-56.5 to 73.5 °C).`);
    process.exit(1);
}

if (P0 >= atm.PRESSURE_INPUT_MAX_PA) {
    console.error(`Error: sea-level pressure must be less than ${atm.PRESSURE_INPUT_MAX_PA} Pa.`);
    process.exit(1);
}

const P    = atm.getPressure(z, P0, T0);
const T    = atm.getTemperature(z, T0);
const rho  = atm.getDensity(P, T, z);
const H    = atm.geometricToGeopotential(z);
const mu   = atm.getViscosity(T);
const M    = atm.getMolarMass(z);
const sos  = atm.getSpeedOfSound(T, M);
const bp   = getBoilingPoint(P);
const mfp  = atm.getMeanFreePath(P, T);
const pa   = atm.getPressureAltitude(P);
const da   = atm.getDensityAltitude(rho);
const H_sc = atm.getScaleHeight(T, M, z);
const n    = atm.getNumberDensity(P, T);
const g    = atm.getGravity(z);

function fmtC(K: number): string {
    const C = K - 273.15;
    const s = formatTemp(C);
    return C >= 0 ? `+${s}` : s;
}

function tempLine(label: string, K: number): void {
    console.log(`${label}: ${formatTemp(K, 3, true)} K  (${fmtC(K)} °C)`);
}

console.log(`Inputs`);
console.log(`  Geometric altitude : ${formatSig(z, 6, true)} m`);
console.log(`  Sea-level temp     : ${formatTemp(T0, 3, true)} K  (${fmtC(T0)} °C)`);
console.log(`  Sea-level pressure : ${formatSig(P0, 6, true)} Pa`);

console.log(`\nOutputs`);
console.log(`  Air pressure       : ${formatSig(P, 6, true)} Pa`);
console.log(`  Altimeter setting  : ${formatSig(atm.getAltimeterSetting(z, P), 6, true)} Pa`);
tempLine(`  Air temperature    `, T);
console.log(`  Air density        : ${formatSig(rho, 6, true)} kg/m³`);
console.log(`  Pressure altitude  : ${formatAltitude(pa, 'm')} m`);
console.log(`  Density altitude   : ${formatAltitude(da, 'm')} m`);
console.log(`  Gravity            : ${formatSig(g, 6, true)} m/s²`);
console.log(`  Geopotential alt   : ${formatSig(H, 6, true)} m`);
console.log(`  Molar mass         : ${formatSig(M, 6, true)} kg/kmol`);
console.log(`  Speed of sound     : ${formatSig(sos, 6, true)} m/s`);
console.log(`  Dynamic viscosity  : ${formatSig(mu, 6, true)} Pa·s`);
console.log(`  Mean free path     : ${formatSig(mfp, 6, true)} m`);
console.log(`  Scale height       : ${formatSig(H_sc, 6, true)} m`);
console.log(`  Number density     : ${formatSig(n, 6, true)} /m³  (${formatSig(convert(n, 'per_m3', 'kmol_m3'), 6, true)} kmol/m³)`);
if (bp !== null) tempLine(`  Boiling point      `, bp);
else console.log(`  Boiling point      : N/A (below triple point)`);

console.log(`\nMole fractions`);
for (const species of MF_SPECIES) {
    const frac = atm.getMoleFraction(species, z);
    if (frac > 1e-20) console.log(`  ${species.padEnd(4)}: ${formatSig(frac, 6, true)}`);
}