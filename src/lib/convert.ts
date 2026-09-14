// Generic unit conversion. Each unit is defined once, as its relationship to
// the SI base/derived unit for its quantity (m, K, Pa, kg/m³, m/s, 1/m³, Pa·s).
// convert() composes those single definitions instead of needing a separate
// hand-written factor for every to/from direction and unit pair.

interface UnitFactor {
	toSI: number;    // multiply a value in this unit by this to get the SI value
	offset?: number; // added after scaling (only nonzero for affine units like °C, °F)
}

const AVOGADRO     = 6.02214076e23;  // 1/mol — SI 2019 exact
const LB_TO_KG	  = 0.45359237;     // kg per lb
const FT3_TO_M3    = 0.028316846592; // m³ per ft³

const UNIT_FACTORS: Record<string, UnitFactor> = {
	// Length (SI: m) — altitude, geopotential/geometric altitude, scale height, mean free path
	'ft': { toSI: 0.3048 },
	'mi': { toSI: 1609.344 },
	'km': { toSI: 1000 },
	'nm': { toSI: 1e-9 },
	'μm': { toSI: 1e-6 },
	'mm': { toSI: 1e-3 },
	'in': { toSI: 0.0254 },

	// Temperature (SI: K)
	'C': { toSI: 1,     offset: 273.15 },
	'F': { toSI: 1/1.8, offset: 273.15 - 32/1.8 },
	'R': { toSI: 1/1.8 },

	// Pressure (SI: Pa)
	'atm':  { toSI: 101325 },
	'hPa':  { toSI: 100 },
	'kPa':  { toSI: 1000 },
	'mmHg': { toSI: 133.322387415 },
	'inHg': { toSI: 3386.388640341 },
	'psi':  { toSI: 6894.75729316836 },

	// Density (SI: kg/m³)
	'lb/ft3':   { toSI: LB_TO_KG / FT3_TO_M3 },
	'slug/ft3': { toSI: 14.5939029372064 / FT3_TO_M3 },
	'sigma':    { toSI: 1.22504083326539 }, // ratio to USSA 1976 sea-level standard density

	// Speed (SI: m/s)
	'km/h': { toSI: 1 / 3.6 },
	'mph':  { toSI: 0.44704 },
	'ft/s': { toSI: 0.3048 },
	'kn':   { toSI: 1852 / 3600 },

	// Number density (SI: 1/m³)
	'per_L':     { toSI: 1000 },
	'per_ft3':   { toSI: 1 / FT3_TO_M3 },
	'mol_m3':    { toSI: AVOGADRO },
	'kmol_m3':   { toSI: AVOGADRO * 1000 },
	'lbmol_ft3': { toSI: AVOGADRO * LB_TO_KG * 1000 / FT3_TO_M3 },

	// Dynamic viscosity (SI: Pa·s)
	'μPa·s':     { toSI: 1e-6 },
	'mPa·s':     { toSI: 1e-3 },
	'lbf·s/ft2': { toSI: 47.8802589803358 },
};

// Converts a value between any two units of the same quantity (e.g. both
// length, or both temperature). Units not found in UNIT_FACTORS — the SI
// base/derived unit itself, e.g. 'm', 'K', 'Pa', 'kg/m3', 'm/s', 'per_m3',
// 'Pa·s' — are treated as identity.
export function convert(value: number, fromUnit: string, toUnit: string): number {
	const from = UNIT_FACTORS[fromUnit] ?? { toSI: 1 };
	const to   = UNIT_FACTORS[toUnit]   ?? { toSI: 1 };
	const si = value * from.toSI + (from.offset ?? 0);
	return (si - (to.offset ?? 0)) / to.toSI;
}
