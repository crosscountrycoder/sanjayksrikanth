import { convert } from './convert.ts';
import { getTemperature, getStandardPressure } from './atmosphere.ts';

// ── Number formatting ────────────────────────────────────────────────────────
// Every displayed number goes through one of these three:
//   formatTemp      temperatures (including dew point and boiling point)
//   formatAltitude  pressure altitude and density altitude
//   formatSig       everything else (including geometric/geopotential altitude)

// toFixed can produce "-0.000" for tiny negative values; show those as unsigned zero.
function fixed(n: number, digits: number): string {
	const s = n.toFixed(digits);
	return /^-0\.?0*$/.test(s) ? s.slice(1) : s;
}

// Temperatures: rounded to `decimals` digits after the decimal point (3 for calculators and
// tables, 1 for graph tick labels), except from 1000 up (999.9995 and up for 3 decimals, i.e.
// wherever that rounding would reach 1000), where it rounds to decimals + 3 significant figures
// instead, so a temperature never shows more digits than that (1340.329 → 1340.33). Output fields
// keep trailing zeros ("15.000"); set `strip` for graph labels, error-message bounds, and values
// written into input boxes ("15").
export function formatTemp(n: number, decimals: number = 3, strip: boolean = false): string {
	let s = fixed(n, decimals);
	if (Math.abs(parseFloat(s)) >= 1000) {
		// Round to significant figures first, so a carry into a new digit (9999.996 → 10000.0)
		// is counted when choosing the number of decimal places.
		const sigFigs = decimals + 3;
		const rounded = Number(n.toPrecision(sigFigs));
		const intDigits = Math.floor(Math.log10(Math.abs(rounded))) + 1;
		s = fixed(rounded, Math.max(0, sigFigs - intDigits));
	}
	return strip ? parseFloat(s).toString() : s;
}

// Pressure/density altitude: nearest integer in m or ft, 3 digits after the decimal point in
// mi or km (~1 m either way); finer precision would just be binary-search noise.
export function formatAltitude(n: number, unit: string): string {
	return (unit === 'm' || unit === 'ft') ? fixed(n, 0) : fixed(n, 3);
}

// Rounds to `sigFigs` significant figures (default 6). Output fields keep trailing zeros
// ("1.00000") and switch to scientific notation below 1e-4 or at 10^sigFigs and above.
// Set `strip` for error-message bounds and values written into input boxes: trailing zeros are
// removed ("1") and JS's own notation rules apply (scientific only below 1e-6 or from 1e21).
export function formatSig(n: number, sigFigs: number = 6, strip: boolean = false): string {
	if (n === 0) return strip ? '0' : (0).toPrecision(sigFigs);
	const rounded = Number(n.toPrecision(sigFigs));
	if (strip) return rounded.toString();
	return Math.abs(rounded) < 1e-4 ? rounded.toExponential(sigFigs - 1) : rounded.toPrecision(sigFigs);
}

// Builds the "Standard conditions" hint shown near the sea-level temperature/pressure
// inputs. When a field is in "Air temperature"/"Air pressure" mode, its standard-atmosphere
// reference value depends on altitude rather than always being 15 °C / 1 atm, so the note
// is phrased accordingly. Returns null if an altitude-dependent value is needed but z is
// invalid or outside the model's range. Lives here rather than atmosphere.ts because it
// only builds a display message — it doesn't calculate any atmospheric property itself.
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
	const T = formatTemp(convert(stdT_K, 'K', tempUnit), 3, true);
	const P = formatSig(convert(stdP_Pa, 'Pa', pressUnit), 6, true);

	if (!tempNeedsAlt && !pressNeedsAlt) {
		if (pressMode === 'altimeter') {
			return `Standard sea-level temperature is ${T} ${tempUnitLabel}; standard altimeter setting is ${P} ${pressUnit}.`;
		}
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
