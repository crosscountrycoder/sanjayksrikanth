import { convert } from './convert.ts';
import { getTemperature, getStandardPressure } from './atmosphere.ts';

// Rounds a number to `sigFigs` significant figures (default 6), switching to scientific
// notation when the *rounded* magnitude falls below `sciThreshold` (default 1e-4, matching
// most user-visible output fields) — rounding first so a value like 9.99999999e-5 that
// rounds up to 1e-4 is still shown as fixed-point rather than scientific. Pass 1e-6 to match
// native toPrecision behavior instead (used for bound-check messages, reset-to-default
// values, and reference-script output, which favor the plainer/native threshold). Set
// `stripTrailingZeros` (default false) to compact the result (e.g. "1" instead of
// "1.00000") for error-message bounds and input values, as opposed to aligned output
// columns, which keep trailing zeros. Note: stripping goes through Number(...).toString(),
// which uses JS's own fixed 1e-6/1e21 notation thresholds — with sciThreshold set to 1e-6
// that exactly reproduces toPrecision's choice, but any other sciThreshold combined with
// stripTrailingZeros can have its notation reverted by that round-trip.
//
// `tempUnit` handles the one case that doesn't fit the above: °C and °F are offset scales,
// so counting sigFigs against the displayed value itself would be meaningless right around
// 0. Pass the temperature unit ('C' or 'F' specifically trigger it) to instead count sigFigs
// against the absolute temperature (Kelvin or Rankine) while still displaying in the original
// unit — e.g. 15 °C (288.15 K) needs 3 decimal places for 6 sig figs, so it renders "15.000"
// rather than 6 digits counted from "15" itself. This bypasses sciThreshold entirely
// (temperatures never use scientific notation here). Leave it '' (default) for every other
// quantity — including K and °R, which are already absolute scales and fall through to the
// generic path below unaffected, so callers can pass a temperature-unit variable straight
// through without narrowing it first.
export function roundSig(
	n: number,
	sigFigs: number = 6,
	sciThreshold: number = 1e-4,
	stripTrailingZeros: boolean = false,
	tempUnit: string = '',
): string {
	if (tempUnit === 'C' || tempUnit === 'F') {
		const T_abs = n + (tempUnit === 'C' ? 273.15 : 459.67);
		const digits = Math.max(0, (sigFigs - 1) - Math.floor(Math.log10(Math.abs(T_abs))));
		const s = n.toFixed(digits);
		return stripTrailingZeros ? parseFloat(s).toString() : s;
	}
	if (n === 0) return stripTrailingZeros ? '0' : (0).toPrecision(sigFigs);
	const rounded = Number(n.toPrecision(sigFigs));
	const s = Math.abs(rounded) < sciThreshold
		? rounded.toExponential(sigFigs - 1)
		: rounded.toPrecision(sigFigs);
	return stripTrailingZeros ? parseFloat(s).toString() : s;
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
	const T = roundSig(convert(stdT_K, 'K', tempUnit), 6, 1e-6, true, tempUnit);
	const P = roundSig(convert(stdP_Pa, 'Pa', pressUnit), 6, 1e-6, true);

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
