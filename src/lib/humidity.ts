// Shared by the dew point graph and table: the six humidity variables, their units, and
// solving for the air state (temperature and dew point) from any two of them.

import * as wp from './water-properties.ts';
import { convert } from './convert.ts';
import { formatTemp, formatSig } from './helpers.ts';

// Water's critical point: the upper limit for temperatures and vapor pressures.
export const T_CRIT_K  = 647.096;
export const P_CRIT_PA = 22064000;

export type Key = 'temp' | 'rh' | 'dew' | 'vp' | 'svp' | 'ah';
export const KEYS: Key[] = ['temp', 'rh', 'dew', 'vp', 'svp', 'ah'];
export const NAMES: Record<Key, string> = {
	temp: 'Temperature', rh: 'Relative humidity', dew: 'Dew point',
	vp: 'Vapor pressure', svp: 'Saturation vapor pressure', ah: 'Absolute humidity',
};

// Each variable's kind sets its units, its SI unit, and how its values are rounded.
// Relative humidity is shown in % and stored as a fraction.
export type Kind = 'temp' | 'press' | 'rh' | 'dens';
export const KIND: Record<Key, Kind> = { temp: 'temp', dew: 'temp', vp: 'press', svp: 'press', rh: 'rh', ah: 'dens' };
const SI_UNIT: Record<Kind, string> = { temp: 'K', press: 'Pa', rh: '%', dens: 'kg/m3' };
const UNITS: Record<Kind, { value: string; label: string }[]> = {
	temp:  [{ value: 'C', label: '°C' }, { value: 'F', label: '°F' }, { value: 'K', label: 'K' }, { value: 'R', label: '°R' }],
	press: [
		{ value: 'atm', label: 'atm' }, { value: 'Pa', label: 'Pa' }, { value: 'hPa', label: 'hPa (mbar)' },
		{ value: 'kPa', label: 'kPa' }, { value: 'mmHg', label: 'mmHg (Torr)' }, { value: 'inHg', label: 'inHg' },
		{ value: 'psi', label: 'psi' },
	],
	rh:    [{ value: '%', label: '%' }],
	dens:  [{ value: 'kg/m3', label: 'kg/m³' }, { value: 'g/m3', label: 'g/m³' }, { value: 'lb/ft3', label: 'lb/ft³' }, { value: 'gr/ft3', label: 'grains/ft³' }],
};
const DEFAULT_UNIT: Record<Kind, string> = { temp: 'C', press: 'hPa', rh: '%', dens: 'g/m3' };

// Short unit labels for headings and messages.
const UNIT_LABELS: Record<string, string> = {
	C: '°C', F: '°F', K: 'K', R: '°R', 'kg/m3': 'kg/m³', 'g/m3': 'g/m³', 'lb/ft3': 'lb/ft³', 'gr/ft3': 'grains/ft³',
};
export const unitLabel = (u: string) => UNIT_LABELS[u] ?? u;
// "15 °C", "50%": no space before a percent sign.
export const withUnit = (v: string, u: string) => u === '%' ? `${v}%` : `${v} ${unitLabel(u)}`;
// "Dew point (°C)"
export const varLabel = (key: Key, u: string) => `${NAMES[key]} (${unitLabel(u)})`;

// Default range and fixed value for each variable, in SI units.
export const DEFAULT_RANGE: Record<Key, [number, number]> = {
	temp: [273.15, 313.15], dew: [263.15, 303.15], rh: [0, 1],
	vp: [0, 5000], svp: [0, 7500], ah: [0, 0.03],
};
export const DEFAULT_FIXED: Record<Key, number> = {
	temp: 298.15, dew: 273.15, rh: 0.5, vp: 1000, svp: 3000, ah: 0.01,
};

// Upper limit of each kind in SI units (the lower limit is always 0).
const MAX_SI: Record<Kind, number> = { temp: T_CRIT_K, press: P_CRIT_PA, rh: Infinity, dens: Infinity };

// Temperature and saturation vapor pressure determine each other, as do dew point and vapor pressure.
export const PARTNER: Partial<Record<Key, Key>> = { temp: 'svp', svp: 'temp', dew: 'vp', vp: 'dew' };

// Variables that can be held fixed while relating y to x: any variable other than x, y, or the
// one-to-one partner of either. A one-to-one pair needs no fixed variable.
export function fixedOptions(x: Key, y: Key): Key[] {
	if (PARTNER[x] === y) return [];
	return KEYS.filter(k => k !== x && k !== y && k !== PARTNER[x] && k !== PARTNER[y]);
}

export const toSI   = (key: Key, v: number, u: string) => KIND[key] === 'rh' ? v / 100 : convert(v, u, SI_UNIT[KIND[key]]);
export const fromSI = (key: Key, v: number, u: string) => KIND[key] === 'rh' ? v * 100 : convert(v, SI_UNIT[KIND[key]], u);
// The unit for a difference between two values in unit u (e.g. a contour interval): the
// "-change" form of °C, °F, or °R, which has no offset; every other unit is unchanged.
export const changeUnit = (u: string) => u === 'C' || u === 'F' || u === 'R' ? `${u}-change` : u;
// Rounds a value for writing into an input box or a message.
export const fmtInput = (key: Key, v: number) => KIND[key] === 'temp' ? formatTemp(v, 3, true) : formatSig(v, 6, true);
// Rounds a value for an output (table cell), keeping trailing zeros.
export const fmtOutput = (key: Key, v: number) => KIND[key] === 'temp' ? formatTemp(v) : formatSig(v);

// Checks a value (SI units) against its variable's limits; returns an error message or ''.
export function rangeError(key: Key, v: number, u: string, name: string): string {
	const kind = KIND[key];
	if (kind === 'temp' || kind === 'press') {
		if (v >= 0 && v <= MAX_SI[kind]) return '';
		const lo = fmtInput(key, fromSI(key, 0, u)), hi = fmtInput(key, fromSI(key, MAX_SI[kind], u));
		return `${name} must be between ${lo} and ${withUnit(hi, u)}.`;
	}
	return v >= 0 ? '' : `${name} must be at least ${withUnit('0', u)}.`;
}

// ── Physics ──────────────────────────────────────────────────────────────────

export interface State { T: number; Td: number } // air temperature and dew point, K

// Finds the air temperature and dew point from two known variables (SI units), or just
// one when it fixes both (no fixed variable). Returns null when no valid state exists.
export function solveState(known: Partial<Record<Key, number>>): State | null {
	const rh = known.rh, ah = known.ah;
	let T  = known.temp ?? (known.svp !== undefined ? wp.getBoilingPoint(known.svp) : NaN);
	let Td = known.dew  ?? (known.vp  !== undefined ? wp.getBoilingPoint(known.vp)  : NaN);

	if (rh === undefined && ah === undefined) {
		// A single variable from a one-to-one pair: saturated air (dew point = temperature).
		if (isNaN(T)) T = Td; else if (isNaN(Td)) Td = T;
	} else if (isNaN(T) && isNaN(Td)) {
		// Relative and absolute humidity: absolute humidity at a fixed relative humidity
		// increases with temperature, so solve for the temperature by bisection.
		if (rh === undefined || ah === undefined || !(rh > 0)) return null;
		const ahAt = (t: number) => wp.getAbsoluteHumidity(rh * wp.waterVaporPressure(t), t);
		if (ahAt(T_CRIT_K) < ah) return null;
		let lo = 0, hi = T_CRIT_K;
		while (hi - lo > 1e-9) {
			const mid = (lo + hi) / 2;
			ahAt(mid) < ah ? lo = mid : hi = mid;
		}
		T = lo;
	}

	if (!isNaN(T) && isNaN(Td)) {
		const e = rh !== undefined ? rh * wp.waterVaporPressure(T) : wp.getVaporPressureFromAbsoluteHumidity(ah!, T);
		Td = e > P_CRIT_PA ? NaN : wp.getBoilingPoint(e);
	} else if (isNaN(T) && !isNaN(Td)) {
		const e = wp.waterVaporPressure(Td);
		if (rh !== undefined) {
			const es = e / rh;
			T = es > P_CRIT_PA ? NaN : wp.getBoilingPoint(es);
		} else {
			// ρ = e·M/(R·T), so T = (e·M/R) / ρ, where e·M/R is the absolute humidity at 1 K.
			T = wp.getAbsoluteHumidity(e, 1) / ah!;
		}
	}

	if (!(T >= 0 && T <= T_CRIT_K && Td >= 0 && Td <= T_CRIT_K)) return null;
	return { T, Td };
}

// A variable's value (SI units) in a given state.
export function valueOf(key: Key, s: State): number {
	switch (key) {
		case 'temp': return s.T;
		case 'dew':  return s.Td;
		case 'vp':   return wp.waterVaporPressure(s.Td);
		case 'svp':  return wp.waterVaporPressure(s.T);
		case 'rh':   return wp.getRelativeHumidity(s.T, s.Td);
		case 'ah':   return wp.getAbsoluteHumidity(wp.waterVaporPressure(s.Td), s.T);
	}
}

// ── Page controls ────────────────────────────────────────────────────────────

// Fills a unit dropdown for a variable, keeping the current unit if it still applies.
export function fillUnits(select: HTMLSelectElement, key: Key) {
	const options = UNITS[KIND[key]];
	const keep = options.some(o => o.value === select.value) ? select.value : DEFAULT_UNIT[KIND[key]];
	select.innerHTML = options.map(o => `<option value="${o.value}">${o.label}</option>`).join('');
	select.value = keep;
	select.dataset.prev = keep;
}

// Fills a variable dropdown with the given variables, keeping the current choice if it's
// still listed and otherwise choosing `fallback` (or the first listed variable).
export function fillVarOptions(select: HTMLSelectElement, keys: Key[], fallback?: Key | null) {
	const prev = select.value as Key;
	select.innerHTML = keys.map(k => `<option value="${k}">${NAMES[k]}</option>`).join('');
	select.value = keys.includes(prev) ? prev : fallback && keys.includes(fallback) ? fallback : keys[0];
}

// Converts the values in a set of inputs when their shared unit dropdown changes, then calls
// `onChange`. `deltas` hold differences (e.g. contour intervals) rather than values.
export function convertOnUnitChange(
	select: HTMLSelectElement, inputs: HTMLInputElement[], key: () => Key | null,
	onChange: () => void, deltas: HTMLInputElement[] = [],
) {
	select.addEventListener('change', () => {
		const k = key();
		if (k) {
			const from = select.dataset.prev!, to = select.value;
			for (const el of inputs) {
				const v = parseFloat(el.value);
				if (!isNaN(v)) el.value = fmtInput(k, fromSI(k, toSI(k, v, from), to));
			}
			for (const el of deltas) {
				const d = parseFloat(el.value);
				if (!isNaN(d)) el.value = formatSig(fromSI(k, toSI(k, d, changeUnit(from)), changeUnit(to)), 6, true);
			}
		}
		select.dataset.prev = select.value;
		onChange();
	});
}
