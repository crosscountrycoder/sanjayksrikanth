// Shared helpers for the SVG graphs (standard atmosphere and dew point): tick placement,
// tick-label rounding, and the color palette for the current light/dark theme.

// Returns n+1 ticks evenly spaced between min and max.
export function linearTicks(min: number, max: number, n: number): number[] {
	const ticks: number[] = [];
	for (let i = 0; i <= n; i++)
		ticks.push(parseFloat((min + (max - min) * i / n).toPrecision(10)));
	return ticks;
}

// Returns n+1 ticks evenly spaced in log10 space between lo and hi (both must be > 0).
export function logTicks(lo: number, hi: number, n: number): number[] {
	const logLo = Math.log10(lo), logHi = Math.log10(hi);
	const ticks: number[] = [];
	for (let i = 0; i <= n; i++)
		ticks.push(Math.pow(10, logLo + (logHi - logLo) * i / n));
	return ticks;
}

// Strips trailing zeros from a rounded number's string without changing its notation
// (Number(...).toString() would flatten a large exponential value back to plain).
function stripZerosKeepNotation(s: string): string {
	const eIdx = s.indexOf('e');
	const mantissa = eIdx === -1 ? s : s.slice(0, eIdx);
	const suffix = eIdx === -1 ? '' : s.slice(eIdx);
	const trimmed = mantissa.includes('.') ? mantissa.replace(/0+$/, '').replace(/\.$/, '') : mantissa;
	return trimmed + suffix;
}

// Tick labels use coarser rounding than the calculators' outputs, with trailing zeros removed:
//   temperatures (incl. dew/boiling point)  formatTemp(v, 1, true) in helpers.ts: 1 digit after the decimal point
//   altitudes                               formatAltitude in helpers.ts
//   everything else                         formatTick: 4 significant figures

// 4 sig figs, scientific notation below 1e-3 or once the value needs more than 4 digits
// before the decimal point, with trailing zeros stripped without disturbing that notation.
export function formatTick(v: number): string {
	if (v === 0) return '0';
	const rounded = Number(v.toPrecision(4));
	return stripZerosKeepNotation(Math.abs(rounded) < 1e-3 ? rounded.toExponential(3) : rounded.toPrecision(4));
}

export function graphColors() {
	const dark = document.documentElement.dataset.theme === 'dark';
	return {
		bg:      dark ? '#000000' : '#ffffff',
		plotBg:  dark ? '#111111' : '#f8f8f8',
		grid:    dark ? '#333333' : '#e0e0e0',
		axis:    dark ? '#cccccc' : '#333333',
		text:    dark ? '#cccccc' : '#333333',
		line:    '#0000ff',
	};
}

// Saves an <svg> element as a standalone .svg file.
export function downloadSvg(svg: SVGElement, filename: string) {
	const str = '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg);
	const url = URL.createObjectURL(new Blob([str], { type: 'image/svg+xml' }));
	const a = Object.assign(document.createElement('a'), { href: url, download: filename });
	document.body.appendChild(a);
	a.click();
	document.body.removeChild(a);
	URL.revokeObjectURL(url);
}
