// Altitude
export function altToM(v: number, unit: string): number {
	return unit === 'ft' ? v * 0.3048
	     : unit === 'mi' ? v * 1609.344
	     : unit === 'km' ? v * 1000
	     : v;
}
export function altFromM(m: number, unit: string): number {
	return unit === 'ft' ? m / 0.3048
	     : unit === 'mi' ? m / 1609.344
	     : unit === 'km' ? m / 1000
	     : m;
}

// Temperature
export function tempToK(v: number, unit: string): number {
	return unit === 'C' ? v + 273.15
	     : unit === 'F' ? (v - 32) / 1.8 + 273.15
	     : unit === 'R' ? v / 1.8
	     : v;
}
export function tempFromK(K: number, unit: string): number {
	return unit === 'C' ? K - 273.15
	     : unit === 'F' ? (K - 273.15) * 1.8 + 32
	     : unit === 'R' ? K * 1.8
	     : K;
}

// Pressure
export function pressureToPa(v: number, unit: string): number {
	return unit === 'atm'  ? v * 101325
	     : unit === 'hPa'  ? v * 100
	     : unit === 'kPa'  ? v * 1000
	     : unit === 'mmHg' ? v * 133.322387415
	     : unit === 'inHg' ? v * 3386.388640341
	     : unit === 'psi'  ? v * 6894.757293168
	     : v;
}
export function pressureFromPa(Pa: number, unit: string): number {
	return unit === 'atm'  ? Pa / 101325
	     : unit === 'hPa'  ? Pa / 100
	     : unit === 'kPa'  ? Pa / 1000
	     : unit === 'mmHg' ? Pa / 133.322387415
	     : unit === 'inHg' ? Pa / 3386.388640341
	     : unit === 'psi'  ? Pa / 6894.757293168
	     : Pa;
}

// Density (from kg/m³)
export function densityFromKgM3(rho: number, unit: string): number {
	return unit === 'lb/ft3'   ? rho * 0.062427960576145
	     : unit === 'slug/ft3' ? rho * 0.0019403203
	     : unit === 'sigma'    ? rho / 1.225034902
	     : rho;
}

// Speed (from m/s)
export function speedFromMS(v: number, unit: string): number {
	return unit === 'km/h' ? v * 3.6
	     : unit === 'mph'  ? v * 2.2369362920544
	     : unit === 'ft/s' ? v / 0.3048
	     : unit === 'kn'   ? v * 1.9438444924406
	     : v;
}

// Mean free path (from m)
export function mfpFromM(m: number, unit: string): number {
	return unit === 'nm' ? m * 1e9
	     : unit === 'μm' ? m * 1e6
	     : unit === 'mm' ? m * 1000
	     : unit === 'km' ? m / 1000
	     : unit === 'in' ? m / 0.0254
	     : unit === 'ft' ? m / 0.3048
	     : unit === 'mi' ? m / 1609.344
	     : m;
}

// Dynamic viscosity (from Pa·s)
export function viscosityFromPaS(mu: number, unit: string): number {
	return unit === 'μPa·s'    ? mu * 1e6
	     : unit === 'P'         ? mu * 10
	     : unit === 'cP'        ? mu * 1000
	     : unit === 'lbf·s/ft2' ? mu * 0.0208854342331501
	     : mu;
}
