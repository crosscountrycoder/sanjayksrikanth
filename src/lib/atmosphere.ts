export const EARTH_RADIUS_M = 6356766;

export const R  = 8.31446261815324; // J/(mol·K)
export const M0 = 0.0289644;        // kg/mol, molar mass of dry air
export const B  = 9.80665 * M0 / R; // g0·M0/R (K/m), from the 1976 Standard Atmosphere

// Layer table: base geopotential altitude (m), base temperature (K), lapse rate (K/m), base pressure (Pa)
export const LAYERS = [
	{ Hb:     0, Tb: 288.15, Lb: -0.0065, Pb: 101325  },
	{ Hb: 11000, Tb: 216.65, Lb:  0,      Pb: 22632.1 },
	{ Hb: 20000, Tb: 216.65, Lb:  0.001,  Pb: 5474.89 },
	{ Hb: 32000, Tb: 228.65, Lb:  0.0028, Pb: 868.019 },
	{ Hb: 47000, Tb: 270.65, Lb:  0,      Pb: 110.906 },
	{ Hb: 51000, Tb: 270.65, Lb: -0.0028, Pb: 66.9389 },
	{ Hb: 71000, Tb: 214.65, Lb: -0.002,  Pb: 3.95642 },
];

export function geometricToGeopotential(h: number): number {
	return (EARTH_RADIUS_M * h) / (EARTH_RADIUS_M + h);
}

export function geopotentialToGeometric(H: number): number {
	return (EARTH_RADIUS_M * H) / (EARTH_RADIUS_M - H);
}

export function getLayerIndex(H: number): number {
	for (let i = LAYERS.length - 1; i >= 0; i--) {
		if (H >= LAYERS[i].Hb) return i;
	}
	return 0;
}

export function getPressure(H: number, P0_Pa: number, T0_K: number = 288.15): number {
	const dT = T0_K - 288.15;
	const layerIdx = getLayerIndex(H);
	let P = P0_Pa;
	for (let i = 0; i <= layerIdx; i++) {
		const { Hb, Lb, Tb } = LAYERS[i];
		const Tbs = Tb + dT;
		const H1 = i < layerIdx ? LAYERS[i + 1].Hb : H;
		if (Lb === 0) {
			P *= Math.exp(-B * (H1 - Hb) / Tbs);
		} else {
			P *= Math.pow((Tbs + Lb * (H1 - Hb)) / Tbs, -B / Lb);
		}
	}
	return P;
}

export function getTemperature(H: number, T0_K: number): number {
	const dT = T0_K - 288.15;
	const { Hb, Tb, Lb } = LAYERS[getLayerIndex(H)];
	return (Tb + dT) + Lb * (H - Hb);
}

export function getSeaLevelTemperature(H: number, T_K: number): number {
	const { Hb, Tb, Lb } = LAYERS[getLayerIndex(H)];
	return T_K - Tb - Lb * (H - Hb) + 288.15;
}

export function getSeaLevelPressure(H: number, P_Pa: number, T0_K: number): number {
	return P_Pa / getPressure(H, 1, T0_K);
}

export function getAltimeterSetting(H: number, P_Pa: number): number {
	return getSeaLevelPressure(H, P_Pa, 288.15);
}

export function getDensity(P_Pa: number, T_K: number): number {
	return (P_Pa * M0) / (R * T_K);
}

// Standard atmosphere values at 86 km geometric (boundary between regimes)
const H86 = geometricToGeopotential(86000);
export const P_AT_86KM   = getPressure(H86, 101325, 288.15);
export const RHO_AT_86KM = getDensity(P_AT_86KM, getTemperature(H86, 288.15));

export function getPressureAltitude(P_Pa: number): number {
	if (P_Pa >= P_AT_86KM) {
		// Lower atmosphere: bisect over geometric -5000 to 86000 m
		let lo = -5000, hi = 86000;
		for (let i = 0; i < 40; i++) {
			const mid = (lo + hi) / 2;
			if (getPressure(geometricToGeopotential(mid), 101325, 288.15) > P_Pa) lo = mid; else hi = mid;
		}
		return hi;
	} else {
		// Upper atmosphere: single ascending scan, stop when pressure crosses target
		return scanUpperForPressure(P_Pa);
	}
}

export function getDensityAltitude(rho: number): number {
	if (rho >= RHO_AT_86KM) {
		// Lower atmosphere: bisect over geometric -5000 to 86000 m
		let lo = -5000, hi = 86000;
		for (let i = 0; i < 40; i++) {
			const mid = (lo + hi) / 2;
			const H   = geometricToGeopotential(mid);
			const T   = getTemperature(H, 288.15);
			const P   = getPressure(H, 101325);
			if (getDensity(P, T) > rho) lo = mid; else hi = mid;
		}
		return hi;
	} else {
		// Upper atmosphere: single ascending scan, stop when density crosses target
		return scanUpperForDensity(rho);
	}
}

// Single-pass upward integration from 86 km, stopping when target is crossed.
function scanUpperForPressure(targetP: number): number {
	return scanUpper((n, T) => n.reduce((sum, ni) => sum + ni, 0) * k_B * T, targetP, P_AT_86KM);
}

function scanUpperForDensity(targetRho: number): number {
	// standard86 must use the same effective scale (P_AT_86KM/P_n7) as getUpperAtmosphere,
	// not RHO_AT_86KM (which uses a different molar mass / temperature).
	return scanUpper((n) => n.reduce((sum, ni, i) => sum + ni * M_[i], 0) / Na, targetRho, P_AT_86KM * rho_n7 / P_n7);
}

function scanUpper(getValue: (n: number[], T: number) => number, target: number, standard86: number): number {
	const STEP = 100;
	// Scale n7 so getValue at Z7 equals the standard 86 km value for this quantity,
	// ensuring the scan is continuous with the lower atmosphere.
	const scale = standard86 / getValue(n7 as number[], T7);
	const n: number[] = n7.map(v => v * scale);
	const fint = [0, 0, 0, 0, 0, 0];
	let T = T7 + dT_std, dT = 0, g = upperGravity(Z7), K = upperEddy(Z7);
	const fprev = [0, 1, 2, 3, 4].map(s => fcalc(Z7, g, K, T, dT, s, n, { val: 0 }));
	let tauint = 0, tauPrev = 0, hLayerStarted = false;
	let Z = Z7;
	let prevVal = getValue(n, T);

	while (Z < Z12) {
		const Znext = Z + STEP;
		const h = STEP;
		const Zmid = Z + h / 2;
		const Tmid = upperTemperature(Zmid) + dT_std, dTmid = upperDT_dZ(Zmid);
		const gmid = upperGravity(Zmid), Kmid = upperEddy(Zmid);
		const Tnext = upperTemperature(Znext) + dT_std, dTnext = upperDT_dZ(Znext);
		const gnext = upperGravity(Znext), Knext = upperEddy(Znext);

		for (let s = 0; s <= 4; s++) {
			const fmid = fcalc(Zmid, gmid, Kmid, Tmid, dTmid, s, n, { val: 0 });
			const fnext = fcalc(Znext, gnext, Knext, Tnext, dTnext, s, n, { val: 0 });
			fint[s] += h * (fprev[s] + 4 * fmid + fnext) / 6;
			fprev[s] = fnext;
		}
		for (let s = 0; s <= 4; s++) n[s] = n7[s] * scale * (T7 / Tnext) * Math.exp(-fint[s]);

		if (Znext >= Zh) {
			if (!hLayerStarted) { tauPrev = (gnext / Tnext) * (M_[5] / Rs); hLayerStarted = true; }
			else { const tauNext = (gnext / Tnext) * (M_[5] / Rs); tauint += h * (tauPrev + tauNext) / 2; tauPrev = tauNext; }
			n[5] = Znext < Z11
				? (n7[5] + Hi - fint[5]) * Math.exp((1 + alpha[5]) * Math.log(T11 / Tnext) + taui - tauint)
				: 0;
		}

		const val = getValue(n, Tnext);
		if (val <= target) {
			// Linearly interpolate between Z and Znext
			return Z + STEP * (prevVal - target) / (prevVal - val);
		}
		prevVal = val;
		T = Tnext; dT = dTnext; g = gnext; K = Knext;
		Z = Znext;
	}
	return Z12; // target not reached within 1000 km
}

// ── Upper atmosphere (86–1000 km geometric) ─────────────────────────────────
// Based on ussa1976.dpr by Steven S. Pietrobon (1999), translated from Pascal.
// All altitudes here are GEOMETRIC (meters).

const k_B  = 1.380622e-23;  // J/K, Boltzmann constant
const Na   = 6.02214076e26;   // 1/kmol, Avogadro constant
const Rs   = 1000 * R;       // J/(kmol·K), gas constant
const g0   = 9.80665;       // m/s²

// Species indices: 0=N2, 1=O, 2=O2, 3=Ar, 4=He, 5=H
const M_: [number,number,number,number,number,number] =
    [28.0134, 15.9994, 31.9988, 39.948, 4.0026, 1.00797]; // kg/kmol

// Boundary conditions at 86 km geometric
const Z7  = 86000;
const T7  = 186.8673;  // K
// Number densities at 86 km (species 0–4) and at 500 km (H, species 5)
const n7: [number,number,number,number,number,number] =
    [1.129794e20, 8.6e16, 3.030898e19, 1.3514e18, 7.58173e14, 8.0e10];
const alpha: [number,number,number,number,number,number] =
    [0.0, 0.0, 0.0, 0.0, -0.40, -0.25]; // thermal-diffusion coefficients

// Pressure and density implied by n7 at T7 (used to normalize scan scales correctly)
const P_n7   = n7.reduce((s, ni) => s + ni, 0) * k_B * T7;
const rho_n7 = n7.reduce((s, ni, i) => s + ni * M_[i], 0) / Na;
// Temperature offset between lower-atm formula and T7 at 86 km (≈ 0.079 K at standard).
// scanUpper uses the same shifted profile so it is the inverse of getUpperAtmosphere at
// standard sea-level conditions, ensuring PA = DA = geometric altitude in the standard atm.
const dT_std = getTemperature(H86, 288.15) - T7;

// H-layer constants
const phi   = 7.2e11;       // 1/(m²·s), vertical flux of H
const taui  = 0.40463343;   // dimensionless, precomputed tau integral at 500 km
const Hi    = 9.8776951e10; // 1/m³, precomputed H number-density integral
const T11   = 999.2356;     // K, kinetic temperature at 500 km
const Z11   = 500000;       // m
const Z12   = 1000000;      // m
const Zh    = 150000;       // m, start of H layer
const Zm    = 100000;       // m, N2 diffusion start
const Zk    = 115000;       // m, end of eddy layer

function upperGravity(Z: number): number {
    return g0 * (EARTH_RADIUS_M / (EARTH_RADIUS_M + Z)) ** 2;
}

function upperEddy(Z: number): number {
    const K7 = 120.0, Ze = 95000, A = 4.00e8;
    if (Z < Ze)       return K7;
    if (Z < Zk - 0.1) return K7 * Math.exp(1 - A / (A - (Z - Ze) ** 2));
    return 0;
}

export function upperTemperature(Z: number): number {
    const T9 = 240.0, T10 = 360.0, Tinf = 1000.0;
    const Z8 = 91000, Z9 = 110000, Z10 = 120000;
    const Tc = 263.1905, A = -76.3232, a_ = 19942.9;
    const Lk9 = 0.012, lambda = 1.875e-5;
    if (Z < Z8)  return T7;
    if (Z < Z9)  { const b = (Z - Z8) / a_; return Tc + A * Math.sqrt(1 - b * b); }
    if (Z < Z10) return T9 + Lk9 * (Z - Z9);
    const b = (EARTH_RADIUS_M + Z10) / (EARTH_RADIUS_M + Z);
    return Tinf - (Tinf - T10) * Math.exp(-lambda * (Z - Z10) * b);
}

function upperDT_dZ(Z: number): number {
    const T9 = 240.0, T10 = 360.0, Tinf = 1000.0;
    const Z8 = 91000, Z9 = 110000, Z10 = 120000;
    const Tc = 263.1905, A = -76.3232, a_ = 19942.9;
    const Lk9 = 0.012, lambda = 1.875e-5;
    if (Z < Z8)  return 0;
    if (Z < Z9)  { const b = (Z - Z8) / a_; return -A * b / (a_ * Math.sqrt(1 - b * b)); }
    if (Z < Z10) return Lk9;
    const b = (EARTH_RADIUS_M + Z10) / (EARTH_RADIUS_M + Z);
    const B = (Tinf - T10) * Math.exp(-lambda * (Z - Z10) * b);
    return lambda * b * b * B;
}

function upperFlux(Z: number, s: number): number {
    // s: species index (1=O, 2=O2, 3=Ar, 4=He only)
    const Q = [-5.809644e-13, 1.366212e-13, 9.434079e-14, -2.457369e-13];
    const U = [56903.11, 86000, 86000, 86000];
    const W = [2.706240e-14, 8.333333e-14, 8.333333e-14, 6.666667e-13];
    const si = s - 1; // 0-based index into Q/U/W (O=0, O2=1, Ar=2, He=3)
    const a_ = Z - U[si], Asq = a_ * a_;
    let x = Q[si] * Asq * Math.exp(-W[si] * Asq * a_);
    if (s === 1 && Z < 97000) {
        const qO = -3.416248e-12, uO = 97000, wO = 5.008765e-13;
        const a2 = uO - Z, A2 = a2 * a2;
        x += qO * A2 * Math.exp(-wO * A2 * a2);
    }
    return x;
}

// Diffusion coefficient constants per species (O, O2, Ar, He, H)
const diffA = [6.986e20, 4.863e20, 4.487e20, 1.700e21, 3.305e21];
const diffB = [0.750,    0.750,    0.870,    0.691,    0.500   ];

function molecularDiffusion(Z: number, T: number, s: number, n: number[]): number {
    // background gas sum: N2 for O/O2, up to O2 for Ar/He, up to He for H
    const maxS = s <= 2 ? 0 : s <= 4 ? 2 : 4;
    let sumn = 0;
    for (let i = 0; i <= maxS; i++) sumn += n[i];
    const si = s - 1; // O=0, O2=1, Ar=2, He=3, H=4
    return diffA[si] * (T / 273.15) ** diffB[si] / sumn;
}

function meanMolarMass(Z: number, n: number[]): number {
    if (Z < Zm) return M0 * 1000; // kg/kmol
    const maxS = 4; // N2..He
    let sumn = 0, sumnm = 0;
    for (let i = 0; i <= maxS; i++) { sumn += n[i]; sumnm += n[i] * M_[i]; }
    return sumnm / sumn;
}

// dnᵢ/dZ integrand for each species (the "fcalc" in the Pascal source)
function fcalc(Z: number, g: number, K: number, T: number, dT: number, s: number, n: number[], tau: { val: number }): number {
    if (s === 0) { // N2
        const M = Z < Zm ? M0 * 1000 : M_[0];
        return M * g / (Rs * T);
    }
    if (s >= 1 && s <= 4) { // O, O2, Ar, He
        let f: number;
        if (Z < Zk - 0.1) {
            const M = meanMolarMass(Z, n);
            const D = molecularDiffusion(Z, T, s, n);
            f = D * (M_[s] + M * K / D + alpha[s] * Rs * dT / g) / (D + K);
        } else {
            f = M_[s] + alpha[s] * Rs * dT / g;
        }
        return g * f / (Rs * T) + upperFlux(Z, s);
    }
    // H (s === 5)
    return 0; // H handled via tau integration separately
}

export interface UpperAtmosphereResult {
    T: number;                   // kinetic temperature, K
    n: [number,number,number,number,number,number]; // number densities, 1/m³
    P: number;                   // total pressure, Pa
    partialP: [number,number,number,number,number,number]; // partial pressures, Pa
    rho: number;                 // mass density, kg/m³
    meanM: number;               // mean molar mass, kg/mol
}

// Integrate from 86 km to target altitude Z using Simpson's rule (10 m steps).
// T_bc and P_bc are the temperature and pressure at 86 km from the lower atmosphere;
// they default to the standard model values. Passing non-standard values ensures
// continuity at the 86 km boundary when sea-level conditions differ from standard.
export function getUpperAtmosphere(
    Ztarget: number,
    T_bc: number = T7,
    P_bc: number = P_AT_86KM,
): UpperAtmosphereResult {
    const STEP = 100;

    // Temperature offset: shift the entire profile so T(86 km) = T_bc
    const dT_offset = T_bc - T7;
    const tempAt = (Z: number) => upperTemperature(Z) + dT_offset;
    const dtempAt = (Z: number) => upperDT_dZ(Z); // derivative unchanged by constant shift

    // Scale all N2..He number densities so that n_total * k_B * T_bc = P_bc at Z7.
    // The 1976 n7 values give P_n7 = n7_total*k_B*T7, which differs slightly from
    // P_AT_86KM (lower-atm formula). We normalize against P_n7 so that the initial
    // pressure exactly equals P_bc, ensuring continuity at the 86 km boundary.
    const n_scale = (P_bc / T_bc) / (P_n7 / T7);

    const n: [number,number,number,number,number,number] = n7.map(v => v * n_scale) as any;
    let T = T_bc, dT = 0, g = upperGravity(Z7), K = upperEddy(Z7);

    // fint[s] accumulates ∫ fcalc dZ for species 0–4
    const fint = [0, 0, 0, 0, 0, 0];

    // Previous fcalc values (for Simpson's rule: need f at Z, Z+h/2, Z+h)
    const fprev = [0, 1, 2, 3, 4].map(s => fcalc(Z7, g, K, T, dT, s, n, { val: 0 }));

    // Tau integration for H
    let tauint = 0;
    let tauPrev = 0;
    let hLayerStarted = false;

    let Z = Z7;

    while (Z < Ztarget - 1e-9) {
        const Znext = Math.min(Z + STEP, Ztarget);
        const h = Znext - Z;

        // Midpoint
        const Zmid = Z + h / 2;
        const Tmid = tempAt(Zmid);
        const dTmid = dtempAt(Zmid);
        const gmid = upperGravity(Zmid);
        const Kmid = upperEddy(Zmid);

        // End point
        const Tnext = tempAt(Znext);
        const dTnext = dtempAt(Znext);
        const gnext = upperGravity(Znext);
        const Knext = upperEddy(Znext);

        // Simpson's rule for species 0–4
        for (let s = 0; s <= 4; s++) {
            const fmid = fcalc(Zmid, gmid, Kmid, Tmid, dTmid, s, n, { val: 0 });
            const fnext = fcalc(Znext, gnext, Knext, Tnext, dTnext, s, n, { val: 0 });
            fint[s] += h * (fprev[s] + 4 * fmid + fnext) / 6;
            fprev[s] = fnext;
        }

        // Update number densities for N2..He
        for (let s = 0; s <= 4; s++) {
            n[s] = n7[s] * n_scale * (T_bc / Tnext) * Math.exp(-fint[s]);
        }

        // H: tau integration (only above Zh)
        if (Znext >= Zh) {
            if (!hLayerStarted) {
                tauPrev = (gnext / Tnext) * (M_[5] / Rs);
                hLayerStarted = true;
            } else {
                const tauNext = (gnext / Tnext) * (M_[5] / Rs);
                tauint += h * (tauPrev + tauNext) / 2;
                tauPrev = tauNext;
            }
            if (Znext < Z11) {
                n[5] = (n7[5] + Hi - fint[5]) * Math.exp((1 + alpha[5]) * Math.log(T11 / Tnext) + taui - tauint);
            } else {
                n[5] = 0;
            }
        }

        T = Tnext;
        dT = dTnext;
        g = gnext;
        K = Knext;
        Z = Znext;
    }

    // Derived quantities
    const partialP = n.map(ni => ni * k_B * T) as [number,number,number,number,number,number];
    const P    = partialP.reduce((a, b) => a + b, 0);
    const rho  = n.reduce((a, ni, i) => a + ni * M_[i], 0) / Na;
    const nTot = n.reduce((a, b) => a + b, 0);
    const meanM = nTot > 0 ? (rho * Na / nTot) / 1000 : M0; // kg/mol

    return { T, n: [...n] as any, P, partialP, rho, meanM };
}

export const SPECIES_NAMES = ['N₂', 'O', 'O₂', 'Ar', 'He', 'H'] as const;

// Speed of sound: a = √(γ · P / ρ), γ = 1.4 for dry air
export function getSpeedOfSound(P_Pa: number, rho: number): number {
    return Math.sqrt(1.4 * P_Pa / rho);
}

// Mean free path: λ = 1 / (√2 · π · σ² · n), n = P / (k_B_SI · T)
// σ from U.S. Standard Atmosphere 1976 Table 2 (effective collision diameter of air)
const SIGMA = 3.65e-10;        // m
const K_B_SI = 1.380649e-23;   // J/K, exact SI value (2019 redefinition)
export function getMeanFreePath(P_Pa: number, T_K: number): number {
    const n = P_Pa / (K_B_SI * T_K);
    return 1 / (Math.SQRT2 * Math.PI * SIGMA * SIGMA * n);
}

// Mole fractions of gases in the lower atmosphere (below 86 km).
// O, He, H are negligible and treated as 0.
export const LOWER_MOLE_FRACTIONS: Record<string, number> = {
    N2: 0.78084, O2: 0.20946, Ar: 0.00934, CO2: 0.00033, O: 0, He: 0, H: 0,
};

// Maps species key to index in upper-atmosphere n[] / partialP[] arrays.
// CO2 is not tracked in the upper atmosphere; it maps to -1.
export const UPPER_SPECIES_INDEX: Record<string, number> = {
    N2: 0, O: 1, O2: 2, Ar: 3, He: 4, H: 5, CO2: -1,
};

export interface AtmosphereSnapshot {
    T: number;
    P: number;
    rho: number;
    partialP: [number, number, number, number, number, number]; // N2, O, O2, Ar, He, H
}

/**
 * Integrates the upper atmosphere ODE from 86 km to Zmax in 100 m steps,
 * storing {T, P, rho} at each step. Returns a Map keyed by Math.round(Z).
 * Includes Z = 86000, 86100, …, and the exact Zmax (final step may be < 100 m).
 * Use buildUpperTable once per graph render and look up individual points in O(1).
 */
export function buildUpperTable(
    Zmax: number,
    T_bc: number,
    P_bc: number,
): Map<number, AtmosphereSnapshot> {
    const table = new Map<number, AtmosphereSnapshot>();
    const STEP = 100;

    const dT_offset = T_bc - T7;
    const tempAt  = (Z: number) => upperTemperature(Z) + dT_offset;
    const n_scale = (P_bc / T_bc) / (P_n7 / T7);

    const n: [number,number,number,number,number,number] = n7.map(v => v * n_scale) as any;
    let T = T_bc, dT = 0;
    const fint = [0, 0, 0, 0, 0, 0];
    const fprev = [0, 1, 2, 3, 4].map(s =>
        fcalc(Z7, upperGravity(Z7), upperEddy(Z7), T, dT, s, n, { val: 0 })
    );
    let tauint = 0, tauPrev = 0, hLayerStarted = false;

    const storeSnap = (Zk: number) => {
        const P        = n.reduce((s, ni) => s + ni * k_B * T, 0);
        const rho      = n.reduce((s, ni, i) => s + ni * M_[i], 0) / Na;
        const partialP = n.map(ni => ni * k_B * T) as [number,number,number,number,number,number];
        table.set(Math.round(Zk), { T, P, rho, partialP });
    };
    storeSnap(Z7);

    let Z = Z7;
    while (Z < Zmax - 1e-9) {
        const h     = Math.min(STEP, Zmax - Z);
        const Znext = Z + h;
        const Zmid  = Z + h / 2;

        const Tmid  = tempAt(Zmid),  dTmid  = upperDT_dZ(Zmid);
        const gmid  = upperGravity(Zmid), Kmid = upperEddy(Zmid);
        const Tnext = tempAt(Znext), dTnext = upperDT_dZ(Znext);
        const gnext = upperGravity(Znext), Knext = upperEddy(Znext);

        for (let s = 0; s <= 4; s++) {
            const fmid  = fcalc(Zmid,  gmid,  Kmid,  Tmid,  dTmid,  s, n, { val: 0 });
            const fnext = fcalc(Znext, gnext, Knext, Tnext, dTnext, s, n, { val: 0 });
            fint[s] += h * (fprev[s] + 4 * fmid + fnext) / 6;
            fprev[s] = fnext;
        }
        for (let s = 0; s <= 4; s++) n[s] = n7[s] * n_scale * (T_bc / Tnext) * Math.exp(-fint[s]);

        if (Znext >= Zh) {
            if (!hLayerStarted) {
                tauPrev = (gnext / Tnext) * (M_[5] / Rs);
                hLayerStarted = true;
            } else {
                const tauNext = (gnext / Tnext) * (M_[5] / Rs);
                tauint += h * (tauPrev + tauNext) / 2;
                tauPrev = tauNext;
            }
            n[5] = Znext < Z11
                ? (n7[5] + Hi - fint[5]) * Math.exp((1 + alpha[5]) * Math.log(T11 / Tnext) + taui - tauint)
                : 0;
        }

        T = Tnext; dT = dTnext; Z = Znext;
        storeSnap(Znext);
    }

    return table;
}

export function waterVaporPressure(T: number): number {
	const Tc = 647.096, Pc = 22064000;
	const theta = 1 - T / Tc;
	const exponent = (Tc / T) * (
		-7.85951783 * theta     + 1.84408259 * theta**1.5 +
		-11.7866497 * theta**3  + 22.6807411 * theta**3.5 +
		-15.9618719 * theta**4  + 1.80122502 * theta**7.5
	);
	return Pc * Math.exp(exponent);
}

export function getBoilingPoint(P_Pa: number): number {
	let lo = 273, hi = 647.096;
	for (let i = 0; i < 30; i++) {
		const mid = (lo + hi) / 2;
		if (waterVaporPressure(mid) < P_Pa) lo = mid; else hi = mid;
	}
	return hi;
}
