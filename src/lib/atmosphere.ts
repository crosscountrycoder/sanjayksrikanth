export const EARTH_RADIUS_M = 6356766;
export const R = 8.31446261815324; // J/(mol·K)

// ── Lower atmosphere composition ─────────────────────────────────────────────
// Mole fractions sum to exactly 1. Molar masses in kg/mol.
// Sources: NOAA GML (CO2, CH4, N2O, 2024 annual means); USSA 1976 (all others).
export const LOWER_COMPOSITION = [
    { key: 'N2',  xi: 0.780859353, Mi: 0.02801340 },
    { key: 'O2',  xi: 0.209351263, Mi: 0.03199880 },
    { key: 'Ar',  xi: 0.009340000, Mi: 0.03994800 },
    { key: 'CO2', xi: 0.000422000, Mi: 0.04400950 },
    { key: 'Ne',  xi: 0.000018180, Mi: 0.02017970 },
    { key: 'He',  xi: 0.000005240, Mi: 0.00400260 },
    { key: 'CH4', xi: 0.000001900, Mi: 0.01604246 },
    { key: 'Kr',  xi: 0.000001140, Mi: 0.08379800 },
    { key: 'H2',  xi: 0.000000500, Mi: 0.00201588 },
    { key: 'N2O', xi: 0.000000337, Mi: 0.04401280 },
    { key: 'Xe',  xi: 0.000000087, Mi: 0.13129300 },
] as const;

export const M0 = LOWER_COMPOSITION.reduce((s, c) => s + c.xi * c.Mi, 0);

// ── Coordinate conversion ────────────────────────────────────────────────────

export function geometricToGeopotential(z: number): number {
    return (EARTH_RADIUS_M * z) / (EARTH_RADIUS_M + z);
}
export function geopotentialToGeometric(H: number): number {
    return (EARTH_RADIUS_M * H) / (EARTH_RADIUS_M - H);
}

// ── Temperature ──────────────────────────────────────────────────────────────
// Unified across all altitudes. Below 86 km: piecewise linear in geometric z,
// with layer boundaries from USSA 1976 geopotential values converted to geometric
// and rounded to 100 m. Above 86 km: USSA 1976 molecular-temperature formula.
// Extrapolates below 0 m with the tropospheric lapse rate.

const TEMP_LAYERS = [
    { z:     0, T: 288.15   },
    { z: 11000, T: 216.65   },
    { z: 20100, T: 216.65   },
    { z: 32200, T: 228.65   },
    { z: 47400, T: 270.65   },
    { z: 51400, T: 270.65   },
    { z: 71800, T: 214.65   },
    { z: 86000, T: 186.8673 },
] as const;

function _upperT(Z: number): number {
    const Z8 = 91000, Z9 = 110000, Z10 = 120000;
    const Tc = 263.1905, A = -76.3232, a_ = 19942.9;
    const T9 = 240.0, T10 = 360.0, Tinf = 1000.0, Lk9 = 0.012, lambda = 1.875e-5;
    if (Z < Z8) return 186.8673;
    if (Z < Z9) { const b = (Z - Z8) / a_; return Tc + A * Math.sqrt(1 - b * b); }
    if (Z < Z10) return T9 + Lk9 * (Z - Z9);
    const b = (EARTH_RADIUS_M + Z10) / (EARTH_RADIUS_M + Z);
    return Tinf - (Tinf - T10) * Math.exp(-lambda * (Z - Z10) * b);
}

function _upperDTdZ(Z: number): number {
    const Z8 = 91000, Z9 = 110000, Z10 = 120000;
    const A = -76.3232, a_ = 19942.9;
    const T10 = 360.0, Tinf = 1000.0, Lk9 = 0.012, lambda = 1.875e-5;
    if (Z < Z8) return 0;
    if (Z < Z9) { const b = (Z - Z8) / a_; return -A * b / (a_ * Math.sqrt(1 - b * b)); }
    if (Z < Z10) return Lk9;
    const b = (EARTH_RADIUS_M + Z10) / (EARTH_RADIUS_M + Z);
    return lambda * b * b * (Tinf - T10) * Math.exp(-lambda * (Z - Z10) * b);
}

export function getTemperature(z: number): number {
    if (z > 86000) return _upperT(z);
    const L = TEMP_LAYERS;
    if (z <= L[0].z) {
        const dTdz = (L[1].T - L[0].T) / (L[1].z - L[0].z);
        return L[0].T + dTdz * (z - L[0].z);
    }
    for (let i = 0; i < L.length - 1; i++) {
        if (z <= L[i + 1].z) {
            const t = (z - L[i].z) / (L[i + 1].z - L[i].z);
            return L[i].T + t * (L[i + 1].T - L[i].T);
        }
    }
    return L[L.length - 1].T;
}

// ── Lower-atmosphere composition ─────────────────────────────────────────────

const H_TROP      = 11000;
const L_CH4       = 17000;
const L_N2O       = 11000;
const H_O_START   = 80000;
const X_O_AT_86KM = 5.942e-4;
const MI_O        = 0.0159994;
const H_86KM = EARTH_RADIUS_M * 86000 / (EARTH_RADIUS_M + 86000);
const _xi_CH4 = LOWER_COMPOSITION[6].xi;
const _xi_N2O = LOWER_COMPOSITION[9].xi;

function _lowerXi(H: number): Record<string, number> {
    const fracs: Record<string, number> = {};
    for (const c of LOWER_COMPOSITION) fracs[c.key] = c.xi;
    fracs['O'] = 0; fracs['H'] = 0;
    if (H > H_TROP) {
        const dH = H - H_TROP;
        const fCH4 = Math.exp(-dH / L_CH4);
        const fN2O = Math.exp(-dH / L_N2O);
        fracs['CH4'] = _xi_CH4 * fCH4;
        fracs['N2O'] = _xi_N2O * fN2O;
        fracs['N2']  = LOWER_COMPOSITION[0].xi + _xi_CH4 * (1 - fCH4) + _xi_N2O * (1 - fN2O);
    }
    if (H >= H_O_START) {
        const t = Math.min((H - H_O_START) / (H_86KM - H_O_START), 1.0);
        const smooth = t * t * (3 - 2 * t);
        fracs['O']  = X_O_AT_86KM * smooth;
        fracs['O2'] = LOWER_COMPOSITION[1].xi - fracs['O'] / 2;
    }
    let total = 0;
    for (const k of Object.keys(fracs)) total += fracs[k];
    const inv = 1 / total;
    for (const k of Object.keys(fracs)) fracs[k] *= inv;
    return fracs;
}

function _lowerMeff(H: number): number {
    const fracs = _lowerXi(H);
    let M = 0;
    for (const c of LOWER_COMPOSITION) M += fracs[c.key] * c.Mi;
    M += fracs['O'] * MI_O;
    return M;
}

// ── Lower-atmosphere pressure (RK4 hydrostatic) ──────────────────────────────

function _g(z: number): number {
    return 9.80665 * (EARTH_RADIUS_M / (EARTH_RADIUS_M + z)) ** 2;
}

function _hydroF(z: number, dT: number): number {
    return _lowerMeff(geometricToGeopotential(z)) * _g(z) / (R * (getTemperature(z) + dT));
}

function _rk4P(z: number, P: number, h: number, dT: number): number {
    const k1 = -P               * _hydroF(z,         dT);
    const k2 = -(P + 0.5*h*k1) * _hydroF(z + 0.5*h, dT);
    const k3 = -(P + 0.5*h*k2) * _hydroF(z + 0.5*h, dT);
    const k4 = -(P + h*k3)     * _hydroF(z + h,      dT);
    return P + (h / 6) * (k1 + 2*k2 + 2*k3 + k4);
}

const Z_STEP = 100;
const Z_MIN  = -6000;
const Z_86K  =  86000;
const _LN    = (Z_86K - Z_MIN) / Z_STEP + 1; // 921
const _LI0   = (0 - Z_MIN) / Z_STEP;          // index of z=0

function _buildLowerTable(dT: number): Float64Array {
    const t = new Float64Array(_LN);
    t[_LI0] = 101325;
    for (let i = _LI0; i < _LN - 1; i++)
        t[i + 1] = _rk4P(Z_MIN + i * Z_STEP, t[i],  Z_STEP, dT);
    for (let i = _LI0; i > 0; i--)
        t[i - 1] = _rk4P(Z_MIN + i * Z_STEP, t[i], -Z_STEP, dT);
    return t;
}

const _LP = _buildLowerTable(0); // standard lower pressure table

function _lowerLookup(t: Float64Array, z: number): number {
    const idx = (z - Z_MIN) / Z_STEP;
    const i   = Math.floor(idx);
    if (idx <= 0)      return t[0];
    if (i >= _LN - 1) return t[_LN - 1];
    return t[i] * Math.pow(t[i + 1] / t[i], idx - i);
}

// Non-standard lower-atmosphere cache (one slot keyed by dT)
let _nsDT = NaN, _nsTable: Float64Array | null = null;

function _lowerNonstd(z: number, P0: number, T0: number): number {
    const dT = T0 - 288.15;
    if (Math.abs(dT) < 1e-9) return _lowerLookup(_LP, z) * (P0 / 101325);
    if (_nsDT !== dT || _nsTable === null) { _nsDT = dT; _nsTable = _buildLowerTable(dT); }
    return _lowerLookup(_nsTable, z) * (P0 / 101325);
}

// ── Upper-atmosphere ODE (86 km – 1000 km) ───────────────────────────────────
// USSA 1976 Part 3 / Pietrobon (1999). Species: 0=N2 1=O 2=O2 3=Ar 4=He 5=H.

const k_B  = 1.380649e-23; // SI 2019 exact value
const Na   = 6.02214076e26;
const Rs   = 1000 * R;
const g0   = 9.80665;

const _uM:   [number,number,number,number,number,number] = [28.0134,15.9994,31.9988,39.948,4.0026,1.00797];
const _uMkg: [number,number,number,number,number,number] = [0.0280134,0.0159994,0.0319988,0.039948,0.0040026,0.00100797];

const _Z7   = 86000;
const _T7   = 186.8673;
const _n7: [number,number,number,number,number,number] =
    [1.129794e20,8.6e16,3.030898e19,1.3514e18,7.58173e14,8.0e10];
const _alpha: [number,number,number,number,number,number] = [0,0,0,0,-0.40,-0.25];
const _Pn7 = _n7.reduce((s, ni) => s + ni, 0) * k_B * _T7;

const _Zh = 150000, _Zm = 100000, _Zk = 115000;
const _Z11 = 500000, _ZMAX = 1000000;
const _taui = 0.40463343, _Hi = 9.8776951e10, _T11 = 999.2356;

function _uGrav(Z: number): number { return g0 * (EARTH_RADIUS_M / (EARTH_RADIUS_M + Z)) ** 2; }

function _uEddy(Z: number): number {
    const K7 = 120, Ze = 95000, A = 4e8;
    if (Z < Ze) return K7;
    if (Z < _Zk - 0.1) return K7 * Math.exp(1 - A / (A - (Z - Ze) ** 2));
    return 0;
}

function _uFlux(Z: number, s: number): number {
    const Q = [-5.809644e-13,1.366212e-13,9.434079e-14,-2.457369e-13];
    const U = [56903.11,86000,86000,86000];
    const W = [2.706240e-14,8.333333e-14,8.333333e-14,6.666667e-13];
    const si = s - 1, a = Z - U[si], A2 = a * a;
    let x = Q[si] * A2 * Math.exp(-W[si] * A2 * a);
    if (s === 1 && Z < 97000) {
        const a2 = 97000 - Z, B = a2 * a2;
        x += -3.416248e-12 * B * Math.exp(-5.008765e-13 * B * a2);
    }
    return x;
}

const _dA = [6.986e20,4.863e20,4.487e20,1.700e21,3.305e21];
const _dB = [0.750,0.750,0.870,0.691,0.500];

function _fcalc(Z: number, g: number, K: number, T: number, dTdz: number, s: number, n: number[]): number {
    if (s === 0) return (Z < _Zm ? M0 * 1000 : _uM[0]) * g / (Rs * T);
    if (s >= 1 && s <= 4) {
        let f: number;
        if (Z < _Zk - 0.1) {
            let sumn = 0; for (let i = 0; i <= (s <= 2 ? 0 : s <= 4 ? 2 : 4); i++) sumn += n[i];
            const D = _dA[s-1] * (T/273.15)**_dB[s-1] / sumn;
            let Mmean = 0, sn = 0;
            for (let i = 0; i <= 4; i++) { sn += n[i]; Mmean += n[i]*_uM[i]; }
            Mmean = Z < _Zm ? M0*1000 : Mmean/sn;
            f = D * (_uM[s] + Mmean*K/D + _alpha[s]*Rs*dTdz/g) / (D + K);
        } else {
            f = _uM[s] + _alpha[s] * Rs * dTdz / g;
        }
        return g * f / (Rs * T) + _uFlux(Z, s);
    }
    return 0;
}

// Precomputed upper atmosphere tables (standard conditions, 86 km to 1000 km).
const _UN  = (_ZMAX - _Z7) / Z_STEP + 1; // 9141
const _uStdP   = new Float64Array(_UN);
const _uStdRho = new Float64Array(_UN);
const _uStdXi: [Float64Array,Float64Array,Float64Array,Float64Array,Float64Array,Float64Array] = [
    new Float64Array(_UN), new Float64Array(_UN), new Float64Array(_UN),
    new Float64Array(_UN), new Float64Array(_UN), new Float64Array(_UN),
];

function _runUpperODE(
    P_bc: number, T_bc: number,
    pOut: Float64Array, rhoOut: Float64Array,
    xiOut: Float64Array[] | null,
): void {
    const nscale = (P_bc / T_bc) / (_Pn7 / _T7);
    const n = _n7.map(v => v * nscale) as [number,number,number,number,number,number];

    function store(i: number, ni: number[], T: number) {
        const P = ni.reduce((s, x) => s + x * k_B * T, 0);
        pOut[i]   = P;
        rhoOut[i] = ni.reduce((s, x, j) => s + x * _uM[j], 0) / Na; // n[1/m³] × M[kg/kmol] / Na[1/kmol] = kg/m³
        if (xiOut) {
            const tot = ni.reduce((s, x) => s + x, 0);
            for (let s = 0; s < 6; s++) xiOut[s][i] = tot > 0 ? ni[s] / tot : 0;
        }
    }
    store(0, n, T_bc);

    const fint = [0,0,0,0,0,0];
    const fprev = [0,1,2,3,4].map(s => _fcalc(_Z7, _uGrav(_Z7), _uEddy(_Z7), T_bc, 0, s, n));
    fprev[5] = 0;
    let tauint = 0, tauPrev = 0, hStarted = false;
    let Z = _Z7;

    for (let step = 0; step < _UN - 1; step++) {
        const Znext = Z + Z_STEP;
        const Zm2   = Z + Z_STEP / 2;
        // Use T_bc offset: above 86 km, temperature = _upperT(z) + (T_bc - _T7)
        const dToff = T_bc - _T7;
        const Tmid  = _upperT(Zm2)  + dToff; const dTm = _upperDTdZ(Zm2);
        const gmid  = _uGrav(Zm2),  Kmid  = _uEddy(Zm2);
        const Tnext = _upperT(Znext) + dToff; const dTn = _upperDTdZ(Znext);
        const gnext = _uGrav(Znext), Knext = _uEddy(Znext);

        for (let s = 0; s <= 4; s++) {
            const fm = _fcalc(Zm2,   gmid,  Kmid,  Tmid,  dTm, s, n);
            const fn = _fcalc(Znext, gnext, Knext, Tnext, dTn, s, n);
            fint[s] += Z_STEP * (fprev[s] + 4*fm + fn) / 6;
            fprev[s] = fn;
        }
        for (let s = 0; s <= 4; s++) n[s] = _n7[s] * nscale * (T_bc / Tnext) * Math.exp(-fint[s]);

        if (Znext >= _Zh) {
            if (!hStarted) { tauPrev = (gnext / Tnext) * (_uM[5] / Rs); hStarted = true; }
            else {
                const tauNext = (gnext / Tnext) * (_uM[5] / Rs);
                tauint += Z_STEP * (tauPrev + tauNext) / 2;
                tauPrev = tauNext;
            }
            n[5] = Znext < _Z11
                ? (_n7[5]*nscale + _Hi - fint[5]) *
                  Math.exp((1 + _alpha[5]) * Math.log(_T11 / Tnext) + _taui - tauint)
                : 0;
        }

        store(step + 1, n, Tnext);
        Z = Znext;
    }
}

// Run standard upper atmosphere table at init (P_bc = standard lower pressure at 86 km).
const _P86std = _LP[_LN - 1];
_runUpperODE(_P86std, _T7, _uStdP, _uStdRho, _uStdXi);

// Non-standard upper atmosphere cache (one slot keyed by P0,T0)
let _nsUpKey = '', _nsUpP: Float64Array | null = null, _nsUpRho: Float64Array | null = null;

function _ensureNsUp(P0: number, T0: number): void {
    const key = `${P0}|${T0}`;
    if (_nsUpKey === key) return;
    const P_bc = _lowerNonstd(Z_86K, P0, T0);
    const T_bc = _T7 + (T0 - 288.15);
    const p = new Float64Array(_UN), rho = new Float64Array(_UN);
    _runUpperODE(P_bc, T_bc, p, rho, null); // xi doesn't change with nscale
    _nsUpP = p; _nsUpRho = rho; _nsUpKey = key;
}

function _upperLookup(arr: Float64Array, z: number): number {
    const idx = (z - _Z7) / Z_STEP;
    const i   = Math.floor(idx);
    if (idx <= 0)      return arr[0];
    if (i >= _UN - 1) return arr[_UN - 1];
    return arr[i] * Math.pow(arr[i+1] / arr[i], idx - i);
}

function _upperXi(s: number, z: number): number {
    const idx = (z - _Z7) / Z_STEP, i = Math.floor(idx);
    const a = _uStdXi[s];
    if (i <= 0) return a[0];
    if (i >= _UN - 1) return a[_UN - 1];
    return a[i] + (idx - i) * (a[i+1] - a[i]);
}

// ── Trace gas decay above 86 km ───────────────────────────────────────────────
// CO2, Ne, CH4, Kr, H2, N2O, Xe are absent from the upper-atmosphere ODE.
// Above 86 km they decay with diffusive scale heights H_s = R·T86 / (M_s·g86).
// Their budget is removed from the 6 ODE species proportionally so fractions sum to 1.

const _TRACE_KEYS = ['CO2', 'Ne', 'CH4', 'Kr', 'H2', 'N2O', 'Xe'] as const;
const _TRACE_MKG  = [0.04400950, 0.02017970, 0.01604246, 0.08379800, 0.00201588, 0.04401280, 0.13129300] as const;
const _g86 = 9.80665 * (EARTH_RADIUS_M / (EARTH_RADIUS_M + Z_86K)) ** 2;
const _TRACE_HS = (_TRACE_MKG as readonly number[]).map(M => R * _T7 / (M * _g86));
// Mole fractions of trace gases at exactly 86 km, from the lower-atmosphere model
const _TRACE_XI86 = (() => {
    const f = _lowerXi(H_86KM);
    return (_TRACE_KEYS as readonly string[]).map(k => f[k] ?? 0);
})();

function _traceDecay(dz: number): number[] {
    return _TRACE_XI86.map((xi0, i) => xi0 * Math.exp(-dz / _TRACE_HS[i]));
}

// Blend zone 86–87 km: lower-model composition (frozen at 86 km) transitions smoothly
// to the ODE upper model. Ensures M is continuous at 86 km so density stays monotonic.
const _BLEND_TOP = 87000;
const _LXI_86 = _lowerXi(H_86KM); // lower-model composition at exactly 86 km

// ── Unified public API ────────────────────────────────────────────────────────
// All functions accept geometric altitude z in meters. No 86 km branch visible to callers.

export function getPressure(z: number): number {
    if (z <= Z_86K) return _lowerLookup(_LP, z);
    return _upperLookup(_uStdP, z);
}

// Pressure with non-standard sea-level conditions.
// Above 86 km: re-integrates the upper atmosphere ODE with the correct boundary
// conditions (T_bc = T7 + dT, P_bc from lower atmosphere at 86 km), cached by (P0, T0).
export function getPressureNonstd(z: number, P0: number, T0: number): number {
    if (z <= Z_86K) return _lowerNonstd(z, P0, T0);
    _ensureNsUp(P0, T0);
    return _upperLookup(_nsUpP!, z);
}

export function getDensity(z: number): number {
    if (z <= Z_86K) {
        const P = _lowerLookup(_LP, z);
        return P * _lowerMeff(geometricToGeopotential(z)) / (R * getTemperature(z));
    }
    // Use P·M/RT so getDensity is consistent with getDensityFromPT and getMolarMass
    // (which now includes trace-gas contributions above 86 km).
    return getPressure(z) * getMolarMass(z) / (R * getTemperature(z));
}

// Density for measured (P, T) at altitude z — uses M(z) from atmospheric composition.
export function getDensityFromPT(P_Pa: number, T_K: number, z: number): number {
    return P_Pa * getMolarMass(z) / (R * T_K);
}

export function getMolarMass(z: number): number {
    if (z <= Z_86K) return _lowerMeff(geometricToGeopotential(z));
    const f = getMoleFractions(z);
    const UPPER_KEYS = ['N2','O','O2','Ar','He','H'] as const;
    let M = 0;
    for (let s = 0; s < 6; s++) M += (f[UPPER_KEYS[s]] ?? 0) * _uMkg[s];
    for (let i = 0; i < _TRACE_KEYS.length; i++) M += (f[_TRACE_KEYS[i]] ?? 0) * _TRACE_MKG[i];
    return M;
}

// Mole fractions at z. Always sums to exactly 1.
// Above 87 km: 6 ODE species + trace gases decaying with diffusive scale heights.
// 86–87 km: smoothly blended from lower model (at 86 km) to ODE+trace model,
//   ensuring molar mass is continuous at 86 km so density remains monotonic.
export function getMoleFractions(z: number): Record<string, number> {
    if (z <= Z_86K) return _lowerXi(geometricToGeopotential(z));
    const UPPER_KEYS = ['N2','O','O2','Ar','He','H'] as const;

    // Build pure upper-model fractions
    const f_upper: Record<string, number> = {};
    for (let s = 0; s < 6; s++) f_upper[UPPER_KEYS[s]] = _upperXi(s, z);
    const trace = _traceDecay(z - Z_86K);
    let traceTot = 0;
    for (let i = 0; i < _TRACE_KEYS.length; i++) {
        f_upper[_TRACE_KEYS[i]] = trace[i];
        traceTot += trace[i];
    }
    const upperScale = 1 - traceTot;
    for (const k of UPPER_KEYS) f_upper[k] *= upperScale;

    if (z >= _BLEND_TOP) return f_upper;

    // Blend zone: mix lower-model composition (frozen at 86 km) with upper model
    const t = (z - Z_86K) / (_BLEND_TOP - Z_86K);
    const alpha = t * t * (3 - 2 * t); // smoothstep
    const f: Record<string, number> = {};
    let tot = 0;
    for (const k of Object.keys(f_upper)) {
        const v = (1 - alpha) * (_LXI_86[k] ?? 0) + alpha * f_upper[k];
        f[k] = v; tot += v;
    }
    const inv = 1 / tot;
    for (const k of Object.keys(f)) f[k] *= inv;
    return f;
}

// ── Derived functions ─────────────────────────────────────────────────────────

export function getSeaLevelTemperature(z: number, T_K: number): number {
    return 288.15 + T_K - getTemperature(z);
}

export function getSeaLevelPressure(z: number, P_Pa: number, T0_K: number): number {
    return P_Pa * 101325 / getPressureNonstd(z, 101325, T0_K);
}

export function getAltimeterSetting(z: number, P_Pa: number): number {
    return P_Pa * 101325 / getPressure(z);
}

// Pressure altitude: geometric z where standard pressure = P_Pa (−6 km to 1000 km).
export function getPressureAltitude(P_Pa: number): number {
    let lo = Z_MIN, hi = _ZMAX;
    for (let i = 0; i < 60; i++) {
        const mid = (lo + hi) / 2;
        getPressure(mid) > P_Pa ? lo = mid : hi = mid;
    }
    return hi;
}

// Density altitude: geometric z where standard density = rho (−6 km to 1000 km).
export function getDensityAltitude(rho: number): number {
    let lo = Z_MIN, hi = _ZMAX;
    for (let i = 0; i < 60; i++) {
        const mid = (lo + hi) / 2;
        getDensity(mid) > rho ? lo = mid : hi = mid;
    }
    return hi;
}

// ── Physical properties ───────────────────────────────────────────────────────

export function getSpeedOfSound(P_Pa: number, rho: number): number | null {
    if (P_Pa < 3e-4) return null; // continuum assumption invalid below ~0.0003 Pa (~160 km)
    return Math.sqrt(1.4 * P_Pa / rho);
}

const SIGMA_MFP = 3.65e-10;
export function getMeanFreePath(P_Pa: number, T_K: number): number {
    const n = P_Pa / (k_B * T_K);
    return 1 / (Math.SQRT2 * Math.PI * SIGMA_MFP * SIGMA_MFP * n);
}

export function waterVaporPressure(T: number): number {
    const Tc = 647.096, Pc = 22064000, theta = 1 - T / Tc;
    return Pc * Math.exp((Tc / T) * (
        -7.85951783*theta + 1.84408259*theta**1.5 - 11.7866497*theta**3 +
        22.6807411*theta**3.5 - 15.9618719*theta**4 + 1.80122502*theta**7.5
    ));
}

export function getBoilingPoint(P_Pa: number): number {
    let lo = 273, hi = 647.096;
    for (let i = 0; i < 30; i++) {
        const mid = (lo + hi) / 2;
        waterVaporPressure(mid) < P_Pa ? lo = mid : hi = mid;
    }
    return hi;
}

// ── Exports for pages ─────────────────────────────────────────────────────────

export const LOWER_MOLE_FRACTIONS: Record<string, number> = {
    ...Object.fromEntries(LOWER_COMPOSITION.map(c => [c.key, c.xi])),
    O: 0, H: 0,
};
