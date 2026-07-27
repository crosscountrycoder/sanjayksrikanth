import {
    R, buildAtmosphereProfile, type AtmosphereProfile,
    getPressureAltitude, getDensityAltitude,
} from './atmosphere';

// Build the sorted altitude array for a graph range: exact endpoints plus
// every 100 m grid point in between.
export function buildGraphAltitudes(zMin: number, zMax: number): number[] {
    const pts: number[] = [];
    const gridStart = Math.ceil((zMin + 1e-9) / 100) * 100;
    if (zMin < gridStart - 1e-9) pts.push(zMin);
    for (let z = gridStart; z < zMax - 1e-9; z += 100) pts.push(z);
    pts.push(zMax);
    return pts;
}

// Compute density at each profile point from P, T, M.
export function profileDensity(profile: AtmosphereProfile): Float64Array {
    const n   = profile.P.length;
    const rho = new Float64Array(n);
    for (let i = 0; i < n; i++) rho[i] = profile.P[i] * profile.M[i] / (R * profile.T[i]);
    return rho;
}

// Compute pressure altitudes for each profile point (fast log-lin mode).
export function profilePressureAltitude(profile: AtmosphereProfile): Float64Array {
    const n  = profile.P.length;
    const pa = new Float64Array(n);
    for (let i = 0; i < n; i++) pa[i] = getPressureAltitude(profile.P[i], true);
    return pa;
}

// Compute density altitudes for each profile point (fast log-lin mode).
export function profileDensityAltitude(profile: AtmosphereProfile, rho: Float64Array): Float64Array {
    const n  = rho.length;
    const da = new Float64Array(n);
    for (let i = 0; i < n; i++) da[i] = getDensityAltitude(rho[i], true);
    return da;
}

export { buildAtmosphereProfile };
