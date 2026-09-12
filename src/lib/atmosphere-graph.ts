import {
    R, buildAtmosphereProfile, type AtmosphereProfile,
} from './atmosphere';

// Build the sorted altitude array for a graph range: 1001 evenly spaced points.
export function buildGraphAltitudes(zMin: number, zMax: number): number[] {
    const N = 1000;
    const pts: number[] = [];
    for (let i = 0; i <= N; i++) pts.push(zMin + (zMax - zMin) * i / N);
    return pts;
}

// Compute density at each profile point from P, T, M.
export function profileDensity(profile: AtmosphereProfile): Float64Array {
    const n   = profile.P.length;
    const rho = new Float64Array(n);
    for (let i = 0; i < n; i++) rho[i] = profile.P[i] * profile.M[i] / (R * profile.T[i]);
    return rho;
}
