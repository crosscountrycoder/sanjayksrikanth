import { G0, LAPSE_RATE, R } from './atmosphere.ts';

// ── Station pressure ─────────────────────────────────────────────────────────
// Barometric formula between two stations, assuming temperature varies linearly
// with geopotential height between them (constant lapse rate). Integrating the
// hydrostatic equation dP/dH = −P · M · g0 / (R · T) gives
//     P_B = P_A · (T_B / T_A) ^ (−(M · g0 / R) · ΔH / (T_B − T_A))
// For an isothermal layer (T_A = T_B) this becomes
//     P_B = P_A · exp(−(M · g0 / R) · ΔH / T)

const M_AIR = 28.9659;
const K_HYDRO = G0 * M_AIR / R; // g0 · M / R (K/m): the "hydrostatic constant" in the exponent.

// Bounds for the standard-lapse-rate station temperature (K): −56.5 °C to 73.5 °C.
export const STATION_TEMP_MIN_K = 216.65;
export const STATION_TEMP_MAX_K = 346.65;

// Geopotential altitude range for station inputs (m): the model's −5 km floor up to the
// standard tropopause, above which a constant lapse rate no longer holds.
export const STATION_ALT_MIN_M = -5000;
export const STATION_ALT_MAX_M = 11000;

// Temperature at station B assuming the standard lapse rate from station A,
// clamped to [STATION_TEMP_MIN_K, STATION_TEMP_MAX_K].
export function getStandardStationTemperature(H_A: number, T_A_K: number, H_B: number): number {
    const T_B_K = T_A_K - LAPSE_RATE * (H_B - H_A);
    return Math.min(Math.max(T_B_K, STATION_TEMP_MIN_K), STATION_TEMP_MAX_K);
}

// Pressure at station B given pressure and temperature at station A.
// Heights are geopotential (m), temperatures in K. Pressure units carry through
// (P_B is returned in the same units as P_A). If T_B is omitted, it is derived
// from T_A using the standard lapse rate.
export function getStationPressure(
    H_A: number, P_A: number, T_A: number,
    H_B: number, T_B: number = getStandardStationTemperature(H_A, T_A, H_B),
): number {
    if (P_A <= 0 || T_A <= 0 || T_B <= 0) return NaN;
    const dH = H_B - H_A;
    const dT = T_B - T_A;
    if (Math.abs(dT) < 1e-4) {
        return P_A * Math.exp(-K_HYDRO * dH / ((T_A + T_B) / 2));
    }
    return P_A * (T_B / T_A) ** (-K_HYDRO * dH / dT);
}

// ── Altimeter setting ────────────────────────────────────────────────────────
// Ratio of station pressure to altimeter setting at geopotential altitude H (m): the
// standard-atmosphere troposphere pressure ratio, (1 − H / (T0/L)) ^ (g0·M / (R·L)),
// with T0 = 288.15 K, L = 0.0065 K/m, and the M and R used above.
export function getAltimeterPressureRatio(H: number): number {
    return (1 - H / 44330.769) ** 5.2560581;
}
