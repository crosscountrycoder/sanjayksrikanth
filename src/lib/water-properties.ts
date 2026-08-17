// ── Water vapor and boiling point ────────────────────────────────────────────
// IAPWS-95 saturation pressure formula, accurate to 0.01°C.
// Below the triple point (611.657 Pa) liquid water cannot exist.

export function waterVaporPressure(T: number): number {
    if (T > 647.096) return NaN; // no vapor pressure above critical temperature
    if (T >= 273.1600117513473) {
        /* Note: 273.1600117513473 is used rather than 273.16 due to rounding, which causes a very small discontinuity
        at the triple point. The IAPWS-95 ice and water formulas intersect at 273.1600117513473 K. */
        const Tc = 647.096, Pc = 22064000, tau = 1 - T / Tc;
        return Pc * Math.exp((Tc / T) * (
            -7.85951783*tau + 1.84408259*tau**1.5 - 11.7866497*tau**3 +
            22.6807411*tau**3.5 - 15.9618719*tau**4 + 1.80122502*tau**7.5
        ));
    }
    const Tt = 273.16, Pt = 611.657, theta = T / Tt;
    return Pt * Math.exp((
        -21.2144006*theta**0.00333333333 +
        27.3203819*theta**1.20666667 -
        6.10598130*theta**1.70333333
    ) / theta);
}

/** Returns the boiling or sublimation point of water in kelvins, at pressure P_Pa pascals.
 * If the value returned is over 273.16 K, this is the boiling point; otherwise, it's the sublimation point.
*/
export function getBoilingPoint(P_Pa: number): number {
    if (P_Pa > 22064000) return NaN; // supercritical fluid, no boiling point
    let lo = 0, hi = 647.096;
    while (hi - lo > 1e-9) {
        const mid = (lo + hi) / 2;
        waterVaporPressure(mid) < P_Pa ? lo = mid : hi = mid;
    }
    return lo;
}

/** Returns relative humidity given temperature and dew point, both in kelvins. 
 * Value is dimensionless and in the range [0, 1] assuming T_d <= T. */
export function getRelativeHumidity(T: number, T_d: number): number {
    return waterVaporPressure(T_d) / waterVaporPressure(T);
}

/** Returns dew point given temperature in kelvins and relative humidity in range [0, 1]. 
 * Value is in kelvins. */
export function getDewPoint(T: number, RH: number): number {
    const e = waterVaporPressure(T) * RH;
    return getBoilingPoint(e);
}

/** Returns temperature given dew point in kelvins and relative humidity in range [0, 1]. 
 * Value is in kelvins. */
export function getTempFromRHAndDewPoint(T_d: number, RH: number): number {
    const e = waterVaporPressure(T_d) / RH;
    return getBoilingPoint(e);
}