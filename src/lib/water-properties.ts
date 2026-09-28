// ── Water vapor and boiling point ────────────────────────────────────────────
// IAPWS saturation pressure formula, accurate to within 72 ppm.
// Formulas from https://iapws.org/public/documents/6dGkr/Supp-sat.pdf (water) and https://iapws.org/public/documents/MdUFK/MeltSub2011.pdf (ice)
// Below the triple point (611.657 Pa) liquid water cannot exist.

export function waterVaporPressure(T: number): number {
    if (T < 0) return NaN; // negative temperature, no vapor pressure
    if (T === 0) return 0; // limit of the ice formula at absolute zero (it evaluates to 0/0 there)
    if (T > 647.096) return Infinity; // no vapor pressure above critical temperature
    if (T >= 273.1600117513473) {
        /* The IAPWS ice and water formulas intersect at precisely 273.1600117513473 K due to rounding. 
        The actual triple point is still 273.16 K.*/
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
 * Formula is accurate to within 0.0024 K, and to within 0.0011 K for pressures below 10 megapascals.
*/
export function getBoilingPoint(P_Pa: number): number {
    if (P_Pa < 0) return NaN; // negative pressure, no boiling point
    if (P_Pa > 22064000) return Infinity; // supercritical fluid, no boiling point
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
export function getTempFromDewPointAndRH(T_d: number, RH: number): number {
    const e = waterVaporPressure(T_d) / RH;
    return getBoilingPoint(e);
}

// ── Absolute humidity ────────────────────────────────────────────────────────
// Water vapor treated as an ideal gas: ρ_w = P_w · M_w / (R · T).

const M_WATER = 18.015268;        // Molar mass of VSMOW water, kg/kmol
const R_UNIV  = 8314.46261815324; // J/(kmol·K)

/** Returns absolute humidity (mass of water vapor per unit volume) in kg/m³,
 * given the water vapor partial pressure P_w in pascals and temperature T in kelvins. */
export function getAbsoluteHumidity(P_w: number, T: number): number {
    return P_w * M_WATER / (R_UNIV * T);
}

/** Returns the water vapor partial pressure in pascals, given absolute humidity rho_w in kg/m³
 * and temperature T in kelvins. Inverse of getAbsoluteHumidity. */
export function getVaporPressureFromAbsoluteHumidity(rho_w: number, T: number): number {
    return rho_w * R_UNIV * T / M_WATER;
}
