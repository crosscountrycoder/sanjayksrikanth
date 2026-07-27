# scripts/

Test and diagnostic scripts for the atmosphere model. Run with `npx tsx <script>` from the repository root.

---

## atm-test.ts

Full atmosphere model test using `src/lib/atmosphere.ts`. Supports a special subcommand and a point-calculator mode.

### Subcommand: verify-mole-fractions

```sh
npx tsx scripts/atm-test.ts verify-mole-fractions
```

Checks that every row in `src/lib/mole-fractions.ts` has species mole fractions summing to 1 within 1 × 10<sup>-13</sup>. Exits 0 on success, 1 if any row fails.

### Point calculator

```sh
npx tsx scripts/atm-test.ts [z_m] [T0_K] [P0_Pa]
```

| Argument | Description | Default |
| :------- | :---------- | :------ |
| `z_m`    | Geometric altitude (m) | `0` |
| `T0_K`   | Sea-level temperature (K) | `288.15` |
| `P0_Pa`  | Sea-level pressure (Pa) | `101325` |

**Outputs** (6 significant figures, temperatures also in °C):

- Air pressure (Pa)
- Altimeter setting (Pa)
- Air temperature (K / °C)
- Air density (kg/m³)
- Pressure altitude (m)
- Density altitude (m)
- Geopotential altitude (m)
- Speed of sound (m/s) — `N/A` above ~150 km where pressure < 5 × 10<sup>-4</sup> Pa
- Dynamic viscosity (Pa·s)
- Mean free path (m)
- Boiling point (K / °C) — `N/A` below the triple point (611.657 Pa)
- Mole fractions of all species present above 10<sup>-20</sup>

**Model notes:** Pressure is integrated via RK4 at 100 m steps using variable molar mass and gravity. Molar mass and species fractions are read from 
a 1 km-step pre-computed table and linearly interpolated. Valid range: −6 km to 1000 km geometric. Non-standard T0 shifts all temperatures by dT = T0 
− 288.15 K.

---

## atm-test-simple.ts

Independent USSA 1976 reference implementation using the standard analytical formulas. Uses a constant molar mass (0.028965742 kg/mol) and 
piecewise-linear temperature layers in geopotential altitude, with no numerical integration. Intended for cross-checking `atm-test.ts` results 
against the USSA 1976 standard.

**Valid range: −5 km to 86 km geometric altitude.**

```sh
npx tsx scripts/atm-test-simple.ts [z_m] [T0_K] [P0_Pa]
```

| Argument | Description | Default |
| :------- | :---------- | :------ |
| `z_m`    | Geometric altitude (m) | `0` |
| `T0_K`   | Sea-level temperature (K) | `288.15` |
| `P0_Pa`  | Sea-level pressure (Pa) | `101325` |

**Outputs** (6 significant figures, temperatures also in °C):

- Air pressure (Pa)
- Altimeter setting (Pa)
- Air temperature (K / °C)
- Air density (kg/m³)
- Pressure altitude (m)
- Density altitude (m)
- Geopotential altitude (m)
- Speed of sound (m/s)
- Dynamic viscosity (Pa·s)
- Mean free path (m)
- Boiling point (K / °C)

**Model notes:** Geometric altitude is converted to geopotential altitude (H = R_E · z / (R_E + z), R_E = 6 356 766 m) before applying the USSA 1976 
geopotential layer formulas. Pressure and density altitude are found by binary search back in geometric altitude. Non-standard T0 shifts all layer 
temperatures by dT = T0 − 288.15 K.

### USSA 1976 temperature layers (geopotential altitude)

| Layer | H_base (m) | H_top (m) | T_base (K) | Lapse rate (K/m) |
| :---- | ---------: | --------: | ---------: | ---------------: |
| Troposphere       |     0 | 11 000 | 288.15 | −0.0065 |
| Tropopause        | 11 000 | 20 000 | 216.65 |  0.0000 |
| Stratosphere 1    | 20 000 | 32 000 | 216.65 | +0.0010 |
| Stratosphere 2    | 32 000 | 47 000 | 228.65 | +0.0028 |
| Stratopause       | 47 000 | 51 000 | 270.65 |  0.0000 |
| Mesosphere 1      | 51 000 | 71 000 | 270.65 | −0.0028 |
| Mesosphere 2      | 71 000 | 86 000 | 214.65 | −0.0020 |

Pressure formulas per layer (T_b includes the dT offset):

- **Gradient layer** (L ≠ 0): P = P_b · (T_b / (T_b + L · ΔH))^(M · g<sub>0</sub> / (R · L))
- **Isothermal layer** (L = 0): P = P_b · exp(−M · g<sub>0</sub> · ΔH / (R · T_b))

where ΔH = H − H_base, g<sub>0</sub> = 9.80665 m/s<sup>2</sup>, M = 0.028965742 kg/mol, R = 8.31446261815324 J/(mol·K).
