# sanjayksrikanth.com

Personal website for Sanjay Srikanth, built with [Astro](https://astro.build) and deployed on Cloudflare Workers. Includes a biography, links, and 
scientific calculators.

## Online tools

### Atmosphere (`/online-tools/atmosphere/`)

A USSA 1976 atmosphere calculator covering −5 km to 1000 km geometric altitude with non-standard temperature and pressure support. Outputs 
pressure, altimeter setting, temperature, density, pressure altitude, density altitude, gravitational acceleration, geopotential altitude, 
molar mass, speed of sound, dynamic viscosity, mean free path, scale height, number density, mole fraction of each gas, and boiling point of water.

The model integrates the hydrostatic equation using RK4 at 50 m steps with variable molar mass and gravity from a pre-computed composition table 
(`src/lib/mole-fractions.ts`), covering all altitudes in a single unified formulation. Four views are available: point calculator, table, graph, and 
technical notes.

## Project structure

```
src/
  lib/
    atmosphere.ts          — core atmosphere model (pressure, temperature, composition, derived quantities, graph profile builder)
    convert.ts             — unit conversion helpers
    water-properties.ts    — boiling point, vapor pressure, humidity/dew point
    mole-fractions.ts      — pre-computed composition table (−5 km to 1000 km, 1 km steps)
  pages/
    index.astro            — home page / biography
    online-tools/
      index.astro          — online tools index
      atmosphere/
        index.astro        — point calculator
        table.astro        — altitude table
        graph.astro        — property-vs-altitude graph
        notes.astro        — technical notes
  layouts/
    Layout.astro
  components/
scripts/
  atm-test.ts             — full-model test script (uses atmosphere.ts)
  atm-test-simple.ts      — USSA 1976 analytical reference script (independent of atmosphere.ts internals)
wrangler.jsonc             — Cloudflare Workers configuration
```

## Development

Requires Node.js ≥ 24.

```sh
npm install       # install dependencies
npm run dev       # start dev server at http://localhost:4321
npm run build     # build to ./dist/
npm run preview   # preview the build locally
```

## Deployment

The site is deployed as a static site served by a Cloudflare Worker to `sanjayksrikanth.com`.

```sh
npm run deploy
```

To test the Worker locally before deploying:

```sh
npm run build
npx wrangler dev
```

## Test scripts

See [`scripts/README.md`](scripts/README.md) for full documentation.

```sh
npx tsx scripts/atm-test.ts [z_m] [T0_K] [P0_Pa]
npx tsx scripts/atm-test-simple.ts [z_m] [T0_K] [P0_Pa]
npx tsx scripts/atm-test.ts verify-mole-fractions
```
