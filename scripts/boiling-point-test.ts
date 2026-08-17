import {getBoilingPoint} from '../src/lib/water-properties.ts';
console.log(getBoilingPoint(0)); // 0
console.log(getBoilingPoint(1e-300)); // close to 0
console.log(getBoilingPoint(611.657)); // 273.16
console.log(getBoilingPoint(101325)); // 373.124
console.log(getBoilingPoint(22064000)); // 647.096