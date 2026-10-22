/**
 * DroneRoutePlanner — Boustrophedon Coverage Path Planning Algorithm
 * ===================================================================
 * Оптимальный маршрут дрона для сплошного покрытия поля (паттерн «вол»).
 *
 * Алгоритм:
 *  1. Вычисляем BOUNDING BOX полигона поля
 *  2. Находим оптимальный угол разворота (перебор 0-90°, минимизируем путь)
 *  3. Поворачиваем полигон в локальную систему координат
 *  4. Нарезаем параллельные полосы с учётом ширины захвата и перекрытия
 *  5. Обрезаем полосы по границе полигона (ray-casting intersection)
 *  6. Соединяем полосы в серпантин (чётные — слева направо, нечётные — справа налево)
 *  7. Поворачиваем маршрут обратно в WGS-84
 *
 * Использование (браузер):
 *   const planner = new DroneRoutePlanner({ swathWidth: 15, overlap: 20 });
 *   const result  = planner.plan(polygonLatLng);
 *
 * Использование (Node.js):
 *   const { DroneRoutePlanner } = require('./drone_route.js');
 */

;(function (global) {
  'use strict';

  const R_EARTH = 6_371_000; // Earth radius, m

  /** Converts degrees → radians */
  const rad = (d) => (d * Math.PI) / 180;

  /** Haversine distance between two [lat, lng] points, metres */
  function haversine([la1, lo1], [la2, lo2]) {
    const dLat = rad(la2 - la1);
    const dLon = rad(lo2 - lo1);
    const a =
      Math.sin(dLat / 2) ** 2 +
      Math.cos(rad(la1)) * Math.cos(rad(la2)) * Math.sin(dLon / 2) ** 2;
    return 2 * R_EARTH * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
  }

  /** Approximate metres-per-degree at a given latitude */
  function mPerDeg(refLat) {
    return {
      lat: 111_132 - 559.8 * Math.cos(rad(2 * refLat)),
      lng: 111_412 * Math.cos(rad(refLat)) - 93.5 * Math.cos(rad(3 * refLat)),
    };
  }

  /** Project [lat,lng] onto a flat plane centered at origin with angle θ (rad) */
  function project([lat, lng], origin, mpd, theta) {
    const x = (lng - origin[1]) * mpd.lng;
    const y = (lat - origin[0]) * mpd.lat;
    return [
      x * Math.cos(theta) - y * Math.sin(theta),
      x * Math.sin(theta) + y * Math.cos(theta),
    ];
  }

  /** Unproject flat-plane point back to [lat,lng] */
  function unproject([px, py], origin, mpd, theta) {
    const x = px * Math.cos(-theta) - py * Math.sin(-theta);
    const y = px * Math.sin(-theta) + py * Math.cos(-theta);
    return [origin[0] + y / mpd.lat, origin[1] + x / mpd.lng];
  }

  /**
   * Find X-coordinate where horizontal line at Y=y crosses the edge p1→p2.
   * Returns null if the edge doesn't cross.
   */
  function edgeIntersectX(y, [x1, y1], [x2, y2]) {
    if ((y1 <= y && y2 > y) || (y2 <= y && y1 > y)) {
      return x1 + ((y - y1) * (x2 - x1)) / (y2 - y1);
    }
    return null;
  }

  /**
   * For a projected polygon, return all x-intersection values at row-latitude y.
   * Sorted ascending so pairs [xs[0],xs[1]], [xs[2],xs[3]], … are strip segments.
   */
  function slicePolygon(y, projPoly) {
    const xs = [];
    const n = projPoly.length;
    for (let i = 0; i < n; i++) {
      const x = edgeIntersectX(y, projPoly[i], projPoly[(i + 1) % n]);
      if (x !== null) xs.push(x);
    }
    xs.sort((a, b) => a - b);
    return xs;
  }

  /**
   * Generate boustrophedon waypoints for one sweep angle θ.
   * Returns { waypoints, totalDist } in projected (flat) space.
   */
  function boustrophedon(projPoly, swathM, overlapFrac) {
    let minY = Infinity, maxY = -Infinity;
    for (const [, y] of projPoly) {
      minY = Math.min(minY, y);
      maxY = Math.max(maxY, y);
    }

    const step = swathM * (1 - overlapFrac);
    const waypoints = [];
    let stripIdx = 0;

    for (let y = minY + step / 2; y < maxY; y += step) {
      const xs = slicePolygon(y, projPoly);
      for (let k = 0; k + 1 < xs.length; k += 2) {
        const goRight = stripIdx % 2 === 0;
        waypoints.push(
          goRight ? [xs[k], y] : [xs[k + 1], y],
          goRight ? [xs[k + 1], y] : [xs[k], y]
        );
        stripIdx++;
      }
    }

    let totalDist = 0;
    for (let i = 1; i < waypoints.length; i++) {
      const dx = waypoints[i][0] - waypoints[i - 1][0];
      const dy = waypoints[i][1] - waypoints[i - 1][1];
      totalDist += Math.sqrt(dx * dx + dy * dy);
    }

    return { waypoints, totalDist, stripCount: stripIdx };
  }

  /**
   * Find the sweep angle that minimises total flight distance.
   * Samples angles 0°–179° in `resolution`-degree increments.
   */
  function optimalAngle(projPoly, swathM, overlapFrac, resolution = 5) {
    let best = { dist: Infinity, theta: 0 };
    for (let deg = 0; deg < 180; deg += resolution) {
      const theta = rad(deg);
      // Rotate polygon by theta and re-project (simple 2D rotation in flat space)
      const rotated = projPoly.map(([x, y]) => [
        x * Math.cos(theta) - y * Math.sin(theta),
        x * Math.sin(theta) + y * Math.cos(theta),
      ]);
      const { totalDist } = boustrophedon(rotated, swathM, overlapFrac);
      if (totalDist < best.dist) best = { dist: totalDist, theta };
    }
    return best.theta;
  }

  // ─── Public class ───────────────────────────────────────────────────────────

  class DroneRoutePlanner {
    /**
     * @param {object} options
     * @param {number} options.swathWidth   Camera swath width in metres (default 15)
     * @param {number} options.overlap      Strip overlap 0–1 fraction (default 0.20)
     * @param {number} options.speed        Cruise speed m/s (default 8)
     * @param {number} options.altitude     Flight altitude m (default 30)
     * @param {boolean} options.optimise    Try to find optimal sweep angle (default true)
     */
    constructor({
      swathWidth = 15,
      overlap = 0.20,
      speed = 8,
      altitude = 30,
      optimise = true,
    } = {}) {
      this.swathWidth = swathWidth;
      this.overlap    = overlap;
      this.speed      = speed;
      this.altitude   = altitude;
      this.optimise   = optimise;
    }

    /**
     * Plan a coverage route for the given field polygon.
     *
     * @param {Array<[number,number]>} polygon  Array of [lat, lng] vertex pairs (WGS-84)
     * @returns {{
     *   waypoints: Array<[number,number]>,  WGS-84 route waypoints
     *   strips: number,                     Number of coverage strips
     *   totalDistanceM: number,             Total flight distance in metres
     *   flightTimeMin: number,              Estimated flight time in minutes
     *   fieldAreaM2: number,                Approximate field area in m²
     *   coveragePercent: number,            Estimated coverage %
     * }}
     */
    plan(polygon) {
      if (polygon.length < 3) throw new Error('Polygon must have at least 3 vertices');

      // Reference point (centroid) for local flat projection
      const refLat = polygon.reduce((s, [la]) => s + la, 0) / polygon.length;
      const refLng = polygon.reduce((s, [, lo]) => s + lo, 0) / polygon.length;
      const origin = [refLat, refLng];
      const mpd    = mPerDeg(refLat);

      // Project polygon to flat metric space
      const projPoly = polygon.map((p) => project(p, origin, mpd, 0));

      // Optionally find the sweep angle that minimises distance
      const theta = this.optimise
        ? optimalAngle(projPoly, this.swathWidth, this.overlap)
        : 0;

      // Rotate polygon by optimal angle in flat space, then plan
      const rotated = projPoly.map(([x, y]) => [
        x * Math.cos(theta) - y * Math.sin(theta),
        x * Math.sin(theta) + y * Math.cos(theta),
      ]);

      const { waypoints: flatWpts, totalDist, stripCount } =
        boustrophedon(rotated, this.swathWidth, this.overlap);

      // Rotate waypoints back and unproject to WGS-84
      const waypoints = flatWpts.map(([px, py]) => {
        const x =  px * Math.cos(-theta) - py * Math.sin(-theta);
        const y =  px * Math.sin(-theta) + py * Math.cos(-theta);
        return unproject([x, y], origin, mpd, 0);
      });

      // Field area approximation (shoelace in projected space)
      let area = 0;
      for (let i = 0; i < projPoly.length; i++) {
        const [x1, y1] = projPoly[i];
        const [x2, y2] = projPoly[(i + 1) % projPoly.length];
        area += x1 * y2 - x2 * y1;
      }
      const fieldAreaM2 = Math.abs(area) / 2;

      // Coverage: swath area / field area
      const swathArea = stripCount * this.swathWidth * (totalDist / Math.max(stripCount, 1));
      const coveragePercent = Math.min(100, Math.round((swathArea / fieldAreaM2) * 100));

      // Real distance via haversine
      let realDistM = 0;
      for (let i = 1; i < waypoints.length; i++) {
        realDistM += haversine(waypoints[i - 1], waypoints[i]);
      }

      return {
        waypoints,
        strips: stripCount,
        totalDistanceM: Math.round(realDistM),
        flightTimeMin: Math.round((realDistM / this.speed / 60) * 10) / 10,
        fieldAreaM2: Math.round(fieldAreaM2),
        coveragePercent,
        sweepAngleDeg: Math.round((theta * 180) / Math.PI),
      };
    }
  }

  // Export
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = { DroneRoutePlanner, haversine, mPerDeg };
  } else {
    global.DroneRoutePlanner = DroneRoutePlanner;
  }
})(typeof window !== 'undefined' ? window : global);
