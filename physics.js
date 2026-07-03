/* CRT deflection-yoke simulator — physics core.
 * Works in browser (window.DeflSim) and under node (module.exports) so the
 * same code the page runs can be unit-tested against analytic results.
 *
 * Model: a toroidal deflection coil pair producing Bx (vertical deflection),
 * wound on a ferrite ring core. Each turn is a CLOSED circuit confined to a
 * single angular position theta: an inner leg along the core's bore radius
 * (gun end -> screen end), a radial crossover at z = L out to the core's
 * outer radius, the outer leg back (screen end -> gun end), and a radial
 * crossover at z = 0 closing the loop. Turns are placed at angles from
 * firstAngle to lastAngle, mirrored to -firstAngle..-lastAngle with reversed
 * current so the coil is symmetric above and below the tube axis; the left
 * coil is the further mirror image through the x = 0 plane, with currents
 * chosen so Bx adds on axis.
 *
 * Field: the SCREENED-TOROID limit of the ferrite core. With mu_r >> 1 the
 * ring is a low-reluctance return path that decouples the outer legs and the
 * radial crossovers from the bore — their flux closes through the ferrite and
 * never reaches the beam. The bore therefore sees only (a) the exact
 * Biot-Savart field of the inner legs and (b) their magnetic images in the
 * core's inner surface. The legs lie ON that surface, so each image coincides
 * with its conductor with strength kappa = (mu_r-1)/(mu_r+1): bore field =
 * (1 + kappa) * B_innerLegs, i.e. very nearly doubled for any real ferrite.
 * This is a boundary-value treatment (in its high-permeability limit), not a
 * scalar boost: it yields the nearly uniform transverse bore field a
 * deflection yoke is designed for, and for WG6100 vertical-coil parameters
 * (46 turns/quadrant, 5.5 A, 40 mm, 19.5 kV) it lands within ~10% of the
 * integral B.dl needed for full-scale deflection with no free parameters.
 * Saturation: the bore flux returns through the ring cross-section crowded by
 * roughly rBore/coreThickness, so the core flux density is estimated as
 * crowd * (1+kappa) * |B_innerLegs| and kappa is derated by the standard
 * tanh knee against coreBSat. At real operating points the core runs far
 * below saturation, so this matters only when driven hard. No hysteresis,
 * gap reluctance, or eddy currents modeled.
 *
 * Electron enters the yoke at full anode energy and drifts — pure magnetic
 * force, relativistic Boris rotation (|v| conserved exactly).
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.DeflSim = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const MU0 = 1.25663706212e-6;
  const QE = 1.602176634e-19;
  const ME = 9.1093837015e-31;
  const C = 2.99792458e8;

  // Defaults are the Wells-Gardner 6100 vertical-coil operating point:
  // 46 turns/quadrant, ~5.5 A full-scale, 40 mm yoke, 19.5 kV anode.
  const defaults = {
    current: 5.5,        // A, full-scale vertical drive
    loops: 46,           // turns per quadrant
    minorRadius: 0.030,  // m, core bore radius at gun end (z = 0)
    majorRadius: 0.0545, // m, core bore radius at screen end (z = yokeLength)
    coreThickness: 0.008,// m, radial thickness of the ferrite ring core
    yokeLength: 0.040,   // m
    firstAngle: 25,      // deg from x-axis
    lastAngle: 80,       // deg
    anodeKV: 19.5,       // beam energy entering the yoke, keV equivalent
    screenZ: 0.198,      // m, measured from yoke entrance; 19" 4:3 half-height ~145 mm
    gunZ: -0.020,        // m, integration start
    coreMuR: 1000,       // small-signal relative permeability of the ferrite core (typical MnZn range ~1000-2500)
    coreBSat: 0.40,      // T, saturation flux density (typical MnZn power ferrite ~0.35-0.5 T)
  };

  // ---- geometry -----------------------------------------------------------
  // Segment soup: flat arrays. Right coil only; the x-mirror is applied
  // analytically inside the field functions (exact, halves the work).
  // `flat` holds every segment of every closed turn (drawing, diagnostics);
  // `flatInner` holds only the bore legs — in the screened-toroid limit the
  // outer legs and crossovers are shunted by the core and contribute nothing
  // to the bore field, so the field functions iterate flatInner alone.
  function buildGeometry(p) {
    const segs = []; // {ax..bz, amp} straight segments, right coil, amp = signed current
    const innerSegs = []; // bore legs only (the screened-toroid field sources)
    const loopsOut = []; // per-turn polyline (right coil) for drawing
    const n = Math.max(1, p.loops | 0);
    const a0 = (p.firstAngle * Math.PI) / 180;
    const a1 = (p.lastAngle * Math.PI) / 180;
    const amp = p.current;
    // bore radius (inner surface, against the tube) tapers with the tube's
    // own funnel; outer radius is the bore plus the core's radial thickness.
    const rIn = (z) => p.minorRadius + ((p.majorRadius - p.minorRadius) * z) / p.yokeLength;
    const rOut = (z) => rIn(z) + p.coreThickness;

    // One turn is a closed quadrilateral at fixed angle theta: bore leg
    // (gun->screen), radial crossover over the core's screen-end face, outer
    // leg (screen->gun), radial crossover over the gun-end face closing the
    // loop. sense flips the leg current direction, used to place a second,
    // oppositely-driven turn at -theta so the coil is symmetric top/bottom
    // (matching the physical winding: turns above and below the tube axis,
    // driven so their fields add rather than cancel).
    function addTurn(th, sense) {
      const ct = Math.cos(th), st = Math.sin(th);
      const a = sense * amp;
      const pts = [
        [rIn(0) * ct, rIn(0) * st, 0],
        [rIn(p.yokeLength) * ct, rIn(p.yokeLength) * st, p.yokeLength],
        [rOut(p.yokeLength) * ct, rOut(p.yokeLength) * st, p.yokeLength],
        [rOut(0) * ct, rOut(0) * st, 0],
      ];
      for (let k = 0; k + 1 < pts.length; k++) {
        const seg = {
          ax: pts[k][0], ay: pts[k][1], az: pts[k][2],
          bx: pts[k + 1][0], by: pts[k + 1][1], bz: pts[k + 1][2],
          amp: a,
        };
        segs.push(seg);
        if (k === 0) innerSegs.push(seg); // bore leg, gun -> screen
      }
      // closing segment back to the start
      segs.push({
        ax: pts[3][0], ay: pts[3][1], az: pts[3][2],
        bx: pts[0][0], by: pts[0][1], bz: pts[0][2],
        amp: a,
      });
      pts.push(pts[0]); // closed polyline for drawing
      loopsOut.push(pts);
    }

    for (let i = 0; i < n; i++) {
      const th = n === 1 ? (a0 + a1) / 2 : a0 + (i * (a1 - a0)) / (n - 1);
      addTurn(th, 1);
      addTurn(-th, -1);
    }
    // Pack into typed arrays for the hot loop: ax ay az bx by bz amp
    function pack(list) {
      const f = new Float64Array(list.length * 7);
      for (let i = 0; i < list.length; i++) {
        const s = list[i], o = i * 7;
        f[o] = s.ax; f[o + 1] = s.ay; f[o + 2] = s.az;
        f[o + 3] = s.bx; f[o + 4] = s.by; f[o + 5] = s.bz;
        f[o + 6] = s.amp;
      }
      return f;
    }
    // Core response: image strength kappa, and the flux-crowding ratio used
    // to estimate core flux density from bore field for the saturation knee
    // (bore return flux ~ B*2*rBore*L squeezed into 2 * thickness * L of ring).
    const kappa = (p.coreMuR - 1) / (p.coreMuR + 1);
    const crowd = ((p.minorRadius + p.majorRadius) / 2) / p.coreThickness;
    return {
      flat: pack(segs), nSegs: segs.length,
      flatInner: pack(innerSegs), nInner: innerSegs.length,
      kappa, crowd,
      loops: loopsOut, params: Object.assign({}, p),
    };
  }

  // Closed-form B of a straight segment a->b carrying current I, at point p:
  // B = (mu0 I / 4pi) * (ab.r1hat - ab.r2hat) * (ab x r1) / |ab x r1|^2
  const KB = MU0 / (4 * Math.PI);

  // Screened-toroid core response: bore field = (1 + kappa_eff) * B_innerLegs,
  // where kappa_eff is the image strength derated by the tanh saturation knee
  // on the estimated core flux density (crowd * linear bore field).
  function coreFactor(geo, bMag) {
    const bCore = geo.crowd * (1 + geo.kappa) * bMag;
    const sat = bCore > 1e-15
      ? (geo.params.coreBSat * Math.tanh(bCore / geo.params.coreBSat)) / bCore
      : 1;
    return 1 + geo.kappa * sat;
  }

  // Full 3-component field of the complete yoke (both coils) at (px,py,pz).
  function fieldB(geo, px, py, pz, out) {
    const f = geo.flatInner, n = geo.nInner;
    let bx = 0, by = 0, bz = 0;
    // mirror = right coil evaluated at (-px, py, pz); a pseudovector reflected
    // through x=0 maps (Bx,By,Bz) -> (Bx,-By,-Bz).
    for (let pass = 0; pass < 2; pass++) {
      const qx = pass === 0 ? px : -px;
      const sy = pass === 0 ? 1 : -1;
      for (let i = 0; i < n; i++) {
        const o = i * 7;
        const ax = f[o], ay = f[o + 1], az = f[o + 2];
        const abx = f[o + 3] - ax, aby = f[o + 4] - ay, abz = f[o + 5] - az;
        const r1x = qx - ax, r1y = py - ay, r1z = pz - az;
        const r2x = r1x - abx, r2y = r1y - aby, r2z = r1z - abz;
        const cx = aby * r1z - abz * r1y;
        const cy = abz * r1x - abx * r1z;
        const cz = abx * r1y - aby * r1x;
        const c2 = cx * cx + cy * cy + cz * cz;
        if (c2 < 1e-20) continue;
        const r1m = Math.sqrt(r1x * r1x + r1y * r1y + r1z * r1z);
        const r2m = Math.sqrt(r2x * r2x + r2y * r2y + r2z * r2z);
        if (r1m < 1e-12 || r2m < 1e-12) continue;
        const k = (KB * f[o + 6] *
          ((abx * r1x + aby * r1y + abz * r1z) / r1m -
           (abx * r2x + aby * r2y + abz * r2z) / r2m)) / c2;
        bx += k * cx;
        by += sy * k * cy;
        bz += sy * k * cz;
      }
    }
    // (bx,by,bz) is the inner-leg Biot-Savart field; images in the core's
    // inner surface scale it by (1 + kappa_eff).
    const cf = coreFactor(geo, Math.sqrt(bx * bx + by * by + bz * bz));
    out.x = bx * cf;
    out.y = by * cf;
    out.z = bz * cf;
    return out;
  }

  // Field of the flat segment array exactly as given (no mirror coil).
  // Used by tests to validate the segment formula against analytic results.
  function rawSegmentsB(flat, nSegs, px, py, pz) {
    let bx = 0, by = 0, bz = 0;
    for (let i = 0; i < nSegs; i++) {
      const o = i * 7;
      const ax = flat[o], ay = flat[o + 1], az = flat[o + 2];
      const abx = flat[o + 3] - ax, aby = flat[o + 4] - ay, abz = flat[o + 5] - az;
      const r1x = px - ax, r1y = py - ay, r1z = pz - az;
      const r2x = r1x - abx, r2y = r1y - aby, r2z = r1z - abz;
      const cx = aby * r1z - abz * r1y;
      const cy = abz * r1x - abx * r1z;
      const cz = abx * r1y - aby * r1x;
      const c2 = cx * cx + cy * cy + cz * cz;
      if (c2 < 1e-20) continue;
      const r1m = Math.sqrt(r1x * r1x + r1y * r1y + r1z * r1z);
      const r2m = Math.sqrt(r2x * r2x + r2y * r2y + r2z * r2z);
      if (r1m < 1e-12 || r2m < 1e-12) continue;
      const k = (KB * flat[o + 6] *
        ((abx * r1x + aby * r1y + abz * r1z) / r1m -
         (abx * r2x + aby * r2y + abz * r2z) / r2m)) / c2;
      bx += k * cx; by += k * cy; bz += k * cz;
    }
    return { x: bx, y: by, z: bz };
  }

  // Bx only, for points on the x=0 plane (where By=Bz=0 by symmetry).
  // Right coil evaluated once, doubled.
  function fieldBxOnPlane(geo, py, pz) {
    const f = geo.flatInner, n = geo.nInner;
    let bx = 0;
    for (let i = 0; i < n; i++) {
      const o = i * 7;
      const ax = f[o], ay = f[o + 1], az = f[o + 2];
      const abx = f[o + 3] - ax, aby = f[o + 4] - ay, abz = f[o + 5] - az;
      const r1x = -ax, r1y = py - ay, r1z = pz - az;
      const r2x = r1x - abx, r2y = r1y - aby, r2z = r1z - abz;
      const cx = aby * r1z - abz * r1y;
      const cy = abz * r1x - abx * r1z;
      const cz = abx * r1y - aby * r1x;
      const c2 = cx * cx + cy * cy + cz * cz;
      if (c2 < 1e-20) continue;
      const r1m = Math.sqrt(r1x * r1x + r1y * r1y + r1z * r1z);
      const r2m = Math.sqrt(r2x * r2x + r2y * r2y + r2z * r2z);
      if (r1m < 1e-12 || r2m < 1e-12) continue;
      const k = (KB * f[o + 6] *
        ((abx * r1x + aby * r1y + abz * r1z) / r1m -
         (abx * r2x + aby * r2y + abz * r2z) / r2m)) / c2;
      bx += k * cx;
    }
    const bInner = 2 * bx; // x-mirror doubles Bx on this plane
    return bInner * coreFactor(geo, Math.abs(bInner));
  }

  // ---- beam ---------------------------------------------------------------
  function beamKinematics(anodeKV) {
    const gamma = 1 + (anodeKV * 1e3 * QE) / (ME * C * C);
    const v = C * Math.sqrt(1 - 1 / (gamma * gamma));
    return { gamma, v };
  }

  // Relativistic Boris push, magnetic field only (E = 0 in the drift region).
  // Beam stays on the x = 0 plane, so B = (Bx, 0, 0) and motion is in y-z.
  function traceTrajectory(geo, opts) {
    const p = geo.params;
    const { gamma, v } = beamKinematics(p.anodeKV);
    const stepM = (opts && opts.stepM) || 0.5e-3; // ~0.5 mm per step
    const dt = stepM / v;
    const halfQdtGm = (-QE * dt) / (2 * gamma * ME); // q = -e
    const maxSteps = 20000;

    let py = 0, pz = p.gunZ;
    let vy = 0, vz = v;
    const pts = [pz, py];
    let maxAbsB = 0, time = 0;

    for (let s = 0; s < maxSteps && pz < p.screenZ; s++) {
      const bx = fieldBxOnPlane(geo, py, pz);
      if (Math.abs(bx) > Math.abs(maxAbsB)) maxAbsB = bx;
      // Boris rotation about x: t = (q dt / 2 gamma m) * B
      const t = halfQdtGm * bx;
      // v' = v + v x t ; with B along x and v in y-z: (v x t)_y = vz*t...
      // v x B for B=(b,0,0): (vy,vz) -> (vy*0... ) compute directly:
      // v x t where t=(t,0,0): (vx,vy,vz)x(t,0,0) = (vy*0-vz*0, vz*t-0, 0-vy*t)
      const wy = vy + vz * t;
      const wz = vz - vy * t;
      const s2 = (2 * t) / (1 + t * t);
      vy = vy + wz * s2;
      vz = vz - wy * s2;
      py += vy * dt;
      pz += vz * dt;
      time += dt;
      pts.push(pz, py);
    }
    const speed = Math.sqrt(vy * vy + vz * vz);
    return {
      points: pts,                 // interleaved z,y pairs
      deflY: py,                   // m at screen plane
      angleDeg: (Math.atan2(vy, vz) * 180) / Math.PI,
      maxBx: maxAbsB,
      transitNs: time * 1e9,
      speedErr: (speed - v) / v,   // Boris should keep this ~0
      gamma, v,
    };
  }

  // Deflection vs current sweep (rebuilds amp only — geometry is linear in I,
  // but trajectory is not, so each point is a full trace).
  function sweepCurrent(params, currents) {
    const out = [];
    for (const I of currents) {
      const geo = buildGeometry(Object.assign({}, params, { current: I }));
      const t = traceTrajectory(geo);
      out.push({ current: I, deflY: t.deflY, angleDeg: t.angleDeg });
    }
    return out;
  }

  return {
    MU0, QE, ME, C, defaults,
    buildGeometry, fieldB, fieldBxOnPlane, rawSegmentsB, coreFactor,
    beamKinematics, traceTrajectory, sweepCurrent,
  };
});
