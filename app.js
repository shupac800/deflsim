/* UI layer: controls, rendering, interaction. Physics lives in physics.js. */
"use strict";
(function () {
  const P = window.DeflSim;
  const state = Object.assign({}, P.defaults);

  // ---- palette from CSS custom properties (re-read on scheme change) ------
  let pal = {};
  function readPalette() {
    const cs = getComputedStyle(document.documentElement);
    const g = (n) => cs.getPropertyValue(n).trim();
    pal = {
      surface: g("--surface-1"), page: g("--page"),
      ink: g("--text-primary"), ink2: g("--text-secondary"), muted: g("--text-muted"),
      grid: g("--grid-line"), baseline: g("--baseline"),
      beam: g("--series-1"), coil: g("--series-2"),
      divMid: g("--div-mid"), divBlue: g("--div-blue"), divRed: g("--div-red"),
    };
  }
  readPalette();

  // sRGB <-> linear-light interpolation for the diverging ramp
  function hex2rgb(h) {
    const s = h.replace("#", "");
    return [0, 1, 2].map((i) => parseInt(s.slice(i * 2, i * 2 + 2), 16) / 255);
  }
  const s2l = (c) => (c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4));
  const l2s = (c) => (c <= 0.0031308 ? c * 12.92 : 1.055 * Math.pow(c, 1 / 2.4) - 0.055);
  function makeRamp(midHex, poleHex, n) {
    const m = hex2rgb(midHex).map(s2l), p = hex2rgb(poleHex).map(s2l);
    const out = [];
    for (let i = 0; i < n; i++) {
      const t = i / (n - 1);
      out.push([0, 1, 2].map((k) => Math.round(255 * l2s(m[k] + (p[k] - m[k]) * t))));
    }
    return out;
  }
  let rampNeg, rampPos; // Bx<0 -> blue arm, Bx>0 -> red arm
  function rebuildRamps() {
    rampNeg = makeRamp(pal.divMid, pal.divBlue, 64);
    rampPos = makeRamp(pal.divMid, pal.divRed, 64);
  }
  rebuildRamps();
  function divColor(t) { // t in [-1,1]
    const a = Math.min(1, Math.pow(Math.abs(t), 0.6));
    const ramp = t < 0 ? rampNeg : rampPos;
    return ramp[Math.round(a * 63)];
  }

  // ---- canvas helpers ------------------------------------------------------
  function setupCanvas(el) {
    const dpr = window.devicePixelRatio || 1;
    const w = el.clientWidth, h = el.clientHeight;
    if (el.width !== Math.round(w * dpr) || el.height !== Math.round(h * dpr)) {
      el.width = Math.round(w * dpr); el.height = Math.round(h * dpr);
    }
    const ctx = el.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h };
  }
  const $ = (id) => document.getElementById(id);
  const sideCv = $("side-canvas"), geoCv = $("geo-canvas"), chartCv = $("chart-canvas");

  // ---- simulation results (cached between renders) -------------------------
  let geo = null, traj = null;
  let grid = null;   // {data: Float32Array, nz, ny, z0, z1, y0, y1, bmax}
  let sweep = null;  // [{current, deflY, angleDeg}]
  let heatCanvas = document.createElement("canvas");

  function computeFast() {
    geo = P.buildGeometry(state);
    traj = P.traceTrajectory(geo);
  }

  function computeGrid() {
    const nz = 150, nyHalf = 32, ny = nyHalf * 2 + 1;
    const zPad = 0.015;
    const z0 = state.gunZ - zPad, z1 = state.screenZ + zPad;
    const ymax = Math.max(state.majorRadius * 1.25, Math.abs(traj.deflY) * 1.15, 0.05);
    const data = new Float32Array(nz * ny);
    let bmax = 0;
    for (let iy = 0; iy <= nyHalf; iy++) {
      const y = (iy / nyHalf) * ymax;
      for (let iz = 0; iz < nz; iz++) {
        const z = z0 + ((z1 - z0) * iz) / (nz - 1);
        const b = P.fieldBxOnPlane(geo, y, z);
        data[(nyHalf + iy) * nz + iz] = b;
        data[(nyHalf - iy) * nz + iz] = b; // Bx is even in y on this plane
        if (Math.abs(b) > bmax) bmax = Math.abs(b);
      }
    }
    grid = { data, nz, ny, z0, z1, y0: -ymax, y1: ymax, bmax };
    // paint offscreen heat image at grid resolution
    heatCanvas.width = nz; heatCanvas.height = ny;
    const hctx = heatCanvas.getContext("2d");
    const img = hctx.createImageData(nz, ny);
    for (let i = 0; i < nz * ny; i++) {
      const c = divColor(bmax ? data[i] / bmax : 0);
      img.data[i * 4] = c[0]; img.data[i * 4 + 1] = c[1]; img.data[i * 4 + 2] = c[2];
      img.data[i * 4 + 3] = 255;
    }
    hctx.putImageData(img, 0, 0);
  }

  function computeSweep() {
    const n = 17, maxI = 6.0;
    const currents = [];
    for (let i = -(n - 1); i <= n - 1; i++) currents.push((maxI * i) / (n - 1));
    sweep = P.sweepCurrent(state, currents);
    const tb = $("sweep-tbody");
    tb.textContent = "";
    for (const s of sweep) {
      const tr = document.createElement("tr");
      for (const v of [s.current.toFixed(2), (s.deflY * 1000).toFixed(1), s.angleDeg.toFixed(2)]) {
        const td = document.createElement("td");
        td.textContent = v;
        tr.appendChild(td);
      }
      tb.appendChild(tr);
    }
  }

  // ---- stat tiles ----------------------------------------------------------
  function fmt(v, d) { return v.toFixed(d); }
  function updateTiles() {
    const dmm = traj.deflY * 1000;
    $("t-defl").textContent = fmt(dmm, 1) + " mm";
    $("t-defl-in").textContent = fmt(dmm / 25.4, 2) + " in at the screen";
    $("t-angle").textContent = fmt(traj.angleDeg, 1) + "°";
    $("t-bmax").textContent = fmt(traj.maxBx * 1000, 2) + " mT";
    const bCore = geo.crowd * Math.abs(traj.maxBx);
    $("t-bmax-sub").textContent =
      "|B" + "ₓ" + "| along beam · core ~" + fmt((bCore / state.coreBSat) * 100, 0) + "% of Bsat";
    $("t-vel").textContent = fmt((traj.v / P.C) * 100, 1) + "% c";
    $("t-vel-sub").textContent = (traj.v / 1e6).toFixed(1) + " Mm/s · transit " + traj.transitNs.toFixed(1) + " ns";
    $("t-sens").textContent = Math.abs(state.current) > 1e-6
      ? fmt(Math.abs(dmm) / Math.abs(state.current), 1) + " mm/A" : "—";
    $("i-driver").textContent = state.current.toFixed(2) + " A";
    $("i-quad").textContent = state.current.toFixed(2) + " A";
    $("i-loop").textContent = state.current.toFixed(2) + " A";
    $("ampturns").textContent =
      (state.loops * state.current).toFixed(2) + " A·turns";
  }

  // ---- side view ------------------------------------------------------------
  const sideMap = {}; // world<->pixel mapping for tooltip
  function renderSide() {
    const { ctx, w, h } = setupCanvas(sideCv);
    ctx.clearRect(0, 0, w, h);
    if (!grid) return;
    const mL = 46, mR = 12, mT = 8, mB = 26;
    const pw = w - mL - mR, ph = h - mT - mB;
    const zpx = (z) => mL + ((z - grid.z0) / (grid.z1 - grid.z0)) * pw;
    const ypx = (y) => mT + ((grid.y1 - y) / (grid.y1 - grid.y0)) * ph;
    Object.assign(sideMap, { zpx, ypx, mL, mT, pw, ph });

    // heat image
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(heatCanvas, mL, mT, pw, ph);

    // axes ticks (z in mm)
    ctx.font = "11px system-ui, sans-serif";
    ctx.fillStyle = pal.muted;
    ctx.strokeStyle = pal.baseline;
    ctx.lineWidth = 1;
    ctx.textAlign = "center";
    for (let zmm = 0; zmm <= grid.z1 * 1000; zmm += 50) {
      const x = zpx(zmm / 1000);
      ctx.fillText(zmm.toString(), x, h - 8);
    }
    ctx.save();
    ctx.textAlign = "right";
    for (const ymm of [-50, 0, 50]) {
      if (ymm / 1000 > grid.y1 || ymm / 1000 < grid.y0) continue;
      ctx.fillText(ymm.toString(), mL - 6, ypx(ymm / 1000) + 4);
    }
    ctx.restore();
    ctx.textAlign = "left";
    ctx.fillText("z (mm)", w - mR - 42, h - 8);
    ctx.save();
    ctx.translate(12, mT + 30); ctx.rotate(-Math.PI / 2);
    ctx.fillText("y (mm)", -30, 0);
    ctx.restore();

    // beam axis
    ctx.strokeStyle = pal.baseline;
    ctx.setLineDash([1, 3]);
    ctx.beginPath(); ctx.moveTo(mL, ypx(0)); ctx.lineTo(mL + pw, ypx(0)); ctx.stroke();
    ctx.setLineDash([]);

    // yoke winding silhouette (projection of the wire band onto the plane)
    const sin0 = Math.sin((state.firstAngle * Math.PI) / 180);
    const sin1 = Math.sin((state.lastAngle * Math.PI) / 180);
    ctx.fillStyle = pal.ink2;
    ctx.globalAlpha = 0.22;
    for (const sgn of [1, -1]) {
      ctx.beginPath();
      ctx.moveTo(zpx(0), ypx(sgn * state.minorRadius * sin0));
      ctx.lineTo(zpx(state.yokeLength), ypx(sgn * state.majorRadius * sin0));
      ctx.lineTo(zpx(state.yokeLength), ypx(sgn * state.majorRadius * sin1));
      ctx.lineTo(zpx(0), ypx(sgn * state.minorRadius * sin1));
      ctx.closePath(); ctx.fill();
    }
    ctx.globalAlpha = 1;

    // screen
    ctx.strokeStyle = pal.ink2;
    ctx.lineWidth = 2.5;
    ctx.beginPath(); ctx.moveTo(zpx(state.screenZ), mT); ctx.lineTo(zpx(state.screenZ), mT + ph); ctx.stroke();

    // trajectory (2px, beam color) with surface halo for legibility on the heat fill
    const pts = traj.points;
    for (const pass of [0, 1]) {
      ctx.strokeStyle = pass === 0 ? pal.surface : pal.beam;
      ctx.lineWidth = pass === 0 ? 4.5 : 2;
      ctx.lineJoin = "round"; ctx.lineCap = "round";
      ctx.beginPath();
      ctx.moveTo(zpx(pts[0]), ypx(pts[1]));
      for (let i = 2; i < pts.length; i += 2) ctx.lineTo(zpx(pts[i]), ypx(pts[i + 1]));
      ctx.stroke();
    }
    // endpoint marker + direct label
    const ex = zpx(pts[pts.length - 2]), ey = ypx(pts[pts.length - 1]);
    ctx.fillStyle = pal.surface;
    ctx.beginPath(); ctx.arc(ex, ey, 6.5, 0, 7); ctx.fill();
    ctx.fillStyle = pal.beam;
    ctx.beginPath(); ctx.arc(ex, ey, 4.5, 0, 7); ctx.fill();
    ctx.fillStyle = pal.ink;
    ctx.font = "600 12px system-ui, sans-serif";
    ctx.textAlign = "right";
    ctx.fillText((traj.deflY * 1000).toFixed(1) + " mm", ex - 10, ey + (traj.deflY < 0 ? 14 : -8));

    // colorbar (top-right corner)
    const cbW = 120, cbH = 8, cbX = mL + pw - cbW - 10, cbY = mT + 10;
    for (let i = 0; i < cbW; i++) {
      const c = divColor((i / (cbW - 1)) * 2 - 1);
      ctx.fillStyle = "rgb(" + c[0] + "," + c[1] + "," + c[2] + ")";
      ctx.fillRect(cbX + i, cbY, 1.2, cbH);
    }
    ctx.strokeStyle = pal.baseline; ctx.lineWidth = 1;
    ctx.strokeRect(cbX - 0.5, cbY - 0.5, cbW + 1, cbH + 1);
    ctx.fillStyle = pal.ink2;
    ctx.font = "10.5px system-ui, sans-serif";
    const bmT = (grid.bmax * 1000).toFixed(grid.bmax * 1000 >= 10 ? 0 : 1);
    ctx.textAlign = "left"; ctx.fillText("−" + bmT, cbX, cbY + cbH + 12);
    ctx.textAlign = "center"; ctx.fillText("0", cbX + cbW / 2, cbY + cbH + 12);
    ctx.textAlign = "right"; ctx.fillText("+" + bmT + " mT", cbX + cbW, cbY + cbH + 12);
    ctx.fillText("signed Bx component", cbX + cbW, cbY + cbH + 25);
    ctx.textAlign = "left";
  }

  // ---- legends --------------------------------------------------------------
  function buildLegends() {
    const mk = (color, label) => {
      const k = document.createElement("span"); k.className = "key";
      const s = document.createElement("span"); s.className = "swatch-line";
      s.style.borderTopColor = color;
      const t = document.createElement("span"); t.textContent = label;
      k.appendChild(s); k.appendChild(t);
      return k;
    };
    const sl = $("side-legend"); sl.textContent = "";
    sl.appendChild(mk(pal.beam, "electron beam"));
    sl.appendChild(mk(pal.ink2, "yoke winding band"));
    const note = document.createElement("span");
    note.textContent = "heat fill: Bx on the x = 0 plane";
    sl.appendChild(note);
    const gl = $("geo-legend"); gl.textContent = "";
    gl.appendChild(mk(pal.coil, "coil windings"));
    gl.appendChild(mk(pal.beam, "electron beam"));
    gl.appendChild(mk(pal.ink2, "screen"));
    const gnote = document.createElement("span"); gnote.id = "geo-note";
    gl.appendChild(gnote);
  }

  // ---- 3D geometry view -------------------------------------------------------
  const view = { yaw: -0.9, pitch: 0.35 };
  function renderGeo() {
    const { ctx, w, h } = setupCanvas(geoCv);
    ctx.clearRect(0, 0, w, h);
    if (!geo) return;
    const cy = Math.cos(view.yaw), sy = Math.sin(view.yaw);
    const cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
    const zMid = state.screenZ / 2;
    const span = Math.max(state.screenZ + 0.03, state.majorRadius * 3);
    const sr = Math.max(Math.abs(traj.deflY) * 1.3, state.majorRadius * 1.2); // screen half-size
    const scale = 0.85 * Math.min(w / span, h / (2.6 * sr));
    // yaw mixes x/z so the tube axis sweeps across the screen; pitch tilts it
    const proj2 = (x, y, z) => {
      const z0 = z - zMid;
      const x1 = cy * x + sy * z0;
      const zr = -sy * x + cy * z0;
      const y1 = cp * y - sp * zr;
      return [w / 2 + x1 * scale, h / 2 - y1 * scale];
    };

    // coil loops (right + mirrored left) — one stroked path per loop, for fidelity to the selected count
    ctx.strokeStyle = pal.coil;
    ctx.lineWidth = 1.2;
    ctx.globalAlpha = 0.85;
    for (let li = 0; li < geo.loops.length; li++) {
      const pts = geo.loops[li];
      for (const mirror of [1, -1]) {
        ctx.beginPath();
        let first = true;
        for (const p of pts) {
          const q = proj2(mirror * p[0], p[1], p[2]);
          if (first) { ctx.moveTo(q[0], q[1]); first = false; } else ctx.lineTo(q[0], q[1]);
        }
        ctx.closePath();
        ctx.stroke();
      }
    }
    ctx.globalAlpha = 1;
    $("geo-note").textContent = geo.loops.length + " loops per coil";

    // screen outline at z = screenZ
    ctx.strokeStyle = pal.ink2;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    const corners = [[-sr, -sr], [sr, -sr], [sr, sr], [-sr, sr]];
    corners.forEach((c, i) => {
      const q = proj2(c[0], c[1], state.screenZ);
      i === 0 ? ctx.moveTo(q[0], q[1]) : ctx.lineTo(q[0], q[1]);
    });
    ctx.closePath(); ctx.stroke();

    // beam trajectory (x = 0 plane)
    const pts = traj.points;
    ctx.strokeStyle = pal.beam;
    ctx.lineWidth = 2;
    ctx.lineJoin = "round"; ctx.lineCap = "round";
    ctx.beginPath();
    let q = proj2(0, pts[1], pts[0]);
    ctx.moveTo(q[0], q[1]);
    for (let i = 2; i < pts.length; i += 2) {
      q = proj2(0, pts[i + 1], pts[i]);
      ctx.lineTo(q[0], q[1]);
    }
    ctx.stroke();
    // spot on the screen
    q = proj2(0, traj.deflY, state.screenZ);
    ctx.fillStyle = pal.surface;
    ctx.beginPath(); ctx.arc(q[0], q[1], 6.5, 0, 7); ctx.fill();
    ctx.fillStyle = pal.beam;
    ctx.beginPath(); ctx.arc(q[0], q[1], 4.5, 0, 7); ctx.fill();

    // z-axis hint
    ctx.strokeStyle = pal.baseline;
    ctx.setLineDash([1, 3]);
    ctx.beginPath();
    q = proj2(0, 0, state.gunZ); ctx.moveTo(q[0], q[1]);
    q = proj2(0, 0, state.screenZ); ctx.lineTo(q[0], q[1]);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = pal.muted;
    ctx.font = "11px system-ui, sans-serif";
    const qg = proj2(0, 0, state.gunZ);
    ctx.fillText("gun", qg[0] - 24, qg[1] + 4);
  }

  // drag to rotate
  let dragging = false, lastX = 0, lastY = 0;
  geoCv.addEventListener("pointerdown", (e) => {
    dragging = true; lastX = e.clientX; lastY = e.clientY;
    geoCv.classList.add("dragging");
    geoCv.setPointerCapture(e.pointerId);
  });
  geoCv.addEventListener("pointermove", (e) => {
    if (!dragging) return;
    view.yaw += (e.clientX - lastX) * 0.008;
    view.pitch = Math.max(-1.4, Math.min(1.4, view.pitch + (e.clientY - lastY) * 0.008));
    lastX = e.clientX; lastY = e.clientY;
    renderGeo();
  });
  geoCv.addEventListener("pointerup", (e) => {
    dragging = false;
    geoCv.classList.remove("dragging");
    geoCv.releasePointerCapture(e.pointerId);
  });

  // ---- deflection vs current chart ---------------------------------------------
  const chartMap = {};
  function renderChart(hoverIdx) {
    const { ctx, w, h } = setupCanvas(chartCv);
    ctx.clearRect(0, 0, w, h);
    if (!sweep) return;
    const mL = 46, mR = 14, mT = 12, mB = 30;
    const pw = w - mL - mR, ph = h - mT - mB;
    const maxI = sweep[sweep.length - 1].current;
    const maxD = Math.max(1, ...sweep.map((s) => Math.abs(s.deflY) * 1000));
    // clean y ticks
    const rawStep = maxD / 4;
    const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
    const yStep = [1, 2, 5, 10].map((m) => m * mag).find((s) => maxD / s <= 5) || 10 * mag;
    const yTop = Math.ceil(maxD / yStep) * yStep;
    const xp = (I) => mL + ((I + maxI) / (2 * maxI)) * pw;
    const yp = (d) => mT + ph / 2 - (d / yTop) * (ph / 2);
    Object.assign(chartMap, { xp, yp, mL, mT, pw, ph, maxI });

    // gridlines + labels
    ctx.font = "11px system-ui, sans-serif";
    ctx.strokeStyle = pal.grid; ctx.lineWidth = 1;
    ctx.fillStyle = pal.muted;
    ctx.textAlign = "right";
    for (let d = -yTop; d <= yTop + 1e-9; d += yStep) {
      const y = Math.round(yp(d)) + 0.5;
      ctx.beginPath(); ctx.moveTo(mL, y); ctx.lineTo(mL + pw, y); ctx.stroke();
      ctx.fillText(Math.round(d).toString(), mL - 7, y + 4);
    }
    ctx.strokeStyle = pal.baseline;
    const zeroY = Math.round(yp(0)) + 0.5;
    ctx.beginPath(); ctx.moveTo(mL, zeroY); ctx.lineTo(mL + pw, zeroY); ctx.stroke();
    const zeroX = Math.round(xp(0)) + 0.5;
    ctx.beginPath(); ctx.moveTo(zeroX, mT); ctx.lineTo(zeroX, mT + ph); ctx.stroke();
    ctx.textAlign = "center";
    const AStep = 1.0;
    for (let I = 0; I <= maxI + 1e-9; I += AStep) {
      ctx.fillText(I.toFixed(0), xp(I), h - 12);
      if (I > 0) ctx.fillText((-I).toFixed(0), xp(-I), h - 12);
    }
    ctx.fillText("coil current (A)", mL + pw / 2, h - 0.5);
    ctx.save();
    ctx.translate(11, mT + ph / 2); ctx.rotate(-Math.PI / 2);
    ctx.fillText("deflection (mm)", 0, 0);
    ctx.restore();

    // series line
    ctx.strokeStyle = pal.beam; ctx.lineWidth = 2;
    ctx.lineJoin = "round"; ctx.lineCap = "round";
    ctx.beginPath();
    sweep.forEach((s, i) => {
      const x = xp(s.current), y = yp(s.deflY * 1000);
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    });
    ctx.stroke();

    // marker at the operating current (ring in surface color)
    const opD = traj.deflY * 1000;
    const opI = state.current;
    if (Math.abs(opI) <= maxI) {
      const x = xp(opI), y = yp(opD);
      ctx.fillStyle = pal.surface;
      ctx.beginPath(); ctx.arc(x, y, 6.5, 0, 7); ctx.fill();
      ctx.fillStyle = pal.beam;
      ctx.beginPath(); ctx.arc(x, y, 4.5, 0, 7); ctx.fill();
      ctx.fillStyle = pal.ink;
      ctx.font = "600 12px system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.fillText(opD.toFixed(1) + " mm", x + 9, y - 8);
    }

    // hover crosshair + lifted point
    if (hoverIdx != null) {
      const s = sweep[hoverIdx];
      const x = Math.round(xp(s.current)) + 0.5;
      ctx.strokeStyle = pal.baseline;
      ctx.beginPath(); ctx.moveTo(x, mT); ctx.lineTo(x, mT + ph); ctx.stroke();
      const y = yp(s.deflY * 1000);
      ctx.fillStyle = pal.surface;
      ctx.beginPath(); ctx.arc(x, y, 6, 0, 7); ctx.fill();
      ctx.fillStyle = pal.beam;
      ctx.beginPath(); ctx.arc(x, y, 4, 0, 7); ctx.fill();
    }
  }

  // chart tooltip (crosshair snaps to nearest current)
  const chartTip = $("chart-tip");
  chartCv.addEventListener("pointermove", (e) => {
    if (!sweep) return;
    const r = chartCv.getBoundingClientRect();
    const frac = (e.clientX - r.left - chartMap.mL) / chartMap.pw;
    const idx = Math.max(0, Math.min(sweep.length - 1, Math.round(frac * (sweep.length - 1))));
    renderChart(idx);
    const s = sweep[idx];
    chartTip.textContent = "";
    const v = document.createElement("div"); v.className = "tt-val";
    v.textContent = (s.deflY * 1000).toFixed(1) + " mm";
    const row = document.createElement("div"); row.className = "tt-row";
    const key = document.createElement("span"); key.className = "tt-key";
    row.appendChild(key);
    row.appendChild(document.createTextNode(
      s.current.toFixed(2) + " A · " + s.angleDeg.toFixed(1) + "°"));
    chartTip.appendChild(v); chartTip.appendChild(row);
    chartTip.style.display = "block";
    const tx = chartMap.xp(s.current);
    chartTip.style.left = Math.min(tx + 12, r.width - 130) + "px";
    chartTip.style.top = "14px";
  });
  chartCv.addEventListener("pointerleave", () => {
    chartTip.style.display = "none";
    renderChart(null);
  });

  // side-view tooltip (bilinear sample of the field grid)
  const sideTip = $("side-tip");
  sideCv.addEventListener("pointermove", (e) => {
    if (!grid) return;
    const r = sideCv.getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    const fz = (px - sideMap.mL) / sideMap.pw, fy = (py - sideMap.mT) / sideMap.ph;
    if (fz < 0 || fz > 1 || fy < 0 || fy > 1) { sideTip.style.display = "none"; return; }
    const z = grid.z0 + fz * (grid.z1 - grid.z0);
    const y = grid.y1 - fy * (grid.y1 - grid.y0);
    const gx = fz * (grid.nz - 1), gy = (1 - (y - grid.y0) / (grid.y1 - grid.y0)) * (grid.ny - 1);
    const x0 = Math.floor(gx), y0 = Math.floor(gy);
    const x1 = Math.min(grid.nz - 1, x0 + 1), y1 = Math.min(grid.ny - 1, y0 + 1);
    const tx = gx - x0, ty = gy - y0;
    const d = grid.data, n = grid.nz;
    const b = d[y0 * n + x0] * (1 - tx) * (1 - ty) + d[y0 * n + x1] * tx * (1 - ty) +
              d[y1 * n + x0] * (1 - tx) * ty + d[y1 * n + x1] * tx * ty;
    sideTip.textContent = "";
    const v = document.createElement("div"); v.className = "tt-val";
    v.textContent = "Bx " + (b * 1000).toFixed(3) + " mT";
    const row = document.createElement("div"); row.className = "tt-row";
    row.textContent = "z " + (z * 1000).toFixed(0) + " mm · y " + (y * 1000).toFixed(0) + " mm";
    sideTip.appendChild(v); sideTip.appendChild(row);
    sideTip.style.display = "block";
    sideTip.style.left = Math.min(px + 14, r.width - 150) + "px";
    sideTip.style.top = Math.max(4, py - 44) + "px";
  });
  sideCv.addEventListener("pointerleave", () => { sideTip.style.display = "none"; });

  // ---- controls -------------------------------------------------------------
  const ctls = document.querySelectorAll(".ctl");
  function refreshCtl(el) {
    const key = el.dataset.key;
    const scale = parseFloat(el.dataset.scale || "1");
    const dec = parseInt(el.dataset.dec, 10);
    const val = (state[key] * scale).toFixed(dec);
    el.querySelector("output").innerHTML = el.dataset.prefix ? el.dataset.prefix + val : val + (el.dataset.unit || "");
    el.querySelector("input").value = state[key];
  }
  ctls.forEach((el) => {
    const inp = el.querySelector("input");
    inp.min = el.dataset.min; inp.max = el.dataset.max; inp.step = el.dataset.step;
    inp.addEventListener("input", () => {
      let v = parseFloat(inp.value);
      const key = el.dataset.key;
      // geometric consistency: keep minor < major, first < last
      if (key === "minorRadius") v = Math.min(v, state.majorRadius - 0.005);
      if (key === "majorRadius") v = Math.max(v, state.minorRadius + 0.005);
      if (key === "firstAngle") v = Math.min(v, state.lastAngle - 5);
      if (key === "lastAngle") v = Math.max(v, state.firstAngle + 5);
      state[key] = v;
      refreshCtl(el);
      update(false);
    });
    inp.addEventListener("change", () => update(true));
  });
  $("reset").addEventListener("click", () => {
    Object.assign(state, P.defaults);
    ctls.forEach(refreshCtl);
    update(true);
  });

  // ---- pipeline ---------------------------------------------------------------
  // fast path every input event; heavy grid + sweep debounced
  let heavyTimer = null;
  function update(immediateHeavy) {
    computeFast();
    updateTiles();
    renderGeo();
    renderSide(); // stale heat image under fresh trajectory until grid lands
    renderChart(null);
    clearTimeout(heavyTimer);
    heavyTimer = setTimeout(() => {
      computeGrid();
      renderSide();
      computeSweep();
      renderChart(null);
    }, immediateHeavy ? 0 : 160);
  }

  window.addEventListener("resize", () => { renderSide(); renderGeo(); renderChart(null); });
  const mq = window.matchMedia("(prefers-color-scheme: dark)");
  mq.addEventListener("change", () => {
    readPalette(); rebuildRamps(); buildLegends();
    if (grid) computeGrid(); // repaint heat image with the new ramp
    renderSide(); renderGeo(); renderChart(null); updateTiles();
  });

  // ---- init ---------------------------------------------------------------------
  buildLegends();
  ctls.forEach(refreshCtl);
  computeFast();
  updateTiles();
  computeGrid();
  computeSweep();
  renderSide(); renderGeo(); renderChart(null);
})();
