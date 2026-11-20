/* Canvas plotting helpers (HiDPI, nice ticks, line / bar / heatmap charts).
   All chart functions take a <canvas> whose CSS width/height define the size. */
(function () {
  'use strict';
  const QM = (window.QM = window.QM || {});

  function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function theme() {
    return {
      ink: css('--ink') || '#1d2433', muted: css('--muted') || '#6b7280', grid: css('--grid') || '#e5e7eb',
      axis: css('--axis') || '#9ca3af', bg: css('--panel') || '#fff', accent: css('--accent') || '#175ec5',
      series: [1, 2, 3, 4, 5, 6, 7, 8].map((i) => css('--s' + i)).filter(Boolean),
    };
  }

  function setup(canvas) {
    const dpr = window.devicePixelRatio || 1;
    const w = Math.max(10, canvas.clientWidth), h = Math.max(10, canvas.clientHeight);
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.font = '12px system-ui, sans-serif';
    return { ctx, w, h };
  }

  function niceStep(range, target = 6) {
    const raw = Math.abs(range) / target || 1, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
    return (f < 1.5 ? 1 : f < 2.25 ? 2 : f < 3.5 ? 2.5 : f < 7.5 ? 5 : 10) * p;
  }
  function tickLabel(v, step) {
    const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9) + (step / Math.pow(10, Math.floor(Math.log10(step))) === 2.5 ? 1 : 0));
    return QM.fmt(+v.toFixed(Math.min(d, 6)), Math.min(d, 6));
  }
  function piLabel(k, den) { // k/den * π
    const g = gcd(Math.abs(k), den), a = k / g, b = den / g;
    if (a === 0) return '0';
    const num = (a === 1 ? '' : a === -1 ? '−' : String(a).replace('-', '−')) + 'π';
    return b === 1 ? num : num + '/' + b;
  }
  function gcd(a, b) { return b ? gcd(b, a % b) : a; }

  /** draws frame/grid/ticks; returns mapping {X(x), Y(y), box} */
  function axes(ctx, w, h, o) {
    const t = theme();
    const pad = Object.assign({ l: 46, r: 12, t: 12, b: o.xlabel ? 34 : 22 }, o.pad || {});
    const box = { x: pad.l, y: pad.t, w: w - pad.l - pad.r, h: h - pad.t - pad.b };
    const X = (x) => box.x + ((x - o.xMin) / (o.xMax - o.xMin)) * box.w;
    const Y = (y) => box.y + box.h - ((y - o.yMin) / (o.yMax - o.yMin)) * box.h;
    ctx.save();
    ctx.strokeStyle = t.grid; ctx.fillStyle = t.muted; ctx.lineWidth = 1;
    // x ticks
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    if (o.piTicks) {
      const span = (o.xMax - o.xMin) / Math.PI; const den = span <= 2.5 ? 4 : span <= 5 ? 2 : 1;
      const stepK = Math.max(1, Math.round((span * den) / 8));
      for (let k = Math.ceil((o.xMin / Math.PI) * den); k <= (o.xMax / Math.PI) * den + 1e-9; k += 1) {
        if (k % stepK) continue; const xp = X((k * Math.PI) / den);
        ctx.beginPath(); ctx.moveTo(xp, box.y); ctx.lineTo(xp, box.y + box.h); ctx.stroke();
        ctx.fillText(piLabel(k, den), xp, box.y + box.h + 4);
      }
    } else if (!o.noXTicks) {
      const s = o.xStep || niceStep(o.xMax - o.xMin, Math.max(3, box.w / 70));
      for (let v = Math.ceil(o.xMin / s) * s; v <= o.xMax + s * 1e-9; v += s) {
        const xp = X(v); ctx.beginPath(); ctx.moveTo(xp, box.y); ctx.lineTo(xp, box.y + box.h); ctx.stroke();
        ctx.fillText(tickLabel(v, s), xp, box.y + box.h + 4);
      }
    }
    // y ticks
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    const sy = o.yStep || niceStep(o.yMax - o.yMin, Math.max(2, box.h / 40));
    for (let v = Math.ceil(o.yMin / sy) * sy; v <= o.yMax + sy * 1e-9; v += sy) {
      const yp = Y(v); ctx.beginPath(); ctx.moveTo(box.x, yp); ctx.lineTo(box.x + box.w, yp); ctx.stroke();
      ctx.fillText(tickLabel(v, sy), box.x - 5, yp);
    }
    // zero lines
    ctx.strokeStyle = t.axis;
    if (o.yMin < 0 && o.yMax > 0) { ctx.beginPath(); ctx.moveTo(box.x, Y(0)); ctx.lineTo(box.x + box.w, Y(0)); ctx.stroke(); }
    ctx.strokeRect(box.x, box.y, box.w, box.h);
    ctx.fillStyle = t.ink;
    if (o.xlabel) { ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; ctx.fillText(o.xlabel, box.x + box.w, h - 2); }
    if (o.ylabel) { ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText(o.ylabel, box.x + 4, box.y + 3); }
    ctx.restore();
    return { X, Y, box };
  }

  /** opts: {xMin,xMax,yMin,yMax, series:[{x,y,color,width,dash,fill,label}], vlines:[{x,color,dash,label}], shade:[{from,to,color}], legend, xlabel, ylabel, piTicks} */
  function line(canvas, o) {
    const { ctx, w, h } = setup(canvas), t = theme();
    if (o.yMin === undefined || o.yMax === undefined) {
      let lo = Infinity, hi = -Infinity; o.series.forEach((s) => s.y.forEach((v) => { if (isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } }));
      if (!isFinite(lo)) { lo = 0; hi = 1; } if (hi - lo < 1e-9) { hi += 0.5; lo -= 0.5; }
      const m = 0.06 * (hi - lo); if (o.yMin === undefined) o.yMin = lo - m; if (o.yMax === undefined) o.yMax = hi + m;
    }
    const { X, Y, box } = axes(ctx, w, h, o);
    ctx.save(); ctx.beginPath(); ctx.rect(box.x, box.y, box.w, box.h); ctx.clip();
    (o.shade || []).forEach((s) => { ctx.fillStyle = s.color || 'rgba(23,94,197,.12)'; const a = X(Math.max(s.from, o.xMin)), b = X(Math.min(s.to, o.xMax)); ctx.fillRect(a, box.y, b - a, box.h); });
    o.series.forEach((s, si) => {
      const col = s.color || t.series[si % t.series.length];
      const xs = s.x, ys = s.y;
      if (s.fill) {
        ctx.beginPath(); ctx.moveTo(X(xs[0]), Y(Math.max(o.yMin, 0)));
        for (let i = 0; i < xs.length; i++) ctx.lineTo(X(xs[i]), Y(ys[i]));
        ctx.lineTo(X(xs[xs.length - 1]), Y(Math.max(o.yMin, 0))); ctx.closePath();
        ctx.fillStyle = s.fill === true ? hexA(col, 0.18) : s.fill; ctx.fill();
      }
      if (s.fillRange) { // fill under curve only for x in [a,b]
        const [a, b] = s.fillRange; ctx.beginPath(); let started = false;
        for (let i = 0; i < xs.length; i++) if (xs[i] >= a && xs[i] <= b) { if (!started) { ctx.moveTo(X(xs[i]), Y(0)); started = true; } ctx.lineTo(X(xs[i]), Y(ys[i])); }
        if (started) { let last = xs.length - 1; while (last > 0 && xs[last] > b) last--; ctx.lineTo(X(xs[last]), Y(0)); ctx.closePath(); ctx.fillStyle = hexA(col, 0.35); ctx.fill(); }
      }
      ctx.beginPath(); ctx.strokeStyle = col; ctx.lineWidth = s.width || 2; ctx.setLineDash(s.dash || []);
      let pen = false;
      for (let i = 0; i < xs.length; i++) { const v = ys[i]; if (!isFinite(v)) { pen = false; continue; } const px = X(xs[i]), py = Y(v); pen ? ctx.lineTo(px, py) : ctx.moveTo(px, py); pen = true; }
      ctx.stroke(); ctx.setLineDash([]);
      if (s.points) { ctx.fillStyle = col; for (let i = 0; i < xs.length; i++) { ctx.beginPath(); ctx.arc(X(xs[i]), Y(ys[i]), 2.5, 0, 7); ctx.fill(); } }
    });
    (o.vlines || []).forEach((v) => { ctx.strokeStyle = v.color || t.ink; ctx.setLineDash(v.dash || [4, 4]); ctx.lineWidth = 1.2; ctx.beginPath(); ctx.moveTo(X(v.x), box.y); ctx.lineTo(X(v.x), box.y + box.h); ctx.stroke(); ctx.setLineDash([]); });
    (o.hlines || []).forEach((v) => { ctx.strokeStyle = v.color || t.muted; ctx.setLineDash(v.dash || [4, 4]); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(box.x, Y(v.y)); ctx.lineTo(box.x + box.w, Y(v.y)); ctx.stroke(); ctx.setLineDash([]); });
    ctx.restore();
    if (o.legend !== false) legend(ctx, box, o.series.filter((s) => s.label).map((s) => ({ label: s.label, color: s.color || t.series[o.series.indexOf(s) % t.series.length] })));
    return { X, Y, box };
  }

  function legend(ctx, box, items) {
    if (!items.length) return; const t = theme();
    ctx.save(); ctx.font = '12px system-ui, sans-serif'; ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
    let x = box.x + box.w - 8, y = box.y + 12;
    const widths = items.map((i) => ctx.measureText(i.label).width + 26);
    const total = widths.reduce((a, b) => a + b, 0);
    if (total < box.w * 0.9) { x = box.x + box.w - total - 4; items.forEach((it, i) => { ctx.fillStyle = it.color; ctx.fillRect(x, y - 2, 14, 4); ctx.fillStyle = t.ink; ctx.fillText(it.label, x + 18, y); x += widths[i]; }); }
    else items.forEach((it, i) => { const yy = y + i * 16, xx = box.x + box.w - Math.max(...widths); ctx.fillStyle = it.color; ctx.fillRect(xx, yy - 2, 14, 4); ctx.fillStyle = t.ink; ctx.fillText(it.label, xx + 18, yy); });
    ctx.restore();
  }

  /** opts: {labels, values, colors?, phases? (radians -> hue), yMax?, yMin?, valueLabels?, xlabel?, ylabel?, highlight? index} */
  function bars(canvas, o) {
    const { ctx, w, h } = setup(canvas), t = theme();
    const n = o.values.length;
    const yMax = o.yMax !== undefined ? o.yMax : Math.max(1e-9, ...o.values) * 1.12;
    const yMin = o.yMin !== undefined ? o.yMin : Math.min(0, ...o.values);
    const rot = n > 10;
    const { X, Y, box } = axes(ctx, w, h, { xMin: 0, xMax: n, yMin, yMax, noXTicks: true, ylabel: o.ylabel, xlabel: o.xlabel, pad: { b: rot ? 48 : (o.xlabel ? 38 : 26) } });
    const bw = box.w / n;
    ctx.save();
    for (let i = 0; i < n; i++) {
      const v = o.values[i];
      ctx.fillStyle = o.phases ? phaseColor(o.phases[i]) : (o.colors ? o.colors[i % o.colors.length] : t.accent);
      if (o.highlight === i) { ctx.fillStyle = t.ink; }
      const y0 = Y(Math.max(0, yMin)), y1 = Y(v);
      ctx.fillRect(X(i) + bw * 0.15, Math.min(y0, y1), bw * 0.7, Math.abs(y1 - y0));
      ctx.fillStyle = t.ink; ctx.textAlign = 'center';
      if (o.labels && (n <= 40 || i % Math.ceil(n / 40) === 0)) {
        ctx.save(); ctx.translate(X(i + 0.5), box.y + box.h + 4);
        if (rot) { ctx.rotate(-Math.PI / 4); ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; } else ctx.textBaseline = 'top';
        ctx.fillText(o.labels[i], 0, 0); ctx.restore();
      }
      if (o.valueLabels && n <= 16) { ctx.textBaseline = 'bottom'; ctx.fillStyle = t.muted; ctx.fillText(QM.fmt(v, 3), X(i + 0.5), Math.min(y0, y1) - 2); }
    }
    ctx.restore();
    return { X, Y, box };
  }

  /** opts: {data (Float32Array nx*ny, row = y index from bottom), nx, ny, xMin,xMax,yMin,yMax, vMax?, xlabel, ylabel, piTicks, marker:{x,y}, hline} */
  function heatmap(canvas, o) {
    const { ctx, w, h } = setup(canvas), t = theme();
    const { X, Y, box } = axes(ctx, w, h, o);
    const off = document.createElement('canvas'); off.width = o.nx; off.height = o.ny;
    const octx = off.getContext('2d'), img = octx.createImageData(o.nx, o.ny);
    let vMax = o.vMax; if (!vMax) { vMax = 0; for (const v of o.data) vMax = Math.max(vMax, v); vMax = vMax || 1; }
    for (let j = 0; j < o.ny; j++) for (let i = 0; i < o.nx; i++) {
      const v = Math.min(1, Math.max(0, o.data[j * o.nx + i] / vMax)), [r, g, b] = ramp(v);
      const k = ((o.ny - 1 - j) * o.nx + i) * 4; img.data[k] = r; img.data[k + 1] = g; img.data[k + 2] = b; img.data[k + 3] = 255;
    }
    octx.putImageData(img, 0, 0);
    ctx.save(); ctx.imageSmoothingEnabled = o.nx < 64 ? false : true;
    ctx.drawImage(off, box.x, box.y, box.w, box.h);
    if (o.hline !== undefined) { ctx.strokeStyle = '#fff'; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(box.x, Y(o.hline)); ctx.lineTo(box.x + box.w, Y(o.hline)); ctx.stroke(); ctx.setLineDash([]); }
    if (o.marker) { ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(X(o.marker.x), Y(o.marker.y), 5, 0, 7); ctx.stroke(); }
    ctx.restore();
    return { X, Y, box };
  }
  // perceptual-ish dark→yellow ramp (similar to rotor sim)
  function ramp(v) {
    v = Math.pow(v, 0.64);
    const stops = [[13, 17, 40], [52, 33, 110], [150, 40, 120], [230, 90, 70], [252, 200, 60], [255, 250, 210]];
    const x = v * (stops.length - 1), i = Math.min(stops.length - 2, Math.floor(x)), f = x - i;
    return stops[i].map((c, k) => Math.round(c + (stops[i + 1][k] - c) * f));
  }
  function phaseColor(arg, light = 50) { const deg = ((arg * 180) / Math.PI + 360) % 360; return `hsl(${deg.toFixed(0)}, 75%, ${light}%)`; }
  function hexA(col, a) {
    if (col.startsWith('#') && (col.length === 7)) { const r = parseInt(col.slice(1, 3), 16), g = parseInt(col.slice(3, 5), 16), b = parseInt(col.slice(5, 7), 16); return `rgba(${r},${g},${b},${a})`; }
    if (col.startsWith('hsl(')) return col.replace('hsl(', 'hsla(').replace(')', `, ${a})`);
    return col;
  }
  /** small phase wheel legend */
  function phaseWheel(ctx, cx, cy, r) {
    for (let a = 0; a < 360; a += 4) { ctx.beginPath(); ctx.moveTo(cx, cy); ctx.arc(cx, cy, r, (a * Math.PI) / 180, ((a + 5) * Math.PI) / 180); ctx.closePath(); ctx.fillStyle = phaseColor((-a * Math.PI) / 180); ctx.fill(); }
  }

  QM.plot = { setup, axes, line, bars, heatmap, niceStep, phaseColor, phaseWheel, ramp, theme, hexA, legend };
})();
