/* Bloch sphere / ball widget (orthographic, drag to rotate).
   const b = QM.bloch(canvas, {labels}); b.draw({vectors:[{v:[x,y,z], color, label}], axis:{v, color, label}, trail:[[x,y,z]...], points:[[x,y,z]...]}) */
(function () {
  'use strict';
  const QM = (window.QM = window.QM || {});

  function bloch(canvas, opts = {}) {
    let yaw = opts.yaw !== undefined ? opts.yaw : -0.55, pitch = opts.pitch !== undefined ? opts.pitch : 0.38;
    let lastState = null, drag = null;
    const labels = Object.assign({ x: ['x', '−x'], y: ['y', '−y'], z: ['|↑⟩ z', '|↓⟩'] }, opts.labels || {});

    canvas.style.cursor = 'grab'; canvas.style.touchAction = 'none';
    canvas.addEventListener('pointerdown', (e) => { drag = { x: e.clientX, y: e.clientY, yaw, pitch }; canvas.setPointerCapture(e.pointerId); canvas.style.cursor = 'grabbing'; });
    canvas.addEventListener('pointermove', (e) => { if (!drag) return; yaw = drag.yaw - (e.clientX - drag.x) * 0.01; pitch = Math.max(-1.4, Math.min(1.4, drag.pitch + (e.clientY - drag.y) * 0.01)); if (lastState) draw(lastState); });
    canvas.addEventListener('pointerup', () => { drag = null; canvas.style.cursor = 'grab'; });

    function proj(p, R, cx, cy) {
      const [x, y, z] = p;
      const x1 = x * Math.cos(yaw) - y * Math.sin(yaw), y1 = x * Math.sin(yaw) + y * Math.cos(yaw);
      const up = z * Math.cos(pitch) - x1 * Math.sin(pitch), depth = x1 * Math.cos(pitch) + z * Math.sin(pitch);
      return { X: cx + R * y1, Y: cy - R * up, d: depth };
    }

    function draw(state) {
      lastState = state;
      const { ctx, w, h } = QM.plot.setup(canvas), t = QM.plot.theme();
      const R = Math.min(w, h) * 0.36, cx = w / 2, cy = h / 2 + 4;
      const P = (p) => proj(p, R, cx, cy);
      // sphere
      const g = ctx.createRadialGradient(cx - R * 0.3, cy - R * 0.3, R * 0.1, cx, cy, R);
      g.addColorStop(0, QM.plot.hexA('#7aa7e8', 0.10)); g.addColorStop(1, QM.plot.hexA('#175ec5', 0.16));
      ctx.fillStyle = g; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 7); ctx.fill();
      ctx.strokeStyle = t.axis; ctx.lineWidth = 1; ctx.stroke();
      // great circles
      const circle = (f) => { let prev = null; for (let k = 0; k <= 96; k++) { const a = (k / 96) * 2 * Math.PI, q = P(f(a)); if (prev) { ctx.setLineDash(q.d < 0 ? [3, 4] : []); ctx.beginPath(); ctx.moveTo(prev.X, prev.Y); ctx.lineTo(q.X, q.Y); ctx.stroke(); } prev = q; } ctx.setLineDash([]); };
      ctx.strokeStyle = QM.plot.hexA('#6b7280', 0.6);
      circle((a) => [Math.cos(a), Math.sin(a), 0]);
      ctx.strokeStyle = QM.plot.hexA('#6b7280', 0.3);
      circle((a) => [Math.cos(a), 0, Math.sin(a)]);
      circle((a) => [0, Math.cos(a), Math.sin(a)]);
      // axes
      ctx.font = '12px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      [['x', [1, 0, 0]], ['y', [0, 1, 0]], ['z', [0, 0, 1]]].forEach(([k, a]) => {
        const p1 = P(a), p2 = P(a.map((c) => -c)), o = P([0, 0, 0]);
        ctx.strokeStyle = t.muted; ctx.setLineDash([2, 3]); ctx.beginPath(); ctx.moveTo(p2.X, p2.Y); ctx.lineTo(p1.X, p1.Y); ctx.stroke(); ctx.setLineDash([]);
        ctx.fillStyle = t.ink; const l1 = P(a.map((c) => c * 1.2)), l2 = P(a.map((c) => -c * 1.2));
        ctx.fillText(labels[k][0], l1.X, l1.Y); ctx.fillStyle = t.muted; ctx.fillText(labels[k][1], l2.X, l2.Y);
        void o;
      });
      // ensemble points
      (state.points || []).forEach((p) => { const q = P(p); ctx.fillStyle = QM.plot.hexA(state.pointColor || '#d9480f', q.d < 0 ? 0.25 : 0.6); ctx.beginPath(); ctx.arc(q.X, q.Y, 2.2, 0, 7); ctx.fill(); });
      // trail
      if (state.trail && state.trail.length > 1) {
        ctx.lineWidth = 1.5;
        for (let i = 1; i < state.trail.length; i++) { const a = P(state.trail[i - 1]), b = P(state.trail[i]); ctx.strokeStyle = QM.plot.hexA(state.trailColor || '#175ec5', b.d < 0 ? 0.25 : 0.7); ctx.beginPath(); ctx.moveTo(a.X, a.Y); ctx.lineTo(b.X, b.Y); ctx.stroke(); }
      }
      // measurement axis
      if (state.axis) {
        const a = state.axis.v, p1 = P(a), p2 = P(a.map((c) => -c));
        ctx.strokeStyle = state.axis.color || '#2f9e44'; ctx.lineWidth = 2; ctx.setLineDash([6, 4]);
        ctx.beginPath(); ctx.moveTo(p2.X, p2.Y); ctx.lineTo(p1.X, p1.Y); ctx.stroke(); ctx.setLineDash([]);
        ctx.fillStyle = state.axis.color || '#2f9e44'; const l = P(a.map((c) => c * 1.12)); ctx.fillText(state.axis.label || 'n', l.X, l.Y);
      }
      // vectors
      (state.vectors || []).forEach((vv) => {
        const o = P([0, 0, 0]), q = P(vv.v), len = Math.hypot(q.X - o.X, q.Y - o.Y);
        ctx.strokeStyle = ctx.fillStyle = vv.color || t.accent; ctx.lineWidth = 3;
        ctx.beginPath(); ctx.moveTo(o.X, o.Y); ctx.lineTo(q.X, q.Y); ctx.stroke();
        if (len > 4) { const ang = Math.atan2(q.Y - o.Y, q.X - o.X); ctx.beginPath(); ctx.moveTo(q.X, q.Y); ctx.lineTo(q.X - 10 * Math.cos(ang - 0.4), q.Y - 10 * Math.sin(ang - 0.4)); ctx.lineTo(q.X - 10 * Math.cos(ang + 0.4), q.Y - 10 * Math.sin(ang + 0.4)); ctx.closePath(); ctx.fill(); }
        ctx.beginPath(); ctx.arc(q.X, q.Y, 3.5, 0, 7); ctx.fill();
        if (vv.label) { ctx.textAlign = 'left'; ctx.fillText(vv.label, q.X + 7, q.Y - 8); ctx.textAlign = 'center'; }
      });
      if (state.title) { ctx.fillStyle = t.ink; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.font = '600 12px system-ui, sans-serif'; ctx.fillText(state.title, 8, 6); }
      if (state.caption) { ctx.fillStyle = t.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom'; ctx.font = '12px system-ui, sans-serif'; ctx.fillText(state.caption, 8, h - 4); }
    }
    return { draw, redraw: () => lastState && draw(lastState) };
  }

  /** Bloch vector of a normalised 2-component state */
  function blochOf(v) {
    const a = QM.V.get(v, 0), b = QM.V.get(v, 1);
    const ab = QM.C.mul(QM.C.conj(a), b);
    return [2 * ab.re, 2 * ab.im, QM.C.abs2(a) - QM.C.abs2(b)];
  }
  QM.bloch = bloch;
  QM.blochOf = blochOf;
})();
