// Minimal canvas plotting (line plots and heatmaps) that follows the page theme.
"use strict";

const Plot = (() => {
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  function setup(canvas) {
    const dpr = window.devicePixelRatio || 1;
    const w = canvas.clientWidth, h = canvas.clientHeight;
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    return { ctx, w, h };
  }

  function niceTicks(lo, hi, count = 5) {
    const span = hi - lo || 1;
    const step0 = span / count, mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => span / s <= count) || 10 * mag;
    const ticks = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) ticks.push(+v.toFixed(10));
    return ticks;
  }

  function fmt(v) {
    const a = Math.abs(v);
    if (a >= 1000) return (v / 1000).toFixed(a >= 10000 ? 0 : 1) + "k";
    return a >= 10 || a === 0 ? String(Math.round(v)) : v.toFixed(1);
  }

  function frame(ctx, box, x0, x1, y0, y1, xlabel, ylabel) {
    const grid = css("--grid"), muted = css("--muted");
    ctx.font = "11px Inter, system-ui, sans-serif";
    ctx.fillStyle = muted;
    ctx.strokeStyle = grid;
    ctx.lineWidth = 1;
    for (const t of niceTicks(x0, x1)) {
      const x = box.l + ((t - x0) / (x1 - x0)) * box.w;
      ctx.beginPath(); ctx.moveTo(x, box.t); ctx.lineTo(x, box.t + box.h); ctx.stroke();
      ctx.textAlign = "center"; ctx.fillText(fmt(t), x, box.t + box.h + 14);
    }
    for (const t of niceTicks(y0, y1, 4)) {
      const y = box.t + box.h - ((t - y0) / (y1 - y0)) * box.h;
      ctx.beginPath(); ctx.moveTo(box.l, y); ctx.lineTo(box.l + box.w, y); ctx.stroke();
      ctx.textAlign = "right"; ctx.fillText(fmt(t), box.l - 6, y + 4);
    }
    ctx.textAlign = "center";
    ctx.fillText(xlabel, box.l + box.w / 2, box.t + box.h + 30);
    ctx.save(); ctx.translate(12, box.t + box.h / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(ylabel, 0, 0); ctx.restore();
  }

  // series: [{x: array, y: array, color, width}]
  function line(canvas, series, opts) {
    const { ctx, w, h } = setup(canvas);
    const box = { l: 52, t: 10, w: w - 64, h: h - 50 };
    let [x0, x1] = opts.xlim || [Infinity, -Infinity], [y0, y1] = opts.ylim || [Infinity, -Infinity];
    if (!opts.xlim) for (const s of series) for (const v of s.x) { x0 = Math.min(x0, v); x1 = Math.max(x1, v); }
    if (!opts.ylim) {
      for (const s of series) for (const v of s.y) { y0 = Math.min(y0, v); y1 = Math.max(y1, v); }
      const pad = (y1 - y0) * 0.08 || 1; y0 -= pad; y1 += pad;
    }
    frame(ctx, box, x0, x1, y0, y1, opts.xlabel || "", opts.ylabel || "");
    ctx.save();
    ctx.beginPath(); ctx.rect(box.l, box.t, box.w, box.h); ctx.clip();
    for (const s of series) {
      ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 1.5; ctx.beginPath();
      const step = Math.max(1, Math.floor(s.x.length / (box.w * 2)));
      for (let i = 0; i < s.x.length; i += step) {
        let yv = s.y[i];
        if (step > 1) for (let j = i; j < Math.min(i + step, s.x.length); j++) yv = Math.max(yv, s.y[j]);  // keep peaks
        const x = box.l + ((s.x[i] - x0) / (x1 - x0)) * box.w;
        const y = box.t + box.h - ((yv - y0) / (y1 - y0)) * box.h;
        i ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
      }
      ctx.stroke();
    }
    ctx.restore();
    return { box, x0, x1, y0, y1 };
  }

  // Perceptual-ish dark-to-bright colour map (viridis stops).
  const STOPS = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];
  function cmap(t) {
    t = Math.min(1, Math.max(0, t)) * (STOPS.length - 1);
    const i = Math.min(STOPS.length - 2, Math.floor(t)), f = t - i;
    return STOPS[i].map((c, j) => Math.round(c + f * (STOPS[i + 1][j] - c)));
  }

  // z[row][col] with rows = y (bottom to top), cols = x. Values already in dB.
  function heatmap(canvas, z, opts) {
    const { ctx, w, h } = setup(canvas);
    const box = { l: 52, t: 10, w: w - 64, h: h - 50 };
    const rows = z.length, cols = z[0].length;
    const img = ctx.createImageData(cols, rows);
    for (let r = 0; r < rows; r++)
      for (let c = 0; c < cols; c++) {
        const [R, G, B] = cmap((z[r][c] - opts.zmin) / (opts.zmax - opts.zmin));
        const o = ((rows - 1 - r) * cols + c) * 4;
        img.data[o] = R; img.data[o + 1] = G; img.data[o + 2] = B; img.data[o + 3] = 255;
      }
    const off = document.createElement("canvas");
    off.width = cols; off.height = rows;
    off.getContext("2d").putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(off, box.l, box.t, box.w, box.h);
    frame(ctx, box, opts.x0, opts.x1, opts.y0, opts.y1, opts.xlabel, opts.ylabel);
    return { ctx, box };
  }

  return { line, heatmap, css };
})();
