// JavaScript port of cpp/include/radar/radar.hpp for the in-browser demo.
// Complex vectors are stored as separate Float64Array re/im pairs for speed.
"use strict";

const C = 299792458;

function nextPow2(n) { let p = 1; while (p < n) p <<= 1; return p; }

// In-place iterative radix-2 FFT on (re, im). inverse=true scales by 1/N.
function fft(re, im, inverse = false) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (2 * Math.PI / len) * (inverse ? 1 : -1);
    const wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      for (let k = 0; k < len / 2; k++) {
        const a = i + k, b = a + len / 2;
        const vr = re[b] * cr - im[b] * ci, vi = re[b] * ci + im[b] * cr;
        re[b] = re[a] - vr; im[b] = im[a] - vi;
        re[a] += vr; im[a] += vi;
        const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
      }
    }
  }
  if (inverse) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
}

function hamming(n) {
  const w = new Float64Array(n);
  for (let i = 0; i < n; i++) w[i] = n < 2 ? 1 : 0.54 - 0.46 * Math.cos((2 * Math.PI * i) / (n - 1));
  return w;
}

// Seeded RNG (mulberry32) so a scenario is reproducible.
function makeRng(seed) {
  let a = seed >>> 0;
  const uni = () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  let spare = null;
  const gauss = () => {  // Box-Muller
    if (spare !== null) { const s = spare; spare = null; return s; }
    const u = Math.max(uni(), 1e-12), v = uni();
    const r = Math.sqrt(-2 * Math.log(u));
    spare = r * Math.sin(2 * Math.PI * v);
    return r * Math.cos(2 * Math.PI * v);
  };
  return { uni, gauss };
}

function drawRcs(mean, model, rng) {
  const exp = (m) => -m * Math.log(Math.max(rng.uni(), 1e-12));
  if (model === 1 || model === 2) return exp(mean);
  if (model === 3 || model === 4) return exp(mean / 2) + exp(mean / 2);
  return mean;
}

function derived(p) {
  const lambda = C / p.fc;
  return {
    lambda,
    nPulse: Math.round(p.T * p.fs),
    ns: Math.round(p.pri * p.fs),
    rangeRes: C / (2 * p.B),
    maxRange: (C * (p.pri - p.T)) / 2,
    vRes: lambda / (2 * p.np * p.pri),
    vMax: lambda / (4 * p.pri),
  };
}

function lfmPulse(p) {
  const n = Math.round(p.T * p.fs), k = p.B / p.T;
  const re = new Float64Array(n), im = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const t = i / p.fs - p.T / 2;
    re[i] = Math.cos(Math.PI * k * t * t);
    im[i] = Math.sin(Math.PI * k * t * t);
  }
  return { re, im };
}

// Returns rx[m] = {re, im} for each pulse, exact delayed chirp + unit-power AWGN.
function simulate(p, targets, seed) {
  const d = derived(p), rng = makeRng(seed);
  const k = p.B / p.T, kSnr = Math.pow(10, p.refSnrDb / 10) * Math.pow(p.refRange, 4);
  const rx = [];
  for (let m = 0; m < p.np; m++) rx.push({ re: new Float64Array(d.ns), im: new Float64Array(d.ns) });
  for (const tg of targets) {
    const tau = (2 * tg.range) / C, fd = (2 * tg.velocity) / d.lambda;
    const scan = drawRcs(tg.rcs, tg.model, rng);
    const p2p = tg.model === 2 || tg.model === 4;
    const first = Math.max(0, Math.ceil(tau * p.fs)), last = Math.min(d.ns - 1, Math.floor((tau + p.T) * p.fs));
    for (let m = 0; m < p.np; m++) {
      const rcs = p2p ? drawRcs(tg.rcs, tg.model, rng) : scan;
      const amp = Math.sqrt((kSnr * rcs) / Math.pow(tg.range, 4));
      const phase = (-4 * Math.PI * tg.range) / d.lambda + 2 * Math.PI * fd * m * p.pri;
      const row = rx[m];
      for (let i = first; i <= last; i++) {
        const t = i / p.fs, u = t - tau - p.T / 2;
        if (u < -p.T / 2 || u >= p.T / 2) continue;
        const ph = Math.PI * k * u * u + phase + 2 * Math.PI * fd * t;
        row.re[i] += amp * Math.cos(ph);
        row.im[i] += amp * Math.sin(ph);
      }
    }
  }
  const s = Math.sqrt(0.5);
  for (const row of rx) for (let i = 0; i < d.ns; i++) { row.re[i] += s * rng.gauss(); row.im[i] += s * rng.gauss(); }
  return rx;
}

function pulseCompress(p, rx, window = true) {
  const d = derived(p), ref = lfmPulse(p);
  if (window) { const w = hamming(ref.re.length); for (let i = 0; i < w.length; i++) { ref.re[i] *= w[i]; ref.im[i] *= w[i]; } }
  const nfft = nextPow2(d.ns + ref.re.length);
  const hr = new Float64Array(nfft), hi = new Float64Array(nfft);
  hr.set(ref.re); hi.set(ref.im);
  fft(hr, hi);
  for (let i = 0; i < nfft; i++) hi[i] = -hi[i];  // conjugate -> matched filter
  return rx.map((row) => {
    const xr = new Float64Array(nfft), xi = new Float64Array(nfft);
    xr.set(row.re); xi.set(row.im);
    fft(xr, xi);
    for (let i = 0; i < nfft; i++) {
      const r = xr[i] * hr[i] - xi[i] * hi[i], im = xr[i] * hi[i] + xi[i] * hr[i];
      xr[i] = r; xi[i] = im;
    }
    fft(xr, xi, true);
    return { re: xr.slice(0, d.ns), im: xi.slice(0, d.ns) };
  });
}

// Power map rd[dopplerBin][rangeBin], zero velocity centred.
function rangeDoppler(pc) {
  const np = pc.length, ns = pc[0].re.length, nd = nextPow2(np), w = hamming(np);
  const rd = Array.from({ length: nd }, () => new Float64Array(ns));
  const cr = new Float64Array(nd), ci = new Float64Array(nd);
  for (let r = 0; r < ns; r++) {
    cr.fill(0); ci.fill(0);
    for (let m = 0; m < np; m++) { cr[m] = pc[m].re[r] * w[m]; ci[m] = pc[m].im[r] * w[m]; }
    fft(cr, ci);
    for (let k = 0; k < nd; k++) rd[(k + nd / 2) % nd][r] = cr[k] * cr[k] + ci[k] * ci[k];
  }
  return rd;
}

function cfar2d(p, rd, c) {
  const d = derived(p), nd = rd.length, nr = rd[0].length;
  const wr = c.gr + c.tr, wd = c.gd + c.td;
  const nTrain = (2 * wr + 1) * (2 * wd + 1) - (2 * c.gr + 1) * (2 * c.gd + 1);
  const alpha = nTrain * (Math.pow(c.pfa, -1 / nTrain) - 1);
  const vRes = d.lambda / (2 * nd * p.pri);
  const dets = [];
  for (let k = 0; k < nd; k++) {
    for (let r = wr; r < nr - wr; r++) {
      const cut = rd[k][r];
      let sum = 0;
      for (let dd = -wd; dd <= wd; dd++) {
        const row = rd[(k + dd + nd) % nd];
        for (let rr = -wr; rr <= wr; rr++) {
          if (Math.abs(dd) <= c.gd && Math.abs(rr) <= c.gr) continue;
          sum += row[r + rr];
        }
      }
      const noise = sum / nTrain;
      if (cut <= alpha * noise) continue;
      let peak = true;
      for (let dd = -1; dd <= 1 && peak; dd++)
        for (let rr = -1; rr <= 1; rr++)
          if ((dd || rr) && rd[(k + dd + nd) % nd][r + rr] > cut) { peak = false; break; }
      if (peak) dets.push({ k, r, range: (r * C) / (2 * p.fs), velocity: (k - nd / 2) * vRes, snrDb: 10 * Math.log10(cut / noise) });
    }
  }
  dets.sort((a, b) => b.snrDb - a.snrDb);
  return { dets, alpha };
}

window.Radar = { C, derived, lfmPulse, simulate, pulseCompress, rangeDoppler, cfar2d, fft, nextPow2 };
