// Pulse-Doppler radar simulation: LFM waveform, Swerling targets, pulse
// compression (matched filter), Doppler processing and 2-D CA-CFAR detection.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include "fft.hpp"

namespace radar {

constexpr double kC = 299792458.0;  // speed of light (m/s)

struct RadarParams {
    double fc = 10e9;         // carrier frequency (Hz), X-band
    double bandwidth = 10e6;  // LFM sweep bandwidth B (Hz)
    double pulse_width = 10e-6;  // pulse duration T (s)
    double fs = 20e6;         // complex baseband sample rate (Hz)
    double pri = 100e-6;      // pulse repetition interval (s)
    int num_pulses = 64;      // pulses per coherent processing interval (CPI)
    double ref_snr_db = -10.0; // single-pulse SNR of a 1 m^2 target at ref_range (integration adds ~41 dB)
    double ref_range = 5000.0;

    double wavelength() const { return kC / fc; }
    int pulse_samples() const { return static_cast<int>(std::lround(pulse_width * fs)); }
    int pri_samples() const { return static_cast<int>(std::lround(pri * fs)); }
    double range_resolution() const { return kC / (2 * bandwidth); }
    double max_range() const { return kC * (pri - pulse_width) / 2; }
    double max_velocity() const { return wavelength() / (4 * pri); }  // unambiguous +/-
    double velocity_resolution() const { return wavelength() / (2 * num_pulses * pri); }
};

// Swerling fluctuation models for radar cross section (RCS).
enum class Swerling { Zero = 0, One = 1, Two = 2, Three = 3, Four = 4 };

struct Target {
    double range;     // m
    double velocity;  // m/s, positive = approaching
    double rcs;       // mean RCS (m^2)
    Swerling model = Swerling::Zero;
};

// Draw an RCS value. I/II: exponential (Rayleigh amplitude). III/IV: chi-square, 4 DOF.
inline double draw_rcs(double mean, Swerling m, std::mt19937_64& rng) {
    switch (m) {
        case Swerling::Zero: return mean;
        case Swerling::One:
        case Swerling::Two: return std::exponential_distribution<double>(1.0 / mean)(rng);
        case Swerling::Three:
        case Swerling::Four: {
            std::exponential_distribution<double> e(2.0 / mean);
            return e(rng) + e(rng);
        }
    }
    return mean;
}

// Baseband LFM (chirp) pulse: s(t) = exp(j*pi*(B/T)*t^2), t in [-T/2, T/2).
inline CVec lfm_pulse(const RadarParams& p) {
    const int n = p.pulse_samples();
    const double k = p.bandwidth / p.pulse_width;  // chirp rate (Hz/s)
    CVec s(n);
    for (int i = 0; i < n; ++i) {
        const double t = i / p.fs - p.pulse_width / 2;
        s[i] = std::polar(1.0, kPi * k * t * t);
    }
    return s;
}

// Fast-time x slow-time data cube (num_pulses rows of pri_samples each).
using Matrix = std::vector<CVec>;

// Simulate received echoes for all pulses in one CPI, with complex AWGN of unit power.
inline Matrix simulate_echoes(const RadarParams& p, const std::vector<Target>& targets,
                              std::uint64_t seed, bool add_noise = true) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> gauss(0.0, std::sqrt(0.5));
    const int ns = p.pri_samples(), np = p.num_pulses;
    const double lambda = p.wavelength();
    const double chirp_rate = p.bandwidth / p.pulse_width;
    const double k_snr = std::pow(10.0, p.ref_snr_db / 10) * std::pow(p.ref_range, 4);  // SNR = k*rcs/R^4

    Matrix rx(np, CVec(ns, cd(0)));
    for (const auto& tgt : targets) {
        const double tau = 2 * tgt.range / kC;
        const double fd = 2 * tgt.velocity / lambda;
        const double scan_rcs = draw_rcs(tgt.rcs, tgt.model, rng);  // Swerling I/III: fixed over CPI
        const bool pulse_to_pulse = tgt.model == Swerling::Two || tgt.model == Swerling::Four;
        for (int m = 0; m < np; ++m) {
            const double rcs = pulse_to_pulse ? draw_rcs(tgt.rcs, tgt.model, rng) : scan_rcs;
            const double amp = std::sqrt(k_snr * rcs / std::pow(tgt.range, 4));
            // Carrier phase from two-way range plus Doppler progression across pulses.
            const double phase = -4 * kPi * tgt.range / lambda + 2 * kPi * fd * m * p.pri;
            // Evaluate the delayed chirp s(t - tau) exactly, so sub-sample delays
            // don't distort the waveform (interpolation would create false echoes).
            const int first = std::max(0, static_cast<int>(std::ceil(tau * p.fs)));
            const int last = std::min(ns - 1, static_cast<int>(std::floor((tau + p.pulse_width) * p.fs)));
            for (int idx = first; idx <= last; ++idx) {
                const double t = idx / p.fs;
                const double u = t - tau - p.pulse_width / 2;  // time within pulse, centred
                if (u < -p.pulse_width / 2 || u >= p.pulse_width / 2) continue;
                rx[m][idx] += amp * std::polar(1.0, kPi * chirp_rate * u * u + phase + 2 * kPi * fd * t);
            }
        }
    }
    if (add_noise)
        for (auto& row : rx)
            for (auto& x : row) x += cd(gauss(rng), gauss(rng));
    return rx;
}

// Pulse compression: correlate each pulse with the (windowed) reference chirp via FFT.
// Output bin r corresponds to a round-trip delay of r samples.
inline Matrix pulse_compress(const RadarParams& p, const Matrix& rx, bool window = true) {
    const int ns = p.pri_samples();
    CVec ref = lfm_pulse(p);
    if (window) {
        const auto w = hamming(ref.size());
        for (std::size_t i = 0; i < ref.size(); ++i) ref[i] *= w[i];
    }
    const std::size_t nfft = next_pow2(ns + ref.size());
    CVec ref_f(nfft, cd(0));
    std::copy(ref.begin(), ref.end(), ref_f.begin());
    fft_inplace(ref_f);
    for (auto& x : ref_f) x = std::conj(x);  // matched filter H(f) = S*(f)

    Matrix out(rx.size(), CVec(ns));
    for (std::size_t m = 0; m < rx.size(); ++m) {
        CVec x(nfft, cd(0));
        std::copy(rx[m].begin(), rx[m].end(), x.begin());
        fft_inplace(x);
        for (std::size_t k = 0; k < nfft; ++k) x[k] *= ref_f[k];
        fft_inplace(x, true);
        std::copy(x.begin(), x.begin() + ns, out[m].begin());
    }
    return out;
}

// Doppler FFT across pulses for each range bin. Returns a power map
// rd[doppler_bin][range_bin] with zero velocity centred (fftshift).
inline std::vector<std::vector<double>> range_doppler_map(const Matrix& compressed, bool window = true) {
    const std::size_t np = compressed.size(), ns = compressed.at(0).size();
    const std::size_t nd = next_pow2(np);
    const auto w = window ? hamming(np) : std::vector<double>(np, 1.0);
    std::vector<std::vector<double>> rd(nd, std::vector<double>(ns));
    CVec col(nd);
    for (std::size_t r = 0; r < ns; ++r) {
        std::fill(col.begin(), col.end(), cd(0));
        for (std::size_t m = 0; m < np; ++m) col[m] = compressed[m][r] * w[m];
        fft_inplace(col);
        for (std::size_t d = 0; d < nd; ++d) rd[(d + nd / 2) % nd][r] = std::norm(col[d]);
    }
    return rd;
}

struct CfarParams {
    int guard_range = 2, guard_doppler = 1;
    int train_range = 8, train_doppler = 4;
    double pfa = 1e-6;
};

struct Detection {
    int doppler_bin, range_bin;
    double power_db, snr_db;  // snr_db = cell power over local noise estimate
    double range_m, velocity_mps;
};

// Cell-averaging CFAR threshold factor for N training cells and target Pfa.
inline double ca_cfar_alpha(int n_train, double pfa) {
    return n_train * (std::pow(pfa, -1.0 / n_train) - 1.0);
}

// 2-D CA-CFAR over the range-Doppler power map, keeping only local maxima.
inline std::vector<Detection> cfar_2d(const RadarParams& p, const std::vector<std::vector<double>>& rd,
                                      const CfarParams& c = {}) {
    const int nd = static_cast<int>(rd.size()), nr = static_cast<int>(rd[0].size());
    const int wr = c.guard_range + c.train_range, wd = c.guard_doppler + c.train_doppler;
    const int n_train = (2 * wr + 1) * (2 * wd + 1) - (2 * c.guard_range + 1) * (2 * c.guard_doppler + 1);
    const double alpha = ca_cfar_alpha(n_train, c.pfa);
    std::vector<Detection> dets;
    const double vres = p.wavelength() / (2 * nd * p.pri);
    for (int d = 0; d < nd; ++d) {
        for (int r = wr; r < nr - wr; ++r) {
            double sum = 0;
            for (int dd = -wd; dd <= wd; ++dd)
                for (int rr = -wr; rr <= wr; ++rr) {
                    if (std::abs(dd) <= c.guard_doppler && std::abs(rr) <= c.guard_range) continue;
                    sum += rd[(d + dd + nd) % nd][r + rr];  // Doppler wraps around
                }
            const double noise = sum / n_train;
            const double cut = rd[d][r];
            if (cut <= alpha * noise) continue;
            bool peak = true;  // suppress neighbours of the same target
            for (int dd = -1; dd <= 1 && peak; ++dd)
                for (int rr = -1; rr <= 1; ++rr)
                    if ((dd || rr) && rd[(d + dd + nd) % nd][r + rr] > cut) { peak = false; break; }
            if (!peak) continue;
            dets.push_back({d, r, pow_db(cut), pow_db(cut / noise), r * kC / (2 * p.fs),
                            (d - nd / 2) * vres});
        }
    }
    std::sort(dets.begin(), dets.end(), [](const Detection& a, const Detection& b) { return a.power_db > b.power_db; });
    return dets;
}

}  // namespace radar
