// Unit tests: each check compares the implementation against theory.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>

#include "radar/radar.hpp"

using namespace radar;

static int failures = 0;
#define CHECK(cond, msg)                                              \
    do {                                                              \
        if (cond) std::printf("  PASS  %s\n", msg);                   \
        else { std::printf("  FAIL  %s  (%s:%d)\n", msg, __FILE__, __LINE__); ++failures; } \
    } while (0)

static void test_fft() {
    std::puts("FFT");
    const int n = 256, k0 = 19;
    CVec x(n);
    for (int i = 0; i < n; ++i) x[i] = std::polar(1.0, 2 * kPi * k0 * i / n);
    const CVec X = fft(x);
    int peak = 0;
    for (int i = 0; i < n; ++i) if (std::abs(X[i]) > std::abs(X[peak])) peak = i;
    CHECK(peak == k0, "complex tone lands in the expected bin");
    CHECK(std::abs(std::abs(X[k0]) - n) < 1e-9, "tone bin magnitude equals N");

    std::mt19937 rng(1);
    std::normal_distribution<double> g;
    CVec y(n);
    for (auto& v : y) v = cd(g(rng), g(rng));
    const CVec back = ifft(fft(y));
    double err = 0, ey = 0, eY = 0;
    for (int i = 0; i < n; ++i) err = std::max(err, std::abs(back[i] - y[i]));
    const CVec Y = fft(y);
    for (int i = 0; i < n; ++i) { ey += std::norm(y[i]); eY += std::norm(Y[i]); }
    CHECK(err < 1e-12, "IFFT(FFT(x)) reconstructs x");
    CHECK(std::abs(ey - eY / n) / ey < 1e-12, "Parseval: energy is preserved");
}

static void test_waveform() {
    std::puts("LFM waveform");
    RadarParams p;
    const CVec s = lfm_pulse(p);
    CHECK(static_cast<int>(s.size()) == 200, "pulse has T*fs = 200 samples");
    bool unit = true;
    for (auto v : s) unit &= std::abs(std::abs(v) - 1) < 1e-12;
    CHECK(unit, "constant envelope (|s| = 1)");
    // Instantaneous frequency sweeps -B/2..+B/2 linearly.
    const double f_start = std::arg(s[1] * std::conj(s[0])) * p.fs / (2 * kPi);
    const double f_end = std::arg(s[199] * std::conj(s[198])) * p.fs / (2 * kPi);
    CHECK(std::abs(f_start + p.bandwidth / 2) < 0.02 * p.bandwidth, "chirp starts at -B/2");
    CHECK(std::abs(f_end - p.bandwidth / 2) < 0.02 * p.bandwidth, "chirp ends at +B/2");
    CHECK(std::abs(p.range_resolution() - 14.99) < 0.01, "range resolution c/2B = 15 m");
}

static void test_pulse_compression() {
    std::puts("Pulse compression");
    RadarParams p;
    p.num_pulses = 1;
    const double R = 3000.0;
    const auto rx = simulate_echoes(p, {{R, 0.0, 1.0}}, 7, false);
    const auto pc = pulse_compress(p, rx, false);
    int peak = 0;
    for (int i = 0; i < p.pri_samples(); ++i) if (std::abs(pc[0][i]) > std::abs(pc[0][peak])) peak = i;
    const double expected = 2 * R / kC * p.fs;
    CHECK(std::abs(peak - expected) <= 1.0, "compressed peak at round-trip delay (3 km)");
    // Unwindowed matched filter: peak amplitude = amp * N (processing gain N = B*T*fs/B).
    const double amp = std::abs(rx[0][static_cast<int>(expected) + 50]);
    const double gain = std::abs(pc[0][peak]) / amp;
    CHECK(gain > 0.9 * p.pulse_samples(), "processing gain ~= pulse length (200)");
    // Mainlobe (-3 dB) width ~ fs/B = 2 samples.
    const double half = std::abs(pc[0][peak]) / std::sqrt(2.0);
    int w = 0;
    for (int i = peak - 10; i <= peak + 10; ++i) w += std::abs(pc[0][i]) >= half;
    CHECK(w >= 1 && w <= 3, "mainlobe width ~= 1/B (compressed from 10 us)");
}

static void test_range_doppler() {
    std::puts("Range-Doppler processing");
    RadarParams p;
    const Target t{4200.0, 25.0, 5.0};
    const auto rd = range_doppler_map(pulse_compress(p, simulate_echoes(p, {t}, 11)));
    int bd = 0, br = 0;
    for (std::size_t d = 0; d < rd.size(); ++d)
        for (std::size_t r = 0; r < rd[d].size(); ++r)
            if (rd[d][r] > rd[bd][br]) { bd = static_cast<int>(d); br = static_cast<int>(r); }
    const double r_m = br * kC / (2 * p.fs);
    const double v = (bd - static_cast<int>(rd.size()) / 2) * p.velocity_resolution();
    std::printf("        peak: %.1f m, %.2f m/s (truth 4200 m, 25 m/s)\n", r_m, v);
    CHECK(std::abs(r_m - t.range) <= p.range_resolution(), "range estimate within one resolution cell");
    CHECK(std::abs(v - t.velocity) <= p.velocity_resolution(), "velocity estimate within one Doppler bin");
}

static void test_cfar() {
    std::puts("CA-CFAR");
    RadarParams p;
    CfarParams c;
    c.pfa = 1e-3;
    // Noise only: count raw threshold crossings (without peak picking) against Pfa.
    const auto rd = range_doppler_map(pulse_compress(p, simulate_echoes(p, {}, 3)), false);
    const int nd = static_cast<int>(rd.size()), nr = static_cast<int>(rd[0].size());
    const int wr = c.guard_range + c.train_range, wd = c.guard_doppler + c.train_doppler;
    const int n_train = (2 * wr + 1) * (2 * wd + 1) - (2 * c.guard_range + 1) * (2 * c.guard_doppler + 1);
    const double alpha = ca_cfar_alpha(n_train, c.pfa);
    long cells = 0, crossings = 0;
    for (int d = 0; d < nd; ++d)
        for (int r = wr + p.pulse_samples(); r < nr - wr; r += 3) {
            double sum = 0;
            for (int dd = -wd; dd <= wd; ++dd)
                for (int rr = -wr; rr <= wr; ++rr)
                    if (std::abs(dd) > c.guard_doppler || std::abs(rr) > c.guard_range)
                        sum += rd[(d + dd + nd) % nd][r + rr];
            ++cells;
            crossings += rd[d][r] > alpha * sum / n_train;
        }
    const double pfa_hat = static_cast<double>(crossings) / cells;
    std::printf("        measured Pfa = %.2e over %ld cells (design 1e-3)\n", pfa_hat, cells);
    CHECK(pfa_hat < 3e-3, "false-alarm rate stays near the design Pfa");

    const std::vector<Target> targets{{2500, -30, 2.0}, {6000, 12, 10.0}};
    const auto dets = cfar_2d(p, range_doppler_map(pulse_compress(p, simulate_echoes(p, targets, 5))));
    int found = 0;
    for (const auto& t : targets)
        for (const auto& d : dets)
            if (std::abs(d.range_m - t.range) < 2 * p.range_resolution() &&
                std::abs(d.velocity_mps - t.velocity) < 2 * p.velocity_resolution()) { ++found; break; }
    CHECK(found == 2, "both targets detected at the right range and velocity");
}

static void test_swerling() {
    std::puts("Swerling RCS models");
    std::mt19937_64 rng(42);
    for (auto m : {Swerling::One, Swerling::Three}) {
        double s = 0, s2 = 0;
        const int n = 200000;
        for (int i = 0; i < n; ++i) { const double x = draw_rcs(4.0, m, rng); s += x; s2 += x * x; }
        const double mean = s / n, var = s2 / n - mean * mean;
        // Exponential: var = mean^2. Chi-square 4 DOF: var = mean^2 / 2.
        const double expected_var = m == Swerling::One ? 16.0 : 8.0;
        CHECK(std::abs(mean - 4.0) < 0.05, m == Swerling::One ? "Swerling I/II mean RCS" : "Swerling III/IV mean RCS");
        CHECK(std::abs(var - expected_var) / expected_var < 0.05,
              m == Swerling::One ? "Swerling I/II variance = mean^2" : "Swerling III/IV variance = mean^2/2");
    }
}

int main() {
    test_fft();
    test_waveform();
    test_pulse_compression();
    test_range_doppler();
    test_cfar();
    test_swerling();
    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL TESTS PASSED", failures, failures == 1 ? "" : "s");
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
