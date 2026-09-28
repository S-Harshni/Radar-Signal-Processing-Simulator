// Command-line radar scenario: simulate a CPI, process it and print detections.
// Usage: radar_sim [--seed N] [--csv rd_map.csv] [--target R V RCS SWERLING]...
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#include "radar/radar.hpp"

using namespace radar;

int main(int argc, char** argv) {
    RadarParams p;
    std::vector<Target> targets;
    std::uint64_t seed = 1;
    std::string csv;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--seed") && i + 1 < argc) seed = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--csv") && i + 1 < argc) csv = argv[++i];
        else if (!std::strcmp(argv[i], "--target") && i + 4 < argc) {
            Target t{std::atof(argv[i + 1]), std::atof(argv[i + 2]), std::atof(argv[i + 3]),
                     static_cast<Swerling>(std::atoi(argv[i + 4]))};
            targets.push_back(t);
            i += 4;
        } else {
            std::fprintf(stderr, "usage: %s [--seed N] [--csv file] [--target R V RCS SWERLING]...\n", argv[0]);
            return 2;
        }
    }
    if (targets.empty())
        targets = {{2500, -30, 2.0, Swerling::One}, {6000, 12, 10.0, Swerling::Zero}, {9000, 45, 20.0, Swerling::Three}};

    std::printf("Radar: fc %.1f GHz, B %.0f MHz, T %.0f us, PRI %.0f us, %d pulses\n", p.fc / 1e9,
                p.bandwidth / 1e6, p.pulse_width * 1e6, p.pri * 1e6, p.num_pulses);
    std::printf("Range res %.1f m, max range %.1f km | velocity res %.2f m/s, max +/-%.1f m/s\n\n",
                p.range_resolution(), p.max_range() / 1e3, p.velocity_resolution(), p.max_velocity());

    const auto rd = range_doppler_map(pulse_compress(p, simulate_echoes(p, targets, seed)));
    const auto dets = cfar_2d(p, rd);

    std::puts("Truth:");
    for (const auto& t : targets)
        std::printf("  %8.1f m  %7.2f m/s  RCS %5.1f m^2  Swerling %d\n", t.range, t.velocity, t.rcs,
                    static_cast<int>(t.model));
    std::printf("\nCFAR detections (Pfa 1e-6): %zu\n", dets.size());
    for (const auto& d : dets)
        std::printf("  %8.1f m  %7.2f m/s  SNR %5.1f dB\n", d.range_m, d.velocity_mps, d.snr_db);

    if (!csv.empty()) {
        std::ofstream f(csv);
        for (const auto& row : rd) {
            for (std::size_t r = 0; r < row.size(); ++r) f << (r ? "," : "") << pow_db(row[r]);
            f << '\n';
        }
        std::printf("\nWrote range-Doppler map (dB) to %s\n", csv.c_str());
    }
    return 0;
}
