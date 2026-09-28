// Radix-2 FFT and small DSP helpers used by the radar pipeline.
#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace radar {

using cd = std::complex<double>;
using CVec = std::vector<cd>;
constexpr double kPi = 3.14159265358979323846;

inline bool is_pow2(std::size_t n) { return n && !(n & (n - 1)); }

inline std::size_t next_pow2(std::size_t n) {
    std::size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

// In-place iterative Cooley-Tukey FFT. inverse=true computes the IFFT (scaled by 1/N).
inline void fft_inplace(CVec& a, bool inverse = false) {
    const std::size_t n = a.size();
    if (!is_pow2(n)) throw std::invalid_argument("FFT size must be a power of two");

    // Bit-reversal permutation.
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    // Butterflies.
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double ang = 2 * kPi / static_cast<double>(len) * (inverse ? 1 : -1);
        const cd wlen(std::cos(ang), std::sin(ang));
        for (std::size_t i = 0; i < n; i += len) {
            cd w(1);
            for (std::size_t k = 0; k < len / 2; ++k) {
                const cd u = a[i + k];
                const cd v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse)
        for (auto& x : a) x /= static_cast<double>(n);
}

inline CVec fft(CVec a) { fft_inplace(a, false); return a; }
inline CVec ifft(CVec a) { fft_inplace(a, true); return a; }

// Hamming window of length n (used to lower range and Doppler sidelobes).
inline std::vector<double> hamming(std::size_t n) {
    std::vector<double> w(n, 1.0);
    if (n < 2) return w;
    for (std::size_t i = 0; i < n; ++i)
        w[i] = 0.54 - 0.46 * std::cos(2 * kPi * static_cast<double>(i) / static_cast<double>(n - 1));
    return w;
}

// Swap halves so zero frequency sits in the middle (like MATLAB fftshift).
template <typename T>
std::vector<T> fftshift(const std::vector<T>& x) {
    const std::size_t n = x.size(), h = n / 2;
    std::vector<T> y(n);
    for (std::size_t i = 0; i < n; ++i) y[(i + h) % n] = x[i];
    return y;
}

inline double pow_db(double p) { return 10.0 * std::log10(p + 1e-30); }

}  // namespace radar
