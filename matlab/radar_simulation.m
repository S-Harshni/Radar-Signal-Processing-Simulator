% RADAR_SIMULATION  Pulse-Doppler radar chain in base MATLAB (no toolboxes).
%   LFM waveform -> target echoes (Swerling RCS) -> pulse compression
%   (matched filter) -> Doppler FFT -> 2-D CA-CFAR detection.
%   Mirrors cpp/include/radar/radar.hpp so results can be cross-checked.
%   Requires MATLAB R2016b+ (local functions in scripts).

clear; close all; rng(3);

%% Radar parameters
p.c   = 299792458;       % speed of light (m/s)
p.fc  = 10e9;            % carrier (Hz)
p.B   = 10e6;            % chirp bandwidth (Hz)
p.T   = 10e-6;           % pulse width (s)
p.fs  = 20e6;            % complex sample rate (Hz)
p.PRI = 100e-6;          % pulse repetition interval (s)
p.Np  = 64;              % pulses per CPI
p.refSNRdB = -10;        % single-pulse SNR of 1 m^2 at refRange
p.refRange = 5000;
lambda = p.c / p.fc;
Ns = round(p.PRI * p.fs);          % fast-time samples per PRI
fprintf('Range res %.1f m, velocity res %.2f m/s, max velocity +/-%.1f m/s\n', ...
    p.c/(2*p.B), lambda/(2*p.Np*p.PRI), lambda/(4*p.PRI));

%% Targets: [range (m), velocity (m/s), mean RCS (m^2), Swerling model]
targets = [2500 -30  2 1;
           6000  12 10 0;
           9000  45 20 3];

%% Transmit waveform (baseband LFM)
n   = round(p.T * p.fs);
t   = (0:n-1).'/p.fs - p.T/2;
k   = p.B / p.T;
ref = exp(1j*pi*k*t.^2);

%% Simulate echoes (fast time x slow time)
rx = zeros(Ns, p.Np);
kSNR = 10^(p.refSNRdB/10) * p.refRange^4;
tFast = (0:Ns-1).'/p.fs;
for i = 1:size(targets, 1)
    R = targets(i,1); v = targets(i,2); sigma = targets(i,3); model = targets(i,4);
    tau = 2*R/p.c;  fd = 2*v/lambda;
    scanRCS = drawRCS(sigma, model);
    u = tFast - tau - p.T/2;                     % time within the echo pulse
    inPulse = (u >= -p.T/2) & (u < p.T/2);
    for m = 1:p.Np
        if model == 2 || model == 4, rcs = drawRCS(sigma, model); else, rcs = scanRCS; end
        amp = sqrt(kSNR * rcs / R^4);
        phase = -4*pi*R/lambda + 2*pi*fd*(m-1)*p.PRI;
        rx(:,m) = rx(:,m) + inPulse .* amp .* exp(1j*(pi*k*u.^2 + phase + 2*pi*fd*tFast));
    end
end
rx = rx + (randn(Ns, p.Np) + 1j*randn(Ns, p.Np)) / sqrt(2);   % unit-power AWGN

%% Pulse compression (FFT-based matched filter, Hamming-weighted reference)
nfft = 2^nextpow2(Ns + n);
H = conj(fft(ref .* hammingWin(n), nfft));
pc = ifft(fft(rx, nfft) .* H);
pc = pc(1:Ns, :);

%% Doppler processing
Nd = 2^nextpow2(p.Np);
rd = fftshift(fft(pc .* hammingWin(p.Np).', Nd, 2), 2);
rdPow = abs(rd).^2;
rangeAxis = (0:Ns-1) * p.c / (2*p.fs);
velAxis   = (-Nd/2:Nd/2-1) * lambda / (2*Nd*p.PRI);

%% 2-D CA-CFAR
[det, alpha] = caCfar2d(rdPow, 2, 1, 8, 4, 1e-6);
[ri, di] = find(det);
fprintf('\nCFAR threshold factor alpha = %.1f (%.1f dB)\n', alpha, 10*log10(alpha));
fprintf('Detections:\n');
for q = 1:numel(ri)
    fprintf('  %8.1f m  %7.2f m/s\n', rangeAxis(ri(q)), velAxis(di(q)));
end

%% Plots
figure('Name', 'Radar signal processing');
subplot(2,2,1); plot(t*1e6, real(ref)); grid on;
xlabel('Time (\mus)'); ylabel('Re\{s(t)\}'); title('Transmitted LFM pulse');
subplot(2,2,2); plot(rangeAxis/1e3, 20*log10(abs(rx(:,1)) + eps), rangeAxis/1e3, 20*log10(abs(pc(:,1)) + eps));
grid on; xlabel('Range (km)'); ylabel('dB'); legend('Raw echo', 'Pulse compressed'); title('Single pulse');
subplot(2,1,2); imagesc(rangeAxis/1e3, velAxis, 10*log10(rdPow.' + eps)); axis xy; colorbar;
hold on; plot(rangeAxis(ri)/1e3, velAxis(di), 'wo', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('Range (km)'); ylabel('Velocity (m/s)'); title('Range-Doppler map with CFAR detections');

%% Local functions
function w = hammingWin(n)
    w = 0.54 - 0.46*cos(2*pi*(0:n-1).'/(n-1));
end

function s = drawRCS(meanRCS, model)
    % Swerling 0: constant. I/II: exponential. III/IV: chi-square with 4 DOF.
    switch model
        case 0,      s = meanRCS;
        case {1, 2}, s = -meanRCS * log(rand);
        otherwise,   s = -(meanRCS/2) * (log(rand) + log(rand));
    end
end

function [det, alpha] = caCfar2d(P, gr, gd, tr, td, pfa)
    % Cell-averaging CFAR on a range x Doppler power map (Doppler wraps).
    [Nr, Nd] = size(P);
    wr = gr + tr; wd = gd + td;
    kernel = ones(2*wr+1, 2*wd+1);
    kernel(tr+1:tr+2*gr+1, td+1:td+2*gd+1) = 0;   % remove guard cells + CUT
    N = nnz(kernel);
    alpha = N * (pfa^(-1/N) - 1);
    Pwrap = [P(:, end-wd+1:end), P, P(:, 1:wd)];   % circular in Doppler
    noise = conv2(Pwrap, kernel, 'same') / N;
    noise = noise(:, wd+1:wd+Nd);
    det = P > alpha * noise;
    det([1:wr, Nr-wr+1:Nr], :) = false;           % ignore range edges
    % Keep local maxima only.
    localMax = P >= max(cat(3, circshift(P,[1 0]), circshift(P,[-1 0]), circshift(P,[0 1]), ...
        circshift(P,[0 -1]), circshift(P,[1 1]), circshift(P,[-1 -1]), circshift(P,[1 -1]), circshift(P,[-1 1])), [], 3);
    det = det & localMax;
end
