% Clear workspace and close all figures
clear;
close all;

% Experiment 1: DSBSC, DSBTC, and Coherent Detection

% Read the audio signal
[S, Fs] = audioread('eric.wav');

% Fourier transform of the original signal
L = length(S);
F = fftshift(fft(S));
f = Fs/2 * linspace(-1, 1, L);

% Plot the original signal spectrum
figure; 
plot(f, abs(F)/L); 
title('Original Signal Spectrum');

% Filter the signal to retain frequencies around 4 kHz
W = 4000;
F(f >= W | f <= -W) = 0;
y = ifft(ifftshift(F));

% Fourier transform of the filtered signal
L = length(y);
F = fftshift(fft(y));
f = Fs/2 * linspace(-1, 1, L);

% Plot the filtered signal spectrum
figure; 
plot(f, abs(F)/L); 
title('Filtered Signal Spectrum');

% Generate time vector for the filtered signal
tstart = 0;
tend = tstart + length(y) / Fs;
t1 = linspace(tstart, tend, length(y));
t1 = t1';

% Plot the filtered signal in the time domain
figure; 
plot(t1, y); 
title('Filtered Signal Time Domain');

% Play the filtered signal
sound(abs(y), Fs);

% Calculate constants for DSBSC and DSBTC modulation
fm = W; % Modulating frequency
fc = 100000; % Carrier frequency
mu = 0.5; % Modulation index
Am = max(y); % Maximum amplitude of the modulating signal
Ac = Am/mu; % Amplitude of the carrier signal

% Resample the signal at 5 times the carrier frequency
y = resample(y, 5*fc, Fs);
Fs = 5*fc;

% DSBSC (Double Sideband Suppressed Carrier) Modulation

% Generate time vector for the resampled signal
tstart = 0;
tend = tstart + length(y) / Fs;
t1 = linspace(tstart, tend, length(y));
t1 = t1';

% DSBSC generation
carrier_signal = Ac .* cos(2*pi*fc*t1);
DSBSC = y .* carrier_signal;

% Fourier transform of DSBSC
L = length(DSBSC);
F = fftshift(fft(DSBSC));
f = Fs/2 * linspace(-1, 1, L);

% Plot DSBSC in the frequency domain
figure; 
plot(f, abs(F) / L); 
title('DSBSC Frequency Domain');

% DSBTC (Double Sideband with Carrier) Modulation

DSBTC = (1 + mu * y / Am) .* carrier_signal;

% Fourier transform of DSBTC
L = length(DSBTC);
F = fftshift(fft(DSBTC));
f = Fs/2 * linspace(-1, 1, L);

% Plot DSBTC in the frequency domain
figure; 
plot(f, abs(F) / L); 
title('DSBTC Frequency Domain');

% Envelope detector for DSBSC modulation
envelopeDSBSC = abs(hilbert(DSBSC));

% Plot DSBSC with envelope detector in the time domain
figure; 
plot(t1, DSBSC);
hold on;
plot(t1, -envelopeDSBSC, 'r-', t1, envelopeDSBSC, '-r', 'LineWidth', 1.5);
hold off;
title('DSBSC Time Domain with Envelope Detector');
ylim([-2 2]);
xlim([2 2.5]);

% Resample the envelope signal for sound playback
envelopeDSBSC = resample(abs(envelopeDSBSC), Fs/5, Fs);
sound(abs(envelopeDSBSC), Fs/5);

% Envelope detector for DSBTC modulation
envelopeDSBTC = abs(hilbert(DSBTC));

% Plot DSBTC with envelope detector in the time domain
figure; 
plot(t1, DSBTC);
hold on;
plot(t1, -envelopeDSBTC, 'r-', t1, envelopeDSBTC, '-r', 'LineWidth', 1.5);
hold off;
title('DSBTC Time Domain with Envelope Detector');
ylim([-5 5]);
xlim([3 3.5]);

% Resample the envelope signal for sound playback
envelopeDSBTC = resample(envelopeDSBTC, Fs/5, Fs);
sound(abs(envelopeDSBTC), Fs/5);

% Coherent Detection with various SNR and frequency/phase errors

% Initialize carrier frequency for coherent detection
fc = 100000;

% Loop through different SNR values
for snr = [0, 10, 30]
    % Generate signal with noise
    noisy_DSBSC = awgn(DSBSC, snr);
    
    % Demodulate using coherent detector
    demodulated = noisy_DSBSC .* cos(2*pi*fc*t1);
    
    % Fourier transform of the demodulated signal
    demodulated_FFT = fftshift(fft(demodulated));
    
    % LPF at modulation frequency
    demodulated_FFT(f >= W | f <= -W) = 0;
    
    % Inverse Fourier transform to get demodulated signal in time domain
    demodulated = ifft(ifftshift(demodulated_FFT));
    
    % Plot demodulated signal in time domain
    figure; 
    plot(t1, demodulated); 
    title(sprintf('%d SNR Demodulated Signal in Time Domain', snr));
    
    % Fourier transform of the demodulated signal
    L = length(demodulated);
    F = fftshift(fft(demodulated));
    f = Fs/2 * linspace(-1, 1, L);
    
    % Plot demodulated signal in frequency domain
    figure; 
    plot(f, abs(F) / L); 
    title(sprintf('%d SNR Demodulated Signal in Frequency Domain', snr));
    
    % Resample to sound the demodulated signal
    demodulated = resample(abs(demodulated), Fs/5, Fs);
    sound(abs(demodulated), Fs/5);
end

% Coherent Detection with frequency error

% Initialize carrier frequency with an error
fc = 100100;

% Loop through different SNR values
for snr = [0, 10, 30]
    % Generate signal with noise
    noisy_DSBSC = awgn(DSBSC, snr);
    
    % Demodulate using coherent detector
    demodulated = noisy_DSBSC .* cos(2*pi*fc*t1);
    
    % Fourier transform of the demodulated signal
    demodulated_FFT = fftshift(fft(demodulated));
    
    % LPF at modulation frequency
    demodulated_FFT(f >= W | f <= -W) = 0;
    
    % Inverse Fourier transform to get demodulated signal in time domain
    demodulated = ifft(ifftshift(demodulated_FFT));
    
    % Plot demodulated signal in time domain
    figure; 
    plot(t1, demodulated); 
    title(sprintf('%d SNR Demodulated Signal with Frequency Error in Time Domain', snr));
    
    % Fourier transform of the demodulated signal
    L = length(demodulated);
    F = fftshift(fft(demodulated));
    f = Fs/2 * linspace(-1, 1, L);
    
    % Plot demodulated signal in frequency domain
    figure; 
    plot(f, abs(F) / L); 
    title(sprintf('%d SNR Demodulated Signal with Frequency Error in Frequency Domain', snr));
    
    % Resample to sound the demodulated signal
    demodulated = resample(abs(demodulated), Fs/5, Fs);
    sound(abs(demodulated), Fs/5);
end

% Coherent Detection with phase error

% Initialize carrier frequency and phase error
fc = 100000;
phase_error = pi/9;

% Loop through different SNR values
for snr = [0, 10, 30]
    % Generate signal with noise
    noisy_DSBSC = awgn(DSBSC, snr);
    
    % Demodulate using coherent detector with phase error
    demodulated = noisy_DSBSC .* cos(2*pi*fc*t1 + phase_error);
    
    % Fourier transform of the demodulated signal
    demodulated_FFT = fftshift(fft(demodulated));
    
    % LPF at modulation frequency
    demodulated_FFT(f >= W | f <= -W) = 0;
    
    % Inverse Fourier transform to get demodulated signal in time domain
    demodulated = ifft(ifftshift(demodulated_FFT));
    
    % Plot demodulated signal in time domain
    figure; 
    plot(t1, demodulated); 
    title(sprintf('%d SNR Demodulated Signal with Phase Error in Time Domain', snr));
    
    % Fourier transform of the demodulated signal
    L = length(demodulated);
    F = fftshift(fft(demodulated));
    f = Fs/2 * linspace(-1, 1, L);
    
    % Plot demodulated signal in frequency domain
    figure; 
    plot(f, abs(F) / L); 
    title(sprintf('%d SNR Demodulated Signal with Phase Error in Frequency Domain', snr));
    
    % Resample to sound the demodulated signal
    demodulated = resample(abs(demodulated), Fs/5, Fs);
    sound(abs(demodulated), Fs/5);
end
