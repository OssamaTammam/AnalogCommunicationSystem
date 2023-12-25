% Clear workspace and close all figures
clear;
close all;

% Experiment 3: Frequency Modulation and Demodulation

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

% Generate time vector
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

% Calculate constants for frequency modulation
kf = 0.2/(2*pi*max(abs(cumsum(y)))./Fs); % Frequency deviation
fc = 100000; % Carrier frequency
Ac = 1; % Amplitude

% Resample the signal at 5 times the carrier frequency
y = resample(y, 5*fc, Fs);
Fs = 5*fc;

% Generate time vector for the resampled signal
tstart = 0;
tend = tstart + length(y) / Fs;
t1 = linspace(tstart, tend, length(y));
t1 = t1';

% FM modulated signal using the frequency modulation equation
y = Ac * cos(2*pi*fc*t1 + 2*pi*kf*cumsum(y)./Fs);

% Fourier transform of the FM modulated signal
L = length(y);
F = fftshift(fft(y));
f = Fs/2 * linspace(-1, 1, L);

% Plot the FM modulated signal in the frequency domain
figure; 
plot(f, F); 
title('FM Modulation in Frequency Domain');

% Frequency Discriminator to calculate frequency deviation
d_y = diff(y);
d_y = [0; d_y];

% Envelope detector for demodulation
envelopeFM = abs(hilbert(d_y)) - mean(abs(hilbert(d_y)));

% Plot the demodulated signal in the time domain
figure; 
plot(t1, envelopeFM); 
title('Demodulated FM in Time Domain using Envelope Detector');
ylim([-2*10^-4 2*10^-4]);

% Resample the demodulated signal for sound playback
envelopeFM = resample(envelopeFM, Fs/5, Fs);
% Play the demodulated signal with increased volume
sound(500.*abs(envelopeFM), Fs/5);
