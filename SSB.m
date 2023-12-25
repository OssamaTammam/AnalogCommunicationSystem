% Clear existing variables and close all figures
clear;
close all;

% Experiment 2

% Read the attached audio file
[S, Fs] = audioread('eric.wav');

% Fourier transform of the original signal
L = length(S);
F = fftshift(fft(S));
f = Fs/2 * linspace(-1, 1, L);

% Plot original signal spectrum
figure; 
plot(f, abs(F)/L); 
title('Original Signal Spectrum');

% Filter the signal to 4kHz
W = 4000;
F(f >= W | f <= -W) = 0;
y = ifft(ifftshift(F));

% Plot filtered signal spectrum
L = length(y);
F = fftshift(fft(y));
f = Fs/2 * linspace(-1, 1, L);
figure; 
plot(f, abs(F)/L); 
title('Filtered Signal Spectrum');

% Time vector calculation
tstart = 0;
tend = tstart + length(y) / Fs;
t1 = linspace(tstart, tend, length(y));
t1 = t1';

% Plot filtered signal in time domain
figure; 
plot(t1, y); 
title('Filtered Signal Time Domain');

% Listen to the filtered signal
sound(abs(y), Fs);

% Calculate constants for modulation
fm = W; % Modulating frequency
fc = 100000; % Carrier frequency
Am = max(y); % Maximum amplitude of the modulating signal
Ac = 2 * Am; % Amplitude of the carrier signal

% Resample at 5 times the carrier frequency
y = resample(y, 5 * fc, Fs);
Fs = 5 * fc;

% SSBSC generation

% Calculate time vector
tstart = 0;
tend = tstart + length(y) / Fs;
t1 = linspace(tstart, tend, length(y));
t1 = t1';

% DSBSC generation
carrier_signal = Ac .* cos(2*pi*fc*t1);
DSBSC = y .* carrier_signal;

% SSBLSB generation by filtering DSBSC
SSBLSB = DSBSC;
L = length(SSBLSB);
f = Fs/2 * linspace(-1, 1, L);
F = fftshift(fft(SSBLSB));

% Filter (SSB in frequency domain)
F(f >= fc | f <= -fc) = 0;
SSBLSB = ifft(ifftshift(F));

% Plot SSBLSB frequency domain
figure; 
plot(f, abs(F) / L); 
title('SSBLSB Frequency Domain');

% Coherent Detection SSBSC using an ideal filter

% Demodulate using coherent detector
demodulated = SSBLSB .* cos(2*pi*fc*t1);

% Fourier transform of the demodulated signal
demodulated_FFT = fftshift(fft(demodulated));

% LPF at modulation frequency
demodulated_FFT(f >= W | f <= -W) = 0;

% Inverse Fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% Plot demodulated signal in time domain
figure; 
plot(t1, demodulated); 
title('Demodulated Signal using ideal filter in Time Domain');

% Fourier transform of the demodulated signal
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2 * linspace(-1, 1, L);

% Plot demodulated signal in frequency domain
figure; 
plot(f, abs(F) / L); 
title('Demodulated Signal using ideal filter in Frequency Domain');

% Resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);

% Coherent Detection SSBSC using Butterworth filter (order = 4)

% Demodulate using coherent detector
demodulated = SSBLSB .* cos(2*pi*fc*t1);

% Design a Butterworth filter
[b, a] = butter(4, W * 2 / Fs);

% Apply the filter to the demodulated signal
demodulated = filtfilt(b, a, demodulated);

% Plot demodulated signal in time domain
figure; 
plot(t1, demodulated); 
title('Demodulated Signal using Butterworth filter in Time Domain');

% Fourier transform of the demodulated signal
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2 * linspace(-1, 1, L);

% Plot demodulated signal in frequency domain
figure; 
plot(f, abs(F) / L); 
title('Demodulated Signal using Butterworth filter in Frequency Domain');

% Resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);

% Coherent Detection SSBSC with different SNR values

for snr = [0, 10, 30]
    % Generate signal with noise
    noisy_SSBLSB = awgn(SSBLSB, snr);
    
    % Demodulate using coherent detector
    demodulated = noisy_SSBLSB .* cos(2*pi*fc*t1);
    
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

% SSBTC generation
SSBTC = carrier_signal + SSBLSB;

% Fourier transform of SSBTC
L = length(SSBTC);
F = fftshift(fft(SSBTC));
f = Fs/2 * linspace(-1, 1, L);

% Plot SSBTC frequency domain
figure; 
plot(f, abs(F) / L); 
title('SSBTC in Frequency Domain');

% Envelope detector for SSBTC
envelopeSSBTC = abs(hilbert(SSBTC));

% Plot SSBTC time domain with envelope detector
figure; 
plot(t1, SSBTC);
hold on;
plot(t1, -envelopeSSBTC, 'r-', t1, envelopeSSBTC, '-r', 'Linewidth', 1.5);
title('SSBTC Time Domain with Envelope Detector (Red)');
hold off;
ylim([-5 5])
xlim([3 3.5])

% Resample to sound the signal
envelopeSSBTC = resample(envelopeSSBTC, Fs/5, Fs);
sound(abs(envelopeSSBTC), Fs/5);
