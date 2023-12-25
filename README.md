# Experiment Overview

- [Experiment Overview](#experiment-overview)
  - [DSB (Double Sideband) - Experiment 1](#dsb-double-sideband---experiment-1)
    - [Overview of Files](#overview-of-files)
    - [Instructions for Running the Code](#instructions-for-running-the-code)
    - [Experiment Steps](#experiment-steps)
  - [SSB (Single Sideband) - Experiment 2](#ssb-single-sideband---experiment-2)
    - [Overview of Files](#overview-of-files-1)
    - [Instructions for Running the Code](#instructions-for-running-the-code-1)
    - [Experiment Steps](#experiment-steps-1)
  - [FM (Frequency Modulation) - Experiment 3](#fm-frequency-modulation---experiment-3)
    - [Overview of Files](#overview-of-files-2)
    - [Instructions for Running the Code](#instructions-for-running-the-code-2)
    - [Experiment Steps](#experiment-steps-2)

## DSB (Double Sideband) - Experiment 1

### Overview of Files

1. **[DSB.m](DSB.m)**: MATLAB script for Experiment 1. Key functionalities include:
   - Reading an audio signal from 'eric.wav'.
   - Filtering the signal to retain frequencies around 4 kHz.
   - Performing DSBSC and DSBTC modulations.
   - Implementing Coherent Detection with various SNRs and introducing frequency/phase errors.

2. **[eric.wav](eric.wav)**: Audio file used in the experiment.

### Instructions for Running the Code

1. Ensure MATLAB is installed on your system.
2. Place 'eric.wav' in the same directory as 'DSB.m'.
3. Run 'DSB.m'.

### Experiment Steps

1. **Original Signal Spectrum**: Displaying the spectrum of the original audio signal.
2. **Filtered Signal Spectrum**: Filtering the signal to 4 kHz and displaying the spectrum.
3. **DSBSC Modulation**: Generating and plotting DSBSC-modulated signals.
4. **DSBTC Modulation**: Generating and plotting DSBTC-modulated signals.
5. **Coherent Detection with SNR and Errors**: Demonstrating coherent detection with different SNR values and introducing frequency/phase errors.

## SSB (Single Sideband) - Experiment 2

### Overview of Files

1. **[SSB.m](SSB.m)**: MATLAB script for Experiment 2. Key functionalities include:
   - Reading an audio signal from 'eric.wav'.
   - Filtering the signal to retain frequencies around 4 kHz.
   - Performing SSBLSB modulation and Coherent Detection using an ideal and Butterworth filter.
   - Demonstrating Coherent Detection with varying SNRs.

2. **[eric.wav](eric.wav)**: Audio file used in the experiment.

### Instructions for Running the Code

1. Ensure MATLAB is installed on your system.
2. Place 'eric.wav' in the same directory as 'SSB.m'.
3. Run 'SSB.m'.

### Experiment Steps

1. **Original Signal Spectrum**: Displaying the spectrum of the original audio signal.
2. **Filtered Signal Spectrum**: Filtering the signal to 4 kHz and displaying the spectrum.
3. **SSBLSB Modulation and Coherent Detection (Ideal Filter)**: Generating and plotting SSBLSB-modulated signals and performing coherent detection using an ideal filter.
4. **SSBLSB Modulation and Coherent Detection (Butterworth Filter)**: Generating and plotting SSBLSB-modulated signals and performing coherent detection using a Butterworth filter.
5. **Coherent Detection with Different SNR Values**: Demonstrating coherent detection with varying SNR values.

## FM (Frequency Modulation) - Experiment 3

### Overview of Files

1. **[FM.m](FM.m)**: MATLAB script for Experiment 3. Key functionalities include:
   - Reading an audio signal from 'eric.wav'.
   - Filtering the signal to retain frequencies around 4 kHz.
   - Performing Frequency Modulation (FM) and Demodulation.

2. **[eric.wav](eric.wav)**: Audio file used in the experiment.

### Instructions for Running the Code

1. Ensure MATLAB is installed on your system.
2. Place 'eric.wav' in the same directory as 'FM.m'.
3. Run 'FM.m'.

### Experiment Steps

1. **Original Signal Spectrum**: Displaying the spectrum of the original audio signal.
2. **Filtered Signal Spectrum**: Filtering the signal to 4 kHz and displaying the spectrum.
3. **FM Modulation**: Calculating constants for FM modulation and displaying the spectrum of the modulated signal.
4. **FM Demodulation**: Implementing frequency discrimination and envelope detection for demodulation. Displaying the demodulated FM signal in the time domain.
5. **Resampling and Playback**: Resampling and playing the demodulated signal.