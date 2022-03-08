##Clearing all the previous figures and windows
clear all;
close all;
clf;
pkg load signal;


[originalSignal, fs] = audioread('Monkey/monkey_sample.wav');
duration = length(originalSignal);
t = linspace(0, length(originalSignal)/fs, length(originalSignal));

##Number of Bits of Filter
BITS = 8;

##For loop describing good enough resolution of sampling
for n = 1:BITS;
  
  x =2^n - 1;
  
  ##Rounding a signal quantized value using "x" parameter
  quantizedSignal = round(originalSignal(1:duration)*x)/x;
  
  ##Noise signal
  noiseSignal = originalSignal(1:duration) - quantizedSignal;
  
  ##Ploting signals: Original, Quantized and Error of Sampling
  ##Blue - Original Signal
  ##Green - Quantized Signal
  ##Red - Noise signal
  
  plot(t, originalSignal(1:duration), '-b; Original signal;',
        t, quantizedSignal, '-g; Quantized signal;',
        t, noiseSignal, '-r; Noise;'); #Error of sampling  
  grid on;
  xlabel ('Time [s]');
  ylabel('Amplitude [-]');
  title({"Singular step of quantization number:", int2str(n)});  
        
  ##(n, :) - all "n" rows will be read in their entirety
  ##std - Standard Deviation
  ##logarythmic scale multiplied by 20 - voltage scale  
  snr(n,:) = 20 * log10(std(originalSignal(1:duration))/std(noiseSignal)); 

  print(["Monkey/Quantized_signal ", int2str(n), ".png"], "-dpng", "-color");
  
  roundedSignal = round(originalSignal*x)/x;
  
  ##Audio save function from folder
  audiowrite(['Monkey/monkey', int2str(n), '.wav'], roundedSignal, fs);
end

################################################################################
#Plotting the Original Signal
figure;
subplot(3,2,1); plot(t, originalSignal(1:duration), '-b');
title('Original Audio Signal in Time Domain');
xlabel('Time [s]'); ylabel('Amplitude [-]');
grid on; hold on;
 
#Plotting the SNR of the Final Quantized Signal
subplot(3,2,2); plot(snr, 'r*');
hold on; plot(snr ,'-b');
hold on; grid on;
title({'Signal-to-Noise Ratio'});
xlabel({'Resolution of quantization'}); ylabel({'SNR value', 'Decibels [dB]'});
 
##Overwriting (copying) final audio file to next one
##one with o number higher of 1
file1 = ['Monkey/monkey' , int2str(BITS), '.wav'];
file2 = ['Monkey/monkey' , int2str(BITS + 1), '.wav'];
copyfile(file1, file2); 

##Spectrum analisys of audio file after quantization and sampling
[finalQuantized, fp] = audioread(file2);

##Subplot of Final Quantized Signal
subplot(3,2,3);
plot(t, finalQuantized(1:duration));
title({'Final Quantized Signal in Time Domain'});
xlabel('Time [s]'); ylabel('Amplitude [-]');
hold on; grid on;

##Apply FFT on finalQuantized
fqFFT = fft(finalQuantized, duration);
fqFFT = fqFFT(1:duration);
fq_magnitudeFFT = fqFFT(1:duration/2);
x_magnitudeFFT = fp.*(0:duration/2-1)/duration;

subplot(3,2,4);
plot(x_magnitudeFFT, abs(fq_magnitudeFFT));
title({'Final Quantized Signal in Frequency Domain'});
xlabel('Frequency [Hz]'); ylabel('Amplitude [-]');
hold on; grid on;

##Plot Sampling Error (Noise)
subplot(3,2,5); plot(t,  noiseSignal(1:duration));
title('Noise in Time Domain'); xlabel('Time [s]'); 
ylabel('Amplitude [-]');
hold on; grid on;    

##Apply FFT on noise
noiseFFT = fft(noiseSignal, duration);
noiseFFT = noiseFFT(1:duration);
noise_magnitudeFFT = noiseFFT(1:duration/2);

subplot(3,2,6);
plot(x_magnitudeFFT, abs(noise_magnitudeFFT));
title({'Noise in Frequency Domain'});
xlabel('Frequency [Hz]'); ylabel('Amplitude [-]');
hold on; grid on;

print(['Monkey/Fig. 1 - Original Plots.png'], "-dpng", '-S1600,900');   

#############################################################################

##Designing FIR Filter
cutOff = [1500 3500]/(fs/2);
order = 510; #511 samples

##Getting the Filter Coefficients or Impulse Response
impulseRes = fir1(order, cutOff, 'pass');
figure;
subplot(2,2,1); stem(impulseRes);
title("Filter Coefficients");
xlabel("Time [s]"); ylabel("Amplitude [-]");
hold on; grid on;

##Apply FFT on Filter Coefficients to get the Frequency Response
freqRes = fft(impulseRes, duration);
freqRes = freqRes(1:duration);
freq_magnitudeFFT = freqRes(1:duration/2);

subplot(2,2,2);
plot(x_magnitudeFFT, abs(freq_magnitudeFFT));
title({'Frequency Response'});
xlabel('Frequency [Hz]'); ylabel('Amplitude [-]');
hold on; grid on;

##Applying the Designed Filter on the finalQuantized
fqFiltered = filter(impulseRes, 1, finalQuantized);
subplot(2,2,3);
plot(t, finalQuantized, '-m');
title("Filtered Audio Signal in Time Domain");
xlabel('Time [s]'); ylabel('Signal [-]');
grid on; hold on;

##Apply FFT on fqFiltered
filteredFFT = fft(fqFiltered, duration);
filteredFFT = filteredFFT(1:duration);
filtered_magnitudeFFT = filteredFFT(1:duration/2);

subplot(2,2,4);
plot(x_magnitudeFFT, abs(filtered_magnitudeFFT));
title({'Filtered Audio Signal in Frequency Domain'});
xlabel('Frequency [Hz]'); ylabel('Amplitude [-]');
hold on; grid on;

print(['Monkey/Fig. 2 - Filtered Plots.png'],"-dpng", '-S1600,900'); 

audiowrite(['Monkey/monkey_filtered.wav'], fqFiltered, fp);