%% EXPERIMENT 10 - SIMULATION OF THE OFDM TRANSMITTER AND RECEIVER

clear all;
nFFT = 64; % FFT size
nDSC = 52; % Number of data subcarriers
nBitPerSym = 52; % Number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nSym = 10^4; % Number of symbols
EbN0dB = 0:10; % Bit to noise ratio
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % Convert to symbol-to-noise ratio

for ii = 1:length(EbN0dB)
    % Transmitter
    ipBit = rand(1, nBitPerSym * nSym) > 0.5; % Random 1's and 0's
    ipMod = 2 * ipBit - 1; % BPSK modulation: 0 -> -1, 1 -> +1
    ipMod = reshape(ipMod, nBitPerSym, nSym).'; % Grouping into multiple symbols
    
    % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
    xF = [zeros(nSym, 6), ipMod(:, 1:nBitPerSym/2), zeros(nSym, 1), ipMod(:, nBitPerSym/2+1:nBitPerSym), zeros(nSym, 5)];
    
    % Taking FFT, normalizing transmit power to 1
    xt = (nFFT/sqrt(nDSC)) * ifft(fftshift(xF.')).';
    
    % Appending cyclic prefix
    xt = [xt(:, 49:64) xt];
    
    % Concatenating multiple symbols to form a long vector
    xt = reshape(xt.', 1, nSym * 80);
    
    % Gaussian noise of unit variance, zero mean
    nt = (1/sqrt(2)) * (randn(1, nSym * 80) + 1j * randn(1, nSym * 80));
    
    % Adding noise, accounting for wasted energy due to cyclic prefix
    yt = sqrt(80/64) * xt + 10^(-EsN0dB(ii)/20) * nt;
    
    % Receiver
    yt = reshape(yt.', 80, nSym).'; % Formatting received vector into symbols
    yt = yt(:, 17:80); % Removing cyclic prefix
    
    % Converting to frequency domain
    yF = (sqrt(nDSC)/nFFT) * fftshift(fft(yt.')).';
    yMod = yF(:, [6+(1:nBitPerSym/2), 7+(nBitPerSym/2+1:nBitPerSym)]);
    
    % BPSK demodulation: +ve value -> 1, -ve value -> -1
    ipModHat = 2 * floor(real(yMod/2)) + 1;
    ipModHat(ipModHat > 1) = 1;
    ipModHat(ipModHat < -1) = -1;
    
    % Converting modulated values into bits
    ipBitHat = (ipModHat + 1) / 2;
    ipBitHat = reshape(ipBitHat.', nBitPerSym * nSym, 1).';
    
    % Counting the errors
    nErr(ii) = sum(ipBitHat ~= ipBit);
end

% BER calculation
simBer = nErr / (nSym * nBitPerSym);
theoryBer = (1/2) * erfc(sqrt(10.^(EbN0dB/10)));

% Plot results
close all;
figure;
semilogy(EbN0dB, theoryBer, 'bs-', 'LineWidth', 2);
hold on;
semilogy(EbN0dB, simBer, 'mx-', 'LineWidth', 2);
axis([0 10 10^-5 1]);
grid on;
legend('Theory', 'Simulation');
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate');
title('Bit Error Probability Curve for BPSK using OFDM');


