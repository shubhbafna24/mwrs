%% EXPERIMENT 9 - PERFORMANCE ANALYSIS OF RAYLEIGH FADING CHANNEL MODEL
%% PART 3 : 
% perform a Monte Carlo simulation to estimate and plot the error probability of a 
% binary orthogonal signaling (BFSK) communication system in Rayleigh fading.

%%Performance of Binary Modulation in Rayleigh Fading Channel
%%transmission through a frequency nonselective channel

clc;
clear;

Eb = 1; % Energy per bit
EbNo_dB = 0:5:35; % Vary the average SNR
No_over_2 = Eb * 10.^(-EbNo_dB / 10); % Noise power
sigma = 1; % Rayleigh parameter
var = sigma^2;
BER = zeros(1, length(EbNo_dB));

% Calculation of error probability using Monte Carlo simulation:
for i = 1:length(EbNo_dB)
    no_errors = 0;
    no_bits = 0;
    
    % Assumption: m = 0 (All zero codeword is transmitted)
    while no_errors <= 10
        u = rand;
        % rand returns a single uniformly distributed random number in (0,1)
        alpha = sigma * sqrt(-2 * log(u)); % Rayleigh distributed with variance set to unity
        noise = sqrt(No_over_2(i)) * randn; % Gaussian noise with zero mean and unit variance
        y = alpha * sqrt(Eb) + noise; % Simulate the input to the detector
        
        if y <= 0
            y_d = 1;
        else
            y_d = 0;
        end
        
        no_bits = no_bits + 1;
        no_errors = no_errors + y_d;
    end
    
    BER(i) = no_errors / no_bits; % Estimated error probability
end

% Calculation of error probability using the theoretical formula:
rho_b = Eb ./ No_over_2 * var;
P2 = 1/2 * (1 - sqrt(rho_b ./ (1 + rho_b))); % Theoretical value

% Plot the results:
semilogy(EbNo_dB, BER, "-*", EbNo_dB, P2, "-o");
title("Monte Carlo Simulation for Performance of BPSK Signal");
xlabel("Average SNR/bit (dB)");
ylabel("Error Probability");
legend("Monte Carlo Simulation", "Theoretical Value");
grid on;
