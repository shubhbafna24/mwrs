%% EXPERIMENT 9 - PERFORMANCE ANALYSIS OF RAYLEIGH FADING CHANNEL
%% PART 1 : 
% Generate a sequence of N=20,000 statistically independent and identically distributed Rayleigh random variables. 
% Plot the histogram for the 20000 symbols and compare it with the corresponding Rayleigh probability density function.

clc;
close all;
N = 20000;
SNR_limit = 35;
SNR_db = -5:0.5:SNR_limit;
SNR = 10.^(SNR_db/10);
u = rand(1, N);
m = floor(2 * rand(1, N));

var = 1; % sigma^2
nstd = sqrt(var);
y = zeros(1, N);

Pe_BPSK_sim = zeros(1, length(SNR));
Pe_BFSK_sim = zeros(1, length(SNR));
Pe_DPSK_sim = zeros(1, length(SNR));

% Generate Rayleigh random variable
r = sqrt(-(2 * var * log(u)));

figure(1);
histogram(r, 100);
title("Rayleigh Random Variable Histogram Plot");
xlabel("Random Variable R");
ylabel("Frequency");

a = 0:0.01:10;
R = (a / var) .* exp(-(a .* a) / (2 * var));

figure(2);
plot(a, R);
title("Rayleigh Probability Density Function (PDF)");
xlabel("Random Variable");
ylabel("Probability");
legend("Variance = 1");

% BPSK simulation
Pe_BPSK_id = 0.5 * (1 - sqrt((var * SNR) ./ (1 + var * SNR)));

% BFSK simulation (coherent)
Pe_BFSK_id = 0.5 * (1 - sqrt(var * SNR ./ (2 + (var * SNR))));

% DPSK simulation
Pe_DPSK_id = 0.5 ./ (1 + var * SNR);

% Comparison of Error Performance for AWGN and Rayleigh Fading Channels
Pe_BPSK_NF = 0.5 * erfc(sqrt(SNR));
Pe_BFSK_NF = 0.5 * erfc(sqrt(SNR / 2));
Pe_DPSK_NF = 0.5 * exp(-SNR); % Non-coherent

figure(3);
semilogy(SNR_db, Pe_BPSK_id, 'r.-', SNR_db, Pe_BFSK_id, 'r*-', SNR_db, Pe_DPSK_id, 'r--', SNR_db, Pe_BPSK_NF, 'b.-', SNR_db, Pe_BFSK_NF, 'b*-', SNR_db, Pe_DPSK_NF, 'b--');
axis([-5 SNR_limit 0.000001 1]);
title('Performance of BPSK, BFSK, and DPSK in Rayleigh Fading');
xlabel('SNR (dB)');
ylabel('Probability of Error');
legend('Pe of BPSK with fading', 'Pe of BFSK with fading', 'Pe of DPSK with fading', 'Pe of BPSK without fading', 'Pe of BFSK without fading', 'Pe of DPSK without fading');
