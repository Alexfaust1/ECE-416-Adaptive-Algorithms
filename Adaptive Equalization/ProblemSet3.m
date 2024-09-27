%% Alexander Faust
%
% ECE 416 - Problem Set 3 - Adaptive Equalization
%
% December 1, 2023
clc; clear; close all;
%% Experiments
% We will examine the cases for alpha = 0.1, 0.2, 0.3 with Pdb = -30 and
% -10. RLS and inverse QRD-RLS will be used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST BENCH
% List experiment parameters:
N0 = 3;
M0 = 5;
K = 1e4;
Mmax = 11;
N_iter = 1;
N_train = Mmax - M0 + N_iter - 1;

% List algorithm parameters:
lambda = 0.9;                   % Step size
delta = 0.01;                   % Forgetting factor

% List experiment parameters for different alphas and PdB's:
Pdb = [-30 -10];
alphas = [0.1, 0.2, 0.3];
var_v = [10^(Pdb(1)/10), 10^(Pdb(2)/10)];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           VARIABLES GLOSSARY
%   x11 - alpha = 0.1, Pdb = -30            x12 - alpha = 0.1, Pdb = -10
%   x21 - alpha = 0.2, Pdb = -30            x22 - alpha = 0.2, Pdb = -10
%   x31 - alpha = 0.3, Pdb = -30            x32 - alpha = 0.3, Pdb = -10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % Generate the channel data for the given experiments:
% 1. Alpha = 0.1, Pdb = -30
[x11, y11, A_11, w0_11] = generateChannel(N_iter, N_train, N0, M0, K, Mmax, alphas(1), Pdb(1));

SNIR_theor11 = -10*log10(4*abs(alphas(1))^2 + var(1));
SNIR_opt11 = -10*log10(var_v(1));

delim;
disp("Theoretical SNIR for a = 0.1, Pdb = -30dB: " + SNIR_theor11);
disp("Optimal SNIR for a = 0.1, Pdb = -30dB: " + SNIR_opt11 + newline);

% Perform the RLS algorithm and examine its performance:
[w_RLS11, xi_RLS11, k_RLS_11, P_11] = doRLS(x11((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_11, lambda, delta, w0_11, N_iter);

% Estimate x and examine the calculated SNIR:
[SNIR_RLS_11_eq, SNIR_RLS_11_raw, x_decoded_RLS_11] = Xest(w_RLS11(:, end), y11, x11, Mmax, ...
                                                      N_iter, N_train, N0, K, "RLS - 11");
% Perform the inverse QRD-RLS algorithm and examine its performance:
[w_iQRD_11, xi_iQRD_11, k_iQRD_11, gamma11, Pch_11] = doinvQRDRLS(x11((Mmax - M0):(Mmax - M0 + N_iter - 1)), ...
                                        A_11, lambda, delta, w0_11, N_iter);
% Estimate x and examine calculated SNIR for iQRD-RLS:
[SNIR_iQRD_11_eq, SNIR_iQRD_11_raw, x_decoded_iQRD_11] = Xest(w_iQRD_11(:, end), y11, x11, Mmax, N_iter, N_train, N0, K, "iQRD RLS - 11");

% Compute the difference between P_11 and Pch * Pch' for computing spectral
% norm:
TEST1 = P_11 - Pch_11*Pch_11';

% Compute the max singular value of this matrix to check spectral norm:
CHECK1 = norm(TEST1);
disp("Spectral norm of P - P_ch*P_ch': " + CHECK1 + newline);
delim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENT 1 COMMENTS:
%
% It was suggested to check that the weights for iQRD-RLS and RLS are the
% same and upon inspection they are the same! This verifies to me that the
% algorithms have to be working as expected. The benefit here is that the
% computation of the w vector becomes easier to compute.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. Run next tests for alpha = 0.1, Pdb = -10dB:
[x12, y12, A_12, w0_12] = generateChannel(N_iter, N_train, N0, M0, K, Mmax, alphas(1), Pdb(2));

SNIR_theor12 = -10*log10(4*abs(alphas(1))^2 + var(2));
SNIR_opt12 = -10*log10(var_v(2));

delim;
disp("Theoretical SNIR for a = 0.1, Pdb = -10dB: " + SNIR_theor12);
disp("Optimal SNIR for a = 0.1, Pdb = -10dB: " + SNIR_opt12 + newline);

% Perform the RLS algorithm and examine its performance:
[w_RLS12, xi_RLS12, k_RLS_12, P_12] = doRLS(x12((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_12, lambda, delta, w0_12, N_iter);
% Estimate x and examine the calculated SNIR:
[SNIR_RLS_12_eq, SNIR_RLS_12_raw, x_decoded_RLS_12] = Xest(w_RLS12(:, end), y12, x12, Mmax, ...
                                                      N_iter, N_train, N0, K, "RLS - 11");
% Perform the inverse QRD-RLS algorithm and examine its performance:
[w_iQRD_12, xi_iQRD_12, k_iQRD_12, gamma12, Pch_12] = doinvQRDRLS(x12((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_12, lambda, delta, w0_12, N_iter);
% Estimate x and examine calculated SNIR for iQRD-RLS:
[SNIR_iQRD_12_eq, SNIR_iQRD_12_raw, x_decoded_iQRD_12] = Xest(w_iQRD_12(:, end), y12, x12, Mmax, N_iter, N_train, N0, K, "iQRD RLS - 12");

% Compute the difference between P_11 and Pch * Pch' for computing spectral
% norm:
TEST2 = P_12 - Pch_12*Pch_12';

% Compute the max singular value of this matrix to check spectral norm:
CHECK2 = norm(TEST2);
disp("Spectral norm of P - P_ch*P_ch': " + CHECK2 + newline);

delim;



% 3. Run tests for alpha = 0.2, Pdb = -30dB:
[x21, y21, A_21, w0_21] = generateChannel(N_iter, N_train, N0, M0, K, Mmax, alphas(2), Pdb(1));

SNIR_theor21 = -10*log10(4*abs(alphas(2))^2 + var(1));
SNIR_opt21 = -10*log10(var_v(1));

delim;
disp("Theoretical SNIR for a = 0.2, Pdb = -30dB: " + SNIR_theor21);
disp("Optimal SNIR for a = 0.2, Pdb = -30dB: " + SNIR_opt21 + newline);

% Perform the RLS algorithm and examine its performance:
[w_RLS21, xi_RLS21, k_RLS_21, P_21] = doRLS(x21((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_21, lambda, delta, w0_21, N_iter);
% Estimate x and examine the calculated SNIR:
[SNIR_RLS_21_eq, SNIR_RLS_21_raw, x_decoded_RLS_21] = Xest(w_RLS21(:, end), y21, x21, Mmax, ...
                                                      N_iter, N_train, N0, K, "RLS - 21");
% Perform the inverse QRD-RLS algorithm and examine its performance:
[w_iQRD_21, xi_iQRD_21, k_iQRD_21, gamma21, Pch_21] = doinvQRDRLS(x21((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_21, lambda, delta, w0_21, N_iter);
% Estimate x and examine calculated SNIR for iQRD-RLS:
[SNIR_iQRD_21_eq, SNIR_iQRD_21_raw, x_decoded_iQRD_21] = Xest(w_iQRD_21(:, end), y21, x21, Mmax, N_iter, N_train, N0, K, "iQRD RLS - 21");

% Compute the difference between P_11 and Pch * Pch' for computing spectral
% norm:
TEST3 = P_21 - Pch_21*Pch_21';

% Compute the max singular value of this matrix to check spectral norm:
CHECK3 = norm(TEST3);
disp("Spectral norm of P - P_ch*P_ch': " + CHECK3 + newline);

delim;




% 4. Run tests for alpha = 0.2, Pdb = -10dB:
[x22, y22, A_22, w0_22] = generateChannel(N_iter, N_train, N0, M0, K, Mmax, alphas(2), Pdb(2));

SNIR_theor22 = -10*log10(4*abs(alphas(2))^2 + var(2));
SNIR_opt22 = -10*log10(var_v(1));

delim;
disp("Theoretical SNIR for a = 0.2, Pdb = -10dB: " + SNIR_theor22);
disp("Optimal SNIR for a = 0.2, Pdb = -10dB: " + SNIR_opt22 + newline);

% Perform the RLS algorithm and examine its performance:
[w_RLS22, xi_RLS22, k_RLS_22, P_22] = doRLS(x22((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_22, lambda, delta, w0_22, N_iter);
% Estimate x and examine the calculated SNIR:
[SNIR_RLS_22_eq, SNIR_RLS_22_raw, x_decoded_RLS_22] = Xest(w_RLS22(:, end), y22, x22, Mmax, ...
                                                      N_iter, N_train, N0, K, "RLS - 22");
% Perform the inverse QRD-RLS algorithm and examine its performance:
[w_iQRD_22, xi_iQRD_22, k_iQRD_22, gamma22, Pch_22] = doinvQRDRLS(x22((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_22, lambda, delta, w0_22, N_iter);
% Estimate x and examine calculated SNIR for iQRD-RLS:
[SNIR_iQRD_22_eq, SNIR_iQRD_22_raw, x_decoded_iQRD_22] = Xest(w_iQRD_22(:, end), y22, x22, Mmax, N_iter, N_train, N0, K, "iQRD RLS - 22");

% Compute the difference between P_11 and Pch * Pch' for computing spectral
% norm:
TEST4 = P_22 - Pch_22*Pch_22';

% Compute the max singular value of this matrix to check spectral norm:
CHECK4 = norm(TEST4);
disp("Spectral norm of P - P_ch*P_ch': " + CHECK4 + newline);

delim;



% 5. Run tests for alpha = 0.3, Pdb = -30dB:
[x31, y31, A_31, w0_31] = generateChannel(N_iter, N_train, N0, M0, K, Mmax, alphas(3), Pdb(1));

SNIR_theor31 = -10*log10(4*abs(alphas(3))^2 + var(1));
SNIR_opt31 = -10*log10(var_v(1));

delim;
disp("Theoretical SNIR for a = 0.3, Pdb = -30dB: " + SNIR_theor31);
disp("Optimal SNIR for a = 0.3, Pdb = -30dB: " + SNIR_opt31 + newline);

% Perform the RLS algorithm and examine its performance:
[w_RLS31, xi_RLS31, k_RLS_31, P_31] = doRLS(x31((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_31, lambda, delta, w0_31, N_iter);
% Estimate x and examine the calculated SNIR:
[SNIR_RLS_31_eq, SNIR_RLS_31_raw, x_decoded_RLS_31] = Xest(w_RLS31(:, end), y31, x31, Mmax, ...
                                                      N_iter, N_train, N0, K, "RLS - 31");
% Perform the inverse QRD-RLS algorithm and examine its performance:
[w_iQRD_31, xi_iQRD_31, k_iQRD_31, gamma31, Pch_31] = doinvQRDRLS(x31((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_31, lambda, delta, w0_31, N_iter);
% Estimate x and examine calculated SNIR for iQRD-RLS:
[SNIR_iQRD_31_eq, SNIR_iQRD_31_raw, x_decoded_iQRD_31] = Xest(w_iQRD_31(:, end), y31, x31, Mmax, N_iter, N_train, N0, K, "iQRD RLS - 31");

% Compute the difference between P_11 and Pch * Pch' for computing spectral
% norm:
TEST5 = P_31 - Pch_31*Pch_31';

% Compute the max singular value of this matrix to check spectral norm:
CHECK5 = norm(TEST5);
disp("Spectral norm of P - P_ch*P_ch': " + CHECK5 + newline);

delim;



% 6. Run final test with alpha = 0.3, Pdb = -10dB:
[x32, y32, A_32, w0_32] = generateChannel(N_iter, N_train, N0, M0, K, Mmax, alphas(3), Pdb(2));

SNIR_theor32 = -10*log10(4*abs(alphas(3))^2 + var(2));
SNIR_opt32 = -10*log10(var_v(2));

delim;
disp("Theoretical SNIR for a = 0.3, Pdb = -10dB: " + SNIR_theor32);
disp("Optimal SNIR for a = 0.3, Pdb = -10dB: " + SNIR_opt32 + newline);

% Perform the RLS algorithm and examine its performance:
[w_RLS32, xi_RLS32, k_RLS_32, P_32] = doRLS(x32((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_32, lambda, delta, w0_32, N_iter);
% Estimate x and examine the calculated SNIR:
[SNIR_RLS_32_eq, SNIR_RLS_32_raw, x_decoded_RLS_32] = Xest(w_RLS32(:, end), y32, x32, Mmax, ...
                                                      N_iter, N_train, N0, K, "RLS - 32");
% Perform the inverse QRD-RLS algorithm and examine its performance:
[w_iQRD_32, xi_iQRD_32, k_iQRD_32, gamma32, Pch_32] = doinvQRDRLS(x32((Mmax - M0) : (Mmax - M0 + N_iter - 1)), ...
                                        A_32, lambda, delta, w0_32, N_iter);
% Estimate x and examine calculated SNIR for iQRD-RLS:
[SNIR_iQRD_32_eq, SNIR_iQRD_32_raw, x_decoded_iQRD_32] = Xest(w_iQRD_32(:, end), y32, x32, Mmax, N_iter, N_train, N0, K, "iQRD RLS - 32");

% Compute the difference between P_11 and Pch * Pch' for computing spectral
% norm:
TEST6 = P_32 - Pch_32*Pch_32';

% Compute the max singular value of this matrix to check spectral norm:
CHECK6 = norm(TEST6);
disp("Spectral norm of P - P_ch*P_ch': " + CHECK6 + newline);

delim;

%% Create SNIR chart for easy comparison:

% Convert numeric arrays to string arrays
SNIR_11 = string([SNIR_RLS_11_raw, SNIR_RLS_11_eq, SNIR_iQRD_11_raw, SNIR_iQRD_11_eq]);
SNIR_21 = string([SNIR_RLS_21_raw, SNIR_RLS_21_eq, SNIR_iQRD_21_raw, SNIR_iQRD_21_eq]);
SNIR_31 = string([SNIR_RLS_31_raw, SNIR_RLS_31_eq, SNIR_iQRD_31_raw, SNIR_iQRD_31_eq]);
SNIR_12 = string([SNIR_RLS_12_raw, SNIR_RLS_12_eq, SNIR_iQRD_12_raw, SNIR_iQRD_12_eq]);
SNIR_22 = string([SNIR_RLS_22_raw, SNIR_RLS_22_eq, SNIR_iQRD_22_raw, SNIR_iQRD_22_eq]);
SNIR_32 = string([SNIR_RLS_32_raw, SNIR_RLS_32_eq, SNIR_iQRD_32_raw, SNIR_iQRD_32_eq]);


Cases = ["SNIR_RLS_raw", "SNIR_RLS_eq", "SNIR_iQRD_raw", "SNIR_iQRD_eq"];

% Create the table
SNIR_table = table(Cases', SNIR_11', SNIR_21', SNIR_31', SNIR_12', SNIR_22', SNIR_32');

SNIR_table.Properties.VariableNames = {'Algorithm', 'Alpha 0.1 | Pdb -30dB', 'Alpha 0.2 | Pdb -30dB', 'Alpha 0.3 | Pdb -30dB', ...
                                       'Alpha 0.1 | Pdb -10dB', 'Alpha 0.2 | Pdb -10dB', 'Alpha 0.3 | Pdb -10dB'};


% Display the table
disp(SNIR_table);

%% One last experiment
% Take ground truth to be SNIR_opt11 and compute the deviation in SNIR up
% to a certain number of iterations:
%
% Test for alpha = 0.1, Pdb = -30dB:

N_test_iterations = 100;

% Pre-allocate:
SNIR_raw_to_plot = zeros(1, length(N_test_iterations));

for iter = 1:N_test_iterations

    N_train_TEST = Mmax - M0 + iter - 1;
    SNIR_raw = 0;
    SNIR_eq = 0;
    % Generate channel on each iteration:
    
    [x_TEST, y_TEST, A_TEST, w0_TEST] = generateChannel(iter, N_train_TEST, N0, M0, K, Mmax, alphas(1), Pdb(1));

    % Perform the RLS algorithm and examine its performance for each 
    % iteration:
    [w_RLS_TEST, ~, ~, ~] = doRLS(x_TEST((Mmax - M0) : (Mmax - M0 + iter - 1)), A_TEST, lambda, delta, w0_TEST, iter);
    
    % Estimate x and examine the calculated SNIR:
    [SNIR_RLS_TEST_eq, SNIR_RLS_TEST_raw, x_decoded_RLS_11] = Xest(w_RLS_TEST(:, end), y_TEST, x_TEST, Mmax, ...
                                                       iter, N_train_TEST, N0, K, "RLS - TEST " + num2str(iter));

    % Compute the deviation - limit - SNIR:
    SNIR_raw_to_plot(iter) = SNIR_opt11 - SNIR_RLS_TEST_raw;
    
end

% Create the plots of the SNIR computed by the RLS algorithm and see how it
% changes with number of iterations of the algorithm that we choose to run:

figure;
stem(1:N_test_iterations, SNIR_raw_to_plot);
title("SNIR_{raw} for RLS Compared to Optimal Channel SNIR");
xlabel("Iteration");
ylabel("|SNIR_{opt} - SNIR_{raw}|^2");

%% Final Experiment Comments
%
% In the quick experiment I ran above, I computed the difference (just a
% simple subtraction from the optimal SNIR i.e. 30dB since this would
% consider no noise or interference in the channel) and the results were
% definitely interesting. The figure suggests that after 10 iterations the
% algorithm begins to perform well and sustains this across further
% iterations.
%
% The ideal performance here would be when the y-axis values are zero
% indicating that the SNIR when using RLS is as good as it gets. i.e. the
% adaptive algorithm can mitigate the effects of noise and interference
% perfectly.






%% Comments on Experiment Runs
%
% Okay so the immediate result is that with PdB = -30dB each algorithm was
% able to estimate x and its delays exactly with zero difference in points
% even when alpha = 0.1, 0.2, and 0.3. However, changing to PdB = -10dB
% hurt the results significantly for each alpha. The trend for this case
% was that increasing alpha led to less deviation in the points so the most
% optimal of these was alpha = 0.3 with PdB = -10dB. 
%   *NOTE that these runs are stochastic so the previous statement might
%         not be true all the time.
%
% Additionally, the RLS and iQRD algorithms must providing the same final
% weight vectors because in the SNIR table they produce the same raw and
% equalized values for each experiment. Although, we can compare this to
% the theoretical and optimal SNIR's in each case to make comparisons. Of
% course, in the optimal case the SNIR will be 10dB or 30dB (depending on
% the Pdb chosen for the experiment - since optimally there will be no
% noise or interference meaning that just the 30dB or 10dB signal power is
% present). Obviously, the trend will be that the SNIR's in the table will
% be bounded below the optimal SNIR when we perform the algorithms with
% added noise. Thus it would be interesting to see what these values become
% if we change the noise variance.
%
% Note also how in just about every case that the SNIR produced from either
% RLS or iQRD RLS algorithms yield higher SNIR than the THEORETICAL
% calculation of SNIR for the current experiment. This is great because it
% shows that using RLS or iQRD for channel equalization gives better SNIR
% values than that predicted theoretically. 
%
% Was strange how my calculation of SNIR for the raw and eq cases were the
% same for both algorithms but I could not figure out why it was doing
% this. It was most likely just an error in the way I calculated this.
%
% However, I believe the point to these two algorithms is that they are
% supposed to produce the same final weight vector except that one of them
% has less computational cost in different scenarios

%% Changing to N = 50 runs: Comments
%
% To avoid re-running everything again in a new file I will just tweak the
% number of iterations to 50 in my testbench and report the following
% analyses here for anything interesting that I observe. 
%
% At first it was clear that increasing the number of iterations kept SNIR
% values practically the same across all runs. One thing to note was that
% the SNIR calculation for RLS or iQRD for the "raw" calculation decreased
% by around 0.2 on average while the "eq" calculation increased by around
% the same margin.
%
% So then I decided to lower the number of iterations to see if I could
% notice some sort of result in that regard. The SNIR for the 
% "eq" calculation stayed the same while the "raw" value decreased. Makes
% sense probably because the algorithm just needs time to run before
% monitoring performance.
%








%% Functions Created

% Function to generate the N_train + M0 + K points x, v
function [x, y, A, w0] = generateChannel(N_iter, N_train, N0, M0, K, Mmax, alpha, Pdb)
    % OUTPUTS:
    %   x   - "symbols" to send i.e. the +-1's equiprobable
    %   y   - output of the channel governed by the channel model
    %   w0  - initial tap weight vector
    %   A   - Matrix [Mmax x N_iter] dimensional s.t. columns of A are the
    %         u(n) vectors (Structure given in the problem set description)

    % Generate the x inputs to the channel model:
    x = 2*randi([0 1], N_train + M0 + K, 1) - 1;

    % Generate the noise 'v' terms:
    var_v = 10^(Pdb/10);
    v = sqrt(var_v) * randn(N_train + M0 + K, 1);

    % Prescribe the filter coefficients:
    h = [zeros(1, N0 - 1), alpha, 1, -alpha];

    % Filter the x inputs with the imulse response of 'y' the channel
    % model:
    y = filter(h, 1, x) + v;
    
    % Finally, construct the A matrix using Toeplitz structure, similar to
    % PS 1:
    %
    % NOTE: 'y' flipped forms the first column, while 'y' forms the first
    % row (with size prescribed in the PS PDF)
    %
    A = toeplitz(flipud(y(1 : Mmax)), (y(Mmax : (N_iter + Mmax - 1))));
    
    % Initialize and construct initial weight vector:
    w0 = zeros(Mmax, 1);
    w0(M0 + 1, 1) = 1;

end

% Function to calculate the x_est(n) vector based on the algorithm output
% at each K time step
function [SNIR_eq, SNIR_raw, x_decoded] = Xest(w_f, y, x, Mmax, N_iter, N_train, ...
                                               N0, K, ALG)
    % OUTPUTS:
    %   SNIR_eq     - Calculates the SNIR based on x_est output produced for 
    %                 each ALG
    %   SNIR_raw    - Calculates the SNIR based on the raw channel output
    %   x_decoded   - The decoded x based on the x_est produced by the ALG

    y_Mmax = toeplitz(flipud(y(1 : Mmax)), (y(Mmax : end)));

    % Compute x_est based on the equation in the PS description:
    % y_Mmax index starts at 
    x_est = w_f' * y_Mmax(:, N_iter + 1 : end);

    % Encode the estimate of x based on a simple sign detection:
    x_decoded = sign(x_est);

    % Create a flag to monitor the delay:
    FLAG = sum(x(N_train + 1 : N_train + K)' ~= x_decoded);
    disp("x[n + N_train] & x_decoded[n] differ by: " + FLAG + " points for " + ALG);

    % Calculate the SNIR for the equalizer output and raw channel output:
    SNIR_raw = -10*log10(mean(abs(x_est' - x(N_train + 1 : N_train + K)).^2 ));
    SNIR_eq = -10*log10(mean(abs(y(N0 + 1 : end) - x(1 : end - N0)).^2 ));
 
end

% Function to perform the RLS algorithm for equalization:
function [w, xi, k, Pfinal] = doRLS(d, A, lambda, delta, w0, N_iter)
        % OUTPUTS:
        %   w       - M x N_iter matrix s.t. each column contains the w tap
        %             weight vector for 1 <= n <= N_iter
        %   xi      - The a-priori estimation error
        %   k       - Matrix s.t. columns are the Kalman gain vectors
        %   Pfinal  - Final computed P(k) matrix i.e. inverse covariance
        %             matrix

        % Initialize RLS vectors:
        k = zeros(size(w0, 1), N_iter);
        P = eye(size(w0, 1)) / delta;
        xi = zeros(1, N_iter);
        W = [w0, zeros(size(w0, 1), N_iter)];

        % Perform N_iter iterations of RLS to determine the tap weights:
        for n = 1:N_iter
            % Note that at n = 1 at the beginning P[n-1] = P[0], and 
            % w[0] = 0, but we initialized it above to satisfy this. So
            % make a mental note that MATLAB "P[1]" = "P[0]" for RLS
            pi_n = P * A(:, n);
            k(:, n) = (lambda + A(:, n)' * pi_n)^-1 * pi_n;
            P = lambda^-1*P - lambda^-1*k(:, n)*A(:, n)' * P;
            xi(1, n) = d(n) - W(:, n)' * A(:, n);
            W(:, n + 1) = W(:, n) + k(:, n) * xi(1, n)';
        end

        % Update the output w vector at the end of RLS based on every entry
        % after the first column (since it's used as the initial weight:
        w = W(:, 2 : end);
        Pfinal = P;

end

% Function to perform inverse QRD-RLS algorithm
function [w, xi, k, gamma, Pch] = doinvQRDRLS(d, A, lambda, delta, w0, N_iter)
    % OUTPUTS:
    %   w, xi   - Same as doRLS
    %   gamma   - Sequence of conversion factors
    %   Pch     - Final value of Cholesky factor for the P matrix

    % Initialize QRD-RLS algorithm parameters similar to RLS algorithm:
    M = size(w0, 1);
    k = zeros(M, N_iter);
    W = [w0, zeros(M, N_iter)];
    xi = zeros(1, N_iter);
    gamma = zeros(1, N_iter);
    u = A;
    Phalf = zeros(M, M, N_iter);
    Phalf(:, :, 1) = delta^(-1/2) * eye(M);

    % Run QRD-RLS for N_iter iterations
    for n = 1:N_iter
        % Form the prearray matrix used in inverse QRD-RLS (in slides):
        prearray = [    1        lambda^(-1/2) * u(:, n)' * Phalf(:, :, n); ...
                    zeros(M, 1),      lambda^(-1/2) * Phalf(:, :, n)     ];
        % Map the prearray to the postarray via the unitary transform:
        R = (qr(prearray'))';
        
        % Create the post array by canceling the e^jtheta factors:
        postarray = R * diag(diag(R) ./ abs(diag(R)));

        % Extract the gamma^-1/2 terms from the computed postarray matrix:
        gamma_half = postarray(1, 1);

        % Now extract the gamma^-1/2 * k term from the bottom left of the
        % lower triangular postarray matrix:
        gamma_half_k = postarray(2 : end, 1);

        % Finally, extract the P^1/2 term from the bottom right of the
        % postarray:
        Phalf(:, :, n + 1) = postarray(2 : end, 2 : end);
        gamma(1, n) = gamma_half^(-2);

        % Can now calculate the Kalman gain and perform the main algorithm
        % iterations:
        k(:, n) = gamma_half^-1 * gamma_half_k;
        xi(n) = d(n) - W(:, n)' * u(:, n);
        W(:, n + 1) = W(:, n) + k(:, n) * xi(n)';

    end

    % Finally, set the output Pch matrix and w weights:
    Pch = Phalf(:, :, end);
    w = W(:, 2 : end);

end

% Function to display a page of #'s to separate tables:
function delim 
    disp("##################################################################");
end