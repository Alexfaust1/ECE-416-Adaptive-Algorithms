%% Alexander Faust
%
% ECE 416 - Question 5
%
% December 1, 2023
clc; clear; close all;

%% Question 5
% In this section we will repeat questions 2-4 but for noise power 10 dB
% above the primary source - thus to save me a lot of time I am going to
% copy and paste my code from the ProblemSet2.m file while changing the
% noise power to +10dB (up from -12 dB) and examine these results. I will
% just summarize my comments at the very end in order to save time and
% simply delete and copy and pasted comments from the prior section:
%% SVD and MUSIC / MVDR Spectra

% Define parameters
d_lambda = 0.5;                         % d/lambda
d = 1;                                  % Assumes unit distance between sensors
N = 100;                                % Number of snapshots
signal_power = [0, -4, -8];             % Signal and noise power levels in dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WE HAVE TWEAKED THIS PARAMETER
noise_power = 10;                       % Noise power in dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AOA = [15, 20, 30 ; 30, 40, -40];       % Source AOA's (in degrees)
num_sources = 3;                        % Assuming 3 sources at 3 AOA's

% Generate sensor locations on the x and y axes in cross formation:
r = d*[(-10:1:10) zeros(1, 20); zeros(1, 21) (-10:1:-1) (1:1:10); zeros(1, 41)];
num_sensors = size(r, 2);

% Compute values of X:
[X, R_est, R_ideal, evalR_est, evalR_ideal, S] = generateX(num_sensors, N, d, d_lambda, r, AOA, signal_power, noise_power);

% Compare the fourth and third singular values of X
[U, sing_vals, V] = svd(X);
sing_vals = diag(sing_vals);

disp("Displaying the ratio between the fourth and third singular values of X: ");
sing_ratio = sing_vals(4) / sing_vals(3);
disp(sing_ratio);
disp("Displaying the ratio between the fourth and third eigen values of X: ");
eig_ratio = abs(evalR_ideal(4)) / abs(evalR_ideal(3)); 
disp(eig_ratio);

% Plot singular values
stem(sing_vals);
xlabel('Singular Value Index');
ylabel('Singular Value Magnitude');
title('Singular Value Plot');
grid on;

% Part b) - Compute eigenvalues of correlation matrix in descending
%           order and plot

figure;
stem(evalR_est);
title("Eigenvalues of correlation matrix in descending order");
xlabel("Position");
ylabel("eigenvalue");

disp("On inspection of the stem plot there are indeed 3 dominant eigenvalues. The value" + newline + ...
     "ratio between the fourth and third largest eigenvalue is:   " + eig_ratio + newline);

% Part c) - Compute projection onto noise subspace from SVD:
P_N = eye(num_sensors) - (U(:, 1:3)*U(:, 1:3)');
% Compute residual values of steering vectors times noise projection
% subspace:
residuals = abs(P_N * S);
disp("Displaying residual values of P_N noise projection times steering vecs: ");
disp(residuals);

% Part d) - Compute MVDR and MUSIC spectra
% Create range of theta and phi AOA's to compute from
thetas = 0:1:90;
phis = -180:1:180;

MUSIC_X = MUSIC(d, d_lambda, P_N, r, thetas, phis);
MVDR_X = MVDR(d, d_lambda, R_est, r, thetas, phis);

plotSpectrum(MUSIC_X, thetas, phis, "MUSIC")
plotSpectrum(MVDR_X, thetas, phis, "MVDR");

% Part e & f) - Take theta = 20 deg slice and obtain 1D plot for phi between
%               -180 <= phi <= 180 deg. THen get plot of 0 <= theta <= 90
%               deg for phi = 30 deg for both MUSIC and MVDR
MUSIC_theta20 = MUSIC(d, d_lambda, P_N, r, AOA(1, 2), phis);
MVDR_theta20 = MVDR(d, d_lambda, R_est, r, AOA(1, 2), phis);

plot1DSpectrum(MUSIC_theta20, phis, "1D MUSIC", "theta", "phi", 20, 40);
plot1DSpectrum(MVDR_theta20, phis, "1D MVDR", "theta", "phi", 20, 40);

MUSIC_phi30 = MUSIC(d, d_lambda, P_N, r, thetas, AOA(2, 1));
MVDR_phi30 = MVDR(d, d_lambda, R_est, r, thetas, AOA(2, 1));

plot1DSpectrum(MUSIC_phi30, thetas, "1D MUSIC", "phi", "theta", 30, 15);
plot1DSpectrum(MVDR_phi30, thetas, "1D MVDR", "phi", "theta", 30, 15);
%
% Part g) - Minimum for each MUSIC and MVDR would be the minimum eigenvalue
%           of their R correlation matrix and noise projection matrix,
%           respectively:
disp("Displaying theoretical minimum of MVDR and MUSIC spectrum, respectively." + newline + ...
     "note that these values should be the minimum eigenvalue of each: ");
disp(min(eig(R_est)));
disp(min(eig(P_N)));

% Part h) - Plug in 3 steering vectors into MUSIC and MVDR to compare lowre
%           bounds from that foudn in part g)
MUSIC_3steer = MUSIC(d, d_lambda, P_N, r, AOA(1, :), AOA(2, :));
MVDR_3steer = MVDR(d, d_lambda, R_est, r, AOA(1, :), AOA(2, :));

disp("Displaying calculated maximum and minimum of MVDR and MUSIC spectrum, respectively." + newline + ...
     "note that these values should be the minimum eigenvalue of each: ");
disp("Max, min MUSIC:" + max(diag(MUSIC_3steer)) + ", " + min(diag(MUSIC_3steer)));
disp("Max, min MVDR:" + max(diag(MVDR_3steer)) + ", " + min(diag(MVDR_3steer)));
%%
% Part i) - Find overall minimum values and compare to lower bounds for the
%           MUSIC and MVDR spectrums generated for all theta and phi's
min_MUSIC = abs(min(MUSIC_X, [], 'all'));
min_MVDR = abs(min(MVDR_X, [], "all"));
theoretical_min_MUSIC = 1/max(eig(P_N), [], "all");
theoretical_min_MVDR = 1/max(eig(R_est^-1), [], "all");

%% Optimal Beamforming: GSC
% List g unit vectors in each direction:
gx = [1 ; 0 ; 0];
gy = [0 ; 1 ; 0];
gz = [0 ; 0 ; 1];

% Compute the ideal GSC weights for the sources:
w_GSC_source1 = GSCWeights(num_sensors, num_sources, gx, S, R_est);
w_GSC_source2 = GSCWeights(num_sensors, num_sources, gy, S, R_est);
w_GSC_source3 = GSCWeights(num_sensors, num_sources, gz, S, R_est);

% Compute the 2-D array response of each source:
GSC_BF1 = arrayResponse(d, d_lambda, r, w_GSC_source1, thetas, phis);
GSC_BF2 = arrayResponse(d, d_lambda, r, w_GSC_source2, thetas, phis);
GSC_BF3 = arrayResponse(d, d_lambda, r, w_GSC_source3, thetas, phis);

% Create the 2-D plots for the GSC beamformers:
plotSpectrum(GSC_BF1, thetas, phis, "GSC Beamformer Source 1");
plotSpectrum(GSC_BF2, thetas, phis, "GSC Beamformer Source 2");
plotSpectrum(GSC_BF3, thetas, phis, "GSC Beamformer Source 3");

% Compute the 1-D array response of each source (take theta = 20 case):
GSC_BF1_theta20 = arrayResponse(d, d_lambda, r, w_GSC_source1, AOA(1, 2), phis);
GSC_BF2_theta20 = arrayResponse(d, d_lambda, r, w_GSC_source2, AOA(1, 2), phis);
GSC_BF3_theta20 = arrayResponse(d, d_lambda, r, w_GSC_source3, AOA(1, 2), phis);

% Compute the 1-D array response for the phi = 30 case:
GSC_BF1_phi30 = arrayResponse(d, d_lambda, r, w_GSC_source1, thetas, AOA(2, 1));
GSC_BF2_phi30 = arrayResponse(d, d_lambda, r, w_GSC_source2, thetas, AOA(2, 1));
GSC_BF3_phi30 = arrayResponse(d, d_lambda, r, w_GSC_source3, thetas, AOA(2, 1));

%% Optimal Beamforming: MVDR
% Calculate ideal weights for the MVDR beamformer:
w_MVDR_source1 = MVDRWeights(S, 1, R_est);
w_MVDR_source2 = MVDRWeights(S, 2, R_est);
w_MVDR_source3 = MVDRWeights(S, 3, R_est);

% Compute the 2-D array response at the AOA's:
MVDR_BF1 = arrayResponse(d, d_lambda, r, w_MVDR_source1, thetas, phis);
MVDR_BF2 = arrayResponse(d, d_lambda, r, w_MVDR_source2, thetas, phis);
MVDR_BF3 = arrayResponse(d, d_lambda, r, w_MVDR_source3, thetas, phis);

% Create the 2-D plots for the MVDR beamformer:
plotSpectrum(MVDR_BF1, thetas, phis, "MVDR Beamformer Source 1");
plotSpectrum(MVDR_BF2, thetas, phis, "MVDR Beamformer Source 2");
plotSpectrum(MVDR_BF3, thetas, phis, "MVDR Beamformer Source 3");

% Compute the 1-D array response for the theta = 20 case and show plots:
MVDR_BF1_theta20 = arrayResponse(d, d_lambda, r, w_MVDR_source1, AOA(1, 2), phis);
MVDR_BF2_theta20 = arrayResponse(d, d_lambda, r, w_MVDR_source2, AOA(1, 2), phis);
MVDR_BF3_theta20 = arrayResponse(d, d_lambda, r, w_MVDR_source3, AOA(1, 2), phis);

% Compute the 1-D array response for the phi = 30 case:
MVDR_BF1_phi30 = arrayResponse(d, d_lambda, r, w_MVDR_source1, thetas, AOA(2, 1));
MVDR_BF2_phi30 = arrayResponse(d, d_lambda, r, w_MVDR_source2, thetas, AOA(2, 1));
MVDR_BF3_phi30 = arrayResponse(d, d_lambda, r, w_MVDR_source3, thetas, AOA(2, 1));


%% MVDR Results Section - Tables and Plots
% Create 2-D plots for optimal MVDR beamformer:
plotArrayResponse(MVDR_BF1, "Optimal MVDR Beamformer", "source 1", AOA, thetas, phis, 2);
plotArrayResponse(MVDR_BF2, "Optimal MVDR Beamformer", "source 2", AOA, thetas, phis, 2);
plotArrayResponse(MVDR_BF3, "Optimal MVDR Beamformer", "source 3", AOA, thetas, phis, 2);

% Creating 1-D plots for the MVDR beamformer when theta = 20 fixed:
plot1DSpectrum(MVDR_BF1_theta20, phis, "1D MVDR", "theta", "phi", 20, 40);
plot1DSpectrum(MVDR_BF2_theta20, phis, "1D MVDR", "theta", "phi", 20, 40);
plot1DSpectrum(MVDR_BF3_theta20, phis, "1D MVDR", "theta", "phi", 20, 40);

% Creating 1-D plots for the MVDR beamformer when phi = 30 fixed:
plot1DSpectrum(MVDR_BF1_phi30, thetas, "1D MVDR", "phi", "theta", 30, 15);
plot1DSpectrum(MVDR_BF2_phi30, thetas, "1D MVDR", "phi", "theta", 30, 15);
plot1DSpectrum(MVDR_BF3_phi30, thetas, "1D MVDR", "phi", "theta", 30, 15);

% Table for steering MVDR beamformer:
MVDRSteerSource1 = SourceShift(MVDR_BF1, AOA);
MVDRSteerSource2 = SourceShift(MVDR_BF2, AOA);
MVDRSteerSource3 = SourceShift(MVDR_BF3, AOA);

Position = ["Array Response (Source 1)"; "Array Response (Source 2)"; "Array Response (Source 3)"];
MVDR_table = table(Position, MVDRSteerSource1, MVDRSteerSource2, MVDRSteerSource3);
disp(MVDR_table);

disp("Based on the table above, we can see that as the beam former is steered to each source, " + newline + ...
     "the array response reaches unity for the current source and highly suppresses the other sources.");


%% GSC Results Section - Tables and Plots
% Create 2-D plots for optimal GSC beamformer:
plotArrayResponse(GSC_BF1, "Optimal GSC Beamformer", "source 1", AOA, thetas, phis, 2);
plotArrayResponse(GSC_BF2, "Optimal MVDR Beamformer", "source 2", AOA, thetas, phis, 2);
plotArrayResponse(GSC_BF3, "Optimal MVDR Beamformer", "source 3", AOA, thetas, phis, 2);

% Creating 1-D plots for the GSC beamformer when theta = 20 fixed:
plot1DSpectrum(GSC_BF1_theta20, phis, "1D GSC", "theta", "phi", 20, 40);
plot1DSpectrum(GSC_BF2_theta20, phis, "1D GSC", "theta", "phi", 20, 40);
plot1DSpectrum(GSC_BF3_theta20, phis, "1D GSC", "theta", "phi", 20, 40);

% Creating 1-D plots for the GSC beamformer when phi = 30 fixed:
plot1DSpectrum(GSC_BF1_phi30, thetas, "1D GSC", "phi", "theta", 30, 15);
plot1DSpectrum(GSC_BF2_phi30, thetas, "1D GSC", "phi", "theta", 30, 15);
plot1DSpectrum(GSC_BF3_phi30, thetas, "1D GSC", "phi", "theta", 30, 15);

% Table for steering GSC beamformer:
GSCSteerSource1 = SourceShift(GSC_BF1, AOA);
GSCSteerSource2 = SourceShift(GSC_BF2, AOA);
GSCSteerSource3 = SourceShift(GSC_BF3, AOA);

GSC_table = table(Position, GSCSteerSource1, GSCSteerSource2, GSCSteerSource3);
disp(GSC_table);

disp("Based on the table above, we can see that as the beam former is steered to each source, " + newline + ...
     "the array response reaches unity for the current source and highly suppresses the other sources.");

%% Adaptive Beamforming - GSC using LMS
% List LMS parameters to be used:
TRIALS = 100;                   % Run 100 experiments/trials of LMS
mu = 0.5;                       % Indicate step size of LMS to run

% Compute the adaptive GSC weights for source 1 through source 3 using LMS:
[w_avg_GSC_LMS_sc1, J_GSC_sc1, w_GSC_LMS_sc1, e_GSC_LMS_sc1] = GSCWeights_LMS(num_sensors, N, num_sources, ...
                d, d_lambda, r, gx, AOA, signal_power, noise_power, TRIALS, mu);
[w_avg_GSC_LMS_sc2, J_GSC_sc2, w_GSC_LMS_sc2, e_GSC_LMS_sc2] = GSCWeights_LMS(num_sensors, N, num_sources, ...
                d, d_lambda, r, gy, AOA, signal_power, noise_power, TRIALS, mu);
[w_avg_GSC_LMS_sc3, J_GSC_sc3, w_GSC_LMS_sc3, e_GSC_LMS_sc3] = GSCWeights_LMS(num_sensors, N, num_sources, ...
                d, d_lambda, r, gz, AOA, signal_power, noise_power, TRIALS, mu);

% Now compute the array responses for each adaptive GSC beamformer:
GSC_Asource1_LMS = arrayResponse(d, d_lambda, r, w_avg_GSC_LMS_sc1(:, N, 1), thetas, phis);
GSC_Asource2_LMS = arrayResponse(d, d_lambda, r, w_avg_GSC_LMS_sc2(:, N, 1), thetas, phis);
GSC_Asource3_LMS = arrayResponse(d, d_lambda, r, w_avg_GSC_LMS_sc3(:, N, 1), thetas, phis);

% Compute array response after just one run of LMS as requested:
GSC_Asource1_LMS_oneRun = arrayResponse(d, d_lambda, r, w_avg_GSC_LMS_sc1(:, 1, 1), thetas, phis);

% Plot the GSC array response for a single run of LMS and for 100 runs of LMS:
plotArrayResponse(GSC_Asource1_LMS_oneRun, "Adaptive LMS GSC", "source 1", AOA, thetas, phis, 2);
plotArrayResponse(GSC_Asource1_LMS, "Adaptive LMS GSC", "source 1", AOA, thetas, phis, 2);

% Plot array responses for the source 2 and source 3 GSC's:
plotArrayResponse(GSC_Asource2_LMS, "Adaptive LMS GSC", "source 2", AOA, thetas, phis, 2);
plotArrayResponse(GSC_Asource3_LMS, "Adaptive LMS GSC", "source 3", AOA, thetas, phis, 2);

% Okay now we will plot the learning curves for each GSC LMS run to see
% where convergence can be observed:
plotLearningCurve(abs(e_GSC_LMS_sc1(:, :, 1)).^2, N, signal_power(1), "ONE-RUN LMS GSC", "Source 1");
plotLearningCurve(J_GSC_sc1, N, signal_power(1), "LMS GSC", "Source 1");

plotLearningCurve(J_GSC_sc2, N, signal_power(2), "LMS GSC", "Source 2");
plotLearningCurve(J_GSC_sc3, N, signal_power(3), "LMS GSC", "Source 3");


% Now compare values between N0 and N0/2 for GSC LMS (i.e. N0 being the
% number of iterations towards the end where LMS seems to already be stable
% and "converged"
N0_GSC_LMS = 10;
[GSC_Asource1_LMS_halfN0] = arrayResponse(d, d_lambda,r, w_avg_GSC_LMS_sc1(:, N0_GSC_LMS/2, 1), ...
                                          thetas, phis);
[GSC_Asource2_LMS_halfN0] = arrayResponse(d, d_lambda,r, w_avg_GSC_LMS_sc2(:, N0_GSC_LMS/2, 1), ...
                                          thetas, phis);
[GSC_Asource3_LMS_halfN0] = arrayResponse(d, d_lambda,r, w_avg_GSC_LMS_sc3(:, N0_GSC_LMS/2, 1), ...
                                          thetas, phis);

% Now calculate the array response for the shifted sources:
GSC_LMS_SteerSource1_halfN0 = SourceShift(GSC_Asource1_LMS_halfN0, AOA);
GSC_LMS_SteerSource2_halfN0 = SourceShift(GSC_Asource2_LMS_halfN0, AOA);
GSC_LMS_SteerSource3_halfN0 = SourceShift(GSC_Asource3_LMS_halfN0, AOA);

disp("Showing table for GSC with LMS for N0 then N0/2:")
table_GSC_LMS_halfN0 = table(Position, GSC_LMS_SteerSource1_halfN0, GSC_LMS_SteerSource2_halfN0, GSC_LMS_SteerSource3_halfN0);
disp(table_GSC_LMS_halfN0);
delim;

% Now calculating array responses for N0 = 10:
[GSC_Asource1_LMS_N0] = arrayResponse(d, d_lambda,r, w_avg_GSC_LMS_sc1(:, N0_GSC_LMS, 1), ...
                                          thetas, phis);
[GSC_Asource2_LMS_N0] = arrayResponse(d, d_lambda,r, w_avg_GSC_LMS_sc2(:, N0_GSC_LMS, 1), ...
                                          thetas, phis);
[GSC_Asource3_LMS_N0] = arrayResponse(d, d_lambda,r, w_avg_GSC_LMS_sc3(:, N0_GSC_LMS, 1), ...
                                          thetas, phis);

% Shift the sources:
GSC_LMS_SteerSource1_N0 = SourceShift(GSC_Asource1_LMS_N0, AOA);
GSC_LMS_SteerSource2_N0 = SourceShift(GSC_Asource2_LMS_N0, AOA);
GSC_LMS_SteerSource3_N0 = SourceShift(GSC_Asource3_LMS_N0, AOA);

table_GSC_LMS_N0 = table(Position, GSC_LMS_SteerSource1_N0, GSC_LMS_SteerSource2_N0, GSC_LMS_SteerSource3_N0);
disp(table_GSC_LMS_N0);
delim;

disp("From the above tables for LMS for GSC, convergence factor doesn't really impact" + newline + ...
     "performance since it converges so fast!");

%% Adaptive Beamforming - MVDR using LMS

% Compute the adaptive MVDR weights for source 1 through source 3 using LMS:
[w_avg_MVDR_LMS_sc1, J_MVDR_sc1, w_MVDR_LMS_sc1, e_MVDR_LMS_sc1] = MVDRWeights_LMS(num_sensors, N, num_sources, ...
                d, d_lambda, r, 1, AOA, signal_power, noise_power, TRIALS, mu);
[w_avg_MVDR_LMS_sc2, J_MVDR_sc2, w_MVDR_LMS_sc2, e_MVDR_LMS_sc2] = MVDRWeights_LMS(num_sensors, N, num_sources, ...
                d, d_lambda, r, 2, AOA, signal_power, noise_power, TRIALS, mu);
[w_avg_MVDR_LMS_sc3, J_MVDR_sc3, w_MVDR_LMS_sc3, e_MVDR_LMS_sc3] = MVDRWeights_LMS(num_sensors, N, num_sources, ...
                d, d_lambda, r, 3, AOA, signal_power, noise_power, TRIALS, mu);

% Now we'll go ahead and compute the array responses for the MVDR
% beamformers for each source. Then for source 1 only we will show the
% beamformer for only one run of LMS
MVDR_Asource1_LMS = arrayResponse(d, d_lambda, r, w_avg_MVDR_LMS_sc1(:, N), thetas, phis);
MVDR_Asource1_LMS_oneRun = arrayResponse(d, d_lambda, r, w_avg_MVDR_LMS_sc1(:, N, 1), thetas, phis);

MVDR_Asource2_LMS = arrayResponse(d, d_lambda, r, w_avg_MVDR_LMS_sc2(:, N), thetas, phis);
MVDR_Asource3_LMS = arrayResponse(d, d_lambda, r, w_avg_MVDR_LMS_sc3(:, N), thetas, phis);

% Alright now show the plots for the results of LMS MVDR for each of the 3
% sources, but also the results for one run of MVDR for source 1
plotArrayResponse(MVDR_Asource1_LMS_oneRun, "Adaptive LMS MVDR - ONE RUN", "source 1", AOA, thetas, phis, 2);
plotArrayResponse(MVDR_Asource1_LMS, "Adaptive LMS MVDR", "source 1", AOA, thetas, phis, 2);

plotArrayResponse(MVDR_Asource2_LMS, "Adaptive LMS MVDR", "source 2", AOA, thetas, phis, 2);
plotArrayResponse(MVDR_Asource3_LMS, "Adaptive LMS MVDR", "source 3", AOA, thetas, phis, 2);

% Plotting Learning curves now...
plotLearningCurve(abs(e_MVDR_LMS_sc1(:, :, 1)).^2, N, signal_power(1), "ONE-RUN LMS MVDR", "Source 1");
plotLearningCurve(J_MVDR_sc1, N, signal_power(1), "LMS MVDR", "Source 1");

plotLearningCurve(J_MVDR_sc2, N, signal_power(2), "LMS MVDR", "Source 2");
plotLearningCurve(J_MVDR_sc3, N, signal_power(3), "LMS MVDR", "Source 3");

% Compute now the array responses after convergence has stabilized:
N0_MVDR_LMS = 26;
[MVDR_Asource1_LMS_halfN0] = arrayResponse(d, d_lambda,r, w_avg_MVDR_LMS_sc1(:, N0_MVDR_LMS/2, 1), ...
                                          thetas, phis);
[MVDR_Asource2_LMS_halfN0] = arrayResponse(d, d_lambda,r, w_avg_MVDR_LMS_sc2(:, N0_MVDR_LMS/2, 1), ...
                                          thetas, phis);
[MVDR_Asource3_LMS_halfN0] = arrayResponse(d, d_lambda,r, w_avg_MVDR_LMS_sc3(:, N0_MVDR_LMS/2, 1), ...
                                          thetas, phis);

% Examine the results for when the beam formers are steered to different
% sources:
MVDR_LMS_SteerSource1_halfN0 = SourceShift(MVDR_Asource1_LMS_halfN0, AOA);
MVDR_LMS_SteerSource2_halfN0 = SourceShift(MVDR_Asource2_LMS_halfN0, AOA);
MVDR_LMS_SteerSource3_halfN0 = SourceShift(MVDR_Asource3_LMS_halfN0, AOA);

disp("Showing table for MVDR with LMS for N0 then N0/2:")
table_MVDR_LMS_halfN0 = table(Position, MVDR_LMS_SteerSource1_halfN0, MVDR_LMS_SteerSource2_halfN0, MVDR_LMS_SteerSource3_halfN0);
disp(table_MVDR_LMS_halfN0);
delim;

[MVDR_Asource1_LMS_N0] = arrayResponse(d, d_lambda,r, w_avg_MVDR_LMS_sc1(:, N0_MVDR_LMS, 1), ...
                                          thetas, phis);
[MVDR_Asource2_LMS_N0] = arrayResponse(d, d_lambda,r, w_avg_MVDR_LMS_sc2(:, N0_MVDR_LMS, 1), ...
                                          thetas, phis);
[MVDR_Asource3_LMS_N0] = arrayResponse(d, d_lambda,r, w_avg_MVDR_LMS_sc3(:, N0_MVDR_LMS, 1), ...
                                          thetas, phis);

% Examine the results for when the beam formers are steered to different
% sources:
MVDR_LMS_SteerSource1_N0 = SourceShift(MVDR_Asource1_LMS_N0, AOA);
MVDR_LMS_SteerSource2_N0 = SourceShift(MVDR_Asource2_LMS_N0, AOA);
MVDR_LMS_SteerSource3_N0 = SourceShift(MVDR_Asource3_LMS_N0, AOA);

table_MVDR_LMS_N0 = table(Position, MVDR_LMS_SteerSource1_N0, MVDR_LMS_SteerSource2_N0, MVDR_LMS_SteerSource3_N0);
disp(table_MVDR_LMS_N0);
delim;

% Finally, compute the average array response for N0-9 <= n <= N0
iterations = 10;
MVDR_Asource1_avg_LMS = avgArrayResponse(d, d_lambda, r, w_avg_MVDR_LMS_sc1, ...
                                        N0_MVDR_LMS, iterations, thetas, phis);
MVDR_Asource2_avg_LMS = avgArrayResponse(d, d_lambda, r, w_avg_MVDR_LMS_sc2, ...
                                        N0_MVDR_LMS, iterations, thetas, phis);
MVDR_Asource3_avg_LMS = avgArrayResponse(d, d_lambda, r, w_avg_MVDR_LMS_sc3, ...
                                        N0_MVDR_LMS, iterations, thetas, phis);
figure;
surf(phis, thetas, MVDR_Asource1_avg_LMS);

%% Adaptive Beamforming - MVDR using RLS
% Like LMS, RLS is an adaptive algorithm we can employ to achieve a faster
% rate of convergence. Let us see how these experiments go!

% LIST RLS PARAMETERS FOR EXPERIMENTS
lamb_MVDR_RLS = 0.95;
del_MVDR = 80;

% First compute the adaptive MVDR weights for source 1 through source 3 using LMS:
[w_avg_MVDR_RLS_sc1, J_MVDR_RLS_sc1, w_MVDR_RLS_sc1, e_MVDR_RLS_sc1] = MVDRWeights_RLS(num_sensors, N, num_sources, ...
                d, d_lambda, r, 1, lamb_MVDR_RLS, AOA, signal_power, noise_power, TRIALS, del_MVDR);
[w_avg_MVDR_RLS_sc2, J_MVDR_RLS_sc2, w_MVDR_RLS_sc2, e_MVDR_RLS_sc2] = MVDRWeights_RLS(num_sensors, N, num_sources, ...
                d, d_lambda, r, 2, lamb_MVDR_RLS, AOA, signal_power, noise_power, TRIALS, del_MVDR);
[w_avg_MVDR_RLS_sc3, J_MVDR_RLS_sc3, w_MVDR_RLS_sc3, e_MVDR_RLS_sc3] = MVDRWeights_RLS(num_sensors, N, num_sources, ...
                d, d_lambda, r, 3, lamb_MVDR_RLS, AOA, signal_power, noise_power, TRIALS, del_MVDR);

% Compute Array responses
MVDR_Asource1_RLS = arrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc1(:, N), thetas, phis);
MVDR_Asource1_RLS_oneRun = arrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc1(:, N, 1), thetas, phis);

MVDR_Asource2_RLS = arrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc2(:, N), thetas, phis);
MVDR_Asource3_RLS = arrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc3(:, N), thetas, phis);

% Now we show the results of performing RLS for MVDR beamforming:
plotArrayResponse(MVDR_Asource1_RLS_oneRun, "Adaptive RLS MVDR - ONE RUN", "source 1", AOA, thetas, phis, 2);
plotArrayResponse(MVDR_Asource1_RLS, "Adaptive LMS MVDR", "source 1", AOA, thetas, phis, 2);

plotArrayResponse(MVDR_Asource2_RLS, "Adaptive LMS MVDR", "source 2", AOA, thetas, phis, 2);
plotArrayResponse(MVDR_Asource3_RLS, "Adaptive LMS MVDR", "source 3", AOA, thetas, phis, 2);

% Next plot the learning curves:
plotLearningCurve(abs(e_MVDR_RLS_sc1(:, :, 1)).^2, N, signal_power(1), "ONE-RUN RLS MVDR", "Source 1");
plotLearningCurve(J_MVDR_RLS_sc1, N, signal_power(1), "RLS MVDR", "Source 1");

plotLearningCurve(J_MVDR_RLS_sc2, N, signal_power(2), "RLS MVDR", "Source 2");
plotLearningCurve(J_MVDR_RLS_sc3, N, signal_power(3), "RLS MVDR", "Source 3");

% choose snapshot
N0_MVDR_RLS = 40;
[MVDR_Asource1_RLS_halfN0] = arrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc1(:, N0_MVDR_RLS/2, 1), ...
                                          thetas, phis);
[MVDR_Asource2_RLS_halfN0] = arrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc2(:, N0_MVDR_RLS/2, 1), ...
                                          thetas, phis);
[MVDR_Asource3_RLS_halfN0] = arrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc3(:, N0_MVDR_RLS/2, 1), ...
                                          thetas, phis);

% Examine the results for when the beam formers are steered to different
% sources with N0/2:
MVDR_RLS_SteerSource1_halfN0 = SourceShift(MVDR_Asource1_RLS_halfN0, AOA);
MVDR_RLS_SteerSource2_halfN0 = SourceShift(MVDR_Asource2_RLS_halfN0, AOA);
MVDR_RLS_SteerSource3_halfN0 = SourceShift(MVDR_Asource3_RLS_halfN0, AOA);

delim;
disp("Showing table for MVDR with RLS for N0 then N0/2:")
table_MVDR_RLS_halfN0 = table(Position, MVDR_RLS_SteerSource1_halfN0, MVDR_RLS_SteerSource2_halfN0, MVDR_RLS_SteerSource3_halfN0);
disp(table_MVDR_RLS_halfN0);
delim;

[MVDR_Asource1_RLS_N0] = arrayResponse(d, d_lambda,r, w_avg_MVDR_RLS_sc1(:, N0_MVDR_LMS, 1), ...
                                          thetas, phis);
[MVDR_Asource2_RLS_N0] = arrayResponse(d, d_lambda,r, w_avg_MVDR_RLS_sc2(:, N0_MVDR_LMS, 1), ...
                                          thetas, phis);
[MVDR_Asource3_RLS_N0] = arrayResponse(d, d_lambda,r, w_avg_MVDR_RLS_sc3(:, N0_MVDR_LMS, 1), ...
                                          thetas, phis);

% Examine the results for when the beam formers are steered to different
% sources with N0:
MVDR_RLS_SteerSource1_N0 = SourceShift(MVDR_Asource1_RLS_N0, AOA);
MVDR_RLS_SteerSource2_N0 = SourceShift(MVDR_Asource2_RLS_N0, AOA);
MVDR_RLS_SteerSource3_N0 = SourceShift(MVDR_Asource3_RLS_N0, AOA);

table_MVDR_RLS_N0 = table(Position, MVDR_RLS_SteerSource1_N0, MVDR_RLS_SteerSource2_N0, MVDR_RLS_SteerSource3_N0);
disp(table_MVDR_RLS_N0);
delim;

% Finally, compute the average array response for last 10 iterations
% N0-9 <= n <= N0
iterations = 10;
MVDR_Asource1_avg_RLS = avgArrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc1, ...
                                        100, iterations, thetas, phis);
MVDR_Asource2_avg_RLS = avgArrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc2, ...
                                        100, iterations, thetas, phis);
MVDR_Asource3_avg_RLS = avgArrayResponse(d, d_lambda, r, w_avg_MVDR_RLS_sc3, ...
                                        100, iterations, thetas, phis);

plotArrayResponse(MVDR_Asource1_avg_RLS, "Adaptive RLS MVDR - average weights", "source 1", AOA, thetas, phis, 2);

%% Adaptive Beamforming - GSC Using RLS
% Lastly, we will examine using RLS for the GSC beamformer and perhaps
% the results will be the best for these cases.

% GSC RLS PARAMETERS:
lamb_GSC_RLS = 0.95;
del_GSC = 50;

% Calculate the RLS weights for GSC (sources 1-3):
[w_avg_GSC_RLS_sc1, J_GSC_RLS_sc1, w_GSC_RLS_sc1, e_GSC_RLS_sc1] = GSCWeights_RLS(num_sensors, N, num_sources, ...
                d, d_lambda, r, gx, lamb_GSC_RLS, AOA, signal_power, noise_power, TRIALS, del_GSC);
[w_avg_GSC_RLS_sc2, J_GSC_RLS_sc2, w_GSC_RLS_sc2, e_GSC_RLS_sc2] = GSCWeights_RLS(num_sensors, N, num_sources, ...
                d, d_lambda, r, gy, lamb_GSC_RLS, AOA, signal_power, noise_power, TRIALS, del_GSC);
[w_avg_GSC_RLS_sc3, J_GSC_RLS_sc3, w_GSC_RLS_sc3, e_GSC_RLS_sc3] = GSCWeights_RLS(num_sensors, N, num_sources, ...
                d, d_lambda, r, gz, lamb_GSC_RLS, AOA, signal_power, noise_power, TRIALS, del_GSC);

% Compute and plot the array responses:
GSC_Asource1_RLS_oneRun = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc1(:, 1, 1), thetas, phis);
GSC_Asource1_RLS = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc1(:, N, 1), thetas, phis);

GSC_Asource2_RLS = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc2(:, N, 1), thetas, phis);
GSC_Asource3_RLS = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc3(:, N, 1), thetas, phis);


% Plot the GSC array response for a single run of RLS and for 100 runs of RLS:
plotArrayResponse(GSC_Asource1_RLS_oneRun, "Adaptive RLS GSC - ONE RUN", "source 1", AOA, thetas, phis, 2);
plotArrayResponse(GSC_Asource1_RLS, "Adaptive RLS GSC", "source 1", AOA, thetas, phis, 2);

plotArrayResponse(GSC_Asource2_RLS, "Adaptive RLS GSC", "source 2", AOA, thetas, phis, 2);
plotArrayResponse(GSC_Asource3_RLS, "Adaptive RLS GSC", "source 3", AOA, thetas, phis, 2);

% Now plot learning curves:
plotLearningCurve(abs(e_GSC_RLS_sc1(:, :, 1)).^2, N, signal_power(1), "ONE-RUN RLS GSC", "Source 1");
plotLearningCurve(J_GSC_RLS_sc1, N, signal_power(1), "RLS GSC", "Source 1");

plotLearningCurve(J_GSC_RLS_sc2, N, signal_power(2), "RLS GSC", "Source 2");
plotLearningCurve(J_GSC_RLS_sc3, N, signal_power(3), "RLS GSC", "Source 3");

% Let's take N0 = 40 here once again to be safe:
N0_GSC_RLS = 40;

% Compute the array response for each of the sources using N0/2:
[GSC_Asource1_RLS_halfN0] = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc1(:, N0_GSC_RLS/2, 1), ...
                                          thetas, phis);
[GSC_Asource2_RLS_halfN0] = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc2(:, N0_GSC_RLS/2, 1), ...
                                          thetas, phis);
[GSC_Asource3_RLS_halfN0] = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc3(:, N0_GSC_RLS/2, 1), ...
                                          thetas, phis);

% Now examine the steered array responses:
GSC_RLS_SteerSource1_halfN0 = SourceShift(GSC_Asource1_RLS_halfN0, AOA);
GSC_RLS_SteerSource2_halfN0 = SourceShift(GSC_Asource2_RLS_halfN0, AOA);
GSC_RLS_SteerSource3_halfN0 = SourceShift(GSC_Asource3_RLS_halfN0, AOA);

disp("Showing table for GSC with RLS for N0 then N0/2:")
table_GSC_RLS_halfN0 = table(Position, GSC_RLS_SteerSource1_halfN0, GSC_RLS_SteerSource2_halfN0, GSC_RLS_SteerSource3_halfN0);
disp(table_GSC_RLS_halfN0);
delim;

% Now N0:
[GSC_Asource1_RLS_N0] = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc1(:, N0_GSC_RLS, 1), ...
                                          thetas, phis);
[GSC_Asource2_RLS_N0] = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc2(:, N0_GSC_RLS, 1), ...
                                          thetas, phis);
[GSC_Asource3_RLS_N0] = arrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc3(:, N0_GSC_RLS, 1), ...
                                          thetas, phis);

GSC_RLS_SteerSource1_N0 = SourceShift(GSC_Asource1_RLS_N0, AOA);
GSC_RLS_SteerSource2_N0 = SourceShift(GSC_Asource2_RLS_N0, AOA);
GSC_RLS_SteerSource3_N0 = SourceShift(GSC_Asource3_RLS_N0, AOA);

table_GSC_RLS_N0 = table(Position, GSC_RLS_SteerSource1_N0, GSC_RLS_SteerSource2_N0, GSC_RLS_SteerSource3_N0);
disp(table_GSC_RLS_N0);
delim;

% The results from the table are once again good!

iterations = 10;
GSC_Asource1_avg_RLS = avgArrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc1, ...
                                        100, iterations, thetas, phis);
GSC_Asource2_avg_RLS = avgArrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc2, ...
                                        100, iterations, thetas, phis);
GSC_Asource3_avg_RLS = avgArrayResponse(d, d_lambda, r, w_avg_GSC_RLS_sc3, ...
                                        100, iterations, thetas, phis);

plotArrayResponse(GSC_Asource1_avg_RLS, "Adaptive RLS GSC - average weights", "source 1", AOA, thetas, phis, 2);



%% CONCLUDING REMARKS - FOR 10dB NOISE POWER ANALYSIS
% So the immediate effect I noticed with increasing the noise power to 10dB
% above the primary source is that MUSIC and MVDR do not do a good job
% picking up on source 2 and source 3. In the surface spectrum plots there
% is only one spike at what I would assume is the primary source which I
% believe is simply due to a lower SNR. The noise power becomes
% overwhelming to the point that it is hard to distinguish "meaningful"
% signals at each of the sources. This is also obvious when studying the
% eigenvalues of the correlation matrix and the singular values of the data
% matrix X. The ratio between the fourth and third singular values, and the
% fourth and third eigenvalues are significantly more large than what they
% were with noise power equal to -12dB below the primary. 
%
% Now, I will talk about the algorithms themselves and what I noticed
% between the two. NOTE: I WILL FIRST DESCRIBE WHAT I NOTICED WITHOUT
% TWEAKING ANY ALGORITHM PARAMETERS (such as step size, forgetting factor,
% etc.)
%
% I do not know if it is a placebo effect or not but right away I actually
% thought I noticed that MVDR was performing BETTER than how it was in the
% other parts. This is very interesting because I thought that there was
% almost no place for MVDR when you could just use GSC and call it a day
% (without considering potential computational costs and similar extraneous
% factors). It seemed as if the overall noise floor was perhaps slightly
% lower than the experiments in the previous section, while "amplifying"
% the array response at the particular source. Thus, it seems as if MVDR is
% more resilient to changes in the noise power.
%
% The GSC algorithms seemed to run basically about the same in terms of the
% array response performance with respect to attenuating the other sources
% far below the noise floor when they were not the "source of interest"
% and "turned on" the source that was of interest. 
%
% As far as learning curves were concerned, I noticed that the algorithms
% were not converging towards J_min and exhibited behavior seen for the
% experiments in the prior section where RLS was employed and divergence
% from J_min was observed. What was strange was that although the learning
% curves were diverging, the overall array response was not suggesting 
% signs of poorer performance due to this phenomenon.
%
% I will now tweak the algorithm parameters and note only the relevant
% changes...
%
% Okay, so changing the step size didn't seem to bring any benefit in terms
% of improving the rate of convergence to J_min for either learning curve.
% I also tried increasing the number of snapshots to N = 1000 in order to
% see if I was "too zoomed in" as in not giving enough snapshots to
% allow the algorithm to begin converging - this did not work either. I
% believe that when the noise power becomes greater than the primary source
% by this factor of 10dB, results can start to perform badly. Only when the
% noise power is below some threshold, are you able to see a benefit in
% terms of learning curve convergence.
%
% Furthermore, the spatial aliasing comment for the previous problem is
% also apparant in this experiment.







%% Functions Created - Copy pasted same functions as in other document
% Function to generate X data:
function [X, R_est, R_ideal, evalR_est, evalR_ideal, S] = generateX(sensors, N, d, d_lambda, r, ...
                                                                AOA, power_signal, power_noise)
    S = Svec(d, d_lambda, r, AOA);

    [A, a] = alphas(power_signal, N);

    % Generate noise vector V:
    [V, v] = generateNoise(sensors, N, power_noise);

    % Compute X from equation in slides:
    X = S*A + V;

    R_est = 1/N * (X*(ctranspose(X)));
    evalR_est = sort(eig(R_est), "descend");

    R_ideal = S * diag(a) * ctranspose(S) + diag(repmat(v, sensors, 1));
    evalR_ideal = sort(eig(R_ideal), "descend");
end

% Function to generate complex valued noise:
function [noise, noiseVariance] = generateNoise(sensors, N, power)
    % Determine variance of the nosie based on samples and noise power
    noiseVariance = (sensors^-1) * 10.^(power/10);
    
    % Generate the output noise:
    noise = (noiseVariance.^0.5) .* (randn(sensors, N) + 1j * randn(sensors, N)) / sqrt(2);

end


% Function to generate alpha complex baseband coefficients for the steering
% vectors:
function [alphas, alphaVariance] = alphas(power_signal, N)
    % Determine number of sources:
    L = size(power_signal, 2);
    alphaVariance = 10.^(power_signal / 10);

    alphas = ((alphaVariance.').^0.5) .* (randn(L, N) + 1j*(randn(L, N)))/sqrt(2);

end

% Function that creates the steering vectors using theta, phi, and d/lambda
% parameters:
function steering_vectors = Svec(d, d_lambda, r, AOA)
    % Convert AOA's to radians:
    AOA_radian = AOA*pi/180;
    
    % Determine lambda from d and d_lambda
    lambda = d/d_lambda;

    % Compute 3D unit vector
    K = (2*pi/lambda) * [sin(AOA_radian(1, :)).*cos(AOA_radian(2, :)) ; ...
                         sin(AOA_radian(1, :)).*sin(AOA_radian(2, :)) ; ...
                         cos(AOA_radian(1, :))];

    % Generate steering vectors
    steering_vectors = (1/sqrt(size(r, 2))) * exp(-1j*(r.') * K);

end

% function to compute MUSIC spectrum
function MUSIC_spectrum = MUSIC(d, d_lambda, Pn, r, phi, theta)

    MUSIC_spectrum = zeros(length(phi), length(theta));

    % COmpute steering vector based on AOA's and MUSIC algorithm
    for i = 1:1:length(phi)
        for j = 1:1:length(theta)
            S = Svec(d, d_lambda, r, [phi(1, i) ; theta(1, j)]);
            MUSIC_spectrum(i, j) = 1 / (S' * Pn * S);
        end
    end
end

% function to compute MVDR spectrum
function MVDR_spectrum = MVDR(d, d_lambda, R, r, phi, theta)

    MVDR_spectrum = zeros(length(phi), length(theta));

    % COmpute steering vector based on AOA's and MUSIC algorithm
    for i = 1:1:length(phi)
        for j = 1:1:length(theta)
            S = Svec(d, d_lambda, r, [phi(1, i) ; theta(1, j)]);
            MVDR_spectrum(i, j) = 1 / (S' * (R^-1) * S);
        end
    end
end

% Function to create GSC weights 
function W_GSC = GSCWeights(num_sensors, num_sources, g, S, R)
    % Construct the M x L steering vector matrix C in the notes along with
    % its quiescent vector:
    C = S;
    quiescent_vector = C * (C' * C)^-1;

    % Calculate vector w that satisfies affine constraints:
    w_q = quiescent_vector * g;

    % Create the Ca matrix whose columns span the orthogonal complement to
    % span{S}:
    [U, ~, ~] = svd(C); 
    Ca = U(:, (num_sources + 1):end);

    % Calculate the ideal weights for the GSC based on equation in notes:
    W_GSC = (eye(num_sensors) - Ca * (inv(Ca'*R*Ca)) * Ca' * R) * w_q;
end

% Function to create the MVDR weights:
function W_MVDR = MVDRWeights(S, source_direction, R)
    % Determine the steering vector for the given source 
    Svec = S(:, source_direction);

    % Calculate the weights for the MVDR beamformer from the equation in
    % notes:
    W_MVDR = (1 / (Svec'*(R^-1)*Svec)) * (R^-1) * Svec;
    
end

% Function that computes the array response at a given theta and phi AOA:
function A_THETA = arrayResponse(d, d_lambda, r, w, theta, phi)
    % Initialize output of array repsonse:
    A_THETA = zeros(length(theta), length(phi));

    for i = 1:length(theta)
        for j = 1:length(phi)
            % Generate steering vectors for each AOA
            S = Svec(d, d_lambda, r, [theta(1, i) ; phi(1, j)]);
            A_THETA(i, j) = abs(w' * S)^2;
        end
    end
end

% Function to average the array response
function A_THETA_avg = avgArrayResponse(d, d_lambda, r, w, N0, iterations, thetas, phis)
    % Initialize output of average array response:
    A_THETA_avg = zeros(length(thetas), length(phis));

    for i = 1:iterations
        A_THETA_avg = A_THETA_avg + arrayResponse(d, d_lambda, r, w(:, (N0 - iterations + 1), 1), thetas, phis);
    end
    % Average the array response:
    A_THETA_avg = A_THETA_avg / iterations;
end

%% Plotting functions and ease of use functions
% Function to plot dB array response surface plot:
function plotArrayResponse(array_response, NAME, SOURCE_NUM, AOA, theta, phi, dB_linear_FLAG)
    % if dB_linear_FLAG = 0:    plot just the dB plot
    % if dB_linear_FLAG = 1:    plot the linear plot only
    % if dB_linear_FLAG = 2:    plot the dB and linear plot
    % srcLocation:              The theta and phi location of the sources
    %                           theta in position 1, phi in position 2

    % Compute dB values of array response:
    array_response_dB = 10*log10(abs(array_response));

    if(dB_linear_FLAG == 0)

        figure;
        surf(phi, theta, array_response_dB);
        title(NAME + " [dB] array response A(\Theta) | " + SOURCE_NUM);
        xlabel("\phi[deg]");
        ylabel("\theta [deg]");
        zlabel("Magnitude [dB]");
        xlim([-180 180]);
        ylim([0 90]);
        
        % Add red dots for 3 source locations (dB)
        hold on;
        scatter3(AOA(2, 1), AOA(1, 1), array_response_dB(find(AOA(1, 1) == theta), find(AOA(2, 1) == phi)), 'r', 'filled');
        scatter3(AOA(2, 2), AOA(1, 2), array_response_dB(find(AOA(1, 2) == theta), find(AOA(2, 2) == phi)), 'r', 'filled');
        scatter3(AOA(2, 3), AOA(1, 3), array_response_dB(find(AOA(1, 3) == theta), find(AOA(2, 3) == phi)), 'r', 'filled');
        hold off;
    end

    if(dB_linear_FLAG == 1)
        figure;
        surf(phi, theta, array_response);
        title(NAME + " [linear] array response A(\Theta) | " + SOURCE_NUM);
        xlabel("\phi[deg]");
        ylabel("\theta [deg]");
        zlabel("Magnitude [linear]");
        xlim([-180 180]);
        ylim([0 90]);

        % Add red dots for source locations (LINEAR)
        hold on;
        scatter3(AOA(2, 1), AOA(1, 1), array_response(find(AOA(1, 1) == theta), find(AOA(2, 1) == phi)), 'r', 'filled');
        scatter3(AOA(2, 2), AOA(1, 2), array_response(find(AOA(1, 2) == theta), find(AOA(2, 2) == phi)), 'r', 'filled');
        scatter3(AOA(2, 3), AOA(1, 3), array_response(find(AOA(1, 3) == theta), find(AOA(2, 3) == phi)), 'r', 'filled');
        hold off;
    end
    
    if(dB_linear_FLAG == 2)
        % Compute dB values of array response:
        array_response_dB = 10*log10(abs(array_response));

        figure;
        subplot(1, 2, 1);
        surf(phi, theta, array_response);
        title(NAME + " [linear] array response A(\Theta) | " + SOURCE_NUM);
        xlabel("\phi[deg]");
        ylabel("\theta [deg]");
        zlabel("Magnitude [linear]");
        xlim([-180 180]);
        ylim([0 90]);

        % Add red dots for source locations (LINEAR)
        hold on;
        scatter3(AOA(2, 1), AOA(1, 1), array_response(find(AOA(1, 1) == theta), find(AOA(2, 1) == phi)), 'r', 'filled');
        scatter3(AOA(2, 2), AOA(1, 2), array_response(find(AOA(1, 2) == theta), find(AOA(2, 2) == phi)), 'r', 'filled');
        scatter3(AOA(2, 3), AOA(1, 3), array_response(find(AOA(1, 3) == theta), find(AOA(2, 3) == phi)), 'r', 'filled');
        hold off;

        subplot(1, 2, 2);
        surf(phi, theta, array_response_dB);
        title(NAME + " [dB] array response A(\Theta) | " + SOURCE_NUM);
        xlabel("\phi[deg]");
        ylabel("\theta [deg]");
        zlabel("Magnitude [dB]");
        xlim([-180 180]);
        ylim([0 90]);
        
        % Add red dots for 3 source locations (dB)
        hold on;
        scatter3(AOA(2, 1), AOA(1, 1), array_response_dB(find(AOA(1, 1) == theta), find(AOA(2, 1) == phi)), 'r', 'filled');
        scatter3(AOA(2, 2), AOA(1, 2), array_response_dB(find(AOA(1, 2) == theta), find(AOA(2, 2) == phi)), 'r', 'filled');
        scatter3(AOA(2, 3), AOA(1, 3), array_response_dB(find(AOA(1, 3) == theta), find(AOA(2, 3) == phi)), 'r', 'filled');
        hold off;
    end
end

% Function to plot 2-D Spectrums
function plotSpectrum(Spectrum, thetas, phis, NAME)
    % Create a figure of the desired MVDR or MUSIC spectrum 
    figure;
    % Plot in dB magnitude ** 10log10 because of power **
    surf(thetas, phis, 10*log10(abs(Spectrum)'));
    title(NAME + " Spectrum of the data");
    % Let x axis plot phis and y axis plot thetas
    xlabel("\phi");
    ylabel("\theta");
    zlabel("Magnitude [dB]");
    xlim([min(thetas) max(thetas)]);
    ylim([min(phis) max(phis)]);
    shading("flat");
    colormap("turbo");
end

% Function to plot 1D spectrums for prescribed algorithm
function plot1DSpectrum(Spectrum, theta_or_phi, SPECTRUM_NAME, CONSTANT_NAME, ...
                        VARYING_VAR, CONSTANT_VAL, LINE_VAL)
    % For MAX_MIN_FLAG = 0:      plot the max of the spectrum as a vertical
    %                            line
    % For MAX_MIN_FLAG = 1:      plot the min of the spectrum as a vertical
    %                            line

    % Create a figure of the desired MVDR or MUSIC spectrum 
    figure;
    % Plot in dB magnitude ** 10log10 because of power **
    Spectrum_dB = 10*log10(abs(Spectrum));
    plot(theta_or_phi, Spectrum_dB);
    title(SPECTRUM_NAME + " Spectrum of the data at \" + CONSTANT_NAME + " = " + num2str(CONSTANT_VAL));
    xlabel("\"+VARYING_VAR+" [deg]");
    ylabel("Magnitude [dB]");
    xline(LINE_VAL, 'r--', "\"+VARYING_VAR+" = "+LINE_VAL);
    xlim([min(theta_or_phi) max(theta_or_phi)]);
end

% Function to plot arbitrary learning curves:
function plotLearningCurve(J, N, power_signal, NAME, SOURCE)
    % Note: power_signal will be used to place J_min line:
    J_min = 10^(power_signal / 10);         % J_min on a LINEAR scale

    figure;
    plot(1:N, J);
    title("Learning curve for " + NAME + " | " + SOURCE);
    xlabel("Snapshot");
    yline(J_min, "--b");
    legend("Learning curve (J)", "J_{min} = " + J_min);
end

% Function to compute shifted array responses for each beamformer so that
% it can be displayed nicely in a table and make my life easier:
function Shifted_Source = SourceShift(array_response, AOA)
    Shifted_Source = [array_response(AOA(1, 1) + 1, AOA(2, 1) + 181); ...
                      array_response(AOA(1, 2) + 1, AOA(2, 2) + 181); ...
                      array_response(AOA(1, 3) + 1, AOA(2, 3) + 181)];
end

% Function to display a page of #'s to separate tables:
function delim 
    disp("##################################################################");
end

%% Adaptive Algorithms - Functions
% Function to perform LMS adaptive weight generation for GSC
function [w_average, J, w, error] = GSCWeights_LMS(num_sensors, N, num_sources, d, d_lambda, r, ...
                                                    g, AOA, power_signal, power_noise, RUNS, mu)
%
% Note that parameters RUNS indicates the number of runs of LMS that we
% wish to perform for GSC weight estimation. FUrthermore, the mu
% parameter which will govern rate of convergence.
%
% OUTPUTS:
%   J - Learning curve
%   w - weights => w_average are the average weights
%   error - estimated error after each run of LMS

    % Preallocation:
    w_a = zeros(num_sensors - num_sources, N, RUNS);
    w = zeros(num_sensors, N, RUNS);
    error = zeros(1, N, RUNS);

    % Run LMS for #RUNS - experiments
    for i = 1:RUNS
        % Compute X data and the steering vectors for each trial:
        [u, ~, ~, ~, ~, S] = generateX(num_sensors, N, d, d_lambda, r, AOA, power_signal, power_noise);

        % Like the non-adaptive GSC weight adjustment, compute the C
        % matrices based on the steering vectors so that adaptive
        % adjustment can be performed:
        C = S;
        [U, ~, ~] = svd(C);
        Ca = U(:, (num_sources + 1):end);
        quiescent_vector = C * (C' * C)^-1;
        % Adjust the weights from the parameters:
        w_q = quiescent_vector * g;

        % Now iterate and perform LMS for N iterations since the number of
        % snapshots determines the total runs of LMS:
        %       Algorithm given in the notes        %
        for n = 1:N
            % Iterate through LMS algorithm
            d_n = w_q' * u(:, n);
            x_n = Ca' * u(:, n);
            error(1, n, i) = d_n - (w_a(:, n, i)' * x_n);
            w_a(:, n + 1, i) = w_a(:, n, i) + mu * x_n * error(1, n, i)';
        end
        % Compute the weight at the end of N iterations of LMS:
        w(:, :, i) = w_q - (Ca * w_a(:, 1:N, i));
    end

    % Compute the learning curve via the mean square error:
    J = mean(abs(error).^2, 3);

    % Finally, compute the average mean across all runs:
    w_average = mean(w, 3);
end

% Function to compute LMS weights for MVDR
function [w_average, J, w, error] = MVDRWeights_LMS(num_sensors, N, num_sources, d, d_lambda, r, ...
                                                    source_number, AOA, power_signal, power_noise, RUNS, mu)
%
% Note that this code will be very similar to LMS for GSC, hence I am
% copying it and adjusting the code inside here
%
% OUTPUTS:
%   J - Learning curve
%   w - weights => w_average are the average weights
%   error - estimated error after each run of LMS

    % Preallocation:
    w_a = zeros(num_sensors - num_sources, N, RUNS);
    w = zeros(num_sensors, N, RUNS);
    error = zeros(1, N, RUNS);

    % Run LMS for #RUNS - experiments
    for i = 1:RUNS
        % Compute X data and the steering vectors for each trial:
        [u, ~, ~, ~, ~, S] = generateX(num_sensors, N, d, d_lambda, r, AOA, power_signal, power_noise);
        % List constraints based on steering vectors
        C = S(:, source_number);
        % Compute the weights:
        [U, ~, ~] = svd(C);
        Ca = U(:, (num_sources + 1):end);
        quiescent_vector = C * (C' * C)^-1;
        % Adjust the weight vector from GSC by removing 'g'
        w_q = quiescent_vector;

        % Now iterate and perform LMS for N iterations since the number of
        % snapshots determines the total runs of LMS:
        %       Algorithm given in the notes        %
        for n = 1:N
            % Iterate through LMS algorithm
            d_n = w_q' * u(:, n);
            x_n = Ca' * u(:, n);
            error(1, n, i) = d_n - (w_a(:, n, i)' * x_n);
            w_a(:, n + 1, i) = w_a(:, n, i) + mu * x_n * error(1, n, i)';
        end
        % Compute the weight at the end of N iterations of LMS:
        w(:, :, i) = w_q - (Ca * w_a(:, 1:N, i));
    end

    % Compute the learning curve via the mean square error:
    J = mean(abs(error).^2, 3);

    % Finally, compute the average mean across all runs:
    w_average = mean(w, 3);
end

% Function to compute weights using RLS for MVDR beamforming
function [w_average, J, w, error] = MVDRWeights_RLS(num_sensors, N, num_sources, d, d_lambda, ...
                            r, source_number, lamb, AOA, power_signal, power_noise, RUNS, del)
    
% Note that this code will be very similar to LMS for MVDR, hence I am
% copying it and adjusting the code inside here
%
% INPUTS: 
%   del - forgetting factor
% OUTPUTS:
%   J - Learning curve
%   w - weights => w_average are the average weights
%   error - estimated error after each run of LMS

    % Preallocation:
    w_a = zeros(num_sensors - num_sources, N, RUNS);
    w = zeros(num_sensors, N, RUNS);
    error = zeros(1, N, RUNS);
    
    % Run RLS for #RUNS - experiments
    for i = 1:RUNS
        % Compute X data and the steering vectors for each trial:
        [u, ~, ~, ~, ~, S] = generateX(num_sensors, N, d, d_lambda, r, AOA, power_signal, power_noise);
        % List constraints based on steering vectors
        C = S(:, source_number);
        % Compute the weights:
        [U, ~, ~] = svd(C);
        Ca = U(:, (num_sources + 1):end);
        quiescent_vector = C * (C' * C)^-1;
        % Adjust the weight vector from GSC by removing 'g'
        w_q = quiescent_vector;
        
        % Construct P_N vector for RLS
        P_n = del^-1 * eye(num_sensors - num_sources);

        % Now iterate and perform RLS for N iterations since the number of
        % snapshots determines the total runs of LMS:
        %       Algorithm given in the notes        %
        % The difference between LMS and RLS is the introduction of the
        % forgetting factor essentially
        for n = 1:N
            % Iterate through RLS algorithm
            d_n = w_q' * u(:, n);
            x_n = Ca' * u(:, n);
            pi_n = P_n * x_n;
            k_n = (lamb + x_n' * pi_n)^-1 * pi_n;
            P_n = lamb^-1 * P_n - lamb^-1 * k_n * x_n' * P_n;
            error(1, n, i) = d_n - (w_a(:, n, i)' * x_n);
            w_a(:, n + 1, i) = w_a(:, n, i) + k_n * error(1, n, i)';
        end
        % Compute the weight at the end of N iterations of LMS:
        w(:, :, i) = w_q - (Ca * w_a(:, 1:N, i));
    end

    % Compute the learning curve via the mean square error:
    J = mean(abs(error).^2, 3);

    % Finally, compute the average mean across all runs:
    w_average = mean(w, 3);

end

% Finally, one more function to perform RLS for GSC since the updated
% weights are dependent on the g vector
function [w_average, J, w, error] = GSCWeights_RLS(num_sensors, N, num_sources, d, d_lambda, r, ...
                                    g, lamb, AOA, power_signal, power_noise, RUNS, del)
    % Preallocation:
    w_a = zeros(num_sensors - num_sources, N, RUNS);
    w = zeros(num_sensors, N, RUNS);
    error = zeros(1, N, RUNS);

    % Run LMS for #RUNS - experiments
    for i = 1:RUNS
        % Compute X data and the steering vectors for each trial:
        [u, ~, ~, ~, ~, S] = generateX(num_sensors, N, d, d_lambda, r, AOA, power_signal, power_noise);

        % Like the non-adaptive GSC weight adjustment, compute the C
        % matrices based on the steering vectors so that adaptive
        % adjustment can be performed:
        C = S;
        [U, ~, ~] = svd(C);
        Ca = U(:, (num_sources + 1):end);
        quiescent_vector = C * (C' * C)^-1;
        % Adjust the weights from the parameters:
        w_q = quiescent_vector * g;
        
        % Compute P_n matrix:
        P_n = del^-1 * eye(num_sensors - num_sources);

        % Now iterate and perform RLS for N iterations since the number of
        % snapshots determines the total runs of RLS:
        %       Algorithm given in the notes        %
        for n = 1:N
            % Iterate through LMS algorithm
            d_n = w_q' * u(:, n);
            x_n = Ca' * u(:, n);
            pi_n = P_n * x_n;
            k_n = (lamb + x_n' * pi_n)^-1 * pi_n;
            P_n = lamb^-1 * P_n - lamb^-1 * k_n * x_n' * P_n;
            error(1, n, i) = d_n - (w_a(:, n, i)' * x_n);
            w_a(:, n + 1, i) = w_a(:, n, i) + k_n * error(1, n, i)';
        end
        % Compute the weight at the end of N iterations of LMS:
        w(:, :, i) = w_q - (Ca * w_a(:, 1:N, i));
    end

    % Compute the learning curve via the mean square error:
    J = mean(abs(error).^2, 3);

    % Finally, compute the average mean across all runs:
    w_average = mean(w, 3);
end