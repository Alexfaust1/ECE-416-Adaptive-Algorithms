%% Alexander Faust
%
% ECE 416 - Problem Set 5 - Unscented Kalman Filter
%
% December 12, 2023
clc; clear; close all;

% List Physical Parameters:
Gm0 = 3.9986*10^5;
Y0 = 6374;
H0 = 13.406;

fs = 10;                % Sample rate of signals [Hz]
dt = 0.1;               % Time step between samples [sec]
%t = n*dt;               % Time variable for discrete formula

% Initial conditions:
iterations = 500;
x0 = [ 6400.4 ;
       349.14 ;
      -1.8093 ;
      -6.7967 ;
       0.6932 ];

%% Theory and Setup
% Computing the Kurtosis
% Referring to slides for the parameters and equations:
w0 = 1/3;
Nx = 5;
M = (Nx - 1)/2;
w_tild = (1 - w0)/(2 * M);
X0 = zeros(1, 2*M);
X00 = zeros(M, 1);

% Iterate over M
for i = 1:M
    X0(i) = 1/sqrt(2*w_tild);
    X0(M + i) = -1/sqrt(2*w_tild);
end

X0 = [X00, diag(X0(1:M)), diag(X0(M + 1 : end))];
sigmas = X0;

w = [w0, repmat(w_tild, 1, 2*M)];

% Compute the empirical kurtosis for each component:
Kemp = zeros(1, M);
for j = 1:M
    Kemp(j) = sum(w .* (X0(j, :).^4));
end

% Display and check if components are 3 as in the assignment:
disp(Kemp);

% Components are both 3 meaning we are good to go

%% Wrapper Unit
% Run the unit 5 times via this wrapper
num_wraps = 5;

for k = 1:num_wraps
    % List disturbance matrices:
    Qx = diag([10^-8, 10^-9, 2.404*10^-5, 2.404*10^-5, 10^-8]);
    Qy = diag([1, 17*10^-3]);
    
    % List UKF initial condition:
    x0_estimate = [ 6400 ;
                     350 ;
                     -2  ;
                     -7  ;
                    0.65 ];
    % List initial condition for Kalman Filter:
    % *** WHEN I RAN IT AS IN THE HOMEWORK IT GOT A POSITIVE DEFINITE ERROR 
    % So I had to change this last row of the initial K0
    K_1_0 = diag([1e-4, 1e-4, 1e-4, 1e-4, 1e-5]);
    
    % Initialize trajectory vectors for the 5 different runs
    X_actual = zeros(5, iterations);
    X_actual(:, 1) = x0;
    X_estimate = zeros(5, iterations);
    X_predicted = zeros(5, iterations);
    y_measured = zeros(2, iterations);
    
    
    % Perform the UKF algorithm:
    y_measured(:, 1) = uhlmeas(X_actual(:, 1), chol(Qy));
    [X_predicted(:, 1), Kp, X_estimate(:, 1), Ke] = ukfilter(x0_estimate, K_1_0, Qx, Qy, y_measured(:, 1), dt, dt);
    
    % Iterate over 500 time steps so as to simulate 50 seconds for the
    % simulation:
    for t = 2:iterations
        % Simulate the actual trajectory using the uhlprocsim function:
        X_actual(:, t) = uhlprocsim(X_actual(:, t - 1), t*dt, dt, 'm', chol(Qx));
        % Measure y from the actual trajectory process:
        y_measured(:, t) = uhlmeas(X_actual(:, t), chol(Qy));
        % Pass the measured y's and the predicted X's to the UKF:
        [X_predicted(:, t), Kp, X_estimate(:, t), Ke] = ukfilter(X_predicted(:, t - 1), Kp, Qx, Qy, y_measured(:, t), t*dt, dt);
    end
    
    % Create subplots for trajectories:
    figure;
    subplot(3, 1, 1)
    % Plot the actual trajectory of x between state 1 and state 2:
    plot(X_actual(1, :), X_actual(2, :), '-X', 'MarkerIndices', [1], 'MarkerFaceColor', ...
        'blue', 'MarkerSize', 10);
    hold on
    % Plot the estimate of x between state 1 and state 2:
    plot(X_estimate(1, :), X_estimate(2, :), '-X', 'MarkerIndices', [1], 'MarkerFaceColor', ...
        'blue', 'MarkerSize', 10);
    % Plot the predicted x between state 1 and state 2:
    plot(X_predicted(1, :), X_predicted(2, :), '-X', 'MarkerIndices', [1], 'MarkerFaceColor', ...
        'blue', 'MarkerSize', 10);
    legend(["Actual Trajectory", "Estimated Trajectory", "Predicted Trajectory"], "Location", "northwest");
    hold off
    
    subplot(3, 1, 2);
    % Plot the actual trajectory of x between state 3 and state 4 now:
    plot(X_actual(3, :), X_actual(4, :), '-X', 'MarkerIndices', [1], 'MarkerFaceColor', ...
        'blue', 'MarkerSize', 10);
    hold on
    % Plot the estimate of x between state 1 and state 2:
    plot(X_estimate(3, :), X_estimate(4, :), '-X', 'MarkerIndices', [1], 'MarkerFaceColor', ...
        'blue', 'MarkerSize', 10);
    % Plot the predicted x between state 1 and state 2:
    plot(X_predicted(3, :), X_predicted(4, :), '-X', 'MarkerIndices', [1], 'MarkerFaceColor', ...
        'blue', 'MarkerSize', 10);
    legend(["Actual Trajectory", "Estimated Trajectory", "Predicted Trajectory"], "Location", "northwest");
    hold off
    
    beta0 = 0.597983;
    beta_estimate = beta0 * exp(X_estimate(5, :));
    beta_predicted = beta0 * exp(X_predicted(5, :));
    beta_actual = beta0 * exp(X_actual(5, :));
    
    % Plot actual beta and the estimated value with the x_5 state variable with
    % superimposed plots.
    subplot(3, 1, 3);
    % Plot the actual trajectory of x between state 1 and state 2:
    plot(beta_actual, '-X', 'MarkerIndices', [1], 'MarkerFaceColor', 'blue', 'MarkerSize', 20);
    hold on
    % Plot the estimate of x between state 1 and state 2:
    plot(beta_predicted, '-X', 'MarkerIndices', [1], 'MarkerFaceColor', 'blue', 'MarkerSize', 20);
    % Plot the predicted x between state 1 and state 2:
    plot(beta_estimate, '-X', 'MarkerIndices', [1], 'MarkerFaceColor', 'blue', 'MarkerSize', 20);
    legend(["Actual Trajectory", "Estimated Trajectory", "Predicted Trajectory"], "Location", "northwest");
    hold off

end

%% Comment Section

% After running the UKF once I proceeded to place it in a wrapper unit and
% ran the experiment 5 times. Across all 5 tests the results were fairly
% the same which implies that UKF does a pretty good job with trajectory
% tracking. 
%
% Note that the closer the initial predicted value of X is to the ground
% truth, the better the tracking will perform. The one irregularity was in
% the graph with the ballistic coefficient. It looked like there was some
% variation between the predicted and estimated trajectory. However, the
% y-axes units were fairly close together so perhaps not as big of a
% deviation in the grand scheme of things. 





%% Functions Created
% Function to perform the unscented Kalman filter:
function [xpnew, Kpnew, xe, Ke] = ukfilter(xp, Kp, Qx, Qy, ymeas, n, dt)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:
    %
    % xp    - Initial predicted state
    % Ke    - Initial predicted error covariance matrix
    % Qx    - Measurement covariance matrix for x
    % Qy    - Measurement covariance matrix for y
    % ymeas - Measured vector 
    % n     - Discrete iteration index that determines simulation time
    % dt    - Discretization parameter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:
    %
    % xe    - Estimated state
    % Ke    - Associated covariance matrix
    % xpnew - Updated predicted value for x
    % Kpnew - Updated predicted value for K
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % According to notes we should do 8 steps for the UKF algorithm:
    % 1. Compute predicted sigma points
    % 2. Apply UT for the measurement equation to obtain predicted
    %    measurements
    % 3. Compute mean predicted measurements
    % 4. Compute cross-covariance between predicted state and predicted
    %    measurements, and covariance of predicted measurements
    % 5. Apply the Kalman correction equations
    % 6. Compute filtered sigma points
    % 7. Apply UT for the process equation to obtain predicted states
    % 8. Compute mean and covariance of predicted state

    % Step 1.
    [wp, sig_xp] = sigmapoints(xp, chol(Kp));

    % Step 2.
    sig_yp = zeros(2, size(sig_xp, 2));
    lengthX = size(sig_xp, 2);

    % Iterate over X to apply the uhlmeas function to obtain the predicted
    % measurements
    for i = 1:lengthX
        sig_yp(:, i) = uhlmeas(sig_xp(:, i), chol(Qy));
    end

    % Step 3.
    yp = sigmamean(wp, sig_yp);

    % Step 4.
    Kyyp = sigmacov(wp, sig_yp, sig_yp) + Qy;
    Kxyp = sigmacov(wp, sig_xp, sig_yp);

    % Step 5.
    alpha = ymeas - yp;
    G =  Kxyp * Kyyp^-1;
    xe = xp + G * alpha;
    Ke = Kp - G * Kyyp * G';

    % Step 6.
    [w, sig_xe] = sigmapoints(xe, chol(Ke));

    % Step 7. 
    new_sig_xp = zeros(size(sig_xe));
    for i = 1:size(sig_xe, 2)
        new_sig_xp(:, i) = uhlprocsim(sig_xe(:, i), n, dt, 'm', chol(Qx));
    end

    % Step 8.
    xpnew = sigmamean(w, new_sig_xp);
    Kpnew = sigmacov(w, new_sig_xp, new_sig_xp) + Qx;


end

% Function to compute sigma points and weights from the Kurtosis of each
% component
function [w, sig] = sigmapoints(mu, Cchol)
    Nx = size(Cchol, 2);
    sig = zeros(Nx, 2*Nx + 1);
    scale = sqrt(Nx / (1 - 1/3));

    % Scale the sigmas accordingly:
    sig(:, 2 : (Nx + 1)) = scale * eye(Nx);
    sig(:, Nx + 2 : end) = -scale * eye(Nx);

    % Compute the weights:
    w = zeros(1, 2*Nx + 1);
    w(:, 1) = 1/3;
    w(:, 2 : end) = (1 - (1/3)) / (2 * Nx);

    sig = Cchol * sig + mu;
    
end


% Function to compute the empirical mean from the sigmas and weights:
function sigmu = sigmamean(w, sig)
    % Just take a simple mean across the lin comb of sigmas and weights:
    sigmu = sum(sig .* w, 2);

end

% Function to compute the empirical cross-covariance for X, Y as given by
% sigma points sigx, sigy with the common weights w.
function sigcov = sigmacov(w, sigx, sigy)
    % Use sigmamean function to compute the mean of each sigma across all
    % tbe weights to pack into the covariance matrix
    sigxMean = sigmamean(w, sigx);
    sigyMean = sigmamean(w, sigy);
    
    % Compute the sum
    summation = w(1) .* ((sigx(:, 1) - sigxMean) * (sigy(:, 1) - sigyMean)');
    % Over all the n's compute the sum
    lengthW = size(w, 2);

    for n = 2:lengthW
        summation = summation + w(n) .* ((sigx(:, n) - sigxMean) * (sigy(:, n) - sigyMean)');
    end

    sigcov = summation;

end












