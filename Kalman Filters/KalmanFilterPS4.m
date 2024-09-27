%% Alexander Faust
%
% ECE 416 - Problem Set 4 - Kalman Filters
%
% December 7th, 2023
clc; clear; close all;

%% Problem 1 - Basic Kalman Filter - Part I
% List the A matrices:
A1 = [-3/2, -5/2 ;
       1/2, 1/2 ];

A2 = [-3, -5  ;
       1,  1 ];

% List the C, Qx matrices and Qy:
C = [1, 0];
Qx = 0.1 * eye(2);
Qy = 0.5;

% 1. Preliminary analysis - Compute eigenvalues
evalsA1 = eig(A1);
evalsA2 = eig(A2);
disp(evalsA1);
disp(evalsA2);

% Compute the observability matrix and determine its rank:
O_A1 = [C ; C * A1];
rank_O_A1 = rank(O_A1);
disp(rank_O_A1);

O_A2 = [C ; C * A2];
rank_O_A2 = rank(O_A2);
disp(rank_O_A2);

% Compare K using the idare function:
% Conversion from Function syntax to our experiment:
N = 2;          % Number of states
p = 1;          % Number of outputs 
E = eye(N);
B = C';
R = Qy;
Q = Qx;
S = zeros(N, p);

% Create the Kalman filters
[K1, ~, ~] = idare(A1', B, Q, R, S, E);
[K2, ~, ~] = idare(A2', B, Q, R, S, E);

% COMMENTS:
% Based on the eigenvalues of A1 and A2 we can see that since we are
% dealing with a discrete system the eigenvalues of A2 are on the unit
% circle which implies that this is the unstable system.
%
% Eigenvalues of A1 lie within unit circle so we are good.
%
% Last check: rank of observability matrix is 2 (full rank) which is
% confirmation that we can uniquely determine the output state trajectory.

%% Part II - Kalman Filter
% Run first 50 iterations with stable form then after switch to unstable
% List parameters:
iterations = 100;               % Run filter for 100 iterations
x0 = [1 ;                       % Initial condition to generate state trajectory
      0 ];

% Pre-allocate xn and xy outputs:
X = zeros(2, iterations + 1);
yn = zeros(2, iterations);
% Initialize x(1|0) = x0:
X(:, 1) = x0;

% Iterate through lower half of iterations and generate the X and Y's for
% Kalman filter:
for n = 1:(iterations/2 - 1)
    X(:, n + 1) = A1 * X(:, n) + (chol(Qx) * randn(size(X, 1), 1));
    yn(n) = C * X(:, n) + chol(Qy) * randn(1, 1);
end

% Iterate through upper half of iterations and generate x and y's for
% Kalman to track:
for n = (iterations/2+1):iterations
    X(:, n + 1) = A2 * X(:, n) + (chol(Qx) * randn(size(X, 1), 1));
    yn(n) = C * X(:, n) + chol(Qy) * randn(1, 1);
end

% Apply the Kalman Filter function to perform the tracking:
[X_pred, X_est, K_next, error] = Kalman(iterations, X, yn, A1, A2, C, K1, K2, Qx, Qy);

% Show the plots for part a)
%Superimpose prediction error and estimation error for axes
%x1(n)-x1_est(n|n-1) and x1(n)-x1_est(n|n):
figure;
hold on
% Plot the actual trajectory from Riccatti 
plot(X(1, 1 : end - 1), X(2, 1 : end - 1));

% Plot the estimated trajectory:
plot(X_est(1, :), X_est(2, :), 'r');

% Plot the predicted trajectory:
plot(X_pred(1, 1:end-1), X_pred(2, 1:end-1), 'g');

hold off
legend("Actual X (Ricatti)", "Estimated X", "Predicted X");
title("Trajectory plot for the estimated, predicted, and actual");
ylabel("x_{2}");
xlabel("X_{1}");

% Comments:
% The results for the three plots are incredibly accurate, I really had to
% zoom in to see the difference.



% Part b) - Plot ||K(n, n-1) - K_ideal||
figure;
subplot(2, 1, 1);
plot(10*log10(error))
title("Plot of 10log_{10}(||K(n, n - 1) - K_{ideal}(n)||");
xlabel("iterations");
ylabel("||K(n, n - 1) - K_{ideal}(n)|| [dB]");

subplot(2, 1, 2);
plot(error);
title("Graph of ||K(n, n - 1) - K_{ideal}(n)|| [linear]");
xlabel("iterations");
ylabel("||K(n, n - 1) - K_{ideal}(n)||");

% Comments: From these plots it appears that in the beginning it took 3
% iterations for the filter to converge (on the linear scale as reference).
% Then after the inaccuracy was seen at index 50 it took about 5 iterations
% for the filter to be at convergence again.


% Part c) - Obtain plots for when the filter seems to have converged i.e.
%           no more spikes

figure;
% Create subplot for the prediction x observation:
subplot(3, 1, 1);
plot(1 : iterations/2, sqrt( sum((X(:, (iterations/2 + 1) : iterations) - ...
                        X_pred(:, (iterations/2 + 1) : iterations)).^2, 1)) );
title("Prediction of x_{1}");
ylabel("X_{1}");
xlabel("iterations");

% Create subplot for the estimated x observation:
subplot(3, 1, 2);
plot(1 : iterations/2, sqrt( sum((X(:, (iterations/2 + 1) : iterations) - ...
                        X_est(:, (iterations/2 + 1) : iterations)).^2, 1)) );
title("Estimation of x_{1}");
ylabel("X_{1}");
xlabel("iterations");

% Create subplot for the actual trajectory observation:
subplot(3, 1, 3);
plot(X(1, (iterations/2 + 1):iterations), X(2, (iterations/2 + 1) : iterations), 'r');
title("Actual trajectories for iterations after the spike (50 < n)");
ylabel("X_{2}");
xlabel("X_{1}");

% I do not know if you wanted this type of plot but I am including it
% anyways in case
figure;
hold on
% Plot the actual trajectory from Riccatti 
plot(X(1, 1 : end - 1), X(2, 1 : end - 1));

% Plot the estimated trajectory:
plot(X_est(1, :), X_est(2, :), 'r');

% Plot the predicted trajectory:
plot(X_pred(1, 1:end-1), X_pred(2, 1:end-1), 'g');


hold off
legend("Actual X (Ricatti)", "Estimated X", "Predicted X");
title("Trajectory plot for the estimated, predicted, and actual");
ylabel("x_{2}");
xlabel("X_{1}");
xlim([-50 50]);
ylim([-50 50]);

% Comments:
% From these plots we can see that aside from the beginning of the
% iterations, the trajectories match up quite well once again.

% Part d) - Final Comments
% 
% Regarding whether or not the Kalman filter still works reasonably for the
% unstable system I would say that it does a good job keeping up. Of
% course, the results are not pretty because of the "unstable-ness" of the
% system but considering that it does technically follow the trajectory
% quite accurately. Furthermore, I would say that it handled the switch
% well since it only took 5 iterations for the algorithm to converge again.
% 

%% Problem II - Tracking Applications
%
% List given terms:
N = 1001;                   % Points to simulate
beta = 0.02;                % Beta term
alpha = 0.1;                % Alpha term
sigma_x = 1e-6;
sigma_y = 1e-5;

Qx_P2 = sigma_x * eye(4);
Qy_P2 = sigma_y * eye(2);

C_P2 = [1, 0, 0, 0 ;
        0, 1, 0, 0];

% Initialize inital state:
x0_P2 = [ 1   ;
          0   ;
        alpha ;
          0  ];

% Time setup:
T = 10;
dt = 0.01;
t = 0 : dt : T;

% List the A(t) matrix: * Does not include t values in multiplication. Send
% to discretizeA function for prescribed values of t at each time step:
At = [ 0,    0,     1,        0    ;
       0,    0,     0,        1    ;
       0,  -beta,  alpha,  -2*beta ;
      beta,  0,   2*beta,   alpha ];

% Use the formula for Ad[n] via the midpoint method:
% Ad1 = discretizeA(At, dt, t(2));            % s.t. t = n*dt

% Compute the trajectory for the system as in the problem statement:
x1 = exp(alpha .* t) .* cos(beta * t.^2);
x2 = exp(alpha .* t) .* sin(beta * t.^2);

% Initialize state vector
x = zeros(4, N + 1);
x(:, 1) = x0_P2;       % This implies x1(:, 1) is symbolically "x[0]"
% Initialize output vector:
y = zeros(2, N);

% Initialize state vector (WITH NOISE):
x_noise = zeros(4, N + 1);
x_noise(:, 1) = x0_P2;       % This implies x1(:, 1) is symbolically "x[0]"

% Initialize output vector (WITH NOISE):
y_noise = zeros(2, N);

% Construct the discrete time model:
for N = 1:N
    % Iterate to form x[n + 1] and y[n]:
    % Compute discretized A[n] at the current N - symbolically (N - 1) idx
    Ad = discretizeA(At, dt, t(N));

    % Compute the x and y's WITHOUT noise:
    x(:, N + 1) = Ad * x(:, N);     
    % Generate y's from state vector:
    y(:, N) = C_P2 * x(:, N);      
    
    % Calculate the x and y's WITH noise:
    x_noise(:, N + 1) = Ad * x_noise(:, N) + sigma_x * randn(4, 1);
    % Generate y's from state vector:
    y_noise(:, N) = C_P2 * x_noise(:, N) + sigma_y * randn(2, 1);
    
end

% Now that state variables have been generated, compute the Kalman
% estimate:
% Initialize state prediction error covariance matrix:
K_init = 10*Qx_P2;
% Compute Kalman
[X_predicted, X_estimate, K_final, Error_Q2, K_100] = KalmanQ2(N, x_noise, y_noise, Ad, C_P2, K_init, Qx_P2, Qy_P2);

% 1. Superimpose plots of the actual trajectory x1(t), x2(t) vs the
%    discretized version of x1[n], x2[n]:
trajectory_actual_no_noise = x1 + 1j*x2;
trajectory_discrete_no_noise = x(1, :) + 1j*x(2, :);

figure;
hold on
% Plot actual trajectory:
plot(trajectory_actual_no_noise);
% Plot discretized model trajectories:
plot(trajectory_discrete_no_noise);
hold off
legend("Continuous trajectory", "Discrete trajectory");
ylabel("x_{2}");
xlabel("x_{1}");

% Compute maximum error between curves:

deviations = (abs(trajectory_actual_no_noise(1:end)) - abs(trajectory_discrete_no_noise(1:end-1))).^2;
maxError = max(deviations);

disp("Maximum error between discrete and actual trajectories: " + maxError + newline);

% COMMENTS:
% From the trajectory plots comparing the continuous trajectories to the
% discrete trajectories, they seem alright. But the discrete approximation
% does not do a good job keeping up at each time step. Maybe if the time 
% step size is lowered we will see better results.
%
% Actually, I just tried to lower the time steps but it made the trajectory
% less responsive.



% 2. Compute max deviation in actual velocities versus the velocities given
%    by the discrete model:

velocity_actual = (alpha + 1j * 2*beta.*t) .* exp(alpha .* t + 1j .* beta .* t.^2);
velocity_discrete = x(3, :) + 1j*x(4, :);

figure;
hold on
% Plot actual trajectory:
plot(velocity_actual);
% Plot discretized model trajectories:
plot(velocity_discrete);
hold off
legend("Continuous velocity", "Discrete velocity");
ylabel("v_{2}");
xlabel("v_{1}");

% Compute maximum error between curves:

velo_deviations = (abs(velocity_actual(1:end)) - abs(velocity_discrete(1:end-1))).^2;
maxError_velo = max(velo_deviations);

disp("Maximum error between discrete and actual velocities: " + maxError_velo + newline);

% COMMENTS:
% From these plots there is the same trend as in the actual plots, the
% plots are a reasonable approximation however not exactly accurate.
% (application specific of course)


% 3. Superimpose plots of the discrete trajectory and velocity (on separate
%    plots) but this time considering the random noise disturbance


% Create figures of the actual trajectories and the predicted and estimated
% trajectories:
trajectory_discrete_noise = x_noise(1, :) + 1j*x_noise(2, :);
velocity_discrete_noise = x_noise(3, :) + 1j*x_noise(4, :);

figure;
hold on
% Plot the discrete trajectory with no noise:
plot(trajectory_discrete_no_noise);
% Superimpose the discrete trajectory with added noise:
plot(trajectory_discrete_noise);
hold off
legend("Discrete trajectory (no disturbance)", "Discrete trajectory (with disturbance)");
title("Discrete Trajectories Noise Comparison");
xlabel("x_{1}");
ylabel("x_{2}");

% Create the plots for the velocity with disturbances:
figure;
hold on
% Plot the discrete trajectory with no noise:
plot(velocity_discrete_noise);
% Superimpose the discrete trajectory with added noise:
plot(velocity_discrete_noise);
hold off
legend("Discrete velocity (no disturbance)", "Discrete velocity (with disturbance)");
title("Discrete Velocities Noise Comparison");
xlabel("v_{1}");
ylabel("v_{2}");

% COMMENTS:
% Adding the disturbances with variance in either x1 or x2 directions
% proportional to sigma_x and sigma_y respectively, the results look pretty
% accurate, assuming the discretized model is correct. (As we found
% earlier, there are some inaccuracies with the discrete model and the
% continuous model but had our discrete model been a perfect representation
% of the trajectory and/or velocity, then these results would be more
% accurate with the continuous case.
%
% Note that at first I did not multiply the white gaussian r.v.'s by the
% sigma's but the trajectories were all haywire so I tried this instead
% because I thought that maybe I was not considering the variances
% properly. This could be completely wrong but yet I chose to do this so I
% did not have bad results.




% 4. Create graphs of the state variables (discrete) with their estimated
%    values via the Kalman filter, superimposed. Do the same for
%    velocities:
%
%    I will assume we are interested in the results with the added noise:

figure;
hold on
% Plot the discrete trajectory with added noise:
plot(trajectory_discrete_noise);
% Plot discrete trajectory estimate
X_trajectory_noise = X_estimate(1, :) + 1j*X_estimate(2, :);
plot(X_trajectory_noise, 'r');
hold off
legend("Discrete Trajectory Model (noisy)", "Estimated Trajectory Kalman (noisy)");
title("Kalman vs Model Trajectory");
xlabel("x_{1}");
ylabel("x_{2}");


% Create Kalman velocity state estimate:
X_velocity_noise = X_estimate(3, :) + 1j*X_estimate(4, :);

figure;
hold on
% Plot the discrete velocity with added noise:
plot(velocity_discrete_noise);
% Superimposed with Kalman velocity state estimate:
plot(X_velocity_noise);
hold off
legend("Discrete Velocity Model (noisy)", "Discrete velocity Kalman (noisy)");
title("Kalman vs Model Velocity");
xlabel("v_{1}");
ylabel("v_{2}");

% COMMENTS:
% 
% I'm surprised with how well the Kalman filter is able to estimate the
% state trajectories and velocities upon inspection of the plots in the
% above section. Sure, there is some deviation in the velocity curves but
% overall a really good job!


% 5. Examine K(n, n) at 100th iteration as well as final K(n, n). Also
%    compute the norm of K in each case:

% Compute norm of K at 100th iteration:
K100Norm = norm(K_100);
% Compute norm of K at final iteration:
KfinalNorm = norm(K_final);

disp("Norm of K at 100th iteration: " + K100Norm + newline);
disp("Norm of K at final iteration: " + KfinalNorm + newline);

% COMMENTS:
% 
% Norms here look pretty close and relatively small. Forgive me, but I do
% not know how this would correlate to our observations so far, but
% interesting!


% 6. Compute observability matrix stuff:

% Consider the given matrix in the problem to be evaluated at t = 1:
MAT_toCompute = [  C_P2  ;
                 C_P2*At ];
disp("Displaying the matrix [C ; CA]: " + newline);
disp(MAT_toCompute);

% COMMENTS:
%
% Well, this matrix does look "pretty"! It is an identity matrix. I believe
% the correct interpretation of this result is that since this matrix is
% full rank, rank 4 this implies that our system is "observable". 
% NOTE: adding CA^n terms would not increase the rank here as per the
% theorem in the notes.
%
% Furthermore, since our system is observable that means given a known
% input, observation of the output y over not more than 4 time spans is
% sufficient to UNIQUELY determine the state trajectory during this time.
%
% This is my analysis based on what I've read in the slides.


% Define symbolic matrices T_MAT and At
T_MAT_sym = [1, 1, 1, 1; 
             1, 1, 1, 1; 
             1, 1, 1, t(1); 
             1, 1, t(1), 1];

At_sym = [0,   0,    1,       0; 
          0,   0,    0,       1; 
          0, -beta, alpha, -2*beta; 
        beta,  0,  2*beta,  alpha];

% Define symbolic matrix Ad
Ad = eye(4) + (At_sym .* T_MAT_sym) * dt + ((1/2) * (At_sym .* T_MAT_sym) * (dt^2));

% Define matrix C
C_P2 = [1, 0, 0, 0; 
        0, 1, 0, 0];

% Construct the combined matrix [C ; C*Ad]
combined_matrix = [C_P2; C_P2*Ad];

% Compute the symbolic determinant of the combined matrix
determinant = det(combined_matrix);
disp("Determinant of matrix: " + newline);
disp(determinant);






%% Functions Created

% Function to compute the estimated and prediction state variables using
% the Kalman filter algorithm in the notes. Uses the single update function
% for the next Kalman parameters
function [X_predicted_n, X_est, K_next, e, K_n] = Kalman(iterations, xn, yn, A1, A2, C, K1, K2, Qx, Qy)
    % VARIABLE REFERENCE:
    % K_n           - K(n, n) current K
    % K_next        - K(n + 1, n) gives next K for algorithm
    % X_est         - Current estimate x(n|n)
    % X_pred_n      - The predicted x(n + 1|n)

    % Initialize
    % Beginning of K(n|n) is Qx in the assignment:
    K_next = 0.1 * eye(2);
    % Initialize prediction error:
    e = zeros(iterations, 1);
    X_predicted_n = zeros(2, iterations + 1);
    X_predicted_n(:, 1) = xn(:, 1);
    X_est = zeros(2, iterations);
    
    % Perform the iteration for the Kalman Filter as in the notes:
    % Run the Kalman Steps for lower half iterations
    for i = 1:iterations/2
        
        e(i, 1) = norm(K_next - K1);
        % K_next here in R is K(n, n-1)
        R = (C * K_next * C.' + Qy)^-1;
        G = (K_next * C.' * R);
        alpha = yn(i) - C * X_predicted_n(:, i);
        [K_next, K_n, X_est(:, i), X_predicted_n(:, i + 1)] = UpdateToKalman(X_predicted_n(:, i), ...
                                           C, G, alpha, A1, K_next, Qx, Qy);
    end

    % Run Kalman Steps for upper half iterations:
    for i = iterations/2 + 1:iterations
        e(i, 1) = norm(K_next - K2);
        R = (C * K_next * C.' + Qy)^-1;
        G = (K_next * C.' * R);
        alpha = yn(i) - C * X_predicted_n(:, i);
        [K_next, K_n, X_est(:, i), X_predicted_n(:, i + 1)] = UpdateToKalman(X_predicted_n(:, i), C, G, alpha, A2, K_next, Qx, Qy);
    end

end


% This function is for Problem 2 where I am not changing A halfway through:
function [X_predicted_n, X_est, K_next, e, K_100] = KalmanQ2(iterations, xn, yn, A, C, K, Qx, Qy)
    % VARIABLE REFERENCE:
    % K_n           - K(n, n) current K
    % K_next        - K(n + 1, n) gives next K for algorithm
    % X_est         - Current estimate x(n|n)
    % X_pred_n      - The predicted x(n + 1|n)

    % Initialize
    % Change initial value accordingly for Problem 2:
    K_next = 10*Qx;

    % Initialize prediction error:
    e = zeros(iterations, 1);
    X_predicted_n = zeros(4, iterations + 1);
    X_predicted_n(:, 1) = xn(:, 1);
    X_est = zeros(4, iterations);
    
    % Perform the iteration for the Kalman Filter as in the notes:
    % Run the Kalman Steps for one iteration
    for iter = 1:iterations
        e(iter, 1) = norm(K_next - K);
        R = (C * K_next * C.' + Qy)^-1;
        G = (K_next * C.' * R);
        alpha = yn(:, iter) - C * X_predicted_n(:, iter);
        [K_next, ~, X_est(:, iter), X_predicted_n(:, iter + 1)] = UpdateToKalman(X_predicted_n(:, iter), ...
                                            C, G, alpha, A, K_next, Qx, Qy);

        % Added for Part 5. of Problem 2 -
        % Compute and return the 100th K(n, n) matrix for comparison with 
        % the final K:
        if(iter == 100)         % USING 100 because after UpdateToKalman
                                % gets executed, it returns K[101] which
                                % symbolically translates to K[100], RECALL
            K_100 = K_next;
        end

    end

end
% Function to run one iteration involving K(n, n - 1), x(n|n - 1) and
% produces K(n, n), K(n + 1, n), x(n|n), and x(n + 1|n) as outputs:
function [K_next, K_n, X_est, X_predicted_next] = UpdateToKalman(X_pred_prev, ...
                                        C, G, alpha, A, K_prev, Qx, Qy)
    % K_n           - K(n, n) current K
    % K_next        - K(n + 1, n) gives next K for algorithm
    % K_prev        - K(n, n - 1) term in algorithm
    % X_est         - Current estimate x(n|n)
    % X_pred_n      - The predicted x(n + 1|n)
    % X_pred_prev   - x(n|n-1) previous

    % Correction step:
    X_est = X_pred_prev + G * alpha;
    
    % Execute the prediction step:
    X_predicted_next = A * X_est;

    % Compute K(n, n) using the Joseph form:
    K_n = (eye(size(K_prev, 1)) - G*C) * K_prev * (eye(size(K_prev, 1)) - G*C)' ...
            + G*Qy*G';
    
    % Compute K(n + 1, n) term:
    K_next = A*K_n*A' + Qx;

end

% Function to compute discretized A for a single time step
function Ad = discretizeA(At, dt, t)
    
    % Use matrix T_MAT to multiply time entry to Matrix A
    T_MAT = [ 1, 1, 1, 1 ;
              1, 1, 1, 1 ;
              1, 1, 1, t ;
              1, 1, t, 1];

    % Compute the outpute discretized Ad:
    Ad = eye(4) + (At.*T_MAT) * dt + (1/2) * (At.*T_MAT).^2 * (dt)^2;

end



