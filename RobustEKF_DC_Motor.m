function [sys, x0, str, ts, simStateCompliance] = RobustEKF_DC_Motor(t, x, u, flag, varargin)
% Robust Extended Kalman Filter for DC Motor Parameter Estimation
% Includes numerical stability and adaptive filtering techniques

% Persistent variables for adaptive filtering
persistent P adaptive_Q adaptive_R;
persistent innovation_history;
persistent iteration_count;

% Default parameters with improved initialization
if nargin < 5 || isempty(varargin{1})
    P0 = diag([10*ones(3,1); 1*ones(6,1)]);  % More conservative initial covariance
else
    P0 = varargin{1};
end

% Adaptive noise covariance initialization
if isempty(adaptive_Q)
    adaptive_Q = diag([0; 0; 0; 0.01*ones(4,1); 0.5; 0.01]);
end

if isempty(adaptive_R)
    adaptive_R = 0.01;
end

% Initialize innovation history for adaptive filtering
if isempty(innovation_history)
    innovation_history = zeros(9, 1000);  % Store last 200 innovations
end

if isempty(iteration_count)
    iteration_count = 0;
end

% Sample time (to be set in Simulink)
Ts = 0.001;  % Default sample time

switch flag
    case 0   % Initialization
        sizes = simsizes;
        sizes.NumContStates  = 0;   % Discrete system
        sizes.NumDiscStates  = 9;   % State vector size
        sizes.NumOutputs     = 9;   % State + Covariance matrix
        sizes.NumInputs      = 2;   % Voltage and Angle inputs
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;
        
        sys = simsizes(sizes);
        
        % Initial state estimates with wider initialization
        x0 = [0;    % θ
              0;    % ω
              0;    % ia
              1;    % Ra
              0.1;  % La
              0.01; % Kb
              0.01; % Kt
              0.1;  % Jeq
              0.1   % Deq
             ];
        
        % Initialize covariance matrix
        P = P0;     
        
        str = [];
        ts  = [Ts 0];   % Discrete sample time
        simStateCompliance = 'UnknownSimState';
        
    case 2   % Update (Discrete-time Robust EKF)
        % Inputs
        Va = max([0;u(1) - 2]);  % Input voltage
        y = u(2);   % Measured angle
        
        % Numerical Stability: Add small regularization to prevent matrix ill-conditioning
        epsilon = 1e-10;
        
        % Nonlinear discrete-time state prediction with improved parameter evolution
        x_pred = robust_discrete_state_update(x, Va, Ts);
        
        % Compute Jacobians with numerical stability
        [F, H] = compute_robust_jacobians(x, Va, Ts);
      
        % Robust Covariance Prediction with Regularization
        P_pred = F * P * F' + adaptive_Q;
        
        % Numerical Stability: Symmetrize and regularize covariance
        P_pred = 0.5 * (P_pred + P_pred') + epsilon * eye(size(P_pred));
        
        % Chi-square test for innovation validation
        innovation = y - x_pred(1);
        innovation_threshold = 3 * sqrt(H * P_pred * H' + adaptive_R);
        
        if abs(innovation) <= innovation_threshold
            % Kalman Gain with Robust Inversion
            K = P_pred * H' / (H * P_pred * H' + adaptive_R + epsilon);
            
            % State Update
            x_updated = x_pred + K * innovation;
            
            % Covariance Update with Joseph form for positive definiteness
            I_KH = eye(9) - K * H;
            P = I_KH * P_pred * I_KH' + K * adaptive_R * K';
            
            % Adaptive Noise Estimation
            iteration_count = iteration_count + 1;
            
            % Update innovation history
            innovation_history(:, mod(iteration_count, 1000) + 1) = innovation;
            
            % Periodically update noise covariances
            if mod(iteration_count, 50) == 0
                [adaptive_Q, adaptive_R] = estimate_noise_covariances(innovation_history);
            end
        else
            % If innovation is too large, use prior estimate
            x_updated = x_pred;
        end
        
        sys = x_updated;
        
    case 3   % Outputs
        % Output state and covariance matrix
        assignin('base', 'current_P_matrix', P);
        assignin('base', 'adaptive_Q', adaptive_Q);
        assignin('base', 'adaptive_R', adaptive_R);
        sys = x(:);
        
    case 9   % Terminate
        sys = [];
        
    otherwise
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
end

% Robust state update with improved parameter evolution
function x_next = robust_discrete_state_update(x, Va, Ts)
    % State vector unpacking
    theta = x(1);
    omega = x(2);
    ia = x(3);
    Ra = x(4);
    La = x(5);
    Kb = x(6);
    Kt = x(7);
    Jeq = x(8);
    Deq = x(9);
    
    % Discrete-time state equations with parameter drift
    x_next = x;
    
    % Kinematic and dynamic equations
    x_next(1) = theta + Ts * omega;  % θ
    x_next(2) = omega + Ts * ((1/Jeq) * (Kt*ia - Deq*omega));  % ω
    x_next(3) = ia + Ts * ((1/La) * (Va - Ra*ia - Kb*omega));  % ia
    
    % Improved parameter evolution with controlled random walk
    param_drift_std = 0.000001;
    x_next(4) = Ra * (1 + param_drift_std * randn);     % Ra with controlled drift
    x_next(5) = La * (1 + param_drift_std * randn);     % La with controlled drift
    x_next(6) = Kb * (1 + param_drift_std * randn);     % Kb with controlled drift
    x_next(7) = Kt * (1 + param_drift_std * randn);     % Kt with controlled drift
    x_next(8) = Jeq * (1 + param_drift_std * randn);    % Jeq with controlled drift
    x_next(9) = Deq * (1 + param_drift_std * randn);    % Deq with controlled drift
end

% Robust Jacobian computation
function [F, H] = compute_robust_jacobians(x, Va, Ts)
    % State vector
    theta = x(1);
    omega = x(2);
    ia = x(3);
    Ra = x(4);
    La = x(5);
    Kb = x(6);
    Kt = x(7);
    Jeq = x(8);
    Deq = x(9);
    
    % Discrete-time state transition Jacobian with improved numerical conditioning
    F = eye(9);
    F(1,2) = Ts;  % dtheta/domega
    F(2,2) = 1 - (Ts*Deq/Jeq);  % domega/domega
    F(2,3) = Ts*(Kt/Jeq);  % domega/dia
    F(2,7) = Ts*(ia/Jeq);
    F(2,8) = -Ts*(Kt*ia - Deq*omega)/(Jeq^2);  % domega/dJeq
    F(2,9) = -Ts*omega/Jeq;  % domega/dDeq
    
    F(3,2) = -Ts*(Kb/La);  % dia/domega
    F(3,3) = 1 - Ts*(Ra/La);  % dia/dia
    F(3,4) = -Ts*(ia/La);  % dia/dRa
    F(3,5) = -Ts*(Va - Ra*ia - Kb*omega)/La^2;  % dia/dLa
    F(3,6) = -Ts*(omega/La);  % dia/dKb

    % Measurement Jacobian (only angle is measured)
    H = [1, 0, 0, 0, 0, 0, 0, 0, 0];
end

% Adaptive noise covariance estimation
function [new_Q, new_R] = estimate_noise_covariances(innovation_history)
    % Compute innovation statistics
    mean_innovation = mean(innovation_history, 2);
    cov_innovation = cov(innovation_history');
    
    % Adaptive Q estimation (process noise)
    new_Q = diag(abs(mean_innovation) + 0.1);
    
    % Adaptive R estimation (measurement noise)
    new_R = mean(diag(cov_innovation));
end