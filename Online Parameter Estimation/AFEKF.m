function [sys, x0, str, ts, simStateCompliance] = AFEKF(t, x, u, flag, varargin)
% Discrete-Time DC Motor Parameter Estimation Adaptive Fading EKF
% Parameters now evolve with process noise

% Default parameters
if nargin < 5 || isempty(varargin{1})
    P0 = diag([1*ones(3,1); 1*ones(2,1)]);  % Initial covariance
else
    P0 = varargin{1};
end

if nargin < 6 || isempty(varargin{2})
    Q = diag([0.01*ones(3,1); 1; 0.1]);   % Process noise covariance

else
    Q = varargin{2};
end

if nargin < 7 || isempty(varargin{3})
    R = 0.01;             % Measurement noise covariance
else
    R = varargin{3};
end

if nargin < 8 || isempty(varargin{4})
    gamma = 0.5;         % Adaptive fading gain
else
    gamma = varargin{4};
end

if nargin < 9 || isempty(varargin{5})
    Ra = 2;
    K = 0.01;
else
    K = varargin{5};
end

% Sample time (to be set in Simulink)
Ts = 0.001;  % Default sample time, should match Simulink setting

persistent P;  % Covariance matrix
persistent lambda;  % Fading factor

switch flag
    case 0   % Initialization
        sizes = simsizes;
        sizes.NumContStates  = 0;   % Discrete system
        sizes.NumDiscStates  = 5;   % State vector size
        sizes.NumOutputs     = 5;  % State + Covariance matrix
        sizes.NumInputs      = 2;   % Voltage and Angle inputs
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;
        
        sys = simsizes(sizes);
        
        % Initial state estimates (include K)
        x0 = [0;   % ω
              0;   % ia
              0.1;   % La
              0.1; % Jeq
              0.0 % Deq
             ];
        
        % Initialize covariance matrix
        P = P0;
        
        % Initialize fading factor
        lambda = 1;
        
        str = [];
        ts  = [Ts 0];   % Discrete sample time
        simStateCompliance = 'UnknownSimState';
        
    case 2   % Update (Discrete-time EKF)
        % Inputs
        Va = max([0;u(1)]);  % Input voltage
        y = u(2);   % Measured angle
        
        % Nonlinear discrete-time state prediction with parameter evolution
        x_pred = discrete_state_update(x, Va, Ts, K, Ra);
        
        % Compute Jacobians
        [F, H] = compute_discrete_jacobians(x, Va, Ts, K, Ra);
        
        % Adaptive Fading Mechanism
        %lambda = max(1, gamma * (innovation'*innovation));

        % Predicted Covariance
        P_pred = F * P * F' +  Q ;        
        
        % Kalman Gain
        K = P_pred * H' / (H * P_pred * H' + R);
        
        % State Update
        x_updated = x_pred + K * (y - x_pred(1));
        
        % Covariance Update
        P = (eye(5) - K*H) * P_pred';
        
        sys = x_updated;
        
    case 3   % Outputs
        % Output state and covariance matrix
        assignin('base', 'current_P_matrix', P);
        sys = x(:);
        
    case 9   % Terminate
        sys = [];
        
    otherwise
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
end

% Separate function for discrete-time state update
function x_next = discrete_state_update(x, Va, Ts, K, Ra)
    % State vector unpacking
    omega = x(1);
    ia = x(2);
    Jeq = x(3);
    Deq = x(4);
    La = x(5);
    % Discrete-time state equations (Euler approximation)
    x_next = x(1:5);
    
    % Kinematic and dynamic equations
    x_next(1) = omega + Ts * ((1/Jeq) * (K*ia - Deq*omega));  % ω
    x_next(2) = ia + Ts * ((1/La) * (Va - Ra*ia - K*omega));  % ia
    
    % Parameter evolution with slight random walk
    % Add small random noise to allow parameter estimation
    x_next(3) = Jeq ; % Jeq with minimal drift 
    x_next(4) = Deq ; % Deq with minimal drift
    x_next(5) = La ; % Deq with minimal drift

end

% Separate function for Jacobian computation
function [F, H] = compute_discrete_jacobians(x, Va, Ts, K, Ra)
    % State vector
    omega = x(1);
    ia = x(2);
    Jeq = x(3);
    Deq = x(4);
    La = x(5);

    % Discrete-time state transition Jacobian
    F = eye(5);
    F(1,1) = 1 - (Ts*Deq/Jeq);  % domega/domega
    F(1,2) = Ts*(K/Jeq);  % domega/dia
    F(1,3) = -Ts*(K*ia + Deq*omega)/(Jeq^2);  % domega/dJeq
    F(1,4) = -Ts*omega/Jeq;  % domega/dDeq
    F(1,5) = 0;  % domega/dLa
  
    F(2,1) = -Ts*(K/La);  % dia/domega
    F(2,2) = 1 - Ts*(Ra/La);  % dia/dia
    F(2,5) = -Ts*(Va - Ra*ia - K*omega)/La^2;  % dia/dLa

    % Measurement Jacobian (only angle is measured)
    H = [1, 0, 0, 0, 0];
end