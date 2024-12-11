function [sys, x0, str, ts, simStateCompliance] = EKF_All_Params(t, x, u, flag, varargin)
% Discrete-Time DC Motor Parameter Estimation EKF
% Parameters now evolve with process noise

% Default parameters
if nargin < 5 || isempty(varargin{1})
    P0 = diag([1*ones(3,1); 0.01*ones(6,1)]);  % Initial covariance
else
    P0 = varargin{1};
end

if nargin < 6 || isempty(varargin{2})
    Q = diag([0.11; 0.11; 0.13; 0.01*ones(4,1); 0.9;0.015]);   % Process noise covariance

else
    Q = varargin{2};
end

if nargin < 7 || isempty(varargin{3})
    R = 0.001;             % Measurement noise covariance
else
    R = varargin{3};
end

if nargin < 8 || isempty(varargin{5})
    Kt = 0.01;           % Torque constant
else
    Kt = varargin{5};
end

% Sample time (to be set in Simulink)
Ts = 0.001;  % Default sample time, should match Simulink setting

persistent P;  % Covariance matrix

switch flag
    case 0   % Initialization
        sizes = simsizes;
        sizes.NumContStates  = 0;   % Discrete system
        sizes.NumDiscStates  = 9;   % State vector size
        sizes.NumOutputs     = 9;  % State + Covariance matrix
        sizes.NumInputs      = 2;   % Voltage and Angle inputs
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;
        
        sys = simsizes(sizes);
        
        % Initial state estimates (include Kt)
        x0 = [0;   % θ
              0;   % ω
              0;   % ia
              1;   % Ra
              0.1; % La
              0.01;% Kb
              0.01;% Kt
              0.1; % Jeq
              0.1  % Deq
             ];
        
        % Initialize covariance matrix
        P = P0;     
        
        str = [];
        ts  = [Ts 0];   % Discrete sample time
        simStateCompliance = 'UnknownSimState';
        
    case 2   % Update (Discrete-time EKF)
        % Inputs
        Va = max([0;u(1)]);  % Input voltage
        y = u(2);   % Measured angle
        
        % Nonlinear discrete-time state prediction with parameter evolution
        x_pred = discrete_state_update(x, Va, Ts);
        
        % Compute Jacobians
        [F, H] = compute_discrete_jacobians(x, Va, Ts);
      
        % Predicted Covariance
        P_pred = F * P * F' +  Q ;        
      
        % Kalman Gain
        K = P_pred * H' / (H * P_pred * H' + R);
        
        % State Update
        x_updated = x_pred + K * (y - x_pred(1));
        
        % Covariance Update
        P = (eye(9) - K*H) * P_pred';
        
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
function x_next = discrete_state_update(x, Va, Ts)
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
    
    % Discrete-time state equations (Euler approximation)
    x_next = x(1:9);
    
    % Kinematic and dynamic equations
    x_next(1) = theta + Ts * omega;  % θ
    x_next(2) = omega + Ts * ((1/Jeq) * (Kt*ia - Deq*omega));  % ω
    x_next(3) = ia + Ts * ((1/La) * (Va - Ra*ia - Kb*omega));  % ia
    
    % Parameter evolution with slight random walk
    % Add small random noise to allow parameter estimation
    x_next(4) = Ra ;  % Ra with minimal drift
    x_next(5) = La ;  % La with minimal drift
    x_next(6) = Kb ;  % Kb with minimal drift
    x_next(7) = Kt ; 
    x_next(8) = Jeq ; % Jeq with minimal drift 
    x_next(9) = Deq ; % Deq with minimal drift
end

% Separate function for Jacobian computation
function [F, H] = compute_discrete_jacobians(x, Va, Ts)
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
    
    % Discrete-time state transition Jacobian
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