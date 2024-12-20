    function [sys, x0, str, ts, simStateCompliance] = EKF_All_Params_Forgetting(t, x, u, flag, varargin)
% Discrete-Time DC Motor Parameter Estimation EKF with Forgetting Factor
% Allows adaptive tracking of parameter changes

% Default parameters
if nargin < 5 || isempty(varargin{1})
    P0 = diag([1*ones(3,1)]);  % Initial covariance
else
    P0 = varargin{1};
end

if nargin < 6 || isempty(varargin{2})
    Q = diag([ 1.2; 5.01; 0.3]);   % Process noise covariance
else
    Q = varargin{2};
end

if nargin < 7 || isempty(varargin{3})
    R = 30.1;             % Measurement noise covariance
else
    R = varargin{3};
end

% Forgetting factor (added parameter)
if nargin < 8 || isempty(varargin{4})
    lambda =1;  % Forgetting factor (between 0 and 1)
else
    lambda = varargin{4};
end

if nargin < 9 || isempty(varargin{5})
    Kt = 0.018; 
    Ra = 5.4;
    La = 0.05;
    Deq = 0.03;
else
    Kt = varargin{5};
end

% Sample time (to be set in Simulink)
Ts = 5e-3;  % Default sample time, should match Simulink setting

persistent P;  % Covariance matrix

switch flag
    case 0   % Initialization
        sizes = simsizes;
        sizes.NumContStates  = 0;   % Discrete system
        sizes.NumDiscStates  = 3;   % State vector size
        sizes.NumOutputs     = 3;   % State + Covariance matrix
        sizes.NumInputs      = 2;   % Voltage and Angle inputs
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;
        
        sys = simsizes(sizes);
        
        % Initial state estimates (include Kt)
        x0 = [0;   % ω
              0;   % ia
              0.1; % Jeq
             ];
        
        % Initialize covariance matrix
        P = P0;     
        
        str = [];
        ts  = [Ts 0];   % Discrete sample time
        simStateCompliance = 'UnknownSimState';
        
    case 2   % Update (Discrete-time EKF with Forgetting Factor)
        % Inputs
        Va = u(1);  % Input voltage
        y = u(2);   % Measured angle
        Params = [Kt,Ra,La,Deq];
        
        % Nonlinear discrete-time state prediction with parameter evolution
        x_pred = discrete_state_update(x, Va, Params, Ts);
        
        % Compute Jacobians
        [F, H] = compute_discrete_jacobians(x, Va, Params, Ts);
      
        % Predicted Covariance with Forgetting Factor
        % Increases uncertainty over time to allow faster adaptation
        P_pred = (F * P * F') ./ lambda + Q;        
      
        % Kalman Gain
        K = P_pred * H' / (H * P_pred * H' + R);
        
        % State Update
        x_updated = x_pred + K * (y - x_pred(1));
        
        % Covariance Update with Forgetting Factor
        P = (eye(3) - K*H) * P_pred' ./ lambda;
        
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

% Separate function for discrete-time state update (unchanged)
function x_next = discrete_state_update(x, Va, Params, Ts)
    % State vector unpacking
    omega = x(1);
    ia = max(0.01,x(2));
    Jeq = min(max(0.01, x(3)), 100);

    Kt = Params(1);
    Ra = Params(2);
    La = Params(3);
    Deq = Params(4);
    
    % Discrete-time state equations (Euler approximation)
    x_next = x(1:3);
    
    % Kinematic and dynamic equations
    x_next(1) = omega + Ts * ((1/Jeq) * (Kt*ia*97.5 - Deq*omega));  % ω
    x_next(2) = ia + Ts * ((1/La) * (Va - Ra*ia - Kt*(omega*97.5)));  % ia
    
    % Parameter evolution with slight random walk
    x_next(3) = Jeq; % Jeq with minimal drift 
end

% Separate function for Jacobian computation (unchanged)
function [F, H] = compute_discrete_jacobians(x, Va, Params, Ts)
    % State vector
    omega = x(1);
    ia = x(2);
    Jeq = x(3);

    Kt = Params(1);
    Ra = Params(2);
    La = Params(3);
    Deq = Params(4);
    
    % Discrete-time state transition Jacobian
    F = eye(3);
    F(1,1) = 1 - (Ts*Deq/Jeq);  % domega/domega
    F(1,2) = Ts*(97.5*Kt/Jeq);  % domega/dia
    F(1,3) = -Ts*(Kt*ia*97.5 - Deq*omega)/(Jeq^2);  % domega/dJeq
    
    F(2,1) = -Ts*(Kt*97.5/La);  % dia/domega
    F(2,2) = 1 - Ts*(Ra/La);  % dia/dia
    
    % Measurement Jacobian (only angle is measured)
    H = [1, 0, 0];
end