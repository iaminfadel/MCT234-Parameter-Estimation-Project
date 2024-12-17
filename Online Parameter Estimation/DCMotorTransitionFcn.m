function [x_next] = DCMotorTransitionFcn(x, u)
% State Transition Function for DC Motor Dual Estimator
% 
% Inputs:
%   x       - Extended state vector 
%   u       - Input vector [voltage]
%   dt      - Time step (discretization interval)
%   params  - Initial parameter estimates
%
% Outputs:
%   x_next  - Next state vector after time update
%
% Extended State Vector Components:
%   x(1)  - Angular position (θ)
%   x(2)  - Angular velocity (ω)
%   x(3)  - Motor current (i)
%   x(4)  - Estimated moment of inertia (J)
%   x(5)  - Estimated viscous friction (b)
%   x(6)  - Estimated back-EMF constant (Ke)
%   x(7)  - Estimated torque constant (Kt)
%   x(8)  - Estimated electrical resistance (R)
%   x(9)  - Estimated electrical inductance (L)
dt = 5e-3;
R = 5.4;
K = 0.018;

% Extract current states and parameter estimates
omega = x(1);
i = x(2);
% Parameter estimates
J = x(3);
B = x(4);
L = x(5);
% Input voltage
V = u(1);

% State Update Equations with Parameter Estimates

% Velocity update (using estimated parameters)
domega = ((1/J) * (K*i - B*omega)) * dt;

% Current update (using estimated parameters)
di = ((1/L) * (V - R*i - K*omega)) * dt;

% Simple parameter evolution model
% These can be adjusted based on expected parameter variability
dJ = 0;      % Moment of inertia typically constant
db = 0;      % Friction coefficient slow-changing
dL = 0;      % Inductance typically stable

% Construct next state vector
x_next = [
    omega + domega;     % Updated velocity
    i + di;             % Updated current
    J + dJ;         % Estimated J (nearly constant)
    B + db;         % Estimated b 
    L + dL;         % Estimated L
];
end