function [x_next] = DCMotorTransitionFcn(x, u, dt)
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

% Extract current states and parameter estimates
theta = x(1);
omega = x(2);
i = x(3);

% Parameter estimates
J_est = x(4);
b_est = x(5);
Ke_est = x(6);
Kt_est = x(7);
R_est = x(8);
L_est = x(9);

% Input voltage
V = u(1);

% State Update Equations with Parameter Estimates
% Position update
dtheta = omega * dt;

% Velocity update (using estimated parameters)
domega = ((1/J_est) * (Kt_est*i - b_est*omega)) * dt;

% Current update (using estimated parameters)
di = ((1/L_est) * (V - R_est*i - Ke_est*omega)) * dt;

% Simple parameter evolution model
% These can be adjusted based on expected parameter variability
dJ = 0;      % Moment of inertia typically constant
db = 0;      % Friction coefficient slow-changing
dKe = 0;     % Back-EMF constant relatively stable
dKt = 0;     % Torque constant relatively stable
dR = 0;      % Resistance can change with temperature
dL = 0;      % Inductance typically stable

% Construct next state vector
x_next = [
    theta + dtheta;     % Updated position
    omega + domega;     % Updated velocity
    i + di;             % Updated current
    J_est + dJ;         % Estimated J (nearly constant)
    b_est + db;         % Estimated b 
    Ke_est + dKe;       % Estimated Ke
    Kt_est + dKt;       % Estimated Kt
    R_est + dR;         % Estimated R
    L_est + dL;         % Estimated L
];
end