function F = DCMotorJacobianFcn(x, Va)
    Ts=0.001;
    Kt=0.01;
    % State vector
    theta = x(1);
    omega = x(2);
    ia = x(3);
    Ra = x(4);
    La = x(5);
    Kb = x(6);
    Jeq = x(7);
    Deq = x(8);
    
    % Discrete-time state transition Jacobian
    F = eye(8);
    F(1,2) = Ts;  % dtheta/domega
    F(2,2) = 1 - (Ts*Deq/Jeq);  % domega/domega
    F(2,3) = Ts*(Kt/Jeq);  % domega/dia
    F(2,7) = -Ts*(Kt*ia + Deq*omega)/(Jeq^2);  % domega/dJeq
    F(2,8) = -Ts*omega/Jeq;  % domega/dDeq
    
    F(3,2) = -Ts*(Kb/La);  % dia/domega
    F(3,3) = 1 - Ts*(Ra/La);  % dia/dia
    F(3,4) = -Ts*(ia/La);  % dia/dRa
    F(3,5) = -Ts*(Va - Ra*ia - Kb*omega)/La^2;  % dia/dLa
    F(3,6) = -Ts*(omega/La);  % dia/dKb

end 