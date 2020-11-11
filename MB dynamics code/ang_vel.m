function [omega1,S_omega1,dRf1_0] = ang_vel(Rf_1_0, dRf1_0, q, dq) 

    % Derivative of rotation matrix
    dRf1_0(:) = [jacobian(Rf_1_0(:,1), q)*dq, jacobian(Rf_1_0(:,2), q)*dq, jacobian(Rf_1_0(:,3), q)*dq];
    dRf1_0 = simplify(dRf1_0,'IgnoreAnalyticConstraints',true);

    % S(omega) = (transpose(R))*dR/dt
    S_omega1 = simplify(transpose(Rf_1_0) * dRf1_0,'IgnoreAnalyticConstraints',true);

    % Angular velocity
    omega1 = [S_omega1(3,2);S_omega1(1,3);S_omega1(2,1)];
end

%{
    
    for i = 1:length(Rf_1_0)
        dRf1_0(:,i) = jacobian(Rf_1_0(:,i), q)*dq; 
    end
    
    %}