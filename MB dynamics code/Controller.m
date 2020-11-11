%% Control

% This function sets up the differential equation as well as the control
% Simple PD control is used
function [out] = Controller(t, y)
    
    %% Variables and default conditions
    %----------------------------------------------------------------------
    global ddpsi2_eqn ddphi1_eqn ddth3_eqn;
    global number_of_links;
    syms dphi1 ddphi1 dpsi2 ddpsi2 dth3 ddth3 'real';
    syms t1 t2 t3 'real';

    % Reference angles and angular velocities
    phi1_ref = (-1/4)*pi; 
    psi2_ref = (1/4)*pi;
    th3_ref = 0*pi;
    dphi1_ref = 0;  
    dpsi2_ref = 0;
    dth3_ref = 0;
    %----------------------------------------------------------------------
    
    %% Second order ODE's
    %----------------------------------------------------------------------
    ycell = num2cell(y);
    if number_of_links == 1
        [phi1, dphi1] = ycell{:};
    elseif number_of_links == 2
        [phi1, dphi1, psi2, dpsi2] = ycell{:};
    else
        [phi1, dphi1, psi2, dpsi2, th3, dth3] = ycell{:};
    end
    %----------------------------------------------------------------------
    
    %% PD Control
    %----------------------------------------------------------------------
    if number_of_links >= 1

        t1 = 10*(phi1_ref - phi1) + 1*(dphi1_ref - dphi1);
        
    end
    
    if number_of_links >= 2

        t2 = 30*(psi2_ref - psi2) + 3*(dpsi2_ref - dpsi2);
        
    end
    
    if number_of_links == 3

        t3 = 50*(th3_ref - th3) + 5*(dth3_ref - dth3);
        
    end    
    %----------------------------------------------------------------------
    
    %% Acceleration
    %----------------------------------------------------------------------
    ddphi1 = (eval(ddphi1_eqn));
    ddpsi2 = (eval(ddpsi2_eqn));
    ddth3 = (eval(ddth3_eqn));
    
    %----------------------------------------------------------------------
    
    %% Output
    %----------------------------------------------------------------------
    out_init = [dphi1; ddphi1; dpsi2; ddpsi2; dth3; ddth3];
    
    if (number_of_links <= 2)
        out = (eval(out_init(1:2*number_of_links,1)));
    else
        out = out_init(1:2*number_of_links,1);
    end
    %----------------------------------------------------------------------
end
