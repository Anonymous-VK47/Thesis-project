% This function pre-allocates the values to plot before the animation
% script is called
function [out,phi_Values,psi_Values,th_Values] = calculation(y,r)
    %% Variables and default conditions
    %----------------------------------------------------------------------
    global number_of_links;% 
    syms phi_Values psi_Values th_Values L1 L2 L3 'real';
    syms phi1 psi2 th3 'real';
    
    r1_0 = r(:,1);
    r2_0 = r(:,2);
    r3_0 = r(:,3);
    
    L1_= 0.3;
    L2_= 0.25;
    L3_= 0.15;
    %----------------------------------------------------------------------
    
    %% Substitution
    %----------------------------------------------------------------------
    
    % Substitute lengths and angles
    pos_link1 = subs(r1_0,[L1,phi1],[L1_,phi_Values]);
    
    if number_of_links >= 2
        
        pos_link2 = subs(r2_0,[L1,L2,phi1,psi2],[L1_,L2_,phi_Values,psi_Values]);
    end
    
    if number_of_links == 3
        
        pos_link3 = subs(r3_0,[L1,L2,L3,phi1,psi2,th3],[L1_,L2_,L3_,phi_Values,psi_Values,th_Values]);
    end
    %----------------------------------------------------------------------
    
    %% Position Vectors
    %----------------------------------------------------------------------
    % Vectorize equation
    [pos_link1_x,pos_link1_y,pos_link1_z] = turn_into_vector(pos_link1);
    
    % Here the simulation values will be used to animate the movement
    % First we create an array of the values of phi for the duration of the
    % simulation    
    phi_Values = y(:,1);
    
    % [x,y,z] positions
    pos_link1_x = eval(pos_link1_x);
    pos_link1_y = eval(pos_link1_y);
    pos_link1_z = eval(pos_link1_z);

    [pos_link1_x_, pos_link1_y_, pos_link1_z_] = size_check(pos_link1_x, pos_link1_y, pos_link1_z, phi_Values);
    out = [pos_link1_x_, pos_link1_y_, pos_link1_z_];

    if number_of_links >= 2
        
        % Vectorize equation
        [pos_link2_x,pos_link2_y,pos_link2_z] = turn_into_vector(pos_link2);
        
        % Here the simulation values will be used to animate the movement
        % First we create an array of the values of psi for the duration of the
        % simulation
        psi_Values = y(:,3);
        
        % [x,y,z] positions
        pos_link2_x = eval(pos_link2_x);
        pos_link2_y = eval(pos_link2_y);
        pos_link2_z = eval(pos_link2_z);
    
        [pos_link2_x_, pos_link2_y_, pos_link2_z_] = size_check(pos_link2_x, pos_link2_y, pos_link2_z, psi_Values);
        out = [pos_link1_x_, pos_link1_y_, pos_link1_z_,pos_link2_x_, pos_link2_y_, pos_link2_z_];
    end
    
    if number_of_links == 3

        % Vectorize equation
        [pos_link3_x,pos_link3_y,pos_link3_z] = turn_into_vector(pos_link3);
        
        % Here the simulation values will be used to animate the movement
        % First we create an array of the values of th for the duration of the
        % simulation
        th_Values = y(:,5);
        
        % [x,y,z] positions
        pos_link3_x = eval(pos_link3_x);
        pos_link3_y = eval(pos_link3_y);
        pos_link3_z = eval(pos_link3_z);
    
        [pos_link3_x_, pos_link3_y_, pos_link3_z_] = size_check(pos_link3_x, pos_link3_y, pos_link3_z, th_Values);
        
        out = [pos_link1_x_, pos_link1_y_, pos_link1_z_,pos_link2_x_, pos_link2_y_, pos_link2_z_,pos_link3_x_, pos_link3_y_, pos_link3_z_];

    end
   %-----------------------------------------------------------------------
end