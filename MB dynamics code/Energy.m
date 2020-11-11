function [Ttot,Vtot,T1_,T2_,T3_,T,V1_,V2_,V3_,V,wrong_Ttot1,wrong_Ttot2,wrong_Ttot3] = Energy(mass,velocity1,velocity2,I,position,variables,wrong_vel,Rf_1_0,Rf_2_0,Rf_3_0)

    %% ENERGY
    
    % Default values and variables
    %----------------------------------------------------------------------
    global number_of_links;
    syms wrong_Ttot1 wrong_Ttot2 wrong_Ttot3 'real';
    T1_='';
    T2_='';
    T3_='';
    T='';
    V1_='';
    V2_= '';
    V3_ = '';
    V = '';
    %----------------------------------------------------------------------
    
    % Kinetic Energy
    %----------------------------------------------------------------------
    
    if number_of_links >= 1
        
        % translational Energy
        T1_t = 0.5*mass(1)*velocity1(:,1)'*velocity1(:,1); 
        
        % rotational Energy
        T1_r = 0.5*(velocity2(:,1))'*I(:,1:3)*velocity2(:,1); 
        
        % Total Kinetic energy of link 1
        T1 = T1_r + T1_t; 
        Ttot = simplify(T1,'IgnoreAnalyticConstraints',true);
        
        % Energies as text for GUI
        T1_ = replace(string(T1),[" - ","- "," + ",],["-","-","+"]);
        T = replace(string(Ttot),[" - ","- "," + ",],["-","-","+"]);
        
        %If student had wrong position
        wT1_t = 0.5*mass(1)*wrong_vel(:,1)'*wrong_vel(:,1);
        wT1 = T1_r + wT1_t;
        
        wrong_Ttot1 = wT1;
        
        % If student forgets to add rotational energy
        wrong_Ttot2 = T1_t;
        
        % If student forgets to add translational energy
        wrong_Ttot3 = T1_r;
    end
    
    if number_of_links >= 2
        % translational Energy
        T2_t = 0.5*mass(2)*(velocity1(:,2))'*velocity1(:,2);
        
        % rotational Energy
        T2_r = 0.5*(velocity2(:,2))'*I(:,4:6)*velocity2(:,2);
        
        % Kinetic energy of link 2
        T2 = T2_r + T2_t;
        
        % Total Kinetic energy
        Ttot = simplify(Ttot + T2,'IgnoreAnalyticConstraints',true);
        
        % Energies as text for GUI
        T2_ = replace(string(T2),[" - ","- "," + ",],["-","-","+"]);
        T = replace(string(Ttot),[" - ","- "," + ",],["-","-","+"]);
        
        %If student had wrong position
        wT2_t = 0.5*mass(1)*wrong_vel(:,1)'*wrong_vel(:,1);
        wT2 = T2_r + wT2_t;
        wTtot = wT2 + wrong_Ttot1;
        
        wrong_Ttot1 = wTtot;
        
        % If student forgets to add rotational energy
        wrong_Ttot2 = T1_t + T2_t;
        
        % If student forgets to add translational energy
        wrong_Ttot3 = T1_r + T2_r;
        
    end
    
    if number_of_links == 3
        
        % translational Energy
        T3_t = 0.5*mass(3)*(velocity1(:,3))'*velocity1(:,3);
        
        % rotational Energy
        T3_r = 0.5*(velocity2(:,3))'*I(:,7:9)*velocity2(:,3); 
        
        % Kinetic energy in link 3
        T3 = T3_t + T3_r;

        % Energies as text for GUI
        Ttot = simplify(Ttot + T3,'Steps',3);
        
        % Energies as text for GUI
        T3_ = replace(string(T3),[" - ","- "," + ",],["-","-","+"]);
        T = replace(string(Ttot),[" - ","- "," + ",],["-","-","+"]);
        
        %If student had wrong position
        wT3_t = 0.5*mass(1)*wrong_vel(:,1)'*wrong_vel(:,1);
        wTtot = T3_r + wT3_t + wrong_Ttot1;
        
        wrong_Ttot1 = wTtot;
        
        % If student forgets to add rotational energy
        wrong_Ttot2 = wrong_Ttot2 + T2_t + T3_t;
        
        % If student forgets to add translational energy
        wrong_Ttot3 = wrong_Ttot3 + T2_r + T3_r;
    end
    
    %----------------------------------------------------------------------

    % Potential Energy
    %----------------------------------------------------------------------

    if number_of_links >= 1
        % Potential Energy of link 1
        V1 = mass(1)*[0 0 variables(4)]*position(:,1); 
        
        % Total potential energy
        Vtot = V1;
        
        % Energies as text for GUI
        V1_ = replace(string(V1),[" - ","- "," + ",],["-","-","+"]);
        V = replace(string(Vtot),[" - ","- "," + ",],["-","-","+"]);
    end
    
    if number_of_links >= 2
        % Potential Energy of link 2
        V2 = mass(2)*[0 0 variables(4)]*position(:,2);
        
        % Total potential energy
        Vtot = Vtot + V2;
        
        % Energies as text for GUI
        V2_ = replace(string(V2),[" - ","- "," + ",],["-","-","+"]);
        V = replace(string(Vtot),[" - ","- "," + ",],["-","-","+"]);
    end
    
    if number_of_links == 3
        % Potential Energy of link 2
        V3 = mass(3)*[0 0 variables(4)]*position(:,3);
        
        % Total potential energy
        Vtot = Vtot + V3;
        
        % Energies as text for GUI
        V3_ = replace(string(V3),[" - ","- "," + ",],["-","-","+"]);
        V = replace(string(Vtot),[" - ","- "," + ",],["-","-","+"]);
    end

    %--------------------------------------------------------------------------
end