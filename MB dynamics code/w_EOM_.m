function [w_EOM,cw1,cw2,cw3] = w_EOM_(Ttot, Vtot, q, dq, ddq, mass, I, variables,wrong_Ttot1,wrong_Ttot2)
    %% Variables and default values
    %----------------------------------------------------------------------
    global number_of_links;
    syms ddphi1_eqn ddpsi2_eqn ddth3_eqn 'real';
    syms wrong_ddphi wrong_ddpsi wrong_ddth 'real';
    syms r1 r2 r3 u v1 v2 v3 'real';
    
    w_EOM = '';

    m1 = mass(1);
    m2 = mass(2);
    m3 = mass(3);
    
    L1 = variables(1);
    L2 = variables(2);
    L3 = variables(3);
    
    g = variables(4);
    t1 = variables(5);
    t2 = variables(6);
    t3 = variables(7);
    
    %Ix1 = I(1,1);
    %Iy1 = I(2,2);
    %Iz1 = I(3,3);
    
    %Ix2 = I(1,4);
    %Iy2 = I(2,5);
    %Iz2 = I(3,6);
    
    %Ix3 = I(1,7);
    %Iy3 = I(2,8);
    %Iz3 = I(3,9);
    %----------------------------------------------------------------------
    
    %% Mass Matrix
    %----------------------------------------------------------------------
    % Check dimensions of matrices before applying hessian
    [wrong_M1,wrong_M2] = check_M(hessian(wrong_Ttot1,dq),hessian(wrong_Ttot2,dq));
    
    wrong_M1 = simplify(wrong_M1,'IgnoreAnalyticConstraints',true);
    wrong_M2 = simplify(wrong_M2,'IgnoreAnalyticConstraints',true);
    
    M = hessian(Ttot,dq);
    %----------------------------------------------------------------------
    
    %% Derivative of Mass Matrix
    %----------------------------------------------------------------------
    wrong_dM1 = simplify(M_derivative(q,dq,wrong_M1),'IgnoreAnalyticConstraints',true);
    wrong_dM2 = simplify(M_derivative(q,dq,wrong_M2),'IgnoreAnalyticConstraints',true);

    dM = M_derivative(q,dq,M);
    %----------------------------------------------------------------------
    
    %% C Matrix
    %----------------------------------------------------------------------
    wrong_C1 = simplify(wrong_dM1* dq - jacobian(wrong_Ttot1,q)','IgnoreAnalyticConstraints',true);
    wrong_C2 = simplify(wrong_dM2 * dq - jacobian(wrong_Ttot2,q)','IgnoreAnalyticConstraints',true);

    C = dM * dq - jacobian(Ttot, q)';
    %----------------------------------------------------------------------
    
    %% G Matrix 
    %----------------------------------------------------------------------
    G = transpose(jacobian(Vtot,q));
    %----------------------------------------------------------------------
    
    %% B input matrix 
    %----------------------------------------------------------------------
    if number_of_links == 1
        B = t1;
    elseif number_of_links == 2
        B = [t1; t2];
    else
        B = [t1; t2; t3];
    end
    %----------------------------------------------------------------------
    
    %% EOM
    %----------------------------------------------------------------------
    if number_of_links == 1
        eqn = M*ddq(1) + C + G - B;
        wrong_eqn1 = simplify(wrong_M1*ddq(1) + wrong_C1 + G - B,'IgnoreAnalyticConstraints',true);
        wrong_eqn2 = simplify(wrong_M2*ddq(1) + wrong_C2 + G - B,'IgnoreAnalyticConstraints',true);
        wrong_eqn3 = simplify(eqn + randi(20),'IgnoreAnalyticConstraints',true);
        
        cw1 = simplify(wrong_M1*v1 + wrong_C1 + G ,'IgnoreAnalyticConstraints',true);
        cw2 = simplify(wrong_M2*v1 + wrong_C2 + G ,'IgnoreAnalyticConstraints',true);
        cw3 = simplify(M*v1 + C + G + randi(20),'IgnoreAnalyticConstraints',true);

    elseif number_of_links == 2
        eqn = M*ddq(1:2,:) + C + G - B;
        wrong_eqn1 = simplify(wrong_M1*ddq(1:2) + wrong_C1 + G - B,'IgnoreAnalyticConstraints',true);
        wrong_eqn2 = simplify(wrong_M2*ddq(1:2) + wrong_C2 + G - B,'IgnoreAnalyticConstraints',true);
        wrong_eqn3 = simplify(eqn + randi(20),'IgnoreAnalyticConstraints',true);
        
        cw1 = simplify(wrong_M1*[v1;v2] + wrong_C1 + G ,'IgnoreAnalyticConstraints',true);
        cw2 = simplify(wrong_M2*[v1;v2] + wrong_C2 + G ,'IgnoreAnalyticConstraints',true);
        cw3 = simplify(M*[v1;v2] + C + G + randi(20),'IgnoreAnalyticConstraints',true);

    else
        eqn = M*ddq(1:3) + C + G - B;
        wrong_eqn1 = simplify(wrong_M1*ddq(1:3) + wrong_C1 + G - B,'IgnoreAnalyticConstraints',true);
        wrong_eqn2 = simplify(wrong_M2*ddq(1:3) + wrong_C2 + G - B,'IgnoreAnalyticConstraints',true);
        wrong_eqn3 = simplify(eqn + randi(20),'IgnoreAnalyticConstraints',true);
        
        cw1 = simplify(wrong_M1*[v1;v2;v3] + wrong_C1 + G ,'IgnoreAnalyticConstraints',true);
        cw2 = simplify(wrong_M2*[v1;v2;v3] + wrong_C2 + G ,'IgnoreAnalyticConstraints',true);
        cw3 = simplify(M*[v1;v2;v3] + C + G + randi(20),'IgnoreAnalyticConstraints',true);
        
    end
    %----------------------------------------------------------------------
    
    %% Substitution
    %----------------------------------------------------------------------
    m1_ = 3;
    m2_ = 2;
    m3_ = 1;
    L1_ = 0.3;
    L2_ = 0.25;
    L3_ = 0.15;
    g_ = 9.81;
    r1_ = 0.3;
    r2_ = 0.25;
    r3_ = 0.15;
    
    %Ix1_ = 1/12*m1*(3*r1_^2+L1^2);
    %Iy1_ = Ix1_;
    %Iz1_ = (1/2)*m1*r1_^2;
    
    %Ix2_ = 1/12*m1*(3*r2_^2+L1^2);
    %Iy2_ = (1/2)*m1*r2_^2;
    %Iz2_ = Ix2_; 
    
    %Ix3_ = 1/12*m1*(3*r3_^2+L1^2);
    %Iy3_ = Ix3_;
    %Iz3_ = (1/2)*m1*r3_^2;

    %wrong_eqn1 = subs(wrong_eqn1,[Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,Ix3,Iy3,Iz3,g, m1, m2, m3, L1, L2, L3],[Ix1_,Iy1_,Iz1_,Ix2_,Iy2_,Iz2_,Ix3_,Iy3_,Iz3_,g_, m1_, m2_, m3_, L1_, L2_, L3_]);
    wrong_eqn1 = subs(wrong_eqn1,[r1,r2,r3,g, m1, m2, m3, L1, L2, L3],[r1_,r2_,r3_,g_, m1_, m2_, m3_, L1_, L2_, L3_]);
    wrong_eqn1 = simplify(wrong_eqn1,'IgnoreAnalyticConstraints',true);

    %wrong_eqn2 = subs(wrong_eqn2,[Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,Ix3,Iy3,Iz3,g, m1, m2, m3, L1, L2, L3],[Ix1_,Iy1_,Iz1_,Ix2_,Iy2_,Iz2_,Ix3_,Iy3_,Iz3_,g_, m1_, m2_, m3_, L1_, L2_, L3_]);
    wrong_eqn2 = subs(wrong_eqn2,[r1,r2,r3,g, m1, m2, m3, L1, L2, L3],[r1_,r2_,r3_,g_, m1_, m2_, m3_, L1_, L2_, L3_]);
    wrong_eqn2 = simplify(wrong_eqn2,'IgnoreAnalyticConstraints',true);

    %wrong_eqn3 = subs(wrong_eqn3,[Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,Ix3,Iy3,Iz3,g, m1, m2, m3, L1, L2, L3],[Ix1_,Iy1_,Iz1_,Ix2_,Iy2_,Iz2_,Ix3_,Iy3_,Iz3_,g_, m1_, m2_, m3_, L1_, L2_, L3_]);
    wrong_eqn3 = subs(wrong_eqn3,[r1,r2,r3,g, m1, m2, m3, L1, L2, L3],[r1_,r2_,r3_,g_, m1_, m2_, m3_, L1_, L2_, L3_]);
    wrong_eqn3 = simplify(wrong_eqn3,'IgnoreAnalyticConstraints',true);
    %----------------------------------------------------------------------
    
    %% Solving for acceleration
    %----------------------------------------------------------------------

    wrong_sol1 = solve(wrong_eqn1, ddq);
    wrong_sol2 = solve(wrong_eqn2, ddq);
    wrong_sol3 = solve(wrong_eqn3, ddq);
 
    if number_of_links == 1
        wrong_ddphi1 = simplify(wrong_sol1,'IgnoreAnalyticConstraints',true);
        wrong_ddphi2 = simplify(wrong_sol2,'IgnoreAnalyticConstraints',true);
        wrong_ddphi3 = simplify(wrong_sol3,'IgnoreAnalyticConstraints',true);

        % Text for GUI
        w_EOM = [replace(string(wrong_ddphi1),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddphi2),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddphi3),[" - ","- "," + ",],["-","-","+"])];
    
    end
    
    if number_of_links >= 2
        wrong_ddphi1 = simplify(wrong_sol1.ddphi1,'IgnoreAnalyticConstraints',true);
        wrong_ddphi2 = simplify(wrong_sol2.ddphi1,'IgnoreAnalyticConstraints',true);
        wrong_ddphi3 = simplify(wrong_sol3.ddphi1,'IgnoreAnalyticConstraints',true);
        
        wrong_ddpsi1 = simplify(wrong_sol1.ddpsi2,'IgnoreAnalyticConstraints',true);
        wrong_ddpsi2 = simplify(wrong_sol2.ddpsi2,'IgnoreAnalyticConstraints',true);
        wrong_ddpsi3 = simplify(wrong_sol3.ddpsi2,'IgnoreAnalyticConstraints',true);
        
        % Text for GUI
        w_EOM = [replace(string(wrong_ddphi1),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddphi2),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddphi3),[" - ","- "," + ",],["-","-","+"]);
                 replace(string(wrong_ddpsi1),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddpsi2),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddpsi3),[" - ","- "," + ",],["-","-","+"])];
             
    end
    
    if number_of_links == 3
        wrong_ddth1 = simplify(wrong_sol1.ddth3,'IgnoreAnalyticConstraints',true);
        wrong_ddth2 = simplify(wrong_sol2.ddth3,'IgnoreAnalyticConstraints',true);
        wrong_ddth3 = simplify(wrong_sol3.ddth3,'IgnoreAnalyticConstraints',true);

        % Text for GUI
        w_EOM = [replace(string(wrong_ddphi1),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddphi2),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddphi3),[" - ","- "," + ",],["-","-","+"]);
                 replace(string(wrong_ddpsi1),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddpsi2),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddpsi3),[" - ","- "," + ",],["-","-","+"]);
                 replace(string(wrong_ddth1),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddth2),[" - ","- "," + ",],["-","-","+"]),replace(string(wrong_ddth3),[" - ","- "," + ",],["-","-","+"])];
        
    end
    %----------------------------------------------------------------------
end