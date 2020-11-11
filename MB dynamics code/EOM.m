function [ddphi1_eqn, ddpsi2_eqn, ddth3_eqn,M1,dM1,C1,G1,e1,e_,M,C,G,c1] = EOM(Ttot, Vtot, q, dq, ddq, mass, I, variables)
  
    %% Variables and default values
    %----------------------------------------------------------------------
    global number_of_links;
    syms ddphi1_eqn ddpsi2_eqn ddth3_eqn 'real';
    syms wrong_ddphi wrong_ddpsi wrong_ddth 'real';
    syms r1 r2 r3 v1 v2 v3 'real';
    
    e_ = '';
    
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
    
    %----------------------------------------------------------------------
    
    %% Mass Matrix
    %----------------------------------------------------------------------
    M = hessian(Ttot,dq);
    M = simplify(M,'Steps',3);
    assignin('base','M',M);

    % Text for GUI
    M1 = "["+ join(replace(string(M),[" - ","- "," + ",],["-","-","+"]))+ "]";
    %----------------------------------------------------------------------
    
    %% Derivative of Mass Matrix
    %----------------------------------------------------------------------
    dM = M_derivative(q,dq,M);
    dM = simplify(dM,'Steps',3);
    
    % Text for GUI
    dM1 = "["+ join(replace(string(dM),[" - ","- "," + ",],["-","-","+"]))+ "]";
    %----------------------------------------------------------------------
    
    %% C Matrix
    %----------------------------------------------------------------------
    C = dM * dq - jacobian(Ttot, q)';
    C = simplify(C,'Steps',3);
    assignin('base','C',C);
   
    % Text for GUI
    C1 = sprintf(replace(join(replace("["+string(C)+"]",[" - ","- "," + ",],["-","-","+"]))," ","\n"));
    %----------------------------------------------------------------------
    
    %% G Matrix 
    %----------------------------------------------------------------------
    G = simplify(transpose(jacobian(Vtot,q)),'Steps',3);
    assignin('base','G',G);
    
    % Text for GUI
    G1 = sprintf(replace(join(replace("["+string(G)+"]",[" - ","- "," + ",],["-","-","+"]))," ","\n"));
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
        eqn = simplify(M*ddq(1) + C + G - B,'Steps',5);
        u = simplify(M*v1 + C + G,'Steps',5);

        % Text for GUI
        c1 = string(u);
        %c1 = split(sprintf(replace(join(replace("["+string(u)+"]",[" - ","- "," + ",],["-","-","+"]))," ","\n")),"\n");
        e1 = join(string(eqn))+"=0";
    elseif number_of_links == 2
        eqn = simplify(M*ddq(1:2,:) + C + G - B,'Steps',5);
        u = simplify(M*[v1;v2] + C + G,'Steps',5);

        % Text for GUI
        c1 = [string(u(1));string(u(2))];
        %c1 = split(sprintf(replace(join(replace("["+string(u)+"]",[" - ","- "," + ",],["-","-","+"]))," ","\n")),"\n");
        e1 = split(sprintf(replace(join(replace("["+string(eqn)+"=0]",[" - ","- "," + ",],["-","-","+"]))," ","\n")),"\n");
    else
        eqn = simplify(M*ddq(1:3) + C + G - B,'Steps',5);
        u = simplify(M*[v1;v2;v3] + C + G,'Steps',5);
        
        % Text for GUI
        c1 = [string(u(1));string(u(2));string(u(3))];
        %c1 = split(sprintf(replace(join(replace("["+string(u)+"]",[" - ","- "," + ",],["-","-","+"]))," ","\n")),"\n");
        e1 = split(sprintf(replace(join(replace("["+string(eqn)+"=0]",[" - ","- "," + ",],["-","-","+"]))," ","\n")),"\n");
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
    
    %eqn = subs(eqn,[Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,Ix3,Iy3,Iz3],[Ix1_,Iy1_,Iz1_,Ix2_,Iy2_,Iz2_,Ix3_,Iy3_,Iz3_]);
    eqn = subs(eqn,[r1,r2,r3],[r1_,r2_,r3_]);
    eqn = subs(eqn, [g, m1, m2, m3, L1, L2, L3], [g_, m1_, m2_, m3_, L1_, L2_, L3_]);
    eqn = simplify(subs(eqn, g, g_),'Steps',5);
    %----------------------------------------------------------------------
    
    %% Solving for acceleration
    %----------------------------------------------------------------------
    
    sol = solve(eqn, ddq);    
    
    if number_of_links == 1
        ddphi1_eqn = simplify(sol(1),'Steps',5);

        % Text for GUI
        e_ = join(string(ddphi1_eqn));
    
    end
    
    if number_of_links >= 2

        ddphi1_eqn = simplify(sol.ddphi1,'Steps',5);
        ddpsi2_eqn = simplify(sol.ddpsi2,'Steps',5);
        
        % File export for testing
        acc1 = fopen('C:\Multibody Dynamics Generator\acc1.txt','wt');
        fprintf(acc1,'%s',join(string(ddphi1_eqn)));
        fclose(acc1);
        
        acc2 = fopen('C:\Multibody Dynamics Generator\acc2.txt','wt');
        fprintf(acc2,'%s',join(string(ddpsi2_eqn)));
        fclose(acc2);
        
        % Text for GUI
        e_ = [replace(string(ddphi1_eqn),[" - ","- "," + ",],["-","-","+"]);replace(string(ddpsi2_eqn),[" - ","- "," + ",],["-","-","+"])];
        
    end
    
    if number_of_links == 3
        ddth3_eqn = simplify(sol.ddth3,'Steps',5);
        
        % File export for testing
        acc3 = fopen('C:\Multibody Dynamics Generator\acc3.txt','wt');
        fprintf(acc3,'%s',join(string(ddth3_eqn)));
        fclose(acc3);
 
        % Text for GUI
        e_ = [replace(string(ddphi1_eqn),[" - ","- "," + ",],["-","-","+"]);replace(string(ddpsi2_eqn),[" - ","- "," + ",],["-","-","+"]);replace(string(ddth3_eqn),[" - ","- "," + ",],["-","-","+"])];
        
    end
 
    %----------------------------------------------------------------------
end