function links(app,number_of_links)
    
    %% Variables to be used
    %----------------------------------------------------------------------
    global ddpsi2_eqn ddphi1_eqn ddth3_eqn;
    global number_of_links l_n;
    
    syms g m1 m2 m3 'real'       
    syms Ix1 Iy1 Iz1 Ix2 Iy2 Iz2 Ix3 Iy3 Iz3 'real' 
    syms phi1 psi2 th3 'real'     
    syms dphi1 dpsi2 dth3 'real'   
    syms ddphi1 ddpsi2 ddth3 ddpsi2_ ddphi1_ ddth3_ 'real'
    syms L1 L2 L3 r1 r2 r3 'real'         
    syms t1 t2 t3 'real'
    syms phi_Values psi_Values th_Values 'real'
    %----------------------------------------------------------------------
    %% Misc 
    %----------------------------------------------------------------------
    %number_of_links = 2;
    
    % Intial Values 
    Init = [0; 0; 0; 0; 0; 0];
    Init = Init(1:2*number_of_links,1);
    
    radii = [r1,r2,r3];
    
    mass = [m1 m2 m3];
    variables = [L1 L2 L3 g t1 t2 t3];
    %----------------------------------------------------------------------
    
    %% Generalised Co-ordinates
    %----------------------------------------------------------------------
    q = [phi1; psi2;th3];
    q = q(1:number_of_links,1);
    
    dq = [dphi1; dpsi2; dth3];
    dq = dq(1:number_of_links,1);

    ddq = [ddphi1; ddpsi2; ddth3];
    ddq = ddq(1:number_of_links,1);
    %----------------------------------------------------------------------
    
    %% Rotation function calls
    %----------------------------------------------------------------------
    % Call function with rotation matrices
    [RotX,RotY,RotZ] = Rotation;

    % Randomize rotations
    [Rot1,Rot2,Rot3,out] = Rot_random(RotX, RotY, RotZ);
    
    % Call kinematics solver
    [Rf_1_0,Rf_2_0,Rf_3_0,velocity1,velocity2,position,R1,R2,R3,R4,R5,rp1,rp2,rp3,rR1,rR2,rR3,drR1,drR2,drR3,dR1,dR2,dR3,S1,v1,S2,v2,S3,v3,wrong_vel,l_n] = Kinematics(variables,q,dq,Rot1,Rot2,Rot3,RotX,RotY,RotZ);
    %----------------------------------------------------------------------
    
    %% Inertia
    
    % Inertias for rotation about each axis in body frame
    % Treat links as cylinders
    
    I11 = diag([Ix1,Iy1,Iz1]);
    
    I1_0 = simplify([Ix1 0 0;
                    0 Iy1 0;
                    0 0 Iz1]*Rf_1_0);
                
    I22 = diag([Ix2,Iy2,Iz2]);
   
    I2_0 = simplify([Ix2 0 0;
                    0 Iy2 0;
                    0 0 Iz2]*Rf_2_0);
                
    I33 = diag([Ix3,Iy3,Iz3]);
      
    I3_0 = simplify([Ix3 0 0;
                    0 Iy3 0;
                    0 0 Iz3]*Rf_3_0);

    
    Ix1 = 1/12*m1*(3*r1^2+L1^2);
    Iy1 = Ix1;
    Iz1 = (1/2)*m1*r1^2;
    
    Ix2 = 1/12*m1*(3*r2^2+L2^2);
    Iy2 = (1/2)*m1*r2^2;
    Iz2 = Ix2; 
    
    Ix3 = 1/12*m1*(3*r3^2+L3^2);
    Iy3 = Ix3;
    Iz3 = (1/2)*m1*r3^2;
    
    I1 = diag([Ix1,Iy1,Iz1]);
    I2 = diag([Ix2,Iy2,Iz2]);
    I3 = diag([Ix3,Iy3,Iz3]);
    
    I = [I1 I2 I3];
    %----------------------------------------------------------------------
    
    %% GUI configs part 1
    %----------------------------------------------------------------------
    % GUI Text
    q1 = "[" + replace(join(string(q))," ",", ") + "]";
    m = "[" + replace(join(string(mass(1,1:number_of_links)))," ",", ") + "]";
    L = "[" + replace(join(string(variables(1,1:number_of_links)))," ",", ") + "]";
    t = "[" + replace(join(string(variables(1,5:4+number_of_links)))," ",", ") + "]";
    
    In1 = "["+ join(string(I11))+ "]";
    In2 = "["+ join(string(I22))+ "]";
    In3 = "["+ join(string(I33))+ "]";
    
    In1_0 = "["+ join(string(I1_0))+ "]";
    In2_0 = "["+ join(string(I2_0))+ "]";
    In3_0 = "["+ join(string(I3_0))+ "]";
    
    dq2 =  create_A_string1(dq);
    ddq3 = create_A_string1(ddq);
    
    if number_of_links == 1
        Link1_text = sprintf("Consider the link as a cylinders with radius, length and mass equal to [%s] (m), %s (m) and %s (kg) respectively. \n\nThe link rotates about the [%s] axis, along the angle %s (degrees). \n\nThe Torque %s (Nm) acts as the input.",radii(1), L , m, out(1), q1, t);
        assignin('base','Link1_text',Link1_text);
    elseif number_of_links == 2
        Link2_text = sprintf("Consider the links as cylinders with radii, lengths and masses equal to [%s,%s] (m), %s (m) and %s (kg) respectively. \n\nLinks 1 and 2 rotate about the [%s,%s] axes along angles %s (degrees) respectively. \n\nTorques %s (Nm) act as inputs to links 1 and 2 respectively.",radii(1:2), L, m, out(1:2), q1, t);
        assignin('base','Link2_text',Link2_text);
    else
        Link3_text = sprintf("Consider the links as cylinders with radii, lengths and masses equal to [%s,%s,%s] (m), %s (m) and %s (kg) respectively. \n\nLinks 1, 2 and 3 rotate about the [%s,%s,%s] axes along angles %s (degrees) respectively. \n\nTorques %s (Nm) act as inputs to links 1,2 and 3 respectively.",radii(1:3), L, m, out(1:3), q1, t);
        assignin('base','Link3_text',Link3_text);
    end
    %----------------------------------------------------------------------
    
    %% Function Calls part 2
    %----------------------------------------------------------------------
    
    % Call energy solver
    [Ttot,Vtot,T1_,T2_,T3_,T,V1_,V2_,V3_,V,wrong_Ttot1,wrong_Ttot2,wrong_Ttot3] = Energy(mass,velocity1,velocity2,I,position,variables,wrong_vel,Rf_1_0,Rf_2_0,Rf_3_0);
    
    % Call EOM solver
    [ddphi1_eqn, ddpsi2_eqn, ddth3_eqn,M1,dM1,C1,G1,e1,e_,M,C,G,c1] = EOM(Ttot,Vtot,q,dq,ddq,mass,I,variables);
    
    % call ode solver (Simulation/Integration)
    [~,y] = simulate_by_integration(@Controller, [0; 1.5], Init);
    %yc = FBCtrl(M,C,G);
    %y = yc;
    %----------------------------------------------------------------------
    %toc
    
    %% GUI text
    %%=--------------------------------------------------------------------
    if number_of_links == 1
        app.TextArea_5.Value = evalin('base', 'Link1_text');
    elseif number_of_links == 2
        app.TextArea_5.Value = evalin('base', 'Link2_text');
    else
        app.TextArea_5.Value = evalin('base', 'Link3_text');
    end
    drawnow();
    %----------------------------------------------------------------------
    
    %% Animation
    %----------------------------------------------------------------------
    
    % Obtain values to plot
    [out,psi_Values]  = calculation(y,position);
    
    %Animate
    animation(out,psi_Values);
    %----------------------------------------------------------------------
    
    %% Gui configs part 2
    %----------------------------------------------------------------------
    % Call function with wrong EOM (For multiple choice)
    [w_EOM,cw1,cw2,cw3] = w_EOM_(Ttot, Vtot, q, dq, ddq, mass, I, variables,wrong_Ttot1,wrong_Ttot2);
    
    % Randomize multiple choice output
    C0 = {"["+c1+"]","["+string(cw1)+"]","["+string(cw2)+"]","["+string(cw3)+"]"};
    C0 = C0(randperm(numel(C0)));
    
    A0 = {{string(wrong_Ttot1);string(V)},{string(wrong_Ttot2);string(V)},{string(wrong_Ttot3);string(V)},{string(Ttot);string(V)}};
    A0 = A0(randperm(numel(A0)));
    
    A = [w_EOM,e_];
    A = A(randperm(numel(A)));
    
    A1 = {w_EOM(:,1),w_EOM(:,2),w_EOM(:,3),e_};
    A1 = A1(randperm(numel(A1)));
    
    % Find position of correct answer 
    Energy_index = find_index(A0,{string(Ttot);string(V)});
    assignin('base','Energy_index',Energy_index);
    
    Control_index = find_index(C0,"["+string(c1)+"]");
    assignin('base','Control_index',Control_index);
    
    if number_of_links == 1
        EOM_index = find_index(A,e_);
    else
        EOM_index = find_index(A1,e_);
    end
    
    assignin('base','EOM_index',EOM_index);
    
    % Multiple choice text (what user sees)
    M_choice_Energy = sprintf("Choose the correct pair of kinetic and potential energies:\n\nA.\nT = %s\n\nV = %s \n\nB.\nT = %s\n\nV = %s \n\nC.\nT = %s\n\nV = %s \n\nD.\nT = %s\n\nV = %s ",A0{1}{:},A0{2}{:},A0{3}{:},A0{4}{:});
    assignin('base','M_choice_Energy',M_choice_Energy);
    
    % Multiple choice text (what user sees)
    if number_of_links == 1
        M_choice_Control = sprintf("Choose the correct pseudo control input for feedback linearization:\n\nA.\nu = \n%s\n\nB.\nu = \n%s\n\nC.\nu = %s\n\nD.\nu = %s\n",C0{1},C0{2},C0{3},C0{4});
    elseif number_of_links == 2
        M_choice_Control = sprintf("Choose the correct pseudo control input for feedback linearization:\n\nA.\nu = \n%s\n%s\n\nB.\nu = \n%s\n%s\n\nC.\nu = \n%s\n%s\n\nD.\nu =\n%s\n%s\n",C0{1},C0{2},C0{3},C0{4});
    else
        M_choice_Control = sprintf("Choose the correct pseudo control input for feedback linearization:\n\nA.\nu = \n%s\n%s\n%s\n\nB.\nu = \n%s\n%s\n%s\n\nC.\nu = \n%s\n%s\n%s\n\nD.\nu = \n%s\n%s\n%s\n",C0{1},C0{2},C0{3},C0{4});
    end
        
    assignin('base','M_choice_Control',M_choice_Control);
    
    %----------------------------------------------------------------------
    
    %% GUI text
    %----------------------------------------------------------------------
    
    if number_of_links == 1
        EOM_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\n\nPosition: \n\nP1 = %s (m)\n\nRotating Link 1 back to inertial frame:  Given by (P1*R1_0)\n\nP1_0 = %s (m)\n\n\nInertia:\n\nI1 = \n%s\n%s\n%s (kg.m^2) \n\nRotating back to inertial frame:\n\nI1_0 = \n%s\n%s\n%s (kg.m^2) \n\nIx1 = %s (kg.m^2)\nIy1 = %s (kg.m^2)\nIz1 = %s (kg.m^2)\n\n\nLinear Velocity: Given by (dP/dt)\n\nV1 = %s (m/s)\n\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\n\nS(ω1) = \n%s\n%s\n%s \n\n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\n\nω1 = %s (rad/s) \n\n\nEnergy: \nGiven by (E = T + V)\n\nKinetic Energy: Given by (T = 0.5*m*dP^2 + 0.5*I* ω^2)\n\nT = %s (J)\n\n\nPotential Energy: Given by (V = m*g*h)\n\nV = %s (J)\n\n\nNow to determine the Equations of motion use the equation:\n\n𝑄= M(𝑞)∙ddq+𝐶(𝑞,dq)+𝐺(𝑞)\n\n\nWhere Q represents the input torques:\n\nQ = %s (Nm)\n\n\nM is the mass matrix given by ( Mij(q) = ∂T/(∂ (dqi)*d(∂qj)) ):\n\nM = %s \n\n\nC contains the centrifugal and coriolis accelerations given by ( C(q,dq) = (dM/dt)*dq – ∂T/∂q ):\n\ndM = %s \n\nC = %s \n\n\nG is the matrix that contains the potential energy given by: ( G(q) = ∂V/∂q )\nG = %s \n\n\nNow the EOM:\n\n%s \n\nSolving for acceleration:\n\nddphi = %s (rad/s^2)",q1,dq2,ddq3,R1,rp1,rR1,In1,In1_0,Ix1,Iy1,Iz1,drR1,dR1,S1,v1,T,V,t,M1,dM1,C1,G1,e1,e_);
        
        Control_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\n\nPosition: \n\nP1 = %s (m)\n\nRotating Link 1 back to inertial frame:  Given by (P1*R1_0)\n\nP1_0 = %s (m)\n\n\nInertia:\n\nI1 = \n%s\n%s\n%s (kg.m^2) \n\nRotating back to inertial frame:\n\nI1_0 = \n%s\n%s\n%s (kg.m^2) \n\nIx1 = %s (kg.m^2)\nIy1 = %s (kg.m^2)\nIz1 = %s (kg.m^2)\n\n\nLinear Velocity: Given by (dP/dt)\n\nV1 = %s (m/s)\n\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\n\nS(ω1) = \n%s\n%s\n%s \n\n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\n\nω1 = %s (rad/s) \n\n\nEnergy: \nGiven by (E = T + V)\n\nKinetic Energy: Given by (T = 0.5*m*dP^2 + 0.5*I* ω^2)\n\nT = %s (J)\n\n\nPotential Energy: Given by (V = m*g*h)\n\nV = %s (J)\n\n\nNow to determine the Equations of motion use the equation:\n\n𝑄= M(𝑞)∙ddq+𝐶(𝑞,dq)+𝐺(𝑞)\n\n\nWhere Q represents the input torques:\n\nQ = %s (Nm)\n\n\nM is the mass matrix given by ( Mij(q) = ∂T/(∂ (dqi)*d(∂qj)) ):\n\nM = %s \n\n\nC contains the centrifugal and coriolis accelerations given by ( C(q,dq) = (dM/dt)*dq – ∂T/∂q ):\n\ndM = %s \n\nC = %s \n\n\nG is the matrix that contains the potential energy given by: ( G(q) = ∂V/∂q )\nG = %s \n\n\nNow the EOM:\n\n%s \n\nFor the pseudo control input, substitute ddphi1 = v1 and u = t1, then solve for u:\n\nu = %s (Nm)",q1,dq2,ddq3,R1,rp1,rR1,In1,In1_0,Ix1,Iy1,Iz1,drR1,dR1,S1,v1,T,V,t,M1,dM1,C1,G1,e1,c1);
        
        M_choice_EOM = sprintf("A.\n\nddphi1 = %s\n\nB.\n\nddphi1 = %s\n\nC.\n\nddphi1 = %s\n\nD.\n\nddphi1 = %s",A);
        
        Energy_sol = sprintf("Solution:\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nPosition: \n\nP1 = %s (m)\n\nRotating Link 1 back to inertial frame:  Given by (P1*R1_0)\n\nP1_0 = %s (m)\n\nInertia:\n\nI1 = \n%s\n%s\n%s (kg.m^2) \n\nRotating back to inertial frame:\n\nI1_0 = \n%s\n%s\n%s (kg.m^2) \n\nIx1 = %s (kg.m^2)\nIy1 = %s (kg.m^2)\nIz1 = %s (kg.m^2)\n\nVelocity: \n\nLinear Velocity: Given by (dP/dt)\n\nV1 = %s (m/s)\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\nω1 = %s (rad/s)\n\nEnergy: \nGiven by (E = T + V)\n\nKinetic Energy: Given by (T = 0.5*m*dP^2 + 0.5*I* ω^2)\n\nT = %s (J)\n\n\nPotential Energy: Given by (V = m*g*h)\n\nV = %s (J)",q1,dq2,ddq3,R1,rp1,rR1,In1,In1_0,Ix1,Iy1,Iz1,drR1,dR1,S1,v1,T,V);

        av_sol = sprintf("Solution:\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\ndR1_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\nω1 = %s (rad/s)",q1,dq2,ddq3,R1,dR1,S1,v1);
        
    elseif number_of_links == 2
        EOM_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 1:\n\nR2_1 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 0:\n\nR2_0 = R1_0*R2_1\n\nR2_0 = \n%s\n%s\n%s \n\n\nPosition: \n\nP1 = %s (m), P2 = %s (m) \n\nRotating Link 1 back to inertial frame:  Given by (P1*R1_0)\n\nP1_0 = %s (m)\n\nRotating Link 2 back to inertial frame:  Given by (P2*R2_0)\n\nP2_0 = %s (m)\n\n\nInertia:\n\nI1 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI1_0 = \n%s\n%s\n%s (kg.m^2) \n\nI2 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI2_0 = \n%s\n%s\n%s (kg.m^2) \n\nIx1 = %s (kg.m^2)\nIy1 = %s (kg.m^2)\nIz1 = %s (kg.m^2)\n\nIx2 = %s (kg.m^2)\nIy2 = %s (kg.m^2)\nIz2 = %s \n\n\nLinear Velocity: Given by (dP/dt)\nV1 = %s (m/s) \n\nV2 = %s (m/s)\n\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\ndR2_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\n\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\n\nω1 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 2:  Given by ( S(ω2) = transpose(R2_0) * dR2_0 )\n\nS(ω2) = \n%s\n%s\n%s \n\nAngular Velocity of Link 2:  Given by ω2 = [S(ω2)(3,2) , S(ω2)(1,3) , S(ω2)(2,1)]\n\nω2 = %s (rad/s)\n\n\nEnergy: \n\nGiven by (E = T + V)\n\n\nKinetic Energy: Given by (T = 0.5*m*dP^2 + 0.5*I* ω^2)\n\nT1 = %s (J)\nT2 = %s (J)\n\nT = T1 + T2 \n\nT = %s (J)\n\n\nPotential Energy: Given by (V = m*g*h)\n\nV1 = %s (J)\nV2 = %s (J)\n\nV = V1+V2\n\nV = %s (J)\n\n\nNow to determine the Equations of motion use the equation:\n\n𝑄= M(𝑞)∙ddq+𝐶(𝑞,dq)+𝐺(𝑞)\n\n\nWhere Q represents the input torques:\n\nQ = %s (Nm)\n\n\nM is the mass matrix given by ( Mij(q) = ∂T/(∂ (dqi)*d(∂qj)) ):\n\nM = \n%s\n%s \n\n\nC contains the centrifugal and coriolis accelerations given by ( C(q,dq) = (dM/dt)*dq – ∂T/∂q ):\n\ndM = \n%s\n%s \n\nC = \n%s \n\n\nG is the matrix that contains the potential energy given by: ( G(q) = ∂V/∂q )\nG = \n%s \n\n\nNow the EOM:\n\n%s \n\nSolving for acceleration:\n\nddphi = %s (rad/s^2)\n\nddpsi = %s (rad/s^2)",q1,dq2,ddq3,R1,R2,R3,rp1,rp2,rR1,rR2,In1,In1_0,In2,In2_0,Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,drR1,drR2,dR1,dR2,S1,v1,S2,v2,T1_,T2_,T,V1_,V2_,V,t,M1,dM1,C1,G1,e1,e_);

        Control_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 1:\n\nR2_1 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 0:\n\nR2_0 = R1_0*R2_1\n\nR2_0 = \n%s\n%s\n%s \n\n\nPosition: \n\nP1 = %s (m), P2 = %s (m) \n\nRotating Link 1 back to inertial frame:  Given by (P1*R1_0)\n\nP1_0 = %s (m)\n\nRotating Link 2 back to inertial frame:  Given by (P2*R2_0)\n\nP2_0 = %s (m)\n\n\nInertia:\n\nI1 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI1_0 = \n%s\n%s\n%s (kg.m^2) \n\nI2 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI2_0 = \n%s\n%s\n%s (kg.m^2) \n\nIx1 = %s (kg.m^2)\nIy1 = %s (kg.m^2)\nIz1 = %s (kg.m^2)\n\nIx2 = %s (kg.m^2)\nIy2 = %s (kg.m^2)\nIz2 = %s \n\n\nLinear Velocity: Given by (dP/dt)\nV1 = %s (m/s) \n\nV2 = %s (m/s)\n\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\ndR2_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\n\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\n\nω1 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 2:  Given by ( S(ω2) = transpose(R2_0) * dR2_0 )\n\nS(ω2) = \n%s\n%s\n%s \n\nAngular Velocity of Link 2:  Given by ω2 = [S(ω2)(3,2) , S(ω2)(1,3) , S(ω2)(2,1)]\n\nω2 = %s (rad/s)\n\n\nEnergy: \n\nGiven by (E = T + V)\n\n\nKinetic Energy: Given by (T = 0.5*m*dP^2 + 0.5*I* ω^2)\n\nT1 = %s (J)\nT2 = %s (J)\n\nT = T1 + T2 \n\nT = %s (J)\n\n\nPotential Energy: Given by (V = m*g*h)\n\nV1 = %s (J)\nV2 = %s (J)\n\nV = V1+V2\n\nV = %s (J)\n\n\nNow to determine the Equations of motion use the equation:\n\n𝑄= M(𝑞)∙ddq+𝐶(𝑞,dq)+𝐺(𝑞)\n\n\nWhere Q represents the input torques:\n\nQ = %s (Nm)\n\n\nM is the mass matrix given by ( Mij(q) = ∂T/(∂ (dqi)*d(∂qj)) ):\n\nM = \n%s\n%s \n\n\nC contains the centrifugal and coriolis accelerations given by ( C(q,dq) = (dM/dt)*dq – ∂T/∂q ):\n\ndM = \n%s\n%s \n\nC = \n%s \n\n\nG is the matrix that contains the potential energy given by: ( G(q) = ∂V/∂q )\nG = \n%s \n\n\nNow the EOM:\n\n%s \n\nFor the pseudo control input, substitute ddphi1 = v1, ddpsi2 = v2 and u = [t1, t2] (u will be a vector of inputs), then solve for u:\n\nu = %s\n\n%s (Nm))",q1,dq2,ddq3,R1,R2,R3,rp1,rp2,rR1,rR2,In1,In1_0,In2,In2_0,Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,drR1,drR2,dR1,dR2,S1,v1,S2,v2,T1_,T2_,T,V1_,V2_,V,t,M1,dM1,C1,G1,e1,c1);
        
        M_choice_EOM = sprintf("A.\n\nddphi1 = %s\n\nddpsi2 = %s\n\nB.\n\nddphi1 = %s\n\nddpsi2 = %s\n\nC.\n\nddphi1 = %s\n\nddpsi2 = %s\n\nD.\n\nddphi1 = %s\n\nddpsi2 = %s",A1{1},A1{2},A1{3},A1{4});
        
        Energy_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 1:\n\nR2_1 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 0:\n\nR2_0 = R1_0*R2_1\n\nR2_0 = \n%s\n%s\n%s \n\n\nPosition: \n\nP1 = %s (m), P2 = %s (m) \n\nRotating Link 1 back to inertial frame:  Given by (P1*R1_0)\n\nP1_0 = %s (m)\n\nRotating Link 2 back to inertial frame:  Given by (P2*R2_0)\n\nP2_0 = %s (m)\n\n\nInertia:\n\nI1 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI1_0 = \n%s\n%s\n%s (kg.m^2) \n\nI2 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI2_0 = \n%s\n%s\n%s (kg.m^2) \n\nIx1 = %s (kg.m^2)\nIy1 = %s (kg.m^2)\nIz1 = %s (kg.m^2)\n\nIx2 = %s (kg.m^2)\nIy2 = %s (kg.m^2)\nIz2 = %s \n\n\nVelocity: \n\n\nLinear Velocity: Given by (dP/dt)\nV1 = %s (m/s) \n\nV2 = %s (m/s)\n\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\ndR2_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\n\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\n\nω1 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 2:  Given by ( S(ω2) = transpose(R2_0) * dR2_0 )\n\nS(ω2) = \n%s\n%s\n%s \n\nAngular Velocity of Link 2:  Given by ω2 = [S(ω2)(3,2) , S(ω2)(1,3) , S(ω2)(2,1)]\n\nω2 = %s (rad/s)\n\n\nEnergy: \n\nGiven by (E = T + V)\n\n\nKinetic Energy: Given by (T = 0.5*m*dP^2 + 0.5*I* ω^2)\n\nT1 = %s (J)\nT2 = %s (J)\n\nT = T1 + T2 \n\nT = %s (J)\n\n\nPotential Energy: Given by (V = m*g*h)\n\nV1 = %s (J)\nV2 = %s (J)\n\nV = V1+V2\n\nV = %s (J)",q1,dq2,ddq3,R1,R2,R3,rp1,rp2,rR1,rR2,In1,In1_0,In2,In2_0,Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,drR1,drR2,dR1,dR2,S1,v1,S2,v2,T1_,T2_,T,V1_,V2_,V);
        
        av_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 1:\n\nR2_1 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 0:\n\nR2_0 = R1_0*R2_1\n\nR2_0 = \n%s\n%s\n%s \n\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\ndR2_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\nω1 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 2:  Given by ( S(ω2) = transpose(R2_0) * dR2_0 ) \nS(ω2) = \n%s\n%s\n%s \n\nAngular Velocity of Link 2:  Given by ω2 = [S(ω2)(3,2) , S(ω2)(1,3) , S(ω2)(2,1)]\nω2 = %s (rad/s)",q1,dq2,ddq3,R1,R2,R3,dR1,dR2,S1,v1,S2,v2);       
        
    else
        EOM_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 1:\n\nR2_1 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 0:\n\nR2_0 = R1_0*R2_1\n\nR2_0 = \n%s\n%s\n%s \n\nRotate from fame 3 to frame 2:\n\nR3_2 = \n%s\n%s\n%s \n\nRotate from fame 3 to frame 0:\n\nR3_0 = R1_0*R2_1*R3_2\n\nR3_0 = \n%s\n%s\n%s \n\n\nPosition: \n\nP1 = %s (m), P2 = %s (m) and P3 = %s (m).\n\nRotating Link 1 back to inertial frame:  Given by (P1*R1_0)\n\nP1_0 = %s (m)\n\nRotating Link 2 back to inertial frame:  Given by (P2*R2_0)\n\nP2_0 = %s (m)\n\nRotating Link 3 back to inertial frame: Given by (P3*R3_0)\n\nP3_0 = %s (m)\n\nInertia:\n\nI1 = \n%s\n%s\n%s (kg.m^2) \n\nRotating back to inertial frame:\n\nI1_0 = \n%s\n%s\n%s (kg.m^2) \n\nI2 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI2_0 = \n%s\n%s\n%s (kg.m^2) \n\nI3 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI3_0 = \n%s\n%s\n%s (kg.m^2)\n\nIx1 = %s (kg.m^2)\nIy1 = %s (kg.m^2)\nIz1 = %s (kg.m^2)\n\nIx2 = %s (kg.m^2)\nIy2 = %s (kg.m^2)\nIz2 = %s (kg.m^2)\n\nIx3 = %s (kg.m^2)\nIy3 = %s (kg.m^2)\nIz3 = %s (kg.m^2)\n\n\nVelocity: \n\nLinear Velocity: Given by (dP/dt)\n\nV1 = %s (m/s) \n\nV2 = %s (m/s) \n\nV3 = %s (m/s).\n\n\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\ndR2_0 = \n%s\n%s\n%s \n\ndR3_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\n\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\n\nω1 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 2:  Given by ( S(ω2) = transpose(R2_0) * dR2_0 )\n\nS(ω2) = \n%s\n%s\n%s \n\nAngular Velocity of Link 2:  Given by ω2 = [S(ω2)(3,2) , S(ω2)(1,3) , S(ω2)(2,1)]\n\nω2 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 3: Given by ( S(ω3) = transpose(R3_0) * dR3_0)\n\nS(ω3) = \n%s\n%s\n%s \n\nAngular Velocity of Link 3:  Given by ω3 = [S(ω3)(3,2) , S(ω3)(1,3) , S(ω3)(2,1)]\n\nω3 = %s (rad/s)\n\n\nEnergy: \nGiven by (E = T + V)\n\n\nKinetic Energy: Given by (T = 0.5*m*dP^2 + 0.5*I* ω^2)\n\nT1 = %s (J)\nT2 = %s (J)\nT3 = %s (J)\n\nT = T1 + T2 + T3\n\nT = %s (J)\n\n\nPotential: Given by (V = m*g*h)\n\nV1 = %s (J)\nV2 = %s (J)\nV3 = %s (J)\n\nV = V1+V2+V3\n\nV = %s (J)\n\n\nNow to determine the Equations of motion use the equation:\n\n𝑄= M(𝑞)∙ddq+𝐶(𝑞,dq)+𝐺(𝑞)\n\n\nWhere Q represents the input torques:\n\nQ = %s (Nm)\n\n\nM is the mass matrix given by ( Mij(q) = ∂T/(∂ (dqi)*d(∂qj)) ):\n\nM = \n%s\n%s\n%s \n\n\nC contains the centrifugal and coriolis accelerations given by ( C(q,dq) = (dM/dt)*dq – ∂T/∂q ):\n\ndM = \n%s\n%s\n%s \n\nC = \n%s \n\n\nG is the matrix that contains the potential energy given by: ( G(q) = ∂V/∂q )\n\nG = \n%s \n\n\nNow the EOM:\n\n%s \n\nSolving for acceleration: \n\nddphi = %s (rad/s^2)\n\nddpsi = %s (rad/s^2)\n\nddth = %s (m/s^2)",q1,dq2,ddq3,R1,R2,R3,R4,R5,rp1,rp2,rp3,rR1,rR2,rR3,In1,In1_0,In2,In2_0,In3,In3_0,Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,Ix3,Iy3,Iz3,drR1,drR2,drR3,dR1,dR2,dR3,S1,v1,S2,v2,S3,v3,T1_,T2_,T3_,T,V1_,V2_,V3_,V,t,M1,dM1,C1,G1,e1,e_);
        
        Control_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 1:\n\nR2_1 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 0:\n\nR2_0 = R1_0*R2_1\n\nR2_0 = \n%s\n%s\n%s \n\nRotate from fame 3 to frame 2:\n\nR3_2 = \n%s\n%s\n%s \n\nRotate from fame 3 to frame 0:\n\nR3_0 = R1_0*R2_1*R3_2\n\nR3_0 = \n%s\n%s\n%s \n\n\nPosition: \n\nP1 = %s (m), P2 = %s (m) and P3 = %s (m).\n\nRotating Link 1 back to inertial frame:  Given by (P1*R1_0)\n\nP1_0 = %s (m)\n\nRotating Link 2 back to inertial frame:  Given by (P2*R2_0)\n\nP2_0 = %s (m)\n\nRotating Link 3 back to inertial frame: Given by (P3*R3_0)\n\nP3_0 = %s (m)\n\nInertia:\n\nI1 = \n%s\n%s\n%s (kg.m^2) \n\nRotating back to inertial frame:\n\nI1_0 = \n%s\n%s\n%s (kg.m^2) \n\nI2 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI2_0 = \n%s\n%s\n%s (kg.m^2) \n\nI3 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI3_0 = \n%s\n%s\n%s (kg.m^2)\n\nIx1 = %s (kg.m^2)\nIy1 = %s (kg.m^2)\nIz1 = %s (kg.m^2)\n\nIx2 = %s (kg.m^2)\nIy2 = %s (kg.m^2)\nIz2 = %s (kg.m^2)\n\nIx3 = %s (kg.m^2)\nIy3 = %s (kg.m^2)\nIz3 = %s (kg.m^2)\n\n\nVelocity: \n\nLinear Velocity: Given by (dP/dt)\n\nV1 = %s (m/s) \n\nV2 = %s (m/s) \n\nV3 = %s (m/s).\n\n\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\ndR2_0 = \n%s\n%s\n%s \n\ndR3_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\n\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\n\nω1 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 2:  Given by ( S(ω2) = transpose(R2_0) * dR2_0 )\n\nS(ω2) = \n%s\n%s\n%s \n\nAngular Velocity of Link 2:  Given by ω2 = [S(ω2)(3,2) , S(ω2)(1,3) , S(ω2)(2,1)]\n\nω2 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 3: Given by ( S(ω3) = transpose(R3_0) * dR3_0)\n\nS(ω3) = \n%s\n%s\n%s \n\nAngular Velocity of Link 3:  Given by ω3 = [S(ω3)(3,2) , S(ω3)(1,3) , S(ω3)(2,1)]\n\nω3 = %s (rad/s)\n\n\nEnergy: \nGiven by (E = T + V)\n\n\nKinetic Energy: Given by (T = 0.5*m*dP^2 + 0.5*I* ω^2)\n\nT1 = %s (J)\nT2 = %s (J)\nT3 = %s (J)\n\nT = T1 + T2 + T3\n\nT = %s (J)\n\n\nPotential: Given by (V = m*g*h)\n\nV1 = %s (J)\nV2 = %s (J)\nV3 = %s (J)\n\nV = V1+V2+V3\n\nV = %s (J)\n\n\nNow to determine the Equations of motion use the equation:\n\n𝑄= M(𝑞)∙ddq+𝐶(𝑞,dq)+𝐺(𝑞)\n\n\nWhere Q represents the input torques:\n\nQ = %s (Nm)\n\n\nM is the mass matrix given by ( Mij(q) = ∂T/(∂ (dqi)*d(∂qj)) ):\n\nM = \n%s\n%s\n%s \n\n\nC contains the centrifugal and coriolis accelerations given by ( C(q,dq) = (dM/dt)*dq – ∂T/∂q ):\n\ndM = \n%s\n%s\n%s \n\nC = \n%s \n\n\nG is the matrix that contains the potential energy given by: ( G(q) = ∂V/∂q )\n\nG = \n%s \n\n\nNow the EOM:\n\n%s \n\nFor the pseudo control input, substitute ddphi1 = v1, ddpsi2 = v2, ddth3 = v3 and u = [t1, t2, t3] (u will be a vector of inputs), then solve for u: \n\nu = %s\n\n%s\n\n%s (Nm)",q1,dq2,ddq3,R1,R2,R3,R4,R5,rp1,rp2,rp3,rR1,rR2,rR3,In1,In1_0,In2,In2_0,In3,In3_0,Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,Ix3,Iy3,Iz3,drR1,drR2,drR3,dR1,dR2,dR3,S1,v1,S2,v2,S3,v3,T1_,T2_,T3_,T,V1_,V2_,V3_,V,t,M1,dM1,C1,G1,e1,c1);
        
        M_choice_EOM = sprintf("A.\n\nddphi1 = %s\n\nddpsi2 = %s\n\nddth3 = %s\n\nB.\n\nddphi1 = %s\n\nddpsi2 = %s\n\nddth3 = %s\n\nC.\n\nddphi1 = %s\n\nddpsi2 = %s\n\nddth3 = %s\n\nD.\n\nddphi1 = %s\n\nddpsi2 = %s\n\nddth3 = %s",A1{1},A1{2},A1{3},A1{4});
        
        Energy_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 1:\n\nR2_1 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 0:\n\nR2_0 = R1_0*R2_1\n\nR2_0 = \n%s\n%s\n%s \n\nRotate from fame 3 to frame 2:\n\nR3_2 = \n%s\n%s\n%s \n\nRotate from fame 3 to frame 0:\n\nR3_0 = R1_0*R2_1*R3_2\n\nR3_0 = \n%s\n%s\n%s \n\n\nPosition: \n\nP1 = %s (m), P2 = %s (m) and P3 = %s (m).\n\nRotating Link 1 back to inertial frame:  Given by (P1*R1_0)\n\nP1_0 = %s (m)\n\nRotating Link 2 back to inertial frame:  Given by (P2*R2_0)\n\nP2_0 = %s (m)\n\nRotating Link 3 back to inertial frame: Given by (P3*R3_0)\n\nP3_0 = %s (m)\n\nInertia:\n\nI1 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI1_0 = \n%s\n%s\n%s (kg.m^2) \n\nI2 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI2_0 = \n%s\n%s\n%s (kg.m^2) \n\nI3 = \n%s\n%s\n%s (kg.m^2)  \n\nRotating back to inertial frame:\n\nI3_0 = \n%s\n%s\n%s (kg.m^2)\n\nIx1 = %s (kg.m^2)\nIy1 = %s (kg.m^2)\nIz1 = %s (kg.m^2)\n\nIx2 = %s (kg.m^2)\nIy2 = %s (kg.m^2)\nIz2 = %s (kg.m^2)\n\nIx3 = %s (kg.m^2)\nIy3 = %s (kg.m^2)\nIz3 = %s (kg.m^2)\n\n\nVelocity: \n\nLinear Velocity: Given by (dP/dt)\n\nV1 = %s (m/s) \n\nV2 = %s (m/s) \n\nV3 = %s (m/s).\n\n\n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\ndR2_0 = \n%s\n%s\n%s \n\ndR3_0 = \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\n\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\n\nω1 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 2:  Given by ( S(ω2) = transpose(R2_0) * dR2_0 )\n\nS(ω2) = \n%s\n%s\n%s \n\nAngular Velocity of Link 2:  Given by ω2 = [S(ω2)(3,2) , S(ω2)(1,3) , S(ω2)(2,1)]\n\nω2 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 3: Given by ( S(ω3) = transpose(R3_0) * dR3_0)\n\nS(ω3) = \n%s\n%s\n%s \n\nAngular Velocity of Link 3:  Given by ω3 = [S(ω3)(3,2) , S(ω3)(1,3) , S(ω3)(2,1)]\n\nω3 = %s (rad/s)\n\n\nEnergy: \nGiven by (E = T + V)\n\n\nKinetic Energy: Given by (T = 0.5*m*dP^2 + 0.5*I* ω^2)\n\nT1 = %s (J)\nT2 = %s (J)\nT3 = %s (J)\n\nT = T1 + T2 + T3\n\nT = %s (J)\n\n\nPotential: Given by (V = m*g*h)\n\nV1 = %s (J)\nV2 = %s (J)\nV3 = %s (J)\n\nV = V1+V2+V3\n\nV = %s (J)",q1,dq2,ddq3,R1,R2,R3,R4,R5,rp1,rp2,rp3,rR1,rR2,rR3,In1,In1_0,In2,In2_0,In3,In3_0,Ix1,Iy1,Iz1,Ix2,Iy2,Iz2,Ix3,Iy3,Iz3,drR1,drR2,drR3,dR1,dR2,dR3,S1,v1,S2,v2,S3,v3,T1_,T2_,T3_,T,V1_,V2_,V3_,V);
        
        av_sol = sprintf("Solution:\n\n\nGeneralised co-ordinates:\n\nq = %s (rad)\ndq = %s (rad/s)\nddq = %s (rad/s^2)\n\n\nRotation Matrices:\n\nRotating from fame 1 to frame 0:\n\nR1_0 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 1:\n\nR2_1 = \n%s\n%s\n%s \n\nRotating from fame 2 to frame 0:\n\nR2_0 = R1_0*R2_1\n\nR2_0 = \n%s\n%s\n%s \n\nRotate from fame 3 to frame 2:\n\nR3_2 = \n%s\n%s\n%s \n\nRotate from fame 3 to frame 0:\n\nR3_0 = R1_0*R2_1*R3_2\n\nR3_0 = \n%s\n%s\n%s \n\nAngular velocity:\n\nDerivative of rotation matrices:  Given by (dR/dt)\n\ndR1_0 = \n%s\n%s\n%s \n\ndR2_0 = \n%s\n%s\n%s \n\ndR3_0 =  \n%s\n%s\n%s \n\nDetermine the skew symmetric matrix for link 1:  Given by ( S(ω1) = transpose(R1_0) * dR1_0 )\n\nS(ω1) = \n%s\n%s\n%s \n\nAngular Velocity of Link 1:  Given by ω1 = [S(ω1)(3,2) , S(ω1)(1,3) , S(ω1)(2,1)]\n\nω1 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 2:  Given by ( S(ω2) = transpose(R2_0) * dR2_0 ) \n\nS(ω2) = \n%s\n%s\n%s \n\nAngular Velocity of Link 2:  Given by ω2 = [S(ω2)(3,2) , S(ω2)(1,3) , S(ω2)(2,1)]\n\nω2 = %s (rad/s)\n\nDetermine the skew symmetric matrix for link 3: Given by ( S(ω3) = transpose(R3_0) * dR3_0)\n\nS(ω3) = \n%s\n%s\n%s \n\nAngular Velocity of Link 3:  Given by ω3 = [S(ω3)(3,2) , S(ω3)(1,3) , S(ω3)(2,1)]\n\nω3 = %s (rad/s)",q1,dq2,ddq3,R1,R2,R3,R4,R5,dR1,dR2,dR3,S1,v1,S2,v2,S3,v3);
        
    end
    
    %+EOM_sol = latex_subs(EOM_sol);
    
    assignin('base','av_sol',av_sol);
    assignin('base','Energy_sol',Energy_sol);
    assignin('base','EOM_sol',EOM_sol);
    assignin('base','Control_sol',Control_sol);
    assignin('base','M_choice_EOM',M_choice_EOM);

    %----------------------------------------------------------------------
 
end
