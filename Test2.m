 %% Variables to be used
%--------------------------------------------------------------------------
syms g m1 m2'real'       
syms Ix1 Iy1 Iz1 Ix2 Iy2 Iz2'real' 
syms phi1 psi2'real'     
syms dphi1 dpsi2'real'   
syms ddphi1 ddpsi2 ddpsi2_ ddphi1_'real'
syms L1 L2 rc1 rc2 'real'         
syms t1 t2'real'
%----------------------------------------------------------------------

%% Generalised Co-ordinates
%--------------------------------------------------------------------------
q = [phi1; psi2];    
dq = [dphi1; dpsi2];
ddq = [ddphi1; ddpsi2];
%----------------------------------------------------------------------

%% Inertia
%--------------------------------------------------------------------------
% Inertias for rotation about each axis in body frame
% Treat links as cylinders
Ix1 = 1/12*m1*(3*rc1^2+L1^2);
Iy1 = Ix1;
Iz1 = (1/2)*m1*rc1^2;

Ix2 = 1/12*m1*(3*rc2^2+L2^2);
Iy2 = (1/2)*m1*rc2^2;
Iz2 = Ix2; 
%----------------------------------------------------------------------

%% KINEMATICS
%--------------------------------------------------------------------------
%  Position
%----------------------------------------------------------------------
% position of link 1 COM in frame 1
r1 = [0;0;L1/2];

% position of link 2 COM in frame 2
r2 = [0;L2/2;0];

% Rotation of link 1 to frame 0 (inertial/base frame)
% Rotation of body and inertia
Rf_1_0 = transpose(RX(q(1)));

% Rotated position of link 1
r1_0 = Rf_1_0 * r1;

% Rotation of link 2 to frame 0 (inertial/base frame)
% Rotation of body and inertia

% Rotation from from 1 to 2
Rf_1_2 = RX(q(2));

% Rotation from frame 2 to 0
Rf_2_0 = simplify(transpose(Rf_1_2 * transpose(Rf_1_0)),'IgnoreAnalyticConstraints',true);

% Rotated position of link 2
r2_0 = r1_0 + Rf_2_0*r2; 

%--------------------------------------------------------------------------

% Velocity (Linear and Angular)
%--------------------------------------------------------------------------
% Linear Velocity in frame 0
dr1 = jacobian(r1_0, q) * dq;
dr2 = simplify(jacobian(r2_0, q) * dq,'IgnoreAnalyticConstraints',true);

% Angular velocity using the skew symmetric matrix
% Derivative of rotation matrix
dRf1_0 = sym(zeros(3, length(Rf_1_0)));

% link 1 angular velocity
[omega1,S_omega1,dRf1_0] = ang_vel(Rf_1_0, dRf1_0, q, dq);

% Derivative of rotation matrix
dRf2_0 = sym(zeros(3, length(Rf_2_0)));

% link 2 angular velocity
[omega2,S_omega2,dRf2_0] = ang_vel(Rf_2_0, dRf2_0, q, dq);

%% ENERGY
%--------------------------------------------------------------------------
I11 = diag([Ix1,Iy1,Iz1]);
I1_0 = simplify(Rf_1_0*I11);

I22 = diag([Ix2,Iy2,Iz2]);
I2_0 = simplify(Rf_2_0*I22);

% Kinetic Energy
%----------------------------------------------------------------------     
% translational Energy
T1_t = 0.5*m1*(dr1)'*dr1; 

% rotational Energy
T1_r = 0.5*omega1'*I1_0*omega1; 

% Total Kinetic energy of link 1
T1 = T1_r + T1_t; 
Ttot = simplify(T1,'IgnoreAnalyticConstraints',true);

% translational Energy
T2_t = 0.5*m2*(dr2)'*dr2;

% rotational Energy
T2_r = 0.5*(omega2)'*I2_0*omega2;

% Kinetic energy of link 2
T2 = T2_r + T2_t;

% Total Kinetic energy
Ttot = simplify(Ttot + T2,'IgnoreAnalyticConstraints',true);

%----------------------------------------------------------------------

% Potential Energy
%----------------------------------------------------------------------
% Potential Energy of link 1
V1 = m1*[0 0 g]*r1_0; 

% Total potential energy
Vtot = V1;

% Potential Energy of link 2
V2 = m2*[0 0 g]*r2_0;

% Total potential energy
Vtot = Vtot + V2;

%% Dynamics
%--------------------------------------------------------------------------
%% Mass Matrix
%----------------------------------------------------------------------
M = hessian(Ttot,dq);
M = simplify(M,'Steps',3);
%----------------------------------------------------------------------

%% Derivative of Mass Matrix
%----------------------------------------------------------------------
dM = M_derivative(q,dq,M);
dM = simplify(dM,'Steps',3);
%----------------------------------------------------------------------

%% C Matrix
%----------------------------------------------------------------------
C = dM * dq - jacobian(Ttot, q)';
C = simplify(C,'Steps',3);
%----------------------------------------------------------------------

%% G Matrix 
%----------------------------------------------------------------------
G = simplify(transpose(jacobian(Vtot,q)),'Steps',3);
%----------------------------------------------------------------------

%% B input matrix 
%----------------------------------------------------------------------
B = [t1; t2];
%----------------------------------------------------------------------

%% EOM
%----------------------------------------------------------------------
eqn = simplify(M*ddq(1:2) + C + G - B,'Steps',5);

%% Substitution
%----------------------------------------------------------------------
m1_ = 3;
m2_ = 2;

L1_ = 0.3;
L2_ = 0.25;

g_ = 9.81;
r1_ = 0.3;
r2_ = 0.25;


eqn = subs(eqn,[rc1,rc2],[r1_,r2_]);
eqn = subs(eqn, [g, m1, m2, L1, L2], [g_, m1_, m2_, L1_, L2_]);
eqn = simplify(subs(eqn, g, g_),'Steps',5);
%----------------------------------------------------------------------

%% Solving for acceleration
%----------------------------------------------------------------------

sol = solve(eqn, ddq);    

ddphi1_eqn = simplify(sol.ddphi1,'Steps',5);
ddpsi2_eqn = simplify(sol.ddpsi2,'Steps',5);

%% Rotation matrices
%----------------------------------------------------------------------
% Rotation about X
function rotX = RX(theta) 
    rotX = [1 0 0;
            0 cos(theta) sin(theta);
            0 -sin(theta) cos(theta)];
end

% Rotation about Y
function rotY = RY(theta)
    rotY = [cos(theta) 0 -sin(theta);
            0 1 0;
            sin(theta) 0 cos(theta)];
end

% Rotation about Z
function rotZ = RZ(theta)
    rotZ = [cos(theta) sin(theta) 0;
            -sin(theta) cos(theta) 0;
            0 0 1];
end