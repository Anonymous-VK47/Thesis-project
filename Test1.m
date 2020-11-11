 %% Variables to be used
%--------------------------------------------------------------------------
syms g m1 'real'       
syms Ix1 Iy1 Iz1 'real' 
syms phi1 'real'     
syms dphi1 'real'   
syms ddphi1 ddphi1_'real'
syms L1 rc1 'real'         
syms t1 'real'
%----------------------------------------------------------------------

%% Generalised Co-ordinates
%--------------------------------------------------------------------------
q = [phi1];    
dq = [dphi1];
ddq = [ddphi1];
%----------------------------------------------------------------------

%% Inertia
%--------------------------------------------------------------------------
% Inertias for rotation about each axis in body frame
% Treat links as cylinders
Ix1 = 1/12*m1*(3*rc1^2+L1^2);
Iy1 = Ix1;
Iz1 = (1/2)*m1*rc1^2;

%----------------------------------------------------------------------

%% KINEMATICS
%--------------------------------------------------------------------------
%  Position
%----------------------------------------------------------------------
% position of link 1 COM in frame 1
r1 = [0;0;L1/2];

% Rotation of link 1 to frame 0 (inertial/base frame)
% Rotation of body and inertia
Rf_1_0 = transpose(RX(phi1));

% Rotated position of link 1
r1_0 = Rf_1_0 * r1;

%--------------------------------------------------------------------------

% Velocity (Linear and Angular)
%--------------------------------------------------------------------------
% Linear Velocity in frame 0
dr1 = jacobian(r1_0, q) * dq;

% Angular velocity using the skew symmetric matrix
% Derivative of rotation matrix
dRf1_0 = sym(zeros(3, length(Rf_1_0)));

% link 1 angular velocity
[omega1,S_omega1,dRf1_0] = ang_vel(Rf_1_0, dRf1_0, q, dq);

%% ENERGY
%--------------------------------------------------------------------------
I11 = diag([Ix1,Iy1,Iz1]);
I1_0 = simplify(Rf_1_0*I11);

% Kinetic Energy
%----------------------------------------------------------------------     
% translational Energy
T1_t = 0.5*m1*(dr1)'*dr1; 

% rotational Energy
T1_r = 0.5*omega1'*I1_0*omega1; 

% Total Kinetic energy of link 1
T1 = T1_r + T1_t; 
Ttot = simplify(T1,'IgnoreAnalyticConstraints',true);
%----------------------------------------------------------------------

% Potential Energy
%----------------------------------------------------------------------
% Potential Energy of link 1
V1 = m1*[0 0 g]*r1_0; 

% Total potential energy
Vtot = V1;
%----------------------------------------------------------------------

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
B = [t1];
%----------------------------------------------------------------------

%% EOM
%----------------------------------------------------------------------
eqn = simplify(M*ddq(1) + C + G - B,'Steps',5);

%% Substitution
%----------------------------------------------------------------------
m1_ = 3;

L1_ = 0.3;

g_ = 9.81;
r1_ = 0.3;

eqn = subs(eqn,[rc1],[r1_]);
eqn = subs(eqn, [g, m1, L1], [g_, m1_, L1_]);
eqn = simplify(subs(eqn, g, g_),'Steps',5);
%----------------------------------------------------------------------

%% Solving for acceleration
%----------------------------------------------------------------------

sol = solve(eqn, ddq);    

ddphi1_eqn = simplify(sol,'Steps',5);

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