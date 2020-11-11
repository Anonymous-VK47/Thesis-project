function [Rf_1_0,Rf_2_0,Rf_3_0,velocity1,velocity2,position,R1,R2,R3,R4,R5,rp1,rp2,rp3,rR1,rR2,rR3,drR1,drR2,drR3,dR1,dR2,dR3,S1,v1,S2,v2,S3,v3,wrong_vel,l_n] = Kinematics(variables, q, dq, Rot1,Rot2,Rot3,RotX,RotY,RotZ)
    % This function deals with the kinematics of the system
    % The positions are established and rotated to their respective
    % reference frames
    % The velocities (linear and angular) are then calculated
    
    %% default values & variables
    %----------------------------------------------------------------------
    global number_of_links;
    global l_n;

    rng shuffle
    l_n = 1;
    assignin('base','l_n',l_n);
    r2_0 = sym(zeros(3,1));
    r3_0 = sym(zeros(3,1));
    R2 = string(zeros(3,3));
    R3 = string(zeros(3,3));
    R4 = string(zeros(3,3));
    R5 = string(zeros(3,3));
    dR2 = string(zeros(3,3));
    dR3 = string(zeros(3,3));
    S2 = string(zeros(3,3));
    v2 = string(zeros(3,1));
    S3 = string(zeros(3,3));
    v3 = string(zeros(3,1));
    Rf_2_0 = zeros(3,3);
    Rf_3_0 = zeros(3,3);

    %----------------------------------------------------------------------
    %% KINEMATICS

    %  Position
    %----------------------------------------------------------------------

    % position of link 1 COM in frame 1
    r1 = [0;0;variables(1)/2];
    rp1 = "[" + replace(join(string(r1))," ",", ") + "]";
    
    % position of link 2 COM in frame 2
    r2 = [0;variables(2)/2;0];
    rp2 = "[" + replace(join(string(r2))," ",", ") + "]";
    
    % position of link 3 COM in frame 3
    r3 = [0;0;-1*variables(3)/2];
    rp3 = create_A_string1(r3);
    
    % Rotation of link 1 to frame 0 (inertial/base frame)
    % Rotation of body and inertia
    Rf_1_0 = transpose(Rot1(q(1)));
    
    % If student uses wrong rotation matrix
    w_R1 = transpose(Rot2(q(1)));
    
    % If student forgets to take the transpose
    w_R2 = transpose(Rf_1_0);
    
    % Add random co-efficient
    w_R3 = Rf_1_0/2;
    
    % Rotation matrix as text for GUI
    R1 = "["+ join(string(Rf_1_0))+ "]";

    % Rotated position of link 1
    r1_0 = Rf_1_0 * r1;
    
    % File export for testing
    pos1 = fopen('C:\Multibody Dynamics Generator\pos1.txt','wt');
    fprintf(pos1,'%s',"["+ join(string(r1_0))+ "]");
    fclose(pos1);
    
    % If student forgets to rotate back to inertial frame
    wrong_pos = r1;

    % Rotation of link 2 to frame 0 (inertial/base frame)
    % Rotation of body and inertia
    if number_of_links >= 2
        l_n = randi(2);
        assignin('base','l_n',l_n);
        
        % Rotation from from 1 to 2
        Rf_1_2 = Rot2(q(2));
        
        % Rotation matrix as text for GUI
        R2 = "["+ join(string(transpose(Rf_1_2)))+ "]";
        
        % Rotation from frame 2 to 0
        Rf_2_0 = simplify(transpose(Rf_1_2 * transpose(Rf_1_0)),'IgnoreAnalyticConstraints',true);
        
        % Rotation matrix as text for GUI
        R3= "["+ join(string(Rf_2_0))+ "]";

        % Rotated position of link 2
        r2_0 = r1_0 + Rf_2_0*r2; 
        
        % File export for testing
        pos2 = fopen('C:\Multibody Dynamics Generator\pos2.txt','wt');
        fprintf(pos2,'%s',"["+ join(string(r2_0))+ "]");
        fclose(pos2);
        
        % If student forgets to rotate back to inertial frame
        wrong_pos = r1_0 + transpose(Rf_1_2)*r2;
        
        % Rotation matrices as text for GUI
        w_R4 = transpose(Rf_1_2 * transpose(w_R1));
        w_R5 = transpose(Rf_1_2 * transpose(w_R2));
        w_R6 = transpose(Rf_1_2 * transpose(w_R3));
    end
    
    % Rotation of link 3 to frame 0 (inertial/base frame)
    % Rotation of body and inertia
    if number_of_links == 3
        l_n = randi(3);
        assignin('base','l_n',l_n);
        
        % Rotation from frame 2 to 3
        Rf_2_3 = Rot3(q(3));
        
        % Rotation matrix as text for GUI
        R4 = "["+ join(string(transpose(Rf_2_3)))+ "]";
        
        % Rotation from frame 3 to 0
        Rf_3_0 = simplify(transpose(Rf_2_3 * Rf_1_2 * transpose(Rf_1_0)),'IgnoreAnalyticConstraints',true);
        
        % Rotation matrix as text for GUI
        R5 = create_A_string1(Rf_3_0);

        % Rotated postion of link 3
        r3_0 = r2_0 + Rf_3_0*r3;
        
        % File export for testing
        pos3 = fopen('C:\Multibody Dynamics Generator\pos3.txt','wt');
        fprintf(pos3,'%s',"["+ join(string(r3_0))+ "]");
        fclose(pos3);
        
        % If student forgets to rotate back to inertial frame
        wrong_pos = r2_0 + transpose(Rf_2_3)*r3;
        
        % Rotation matrices as text for GUI
        w_R7 = transpose(Rf_2_3 * w_R4);
        w_R8 = transpose(Rf_2_3 *w_R5);
        w_R9 = transpose(Rf_2_3 *w_R6);
    end
    
    % Output Vector containing rotated positions of each link
    position = [r1_0, r2_0, simplify(r3_0,'IgnoreAnalyticConstraints',true)];
    
    % Rotated positions of each link as text for GUI
    rR1 = create_A_string1(r1_0);
    rR2 = create_A_string1(r2_0);
    rR3 = create_A_string1(r3_0);

    %--------------------------------------------------------------------------

    % Velocity (Linear and Angular)
    %--------------------------------------------------------------------------

    % Linear
    
    % Linear Velocity in frame 0
    dr1 = jacobian(r1_0, q) * dq;
    dr2 = simplify(jacobian(r2_0, q) * dq,'IgnoreAnalyticConstraints',true);
    dr3 = simplify(jacobian(r3_0, q) * dq,'IgnoreAnalyticConstraints',true);
    
    % File export for testing
    vel1 = fopen('C:\Multibody Dynamics Generator\vel1.txt','wt');
    fprintf(vel1,'%s',"["+ join(string(dr1))+ "]");
    fclose(vel1);
    vel2 = fopen('C:\Multibody Dynamics Generator\vel2.txt','wt');
    fprintf(vel2,'%s',"["+ join(string(dr1))+ "]");
    fclose(vel2);
    vel3 = fopen('C:\Multibody Dynamics Generator\vel3.txt','wt');
    fprintf(vel3,'%s',"["+ join(string(dr3))+ "]");
    fclose(vel3);
    
    % Wrong linear velocity
    wrong_vel = simplify(jacobian(wrong_pos, q) * dq,'IgnoreAnalyticConstraints',true);
    
    % Output vector containing linear velocities of each link
    velocity1 = [dr1, dr2, dr3];
    
    % Wong Velocities of each link as text for GUI
    drR1 = create_A_string1(dr1);
    drR2 = create_A_string1(dr2);
    drR3 = create_A_string1(dr3);

    % Angular velocity using the skew symmetric matrix
    
    % Derivative of rotation matrix
    dRf1_0 = sym(zeros(3, length(Rf_1_0)));
    
    % link 1 angular velocity
    [omega1,S_omega1,dRf1_0] = ang_vel(Rf_1_0, dRf1_0, q, dq);
    
    % File export for testing
    ome1 = fopen('C:\Multibody Dynamics Generator\ome1.txt','wt');
    fprintf(ome1,'%s',"["+ join(string(omega1))+ "]");
    fclose(ome1);
    
    % Wrong anular velocities of link 1
    [w_o1,~,~] = ang_vel(w_R1, dRf1_0, q, dq);
    [w_o2,~,~] = ang_vel(w_R2, dRf1_0, q, dq);
    [w_o3,~,~] = ang_vel(w_R3, dRf1_0, q, dq);
    
    % Randomize angular velocities for multiple choice
    A2 = {w_o1,w_o2,w_o3,omega1};
    A2 = A2(randperm(numel(A2)));
    
    % Find position of correct angular velocity
    av_index1 = find_index(A2,omega1);
    assignin('base','av_index1',av_index1);
    
    % Gui text (What user sees)
    M_choice_av1 = sprintf("Find the angular velocity of link 1:\n\nA.\nω = [%s,%s,%s]\n\nB.\nω = [%s,%s,%s]\n\nC.\nω = [%s,%s,%s]\n\nD.\nω = [%s,%s,%s]",A2{1},A2{2},A2{3},A2{4});
    assignin('base','M_choice_av1',M_choice_av1);
    
    % Derivative of rotation marix as text for GUI
    dR1 = "["+ join(replace(string(dRf1_0),[" - ","- "," + ",],["-","-","+"]))+ "]";
    %dR1 = replace(string(dR1),["-","+"],[" - "," + ",]);
    
    % Skew symmetric matrix as text for GUI
    S1 = "["+ join(string(S_omega1))+ "]";
    
    if number_of_links >= 2
        % Derivative of rotation matrix
        dRf2_0 = sym(zeros(3, length(Rf_2_0)));
        
        % link 2 angular velocity
        [omega2,S_omega2,dRf2_0] = ang_vel(Rf_2_0, dRf2_0, q, dq);
        
        % File export for testing
        ome2 = fopen('C:\Multibody Dynamics Generator\ome2.txt','wt');
        fprintf(ome2,'%s',"["+ join(string(omega2))+ "]");
        fclose(ome2);
        
        % Skew symmetric matrix and derivative of rotation matrix as text for GUI
        S2 = "["+ join(string(S_omega2))+ "]";
        dR2 = "["+ join(replace(string(dRf2_0),[" - ","- "," + ",],["-","-","+"]))+ "]";
        %dR2 = replace(string(dR2),["-","+"],[" - "," + ",]);
        
        % Wrong anular velocities of link 2
        [w_o4,~,~] = ang_vel(w_R4, dRf1_0, q, dq);
        [w_o5,~,~] = ang_vel(w_R5, dRf1_0, q, dq);
        [w_o6,~,~] = ang_vel(w_R6, dRf1_0, q, dq);
        
        % Randomize angular velocities for multiple choice
        A3 = {w_o4,w_o5,w_o6,omega2};
        A3 = A3(randperm(numel(A3)));
        
        % Find position of correct angular velocity
        av_index2 = find_index(A3,omega2);
        assignin('base','av_index2',av_index2);
        
        % Gui text (What user sees)
        M_choice_av2 = sprintf("Find the angular velocity of link %s:\n\nA.\nω = [%s,%s,%s]\n\nB.\nω = [%s,%s,%s]\n\nC.\nω = [%s,%s,%s]\n\nD.\nω = [%s,%s,%s]",l_n,A3{1},A3{2},A3{3},A3{4});
        assignin('base','M_choice_av2',M_choice_av2);
        
    end
    
    if number_of_links == 3
        % Derivative of rotation matrix
        dRf3_0 = sym(zeros(3, length(Rf_3_0)));
        
        % link 3 angular velocity
        [omega3,S_omega3,dRf3_0] = ang_vel(Rf_3_0, dRf3_0, q, dq);
        
        % File export for testing
        ome3 = fopen('C:\Multibody Dynamics Generator\ome3.txt','wt');
        fprintf(ome3,'%s',"["+ join(string(omega3))+ "]");
        fclose(ome3);
        
        % Skew symmetric matrix and derivative of rotation matrix as text for GUI
        dR3 = "["+ join(replace(string(dRf3_0),[" - ","- "," + ",],["-","-","+"]))+ "]";
        %dR3 = replace(string(dR3),["-","+"],[" - "," + ",]);
        S3 = "["+ join(string(S_omega3))+ "]";
        
        % Wrong anular velocities of link 3
        [w_o7,~,~] = ang_vel(w_R7, dRf1_0, q, dq);
        [w_o8,~,~] = ang_vel(w_R8, dRf1_0, q, dq);
        [w_o9,~,~] = ang_vel(w_R9, dRf1_0, q, dq);
        
        % Randomize angular velocities for multiple choice
        A4 = {w_o7,w_o8,w_o9,omega3};
        A4 = A4(randperm(numel(A4)));
        
        % Find position of correct angular velocity
        av_index3 = find_index(A4,omega3);
        assignin('base','av_index3',av_index3);
        
        % Gui text (What user sees)
        M_choice_av3 = sprintf("Find the angular velocity of link %s:\n\nA.\nω = [%s,%s,%s]\n\nB.\nω = [%s,%s,%s]\n\nC.\nω = [%s,%s,%s]\n\nD.\nω = [%s,%s,%s]",l_n,A4{1},A4{2},A4{3},A4{4});
        assignin('base','M_choice_av3',M_choice_av3);
    end
    
    % Output vectors with angular velocity and text versions for GUI
    if number_of_links == 1
        velocity2 = omega1;
        v1 = create_A_string1(omega1);
    elseif number_of_links == 2    
        velocity2 = [omega1, omega2];
        v1 = create_A_string1(omega1);
        v2 = create_A_string1(omega2);
        
    else    
        velocity2 = [omega1, omega2, omega3];
        v1 = create_A_string1(omega1);
        v2 = create_A_string1(omega2);
        v3 = create_A_string1(omega3);
    end

    %--------------------------------------------------------------------------
end