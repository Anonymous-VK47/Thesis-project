function [Rot1,Rot2,Rot3,out] = Rot_random(RotX, RotY, RotZ)
    % Randomize rotations
    global number_of_links;

    % Matrix with function handles of rotation matrices for 2 or more links
    A = {RotX, RotY, RotZ};
    
    % Matrix with function handles of rotation matrices for 1 link
    % (To avoid single link rotating about its own axis, appearing
    % stationary)
    A1 = {RotX, RotY};
    
    % Default values
    Rot2 = A{2};
    Rot3 = A{3};

    %% Randomize matrices
    %----------------------------------------------------------------------
    
    % For single link use A1
    if number_of_links == 1
        Rot_list = A1(randperm(numel(A1)));
        Rot1 = Rot_list{1};
        
    % For 2 or more links use A
    else
        Rot_list = A(randperm(numel(A)));
        Rot1 = Rot_list{1};
        Rot2 = Rot_list{2};
        Rot3 = Rot_list{3};
    end
    %----------------------------------------------------------------------
    
    %% Order of rotations
    %----------------------------------------------------------------------
    % Create a vector with order or random rotations
    if isequal(Rot1,A{1})
        x = "X";
    elseif isequal(Rot1,A{2})
        x = "Y";
    else
        x = "Z";
    end
    
    if isequal(Rot2,A{1})
        y = "X";
    elseif isequal(Rot2,A{2})
        y = "Y";
    else
        y = "Z";
    end

    if isequal(Rot3,A{1})
        z = "X";
    elseif isequal(Rot3,A{2})
        z = "Y";
    else
        z = "Z";
    end
    %----------------------------------------------------------------------
    
    % Rotation order vector
    out = [x,y,z];

end