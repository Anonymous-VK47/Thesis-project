function [RotX,RotY,RotZ] = Rotation

% Return function handles to each rotation function
    RotX = @RX;
    RotY = @RY;
    RotZ = @RZ;
    
end

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