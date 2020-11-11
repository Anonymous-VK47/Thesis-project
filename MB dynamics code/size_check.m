function [pos_link1_x, pos_link1_y, pos_link1_z] = size_check(p1,p2,p3,temppsi)  

    if size(p1) == size(1)
        pos_link1_x = p1.*ones(length(temppsi),1);
    else
        pos_link1_x = p1;
    end
    if size(p2) == size(1)
        pos_link1_y = p2.*ones(length(temppsi),1);
    else
        pos_link1_y = p2;
    end
    if size(p3) == size(1)
        pos_link1_z = p3.*ones(length(temppsi),1);
    else
        pos_link1_z = p3;
    end
    
end