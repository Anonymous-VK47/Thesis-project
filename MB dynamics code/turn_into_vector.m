function [pos_link1_x,pos_link1_y,pos_link1_z] = turn_into_vector(p1)    

    pos_link1_x = vectorize(p1(1));
    pos_link1_y = vectorize(p1(2));
    pos_link1_z = vectorize(p1(3));
    
end
    