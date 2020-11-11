function [wrong_M1,wrong_M2] = check_M(a,b)
    global number_of_links;
    
    wrong_M1 = a;
    wrong_M2 = b;
    if (number_of_links==1) && isequal(a(1,1),sym(0)) 
        wrong_M1(1,1) = 1;
    elseif (number_of_links==2) && isequal(a(2,2),sym(0)) 
        wrong_M1(2,2) = 1;
    elseif (number_of_links==3) && isequal(a(3,3),sym(0)) 
        wrong_M1(3,3) = 1;
    elseif (number_of_links==1) && isequal(b(1,1),sym(0)) 
        wrong_M2(1,1) = 1;
    elseif (number_of_links==2) && isequal(b(2,2),sym(0))
        wrong_M2(2,2) = 1;
    elseif (number_of_links==3) && isequal(b(3,3),sym(0))
        wrong_M2(3,3) = 1; 
    end
end