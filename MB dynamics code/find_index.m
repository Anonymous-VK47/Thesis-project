function [index] = find_index(A,e_)
    index = 4;
    for i=1:4
        if isequal(A{i},e_)
            index = i;
        end
    end
end