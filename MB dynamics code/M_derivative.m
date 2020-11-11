function [dM] = M_derivative(q,dq,M)
    dim = size(M);
    dM = sym(ones(dim(1),dim(2)));
    
    for i = 1:length(M)
        for j = 1:length(M)
            dM(i,j) = jacobian(M(i,j),q)*dq;
        end
    end
    
    dM = simplify(dM,'IgnoreAnalyticConstraints',true);
end