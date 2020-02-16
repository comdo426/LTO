function [f g h] = final_mass_fcn(X, s)
    f = -X(21*s + 7);
    
    if nargout > 1
    
    g = zeros(length(X), 1);
    g(21*s+7) = -1;
    
    if nargout > 2
    
    h = zeros(length(X), length(X));
    end
    end
end