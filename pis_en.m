function [En,ABknl] = pis_en(V,n,l)
%PIS_EN returns energies for a given set of n and l

% output arguments:
% En: energies for each pair of n and l
% ABknl: the zero of j for n and l, 
%        or a vector of expansion coefficients (for multistep case)

    if isscalar(V)
        % particle in an infinite step
        % pretabulated roots - index is l,n
        % Roots of j_l are roots of J_{l+1/2}
        roots = besselzero((0:max(l)).'+0.5,max(n));
        ind = sub2ind(size(roots),l+1,n);
        ABknl(:,1) = roots(ind);
        % energy is just k^2; stay unitless, with energy scale 1/(2*m*r^2)
        En = ABknl.^2;
        % also return normalization (unit sphere, R = 1)
        ABknl(:,2) = sqrt(2) ./ abs(sqrt(pi./(2*ABknl(:,1))) .* besselj(l+1+1/2,ABknl(:,1)));
    else
        % functionality for multistep not integrated yet
        
        % vector of length n / l
        En = 0;
        
        % matrix of row length n / l, column length of 2*steps - 2
        ABknl = 0;
    end
end