function U = Ulmu(u,m)
% return elements of the transformation matrix U for complex -> real
% spherical harmonics: U_{lm}^u, where u is row and m is col index.
% vectorized.
    
    % commented out line returns negative real spherical harmonics for m
    % odd
    % U = (m==0).*(u==0) 
    %   + (-u>0).*1i.*(-1).^m.*(m==u)
    %   + ((u>0).*(m==u)
    %   - 1i*(-u>0).*(m==-u) 
    %   + (u>0).*(-1).^m.*(m==-u))./sqrt(2);
    
    % fixed line returns positive real spherical harmonics
    U = (m==0).*(u==0) ...
        + ((u>0).*(-1).^m.*(m==u) ...
        + 1i*(u<0).*(m==u) ...
        - 1i*(u<0).*(-1).^m.*(m==-u) ...
        + (u>0).*(m==-u))./sqrt(2);
end