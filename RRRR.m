function pppp = RRRR(k,l,knl,N1234,r1,r2)
% Returns a vectorized output of the product of four spherical Bessel
% functions.
% integral2: "The function fun must accept two arrays of the same size and 
% return an array of corresponding values. It must perform element-wise
% operations."
% integral2 can handle singularities at the boundary (r = 0) - but slows
% down the code 
    
    Ra = sqrt(pi./(2*knl(1)*r1)) .* besselj(l(1) + 1/2,knl(1)*r1);
    Rb = sqrt(pi./(2*knl(2)*r1)) .* besselj(l(2) + 1/2,knl(2)*r1);
    Rc = sqrt(pi./(2*knl(3)*r2)) .* besselj(l(3) + 1/2,knl(3)*r2);
    Rd = sqrt(pi./(2*knl(4)*r2)) .* besselj(l(4) + 1/2,knl(4)*r2);

    % special value at 0, r1
    ind = r1 == 0;
    if any(ind)
        if l(1) == 0
            Ra(ind) = 1;
        else
            Ra(ind) = 0;
        end
        if l(2) == 0
            Rb(ind) = 1;
        else
            Rb(ind) = 0;
        end
    end
    % repeat for r2
    ind = r2 == 0;
    if any(ind)
        if l(3) == 0
            Rc(ind) = 1;
        else
            Rc(ind) = 0;
        end
        if l(4) == 0
            Rd(ind) = 1;
        else
            Rd(ind) = 0;
        end
    end

    pppp = N1234 * r1.^(1-k).*r2.^(k+2) .* Ra.*Rb.*Rc.*Rd;
end