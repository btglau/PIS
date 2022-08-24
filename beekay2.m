function bk = beekay2(l,u,U)
% calculate the output of b^k for angular momentum integrals in
% the uncoupled angular momentum representation, with real spherical
% harmonics

% this upgraded function uses a precomputed U

% Ref: Brink & Satchler, Angular momentum, Section 6.3
% and: 10.1016/S0166-1280(96)90531-X or J MOL STRUC-THEOCHEM 368 31-37 1996
    
    % move 0 to right
    l_save = l(2:3);
    [~,ix] = sort(abs(u),'descend');
    u = u(ix);
    l = l(ix);
    li = l+1;

    bk = 0;
    % calculate the U part
    switch nnz(u(2:3))
        case 2
            if abs(u(1)) == abs(u(2) - u(3))
                bk = 2*real(conj(U{li(1)}(u(1)+li(1),u(2)-u(3)+li(1)))*...
                                 U{li(2)}(u(2)+li(2),u(2)+li(2))*...
                                 U{li(3)}(u(3)+li(3),-u(3)+li(3)));
                % for calculating the Gaunt coeff
                u(1) = u(2) - u(3);
                u(3) = -u(3);
            else % abs(u(1)) == abs(u(2) + u(3))
                bk = 2*real(conj(U{li(1)}(u(1)+li(1),u(2)+u(3)+li(1)))*...
                                 U{li(2)}(u(2)+li(2),u(2)+li(2))*...
                                 U{li(3)}(u(3)+li(3),u(3)+li(3)));
                u(1) = u(2) + u(3);
            end
        case 1
            bk = 2*real(conj(U{li(1)}(u(1)+li(1),u(2)+li(1)))*...
                             U{li(2)}(u(2)+li(2),u(2)+li(2)));
            u(1) = u(2);
        case 0
            if u(1) == 0
                bk = 1;
            end
    end 
    
    if bk ~= 0
        % calculate Gaunt coefficient
        bk = bk*(-1)^u(1);
        u(1) = -u(1);
        bk = bk*sqrt((2*l_save(1)+1)*(2*l_save(2)+1))*Wigner3j(l,[0;0;0])*Wigner3j(l,u);
    end
end