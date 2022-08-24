function [n,l,m,m2,En,ABknl] = pis_AO(lmax)
%PIS_AO Return a sequence of n, l, m, and matlab-indexed m2 for later use,
%based upon a given lmax

% the energy cutoff is defined as E_cut = 1/2 k_lmax,1^2, where k_lmax,1 is
% the first zero k of lmax. See J Chem Phys, 118, 24, 2003

    if isscalar(lmax)
        sphjz = besselzero(((0:lmax)+0.5)',20); % row is l, col is n
        [l,n] = find(sphjz.^2<=sphjz(end,1).^2);
        l = l-1; % convert matlax index to quantum numbers

        % make a nmax and lmax matrix to feed to existing pis AO functions
        nl = zeros(2,lmax+1);
        nl(2,:) = 0:lmax;
        for a = 0:lmax
            nl(1,a+1) = numel(n(l==a));
        end
    else
        % vector lmax = number of n's for each l
        nl = zeros(2,length(lmax));
        nl(1,:) = lmax;
        nl(2,:) = 0:length(lmax)-1;
    end

    [n,l,m,m2] = pis_nlm(nl(1,:),nl(2,:));
    [En,ABknl] = pis_en(0,n,l);
end