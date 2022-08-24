function [n,l,m,m2] = pis_nlm(nmax,lmax)
%PIS_NLM Return a sequence of n, l, m, and matlab-indexed m2 for later use,
%based upon input to the program.
    
    % allowed values of m
    if ~isscalar(nmax)
        if isscalar(lmax)
            assert(length(nmax) == lmax+1,'# of ef per not equal to # of l''s!')
        else
            assert(length(nmax)<=length(lmax),'Must specify # of ef for each l!')
        end
    end
    if isscalar(lmax)
        ind = (0:lmax)*2 + 1;
        ind2 = 0:lmax;
    else
        ind = lmax*2 + 1;
        ind2 = lmax;
    end
    
    % total number of ef's
    ind = ind.*nmax;

    % preallocate
    n = zeros(sum(ind),1);
    l = repelem(ind2,ind).';
    m = zeros(size(n));
    % m2 is a matlab compatible index for the m's
    m2 = zeros(size(n));

    % add index 1
    ind = cumsum([1 ind]);

    nl_max = nmax;
    for a = 1:length(ind)-1
        if ~isscalar(nmax)
            nl_max = nmax(a);
        end
        n(ind(a):ind(a+1)-1) = repmat((1:nl_max).',2*ind2(a) + 1,1);
        m(ind(a):ind(a+1)-1) = repelem(-ind2(a):ind2(a),nl_max);
        m2(ind(a):ind(a+1)-1) = m(ind(a):ind(a+1)-1) + ind2(a) + 1;
    end
end

