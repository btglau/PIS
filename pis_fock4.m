function pis_fock4(varargin)
% a function that will calculate and save integrals for pyscf
% computes the entire eri matrix in a SPARSE representation,
% using real spherical harmonics to achieve 8-fold symmetry
% version 4: precomputes radial integrals by using magnetic angular momentum
% degeneracy
% optimizations: replace 'intersect' with ismembc

% see: 10.1016/S0166-1280(96)90531-X
% or: J MOL STRUC-THEOCHEM 368 31-37 1996

    % for compiled function param, need to convert from string to numbers
    if isdeployed
        for a = 1:nargin
            temp = [];
            if ischar(varargin{a})
                % str2num only converts if input is valid matlab syntax for a
                % number or matrix, leaving the name of name-value pairs as
                % strings
                temp = str2num(varargin{a});
            end
            if ~isempty(temp)
                varargin{a} = temp;
            end
        end
    end
    p = inputParser;
    p.FunctionName = 'input parser for pis_fock';
    p.CaseSensitive = true;
    
        % basis set parameters
        addParameter(p,'V',0);
            % potential steps
        addParameter(p,'r',1);
            % radius steps
        addParameter(p,'me',1);
            % effective mass steps
        addParameter(p,'er',1);
            % dielectric constant (steps not implemented yet)
        addParameter(p,'lmax',2,@(x) all(x>=0));
            % which angular momentum quantum numbers to calculate
        addParameter(p,'noint',false)
            % if true, don't do the integrals (for checking outputs)
        addParameter(p,'eris_per_piece',2E8);
            % NOTE: large numbers helps because it encompasses more of the
            % magnetic angular momentum degeneracy!!
            % todo: adaptive refinement of eris_per_piece, by sampling
            % random values of next and increasing epp until no change
            %{
            [n,l,u,~,En,ABknl] = pis_AO([10,9,9,8,8,7,7,7]);
            next = uint64([1 1 1 1 1 1]);
            N = uint64(length(n));
            
            [csec,~] = pis_8f_pw2_mex(N,uint64(1E7),next);C = unique([n(csec) l(csec)],'rows');size(C)
            ans =
                  362024           8
            
            [csec,~] = pis_8f_pw2_mex(N,uint64(5E7),next);C = unique([n(csec) l(csec)],'rows');size(C)
            ans =
                  472792           8
            
            [csec,~] = pis_8f_pw2_mex(N,uint64(7.5E7),next);C = unique([n(csec) l(csec)],'rows');size(C)
            ans =
                  773358           8 <- unlucky, hitting a patch of 'new' angular momentum 

            [csec,~] = pis_8f_pw2_mex(N,uint64(1E8),next);C = unique([n(csec) l(csec)],'rows');size(C)
            ans =
                  774786           8
            
            tic;[csec,~] = pis_8f_pw2_mex(N,uint64(2E8),next);toc;C = unique([n(csec) l(csec)],'rows');size(C)
            ans =

                  774786           8
            
            tic;[csec,~] = pis_8f_pw2_mex(N,uint64(3E8),next);toc;C = unique([n(csec) l(csec)],'rows');size(C)  
            ans =

                 1275344           8

            %}
        addParameter(p,'resume','')
            % string to resume a previous calculation.
            
    % parse the params and insert them into the param structure
    parse(p,varargin{:});
    args = p.Results;
    epp = args.eris_per_piece;
    
    % fresh start, set up file name
    job_title = getenv('SLURM_JOB_ID');
    tor = datestr(now,'dd-mm-yy');
    
    if ~isempty(args.resume)
        % load in previous file
        if ~isempty(job_title)
            % if running on cluster (assuming that the matlab and basis set
            % folder are sitting at the same level)
            file_name = fullfile('..','basissets_8fold',args.resume);
        else
            file_name = args.resume;
        end
        %load(file_name);
        m = matfile(file_name,'Writable',true);
        
        if ~isempty(who(m,'Hcore'))
            disp(m.args)
            disp('This basis set has been calculated to completion!')
            return
        end
        args = m.args;        
        args.resume = 'UNFINISHED';
        args.eris_per_piece = epp;
    end    

    % atomic units
    Eh = 27.211386; % eV/Ha
    a0 = 5.291E-11; % a0/m

    % potential profile - last entry runs from r_(N-1) -> infty
    % if V = 0, assume an infinite step
    V = args.V;
    if isscalar(V);V = 0;end
    r = args.r; % last r is infinity (always one less element than V)
    me = args.me; % effective mass in each region
    er = args.er; % relative dielectric constant
    
    % eigenfunctions to calculate
    lmax = args.lmax;
    
    % some sanity checks
    assert(length(V) == length(me),'Number of potential steps must equal number of masses!')
    % turn off integrals
    noint = args.noint;

    % atomic units scaling
    V = V/Eh;
    r = cumsum(r)*1E-9/a0;
    
    if isempty(args.resume)
        % generate file name
        if isscalar(V)
            % infinite well
            file_name = ['8f' ...
                    '_l' strrep(num2str(lmax),' ','') ...
                    '_' job_title '_' tor];
        else
            % stepped potential
            file_name = ['sp' ...
                    '_l' regexprep(num2str(lmax),'\s+','') ...
                    '_V' regexprep(num2str(V),'\s+','')...
                    '_r' regexprep(num2str(r*a0/1E-9),'\s+','')...
                    '_me' regexprep(num2str(me),'\s+','')...
                    '_er' regexprep(num2str(er),'\s+','')...
                    '_' job_title '_' tor];
        end
        if ~isempty(job_title)
            % running on cluster, save to folder, else save it to current
            % working directory
            file_name = fullfile('..','basissets_8fold',file_name);
        end
    end
    
    % try to get walltime assigned format: days-hrs
    walltime = getenv('SLURM_TIME');
    if isempty(walltime)
        walltime = '1-0';
    end
    walltime = sum(cellfun(@str2num,strsplit(walltime,'-')).*[24 1]*3600);

    %% nlm - create three vectors for n, l, and m
    [n,l,u,~,En,ABknl] = pis_AO(lmax);
    n = n.';
    l = l.';
    u = u.';
    knl = ABknl(:,1).';
    Nnl = ABknl(:,2).';

    % how many atomic orbitals
    N = uint64(numel(n));
    fprintf('%g basis functions!\n',N)
    args.N = N;

    N8 = nchoosek(nchoosek(N+1,2)+1,2);
    fprintf('%g 1e- integrals\n',nchoosek(N+1,2));
    fprintf('%g 2e- integrals (8 fold symmetry in indices)\n',N8);
    fprintf('%g possible 2e- integrals\n',N^4);
    
    % complex -> real transformation matrices
    % lmax = l2+l3
    Umu = cell(2*max(l)+1,1);
    for a = 0:2*max(l)
        Umu{a+1} = Ulmu((-a:a).',-a:a);
    end
    
    %% calculate eri of form (12|34), or in coordinates, (1*1|2*2)
    eris_per_piece = uint64(min([N8/10 args.eris_per_piece]));
    
    if isempty(args.resume)
        fprintf('Fresh calculation\n')
        % if not resuming, start from the beginning
        eri_done = 0;
        integrals_done = 0;
        unique_integrals_done = 0;
        N8_ind = 1;
        % initialize piecewise indices generator
        next = uint64([1 1 1 1 1 1]);
        
        % create matfile to enable saving eri's directly to disk
        save_time = tic;
        m = matfile(file_name,'Writable',true);
        m.eri = zeros(N8,1);
        m.args = args;
        m.N8_ind = N8_ind;
        fprintf('Creating eri''s on disk... ');
        toc(save_time);
    else
        % resume a calculation
        eri_done = m.eri_done;
        integrals_done = m.integrals_done;
        unique_integrals_done = m.unique_integrals_done;
        next = m.next;
        N8_ind = m.N8_ind;
        fprintf(['Resuming prior calculation at N8_ind: %u, next: ' sprintf('%u ',next) '\n'],N8_ind)
    end
    
    % dump flags
    disp(file_name)
    fprintf('eris_per_piece: %g, size (MB): %g \n\n',eris_per_piece,double(eris_per_piece)*8/1024^2);
    disp(args)
    
    parallelOpen
    run_time = tic;
    last_time = 0;
    while true
        % save state of calculation in case time runs out
        m.eri_done = eri_done;
        m.integrals_done = integrals_done;
        m.unique_integrals_done = unique_integrals_done;
        m.next = next;
        m.N8_ind = N8_ind;
        % Gracefully exit if there isn't enough time
        if walltime - ceil(toc(run_time)) < 1.5*last_time
            fprintf(['Out of time, graceful exit at N8_ind: %u, next: ' sprintf('%u ',next) '\n'],N8_ind);
            parallelClose
            return
        end
        
        % some counters
        unique_integrals_sec = 0;
        integrals_sec = 0;
        
        % generate a chunk of indices
        ind_time = tic;
        [csec,next] = pis_8f_pw2_mex(N,eris_per_piece,next);
        u_csec = u(csec);
        l_csec = l(csec);
        
        % pre-evaluate radial integrals - they do not depend on magnetic
        % angular momentum, which contain a lot of degeneracy
        % C = A(IA,:) and A = C(IC,:)
        % IA grabs the unique values for the radial integral
        % IC maps back onto the full eri w/ degeneracy
        [~,IA,IC] = unique([n(csec) l_csec],'rows');
        csec_rad = csec(IA,:);
        rad_int = cell(numel(IA),1);
        fprintf('\n  Make a chunk of indices and find uniques (%.4g%% density): ',...
            double(length(IA))/double(eris_per_piece)*100)
        toc(ind_time)
        
        % First pass: evaluate unique radial integrals assuming the lower
        % limit of k is always l(1)-l(2) - ignore the u dependence
        rad_time = tic;
        knl_rad = knl(csec_rad);
        Nnl_rad = Nnl(csec_rad);
        l_rad = l(csec_rad);
        parfor b = 1:numel(rad_int)
            l_ = l_rad(b,:);
            knl_ = knl_rad(b,:);
            Nnl_ = Nnl_rad(b,:);
            N1234 = Nnl_(1)*Nnl_(2)*Nnl_(3)*Nnl_(4);
            % parity of l over b^k(lab) and b^k(lcd)
            lmax_12 = l_(1)+l_(2);
            lmax_34 = l_(3)+l_(4);
            if bitget(lmax_12,1) == bitget(lmax_34,1)
                % k is l in R^l and b^l
                lmin_12 = abs(l_(1)-l_(2));
                lmin_34 = abs(l_(3)-l_(4));
                %k = intersect(lmin_12:2:lmax_12,...
                %              lmin_34:2:lmax_34).';
                lrange_12_full = lmin_12:2:lmax_12;
                k = lrange_12_full(ismembc(lrange_12_full,lmin_34:2:lmax_34));
                if ~isempty(k)
                    r12 = zeros(numel(k),1);
                    for c = 1:numel(k)
                        % calculate R^k
                        if ~noint
                            r12(c) = integral2(@(r1,r2) RRRR(k(c),l_,knl_,N1234,r1,r2),0,1,0,@(r1)r1) + ...
                                     integral2(@(r1,r2) RRRR(k(c),flip(l_),flip(knl_),N1234,r1,r2),0,1,0,@(r1)r1);
                        else
                            r12(c) = 1;
                        end
                        unique_integrals_done = unique_integrals_done + 1;
                        unique_integrals_sec = unique_integrals_sec + 1;
                    end
                    rad_int{b} = r12;
                end
            end 
        end
        fprintf('  Integrate unique radial integrals, ')
        toc(rad_time)
        
        % Second pass: evaluate magnetic angular momentum
        ang_time = tic;
        rad_int = rad_int(IC);
        eri_csec = zeros(size(csec,1),1);
        parfor b = 1:size(csec,1)
            rad_int_ = rad_int{b};
            if any(rad_int_)
                u_ = u_csec(b,:).';
                l_ = l_csec(b,:).';
                % parity of l over b^k(lab) and b^k(lcd)
                lmax_12 = l_(1)+l_(2);
                lmax_34 = l_(3)+l_(4);
                if bitget(lmax_12,1) == bitget(lmax_34,1)
                    lmin_12 = max(abs(l_(1)-l_(2)),min(abs(u_(1)+u_(2)),abs(u_(1)-u_(2))));
                    lmin_12 = lmin_12 + bitget(lmin_12+lmax_12,1);
                    lmin_34 = max(abs(l_(3)-l_(4)),min(abs(u_(3)+u_(4)),abs(u_(3)-u_(4))));
                    lmin_34 = lmin_34 + bitget(lmin_34+lmax_34,1);
                    % k is l in R^l and b^l
                    %k = intersect(lmin_12:2:lmax_12,...
                    %              lmin_34:2:lmax_34).';
                    lrange_12 = lmin_12:2:lmax_12;
                    k = lrange_12(ismembc(lrange_12,lmin_34:2:lmax_34)).'; % needs to be a column vector
                    if ~isempty(k)
                        % previous loop calculated radial integral for all
                        % possible values of k. subset them out here.
                        %k_full = intersect(abs(l_(1)-l_(2)):2:lmax_12,...
                        %                   abs(l_(3)-l_(4)):2:lmax_34).';
                        lrange_12_full = abs(l_(1)-l_(2)):2:lmax_12;
                        % k_full is the naive range assumed in precalculating R
                        k_full = lrange_12_full(ismembc(lrange_12_full,abs(l_(3)-l_(4)):2:lmax_34));
                        k_ind = ismembc(k_full,k);
                        R = rad_int_(k_ind);
                        
                        % precalculate all allowed values of u assuming the
                        % maximum value of -u:u
                        u_12 = [u_(1)+u_(2) -(u_(1)+u_(2)) u_(1)-u_(2) -(u_(1)-u_(2))];
                        u_34 = [u_(3)+u_(4) -(u_(3)+u_(4)) u_(3)-u_(4) -(u_(3)-u_(4))];
                        ku_range = -k(end):k(end);
                        u_12i = ismembc(u_12,ku_range);
                        u_34i = ismembc(u_34,ku_range);
                        u_full = u_12(u_12i);
                        u_full = u_full(ismembc(u_full,sort(u_34(u_34i))));
                        u_full = unique(u_full);
                        
                        U = zeros(numel(k),1);
                        if ~isempty(u_full)
                            for c = 1:numel(k)
                                % for each l, check which \mu = -l,...,l are common
                                % with b^k(lab) and b^k(lcd)
                                %u_i = intersect(intersect(-k(c):k(c),u_12),intersect(-k(c):k(c),u_34));

                                % subset u_i's
                                u_i = u_full(abs(u_full) <= k(c));

                                if ~isempty(u_i)
                                    % for each u_i(ntersect)
                                    for d = 1:numel(u_i)
                                        % calculate b^k(uab)
                                        l_s = [k(c);l_(1:2)];
                                        u_s = [u_i(d);u_(1:2)];
                                        %bkab = beekay(l_s,u_s);
                                        bkab = beekay2(l_s,u_s,Umu);
                                        % b^k(ucd)
                                        if bkab ~= 0
                                            l_s(2:3) = l_(3:4);
                                            u_s(2:3) = u_(3:4);
                                            %bkcd = beekay(l_s,u_s);
                                            bkcd = beekay2(l_s,u_s,Umu);
                                            U(c) = U(c) + bkab*bkcd;
                                        end
                                    end
                                    if U(c) ~= 0
                                        % i.e. if you did a radial integration for 
                                        % each angular integration
                                        integrals_done = integrals_done + 1;
                                        integrals_sec = integrals_sec + 1;
                                    end
                                end
                            end
                        end
                        eri_csec(b) = sum(R.*U);
                        eri_done = eri_done + 1;
                    end
                end
            end
        end
        fprintf('  Angular integrations: ');
        toc(ang_time)
        % insert them into vector
        save_time = tic;
        if N8_ind + eris_per_piece > N8
            fprintf('*** Last Segment! ***\n')
            m.eri(N8_ind:N8,1) = eri_csec;
        else
            m.eri(N8_ind:N8_ind+eris_per_piece-1,1) = eri_csec;
        end
        fprintf('  Save chunk of eris: ')
        toc(save_time)
        fprintf('  %u unique integrals and %u lazy integrals in this section\n',unique_integrals_sec,integrals_sec)
        fprintf('  average time (us) per eri in this chunk: %g\n',toc(ind_time)/double(eris_per_piece)*1E6)
        N8_ind = N8_ind + eris_per_piece;
        if N8_ind > N8 || all(next==0)
            fprintf('100%%! ');
            toc(ind_time);
            break
        end
        fprintf(['%.4g%% done, next: ' sprintf('%u ',next) ', '],double(N8_ind)/double(N8)*100);
        toc(ind_time)
        % use the time of this run to guess the time of the next run (to
        % see if we need to terminate gracefully)
        last_time = toc(ind_time);
    end
    parallelClose
    
    % diagnostics
    fprintf('\n\neris done: %u\n',eri_done)
    fprintf('lazy, full 4 function, 1d integrals: %u\n',integrals_done);
    fprintf('unique 4 function, 1d integrals done: %u\n',unique_integrals_done);
    fprintf('fraction of eris done: %.4g\n',eri_done/double(N8))

    %% save it to hdf5 file
    % only bother with the 1 electron matrices if the eri is finished
    m.unique_integrals_done = unique_integrals_done;
    m.integrals_done = integrals_done;
    m.eri_done = eri_done;
    m.next = next;
    m.Hcore = diag(En);
    m.ovlp = eye(N);
    % commented out bottom line .. most likely because code up top takes
    % its place?
    %save(file_name,'args','Hcore','ovlp','-v7.3')
end