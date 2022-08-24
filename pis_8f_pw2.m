function [C2,next] = pis_8f_pw2(norb,count,curr)
% piecewise evaluation of the 8f index according to pyscf convention
% this function returns a list of eri's that are indexed in the pyscf
% scheme: (ij|kl), i>=j, k>=l, ij>=kl
%
% this function is written for mex'ing by explicitly defining the variable
% types
%
% arguments: (norb,count,curr)
% norb = integer, number of orbitals, uint64
% count = integer, number of indices to generate in this function call, uint64
% curr = 1x6 vector containing the previous exit state of the function, uint64
% curr = [i j k l ij kl]
    
    stop = false;
    cont_jkl = true(3,1);
    cont_kl = true;
    
    % stop condition (overwritten if ijkl > count)
    next = uint64(zeros(1,6));
    
    C2 = uint16(zeros(count,4));
    ijkl = uint64(1);
    ij = curr(5);
    for i = curr(1):norb
        for j = curr(2):i
            if cont_kl
                kl = curr(6);
                cont_kl = false;
            else
                kl = uint64(1);
            end
            for k = curr(3):i
                for l = curr(4):k
                    if ij >= kl
                        C2(ijkl,:) = [i j k l];
                        ijkl = ijkl + uint64(1);
                    end
                    kl = kl + uint64(1);
                    if ijkl > count
                        next = [i j k l+uint64(1) ij kl];
                        stop = true;
                        break
                    end
                end
                if cont_jkl(3)
                    curr(4) = uint64(1);
                    cont_jkl(3) = false;
                end
                if stop
                    break
                end
            end
            if cont_jkl(2)
                curr(3) = uint64(1);
                cont_jkl(2) = false;
            end
            if stop
                break
            end
            ij = ij + uint64(1);
        end
        if stop
            break
        end
        if cont_jkl(1)
            curr(2) = uint64(1);
            cont_jkl(1) = false;
        end
    end
    
    if ~stop
        % if the loop finished, prune zeros - ijkl increments by 1 even
        % after the loop is done in principle
        C2(ijkl:end,:) = [];
    end
end