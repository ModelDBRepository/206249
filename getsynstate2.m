function [brweights, branch_syns, nrnweights, nrn_syns, brstrengths, brsynratio] = getsynstate2(fn, npyrs, nbranches, ninputs)
        %defaults
        brweights = zeros(ninputs, npyrs*nbranches);
        nrnweights = zeros(ninputs, npyrs);
        branch_syns = zeros(ninputs, npyrs*nbranches);
        nrn_syns = zeros(ninputs, npyrs);
        %ff = sprintf('./data/%s_%d_%d/synstate.dat', CONDITION, ncase, run-1)
        ss = load(fn);

        for i=1:size(ss,1)
            bid=ss(i,2);
            nid=ss(i,3);
            srcid=ss(i,5);
            bstrength = ss(i,6);
            
            w=ss(i,7);
            

            if ((srcid >=0) && (bid <= npyrs*nbranches))
                brweights(srcid+1, bid+1) = brweights(srcid+1,  bid+1) + w;
                brstrengths(srcid+1, bid+1)=bstrength;
                nrnweights(srcid+1,   nid+1) = nrnweights(srcid+1,  nid+1) + w;
            end
            if (srcid >=0 && bid <= npyrs*nbranches &&  w > 0.7)
                branch_syns(srcid+1,  bid+1) = branch_syns(srcid+1, bid+1)+1;
                nrn_syns(srcid+1,  nid+1) = nrn_syns(srcid+1, nid+1)+1;
            end

        end

end
