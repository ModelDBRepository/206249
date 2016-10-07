defaults;
npatterns=1;
actPpre= [];
actPpost = [];

histPpre = [];
histPpost = [];
histIpre = [];
histIpost = [];

prePop = [];;
postPop = [];


nruns=10;

brcase = [];
multiruns=1;
for ncase=1:4
    nbranches = ncase*10;
    brsynsH = zeros(nruns,10);
    brsynsK = zeros(nruns,10);
    nrnsyns = zeros(nruns,npyrs);

    nrnbins = [0:80];
    nrnsynsH = zeros(nruns, length(nrnbins));
    sumweights = zeros(nruns, 2);

    %%CONDITION='brtestL'

    brws = [];
    brsyns = zeros(nruns,npyrs*nbranches);
    nsbefore = [];
    nsafter= [];


    for run=1:nruns
        
 
        sfn=sprintf('./data/%s_%d_%d/spikes.dat', CONDITION, nbranches, run-1);
        fspk = fopen(sfn);
      
        nn =1;
        while ~feof(fspk)
           line = fgets(fspk);
           c = sscanf(line, '%d');

           for ns=1:multiruns
              tt= (ns-1)*stimduration;
              nsbefore(run, ns, nn) = length(find(c>tt & c<tt+stimduration));

              tt= (multiruns+ns)*stimduration;
              nsafter(run,  ns, nn)  = length(find(c>tt & c<tt+stimduration));
           end

           nn = nn + 1;
           if (nn>nneurons) 
               break;
           end
        end
        fclose(fspk);

        npat=0;
        %tend = npat*stimduration;
        %actpre  = sum( raster(1:npyrs, tend+1:tend+stimduration),2)/2; % for 2000 msec

        actpre = nsbefore(run, 1, 1:npyrs)/(stimduration/1000);

        %tend = 2*npatterns*stimduration+npat*stimduration;
        %actpost = sum( raster(1:npyrs, tend+1:tend+stimduration),2)/2; % for 2000 msec
        actpost = nsafter(run, 1, 1:npyrs)/(stimduration/1000);

        actPpre(ncase, run,  :) = actpre;
        actPpost(ncase, run, :) = actpost;

        [bw, bs, nw, ns] = getsynstate2(sprintf('./data/%s_%d_%d/synstate.dat', CONDITION, nbranches, run-1), npyrs, nbranches, ninputs);

        brws(run,:) = bw(1,:);
        brsyns(run,:) = bs(1,:);
       
        brsynsH(run,:) = hist(bs(1,:), [0:9]);
        brsynsK(run,:) = hist(bs(1,:), [0:9]);
        brsynsH(run,:) = brsynsH(run,:)/sum(brsynsH(run,:));

        nrnsyns(run,:) = ns(1,:);
        nrnsynsH(run,:) = hist(ns(1,:), nrnbins);
        nrnsynsH(run,:) = nrnsynsH(run,:)/sum(nrnsynsH(run,:));
        
    end
    brcase(ncase) = mean(brsyns(find(brsyns>0)));
    brcases(ncase) = std(brsyns(find(brsyns>0)))/sqrt(10);
end


ba = (mean(actPpre>CUTOFF,3));
bb = (mean(actPpost>CUTOFF,3));

figure
errorbar(mean(ba,2),std(ba,0,2));

hold on
errorbar(mean(bb,2),std(bb,0,2));
title(sprintf('Active population - %s',CONDITION))
legend('Pre', 'Post');



