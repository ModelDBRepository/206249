defaults
nruns=10;

totfiring = zeros(nruns, 4);
totactive = zeros(nruns, 4);

coract = zeros(nruns, 9);

%CUTOFF=5; % Hz

npatterns=1;
ndays  = 4;

brweights = zeros(ndays, npyrs*nbranches, nruns);
nrnweights = zeros(ndays, npyrs, nruns);
branch_syns = zeros(ndays, npyrs*nbranches, nruns);
brsynratio= zeros(ndays, nruns);



for run=1:nruns
        for ncase=1:ndays
            
            sfn=sprintf('./data/%s_%d_%d/spikesperpattern.dat', CONDITION, ncase, run-1);
            spk = load( sfn);
            spk = spk(1, 1:npyrs)/(stimduration/1000);
            
            %pop = spk(spk>=CUTOFF);
            pop = spk(spk>=CUTOFF);
            totfiring(run, ncase) = mean(spk, 2);
            totactive(run, ncase) = sum(spk>CUTOFF,2)/npyrs;



            ff = sprintf('./data/%s_%d_%d/synstate.dat', CONDITION, ncase, run-1);
            ss = load(ff);

            for i=1:size(ss,1)
                bid=ss(i,2);
                nid=ss(i,3);
                srcid=ss(i,5);
                bstrength = ss(i,6);
                w=ss(i,7);
                if (srcid ==0 && bid <= npyrs*nbranches)
                    brweights( ncase, bid+1, run) = brweights(ncase, bid+1, run) + w;
                    brstrengths(ncase, bid+1)=bstrength;
                    nrnweights( ncase, nid+1,run) = nrnweights(ncase, nid+1,run) + w;
                end
                if (srcid ==0 && bid <= npyrs*nbranches &&  w > 0.7)
                    branch_syns(ncase, bid+1, run) = branch_syns(ncase, bid+1, run)+1;
                end
            end
            
            brsynratio(ncase,run) = sum(branch_syns(ncase,:, run)>3)/(nbranches*npyrs);
        end
end

close all
mf = mean(totfiring,1);
sf = std(totfiring,0,1)/sqrt(nruns);

mact = mean(totactive,1);
sact = std(totactive,0,1)/sqrt(nruns);

barwitherr(100.* sact(1,:,1), 100.* mact(1,:,1));
title('% coding neurons')
%ylabel('% Active neurons')
xlabel('Day')
ylim([0,80]);

export_fig(sprintf('./figs/%s_pops.pdf',CONDITION), '-transparent')

figure
barwitherr(sf(1,:,1), mf(1,:,1));

title('Average firing rate [Hz]')
%ylabel('Avg firing rate [Hz]')
%xlabel('Number of trainings')
%ylim([0,70]);
export_fig(sprintf('./figs/%s_rates.pdf',CONDITION), '-transparent')

figure
hs=hist(branch_syns(1, :), [0:8]);

bar(hs(:,2:end)/sum(hs(:)))
title('1st day')
ylim([0,0.2]);
export_fig(sprintf('./figs/%s_brsyns_1day.pdf',CONDITION), '-transparent')

figure
hs=hist(branch_syns(4, :), [0:8]);
bar(hs(:,2:end)/sum(hs(:)))
ylim([0,0.2]);
title('4th day')
export_fig(sprintf('./figs/%s_brsyns_4day.pdf',CONDITION), '-transparent')


figure
barwitherr(100.0*std(brsynratio,0,2)/sqrt(nruns), 100.0*mean(brsynratio,2))
title('Branches with >2 potentiated synapses')
ylabel('Percentage')
xlabel('Day')
ylim([0,16]);
export_fig(sprintf('./figs/%s_brsyns.pdf',CONDITION), '-transparent')


figure
aa = mean(branch_syns(:,:)');
ss = std(branch_syns(:,:)');
errorbar(ss, aa, 'o')
title('Potentiated synapses per branch');
ylabel('Number of synapses')
xlabel('Day')
ylim([0,6]);
export_fig(sprintf('./figs/%s_syn_per_branch.pdf',CONDITION), '-transparent')



figure
tweights = (squeeze(sum(nrnweights,2)));

barwitherr(std(tweights,0,2), mean(tweights,2))

%barwitherr(100.0*std(brsynratio,0,2)/sqrt(nruns), 100.0*mean(brsynratio,2))
title('Total Syn Weights')
ylabel('Total Syn Weight')
xlabel('Day')
ylim([0,12000]);
export_fig(sprintf('./figs/%s_tweights.pdf',CONDITION), '-transparent')


CONDITION
b1 = branch_syns(1,:);
b1 = b1(b1>0);
b4 = branch_syns(4,:);
b4 = b4(b4>0);
mean(b1)
std(b1)
mean(b4)
std(b4)
[h,p] = ttest2(b1, b4)