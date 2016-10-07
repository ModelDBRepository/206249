    
defaults
%set(0,'DefaultAxesFontSize', 8);
%set(0,'DefaultTextFontSize', 8);


nruns=10;
totfiring = zeros(nruns, 8, 2);
totactive = zeros(nruns, 8, 2);

coract = zeros(nruns, 8);
coractPi = zeros(nruns, 8);
overlap = zeros(nruns, 8);


%CUTOFF=5; % Hz

diffs = [ -1440,  -300,  -120, -60, 60,  120, 300, 1440];


difflabels = { '-24H', '-5H', '-2H', '-1H', '+1H','+2H', '+5H', '+24H'};
%CONDITION='weakstrongL'


brsyns = zeros(nruns,npyrs*nbranches);
brsynsH = zeros(nruns,10);

nrnsyns = zeros(nruns,npyrs);

nrnbins = [0:40];
nrnsynsH = zeros(nruns, length(nrnbins));

brcors = zeros(nruns,7);
nrncors = zeros(nruns, 7);

brcommon = [];
brtsyns = [];
nrntsyns = [];
%nruns=10;
for run=1:nruns
        for ncase=1:length(diffs)
            
            diff = diffs(ncase)
            
            if (diff < 0)
                p1 = 0;
                p2 = -diff;
            else
                p1 = diff;
                p2 = 0;
            end
            sfn=sprintf('./data/%s_%d_%d_%d/spikesperpattern.dat', CONDITION, p2, p1,run-1);
            
            spk = load( sfn);

            
            weakmem = 2;
            if (diff>0)
                spk = flipud(spk);
                weakmem = 1;
            end
            spk = spk(:, 1:npyrs)/(stimduration/1000);
            
            totfiring(run, ncase, 1) = mean(spk(1, spk(1,:)>0), 2);
            totfiring(run, ncase, 2) = mean(spk(2, spk(2,:)>0), 2);
            totactive(run, ncase, :) = sum(spk>CUTOFF,2);
            
            %s = corrcoef(spk');
            
            
            [s,p] = corrcoef(spk');
            
            coract(run, ncase) = s(1,2);
            coractPi(run, ncase) = p(1,2);

            spcut = spk >CUTOFF;
            overlap(run, ncase) = sum(spcut(1,:)&spcut(2,:) )/( sum(spcut(:))/2 );

            [bw, bs, nw, ns] = getsynstate(sprintf('./data/%s_%d_%d_%d/synstate.dat', CONDITION, p2, p1,run-1));
            
            brcors(run,ncase) = corr(bw(1,:)', bw(2,:)');
            nrncors(run,ncase) = corr(nw(1,:)', nw(2,:)');
            brsyns = bs(1,:);
            brsyns2 = bs(2,:);
            
            brcommon(ncase, run) =  sum( (brsyns>1) & (brsyns2>1) ,2)/sum( (brsyns>1) | (brsyns2>1) ,2);
            
            brtsyns( ncase, run) = sum(bs(weakmem,:)>0);
            nrntsyns(ncase, run) = sum(ns(weakmem,:)>0);


            % brsynsH(run,:) = hist(bs(1,:), [0:9]);
            % brsynsH(run,:) = brsynsH(run,:)/sum(brsynsH(run,:));

            % nrnsyns(run,:) = ns(1,:);
            % nrnsynsH(run,:) = hist(ns(1,:), nrnbins);
            % nrnsynsH(run,:) = nrnsynsH(run,:)/sum(nrnsynsH(run,:));

        end
end



close all
mf = mean(totfiring,1);
sf = std(totfiring,0,1);

figure
barwitherr(sf(1,:,1), mf(1,:,1));
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', diffs)
title('Strong Memory F')
%ylabel('Avg firing rate [Hz]')
%xlabel('Weak-Strong Interval [minutes]')
ylim([0,50])
xlim([0,13])
savefig(sprintf('./figs/%s_fstrong',CONDITION), 'pdf')



figure
barwitherr(sf(1,:,2), mf(1,:,2));
results([CONDITION, 'ff']) = [ sf(1,:,:); mf(1,:,:) ];
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', diffs)
title('Weak Memory F')
%ylabel('Avg firing rate [Hz]')
%xlabel('Weak-Strong Interval [minutes]')
ylim([0,50]);
xlim([0,13]);
savefig(sprintf('./figs/%s_fweak',CONDITION), 'pdf')

figure
mact = 100.*mean(totactive,1)/npyrs
sact = (100.*std(totactive,0,1)/npyrs)/sqrt(nruns)


results([CONDITION, 'act']) = [ sact; mact];

results([CONDITION, 'act_d']) = [ 100.*(totactive)/npyrs ];

barwitherr(sact(1,:,1), mact(1,:,1));
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', diffs)
title('Strong Memory')
%ylabel('% Active neurons')
%xlabel('Weak-Strong Interval [minutes]')
ylim([0,35])
xlim([0,13])
savefig(sprintf('./figs/%s_actstrong',CONDITION), 'pdf')

figure
barwitherr(sact(1,:,2)/sqrt(nruns), mact(1,:,2)); % weak
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', diffs)
title('Weak Memory')
ylabel('% Active neurons')
%xlabel('Weak-Strong Interval [minutes]')
ylim([0,35])
xlim([0,13])
savefig(sprintf('./figs/%s_actweak',CONDITION), 'pdf')


results([ CONDITION, 'firingcor']) = [ 100.*std(coract, 0,1) ; 100.*mean(coract) ];
results([ CONDITION, 'firingcor_d']) = [  100.*(coract) ];
results([ CONDITION, 'firingcor_pi']) = [  100.*(coractPi) ];


figure
results([ CONDITION, 'common']) = [ 100.*std(overlap, 0,1) ; 100.*mean(overlap) ];
results([ CONDITION, 'common_d']) = [  100.*(overlap) ];


barwitherr(std(overlap, 0,1)/sqrt(nruns) , mean(overlap));
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', diffs)
title('Correlation between population firing rates')
%ylabel('Correlation coefficient')
%xlabel('Weak-Strong Interval [minutes]')
%ylim([0,1.0])
%xlim([0,13])

set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', difflabels)
rotateXLabels(gca, 60)
xlim([0,length(diffs)+1])    
ylim([0,1]);

savefig(sprintf('./figs/%s_corr',CONDITION), 'pdf');
figure

results([ CONDITION, 'brcors']) = [ std(brcors, 0,1); mean(brcors) ];
results([ CONDITION, 'brcors_d']) = (brcors) ;


barwitherr(std(brcors, 0,1)/sqrt(nruns), mean(brcors));
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', diffs)
title('Correlated synapses per branch')
%ylabel('Correlation')
%xlabel('Weak-Strong Interval [minutes]')
ylim([0,1.0])
xlim([0,13])
savefig(sprintf('./figs/%s_brcorr',CONDITION), 'pdf')


figure
results([CONDITION 'nrncors']) = [std(nrncors, 0,1); mean(nrncors)];
results([CONDITION 'nrncors_d']) = nrncors;


barwitherr(std(nrncors, 0,1), mean(nrncors));
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', diffs)
title('Correlated synapses per neuron')
%ylabel('Correlation')
%xlabel('Weak-Strong Interval [minutes]')
ylim([0,1.0])
xlim([0,13])
savefig(sprintf('./figs/%s_nrnrcorr',CONDITION), 'pdf')



