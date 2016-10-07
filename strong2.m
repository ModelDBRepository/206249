defaults
close all
nruns=10;

diffs = [ 60, 120, 300, 1440 ];

difflabels = { '1H','2H', '5H', '24H'};

ndiffs = length(diffs);

totcommon = zeros(nruns, ndiffs);

totfiring = zeros(nruns, ndiffs, 2);
totactive = zeros(nruns, ndiffs, 2);

overlap = zeros(nruns, ndiffs);
coract = zeros(nruns, ndiffs);

brsyns = zeros(nruns,npyrs*nbranches);
brsynsH = zeros(nruns,10);

nrnsyns = zeros(nruns,npyrs);

nrnbins = [0:40]
nrnsynsH = zeros(nruns, length(nrnbins));

brcors = zeros( nruns, ndiffs);
nrncors = zeros(nruns, ndiffs);


histCSUSrange = [1:50];
histCSUS = zeros(nruns, length(diffs), length(histCSUSrange));

brcommon = [];

for ncase=1:length(diffs)
        for run=1:nruns
            diff = diffs(ncase);
            
            if (diff < 0)
                p1 = 0;
                p2 = -diff;
            else
                p1 = diff;
                p2 = 0;
            end
            
            sfn=sprintf('./data/%s_%d_%d/spikesperpattern.dat', CONDITION, p1,run-1)
            
            spk = load( sfn);

            if (diff<0)
                spk = flipud(spk);
            end
            
            spk = spk(:, 1:npyrs)/(stimduration/1000); %duration
            
            
            totfiring(run, ncase, 1) = mean(spk(1, :), 2);
            totfiring(run, ncase, 2) = mean(spk(2, :), 2);
            totactive(run, ncase, :) = sum(spk>CUTOFF,2)/npyrs;
            
            
            ffactive(run, ncase, 1) = mean(spk(1, spk(1,:)>CUTOFF), 2);
            ffactive(run, ncase, 2) = mean(spk(2, spk(2,:)>CUTOFF), 2);
            
            
            spk1 = (spk(1,:)>=CUTOFF);
            spk2 = (spk(2,:)>=CUTOFF);
            totcommon(run, ncase) = sum(spk1 & spk2)/( sum(spk1)+sum(spk2));
            
            
            s = corrcoef(spk');
            coract(run, ncase) = s(1,2);

            overlap(run, ncase) = corr( (spk(1,:)>=CUTOFF)' ,  (spk(2,:)>=CUTOFF)' );
            
            if (1)
                [bw, bs, nw, ns] = getsynstate(sprintf('./data/%s_%d_%d/synstate.dat', CONDITION, p1,run-1));

                %brcors(run,ncase) = corr(bw(1,:)', bw(2,:)');
                %nrncors(run,ncase) = corr(nw(1,:)', nw(2,:)');


                brsyns = bs(1,:);
                brsyns2 = bs(2,:);

                brcommon(run, ncase) =  sum( (brsyns>1) & (brsyns2>1) ,2)/ sum( (brsyns>1) | (brsyns2>1) ,2);
                
                brtot = brsyns+brsyns2;
                histCSUS(run, ncase, :)  = histc(brtot(find((brsyns>0) & (brsyns2>0))), [1:50]);
                
                
               % brsyns(run,:) = bs(1,:);
               % brsynsH(run,:) = hist(bs(1,:), [0:9]);
               % brsynsH(run,:) = brsynsH(run,:)/sum(brsynsH(run,:));

               % nrnsyns(run,:) = ns(1,:);
               % nrnsynsH(run,:) = hist(ns(1,:), nrnbins);
               % nrnsynsH(run,:) = nrnsynsH(run,:)/sum(nrnsynsH(run,:));
            end
        end
        
        
end

if (strcmp(CONDITION,'strong2L') || strcmp(CONDITION,'strong2NL'))
    COL='r';
else
    COL='b';
end

%close all
mf = mean(totfiring,1)
sf = std(totfiring,0,1)/sqrt(nruns)
figure

hold on
errorbar(mf(:,:,1), sf(:,:,1));
errorbar(mf(:,:,2), sf(:,:,2), 'g');
hold off

%barwitherr(sf(1,:,1), mf(1,:,1));
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', difflabels)
title(sprintf('Avg firing rate \nAll Neurons [Hz] - %s', CONDITION))
xlabel('Interval Between Memories [minutes]')

%ylim([0,10])

export_fig(sprintf('./figs/%s_ffirst.pdf',CONDITION), '-transparent')


figure
mact = 100.0*mean(totactive,1)
sact = 100.0*std(totactive,0,1)/(sqrt(nruns))


hold on
errorbar(mact(:,:,1), sact(:,:,1), COL);
errorbar(mact(:,:,2), sact(:,:,2), 'g');

hold off

%barwitherr(sact(1,:,1), mact(1,:,1));
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', difflabels)
title('% Active neurons');
%xlabel('Weak-Strong Interval [minutes]')
ylim([10,50])
export_fig(sprintf('./figs/%s_actfirst.pdf',CONDITION), '-transparent')



figure
barwitherr(std(coract, 0,1)/sqrt(nruns), mean(coract), COL);
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', difflabels)
title(sprintf('Similarity between population firing patterns %s ', CONDITION))
ylabel('Similarity')
%xlabel('Weak-Strong Interval [minutes]')
%ylim([0,.8])
export_fig(sprintf('./figs/%s_corr.pdf',CONDITION), '-transparent')





figure
barwitherr(std(overlap, 0,1)/sqrt(nruns), mean(overlap), COL);
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', difflabels)
title(sprintf('Active population similarity %s ', CONDITION))
ylabel('Similarity')
%xlabel('Weak-Strong Interval [minutes]')
ylim([-.1,1.0])
export_fig(sprintf('./figs/%s_overlap.pdf',CONDITION), '-transparent')

%close all
mf = mean(ffactive,1)
sf = std(ffactive,0,1)/sqrt(nruns)
figure


hold on
errorbar(mf(:,:,1), sf(:,:,1));
errorbar(mf(:,:,2), sf(:,:,2), 'g');
hold off
%barwitherr(sf(1,:,1), mf(1,:,1));
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', difflabels)
title('Active neurons firing rate [Hz]')
xlabel('Interval Between Memories [minutes]')
ylim([0,100])
export_fig(sprintf('./figs/%s_ffactive.pdf',CONDITION), '-transparent')


figure
barwitherr(100.*std(totcommon)/(sqrt(nruns)), 100.*mean(totcommon), COL);
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', difflabels)
title(sprintf('Common Active Neurons %s', CONDITION))
ylabel('% common neurons')
%xlabel('Weak-Strong Interval [minutes]')
ylim([0,50])
export_fig(sprintf('./figs/%s_totcommon.pdf',CONDITION), '-transparent')


figure
barwitherr(100.*std(brcommon)/(sqrt(nruns)), 100.*mean(brcommon), COL);
set(gca, 'XTick', [1:length(diffs)])
set(gca, 'XTickLabel', difflabels)
title('% Branches with clusters of both memories')
ylabel('% branches')
%xlabel('Weak-Strong Interval [minutes]')
ylim([0,50])
export_fig(sprintf('./figs/%s_brcommon.pdf',CONDITION), '-transparent')


figure;
cc=winter(6);
nn = [1,3,4];
color = [0,0,1];

for i=1:length(nn)
    ncase = nn(i);
    
    diff = diffs(ncase);
    
    aa  = mean(histCSUS(:,ncase,:));
    bb = stderr(histCSUS(:,ncase,:));
    aa = aa(:);
    
    bb =bb(:);
    errorbar(aa,bb, 'Color',  cc(ncase,:));
    %set(gca, 'Color', 'o');
    hold on
    xlim([0,15]);
end

ylabel('Number of clusters of both memories');
xlabel('Synapses per  cluster');

legend({'1 hour', '3 hours', '5 hours'});
export_fig(sprintf('./figs/%s_clustering2.pdf',CONDITION), '-transparent')
hold off;
