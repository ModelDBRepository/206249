defaults

npatterns=1
actPpre= []
actPpost = []

histPpre = []
histPpost = []
histIpre = []
histIpost = []

prePop = []
postPop = []
brws = [];


nruns=10

brsyns = zeros(nruns,npyrs*nbranches);
brsynsH = zeros(nruns,10);
brsynsK = zeros(nruns,10);
nrnsyns = zeros(nruns,npyrs);
csusbefore = zeros(nruns,1);
csusafter = zeros(nruns,1);

nrnbins = [0:80];
nrnsynsH = zeros(nruns, length(nrnbins));
sumweights = zeros(nruns, 2);

iactPpre = [];
iactPpost = [];
%CONDITION='sparseL'

%close all

nsbefore = [];
nsafter= [];

multiruns = 1;

synapsesCSUS = zeros(nruns, 10);


for run=1:nruns
    
    raster = zeros(npyrs, stimduration*npatterns*3);
    avraster = zeros(npyrs, stimduration*npatterns*3);
    sfn=sprintf('./data/%s_%d/spikes.dat', CONDITION, run-1);
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


    if (0)
        nfilts=40;
        filt = ones(1,nfilts);
        filt=filt/nfilts;

        nn=1;
        while ~feof(fspk)
           line = fgets(fspk);
           c = sscanf(line, '%d');
           raster(nn, c') = 1;
           ab=filter( filt, [1], raster(nn,:));
           avraster(nn,:) = 1000*ab;
           nn = nn + 1;
           if (nn>nneurons) 
               break;
           end
        end
        fclose(fspk);

        figure;
        imagesc(avraster(1:160,1:stimduration), [0,200]);
        xlabel('Time [msec]')
        ylabel('Neuron #')
        colorbar();
        savefig(sprintf('./figs/ENC_%s_rates_pre', CONDITION), 'pdf')
        figure();
        imagesc(avraster(1:160,2*stimduration:2*stimduration+stimduration) , [0,200]);
        colorbar()
        xlabel('Time [msec]')
        ylabel('Neuron #')
        savefig(sprintf('./figs/ENC_%s_rates_post', CONDITION), 'pdf')

        myit = 50;
        bins = [];

        for n=1:myit:size(raster,2)
            ss = sum(sum(raster(1:160,n:n+myit-1)))
            bins(end+1) =ss;
        end

        figure;
        bar(bins(4000/myit:6000/myit), 'r');
        hold on
        bar(bins(1:2000/myit));
        hold off
        ivv = 8;
        set(gca,'Xtick', 0:ivv:2000/(myit), 'XTickLabel',[0:ivv*myit:2000]);
        xlim([0,2000/myit+1])
        %ylim([0,85])
        ylabel('Spike count')
        xlabel('Time [msec]')
        savefig(sprintf('./figs/ENC_%s_inst_rates', CONDITION), 'pdf')
   
    end
    
    npat=0;
    %tend = npat*stimduration;
    %actpre  = sum( raster(1:npyrs, tend+1:tend+stimduration),2)/2; % for 2000 msec
    actpre = nsbefore(run, 1, 1:npyrs)/(stimduration/1000);
    %tend = 2*npatterns*stimduration+npat*stimduration;
    %actpost = sum( raster(1:npyrs, tend+1:tend+stimduration),2)/2; % for 2000 msec
    actpost = nsafter(run, 1, 1:npyrs)/(stimduration/1000);

    actPpre(run, :) = actpre(:);
    actPpost(run, :) = actpost(:);
        
       
    %tend = npat*stimduration;
    %iactpre  = sum( raster(npyrs:nneurons, tend+1:tend+stimduration),2)/2;    
    iactpre = nsbefore(run, 1, npyrs:nneurons);
    
    
    %tend = 2*npatterns*stimduration+npat*stimduration;
    %iactpost = sum( raster(npyrs:nneurons, tend+1:tend+stimduration),2)/2;
    iactpost = nsbefore(run, 1, npyrs:nneurons);

    iactPpre(run, :) = iactpre;
    iactPpost(run, :) = iactpost;

    [bw, bs, nw, ns] = getsynstate(sprintf('./data/%s_%d/synstate.dat', CONDITION, run-1));
    
    
    brws(run,:) = bw(1,:);
    brsyns(run,:) = bs(1,:);
    brsynsH(run,:) = hist(bs(1,:), [0:9]);
    brsynsK(run,:) = hist(bs(1,:), [0:9]);
    brsynsH(run,:) = brsynsH(run,:)/sum(brsynsH(run,:));
    
    
    nn = ns(1,:);
    nrnsyns(run,:) = ns(1,:);
    nrnsynsH(run,:) = hist(ns(1,:), nrnbins);
    nrnsynsH(run,:) = nrnsynsH(run,:)/sum(nrnsynsH(run,:));

    l1 = load(sprintf('./data/%s_%d/sum-weights.txt', CONDITION, run-1));
    sumweights(run, 1) = l1(1,1);
    sumweights(run, 2) = l1(2,1);

    
    dta = load(sprintf('./data/%s_%d/syn-pre.txt', CONDITION, run-1));
    dta= dta(find(dta(:,1)==0), :);
    bwCS = dta(find(dta(:,2)<3), 3);
    bwUS = dta(find(dta(:,2)>=3), 3);
    csusbefore(run) = length(intersect(bwCS, bwUS))
    
    
    dta = load(sprintf('./data/%s_%d/syn-post.txt', CONDITION, run-1));
    dta= dta(find(dta(:,1)==0), :);
    dta= dta(find(dta(:,5)>0.7), :);
    bwCS = dta(find(dta(:,2)<3), 3);
    

    csusafter(run) = length(intersect(bwCS, bwUS));

    bwUS = dta(find(dta(:,2)>=3), 3);
     
    brstim = zeros(3200, 5);
    unCS = unique(bwCS);
    unUS = unique(bwUS);
    brstim(unCS+1, 1) = histc(bwCS, unCS);
    brstim(unUS+1, 2) = histc(bwUS, unUS);
    
    brstim = brstim(find(brstim(:,1)>0),:);
    brstim = brstim(find(brstim(:,2)>0),:);
    synapsesCSUS(run, :) = histc(brstim(:,1)+brstim(:,2), [1:10]);   % total number of CS+US synapses in  branches that contain at least 1 CS and 1 US
end


nextplot(2,2)
b = hist(actPpre(:), [0:19]);
h=bar(b/sum(b(:)));
set(h,'edgecolor' , 'none')
title('Pre - Excitatory')
%ylim([0,.8]);
%xlabel('Avg Firing Rate [Hz]')
%ylabel('Probability')

nextplot
b = hist(actPpost(:), [0:19]);
h=bar(b/sum(b(:)))
set(h,'edgecolor' , 'none')
%ylim([0,.8]);
title('Post - Excitatory')

nextplot
b = hist(iactPpre(:), [0:40]);
h=bar(b/sum(b(:)))
set(h,'edgecolor' , 'none')
%ylim([0,0.2])
title('Pre - Inhibitory')


nextplot
b = hist(iactPpost(:), [0:40]);
h=bar(b/sum(b(:)))
set(h,'edgecolor' , 'none')
%ylim([0,0.2])
title('Post - Inhibitory')

export_fig(sprintf('./figs/%s_rates.pdf',CONDITION), '-transparent')


nextplot(1,1)
set(gcf, 'Position', [0,0, 440,300])
dta = actPpre(:);
ad =sort(dta,1,'descend');
bar(ad)
kpre = trevrolls(ad);
title('Pre-training firing rates')

ylabel('Average Firing Rate [Hz]');
xlabel('Excitatory neuron')
ylim([0,40]);
%set(gcf, 'Position', [0 0 15 10])
export_fig(sprintf('./figs/%s_sortedpre.pdf',CONDITION), '-transparent')
nextplot(1,1)

set(gcf, 'Position', [0,0, 440,300])

dta = actPpost(1,:);
dta = dta(:);
ad =sort(dta,1,'descend');
bar(ad)
title('Post-training firing rates')
kpost  = trevrolls(ad);
xlabel('Excitatory neuron');
ylabel('Average Firing Rate [Hz]');
ylim([0,40]);
export_fig(sprintf('./figs/%s_sortedpost.pdf',CONDITION), '-transparent')



nextplot(1,3)

kpre  = [];
kpost = [];
for i=1:nruns
    
    kpre(i) = trevrolls(actPpre(i,:));
    kpost(i) = trevrolls(actPpost(i,:));
    %actPpre(i,find(actPpre(i,:)<5)) = 0;
    %actPpost(i,find(actPpost(i,:)<5)) = 0;
    
    kpre(i) = trevrolls(actPpre(i,: ));
    kpost(i) = trevrolls(actPpost(i,:));
end





eb_sp_errs = [stderr(kpre), stderr(kpost)];
eb_sp_mean = [mean(kpre), mean(kpost) ];


%set(gcf, 'Position', [0,0, 440,300])
%barwitherr( eb_sp_errs, eb_sp_mean);

%bar([kpre, kpost]);

%title(sprintf('TrevRolls- %s',CONDITION))
%set(gca, 'XTickLabel', {'Pre', 'Post'})
%savefig(sprintf('./figs/%s_sparsed',CONDITION), 'pdf')




%hold on
%plot(cumsum(hist(actPpost(:), [0:29])), 'r');
%legend('After');
%hold off

%savefig(sprintf('./figs/%s_csum',CONDITION), 'pdf')
% nextplot(1,2)



nextplot
z = []


z(:,1) = mean(actPpre,2);  %  sum(actPpre,2)./sum(actPpre>CUTOFF,2); %mean(actPpre,2)
z(:,2) = mean(actPpost,2); %  sum(actPpost,2)./sum(actPpost>CUTOFF,2);;

eb_firing_errs = std(z);
eb_firing_mean = mean(z)

%barwitherr(eb_firing_errs, eb_firing_mean);
%title('Average Firing Rate')
%ylabel('Avg Firing Rate [Hz]')
%set(gca, 'XTickLabel', {'Pre', 'Test'})

%export_fig(sprintf('./figs/%s_sizes.pdf',CONDITION), '-transparent')


nextplot
 z = [];

 z(:,1) = mean(actPpre>CUTOFF,2);
 z(:,2) = mean(actPpost>CUTOFF,2);
% 
 A1=std(z.*100)
 A2=mean(z.*100)
% 


eb_act_errs = stderr(z.*100);
eb_act_mean = mean(z.*100);

 %barwitherr(eb_act_errs, eb_act_mean);
 %title('Active neurons')
 %ylabel('% of neurons')
 %ylim([0,35]);
%set(gca, 'XTickLabel', {'Pre', 'Test'})

 
 
 
nextplot(1,1)
%hst = brsynsH(:,2:end)
%h=barwitherr(std(hst,0,1)/sqrt(nruns), mean(hst,1))

h = hist_percents(brsyns, npyrbranches)
xlabel('Potentiated synapses per branch')
ylabel('% Branches')
%ylim([0,100])

%set(gcf, 'Position', [0,0, 490,300])
%export_fig(sprintf('./figs/%s_brsyns.pdf', CONDITION), '-transparent')


nextplot(1,1)

h = hist_percents(nrnsyns, npyrs)

ylim([0,15]);
xlabel('Potentiated synapses per neuron')
ylabel('% Excitatory neurons')
set(gcf, 'Position', [0,0, 440,300])
export_fig(sprintf('./figs/%s_nrnsyns.pdf', CONDITION), '-transparent')


vals = sum(brsyns>2,2);
eb_clu_errs = stderr(vals(:));
eb_clu_mean = mean(vals(:));




%nextplot(1,1)
%set(gcf, 'Position', [0,0, 440,300])
%bar(sort(brws(:), 1, 'descend'))
%set(gca, 'XTickLabel', num2str(get(gca, 'XTick')));

%ylabel('Synaptic potentiation');
%xlabel('Dendritic branch');

%export_fig(sprintf('./figs/%s_brwsorted.pdf', CONDITION), '-transparent')



%nextplot(1,1)

we= stderr(synapsesCSUS);
wm= mean(synapsesCSUS );

%set(gcf, 'Position', [0,0, 440,300])
%barwitherr( we, wm, 'r');
%title(sprintf('Size of CS+US clusters  %s',CONDITION))
%set(gca, 'XTickLabel', [1:10])

%savefig(sprintf('./figs/%s_CSUSsyn',CONDITION), 'pdf');
