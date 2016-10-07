conds = {'sparseS2', 'sparseS2G', 'sparseS2L'};

for ncondition=1:length(conds) 
       CONDITION =conds{ncondition}
    

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
       
      ns=1;
      tt= (ns-1)*stimduration;
      nsbefore(run, ns, nn) = length(find(c>tt & c<tt+stimduration));
      
      tt= (ns+1)*stimduration;
      nsafter(run,  ns, nn)  = length(find(c>tt & c<tt+stimduration));
       
      tt= (ns+2)*stimduration;
      nsafterUS(run,  ns, nn)  = length(find(c>tt & c<tt+stimduration));
        
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
    actpostUS = nsafterUS(run, 1, 1:npyrs)/(stimduration/1000);


    actPpre(run, :) = actpre(:);
    actPpost(run, :) = actpost(:);
    actPpostUS(run, :) = actpostUS(:);        
       
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


nextplot(2,2);
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

savefig(sprintf('./figs/%s_rates',CONDITION), 'pdf')


 



bovs = sum(actPpost>CUTOFF & actPpostUS>CUTOFF,2)/npyrs;
bara = 100.0*[mean(sum(actPpost>CUTOFF,2)/npyrs), mean(sum(actPpostUS>CUTOFF,2)/npyrs), mean(bovs)];
erra = 100.0*[stderr(sum(actPpost>CUTOFF,2)/npyrs), stderr(sum(actPpostUS>CUTOFF,2)/npyrs), stderr(bovs)];

nextplot(1,1)
set(gcf, 'Position', [0,0, 440,300])
barwitherr( erra, bara);
ylim([0,50]);
title(('Coding population (%) '))
set(gca, 'XTickLabel', {'S1', 'S2', 'S1 & S2'});
export_fig(sprintf('./figs/%s_S1S2.pdf', CONDITION), '-transparent')

ovs = sum((actPpost>=CUTOFF) & (actPpostUS>=CUTOFF),2)./(sum(actPpost>=CUTOFF | actPpostUS>=CUTOFF,2))

results(sprintf('bara_%s', CONDITION)) = bara(1);
results(sprintf('movs_%s', CONDITION)) = mean(ovs);
results(sprintf('eovs_%s', CONDITION)) = stderr(ovs);


end


pover = @(x) (x*x)/(x+x-(x*x))

if (1)
    nextplot(1,1)
    set(gcf, 'Position', [0,0, 440,300])
    erra = 100.*[ results('eovs_sparseS2G'), results('eovs_sparseS2L'), results('eovs_sparseS2')];
    mova = 100.*[ results('movs_sparseS2G'), results('movs_sparseS2L'), results('movs_sparseS2')];

    barwitherr( erra, mova);
    ylim([0,110]);
    title(('Overlapping population S1&S2 (%) '))
    set(gca, 'XTickLabel', {'Global', 'Local', 'Both'});
    export_fig(sprintf('./figs/%s_OVS1S2.pdf', CONDITION), '-transparent')
end;

