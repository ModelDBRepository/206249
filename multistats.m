defaults

close all



npatterns=10
nruns =10


%CONDITION='multiL';
%ISI=60


spks = zeros(npatterns, npyrs, nruns);
pops = zeros(npatterns, npyrs, nruns);
corr_spk = zeros(npatterns, npatterns, nruns);
corr_pop = zeros(npatterns, npatterns, nruns);

branch_syns = zeros(ninputs, npyrs*nbranches, nruns);
br_hists = zeros(ninputs, 12, nruns);
clustering = zeros(ninputs, nruns);
brstrengths = zeros(ninputs, npyrs*nbranches);

brweights = zeros(ninputs, npyrs*nbranches, nruns);
nrnweights = zeros(ninputs, npyrs, nruns);
brweightcors = zeros(ninputs, ninputs, nruns);
brsyncors= zeros(ninputs, ninputs, nruns);
nrnweightcors = zeros(ninputs, ninputs, nruns);

brcommon = zeros(ninputs, ninputs, nruns);


clust_all = {};
clust_all = cell(9,1);

for run = 1:nruns
    sprintf('./data/%s_%d_%d/spikesperpattern.dat', CONDITION, ISI,run-1)
    spk = load( sprintf('./data/%s_%d_%d/spikesperpattern.dat', CONDITION, ISI,run-1));
    
    recallspikes = spk(:, 1:npyrs)/(stimduration/1000); 
% 
%     figure()
%     imagesc(recallspikes');
%     rsc = recallspikes'
%     rsc = diag(1./sum(rsc,2))*rsc
%     rsc(isnan(rsc)) =0
%     
%     return;
%hist(recallspikes(:),20);
    
    pop = recallspikes>CUTOFF; %Hz
    spks(:, :, run) = recallspikes;
    pops(:, :, run) = pop;
    
    corr_spk(:,:, run) = corrcoef(recallspikes');
    %corr_pop(:,:, run) = corrcoef(pop');
    
    %corr_pop(:,:, run) = corrcoef(pop');
    
    for nk=1:10
        for nl=1:10
            corr_pop(nk,nl,run) = sum(pop(nk,:)&pop(nl,:)) / ((sum(pop(nk, :))+sum(pop(nl,:)) )/2);
        end
    end
    
    ff = sprintf('./data/%s_%d_%d/synstate.dat', CONDITION, ISI,run-1);   
    ss = load(ff);
    
    for i=1:size(ss,1)
        bid=ss(i,2);
        nid=ss(i,3);
        srcid=ss(i,5);
        bstrength = ss(i,6);
        w=ss(i,7);
        if (srcid >= 0 && bid <= npyrs*nbranches)
            brweights(srcid+1, bid+1, run) = brweights(srcid+1, bid+1, run) + w;
            brstrengths(srcid+1, bid+1)=bstrength;
            nrnweights(srcid+1, nid+1,run) = nrnweights(srcid+1, nid+1,run) + w;
        end
        if (srcid >= 0 && bid <= npyrs*nbranches &&  w > 0.7)
            branch_syns(srcid+1, bid+1, run) = branch_syns(srcid+1, bid+1, run)+1;
        end
    end

    for i=1:npatterns
        bs = branch_syns(i, :, run);
        b = bs(bs>0);
        x = 1:12;
        [d, h] = hist(b, x);
        br_hists(i, :, run) = d;
        ss = sum(d(1:end));
        if (ss>0)
            clustering(i,run) = sum(d(2:end))/ss;
        end
    end
    brweightcors(:, :, run) = corrcoef(brweights(:,:, run)'); 
    brsyncors(:, :, run) = corrcoef(branch_syns(:,:, run)'); 
    nrnweightcors(:, :, run) = corrcoef(nrnweights(:,:, run)'); 
    
    brclust = branch_syns(:,:,run)>1;
    for nk=1:10
        for nl=1:10
            brcommon(nk,nl, run) = sum(brclust(nk,:)&brclust(nl,:))/(sum(brclust(nk,:)|(brclust(nl,:))));
        end
    end
    
      
    brclust = branch_syns(:,:,run)>0;
    for nk=1:10
        for nl=1:nk-1
            brtot = branch_syns(nk, :, run) + branch_syns(nl, :, run);
            brl = brtot(find(brclust(nk,:)&brclust(nl,:)));
            clust_all{nk-nl} = [clust_all{nk-nl} brl];
            %brcommon(nk,nl, run) = sum(brclust(nk,:)&brclust(nl,:))/(sum(brclust(nk,:)|(brclust(nl,:))));
        end
    end
end



% figure();
% imagesc(pops(:,:,1)');
% colorbar()
% %title('Firing rates per memory recall (Hz)');
% xlabel('Memory #');
% ylabel('Pyramidal Neuron #');


m = mean(corr_spk, 3);
m_corr = m;
s_corr = std(corr_spk, 0,3);
figure();


imagesc(m, [-0.2, 1.0]);
%colorbar()
title('Similarity between population firing patterns');
axis off
%xlabel('Memory #')
%ylabel('Memory #')


export_fig(sprintf('./figs/%s_ffsim.pdf',CONDITION), '-transparent')



m = mean(corr_pop, 3)';
figure();
imagesc(m, [0, 1.0]);
%colorbar();
axis off
title('Population overlap between memories');
xlabel('Memory #')
ylabel('Memory #')
%colorbar()

export_fig(sprintf('./figs/%s_popoverlap.pdf',CONDITION), '-transparent')


b = zeros(npatterns, 2);
for i=1:npatterns
    cors =  sum(m .* tril(circshift(eye(npatterns), i)));
    cors = cors(1:npatterns-i);
    b(i,1) = mean(cors);
    b(i,2) = stderr(cors);
end

figure()
bars_ovl = 100.0*b(1:end-2,1);
bars_err = 100.0*b(1:end-2,2);
barwitherr(100.0*b(1:end-2,2), 100.0*b(1:end-2,1), COL)
xlabel('Hours between memories')
ylabel('% Overlapping population')
title('Overlap between populations')
ylim([0,100])
export_fig(sprintf('./figs/%s_popsimbar.pdf',CONDITION), '-transparent')



figure()
m = mean(nrnweightcors, 3)';
imagesc(m, [-0.2, 1.0]);
title('Similarity of synaptic projection patterns per neuron')
axis off
%colorbar()
export_fig(sprintf('./figs/%s_nrnw.pdf',CONDITION), '-transparent')


b = zeros(npatterns, 2);
for i=1:npatterns
    cors =  sum(m .* tril(circshift(eye(npatterns), i)));
    cors = cors(1:npatterns-i);
    b(i,1) = mean(cors);
    b(i,2) = stderr(cors);
end

figure()
barwitherr(b(1:end-2,2), b(1:end-2,1), COL)
xlabel('Hours between memories')
ylabel('Average similarity')
title('Similarity of synaptic projection patterns per neuron')

%ylim([,0.8])
export_fig(sprintf('./figs/%s_nrnwbar.pdf',CONDITION), '-transparent')



figure()
m = mean(brweightcors, 3)';
imagesc(m, [-0.2, 1.0]);
title('Similarity of synaptic projection patterns per branch')
axis off
%colorbar()

export_fig(sprintf('./figs/%s_brw.pdf',CONDITION), '-transparent')


b = zeros(npatterns, 2);
for i=1:npatterns
    cors =  sum(m .* tril(circshift(eye(npatterns), i)));
    cors = cors(1:npatterns-i);
    b(i,1) = mean(cors);
    b(i,2) = stderr(cors);
end

figure()
barwitherr(b(1:end-2,2), b(1:end-2,1), COL)
xlabel('Hours between memories')
ylabel('Average similarity')
title('Similarity of synaptic projection patterns per branch')
%ylim([0,0.4])
export_fig(sprintf('./figs/%s_brwbar.pdf',CONDITION), '-transparent')





figure()
tp = sum(pops, 2)*100.0/npyrs;
m_p = mean(tp, 3)
s_p = std(tp, 0, 3)
bar(m_p);
hold on
h=errorbar(m_p', s_p')
set(h(1), 'color', 'red');set(h(1), 'LineStyle', 'None');
hold off;
title('Active population per memory')
%ylabel('Active Pyr. Neurons (%)')
%xlabel('Memory #')





figure()
tp = sum(spks, 2)/npyrs;
m_p = mean(tp, 3)
s_p = std(tp, 0, 3)
bar(m_p);
hold on
h=errorbar(m_p', s_p')
set(h(1), 'color', 'red');set(h(1), 'LineStyle', 'None');
hold off;
title('Avg Firing rate of pyramidal neurons')
%ylabel('Firing Rate [Hz]')
%xlabel('Memory #')

if (0)
    figure()
    nmem =2
    sb = mean(br_hists, 3);
    sbd = std(br_hists,0,3);
    bar(sb(nmem, :))
    hold on
    h = errorbar(sb(nmem, :), sbd(nmem,:));

    set(h(1), 'LineStyle', 'None');
    title('Distribution of potentiated synapses per branch')
    xlabel('Number of potentiated synapses')
    ylabel('Number of branches')
    yl = ylim(); yl(1) = 0; ylim(yl);
    
    %saveas(gcf,'./figs/norep4.eps', 'epsc');
    %imagesc(corrcoef(mp'));

    figure()
    m_clust = mean(clustering, 2)
    s_clust = std(clustering, 0, 2)
    bar(m_clust)
    hold on
    h = errorbar(m_clust, s_clust);
    set(h(1), 'LineStyle', 'None');
    title('Clustered synapses per memory')
    xlabel('Memory number')
    ylabel('Percentage of clustered synapses')



    % figure()
    % m = mean(brweightcors, 3)';
    % imagesc(m, [-0.2, 1.0]);
    % title('Correlation between branches')
    % colorbar()
    % %imagesc(corrcoef(mp'));

    figure()
    m = mean(brsyncors, 3)';
    imagesc(m, [-0.2, 1.0]);
    %title('Correlation between branches')
    %colorbar()

end




figure()
mm=sum(branch_syns(:,1:npyrs*nbranches,:)>0, 1);
mean(mm(:))
std(mm(:))

xbins = [0:10];
mh = zeros(nruns,size(xbins,2))
for i=1:nruns
    [d,h] = hist(mm(:,:,i), xbins);
    mh(i, :) = d;
end

mh = mh/(npyrs*nbranches);
barwitherr(std(mh,0,1), mean(mh, 1), COL)


title(sprintf('Memories represented per branch'))
ylabel('Probability')
xlabel( 'Number of memories represented');
set(gca,'Xtick', [0:11], 'XTickLabel',[0 0:11]);

export_fig(sprintf('./figs/%s_brr.pdf',CONDITION), '-transparent')





figure()
m = (mean(brcommon, 3)') ;
imagesc(m, [0, 1.0]);
title('% Branches with clusters of both memories')
axis off
%colorbar()
export_fig(sprintf('./figs/%s_brcommon.pdf',CONDITION), '-transparent')


b = zeros(npatterns, 2);
for i=1:npatterns
    cors =  sum(m .* tril(circshift(eye(npatterns), i)));
    cors = cors(1:npatterns-i);
    b(i,1) = mean(cors);
    b(i,2) = stderr(cors);
end

figure()
barwitherr(100.*b(1:end-2,2), 100.*b(1:end-2,1), COL)
xlabel('Hours between memories')
ylabel('% branches with clusters')
title('% Branches with clusters of both memories')

%ylim([,0.8])
export_fig(sprintf('./figs/%s_brcommon_bar.pdf',CONDITION), '-transparent')

brovl_err = 100.*b(1:end-2,2);
brovl_mean = 100.*b(1:end-2,1);



figure;
cc=winter(12);
nn = [1,2,3,4,5,6,7,8,9];
color = [0,0,1];

for i=1:length(nn)
    ncase = nn(i);
    
    %diff = diffs(ncase);
    
    
    %aa  = mean(histCSUS(:,ncase,:));
    %bb = stderr(histCSUS(:,ncase,:));
    %aa = aa(:);
    
    %bb =bb(:);
    aa = histc(clust_all{i}, [1:20])/(10*(10-i));
    plot(aa, 'Color',  cc(ncase,:));
    hold on
end
ylim([0,2000]);
xlim([0,15]);
ylabel('Number of clusters of both memories');
xlabel('Synapses per  cluster');

legend({'1 hour', '2 hours', '3 hours', '4', '5','6','7','8', '9'});
export_fig(sprintf('./figs/%s_clustering2.pdf',CONDITION), '-transparent')
hold off;

