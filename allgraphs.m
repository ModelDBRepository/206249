
figsdir = './figs/'
results = containers.Map()

% strong2 weakstrong sparse repeated multi

close all
%analyze='strong2';

if (strcmp(analyze,'sparse'))
    
    bars = [];
    errs = [];
    bars_sp = [];
    brssum = [];
    local_brssum = [];
    nrnssum = [];
    local_nrnssum = [];
    bars_ff = [];
    bars_act = [];

    CONDITION='sparse';
    %CONDITION='Nsparse';
    an_stats
    
    brssum = brsyns;
    nrnssum = nrnsyns;
    
    both_actppre = actPpre;
    both_actppost = actPpost;
    
    both_kpre = kpre;
    both_kpost = kpost;
    
    d = mean(sumweights);
    bars(1) = d(1);
    bars(4) = d(2);
    d = std(sumweights);
    errs(1) = d(1);
    errs(4) = d(2);
    
    bars_sp(1,:) = [ eb_sp_mean(1), eb_sp_errs(1)]; % Pre
    bars_sp(4,:) = [ eb_sp_mean(2), eb_sp_errs(2)]; % both 
    
    bars_ff(1,:) = [ eb_firing_mean(1), eb_firing_errs(1)]; % Pre
    bars_ff(4,:) = [ eb_firing_mean(2), eb_firing_errs(2)]; % both
    
    bars_act(1,:) = [ eb_act_mean(1), eb_act_errs(1)]; % Pre
    bars_act(4,:) = [ eb_act_mean(2), eb_act_errs(2)]; % both
    
    tr = []
    for i=1:nruns
        tr(1, i) = trevrolls(brws(i,:))
    end
    csus = csusafter;

    CONDITION='sparseL';
    %CONDITION='NsparseL';
    an_stats
    local_csus = csusafter;
    local_brssum = brsyns;
    local_nrnssum = nrnsyns;
    
    local_actppre = actPpre;
    local_actppost = actPpost;
        
    local_kpre = kpre;
    local_kpost = kpost;
    
    d = mean(sumweights);
    bars(3) = d(2);
    d = std(sumweights);
    errs(3) = d(2);
    
    bars_sp(3,:)  = [ eb_sp_mean(2), eb_sp_errs(2)]; % Post Global
    bars_ff(3,:)  = [ eb_firing_mean(2), eb_firing_errs(2)]; % Post Global
    bars_act(3,:) = [ eb_act_mean(2), eb_act_errs(2)]; % Post Global

    for i=1:nruns
        tr(2, i) = trevrolls(brws(i,:))
    end
    

    CONDITION='sparseG';
    %CONDITION='NsparseG';
    an_stats
    global_csus = csusafter;
    global_brssum = brsyns;
    global_nrnssum = nrnsyns;
    global_actppre = actPpre;
    global_actppost = actPpost;
    global_kpre = kpre;
    global_kpost = kpost;
    
    d = mean(sumweights);
    bars(2) = d(2);
    d = std(sumweights);
    errs(2) = d(2);
    
    bars_sp(2,:) = [ eb_sp_mean(2), eb_sp_errs(2)]; % Post Global
    bars_ff(2,:) = [ eb_firing_mean(2), eb_firing_errs(2)]; % Post Global
    bars_act(2,:) = [ eb_act_mean(2), eb_act_errs(2)]; % Post Global

    for i=1:nruns
        tr(2, i) = trevrolls(brws(i,:))
    end



    nextplot(1,1)
    barwitherr(errs, bars);
    %set(h(3), 'facecolor', 'g');
    set(gcf, 'Position', [0,0, 400,300])
    title('Total synaptic weight')
    %ylabel('Total Synaptic Weight')
    %ylim([0,3000]);
    set(gca, 'XTickLabel', {'Pre', 'Somatic', 'Local','S&L' })
    export_fig(sprintf('%s/sparisty_sumweights.pdf', figsdir), '-transparent')


    nextplot(1,1);
    set(gcf, 'Position', [0,0, 400, 300])
    barwitherr( bars_sp(:,2)', bars_sp(:,1)');
    %bar([kpre, kpost]);
    title(sprintf('Population firing sparseness'))
    %ylim([0.4,0.8]);
    set(gca, 'XTickLabel', {'Pre', 'Somatic', 'Local','S&L'})
    export_fig(sprintf('%s/sparse_sp.pdf', figsdir), '-transparent')
    ACT_CUTOFF=5;
    
    nextplot(1,1);
    set(gcf, 'Position', [0,0, 400, 300])
    bars = [];
    %ffs = [];
    for i=1:nruns
        bars(i,1) = mean(both_actppre(i,both_actppre(i,:)>ACT_CUTOFF));
        bars(i,4) = mean(both_actppost(i, both_actppost(i,:)>ACT_CUTOFF));
        bars(i,3) = mean(local_actppost(i, local_actppost(i,:)>ACT_CUTOFF));
        bars(i,2) = mean(global_actppost(i, global_actppost(i,:)>ACT_CUTOFF));
        
    end
    
    cbars = [];
    for i=1:nruns
        cbars(i,1) = mean(both_actppre(i,both_actppre(i,:)>CUTOFF));
        cbars(i,4) = mean(both_actppost(i, both_actppost(i,:)>CUTOFF));
        cbars(i,3) = mean(local_actppost(i, local_actppost(i,:)>CUTOFF));
        cbars(i,2) = mean(global_actppost(i, global_actppost(i,:)>CUTOFF));
    end

    %h = barwitherr( [std(cbars); std(bars)]', [mean(cbars); mean(bars)]');
    h = barwitherr( std(cbars), mean(cbars));
    title(sprintf('Avg Firing Rate (Hz)'))
    %ylim([0.4,0.8]);
    %set(h(2), 'facecolor', [0.7, 0.82,0.9]); 
    %bar([kpre, kpost]);
    %ylim([0,40]);
    %legend('Coding neurons', 'Active neurons')
    set(gca, 'XTickLabel', {'Pre', 'Somatic', 'Local','S&L'})
    export_fig(sprintf('%s/sparse_ff.pdf', figsdir), '-transparent')

    
    ant = bars; % [anpre anpostG anpostL anpostB];
    %%[p, tbl, stats] = anova1(ant)
    %%[res,means] = multcompare(stats,'CType','bonferroni')
    

    nextplot(1,1);
    set(gcf, 'Position', [0,0, 400, 300])
    
    vl = 100*sum((local_actppost>ACT_CUTOFF)')/npyrs;
    vg = 100*sum((global_actppost>ACT_CUTOFF)')/npyrs;
    vb = 100*sum((both_actppost>ACT_CUTOFF)')/npyrs;
    
    vp = 100*sum((both_actppre>ACT_CUTOFF)')/npyrs;
    mbars = [mean(vp), mean(vg), mean(vl), mean(vb)];
    ebars = [std(vp), std(vg), std(vl), std(vb)];
    
    bars = bars_act(1:4,1);
    errs = bars_act(1:4,2);
    h=barwitherr( errs, bars) ; %bars_act(1:4,2)', bars_act(1:4,1)');
    %set(h(2), 'facecolor', [0.7, 0.82,0.9]); 
    %bar([kpre, kpost]);
    title('Coding neurons (%)')
    %ylim([0,40]);
    %legend('Coding neurons', 'Active neurons')
    set(gca, 'XTickLabel', { 'Pre', 'Somatic', 'Local','S&L'})
    export_fig(sprintf('%s/sparse_act.pdf', figsdir), '-transparent')

    
    anpre = 100*sum(global_actppre>10,2)/400;
    anpost = 100*sum(global_actppost>10,2)/400;
    anpostL = 100*sum(local_actppost>10,2)/400;
    anpostB = 100*sum(both_actppost>10,2)/400;
    ant = [anpre anpost anpostL anpostB];
    %%[p, tbl, stats] = anova1(ant);
    %%[res,means] = multcompare(stats,'CType','bonferroni');
    

    
    nextplot(1,1)
    
    set(gcf, 'Position', [0,0, 360, 300])
    av = [];
    
    me = [ std(global_brssum(global_brssum>0)), std(local_brssum(local_brssum>0)), std(brssum(brssum>0))]/10.;
    mm = [ mean(global_brssum(global_brssum>0)), mean(local_brssum(local_brssum>0)), mean(brssum(brssum>0))];
    
   
    barwitherr( me, mm);
    title(('Avg potentiated synapses per branch'))
   %ylim([0,5]);
    ylabel('# Synapses');
    set(gca, 'XTickLabel', { 'Somatic', 'Local', 'S&L'})
    export_fig(sprintf('%s/sparse_syn_per_branch.pdf', figsdir), '-transparent')
    
    
    
    sg = global_brssum(global_brssum>0);
    sl= local_brssum(local_brssum>0);
    sb= brssum(brssum>0);
    
    sg = sg(1:500);
    sl = sl(1:500);
    sb = sb(1:500);
    
    sall = [sg' sl' sb'];
    
    %sall = [sb'  sg'];
    lab = [repmat(1,1,length(sb)) repmat(2,1,length(sg)) repmat(3,1,length(sb)) ];
    %[p, tbl, stats] = kruskalwallis(sall, lab)
    
    boxplot(sall, lab, 'Labels', {'Somatic', 'Local', 'S&L'}, 'colorgroup', lab)
    ylabel('# potentiated synapses per branch')
    %export_fig(sprintf('%s/BOX_BRALLOC.pdf', figsdir), '-transparent')
  
    
    %return;
    
    %%[res,means] = multcompare(stats,'CType','bonferroni')
      
    nextplot(1,1)
    
    set(gcf, 'Position', [0,0, 360, 300])
    av = [];
    
    me = [ std(global_nrnssum(global_nrnssum>0)), std(local_nrnssum(local_nrnssum>0)), std(nrnssum(nrnssum>0))]/10.;
    mm = [ mean(global_nrnssum(global_nrnssum>0)), mean(local_nrnssum(local_nrnssum>0)),  mean(nrnssum(nrnssum>0))];

    
    h=barwitherr( me, mm);
    title(('Avg potentiated synapses per neuron'))
    %ylim([0.2,1.0]);
    ylabel('# Synapses');
    set(gca, 'XTickLabel', { 'Somatic', 'Local', 'S&L'})
    export_fig(sprintf('%s/sparse_syn_per_neuron.pdf', figsdir), '-transparent')
    
      
    sall = [];
    sg = global_nrnssum(global_nrnssum>0);
    sl= local_nrnssum(local_nrnssum>0);
    sb= nrnssum(nrnssum>0);
    %sall = [sg' sl' sb'];
    %lab = [repmat(1,1,length(sg)) repmat(2,1,length(sl)) repmat(3,1,length(sb)) ];
     
    %sg = sg(1:20);
    %sl = sl(1:20);
    %sb = sb(1:20);
    
    sall = [sg' sl' sb'];
    %sall = [sg'  sl'];
    lab = [repmat(1,1,length(sg)) repmat(2,1,length(sl)) repmat(3,1,length(sb)) ];
    [p, tbl, stats] = kruskalwallis(sall, lab)
     
    boxplot(sall, lab, 'Labels', {'Somatic', 'Local', 'S&L'}, 'colorgroup', lab)
    ylabel('# potentiated synapses per neuron')
    export_fig(sprintf('%s/BOX_NRNALLOC.pdf', figsdir), '-transparent')

    
    [p, tbl, stats] = kruskalwallis(sall, lab)
    [res,means] = multcompare(stats,'CType','bonferroni')
    
    nextplot(1,1)
    set(gcf, 'Position', [0,0, 400, 300])
    ms = [  mean(sum(global_nrnssum>0, 2)) ;    mean(sum(local_nrnssum>0, 2)); mean(sum(nrnssum>0, 2))   ];
    ss = [ std(sum(global_nrnssum>0, 2)) ;  std(sum(local_nrnssum>0, 2)); std(sum(nrnssum>0, 2))];
    ms = 100*ms/(npyrs);
    ss = 100*ss/(npyrs);
    barwitherr( ss, ms);
    title(('Neurons with at least 1 synapse'))

    ylabel('% Neurons');
    set(gca, 'XTickLabel', { 'Somatic', 'Local', 'S&L'})
    export_fig(sprintf('%s/sparse_nrn_alloc.pdf', figsdir), '-transparent')
    
    sg = sum(global_nrnssum>0,2);
    sl= sum(local_nrnssum>0,2);
    sb= sum(nrnssum>0,2);
    sall = [sg sl sb];
    %lab = [repmat(1,1,length(sg)) repmat(2,1,length(sl)) repmat(3,1,length(sb)) ];
    %%[p, tbl, stats] = anova1(sall)
    [res,means] = multcompare(stats,'CType','bonferroni')
    
    
    nextplot(1,1)
    set(gcf, 'Position', [0,0, 400, 300])
    ms = [mean(sum(global_nrnssum, 2));  mean(sum(local_nrnssum, 2)); mean(sum(nrnssum, 2))];
    ss = [std(sum(global_nrnssum, 2));  std(sum(local_nrnssum, 2));  std(sum(nrnssum, 2)) ];
    %ms = 100*ms/(npyrs);
    %ss = 100*ss/(npyrs);
    barwitherr( ss, ms);
    title(('Total potentiated synapses'))
    %ylim([0.2,1.0]);
    ylabel('# Potentiated Synapses');
    set(gca, 'XTickLabel', { 'Somatic', 'Local', 'S&L'})
    export_fig(sprintf('%s/sparse_total_syns.pdf', figsdir), '-transparent')
    
    
    
    nextplot(1,1)
    set(gcf, 'Position', [0,0, 400, 300])
    ms = [mean(csus);  mean(local_csus)];
    ss = [stderr(csus);  stderr(local_csus)];
    %ms = 100*ms/(npyrs);
    %ss = 100*ss/(npyrs);
    barwitherr( ss, ms);
    
    title(('CS-US clustering'))
    %ylim([0.2,1.0]);
    ylabel('# CS/US clusters');
    set(gca, 'XTickLabel', { 'S&L', 'Local'})
    export_fig(sprintf('%s/sparse_csus.pdf', figsdir), '-transparent')
end


if (0)
    bspar = [];
    bpop = [];
    bfir = [];
    bclu = [];
    ni=1;

    brov = 'nbrov'
    for i=0.1:0.1:1.0
        CONDITION = sprintf('%s%.1f', brov, i);
        an_stats
        bspar(ni,: ) = [eb_sp_mean(2), eb_sp_errs(2)];
        bfir(ni, :)  = [eb_firing_mean(2), eb_firing_errs(2)];
        bpop(ni, :)  = [eb_act_mean(2), eb_act_errs(2)];
        bclu(ni, :)  = [eb_clu_mean, eb_clu_errs];
	ni = ni+1
    end

    figure
    
	%barwitherr(bspar(:,2), bspar(:,1)); title('Sparseness');
    %set(gca, 'XTick', [1:length(diffs)])
    %set(gca, 'XTickLabel', [0.1:0.1:1.0])
    %export_fig(sprintf('./figs/%s_sp.pdf', brov), '-transparent')

   
    figure
	barwitherr(bfir(:,2), bfir(:,1)); title('Firing rate');
    set(gca, 'XTickLabel', [0.1:0.1:1.0])
    ylim([0,5]);
    export_fig(sprintf('./figs/%s_fir.pdf', brov), '-transparent')
    
    figure
	barwitherr(bpop(:,2), bpop(:,1)); title('Population');
    set(gca, 'XTickLabel', [0.1:0.1:1.0])
    ylim([0,50]);
    export_fig(sprintf('./figs/%s_act.pdf', brov), '-transparent')
    
    figure
	barwitherr(bclu(:,2)/32, bclu(:,1)/32); title('% Branches with > 2 potentiated synapses');
    set(gca, 'XTickLabel', [0.1:0.1:1.0])
    ylim([0,20]);
    export_fig(sprintf('./figs/%s_clu.pdf', brov), '-transparent')

end


if (strcmp(analyze,'weakstrong'))
    
    results = containers.Map();
CONDAR={'weakstrongN','weakstrong'}
for kk=1:2
     COND= CONDAR{kk}  %'weakstrongN';

     CONDITION=[COND 'L'];
     weastrong
     local_brcommon = brcommon;
     local_brtsyns = brtsyns;
     local_nrntsyns = nrntsyns;

 
     CONDITION=[COND 'G'];
     weastrong
     global_brcommon = brcommon;
     global_brtsyns = brtsyns;
     global_nrntsyns = nrntsyns;
     
     CONDITION=COND; %'weakstrong';
     weastrong

     
   
    close all
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b1 = results([COND 'act'])
    %errorbar( b(2,:,2), b(1,:,2));
    %hold on
    b2 = results([COND 'Lact'])
    %errorbar( b(2,:,2), b(1,:,2), 'g');
    %hold off
    
    b22 = results([COND 'Gact'])
    
    b3 = [ b1(1,:,2); b2(1,:,2) ; b22(1,:,2)];
    b4 = [ b1(2,:,2); b2(2,:,2) ; b22(2,:,2)];
    
    h = barwitherr(b3', b4', 'LineStyle', 'none');
    set(h(3), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0.0,0.5,0]);
    
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    rotateXLabels(gca, 60)
    
    ylim([0,30])
    title(sprintf('Coding Neurons (%%)', CONDITION));
    
    xlim([0,9])
    export_fig(sprintf('./figs/WS_act_%s.pdf', COND), '-transparent')

    
     c3 = results([COND 'act_d']);
     c2 = results([COND 'Lact_d']);
     c1 = results([COND 'Gact_d']);
     c3 = c3(:,:,2);
     c2 = c2(:,:,2);
     c1 = c1(:,:,2);
     [p, tbl, stats] = anova1(c1)
     [res,means] = multcompare(stats, 'CType', 'bonferroni')
%     
   
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b1 = results([COND 'act'])
    b2 = results([COND 'Lact'])
    b22 = results([COND 'Gact'])  
    
    b3 = [ b22(1,:,1);  b2(1,:,1); b1(1,:,1) ];
    b4 = [ b22(2,:,1);  b2(2,:,1); b1(2,:,1)  ];
    h = barwitherr(b3', b4', 'LineStyle', 'none');
      set(h(3), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0,0.5,0.0]);
    
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    rotateXLabels(gca, 60)

    ylim([0,40])
    ylabel('Coding Neurons (%)'); 
    title('Strong memory')
    xlim([0,9])
    export_fig(sprintf('./figs/WS_act_strong_%s.pdf', COND), '-transparent')
        
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b1 = results([COND 'ff'])
    b2 = results([COND 'Lff'])
    b22 = results([COND 'Gff'])

    b3 = [ b22(1,:,2); b2(1,:,2) ; b1(1,:,2) ];
    b4 = [ b22(2,:,2) ; b2(2,:,2); ; b1(2,:,2)];
    
    h = barwitherr(b3', b4', 'LineStyle', 'none');
    set(h(3), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0.0,0.5,0]);
    
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    rotateXLabels(gca, 60)

    title('Avg Firing Frequency (Hz)'); 
    xlim([0,9])
    export_fig(sprintf('./figs/WS_ff_%s.pdf', COND), '-transparent')


     c3 = results([COND 'act_d']);
     c2 = results([COND 'Lact_d']);
     c1 = results([COND 'Gact_d']);
     c3 = c3(:,:,2);
     c2 = c2(:,:,2);
     c1 = c1(:,:,2);
     [p, tbl, stats] = anova1(c2)
     [res,means] = multcompare(stats, 'CType', 'bonferroni')
    
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b1 = results([COND 'firingcor'])
    b2 = results([COND 'Lfiringcor'])
    b22 = results([COND 'Gfiringcor'])
        
    b3 = [ b22(1,:);  b2(1,:); b1(1,:) ];
    b4 = [ b22(2,:);  b2(2,:); b1(2,:) ];
    
    h = barwitherr(b3', b4', 'LineStyle', 'none');
    set(h(3), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0,0.5,0.0]);
    
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    rotateXLabels(gca, 60)
    ylim([0,100]);

    ylabel('Correlation'); 
    title('Firing rate vector correlation')
    xlim([0,9])
    export_fig(sprintf('./figs/WS_firingcor_%s.pdf', COND), '-transparent')
    
    
    
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b1 = results([COND 'common'])
    b2 = results([COND 'Lcommon'])
    b22 = results([COND 'Gcommon'])
        
    b3 = [ b22(1,:);  b2(1,:); b1(1,:) ];
    b4 = [ b22(2,:);  b2(2,:); b1(2,:) ];
    
    h = barwitherr(b3', b4', 'LineStyle', 'none');
    set(h(3), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0,0.5,0.0]);
    
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    rotateXLabels(gca, 60)
    ylim([0,100]);

    ylabel('Common neurons (%)'); 
    title('Common recruited neurons ')
    xlim([0,9])
    export_fig(sprintf('./figs/WS_coract_%s.pdf', COND), '-transparent')

 
    
    
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b1 = results([COND 'brcors'])
    b2 = results([COND 'Lbrcors'])
    b22 = results([COND 'Gbrcors'])
     
    b3 = [b22(1,:); b1(1,:); b2(1,:) ];
    b4 = [b22(2,:); b1(2,:); b2(2,:) ];
    
    h = barwitherr(b3', b4', 'LineStyle', 'none');
    set(h(3), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0,0.5,0.0]);
    
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    rotateXLabels(gca, 60);
    ylabel('Similarity');
    title('Similarity of synaptic projection patterns per branch')
    ylim([0,.6]);
    xlim([0,9])
    export_fig(sprintf('./figs/WS_brcor_%s.pdf', COND), '-transparent')

    
    
    
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b1 = results([COND 'nrncors'])
    b2 = results([COND 'Lnrncors'])
    b22 = results([COND 'Gnrncors'])
     
    b3 = [b22(1,:); b1(1,:); b2(1,:) ];
    b4 = [b22(2,:); b1(2,:); b2(2,:) ];
    
    h = barwitherr(b3', b4', 'LineStyle', 'none');
    set(h(3), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0,0.5,0.0]);
    
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    rotateXLabels(gca, 60);
    ylabel('Similarity');
    title('Similarity of synaptic projection patterns per neuron')
    %ylim([=,.8]);
    xlim([0,9])
    export_fig(sprintf('./figs/WS_nrncor_%s.pdf', COND), '-transparent')

    
    
    
    
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b3 = 100.*[ mean(global_brcommon,2), mean(local_brcommon,2), mean(brcommon,2)];
    b4 = 100.*[  std(global_brcommon,0,2), std(local_brcommon,0,2), std(brcommon,0,2)];
    
    h = barwitherr(b4, b3, 'LineStyle', 'none');
    set(h(3), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0.0,0.5,0]);
    
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    
    rotateXLabels(gca, 60)
    ylim([0,50]);
    title('% Branches with clusters of both memories'); 
    ylabel('% Branches')
    xlim([0,9])
    export_fig(sprintf('./figs/WS_brcommon_%s.pdf', COND), '-transparent')
    
    
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b3 = [ mean(brtsyns,2), mean(local_brtsyns,2)];
    b4 = [ std(brtsyns,0,2), std(local_brtsyns,0,2)];
    h = barwitherr(b4, b3, 'LineStyle', 'none');
    set(h(1), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', diffs)
    rotateXLabels(gca, 60)
    ylim([0,900]);
    xlim([0,13])
    title('Number of branches containing the weak memory'); 
    ylabel('# branches')
    export_fig(sprintf('./figs/WS_brtsyns_%s.pdf', COND), '-transparent')
    
    
    figure;
    set(gcf, 'Position', [0,0, 440,300])
    b3 = [ mean(nrntsyns,2), mean(local_nrntsyns,2)];
    b4 = [ std(nrntsyns,0,2), std(local_nrntsyns,0,2)];
    h = barwitherr(b4, b3, 'LineStyle', 'none');
    set(h(1), 'facecolor', [0,0,0.5]);
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', diffs)
    rotateXLabels(gca, 60)
    ylim([0,600]);
    xlim([0,13])
    title('Number of neurons containing the weak memory'); 
    ylabel('# neurons')
    export_fig(sprintf('./figs/WS_nrntsyns_%s.pdf', COND), '-transparent')
    

    
     [p, tbl, stats] = anova1(results([COND 'common_d']))
     [res,means] = multcompare(stats, 'CType', 'bonferroni')
    
     [p, tbl, stats] = anova1(results([COND 'firingcor_d']))
     [res,means] = multcompare(stats, 'CType', 'bonferroni')
     
     [p, tbl, stats] = anova1(results([COND 'brcors_d']))
     [res,means] = multcompare(stats, 'CType', 'bonferroni')
     
     [p, tbl, stats] = anova1(results([COND 'nrncors_d']))
     [res,means] = multcompare(stats, 'CType', 'bonferroni')
     
     
    end
end



if (strcmp(analyze,'strong2'))
    close all

    CONDITION='strong2L'
    strong2
    local_brcommon = brcommon;
    local_totactive= totactive;
    local_totcommon = totcommon;

    CONDITION='strong2G'
    strong2
    global_brcommon = brcommon;
    global_totactive= totactive;
    global_totcommon = totcommon;

    CONDITION='strong2'
    strong2

    
    close all
    
    for NMEM=1:2
    
        figure;
        hold on
        mact = 100.0*mean(totactive,1)
        sact = 100.0*std(totactive,0,1)/(sqrt(nruns))
        %errorbar(mact(:,:,1), sact(:,:,1), 'b.');
        errorbar(mact(:,:,NMEM), sact(:,:,NMEM), 'b');

        mact = 100.0*mean(local_totactive,1)
        sact = 100.0*std(local_totactive,0,1)/(sqrt(nruns))
        %errorbar(mact(:,:,1), sact(:,:,1), COL);
        errorbar(mact(:,:,NMEM), sact(:,:,NMEM), 'r');

        mact = 100.0*mean(global_totactive,1)
        sact = 100.0*std(global_totactive,0,1)/(sqrt(nruns))
        %errorbar(mact(:,:,1), sact(:,:,1), COL);
        errorbar(mact(:,:,NMEM), sact(:,:,NMEM), 'g');
        hold off

        set(gca, 'XTick', [1:length(diffs)])
        set(gca, 'XTickLabel', difflabels)
        title('% Active neurons');
        %xlabel('Weak-Strong Interval [minutes]')
        ylim([10,50])
        export_fig(sprintf('./figs/%s_ACT%d.pdf',CONDITION, NMEM), '-transparent')
    end

    

    for NMEM=1:2
        figure
        hold on
        tt = totactive(:,:,NMEM)
        dt = bsxfun(@rdivide, bsxfun(@minus, tt, tt(:,4)), tt(:,4))
        errorbar(100*mean(dt), NMEM*100*stderr(dt), 'b')

        tt = local_totactive(:,:,NMEM)
        dt = bsxfun(@rdivide, bsxfun(@minus, tt, tt(:,4)), tt(:,4))
        errorbar(100*mean(dt), NMEM*100*stderr(dt), 'r')

        tt = global_totactive(:,:,NMEM)
        dt = bsxfun(@rdivide, bsxfun(@minus, tt, tt(:,4)), tt(:,4))
        errorbar(100*mean(dt), NMEM*100*stderr(dt), 'g')
        hold off

        set(gca, 'XTick', [1:length(diffs)])
        set(gca, 'XTickLabel', difflabels)
        title('Increase in Coding Neurons (%)');
        %xlabel('Weak-Strong Interval [minutes]')
        ylim([0,100])
        export_fig(sprintf('./figs/%s_increase%d.pdf',CONDITION, NMEM), '-transparent')
    end

    figure
    bm = [100.*std(global_totcommon)/(sqrt(nruns)); 100.*std(local_totcommon)/(sqrt(nruns)); 100.*std(totcommon)/(sqrt(nruns)) ];
    be = [100.*mean(global_totcommon); 100.*mean(local_totcommon) ; 100.*mean(totcommon)];

    h=barwitherr(bm', be');
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0,0.5,0.0]);    
    set(h(3), 'facecolor', [0,0.0,0.5]);
    
    title(sprintf('Common Coding Neurons (%%) %s', CONDITION))
    ylabel('% common neurons')
    %xlabel('Weak-Strong Interval [minutes]')
    ylim([0,50])
    export_fig(sprintf('./figs/%s_common.pdf',CONDITION), '-transparent')

    
    figure
    bm = [100.*std(global_brcommon)/(sqrt(nruns)); 100.*std(local_brcommon)/(sqrt(nruns));; 100.*std(brcommon)/(sqrt(nruns)) ];
    be = [100.*mean(global_brcommon);  100.*mean(local_brcommon); 100.*mean(brcommon)];
     
    h=barwitherr(bm', be');
    set(h(2), 'facecolor', [0.5,0.0,0]);
    set(h(1), 'facecolor', [0,0.5,0.0]);    
    set(h(3), 'facecolor', [0,0.0,0.5]);    
    
    set(gca, 'XTick', [1:length(diffs)])
    set(gca, 'XTickLabel', difflabels)
    title('% Branches with clusters of both memories')
    ylabel('% branches')
    %xlabel('Weak-Strong Interval [minutes]')
    ylim([0,50])
    export_fig(sprintf('./figs/%s_brcommon.pdf',CONDITION), '-transparent')

end



if (strcmp(analyze,'dir'))
    close all
    
    DIRS = {'dir1', 'dir2'};
    for ddir =1:length(DIRS)
        dir = DIRS{ddir};
        
        CONDITION=[ dir 'L']
        dir2
        local_brcommon = brcommon;
        local_totactive= totactive;
        local_totcommon = totcommon;
        
        CONDITION= [dir 'G']
        dir2
        global_brcommon = brcommon;
        global_totactive= totactive;
        global_totcommon = totcommon;
        
        CONDITION= dir
        dir2
        
        close all
        
        for NMEM=1:2
            
            figure;
            hold on
            mact = 100.0*mean(totactive,1)
            sact = 100.0*std(totactive,0,1)/(sqrt(nruns))
            %errorbar(mact(:,:,1), sact(:,:,1), 'b.');
            errorbar(mact(:,:,NMEM), sact(:,:,NMEM), 'b');
            
            mact = 100.0*mean(local_totactive,1)
            sact = 100.0*std(local_totactive,0,1)/(sqrt(nruns))
            %errorbar(mact(:,:,1), sact(:,:,1), COL);
            errorbar(mact(:,:,NMEM), sact(:,:,NMEM), 'r');
            
            mact = 100.0*mean(global_totactive,1)
            sact = 100.0*std(global_totactive,0,1)/(sqrt(nruns))
            %errorbar(mact(:,:,1), sact(:,:,1), COL);
            errorbar(mact(:,:,NMEM), sact(:,:,NMEM), 'g');
            hold off
            
            set(gca, 'XTick', [1:length(diffs)])
            set(gca, 'XTickLabel', difflabels)
            title('% Active neurons');
            %xlabel('Weak-Strong Interval [minutes]')
            ylim([0,70])
            export_fig(sprintf('./figs/%s_ACT%d.pdf',CONDITION, NMEM), '-transparent')
        end
        
        for NMEM=1:2
            figure
            hold on
            tt = totactive(:,:,NMEM)
            dt = bsxfun(@rdivide, bsxfun(@minus, tt, tt(:,4)), tt(:,4))
            errorbar(100*mean(dt), 100*stderr(dt), 'b')
            
            tt = local_totactive(:,:,NMEM)
            dt = bsxfun(@rdivide, bsxfun(@minus, tt, tt(:,4)), tt(:,4))
            errorbar(100*mean(dt), 100*stderr(dt), 'r')
            
            tt = global_totactive(:,:,NMEM)
            dt = bsxfun(@rdivide, bsxfun(@minus, tt, tt(:,4)), tt(:,4))
            errorbar(100*mean(dt), 100*stderr(dt), 'g')
            hold off
            
            set(gca, 'XTick', [1:length(diffs)])
            set(gca, 'XTickLabel', difflabels)
            title('Increase in Coding Neurons (%)');
            %xlabel('Weak-Strong Interval [minutes]')
            ylim([0,100])
            export_fig(sprintf('./figs/%s_increase%d.pdf',CONDITION, NMEM), '-transparent')
        end
        
        figure
        bm = [100.*std(global_totcommon)/(sqrt(nruns)); 100.*std(local_totcommon)/(sqrt(nruns)); 100.*std(totcommon)/(sqrt(nruns)) ];
        be = [100.*mean(global_totcommon); 100.*mean(local_totcommon) ; 100.*mean(totcommon)];
        
        h=barwitherr(bm', be');
        set(gca, 'XTick', [1:length(diffs)])
        set(gca, 'XTickLabel', difflabels)
        set(h(2), 'facecolor', [0.5,0.0,0]);
        set(h(1), 'facecolor', [0,0.5,0.0]);
        
        title(sprintf('Common Coding Neurons (%%) %s', CONDITION))
        ylabel('% common neurons')
        %xlabel('Weak-Strong Interval [minutes]')
        ylim([0,80])
        export_fig(sprintf('./figs/%s_common.pdf',CONDITION), '-transparent')
        
        
        figure
        bm = [100.*std(global_brcommon)/(sqrt(nruns)); 100.*std(local_brcommon)/(sqrt(nruns));; 100.*std(brcommon)/(sqrt(nruns)) ];
        be = [100.*mean(global_brcommon);  100.*mean(local_brcommon); 100.*mean(brcommon)];
        
        h=barwitherr(bm', be');
        set(h(2), 'facecolor', [0.5,0.0,0]);
        set(h(1), 'facecolor', [0,0.5,0.0]);
        
        set(gca, 'XTick', [1:length(diffs)])
        set(gca, 'XTickLabel', difflabels)
        title('% Branches with clusters of both memories')
        ylabel('% branches')
        %xlabel('Weak-Strong Interval [minutes]')
        ylim([0,50])
        export_fig(sprintf('./figs/%s_brcommon.pdf',CONDITION), '-transparent')
    end
end





if (strcmp(analyze, 'three'))   
    %CONDITION='threewww'
    %three
    CONDITION='threesws'
    three
    CONDITION='threewsw'
    three
    
end


if (strcmp(analyze, 'repeated'))
    %CONDITION='repeated'
    %repeated
    %global_brsynratio = brsynratio;
    
    CONDITION='repeatedUL'
    repeated
    local_brsynratio = brsynratio;
    
    
    ta = local_brsynratio';
    t = table([1 2 3 4 5 6 7 8 9 10]' , ta(:,1), ta(:,2), ta(:,3),ta(:,4));
    
    CONDITION='repeatedUG'
    repeated
    global_brsynratio = brsynratio;
    
    
end


if (strcmp(analyze, 'multi'))
    CONDITION='multiG'
    COL='g'
    multistats
    global_mm = mm;
    i=3
    cors = sum(m .* tril(circshift(eye(npatterns), i)));
    global_cors = cors(1:npatterns-i);
    
    
    COL='r'
    CONDITION='multiL'
    multistats
    
    local_mm = mm;
    i=3
    cors = sum(m .* tril(circshift(eye(npatterns), i)));
    local_cors = cors(1:npatterns-i);
    
    aa = [global_cors' local_cors']
    anova1(aa)
    
    COL='b'
    CONDITION='multiGN'
    multistats
    globalN_mm = mm;
        
    COL='r'
    CONDITION='multiLN'
    localN_mm = mm;
    multistats
end


if (strcmp(analyze,'brtest'))
    
    CONDITION='brtestL'
    an_brtest
    local_bb = bb;
    local_ba = ba;
    local_brcase = brcase;
    local_brcases = brcases;
    local_actPpre = actPpre;
    local_actPost = actPpost;
    
    CONDITION='brtestG'
    an_brtest
    global_bb = bb;
    global_ba = ba;
    global_brcase = brcase;
    global_brcases = brcases;
    global_actPpre = actPpre;
    global_actPost = actPpost;
    
    CONDITION='brtest'
    an_brtest
   
    close all
    
    
    figure
    %errorbar(100.*mean(ba,2),100.*std(ba,0,2), 'c')
    %hold on
    errorbar(100.*mean(bb,2),100.*std(bb,0,2), 'b')
    hold on
    errorbar(100.*mean(local_bb,2),100.*std(local_bb,0,2), 'r')
    errorbar(100.*mean(global_bb,2),100.*std(global_bb,0,2), 'g')

    set(gca, 'XTick', [1:4])
    set(gca, 'XTickLabel', [10:10:40])
    title(sprintf('Active population'))
    legend( 'S&L', 'Local', 'Somatic');
    ylabel('% Active neurons')
    xlabel('Number of branches per neuron')
    hold off
    export_fig(sprintf('./figs/brtest_act.pdf'), '-transparent')
  
    
    figure
    
    
    errorbar(brcase,brcases, 'b')
    hold on
    errorbar(local_brcase,local_brcases, 'r')
    errorbar(global_brcase,global_brcases, 'g')
    set(gca, 'XTick', [1:4])
    set(gca, 'XTickLabel', [10:10:40])
    xlabel('Number of branches per neuron')
    ylabel('Potentiated synapses per branch');
    title(sprintf('Avg. Potentiated synapses per branch'))
    legend( 'S&L', 'Local', 'Somatic');
    export_fig(sprintf('./figs/brtest_clu.pdf'), '-transparent')

    
    
    figure
    
    errorbar(mean(mean(actPpre,3),2),std(mean(actPpre,3),0,2), 'b')
    hold on
    
    errorbar(mean(mean(local_actPpre,3),2),std(mean(local_actPpre,3),0,2), 'r')
    
    set(gca, 'XTick', [1:4])
    set(gca, 'XTickLabel', [10:10:40])
    xlabel('Number of branches per neuron')
    ylabel('Average firing rate [Hz]');
    title(sprintf('Average firing rate'))
    legend('Somatic', 'Local');
    export_fig(sprintf('./figs/brtest_ff.pdf'), '-transparent')
        
end


