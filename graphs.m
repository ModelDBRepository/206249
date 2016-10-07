clear all
close all
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontSize', 18)



mm = zeros(2, 2)
ee = zeros(2, 2)

mm_clust = zeros(2, 2)
ee_clust = zeros(2, 2)

mm_corr = zeros(2, 1);
ee_corr = zeros(2, 1);

ng = 1

%fig_prefix='2memC_';
%fbase='N100.B40.I10.i6.P2.p1.T%d.S1980.w0c'

for mm_mins=60:60:240
    fn = sprintf(fbase, mm_mins);
    multistats
    mm(ng,:) =m_p';
    ee(ng,:) =s_p';
    mm_clust(ng,:) =m_clust';
    ee_clust(ng,:) =s_clust';
    mm_corr(ng,:) = m_corr(1,2);
    ee_corr(ng,:) = s_corr(1,2);
    ng = ng+1
end

figure()
barwitherr(ee, mm)
ylim([0,40])
set(gca,'XTickLabel',{'1 hour', '2 hours', '3 hours', '4 hours'})

xlabel('Interval')
ylabel('% recruited excitatory neurons')
legend('Memory 1', 'Memory 2')
saveas(gcf, sprintf('./figs/%s_pop.eps', fig_prefix));

figure()
barwitherr(ee_clust*100, mm_clust*100)
ylim([0,100])
set(gca,'XTickLabel',{'1 hour', '2 hours', '3 hours', '4 hours'})

xlabel('Interval')
ylabel('% clustered synapses')
legend('Memory 1', 'Memory 2')
saveas(gcf, sprintf('./figs/%s_clust.eps', fig_prefix));

figure()
barwitherr(ee_corr, mm_corr)
ylim([-0.1,1])
set(gca,'XTickLabel',{'1 hour', '2 hours', '3 hours', '4 hours'})
xlabel('Interval')
ylabel('% correlation')
saveas(gcf, sprintf('./figs/%s_corr.eps', fig_prefix));
