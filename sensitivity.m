defaults
close all
nruns=10;
ndiffs = 7;


DATADIR='./data';

totcommon = zeros(nruns, ndiffs, 3);
totfiring = zeros(nruns, ndiffs, 3, 2);
totactive = zeros(nruns, ndiffs, 3, 2);
totactive2 = zeros(nruns, ndiffs, 3, 2);

overlap = zeros(nruns, ndiffs, 3);
overlap2 = zeros(nruns, ndiffs, 3);
coract = zeros(nruns, ndiffs);

brsyns = zeros(nruns,npyrs*nbranches);
brsynsH = zeros(nruns,10);

nrnsyns = zeros(nruns,npyrs);

nrnbins = [0:40];
nrnsynsH = zeros(nruns, length(nrnbins));

brcors = zeros( nruns, ndiffs);
nrncors = zeros(nruns, ndiffs);


brcommon = [];

CONDITION='GPL';

ncase=1;

conds = {'connectivityParam' ,'homeostasisTimeParam' ,'CREBTimeParam' ,'inhibitionParam' ,'globalPRPThresh' ,'localPRPThresh' ,'dendSpikeThresh' ,'initWeight' ,'maxWeight', 'nBranchesParam', 'nNeuronsParam'};

condLabels = {'Afferent connections'  ,'Homeostasis time constant' ,'CREB time constant' ,'Feedback inhibition connections' ,'Somatic PRP threshold' ,'Local PRP Threshold' ,'Dendritic spike threshold' ,'Initial synapse weight' ,'Maximum synapse weight', 'Number of branches per neuron', 'Number of neurons'};

%conds = {'connectivityParam'};
%condLabels = {'Number of branches per neuron'};

PRPs = {'NPERTG', 'NPERTL', 'NPERT'};

boxes = [];
boxesWS = [];


for i=1:length(conds)
    CONDITION=conds{i}
    
    if (1)
        
        PRPs = {'PERT', 'PERTL', 'PERTG'};
        
        for prp=1:length(PRPs)
            PRP=PRPs{prp}
            
            gpvalues={};
            
            ncase=1;
            
            for gpv=0.7:0.1:1.3
                
                if (strcmp(CONDITION,'nNeuronsParam'))
                    npyrs = ceil(gpv*400);
                end
                
                for run=1:nruns
                    sfn=sprintf('%s/%s_%s_%.1f_%d/spikesperpattern.dat', DATADIR, PRP, CONDITION, gpv, run-1)
                    spk = load( sfn);
                    spk = spk(:, 1:npyrs)/(stimduration/1000); %duration
                    totactive(run, ncase, prp, :) = sum(spk>CUTOFF,2)/npyrs;
                    totfiring(run, ncase, prp, 1) = mean(spk(1, :), 2);
                    ffactive(run, ncase, prp, 1) = mean(spk(1, spk(1,:)>CUTOFF), 2);
                    
                end
                
                for run=1:nruns
                    sfn=sprintf('%s/N%s_%s_%.1f_%d/spikesperpattern.dat', DATADIR, PRP, CONDITION, gpv, run-1)
                    spk = load( sfn);
                    spk = spk(:, 1:npyrs)/(stimduration/1000); %duration
                    totactive2(run, ncase, prp, :) = sum(spk>CUTOFF,2)/npyrs;
                end
                gpvalues{ncase} = sprintf('%.1f', gpv);
                ncase=ncase+1;
                
            end
        end
        figure
        
        mact = 100.0*mean(totactive,1);
        sact = 100.0*std(totactive,0,1)/(sqrt(nruns));
        
        boxes(i,:) = mact(:,:,2)';
        
        hold on
        
        errorbar(mact(:,:,1,2), sact(:,:,1,2));
        errorbar(mact(:,:,2,2), sact(:,:,2,2), 'r');
        errorbar(mact(:,:,3,2), sact(:,:,3,2), 'g');
        
        
        if (~strcmp(CONDITION,'connectivityParam') && ~strcmp(CONDITION,'dendSpikeThresh'))
            mact2 = 100.0*mean(totactive2,1);
            sact2 = 100.0*std(totactive2,0,1)/(sqrt(nruns));

            errorbar(mact2(:,:,1,2), sact2(:,:,1,2), '--');
            errorbar(mact2(:,:,2,2), sact2(:,:,2,2), 'r--');
            errorbar(mact2(:,:,3,2), sact2(:,:,3,2), 'g--');
        
        end
        
        hold off
        
        set(gca, 'XTick', [1,4,7], 'XTickLabel',{'70%', '100%', '130%'});
        xlabel(condLabels{i});
        ylim([0,100]);
        set(gca, 'YTick', [0, 20,40,60,80]);
        
        %return
        
        ylabel('% coding neurons');
        export_fig(sprintf('./Nfigs/ALLGP_%s.pdf',CONDITION), '-transparent')
        
    end
    
    if (1)
        
        PRPs = {'PERTWS', 'PERTWSL', 'PERTWSG'};
        
        for prp =1:length(PRPs)
            PRP=PRPs{prp}
            
            gpvalues={};
            ncase=1;
            for gpv=0.7:0.1:1.3
                if (strcmp(CONDITION,'nNeuronsParam'))
                    npyrs = ceil(gpv*400);
                end
                
                for run=1:nruns
                    
                    sfn=sprintf('%s/%s_%s_%.1f_%d/spikesperpattern.dat',  DATADIR, PRP, CONDITION, gpv, run-1)
                    spk = load( sfn);
                    spk = spk(:, 1:npyrs)/(stimduration/1000); %duration
                    spcut = spk>CUTOFF;
                    ftot = sum(spcut(:));
                    if (ftot ==0)
                        ftot = 1;
                    end
                    overlap(run, ncase, prp) = sum(spcut(1,:)&spcut(2,:) )/( ftot/2 );
                end
                
                for run=1:nruns
                    
                    sfn=sprintf('%s/N%s_%s_%.1f_%d/spikesperpattern.dat',  DATADIR, PRP, CONDITION, gpv, run-1)
                    spk = load( sfn);
                    spk = spk(:, 1:npyrs)/(stimduration/1000); %duration
                    spcut = spk>CUTOFF;
                    ftot = sum(spcut(:));
                    if (ftot ==0)
                        ftot = 1;
                    end
                    overlap2(run, ncase, prp) = sum(spcut(1,:)&spcut(2,:) )/( ftot/2 );
                end
                gpvalues{ncase} = sprintf('%.1f', gpv);
                ncase=ncase+1;
                
            end
        end

        figure
        mact = 100.0*mean(overlap,1);
        sact = 100.0*std(overlap,0,1)/(sqrt(nruns));

        %boxesWS(i,:) = mact';

        hold on
        errorbar(mact(:,:,1), sact(:,:,1));
        errorbar(mact(:,:,2), sact(:,:,2), 'r');
        errorbar(mact(:,:,3), sact(:,:,3), 'g');

        if (~strcmp(CONDITION,'connectivityParam') && ~strcmp(CONDITION,'dendSpikeThresh'))
            mact2 = 100.0*mean(overlap2,1);
            sact2 = 100.0*std(overlap2,0,1)/(sqrt(nruns));

            errorbar(mact2(:,:,1), sact2(:,:,1), '--');
            errorbar(mact2(:,:,2), sact2(:,:,2), 'r--');
            errorbar(mact2(:,:,3), sact2(:,:,3), 'g--');
        end

        
        hold off
        %%errorbar(mact, sact);

        %set(gca, 'XTick', [1,6,11], 'XTickLabel',{'50%', '100%', '150%'});
        set(gca, 'XTick', [1,4,7], 'XTickLabel',{'70%', '100%', '130%'});

        xlabel(condLabels{i});
        ylim([0,100]);
        %xlim([3,9]);
        set(gca, 'YTick', [0, 20,40,60,80]);
        ylabel('% overlapping neurons');
        export_fig(sprintf('./Nfigs/ALLGPWS_%s.pdf',CONDITION), '-transparent')

        
        
        if (0)
            close all
            condLabels2 = {'Nin'  ,'tH','tCREB' ,'Ninh' ,'soma' ,'local' ,'dspike' ,'winit' ,'wmax', 'Nbr', 'Nn'};
            figure
            boxplot(boxes')
            ylabel('% coding neurons');
            set(gca, 'XTick', [1:11]);
            set(gca, 'XTickLabel',condLabels2);
            rotateXLabels(gca, 90)
            export_fig(sprintf('./figs/GPWS_BOXES.pdf'), '-transparent')
            
            figure
            boxplot(boxesWS')
            ylabel('% overlapping neurons');
            set(gca, 'XTick', [1:11]);
            set(gca, 'XTickLabel',condLabels2);
            rotateXLabels(gca, 90)
            export_fig(sprintf('./figs/GPWS_BOXESWS.pdf'), '-transparent')
            
        end
    end
end  %sensitivity conds




if (0)
    
    X = [];
    Y=[];
    nrow=1;
    for  CREBTimeParam= 0.5:0.5:1.5
        for  connectivityParam= 0.5:0.5:1.5
            for  globalPRPThresh = 0.5:0.5:1.5
                for  localPRPThresh = 0.5:0.5:1.5
                    for  dendSpikeThresh = 0.5:0.5:1.5
                        for  initWeight = 0.5:0.5:1.5
                            for run = 0:4
                                
                                
                                fn = sprintf('%.1f_%.1f_%.1f_%.1f_%.1f_%.1f_%d',  CREBTimeParam, connectivityParam, globalPRPThresh, localPRPThresh, dendSpikeThresh, initWeight, run)
                                
                                %${CREBTimeParam}_${connectivityParam}_${globalPRPThresh}_${localPRPThresh}_${dendSpikeThresh}_${initWeight}_${run}
                                X(nrow, :) = [ CREBTimeParam, connectivityParam, globalPRPThresh, localPRPThresh, dendSpikeThresh, initWeight, run];
                                
                                sfn=sprintf('./data/SENS_%s/spikesperpattern.dat', fn);
                                spk = load( sfn);
                                spk = spk(:, 1:npyrs)/(stimduration/1000); %duration
                                totactive = sum(spk>CUTOFF,2)/npyrs;
                                Y(nrow,1) = totactive;
                                
                                sfn=sprintf('./data/SENSWS_%s/spikesperpattern.dat', fn);
                                spk = load( sfn);
                                spk = spk(:, 1:npyrs)/(stimduration/1000); %duration
                                spcut = spk>CUTOFF;
                                ftot = sum(spcut(:));
                                if (ftot ==0)
                                    ftot = 1;
                                end
                                Y(nrow,2) = sum(spcut(1,:)&spcut(2,:) )/( ftot/2 );
                                nrow = nrow+1;
                                
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    % ${CREBTimeParam}_${connectivityParam}_${globalPRPThresh}_${localPRPThresh}_${dendSpikeThresh}_${initWeight}_${run}
    
    condLabels = {'tCREB' ,'Nin' ,'soma' ,'local' ,'dspike' ,'winit', 'run'};
    
    x = array2table(X);
    y = array2table(Y);
    
    opt = sdo.AnalyzeOptions;
    opt.Method = 'StandardizedRegression';
    opt.MethodOptions = 'Ranked'
    r = sdo.analyze(x,y, opt)
    
    res =  table2array(r);
    res = res(1:6, :);
    bar(abs(res))
    ylabel('Correlation')
    xlabel('Parameter');
    set(gca, 'XTickLabel',condLabels);
    %rotateXLabels(gca, 60)
    
end
