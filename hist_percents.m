function [hh, dm] = hist_percents(data, ytotal)

    dd = [];
    de = [];
    
    flat = data(:);
    flat = flat(:);
    
    mmin = floor(min(flat));
    mmax = ceil(max(flat));
    

    xbins = [1:1:mmax]
    
    
    for i=1:size(data,1)
        ddd = data(i,:);
        dd(i,:) = hist(ddd(ddd>0), xbins)/(ytotal/100);
    end
    
    dm = mean(dd,1);
    ds = std(dd,0,1)/sqrt(size(data,1));
    
    xmin=mmin;
    if (mmin==0)
        %dm(1) =0;
        %ds(1) =0;
        %xmin=xmin+1;
    end

    hh  = barwitherr(ds, dm);
    %set(gca , 'XTick', linspace(0,length(xbins)) );
    %set(gca , 'XTickLabel', xbins);
    
    %al = get(gca , 'XTickLabel');
    %al = str2num(al);
    %al = al + mmin-1
    %al(al<0) = '';
    %set(gca , 'XTickLabel', al);
    %xmin = floor(min(data(data>0)))
    %xlim([xmin, mmax]);
end