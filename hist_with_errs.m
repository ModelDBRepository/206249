function hh = hist_with_errs(data, bins)

    dd = [];
    de = [];
    
    
    mmin = (min(data(:)))
    mmax = ceil(max(data(:)))
       
    xbins = [mmin:1:mmax];
    
    for i=1:size(data,1)
        ddd = data(i,:);
        dd(i,:) = hist(ddd(:), xbins);
    end
    
    dm = mean(dd,1);
    ds = std(dd,0,1)/sqrt(size(data,1));
    if (mmin==0)
        dm(1) =0;
        
    hh  = barwitherr(ds, dm);
    %set(gca , 'XTick', linspace(0,length(xbins)) );
    %al = get(gca , 'XTickLabel');
    %al = str2num(al);
    %al = al + mmin-1
    %set(gca , 'XTickLabel', al);
    xmin = floor(min(data(data>0)))
    %xlim([xmin, mmax+1]);
end