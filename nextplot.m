function nextplot(varargin)

    persistent plotno
    persistent nx
    persistent ny
    if (nargin>0)
        figure;
        

        plotno=1;
        nx = varargin{1};
        ny = varargin{2};
        subplot(nx, ny,plotno);
    else
        plotno=plotno+1;
        subplot(nx, ny, plotno);
    end

end