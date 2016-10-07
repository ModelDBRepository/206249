xx = [0:240]
y1 = exp(-xx/(1*60))
%y1 = 1 - xx./(z*60)
%y1(find(y1<0)) = 0.

y2 = (xx/30) .*exp(1-xx/(30))
y4 = [zeros(1,240), y2, zeros(1,240)]

close all
figure; 
plot(y1);hold on; plot(y2); hold off;
figure
d = [];
for i=1:480
    y3 = zeros(size(y4));
    y3(i:i+241-1) = y1;

    s = sum(y3.*y4);
    d(i) = s;
end

plot(d)
