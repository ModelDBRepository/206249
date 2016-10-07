function kurt = trevrolls(actPpost)

mpost   = mean( actPpost(:) );
stdpost = std(  actPpost(:) );

s1 =0;
s2 =0;
n = length(actPpost);

for i=1:length(actPpost)
    s1 = s1 + (actPpost(i)/n);
    s2 = s2 + (actPpost(i)^2)/n;    
end

kurt = (s1^2)/s2;
kurt = 1-kurt;