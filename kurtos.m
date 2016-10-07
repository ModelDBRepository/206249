function kurt = kurtos(actPpost)

mpost   = mean( actPpost(:) );
stdpost = std(  actPpost(:) );

kpost = 0;
for i=1:length(actPpost)
    kpost  = kpost  + ( (actPpost(i) -mpost )/stdpost )^4;
end

kurt = (kpost / length(actPpost)) - 3.0;
