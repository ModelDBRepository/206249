
import sys
import matplotlib
import pylab
import numpy


bdir = sys.argv[1]

ar = numpy.genfromtxt("./data/%s/memorytest.dat"%(bdir)) 
nins =  len(ar)

pops  = numpy.zeros( (nins,nins));
coss  = numpy.zeros( (nins,nins));
coact = numpy.zeros( (nins,nins));



for i in range(nins):
	for k in range(len(ar[i])):
		if (ar[i][k] < 7.0): ar[i][k] =  0.0;
			#else: ar[i][k]  =0.0

print ar

for i in range(nins):
	for j in range(nins):
		pops[i][j] =  numpy.dot(ar[i]  , ar[j] ) 
		s1 = numpy.sqrt(numpy.dot(ar[i], ar[i]))
		s2 = numpy.sqrt(numpy.dot(ar[j], ar[j]))
		coss[i][j] = numpy.dot(ar[i], ar[j]) / (s1*s2)

pylab.figure();
pylab.imshow(pops, interpolation='nearest', aspect='auto',cmap='hot')
pylab.colorbar()

pylab.figure();
pylab.imshow(coss, interpolation='nearest', aspect='auto',cmap='hot')
pylab.colorbar()


dr = numpy.genfromtxt("./data/%s/synweights.txt"%(bdir)) 
maxbid = max(dr[:,1]);
maxinpid = max(dr[:,3]);

myar =  numpy.zeros((maxinpid+1,maxbid+1));

for i in range(len(dr)):
	if (dr[i,3] >=0):
		v = dr[i,4]
		if (v <= 0.4): v = 0.0
		myar[dr[i,3], dr[i,1]] += v

maxinpid += 1
coss  = numpy.zeros( (maxinpid,maxinpid));


for i in range(int(maxinpid) ):
	for j in range(int(maxinpid)):
		s1 = numpy.sqrt(numpy.dot(myar[i], myar[i]))
		s2 = numpy.sqrt(numpy.dot(myar[j], myar[j]))
		#if (i==j): coss[i][j] = 0.0
		#else: 
		coss[i][j] = numpy.dot(myar[i], myar[j]) / (s1*s2)


pylab.figure();
pylab.imshow(coss, interpolation='nearest', aspect='auto',cmap='hot')
pylab.colorbar()



pylab.figure();
for i in [4]: pylab.plot(coss[i, :], label='Mem #4 - Overlap');

pylab.figure();
for i in [1,5,9]: pylab.plot(coss[i, :], label='3 Mems / Overlap');

pylab.show()


