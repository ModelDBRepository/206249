import matplotlib
import pylab
import numpy
import NeuroTools.signals as signals

def visalloc():
	for it in [60]:
		for nm in [10,19]:
			for nbranches in range(10, 60, 16):

				coss = numpy.zeros((nm, nm))
				pcos = numpy.zeros((nm, nm))

				pat = numpy.genfromtxt("./data/alloc/p_m%di%db%d.dat"%(nm, it, nbranches))

				ar = numpy.genfromtxt("./data/alloc/r_m%di%db%d.dat"%(nm, it, nbranches))
				ar = ar[:, 0:64];


				scx = []
				scy = []
				for i in range(nm):
					for j in range(i+1):
						s1 = numpy.sqrt(numpy.dot(ar[i], ar[i]))
						s2 = numpy.sqrt(numpy.dot(ar[j], ar[j]))
						coss[i][j] = numpy.dot(ar[i], ar[j]) / (s1*s2)

						s1 = numpy.sqrt(numpy.dot(pat[i], pat[i]))
						s2 = numpy.sqrt(numpy.dot(pat[j], pat[j]))
						pcos[i][j] = numpy.dot(pat[i], pat[j]) / (s1*s2)

						scx.append( pcos[i][j])
						scy.append( coss[i][j])

						

				pylab.figure()
				pylab.title(" Mems %d intvl %d branches %d" %(nm, it, nbranches));
				pylab.imshow(coss , interpolation='nearest', aspect='auto',cmap='hot')
				pylab.imshow(pcos , interpolation='nearest', aspect='auto',cmap='hot')
				pylab.colorbar();
				pylab.savefig("cos_%di%db%d.png"%(nm, it, nbranches));

				#pylab.scatter(scx, scy)
				#pylab.title("m %d i %d B%d"%(nm, it, nbranches));
				#pylab.xlabel("Input Pattern similarity ")
				#pylab.ylabel("Firing pattern similarity ")
				#pylab.savefig("bn_%di%db%d.png"%(nm, it, nbranches));
	pylab.show()

			
			




def spikestats():
	
	ff = open('./data/0/spikes.dat', 'r') 
	fdata = ff.readlines()

	lst = []
	nid=0
	for l in fdata:
		ar = numpy.fromstring(l, sep=' ')
		for a in ar:
			lst.append( (nid, a))
		nid += 1

	slist = signals.SpikeList(lst, range(0,nid))
	slist.raster_plot()

	#pylab.figure()
	#pylab.imshow(spraster , interpolation='nearest', aspect='auto',cmap='hot')


	#pp = numpy.correlate(spraster[100, :], spraster[102,:], 'full');
	#print pp;
	#pylab.plot(pp)

	pylab.show()





def printltp(): 
	lar = numpy.genfromtxt('/tmp/eltp.dat') 
	lar = lar.transpose()
	pylab.figure();
	print len(lar);

	labels = ['Induction', 'ELTP', 'Rb', 'Pb', 'Rn', 'Pn','iltp', 'W', 'T']

	for i in range(len(lar)):
		if ( i != 6):
			pylab.plot(lar[i, :], label=labels[i])


	pylab.legend()

	mmax = len(lar[0, :]) 
	tstep = 50
	th = 3600/tstep


	locs, labels = pylab.xticks( numpy.arange(0, mmax, th), numpy.arange(0, mmax/th, 1))

	#pylab.plot(lar[1, :] + lar[3, :])

	pylab.show()




spikestats()
