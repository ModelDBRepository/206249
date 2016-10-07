import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
mfs=22
matplotlib.rc('xtick', labelsize=mfs) 
matplotlib.rc('ytick', labelsize=mfs) 
matplotlib.rc('axes', labelsize=mfs) 


from pylab import *

import numpy as np
import re
import os

datadir = sys.argv[1]

RUNSTATS = int(sys.argv[2])


datadir = datadir.replace('./data/', '');
datadir = datadir.replace('data/', '');
figprefix = datadir


p = re.match(r'N(\d+).B(\d+).I(\d+).i(\d+).P(\d+).p(\d+).T(\d+).S(\d+).(.*)', datadir)
RSEED  = int(p.group(8))
NTOTAL       = int(p.group(1))   #inh + pyr neurons
NBRANCHES    = int(p.group(2))
NINPUTS      = int(p.group(3))
NPERINPUT    = int(p.group(4))
NPATTERNS    = int(p.group(5))
NPERPATTERN    = int(p.group(6))
INTERSTIM    = int(p.group(7))
NPYR         = int(0.8*NTOTAL)
PYR_IDS      = range(0 , NPYR)
IN_IDS       = range(NPYR, NTOTAL)

suff = p.group(9)

NRUNS = 10


def pltbar(fn, ttitle, xtitle, ytitle, mylim=-1):
	bp = None
	for i in range(0,NRUNS):
		ddir = "N%d.B%d.I%d.i%d.P%d.p%d.T%d.S%d.%s"%(NTOTAL, NBRANCHES, NINPUTS, NPERINPUT, NPATTERNS, NPERPATTERN, INTERSTIM, RSEED+i, suff)
		dd = np.load("./data/%s/%s.npy"%(ddir, fn))
		if (bp == None):
			bp = np.zeros((NRUNS, len(dd)))
		bp[i,:] = dd

	figure()
	title(ttitle)
	xlabel(xtitle)
	ylabel(ytitle)
	m = np.average(bp, axis=0)
	s = np.std(bp, axis=0)
	#print m
	#print s
	bar(range(0,len(m)), m, yerr=s, ecolor='0.3')
	if (mylim == -1):
		ylim([0,max(m) + max(s)*1.2])
	else:
		ylim([0,mylim])

	print fn, "= ", np.average(bp), " +- ", np.std(bp)

	savefig("./data/%s/sum_%s.pdf"%(figprefix, fn))


def plthist(fn, ttitle, xtitle, ytitle, mylim=-1):
	bp = []
	for i in range(0,NRUNS):
		ddir = "N%d.B%d.I%d.i%d.P%d.p%d.T%d.S%d.%s"%(NTOTAL, NBRANCHES, NINPUTS, NPERINPUT, NPATTERNS, NPERPATTERN, INTERSTIM, RSEED+i, suff)
		dd = np.load("./data/%s/%s.npy"%(ddir, fn))
		bp.append(dd)

	figure()
	#title(ttitle)
	xlabel(xtitle)
	ylabel(ytitle)
	#pp = np.array(bp).flatten()
	pp = bp
	bars = np.zeros((len(pp), 20))
	i=0
	for row in pp:
		y,be = np.histogram(row, bins=20, range=(0.,20.))
		bars[i, :] = y
		i = i+1
	print bars

	#y,be = np.histogram(pp, bins=20, range=(0,20))
	#mstd = np.sqrt(y)
	#hist(pp, bins=20, range=(0,20))
	n = np.average(bars,axis=0)
	st = np.std(bars,axis=0)
	#bincenters = 0.5*(be[1:]+be[:-1])
	bincenters = be[:-1]
	
	clc = '#2060e0'
	bar(bincenters, n,  yerr=st, ecolor=clc,facecolor=clc, edgecolor='white' )
	ylim([0,max(n) + max(st)*1.2])
	xlim([0,20])



	#dpp = pp[pp>0]
	dpp = pp
	print fn, "= ", np.average(dpp), " +- ", np.std(dpp)
	savefig("./data/%s/sum_%s.pdf"%(figprefix, fn))



p = np.zeros((NPATTERNS, NPATTERNS))
ovl = np.zeros((NPATTERNS, NPATTERNS))
brpatterns = np.zeros((NRUNS, 12))



for i in range(0,NRUNS):
	ddir = "N%d.B%d.I%d.i%d.P%d.p%d.T%d.S%d.%s"%(NTOTAL, NBRANCHES, NINPUTS, NPERINPUT, NPATTERNS, NPERPATTERN, INTERSTIM, RSEED+i, suff)
	if (RUNSTATS):
		print("python stats.py %s NO"%( ddir))
		os.system("python stats.py %s NO"%( ddir))

	dd = np.load("./data/%s/spcountscorr.npy"%(ddir))
	p += dd

	dd = np.load("./data/%s/spoverlap.npy"%(ddir))
	ovl += dd
	

p /= NRUNS
ovl /= NRUNS


MDISTANCE=1

figure()
#title("Firing rate correlation (interval=%d min)"%(INTERSTIM))
xlabel("Event No")
ylabel("Event No")

av =[]
nav = []
for i in range(1, NRUNS): 
	av.append(p[i-1,i])

	for j in range(0, i): 
		if (j < i-MDISTANCE):
			nav.append(p[i,j])

print "ScCorr= ", np.average(av), " +- " , np.std(av)
print "NScCorr= ", np.average(nav), " +- " , np.std(nav)

imshow(p , interpolation='nearest', aspect='auto',cmap='jet', vmin=-0.0, vmax=1.0)
colorbar()
savefig("./data/%s/sum_sc_corr.pdf"%(figprefix))


figure()
#title("Active populations correlation (interval=%d min)"%(INTERSTIM))
xlabel("Event No")
ylabel("Event No")

av =[]
nav = []
for i in range(1, NRUNS): 
	av.append(ovl[i-1,i])
	for j in range(0, i): 
		if (j < i-MDISTANCE):
			nav.append(ovl[i,j])


print "Overlap= ", np.average(av), " +- " , np.std(av)
print "Nonoverlap= ", np.average(nav), " +- " , np.std(nav)

imshow(ovl , interpolation='nearest', aspect='auto',cmap='jet',vmin=-0.0, vmax=1.)
colorbar()
savefig("./data/%s/sum_over.pdf"%(figprefix))



#plthist('above_30_dt', "Percentage of neurons active during recall (interval=%d min)"%(INTERSTIM), "Event No", "Percentage of pyr. neurons", 100)

plthist('br_hits_dt', "Potentiated synapses per branch (interval=%d min)"%(INTERSTIM), "Number of potentiated synapses", "Number of branches")

plthist('br_patterns_dt', "Patterns represented per branch (interval=%d min)"%(INTERSTIM), "Number of events", "Number of branches")

plthist('pat_neuron_dt', "Number of neurons responsive to N patterns (interval=%d min)"%(INTERSTIM), "Number of events", "Number of neurons")


if (len(sys.argv)  < 4):
	show()


