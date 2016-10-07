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

datadir = 'N400.B20.I12.i1.P2.p6.T60.S1980.sc_sim0'
RUNSTATS = 0


datadir = datadir.replace('./data/', '');
datadir = datadir.replace('data/', '');
savedir = datadir


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


if (1):
	bp = zeros((6, 4, NRUNS) )
	n1 = 0
	figure()
	for sim in [1,2,3,4,5,6]:
		n2 = 0
		pp = []
		pe = []
		for itvl in [60, 120, 180, 300]: 
			for i in range(0,NRUNS):
				ddir = "N%d.B%d.I%d.i%d.P%d.p%d.T%d.S%d.sn_sim%d"%(NTOTAL, NBRANCHES, NINPUTS, NPERINPUT, NPATTERNS, NPERPATTERN, itvl, RSEED+i, sim)
				dd = np.load("./data/%s/spcountscorr.npy"%(ddir))
				bp[n1, n2, i] = dd[0,1]
			pp.append( np.average(bp[n1,n2]) )
			pe.append( np.std(bp[n1,n2]))
			n2 += 1
		h = errorbar(range(len(pp)), pp, yerr=pe, label="%d %%"%(100*(6-sim)/6))
		ylim(-0.1, 1.2)
		legend()
		n1 += 1
		
	lab = ["1 hour","2 hours", "3 hours", "5 hours"]
	xticks(range(len(pp)), lab)
	legend(loc="center right")
	xlim(-0.2, 4.2)


show()


