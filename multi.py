#!/usr/bin/python  -i
import matplotlib
from pylab import *
import numpy as np
import re
import os


NPYR = 160

if (len(sys.argv) < 2):
	datadir = os.popen("ls -t ./data | head -1").read().strip();
else:
	datadir = sys.argv[1]
	datadir = datadir.replace('./data/', '');
	datadir = datadir.replace('data/', '');


#N100B.20.I10.i6.P10.p3.T160.S1980
def loadtimings(filename, tduration):
	ff = open(filename, 'r') 
	fdata = ff.readlines()
	sx = len(fdata)
	sy = tduration;
	raster = np.zeros( (sx, sy) );
	nid=0
	for l in fdata:
		ar = np.fromstring(l, sep=' ' , dtype=int)
		raster[nid, ar] = 1
		nid += 1

	return raster


def getdatafromrun(datadir):
	p = re.match(r'N(\d+)B.(\d+).I(\d+).i(\d+).P(\d+).p(\d+).T(\d+).S(\d+)', datadir)
	print datadir
	STIMDURATION = 1600 + 200 + 200 # Taken from constructs.cpp
	NTOTAL       = int(p.group(1))   #inh + pyr neurons
	NBRANCHES    = int(p.group(2))
	NINPUTS      = int(p.group(3))
	NPERINPUT    = int(p.group(4))
	NPATTERNS    = int(p.group(5))
	INTERSTIM    = int(p.group(7))
	NPYR         = int(0.8*NTOTAL)
	PYR_IDS      = range(0 , NPYR)
	IN_IDS       = range(NPYR, NTOTAL)

	ry = STIMDURATION*NINPUTS + (NINPUTS*(NINPUTS-1))*STIMDURATION

	raster = loadtimings("./data/%s/spikes.dat"%(datadir), ry)

	spcounts = np.zeros( (len(PYR_IDS), NPATTERNS) )
	spcoeff = np.zeros( (len(PYR_IDS), NPATTERNS) )
	cors = np.zeros( (NPATTERNS,len(PYR_IDS)) )

	for n in range(NPATTERNS):
		tstart = NPATTERNS*STIMDURATION*1 + STIMDURATION*n
		tend = tstart +STIMDURATION 
		r = raster[0:320, tstart:tend]
		sums = np.sum(r,1)
		#corrs = np.nan_to_num(np.corrcoef(av)) # discard NaN values
		for nid in PYR_IDS:
			#s =0;
			#s += cc[nid, NTOTAL + n*NPERINPUT+inpid];
			#cors[n, nid] = s / len(PYR_IDS)
			spcounts[nid, n] = sums[nid]
			#for inpnid in range(NPERINPUT):
			#	spcoeff[nid, n] += corrs[nid, NTOTAL + n*NPERINPUT +inpnid]/NPERINPUT

	spc = np.zeros(NPATTERNS)
	for n in range(NPATTERNS):
		spc[n] = (spcounts[:, n]>16).sum()
	
	pp = np.corrcoef(spcounts.transpose())

	return [spcounts, spc, pp]



def runweaks():
	nruns = 5
	ptns = 2
	delay = 60
	doruns =0
	#creb = '-c'
	creb = ''

	if (len(sys.argv) >1):
		doruns = 1

	runids = range(1981, 1981+nruns)

	if (doruns):
		for rid in runids:
			os.system("./lamodel -P %d -T %d -S %d %s"%(ptns, delay, rid, creb))

	spcnts = np.zeros((len(runids), ptns))
	sphi = np.zeros((len(runids), ptns))

	cormat1 = np.zeros((len(runids), ptns))
	i=0
	for rid in runids:
		#system("./lamodel -P 2 -T 30")
		res = getdatafromrun("N200B.40.I10.i6.P%d.p1.T%d.S%d"%( ptns, delay, rid))
		spcnts[i,:] = res[0].sum(axis=0)
		sphi[i,:] = res[1]
		cormat1[i,:] = res[2][0]
		i+= 1

	spcnts /= (320.*1.8)
	pp   =  np.average(spcnts, axis=0)
	perr =  np.std(spcnts, axis=0)

	figure()
	subplot(131)
	xlabel("spike rate mem1, mem2")
	ylabel("avg spike rate")
	bar(range(ptns), pp, yerr =perr, facecolor='#777777')
	ylim((0,30))


	sphi /= 320.
	subplot(132)
	xlabel("recruited neurons mem1, mem2")
	ylabel("avg recruited neurons")
	pp   =  np.average(sphi, axis=0)
	perr =  np.std(sphi, axis=0)
	bar(range(ptns), pp, yerr =perr, facecolor='#777777')
	ylim((0,0.6))

	subplot(133)
	xlabel("correlation value")
	ylabel("overlap mem1, mem2")
	pp   =  np.average(cormat1, axis=0)
	perr =  np.std(cormat1, axis=0)
	bar(range(ptns-1), pp[1:], yerr =perr[1:], facecolor='#777777')
	ylim((-1,1))

	show()

runweaks()


