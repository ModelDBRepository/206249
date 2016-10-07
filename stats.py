#!/usr/bin/python -i

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

np.set_printoptions( threshold=999999999999999)

if (len(sys.argv) < 2):
	datadir = os.popen("ls -t ./data | head -1").read().strip();
else:
	datadir = sys.argv[1]
	datadir = datadir.replace('./data/', '');
	datadir = datadir.replace('data/', '');

IMODE=1
if (len(sys.argv) > 2):
	IMODE =0
	

ACTIVE_CUTOFF = 10

#N400B.20.I10.i6.P1.p1.T30.S1980.s
p = re.match(r'N(\d+).B(\d+).I(\d+).i(\d+).P(\d+).p(\d+).T(\d+).S(\d+).(\w+)', datadir)

print datadir
#os.system("cp constructs.cpp ./data/%s/"%( datadir))

STIMDURATION = 1800 + 100 + 100 # Taken from constructs.cpp
NTOTAL       = int(p.group(1))   #inh + pyr neurons
NBRANCHES    = int(p.group(2))
NINPUTS      = int(p.group(3))
NPERINPUT    = int(p.group(4))
NPATTERNS    = int(p.group(5))
INTERSTIM    = int(p.group(7))
SUFFIX    = p.group(9)
NPYR         = int(0.8*NTOTAL)
PYR_IDS      = range(0 , NPYR)
IN_IDS       = range(NPYR, NTOTAL)


def getoverlaps(patarr, threshold):
	ta = patarr > threshold
	cors = np.zeros((ta.shape[0], ta.shape[0]))
	s = ta[i,:].sum() + ta[j,:].sum();
	for i in range(ta.shape[0]):
		for j in range(ta.shape[0]):
			cors[i][j] = 2.0*float((ta[i,:] & ta[j,:]).sum()) / float(ta[i,:].sum() + ta[j,:].sum())
	return cors


spcoeff = None

def getrates(period, patternid=0, rates=1):
	if (period == 'pre'):
		tstart = 0
		tend = tstart + NPATTERNS*STIMDURATION
	elif (period ==  'during'):
		tstart = NPATTERNS*STIMDURATION*0
		tend = tstart + NPATTERNS*STIMDURATION
	elif (period == 'post'):
		tstart = NPATTERNS*STIMDURATION*1
		tend = tstart + NPATTERNS*STIMDURATION
	elif (period=='single'):
		tstart = NPATTERNS*STIMDURATION*3
		tend = tstart + NINPUTS*STIMDURATION

	if (patternid>=0 and period != 'single'):
		tstart = tstart + STIMDURATION*patternid
		tend = tstart + STIMDURATION

	if (rates==1): return  averaged[:, tstart:tend]
	else: return  raster[:, tstart:tend]




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
		raster[nid,0] =0 # XXX bug
		nid += 1

	return raster




def plotdistancesspikes():
	global spcoeff
	global spcounts


	ry = STIMDURATION*NINPUTS + (NINPUTS*(NINPUTS-1))*STIMDURATION
	NMEMS = NINPUTS

	raster = loadtimings("./data/%s/spikes.dat"%(datadir), ry)

	windowsize=5.
	WND = np.repeat(1.0, windowsize) / windowsize
	#averaged = np.zeros(raster.shape)


	counts = np.zeros((NMEMS,NMEMS))
	tot =0
	for i in range(NMEMS):
		for j in range(NMEMS):
			if (i < j):
				tstart = NPATTERNS*STIMDURATION + STIMDURATION*(tot)
				tend = tstart +STIMDURATION 
				r = raster[:, tstart:tend]
				sums = np.sum(r)
				tot += 1
				counts[i][j] = sums

	figure()
	imshow(raster , interpolation='nearest', aspect='auto',cmap='jet')
	colorbar()
	#savefig("./data/%s/raster.png"%(datadir))

	figure()
	imshow(counts , interpolation='nearest', aspect='auto',cmap='jet')
	#colorbar()
	#tight_layout()
	savefig("./data/%s/distances.png"%(datadir))

	return





def plotspikes():
	global spcoeff
	global spcounts
	ry = STIMDURATION*NPATTERNS*2 + 0*NINPUTS*STIMDURATION
	raster = loadtimings("./data/%s/spikes.dat"%(datadir), ry)

	windowsize=50.
	WND = np.repeat(1.0, windowsize) / windowsize
	averaged = np.zeros(raster.shape)
	for nid in range(raster.shape[0]):
		av = np.convolve(raster[nid, :], WND)
		averaged[nid, : ] = av[windowsize/2:-windowsize/2+1]

	figure()
	imshow(averaged *(500.), interpolation='nearest', aspect='auto', cmap='jet')
	colorbar()
	#tight_layout()
	#savefig("./data/%s/raster.png"%(datadir));

	spcounts = np.zeros( (len(PYR_IDS), NPATTERNS) )
	spcoeff = np.zeros( (len(PYR_IDS), NPATTERNS) )
	cors = np.zeros( (NPATTERNS,len(PYR_IDS)) )
	patresp = np.zeros((len(PYR_IDS), NPATTERNS))

	for n in range(NPATTERNS):
		tstart = NPATTERNS*STIMDURATION*1 + STIMDURATION*n
		tend = tstart +STIMDURATION 
		r = raster[:, tstart:tend]
		av = averaged[:, tstart:tend]
		sums = np.sum(r,1)
		corrs = np.nan_to_num(np.corrcoef(av)) # discard NaN values

		for nid in PYR_IDS:
			#s =0;
			#s += cc[nid, NTOTAL + n*NPERINPUT+inpid];
			#cors[n, nid] = s / len(PYR_IDS)
			stt = sums[nid]
			spcounts[nid, n] = stt
			
			#for inpnid in range(NPERINPUT):
				#print "%d %d %d\n"%( NTOTAL, NPERINPUT, NTOTAL + n*NPERINPUT +inpnid)
				#spcoeff[nid, n] += corrs[nid, NTOTAL + n*NPERINPUT +inpnid]/NPERINPUT
	spcounts /= 2.0

	fircorr = np.zeros((len(PYR_IDS), NPATTERNS))

	if (0):
		for nid in [1, 2, 3]:
			nn = np.zeros((NPATTERNS, STIMDURATION))
			for sp in range(NPATTERNS):
				tstart = NPATTERNS*STIMDURATION*1 + STIMDURATION*sp
				tend = tstart +STIMDURATION 
				nn[sp, :] = raster[nid, tstart:tend]

			gg = np.nan_to_num(np.corrcoef(nn))
			figure()
			imshow(gg , interpolation='nearest', aspect='auto', cmap='jet')
			colorbar()
		return;



	spcountspre = np.zeros( (len(PYR_IDS), NPATTERNS) )
	for n in range(NPATTERNS):
		tstart = NPATTERNS*STIMDURATION*0 + STIMDURATION*n
		tend = tstart +STIMDURATION 
		r = raster[:, tstart:tend]
		#av = averaged[:, tstart:tend]
		sums = np.sum(r,1)

		for nid in PYR_IDS:
			#s =0;
			#s += cc[nid, NTOTAL + n*NPERINPUT+inpid];
			#cors[n, nid] = s / len(PYR_IDS)
			spcountspre[nid, n] = sums[nid]




	spcoeff = spcoeff.transpose()

	"""
	figure()
	title("Average firing rate per pattern recall")
	ylabel("number of spikes PRE")
	xlabel("#pattern")
	imshow(spcountspre, interpolation='nearest', aspect='auto',cmap='jet')
	np.save("./data/%s/spikespre"%(datadir), spcountspre);

	colorbar()
	"""

	figure()
	title("Average firing rate (Hz) during pattern recall")
	xlabel("Memory No.")
	ylabel("Pyramidal Neuron No.")
	imshow(spcounts , interpolation='nearest', aspect='auto',cmap='jet')
	np.save("./data/%s/spikespost"%(datadir), spcounts);
	savefig("./data/%s/spikespost.pdf"%(datadir))
	colorbar()


	figure()
	if (NPATTERNS>1):
		pp = np.corrcoef((spcounts).transpose())
		xlabel("Memory #")
		ylabel("Memory #")
		title("Firing rates Correlation")
		imshow(pp.transpose() , interpolation='nearest', aspect='auto',cmap='jet', vmin=-0.2, vmax=1.)
		colorbar()
		np.save("./data/%s/spcountscorr"%(datadir), pp.transpose());
		xticks(arange(NPATTERNS))

		figure()
		activeonly = np.corrcoef((spcounts>ACTIVE_CUTOFF).transpose())
		xlabel("Memory #")
		ylabel("Memory #")
		title("Population overlap")
		imshow(activeonly , interpolation='nearest', aspect='auto',cmap='jet', vmin=-0.2, vmax=1.)
		colorbar()
		#np.save("./data/%s/spcountscorr"%(datadir), pp.transpose());
		xticks(arange(NPATTERNS))




	figure()
	ylabel("% neurons")
	xlabel("# mem")
	title("Percentage of neurons > %d spikes"%(ACTIVE_CUTOFF))
	spc = np.zeros(NPATTERNS)
	for n in range(NPATTERNS):
		spc[n] = (spcounts[:,n] > ACTIVE_CUTOFF).sum()

	vals = 100*spc/NPYR
	print "VV"
	print vals
	bar( arange(NPATTERNS), vals )
	#ylim((0,60))
	np.save("./data/%s/above_%d_dt"%(datadir, ACTIVE_CUTOFF), vals)
	savefig("./data/o_%s_%d.pdf"%( SUFFIX, INTERSTIM));
	xticks(arange(NPATTERNS))


	figure()
	title("avg firing rate per pattern (interval=%d min)"%(INTERSTIM))
	ylabel("avg firing rate (Hz)")
	xlabel("# mem")
	nbins = 20
	pp = zeros( (NPATTERNS, nbins));
	perr = np.std(spcounts.transpose(), 1)
	vv = np.sum(spcounts.transpose(), 1)/(NPYR*1.8) ;

	bar(arange(NPATTERNS), vv)
	ylim((0,50))

	np.save("./data/%s/avgrate"%(datadir), vv)
	xticks(arange(NPATTERNS))

	#	hist = np.histogram(st[p, :], bins=nbins)
	#	pp[p, :] = hist[0]
	#imshow(pp , interpolation='nearest', aspect='auto',cmap='jet')
	#colorbar()


	"""
	figure()
	title("spikes overlap per mem")
	ylabel("#mem")
	xlabel("#mem")
	if 0:
		pp = getoverlaps(spcounts.transpose(), 20)
		bar([0], [pp[0,1]])
		ylim((0.,1.))
	else:
		if (NPATTERNS > 1):
			#pp = np.corrcoef((spcounts).transpose())
			pp = getoverlaps(spcounts.transpose(), ACTIVE_CUTOFF)
			pp = pp.transpose();
			imshow(pp , interpolation='nearest', aspect='auto',cmap='jet', vmin=0, vmax=1.)
			np.save("./data/%s/spoverlap"%(datadir), pp)
			colorbar()
			xticks(arange(NPATTERNS))
			yticks(arange(NPATTERNS))

	#tight_layout()
	savefig("./data/%s/firing.png"%(datadir));
	"""

	#figure()
	#title("Dist. of frequences during recall")
	#hist(spcounts.flatten()/1.800, bins=100)


	"""
	figure()
	title("#patterns / neuron > 30 spikes")
	xlabel("#patterns")
	ylabel("#neurons")
	n, b, l = hist( (spcounts>ACTIVE_CUTOFF).sum(axis=1), bins=12, range =(0,12))
	print n
	print b
	np.save("./data/%s/pat_neuron"%(datadir), n)
	np.save("./data/%s/pat_neuron_dt"%(datadir), (spcounts>ACTIVE_CUTOFF).sum(axis=1))


	figure(figsize=(15,15))
	subplot(221)
	title("#firing correlations coefficients")
	pp = np.zeros( (len(PYR_IDS), NPATTERNS) )
	ylabel("#pattern");
	xlabel("#neuron");
	imshow(spcoeff.transpose(), interpolation='nearest', aspect='auto',cmap='jet', vmin=-.2, vmax=1)
	colorbar()

	subplot(222)
	title("#corr. coeff.")
	ylabel("#pattern");
	xlabel("#pattern");
	if (NPATTERNS>1):
		imshow(np.corrcoef(spcoeff), interpolation='nearest', aspect='auto',cmap='jet',  vmin=-.2, vmax=1)
		colorbar()

	subplot(223)
	title("recruited neurons")
	ylabel("% recruited pyr. neurons");
	xlabel("#pattern");

	n1 = []
	n2 = []
	n3 = []
	for i in range(NPATTERNS):
		n1.append((spcoeff[i,:]>0.2).sum()/float(NPYR))
		n2.append((spcoeff[i,:]>0.3).sum()/float(NPYR))
		n3.append((spcoeff[i,:]>0.5).sum()/float(NPYR))
	
	#plot(n1, label='0.2', marker='o')
	plot(n2, label='0.3', marker='o')
	#plot(n3, label='0.6', marker='o')
	legend()

	subplot(224)
	hst = []
	thresh = 0.3
	for i in range(NPYR):
		hst.append((spcoeff[:,i]>thresh).sum())

	title('# of patterns per neuron (thresh=%f)'%(thresh))
	hist(hst, bins=20)
	xlabel('# of patterns')
	ylabel('# of neurons')

	#tight_layout()
	savefig("./data/%s/firingcorr.png"%(datadir), dpi=200);
	"""



def plotsyns():
	global brhits
	ss = np.loadtxt("./data/%s/synstate.dat"%(datadir))

	brweights= np.zeros((NINPUTS, NPYR*NBRANCHES))
	nrnweights= np.zeros((NINPUTS, NPYR))
	brhits = np.zeros((NPYR*NBRANCHES))

	brpatterns = np.zeros((NBRANCHES*NPYR, NINPUTS))

	learnedhits= np.zeros(NPYR*NBRANCHES)
	unlearnedhits= np.zeros(NPYR*NBRANCHES)

	mem_branches = np.zeros((NPATTERNS, NPYR*NBRANCHES));

	for l in ss:
		bid = l[1]
		nrnid = l[2]
		inputid = l[4]
		srcnrnid=l[3]
		w = l[6]

		if (nrnid < NPYR and inputid >= 0):
			if (w > 0.999):
				brhits[bid] += 1
				brpatterns[bid, inputid] = 1
				if (spcounts[nrnid, inputid] > ACTIVE_CUTOFF):
					learnedhits[bid] += 1
				else:
					unlearnedhits[bid] += 1
			if (w > 1.0):
				mem_branches[inputid, bid] += 1

			brweights[inputid, bid] += w
			nrnweights[inputid, nrnid] += w
	
	wmax =  nrnweights.max()
	brmax =  brweights.max()
	nrnweights /= wmax
	brweights /= brmax

	#figure()
	#title("Number of potentiated synapses per branch");
	#hist(mem_branches[0, :],  range=(1,19), bins=19)

	#m = np.average(mem_branches, axis=1)
	#st = np.std(mem_branches, axis=1)

	#figure()
	#bar(range(0,len(m)), m, yerr=st)

	#figure()
	#title("Potentiated  synapses per branch")
	#imshow(mem_branches.transpose(), interpolation='nearest', aspect='auto',cmap='jet')
	#colorbar()

	if (NPATTERNS>1):
		figure()
		title("Correlation of Potentiated  synapses per neuron")
		pp = np.corrcoef(nrnweights) 
		imshow(pp , interpolation='nearest', aspect='auto',cmap='jet');
		colorbar()

		figure()
		title("Correlation of Potentiated  synapses per branch")
		pp = np.corrcoef(mem_branches) 
		imshow(pp , interpolation='nearest', aspect='auto',cmap='jet', vmin=-0.2, vmax=1.)
		colorbar()







def plotdspikes():
	ry = STIMDURATION*NPATTERNS*2 

	filename = ("./data/%s/branchspikes.dat"%(datadir))
	ff = open(filename, 'r') 
	fdata = ff.readlines()
	sx = len(fdata)

	branch_dspikes = np.zeros((NPATTERNS, NBRANCHES*NPYR))
	bid =0
	for l in fdata:
		ar = np.fromstring(l, sep=' ' , dtype=int)
		for k in range(NPATTERNS):
			print bid
			print (len(ar))
			branch_dspikes[k, bid] += ar[NPATTERNS + k]
		bid += 1
		if (bid >= NPYR*NBRANCHES): break
			
	branch_dspikes = branch_dspikes > 1

	figure()
	xlabel("#pattern")
	ylabel("#branch")
	title("dend. spikes per memory")
	imshow(branch_dspikes.transpose() , interpolation='nearest', aspect='auto',cmap='jet')
	colorbar()

	figure()
	title("corr. coeff ")
	ylabel("#pattern")
	xlabel("#pattern")
	pp = np.corrcoef(branch_dspikes) 
	imshow(pp , interpolation='nearest', aspect='auto',cmap='jet', vmin=-0.2, vmax=1.)
	colorbar()
	#savefig("./data/%s/dspikes.png"%(datadir));

	"""
	figure(figsize=(15,6))
	subplot(121)
	ylabel("#pattern")
	xlabel("#pyr neyron")
	title("dend. spike rate per neuron")
	imshow(nrnspcounts , interpolation='nearest', aspect='auto',cmap='jet')
	colorbar()

	subplot(122)
	title("corr. coeff ")
	ylabel("#pattern")
	xlabel("#pattern")
	pp = np.corrcoef(nrnspcounts) 
	imshow(pp , interpolation='nearest', aspect='auto',cmap='jet', vmin=-.2, vmax=1.)
	colorbar()
	#tight_layout()
	savefig("./data/%s/nrndspikes.png"%(datadir));
	"""




def plotproteins():
	""" ss = np.loadtxt("./data/%s/branchproteins.dat"%(datadir))
	figure(figsize=(15,6))
	ylabel("#branch")
	xlabel("protein")
	title("protein levels per branch ")
	imshow( ss, interpolation='nearest', aspect='auto',cmap='jet')
	colorbar()
	savefig("./data/%s/proteins.png"%(datadir));
	"""

	ss = np.loadtxt("./data/%s/synstate.dat"%(datadir))
	brweights= np.zeros((NPATTERNS, NPYR*NBRANCHES))
	nrnweights= np.zeros((NPATTERNS, NPYR))
	synToPattern = []
	for l in ss:
		bid = l[1]
		nrnid = l[2]
		patternid = l[4]
		srcnrnid=l[3]
		w = l[6]
		synToPattern.append(patternid)

	synsort = []
	for n in range(NPATTERNS):
		k =0
		for l in synToPattern:
			if (n == synToPattern[int(l)]):
				synsort.append( k)
			k+= 1


	
	ss = np.loadtxt("./data/%s/tags.dat"%(datadir))
	print len(synsort)
	figure(figsize=(15,6))
	ylabel("#syn")
	xlabel("tag")
	title("tags per synapse ")

	ab = ss[ss[:,0].argsort()]
	imshow(ab[:,1:], interpolation='nearest', aspect='auto',cmap='jet')
	colorbar()

	savefig("./data/%s/tags.png"%(datadir));

	ss = np.loadtxt("./data/%s/nrnprotein.dat"%(datadir))
	figure(figsize=(15,6))
	ylabel("#nrn")
	xlabel("protein")
	title("proteins per neuron")
	imshow( ss, interpolation='nearest', aspect='auto',cmap='jet')
	colorbar()

	savefig("./data/%s/nrnprotein.png"%(datadir));


	ss = np.loadtxt("./data/%s/weighthistory.dat"%(datadir))
	figure(figsize=(15,6))
	ylabel("#synapse")
	xlabel("weight")
	title("weight per synapse")
	ab = ss[ss[:,0].argsort()]
	imshow( ab[:,1:], interpolation='nearest', aspect='auto',cmap='jet')
	colorbar()

	savefig("./data/%s/nrnweights.png"%(datadir));






def plotbstrengths():
	ss = np.loadtxt("./data/%s/branchstrengths.dat"%(datadir))
	figure(figsize=(15,6))
	#subplot(121)
	ylabel("#branch")
	xlabel("branch strength")
	title("branch strengths ")
	imshow( ss, interpolation='nearest', aspect='auto',cmap='jet')
	colorbar()

	#subplot(122)
	#title("corr. coeff ")
	#ylabel("#branch")
	#xlabel("#branch")
	#pp = np.corrcoef(ss) 
	#imshow( pp , interpolation='nearest', aspect='auto',cmap='jet')
	#colorbar()

	tight_layout()
	savefig("./data/%s/bstrengths.png"%(datadir));



def ansims():
	for sim in range(1,6):
		for i in range(0,NRUNS):
			ddir = "N%d.B%d.I%d.i%d.P%d.p%d.T%d.S%d.sn_sim%d"%(NTOTAL, NBRANCHES, NINPUTS, NPERINPUT, NPATTERNS, NPERPATTERN, INTERSTIM, RSEED+i, sim)
			dd = np.load("./data/%s/spcountscorr.npy"%(ddir))
			p += dd



if (0):
	NRUNS = 10
	p = np.zeros((NPATTERNS, NPATTERNS))
	for i in range(0,NRUNS):
		ddir = "N%d.B%d.I%d.i%d.P%d.p%d.T%d.S%d.%s"%(NTOTAL, NBRANCHES, NINPUTS, NPERINPUT, NPATTERNS, NPERPATTERN, INTERSTIM, RSEED+i, suff)
		
		dd = np.load("./data/%s/spcountscorr.npy"%(ddir))
		p += dd
	p /= NRUNS
	figure()
	imshow(p , interpolation='nearest', aspect='auto',cmap='jet')

else:
	if (IMODE): ion()
	plotspikes()
	plotsyns()
	#plotdspikes()
	#plotproteins()
	#plotbstrengths()
	if (IMODE): show()


