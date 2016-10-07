/* 
    Version: $Id: constructs.cpp 172 2014-02-12 10:06:07Z gk $
    LAModel main imlementation file

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "constructs.h"
//#include "wxglmodel.h"
#include <iterator>

const int DEBUG_NID = 91;
const int DEBUG_BID = 2270;
const int DEBUG_SID = 3024;
const int CUTOFF = 10.0;

int LANetwork::RSEED = 1980;


static void updminmax(float&  cmin, float& cmax, float val)
{
	if (val < cmin )
		cmin = val;
	if (val > cmax )
		cmax = val;
}





#define VEC_REMOVE(vec, t) (vec).erase(std::remove((vec).begin(), (vec).end(), t), (vec).end())

inline float sigmoid(float x, float x0, float broad)
{
	return (1.0 / (1.0 + exp(-x - x0)/broad));
}




inline double nuclearproteinalpha(float x)
{
	return (x>20.)*((x-20.*60.)/(30.*60)) * exp(1. - (x-20.*60. )/(30.*60.));

	//double d =  (x-60.*60)/(40.0*60.);
	//return (exp(-d*d));

	//double d =  ((x-3.*60)/(23.*60)) * exp(1. - (x-(3.*60.) )/(23.*60.)); // 4. or 6.
	//if (d >0.) return d;
	//else return 0.;
	
}




inline double branchproteinalpha(float x)
{
	return ((x)/(15.*60)) * exp(1. - (x )/(15.*60.));

	
}


inline void adjust_with_bounds(float& val, float direction, float max, float min)
{
	if (direction>=0.0) val += direction*(max-val);
	else val += direction*(val - min);
}


inline float step_with_bounds(float cur_val, float dir, float max, float min)
{
	if (dir >0) return dir*(max - cur_val);
	else  return dir*(cur_val - min);
}




inline float caDP(float x)
{
	//return x >0.2 ? 1.0 : 0.0;
	//float f =  (2.0/(1.+exp(-(x*10.-3.5)*10.))) - (1.0/(1.+exp(-(x*10.0-0.5)*19.)));


	float f = (1.3/(1.+exp(-(x*10.-3.5)*10.))) - (0.3/(1.+exp(-(x*10.0-2.0)*19.)));
	return f;
	//return f;//(2.0/(1.+exp(-(x*10.-0.7)*10.))) - (1.0/(1.+exp(-(x*10.0-0.2)*19.)));
	//return (2.0/(1.+exp(-(x*10.-3.1)*8.))) - (1.0/(1.+exp(-(x*10.0-2.1)*6.)));
}


inline void program_input(nrn_list lst, int tstart, int duration, float freq, float randomness, int limitActive = -1)
{
	int tot =0;
	int skip = 0;
	if (limitActive == -999)
		skip = 3;
	for (nrn_iter n = lst.begin(); n != lst.end(); ++n)
	{
		if (skip>0)
		{
			skip--;
			continue;
		}

		LAInput* in = (LAInput*)(*n);
		in->Program(tstart, duration, freq, randomness);
		if (limitActive >0 && ++tot >= limitActive)
			return;
	}
}


void LANetwork::CreateNeurons(int number, int n_branches_per_neuron, char type, vector<LANeuron*>* appendTo = 0, int inputId =-1, bool sepTypes = 0)
{
	for (int i =0 ;i < number; i++)
	{
		LANeuron* n ;
		if (type == 'S')
		{
			n = new LAInput();
			n->network = this;

			n->input_id = inputId;
			((LAInput*)n)->groupIdx = i;
		}
		else
		{
			n = new LANeuron();
			n->network = this;
			for (int tt =0; tt < n_branches_per_neuron; tt++)
			{
				LABranch* bb  = new LABranch;
				bb->bid = this->branches.size();
				bb->neuron = n;

				this->branches.push_back(bb);
				n->branches.push_back(bb);
			}
		}
		n->type = type;

		n->V = param_V_rest +1.0;
		n->nid = this->neurons.size();

		n->glx = 1.9*(n->nid%50)/50. - 0.9;
		n->gly = 1.9*(n->nid/50)/50. - 0.9;

		this->neurons.push_back(n);
		if (appendTo)
			appendTo->push_back(n);
	}
}



void LANetwork::CalculateDistances()
{
	for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
		for (nrn_iter nb = neurons.begin(); nb != neurons.end(); ++nb)
			if (nb != na)
			{
				LANeuron* a = *na, *b=*nb;
				
				float dx = b->pos_x - a->pos_x;
				float dy = b->pos_y - a->pos_y;
				float dz = b->pos_z - a->pos_z;
				distances[pair<int,int>(a->nid,b->nid)] = (double)sqrt(dx*dx + dy*dy + dz*dz);
			}
}



int LANetwork::ConnectByDistance(vector<LANeuron*> fromList, vector<LANeuron*> toList, float fromDistance, float toDistance, int nNeuronPairs, int nSynapsesPerNeuron, float weight, bool isPlastic= false, bool randomizeweight = false, float overlap =-1.0)
{
	int tpairs =0;
	while(true)
	{
		LANeuron* a = fromList.at(int(rgen()*(float)fromList.size()));
		LANeuron* b = toList.at(int(rgen()*(float)toList.size()));

		if (1) 
		{
			for (int i =0; i < nSynapsesPerNeuron; i++)
			{
				LASynapse* syn = new LASynapse();
				syn->sid  = this->synapsesCounter++;
				syn->source_nrn = a;
				syn->target_nrn = b;
				syn->isPlastic = isPlastic;

				if (randomizeweight)
					syn->weight = rgen()* 0.5 + 0.2;
				else
					syn->weight  = weight;


				float rval;
				if (0 && overlap >=0.0) // ONLY FOR INPUT NEURONS
				{
					LAInput* in = (LAInput*) a;
					if (in->groupIdx<3) // CS
						rval = rgen()*overlap*float(syn->target_nrn->branches.size());
					else
						rval = (1-rgen()*overlap)*float(syn->target_nrn->branches.size());
				}
				else
				{
					rval = rgen()*float(syn->target_nrn->branches.size());
				}

				syn->target_branch = syn->target_nrn->branches[(int)rval];

				syn->source_nrn->outgoing.push_back(syn);
				syn->target_nrn->incoming.push_back(syn);
				syn->target_branch->synapses.push_back(syn);
				syn->pos = rgen();
				this->synapses.push_back(syn);
			}
		}

		if (++tpairs >= nNeuronPairs) break;
	}
	return tpairs;
}




inline int randomindex(int max)
{
	return (int)(rgen()*float(max));
}

int LANetwork::PurgeInputSynapses(int totalToRemove,float weightLimit)
{

	int totalFound =0;
	int totalTries = totalToRemove*3;
	while (totalFound < totalToRemove && totalTries-->0)
	{
		nrn_list lst = this->inputs_cs.at(randomindex(this->inputs_cs.size()));
		LANeuron* n = lst.at(randomindex(lst.size()));
		if (n->outgoing.size())
		{
			LASynapse* s = n->outgoing.at(randomindex(n->outgoing.size()));
			if (s->target_nrn->type == 'P'  && s->weight <= weightLimit)
			{
				//candidate for deletion

				//pthread_mutex_lock(&this->synapses_mutex);

				VEC_REMOVE(s->source_nrn->outgoing, s);
				VEC_REMOVE(s->target_nrn->incoming, s);
				VEC_REMOVE(s->target_branch->synapses, s);
				VEC_REMOVE(this->synapses, s);
				//cout << " Removing " << s->sid << endl;
				delete s; 

				//pthread_mutex_unlock(&this->synapses_mutex);

				totalFound++;
			}
		}
	}

	if (totalFound < totalToRemove)
	{
		cout << " Warning: not enougn synapses to remove: " << totalToRemove << endl;  
	}

	return totalFound;
}



int LANetwork::CreateInputSynapses( int totalToAdd)
{
	
	int totalFound = 0;
	while (totalFound < totalToAdd)
	{
		this->ConnectByDistance(this->inputs_cs[randomindex(this->inputs_cs.size())], this->pyr_list, 0.0, 10.0,  1, 1, initWeight, true);
		totalFound++;
	}
	return totalFound;
}



void LANetwork::CreateFearNet(int nneurons, int nbranches, int ninputs, int nneuronsperinput)
{
	this->n_neurons = nneurons;
	this->n_inputs = ninputs;
	this->n_neurons_per_input = nneuronsperinput;
	this->n_branches_per_neuron = nbranches;

	this->CreateNeurons(this->n_neurons*0.8, this->n_branches_per_neuron, 'P', &this->pyr_list);

	this->CreateNeurons(this->n_neurons*0.2, 1 , 'I', &this->in_list);
	
	for (int i=0;i < n_inputs;i++)
	{
		vector<LANeuron*> nrnsCS, nrnsNoise;
		CreateNeurons(this->n_neurons_per_input, 0, 'S', &nrnsCS, i, true);
		this->inputs_cs.push_back(nrnsCS);
	}


	CreateNeurons(10, 0, 'S', &this->bg_list, -1, true);

	this->CalculateDistances();

	if (1) 
	{
	/*
		int nconn = 1.0*float(this->pyr_list.size()*this->in_list.size());
		nrn_list pl, nl;
		for (int i=0; i < nconn; i++)
		{
			pl.clear();
			nl.clear();
			pl.append(this->pyr_list.at(randomindex(this->pyr_list.size())));
			nl.append(this->in_list.at(randomindex(this->pyr_list.size())));
		}
		*/

		//ConnectByDistance(this->pyr_list, this->pyr_list, 0., 10., 0.5*float(this->pyr_list.size()*this->pyr_list.size()), 1, 1.0, false);

		ConnectByDistance(this->pyr_list, this->in_list, 0., 10., this->inhibitionParam*0.4*float(this->pyr_list.size()*20), 1, 0.5, false);

		ConnectByDistance(this->in_list, this->pyr_list, 0.0, 10.,this->inhibitionParam*0.6*float(this->pyr_list.size()*20), 1, 0.5, false);

	}

	// CS inputs  projections
	for (int i =0; i < this->n_inputs ; i++)
	{
		// For clustered input data see Takahashi_2012
		//if (this->altConnectivity) this->ConnectByDistance(this->inputs_cs[i], this->pyr_list, 0.0, 10.0, 1.6* 1.06*float(this->branches.size()), 1, INITSTRENGTH, true);
		//else

		//1.8
		// 
		
		this->ConnectByDistance(this->inputs_cs[i], this->pyr_list, 0.0, 10.0, this->connectivityParam*this->synapseMult*1.6*float(this->pyr_list.size()*20), 1, initWeight, true, false, this->branchOverlap);
		//
		
		//this->ConnectByDistance(this->inputs_cs[i], this->pyr_list, 0.0, 10.0,1.3* 1.06*float(this->branches.size()), 1, INITSTRENGTH, true, false, this->branchOverlap);

		/* 
		Quote:
		-Synaptic inputs are not clustered in fast-spiking parvalbumin-positive interneurons
		-Inhibitory interneurons exhibited highly localized alcium activity in their aspiny dendrites, 
			as reported previously (S10, 11). Because they are highly excitable and can fire in response to a 
			single excitatory input (S12), they may not require dendritic integration via synaptic clustering. 
		For another type of interneuron, see ref. S10.
		Takahashi et al. Science 2012 (supplement figure S6) 
		*/

	}


	ConnectByDistance(this->bg_list, this->pyr_list, 0.0, 10., 0.2*float(this->pyr_list.size()*20), 1, 0.6, false);

	//this->spikesFile = fopen("./data/0/spikes.dat", "w");
	
	/*&
	for (int i=0; i < 0.20*this->pyr_list.size(); i++)
	{
		//LANeuron* a = this->pyr_list.at(int(rgen()*(float)this->pyr_list.size()));
		LANeuron* a = this->pyr_list.at(i);
		if (!this->disableCreb)
		{
			//a->crebLevel = 1.0;
		}
		//a->prpTransients.push_back( pair<float,float>(0, 1.) );
	}	
	*/

	this->RecordInitialWeightSums();
}



void LANetwork::RunPattern(vector<int>& pattern, float hifreq,  float lowfreq, int duration, int limitActive )
{

	int pad = (4000 - duration)/2;
	//printf("Running pattern: [ ");


	for (size_t j=0; j < pattern.size(); j++)
	{
		if (pattern[j]>0)
		{
			
			double ts =0;
			for (syn_iter si=  this->synapses.begin() ; 
				 si != this->synapses.end() ;
				si++)
			{
				LASynapse* s = *si;
				if (s->source_nrn->input_id == (int)j)
				{
					ts+= s->weight;
				}
			}
			printf("Total synaptic drive: %f\n", ts);


			program_input(this->inputs_cs[j], pad, duration, hifreq, 0.5, limitActive);
		}
		else
			program_input(this->inputs_cs[j], pad, duration, lowfreq, 0.5, limitActive);

		//printf(" %d", pattern[j] ? 1 : 0);
	}

	program_input(this->bg_list, pad, duration, 0.5, 0.5, -1);

	//printf("] @%f Hz dur=%d \n", hifreq, duration);

	this->ResetSpikeCounters();
	this->StimDynamics(duration+pad+pad);
	
	int tActive =0;
	int tSpikes =0;
	
	for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
	{
		LANeuron* n = *ni;
		if (float(n->total_spikes) /4000.> CUTOFF )
			tActive++;
		tSpikes += n->total_spikes;
	}
	
	printf("Active: %d (%.2f%%), avgF: %.2f\n", tActive, 100.0*float(tActive)/float(this->pyr_list.size()), tSpikes/(float(this->pyr_list.size())*2));
}



void LANetwork::RunStore2(int n_patterns, int activePerPattern, int delayMinutes, int testMode, int patternsOverlapping)
{
	vector<int> pattern;

	int inpfrequency=30;
	int lowfreq = 0;
	int multiruns = 1;


	int n_ones = activePerPattern;
	printf("Running net with %d pyr. neurons, %d branches, %d synapses [%s,%s]\n", (int)this->pyr_list.size(),  (int)branches.size(),  (int)synapses.size(), localProteins ? "Local" : "Global", disableCreb ? "-CREB" : "+CREB");
	this->mfile.open("./data/bdata.dat");
	this->vfile.open("./data/bvoltage.dat");


	fill(pattern.begin(), pattern.end(), 0);
	for (int i =0; i < this->n_inputs; i++)
		pattern.push_back( (i < n_ones) ? 1 : 0);


	for (int i =0; i < n_patterns; i++)
	{
		if (patternsOverlapping>=0)
		{
			fill(pattern.begin(), pattern.end(), 0);
			for (int k =0;  k < n_ones; k++)
				pattern[i*patternsOverlapping+k] =1;
		}
		else
		{
			fill(pattern.begin(), pattern.end(), 0);
			if (repeatedLearning)
				pattern[0] = 1;
			else
				pattern[i] = 1; 
		}
		this->patterns.push_back( pattern );

		cout<< "Pattern  " << i<< " : [";
		copy(pattern.begin(), pattern.end(), ostream_iterator<int>(cout, " "));
		cout << "]" << endl;

	}


	int np =0;

	this->Interstim(1*60);

	PrintSynapsesSnapshot( datadir + "/syn-pre.txt");

	this->ReportSumWeights();

	if (this->pretraining)
	{
		this->runningMode = RUN_PRE;
		this->enablePlasticity = false;
		for (int nr=0; nr < multiruns; nr++)
		{
			np =0;
			for (vector<vector<int> >::iterator it = this->patterns.begin(); it != this->patterns.end(); it++)
			{
				pattern = *it;
				cout<< "Pretraining" << np<< endl;
				this->runningPatternNo = np;

				RunPattern(pattern, inpfrequency, lowfreq, stimDurationParam*3800., 3);

				cout << "tiny interstim  ..." << endl;
				this->Interstim(5*60);
				cout << "done" << endl;
				np++;
			}
		}
		
	}

	cout << "Training .. " << endl;
	this->enablePlasticity = true;
	np =0;
	this->runningMode = RUN_TRAINING;
	for (vector<vector<int> >::iterator it = this->patterns.begin(); it != this->patterns.end(); it++)
	{
		pattern = *it;
		cout<< "Training" << np<< endl;
		this->runningPatternNo = np;

		if (std::find(isWeakMem.begin(), isWeakMem.end(), np) != isWeakMem.end())
		{
			printf("Weak: %d\n", np);
			RunPattern(pattern, inpfrequency, lowfreq, 2700, -1);
		}
		else
			RunPattern(pattern, inpfrequency, lowfreq,  3800, -1);

		cout << "Interstim " << delayMinutes << " mins ..." << endl;
		//this->SaveCalcs();
		this->Interstim(delayMinutes*60);
		cout << "done" << endl;
		np++;
		this->ReportSumWeights();
	}


	this->runningMode = RUN_TEST;
	cout << "Large interstim ..." << endl;
	this->enablePlasticity = false;
	for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
	{
		LANeuron* nrn = *na;
		nrn->crebLevel = 0.0;
	}
	this->Interstim((int)(this->homeostasisTime*3600.));
	cout << " done" << endl;


	this->ReportSumWeights();

	this->isRecall = true;
	cout << "Recall .. " << endl;
	np =0;
	this->enablePlasticity = false;
	int n =0;
	


	PrintSynapsesSnapshot(datadir + "/syn-post.txt");

	this->spikesPerPattern.resize(n_patterns);
	for (int i=0; i < n_patterns; i++)
		this->spikesPerPattern[i].resize(this->neurons.size());


	this->enablePlasticity = false;

	for ( vector<vector<int> >::iterator it = patterns.begin(); it != patterns.end(); it++ )
	{
		pattern = *it;
		n++ ;
		for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
		{
			LANeuron* nrn = *na;
			nrn->crebLevel = 0.0;
		}

		cout<< endl<< "Testing " << np<< endl;
		this->runningPatternNo = np;
		RunPattern(pattern, inpfrequency, lowfreq, 3800, 3);
		this->Interstim(60);

		for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
		{
			LANeuron* nrn = *na;
			this->spikesPerPattern[np][nrn->nid] = nrn->total_spikes;
		}
		np++;
	}

	/*

	np =0;
	n=0;


	for ( vector<vector<int> >::iterator it = patterns.begin(); it != patterns.end(); it++ )
	{
		pattern = *it;
		n++ ;

		for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
		{
			LANeuron* nrn = *na;
			nrn->crebLevel = 0.0;
		}


		cout<< endl<< "Testing US " << np<< endl;
		this->runningPatternNo = np;
		RunPattern(pattern, inpfrequency, lowfreq, 3800, -999);

		this->Interstim(60);
		for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
		{
			LANeuron* nrn = *na;
			this->spikesPerPattern[np][nrn->nid] = nrn->total_spikes;
		}
		np++;
	}

	*/


	this->mfile.close();

}


template <typename T> static void PrintVector( vector<T>&  ar, ostream& outfile) 
{
	for (typename vector<T>::iterator it = ar.begin(); it != ar.end(); it++)
	{
		outfile << *it << ' ';
	}
	outfile << std::endl;
}





void LANetwork::StoreDataFiles( bool extras = false )
{
	vector<int> pattern;

	string dirname = this->datadir;
	ofstream paramsdat((dirname + "/parameters.txt").c_str());
	paramsdat <<"total_neurons="<< this->neurons.size() << endl;
	paramsdat <<"total_pyramidals="<< this->pyr_list.size() << endl;
	paramsdat <<"branches_per_neuron="<< this->n_branches_per_neuron << endl;
	paramsdat <<"number_inputs="<< this->inputs_cs.size() << endl;
	paramsdat <<"neurons_per_input="<< this->n_neurons_per_input << endl;
	paramsdat <<"rseed="<< RSEED << endl;


	ofstream patternsdat((dirname + "/patterns.txt").c_str());
	for (vector<vector<int> >::iterator it = this->patterns.begin(); it != this->patterns.end(); it++)
	{
		pattern = *it;
		copy(pattern.begin(), pattern.end(), ostream_iterator<int>(patternsdat, " "));
		patternsdat << endl;
	}

	ofstream synstatedat((dirname + "/synstate.dat").c_str());

	ofstream spikesdat((dirname + "/spikes.dat").c_str());
	ofstream crebdat((dirname + "/creb.dat").c_str());
	ofstream voltagedat((dirname + "/voltages.dat").c_str());
	ofstream branchspikesdat((dirname +"/branchspikes.dat").c_str());
	ofstream branchcalcium((dirname + "/branchcalcium.dat").c_str());
	ofstream weightsdat((dirname + "/weights.dat").c_str());
	ofstream branchproteins((dirname + "/branchproteins.dat").c_str());
	ofstream branchstrengths((dirname + "/branchstrengths.dat").c_str());
	ofstream tagsdat((dirname + "/tags.dat").c_str());
	ofstream nrnproteindat((dirname + "/nrnprotein.dat").c_str());
	ofstream weighthistorydat((dirname + "/weighthistory.dat").c_str());
	ofstream dbgneuron((dirname + "/dbgneuron.dat").c_str());

	PrintVector<float>( dbgNeuron, dbgneuron);

	for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
	{
		LANeuron* nrn = *na;

		PrintVector<int>( nrn->spikeTimings, spikesdat);
		//PrintVector<float>( nrn->crebHistory, crebdat);
		PrintVector<float>( nrn->proteinHistory, nrnproteindat);

		for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
		{
			LABranch* b = *bi;

			PrintVector<float>(b->branchSpikesHistory, branchspikesdat);

			for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
			{
				LASynapse* s = *si;
				//PrintVector<float>(s->weightHistory, weightsdat);
				synstatedat << s->sid<<" " 
					<< b->bid<<" "
					<< nrn->nid  << " "
					<< s->source_nrn->nid <<" " 
					<< s->source_nrn->input_id<< " "
					<< b->strength  << " "
					<< s->weight << " " 
					<<endl;

				if (s->tagHistory.size())
				{
					tagsdat << s->source_nrn->input_id<< " ";
					PrintVector<float>(s->tagHistory, tagsdat);
				}

				if (s->weightHistory.size())
				{
					weighthistorydat << s->source_nrn->input_id<< " ";
					PrintVector<float>(s->weightHistory, weighthistorydat);
				}
			}
		}
	}

	ofstream sppdat((dirname + "/spikesperpattern.dat").c_str());
	for (uint i =0; i < this->spikesPerPattern.size(); i++)
	{
		for (uint j =0; j < this->spikesPerPattern[i].size(); j++) sppdat << this->spikesPerPattern[i][j] << " ";
		sppdat << endl;
	}
}



void LANetwork::StimDynamics(int duration) // stimulation dynamics, with dt=1msec 
{
	int t = 0;
	bool spikeState[this->neurons.size()+1];
	int lastSpikeT[this->neurons.size()+1];

	fill_n(spikeState, this->neurons.size()+1, 0);
	fill_n(lastSpikeT, this->neurons.size()+1, 0);

	//FILE* outfile 	= fopen("./data/0/wdata.dat", "w");

	for (syn_iter si = this->synapses.begin(); si != this->synapses.end(); ++si)
	{
		(*si)->calcium = 0.0;
	}


	for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ++ni)
	{
		LANeuron* n = *ni;
		n->wadapt = 0.0;
		n->vspike =0.;
		n->vreset =0.;
		n->V =0.;
	}



	for (branch_iter bi=this->branches.begin(); bi != this->branches.end(); ++bi)
	{
		(*bi)->totcalc = 0.0;
		(*bi)->depol = 0.0;
		(*bi)->dspike = 0.0;
		(*bi)->dspikeT = -1;
	}

	for (t=0; t < duration; t++)
	{
		for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ++ni)
		{
			LANeuron* n = *ni;
			double s_inh = 0.0;
			double sv =0.0 ; //, ss =0;
			
			for (branch_iter bi=n->branches.begin(); bi != n->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				float sumepsp =0.;
				float sumipsp =0.;

				//if (this->debugMode &&  b->bid==DEBUG_BID) vfile<< b->depol << " " << n->V << endl;

				for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s = *si;
					if (spikeState[s->source_nrn->nid])
					{

						if (s->source_nrn->type == 'I') sumipsp += ( s->weight);
						else sumepsp += (  s->weight);
					}
				}


				b->depol +=  (4.0*sumepsp) -  b->depol/20.0;

				s_inh += sumipsp;

				//if (b->bid == 40) printf("BID %d f=%f epsp=%f\n", b->bid, b->depol, sumepsp);
				// http://www.cell.com/neuron/abstract/S0896-6273(12)00761-1: Inhibition blocks only weak dendritic spikes

				/*
				if ( b->dspikeT > 0)
				{
					float tds = t - b->dspikeT;
					if (tds < 40)
						b->depol = 40.;
					else if (tds < 60)
					{
					}
					else
					{
						b->dspikeT = -1; // Allow next dend spike
						b->depol =0.;// Reset dspike
					}
				}
				*/

			

				if (n->type == 'P' && b->dspikeT < t-100 && (n->vspike + b->depol) > 30.*dendSpikeThresh) // Generate a dendritic branch spike
				{
					//b->depol = 60.0;
					b->depol = 50;
					b->dspikeT =t;
					b->branch_spikes++;
					n->dend_spikes++;
					if (this->enablePlasticity && b->strengthTag < 1.0) b->strengthTag +=  (1.0 - b->strengthTag)*0.20;
				}

			//	if (b->dspikeT >0 && t-b->dspikeT < 20) b->depol = 50.;



				for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s = *si;
					if (spikeState[s->source_nrn->nid])
					{
						if (this->enablePlasticity && s->isPlastic) // && t > 50b)
						{
							float depol =  b->depol + n->vspike;
							//float depol =  n->vspike;
							//NMDA dynamics: Jahr_1993
							//if (depol > 10.0)
							if (depol > 1.0) //n->vspike > 10)
							{
								//float ff = exp(-dt/70.0 ) * (1.0/(1.+exp( (-(depol-30.0)/7.0))));
								float ff =  (1.0/(1.+exp( (-(depol-30.0)/5.0))));

								//s->calcium +=  ( 0.5*ff)*(1.0 - s->calcium) ;
								//s->calcium +=  2.0*(ff);
								//cout << " Calc " << s->calcium<< endl;
								s->calcium +=  ff/10.;
							}
						}
					}
				}




				sv +=  (b->depol)*b->strength ;

				//if (b->dspike > 0.0) b->dspike -= (b->dspike*b->dspike)/40.0; // see Losoczy 2008  for epsp durations
				//b->branchVoltageHistory.push_back(b->depol1 + n->vspike + b->dspike);
				//b->branchCalciumHistory.push_back(bcalc);
			}


			n->inh_cur += 0.18*s_inh; 

			if (n->type == 'S') // Source neuron
			{
				LAInput* in = (LAInput*)n;
				//if (t ==0) printf ("nid=%d , T=%d, next=%d\n", in->nid, t, in->nextSpikeT);
				if (in->spikeTimes && t >= in->nextSpikeT && in->nextSpikeT >0)
				{
					//printf("ID %d spike! t=%d\n", in->nid, t );
					if (in->curSpike < in->totalSpikes)
						in->nextSpikeT = in->spikeTimes[in->curSpike++];
					else 
						in->nextSpikeT = -1;
					spikeState[in->nid] = 1;
					n->isSpiking = true;
					lastSpikeT[in->nid] = t;
					in->total_spikes++;
					in->V = param_thresh;

				}
				else
				{
					spikeState[in->nid] = 0;
					n->isSpiking = false;
					in->V = param_E_L;
				}
			}
			else /// PYR OR IN NEURON
			{

				if (spikeState[n->nid])
				{
					//spikeState[n->nid] =0;
					n->V = 0.0;
					n->wadapt += 0.18;
				}
				
				n->exc_cur = sv*0.12;

				if (n->type == 'I') // Interneuron
	 			{
					n->V +=  (n->exc_cur - n->inh_cur) - (n->V)/20. - n->wadapt*(n->V+10.0);
					n->wadapt -= n->wadapt/70.;
				}
				else
				{
					n->V +=  n->exc_cur - 3.0*s_inh - (n->V)/30. -   n->wadapt*(n->V+10.0) ;

					if (this->disableCreb)
						n->wadapt -= n->wadapt/180.;
					else
						n->wadapt -= n->wadapt/((180. - 70.0*(n->crebLevel>0.2 ? 1. : 0.)));
				}


				if ( lastSpikeT[n->nid] < t-2 && n->V > (20.0 - (n->crebLevel>100.2 ? 2.0 : 0.) ))
				{
					spikeState[n->nid] = 1;
					lastSpikeT[n->nid] = t;
					n->total_spikes++;
					n->isSpiking = true;
					n->vspike = 30.0; 
					n->V = 70.;

				}
				else
				{
					spikeState[n->nid] = 0;
					n->isSpiking = false;
				}

				if (n->vspike > 0) n->vspike -= n->vspike / 17.0; 
				// Legenstein&maas eta_reset = 20msec
				

				if (n->wadapt <-0.) n->wadapt =0.;
				//if (n->wadapt < 0.) n->wadapt =0.;

			}

			if (spikeState[n->nid])
			{
				n->spikeTimings.push_back(t+ this->Tstimulation);
			}
		}

		#ifdef WXGLMODEL_H
		if (this->wx_window && t%10==0)
		{
			this->wx_window->UpdatePanel();
		}
		#endif
		

		//usleep(10000);
	}



	for (branch_iter bi=this->branches.begin(); bi != this->branches.end(); ++bi)
	{
		LABranch* b = *bi;
		b->branchSpikesHistory.push_back(b->branch_spikes);
	}


	//fflush(stdout);
	//fflush(spikesFile);

	this->Tstimulation += duration;

}


void LANetwork::SaveSnapshot(char* filename)
{
	ofstream of(filename);
	for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
	{
		LANeuron* nrn = *ni;
		if (nrn->type == 'P')
		{
			float totTag =0;
			float totProtein =0;
			for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s = *si;
					if (s->isPlastic)
					{
						totTag += s->stag;
						//of<<s->sid << " "<< s->source_nrn->input_id << " "<< b->bid << ' ' << nrn->nid <<' '<< b->protein  <<' '<< s->stag <<' '<< endl;
					}
				}
				totProtein += b->protein;
			}
			of << nrn->nid << ' ' << totProtein << ' '<<  totTag << endl;
		}
	}

}



void LANetwork::Interstim(int durationSecs)
{
	int tstop = T + durationSecs;
	this->isInterstim = true;

	printf("Interstim %d seconds (T=%d) plast=%d G=%d, L=%d ... \n", durationSecs, T, this->enablePlasticity, this->globalProteins, this->localProteins);

	float tstep = 60.0;
	int totalWeak =0;
	float weightLimit =  initWeight + 0.0;

	// Count the total weak afferent synapses
	for (input_iter ii = this->inputs_cs.begin(); ii != this->inputs_cs.end(); ++ii)
	{
		nrn_list lst = *ii;
		for (nrn_iter n = lst.begin(); n != lst.end(); ++n)
			for (syn_iter si=(*n)->outgoing.begin(); si != (*n)->outgoing.end(); ++si)
				if ((*si)->weight <= weightLimit)
					totalWeak++;
	}

	int trec =0, tactrec=0;
	int totTagged=0, totTaggedMinus=0, totBProd =0;
	float totLTP=0., totLTD=0.;
	int totact =0, totSact =0;
	int totbspikes =0, totSpikes=0;
	float maxSpk =0;
	
		float actmin=9999, actmax=0;
		float nactmin=9999, nactmax=0;

	for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
	{
		LANeuron*  nrn = *ni;
		float nrnCalc =0.0;
		if (nrn->total_spikes > 4)
			totSact++;


		for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
		{
			LABranch* b = *bi;
			totbspikes += b->branch_spikes;

			if (!this->enablePlasticity)
				continue;

			for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
			{
				LASynapse* s =*si;
				float ctag = caDP(s->calcium);

				if (fabs(s->stag) < 0.1) /// Avoid erasing existing tags
				{
					s->stag = ctag;
					
					if (s->stag > 0.1)
					{
						totTagged++;
						totLTP += s->stag;
						//if (s->stag > 0.5 && b->prpTransients.size()>0) cout <<"Candidate is "<<s->sid<< " " << b->bid << " "<< nrn->nid<<endl;
					}
					if (s->stag < -0.1)
					{
						totTaggedMinus++;
						totLTD += s->stag;
					}	
				}

				/* */
				//printf("SID=%d Calc=%f tag=%f %f\n", s->sid, s->calcium, ctag, s->stag);
				//s->stag = ctag ;
				//* (1.0 - s->stag) ; // * ((1.- 1. / (1.+ exp(-(totw-2.5)*1.5))));
				//else if (ctag < -0.05) s->stag =  ctag  ; // * (s->stag - (-1.0)) ;  // * ((1.- 1. / (1.+ exp(-(totw-2.5)*1.5))));
				//
				//if (s->stag > 1.5 && nrn->prpTransients.size()==1) printf("SID=%d tag=%f PR=%f\n", s->sid, s->stag, nrn->proteinRate);

				//if (nrn->nid  == 65) printf("ST=%f %d, pr=%f\n", s->stag, s->sid, b->proteinRate);
				//if (s->stag > 0.2) printf("ID=%d %f\n", s->sid, s->stag);



				//s->eltp = caDP(s->calcium);
				//if (s->calcium >0.0) printf("NID=%d sid=%d calc=%f tag=%f\n", nrn->nid, s->sid, s->calcium, s->stag);
				b->totcalc += s->calcium;
				s->calcium = 0.;
			}


			
			if ( b->totcalc  > this->localPRPThresh*BPROD_CUTOFF) // This branch should produce PRPs now BPROD
			{
				b->prpTransients.push_back( pair<float,float>(T, b->totcalc));
				//printf ("BID %d  NID %d BCALC=%f NCALC=%f S=%d\n", b->bid, nrn->nid, b->totcalc, nrnCalc, (int)b->prpTransients.size());
				totBProd++;
			}

			nrnCalc +=  b->totcalc;

		}


		if (nrn->total_spikes > CUTOFF*4.0)
		{
			totSpikes += nrn->total_spikes;
			totact++;
			//printf("ACT %d\t %.2f\t%d\t%.2f\n",  nrn->nid  , nrn->crebLevel  , (nrn->total_spikes) , nrnCalc);

			updminmax(actmin, actmax, nrnCalc);

		}
		else
			updminmax(nactmin, nactmax, nrnCalc);



		if (maxSpk < nrn->total_spikes)
			maxSpk = nrn->total_spikes;


		if (this->enablePlasticity)
		{
			if (nrnCalc > this->globalPRPThresh*GPROD_CUTOFF) /// GPROD
			{
				nrn->prpTransients.push_back( pair<float,float>(T, nrnCalc) );

				if (!this->disableCreb ) 
					nrn->crebLevel=1.0;

				if (nrn->total_spikes > CUTOFF*4)
					tactrec ++;
				trec ++;
				//printf("id %d nrnCalc=%f spk=%d, dspk=%d CREB=%f trans=%d\n", nrn->nid, nrnCalc, nrn->total_spikes, nrn->dend_spikes, nrn->crebLevel,(int)nrn->prpTransients.size());
			}
		}



		nrn->totcalc  = nrnCalc;

	}



	printf("\n\nAct=[%f,%f] nonact=[%f,%f] \n", actmin, actmax, nactmin, nactmax);

	
	if (this->runningMode == RUN_TRAINING)
	{
		char buf[256];
		sprintf(buf,  "./data/r%d.dat", this->runningPatternNo);
		//this->SaveSnapshot(buf);
	}
	
	printf("TAGGED: %d/%d +%.1f/%.1f GPROD:%d (%.1f%%) BPROD:%d, ACT:%d (%.1f%% ), act+gprod:%d FREQ:%.1f max %.1f DSPK:%d sact:%d\n", totTagged, totTaggedMinus, totLTP, totLTD, trec, 100.*(float)trec/(float(this->pyr_list.size())), totBProd, totact, 100.*(float)totact/(float(this->pyr_list.size())), tactrec, (float)totSpikes/((float)this->pyr_list.size()*4.0), float(maxSpk)/4.0, totbspikes, totSact);

	for (; T < tstop; T+= tstep)
	{
		for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
		{
			LANeuron*  nrn = *ni;

			float totalSynapseWeight =0.0;
			float totalBranchStrength =0.0;
			int totalSynapses =0;
			nrn->protein =0.0;

			// "Whole-neuron" distribution of proteins
			nrn->proteinRate =0;
			for (pair_iter ii = nrn->prpTransients.begin(); ii != nrn->prpTransients.end(); ++ii)
			{
				pair<float, float> p = *ii;
				int td= (T - p.first);
				float al = (nuclearproteinalpha(td));
				if (nrn->proteinRate < al)
					nrn->proteinRate =  al;
			}

			
			
			for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				b->proteinRate =0.;

				for (pair_iter ii = b->prpTransients.begin();ii != b->prpTransients.end(); ++ii)
				{
					pair<float, float> p = *ii;
					float td = float(T - p.first);
					float al = (branchproteinalpha(td));
					if (b->proteinRate < al)
						b->proteinRate = al;
				}


				float  f =0.;

				if (this->localProteins)
					f = 1.0*b->proteinRate; 
				else if (this->globalProteins)
					f =  1.0* nrn->proteinRate;
				else
				{
					f = 1.0*b->proteinRate + 1.0* nrn->proteinRate;
					if (f>1.0) f = 1.0;
				}



				b->protein = f; // += tstep*(f/2000. - b->protein/3600.);


				vector<LASynapse*> candidates;


				for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s =*si;

					if (s->stag != 0.0) 
					{
						s->stag -= (tstep/3600.)* s->stag;
						//s->stag -= (s->stag>0. ? +1. : -1.)*tstep / (3.4*3600);


						//if (s->sid == 3024) printf("STAG=%f\n", s->stag);

						if (b->protein > 0.1 && (s->stag >0.1 || s->stag < 0.1))
						{
							//candidates.push_back(s);
							float fw = s->stag* b->protein;
							//if (s->stag >0.)
							s->weight += tstep * fw/400.;
						}
					}

					if(1)
					{
						if (s->weight  > maxWeight)
						{
							s->weight = maxWeight;
							//s->stag =0.;
						}
						else if (s->weight  < 0.)
						{
							s->weight = 0.;
							//s->stag =0.;
						}
					}


					totalSynapseWeight += s->weight;
					totalSynapses++;

					//Homeostasis
			
					//if (1) //this->runningMode  != RUN_TRAINING)
					//if (this->runningMode  == RUN_TEST)
					if (1)
					{
						s->weight += s->weight * (1.0 - nrn->synScaling )*tstep/(7.*24.0* 3600.*homeostasisTimeParam); // Synaptic scaling (Ref: http://www.sciencedirect.com/science/article/pii/S0896627308002134)
					}
					//if (s->weight > 1.2) printf("SID =  %d %f\n", s->sid, s->weight);


					if (this->debugMode && b->bid == DEBUG_BID)
					{
						//cout << b->proteinRate << " " << nrn->proteinRate << " " << b->protein << " " << s->stag << " " << s->weight << " " << endl;
						//this->mfile << " "  <<   s->stag << " " << s->weight ;
					}
				}

				/*
				if (candidates.size()>0 && b->protein >0.1)
				{
					LASynapse* s = candidates.at(int(rgen()*(float)candidates.size()));

					// Consolidation
					float q = b->protein* (s->stag>0 ? 1.0 : -1.0) *tstep/600.0;
					float dw;

					if (q>=0)
						dw =  q; //(MAXSTRENGTH - s->weight);
					else
						dw = q; // *(s->weight - 0.0);
					s->weight += dw;
					s->stag -= q*1.0;

					totCons += dw;

					if (s->weight  > MAXSTRENGTH)
					{
						s->weight = MAXSTRENGTH;
						//s->stag =0.;
					}
					else if (s->weight  < 0.)
					{
						s->weight = 0.;
						//s->stag =0.;
					}
				}
				*/



				//if (this->debugMode && b->bid == DEBUG_BID) this->mfile << endl;



				if (0)  // Branch strength potentiation 
				{
					if (b->strengthTag>0.01)
					{
						b->strength += tstep*b->strengthTag*(1.5- b->strength) / (2.0*60.*60.*BSPTimeParam); // Max branch strength=2, time to get there: 40 minutes losonczy 2008 figure 3 http://www.nature.com/nature/journal/v452/n7186/abs/nature06725.html
						b->strengthTag -= tstep*b->strengthTag / (10.*60.);
					}
					totalBranchStrength += b->strength;
					b->strength += tstep*b->strength*(1.0 - nrn->branch_scaling) / (3.0*3600.); // Branch homeostatic scaling 
				}


				if (T%800 ==0)
				{
					b->branchProteinHistory.push_back(b->protein);
					//b->branchStrengthHistory .push_back(b->strength);
				}
			}

			//if (T%500==0) nrn->proteinHistory.push_back(nrn->proteinRate);

			// Synaptic homeostasis / synaptic scaling
			//if (nrn->synapticWeightsInitialSum != 0.0 )
			if (totalSynapses>0)
				nrn->synScaling = totalSynapseWeight / (initWeight*float(totalSynapses));
		
			else
				nrn->synScaling = 1.0;

			// Branch plasticity homeostasis
			nrn->branch_scaling = totalBranchStrength/((float(1.0) * float(nrn->branches.size())));


			//creb Drop
			if (nrn->crebLevel >0.0)
			{
				nrn->crebLevel -= tstep/(3600.*8.*CREBTimeParam );
				//if (nrn->nid == DEBUG_NID) printf("CREBlev=%f\n", nrn->crebLevel);
				//nrn->crebLevel -= nrn->crebLevel * tstep / (3600.0*4.);
			}


			if ( this->debugMode &&  nrn->nid == DEBUG_NID) this->mfile << nrn->V <<endl;
		
		}


		#ifdef WXGLMODEL_H
		if (this->wx_window && T%20 ==0)
		{
			this->wx_window->UpdatePanel();
		}
		#endif
	}
	


	//printf("CONSOLIDATED: %f\n", totCons);
	

	//this->mfile.flush();

	// zuo et al. 2005
	// We assume 100% of  filopodia turn over over 24 hours. 
	// We assume that weak synapses (W < 0.3) are filopodia
	
	if (this->enablePruning)
	if (this->enablePlasticity)
	{
		
		this->isPruningSynapses = true;
		int synapsesToPrune = 0.005*float(totalWeak) * float(durationSecs)/3600.0;
		cout << "removing "<< synapsesToPrune << "synapses" << endl;
		int found = this->PurgeInputSynapses(synapsesToPrune, weightLimit);
		if (found > 0)
			this->CreateInputSynapses(found);
		cout << "Renewed " << found <<" synapses"  << endl;
		this->isPruningSynapses = false;

	}
	
	this->isInterstim = false;
}




void LANetwork::Begin()
{
	
}




