#include "constructs.h"
#include "wxmodel.h"
//#include <gluplot.h>
//using namespace glp;

#define CIMG  0
#define STDP 0



#if  CIMG
#include "CImg.h"
using namespace cimg_library;
#endif

//static FILE*  eltpfile = 0;
static FILE* eltpfile =  NULL;
MTRand rgen(1980); 



inline float sigmoid(float x, float x0, float broad)
{
	return (1.0 / (1.0 + exp(-x - x0)/broad));
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

inline float branch_nonlinearity(float V)
{
	if (V < 18.0) return V;
	else return (18.0+ 30.0 / (1.0 + exp(-(V-22.0))));
	//if (V <50.0) return V; 
	//else return 50.0;

}


inline float calcium_step(float post_depolarization)
{
	//return 1.0/(1.+exp( (-(post_depolarization-10.0))*2.0));
	return 1.0/(1.+exp( (-(post_depolarization-23.0))*1.0));
}


inline float caDP(float ca)
{
	return (1.4/(1.+exp(-(ca*10.-5.0)*10.))) - (0.4/(1+exp(-(ca*10.0-4.0)*10.0)));
}


void LANetwork::CreateNeurons(int number, char type, vector<LANeuron*>* appendTo = 0, int inputId =0)
{
	for (int i =0 ;i < number; i++)
	{
		LANeuron* n ;
		if (type == 'S')
			n = new LAInput();
		else
		{
			n = new LANeuron();

			for (int tt =0; tt < N_BRANCHES_PER_NEURON; tt++)
			{
				LABranch* bb  = new LABranch;
				bb->bid = this->branches.size();
				bb->neuron = n;

				this->branches.push_back(bb);
				n->branches.push_back(bb);
			}
		}
		n->type = type;
		n->input_id = inputId;
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


int LANetwork::ConnectByDistance(vector<LANeuron*> fromList, vector<LANeuron*> toList, float fromDistance, float toDistance, int nNeuronPairs, int nSynapsesPerNeuron, bool isPlastic= false)
{
	int tpairs =0;
	while(true)
	{
		LANeuron* a = fromList.at(int(rgen()*(float)fromList.size()));
		LANeuron* b = toList.at(int(rgen()*(float)toList.size()));

		if (distances[pair<int,int>(a->nid, b->nid)] >= fromDistance && distances[pair<int,int>(a->nid, b->nid)] <= toDistance) 
		{
			for (int i =0; i < nSynapsesPerNeuron; i++)
			{
				LASynapse* syn = new LASynapse();
				syn->sid  = this->synapses.size();
				syn->source_nrn = a;
				syn->target_nrn = b;
				syn->isPlastic = isPlastic;
				syn->target_branch = syn->target_nrn->branches[rgen()*float(syn->target_nrn->branches.size())];

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


void LANetwork::Create()
{
	int n_neurons = N_NEURONS;
	this->nNeurons = n_neurons;

	CreateNeurons(n_neurons*0.8, 'P', &this->pyr_list);
	CreateNeurons(n_neurons*0.2, 'I', &this->in_list);
	

	int n_inputs = N_INPUTS ;
	int n_per_input = 6;
	for (int i=0;i < n_inputs;i++)
	{
		vector<LANeuron*> nrnsCS, nrnsUS;

		CreateNeurons(n_per_input, 'S', &nrnsCS, i);
		CreateNeurons(n_per_input, 'S', &nrnsUS, i);


		this->inputs_cs.push_back(nrnsCS);
		this->inputs_us.push_back(nrnsUS);

	}

	this->CalculateDistances();


	if (0) { 
		ConnectByDistance(this->pyr_list, this->pyr_list, 0.0, 3.5, 0.02*float(this->pyr_list.size()*this->pyr_list.size()), 1, false);
		ConnectByDistance(this->pyr_list, this->in_list, 0.0, 3.6, 0.02*float(this->pyr_list.size()*this->in_list.size()), 1, false);
		ConnectByDistance(this->in_list, this->pyr_list, 0.0, 3.6, 0.02*float(this->pyr_list.size()*this->in_list.size()), 1, false);
	}



	int n_per_neuron = 1;

	for (unsigned int i =0; i < this->inputs_cs.size(); i++)
	{
		ConnectByDistance(this->inputs_cs[i], this->pyr_list, 0.0, 3.0, 1.0*float(n_per_input*this->pyr_list.size()), n_per_neuron, true);
		//ConnectByDistance(this->inputs_cs[i], this->in_list, 0.0, 3.0, 0.4*float(n_per_input*this->in_list.size()), n_per_neuron, false);

		ConnectByDistance(this->inputs_us[i], this->pyr_list, 0.0, 3.0, 1.0*float(n_per_input*this->pyr_list.size()), n_per_neuron, false);
		//ConnectByDistance(this->inputs_us[i], this->in_list, 0.0, 3.0, 0.4*float(n_per_input*this->in_list.size()), n_per_neuron, false);

	}

	eltpfile = fopen("/tmp/eltp.dat", "w");

}



void LANetwork::StimDynamics(int duration) // stimulation dynamics, with dt=1msec 
{
	int t = 0;
	bool spikeState[this->neurons.size()+1];
	int lastSpikeT[this->neurons.size()+1];
	fill_n(spikeState, this->neurons.size()+1, 0);
	fill_n(lastSpikeT, this->neurons.size()+1, 0);

	this->voltageData = new Arr2D(duration, this->neurons.size());



#if CIMG
	CImg<float> vimg(duration, this->neurons.size());
	CImg<float> synimg(duration, this->neurons.size());
	vimg.fill(0);
	synimg.fill(0);
#endif

	FILE* outfile = fopen("/tmp/wdata.dat", "w");
	FILE* spfile = fopen("/tmp/spdata.dat", "w");

	for (t=0; t < duration; t++)
	{
		for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ++ni)
		{
			LANeuron* n = *ni;
			double sv =0.0;

			for (branch_iter bi=n->branches.begin(); bi != n->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				if (spikeState[n->nid]) // Back propagating action potential
					b->depol =  40.0;

				for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s = *si;
					if (spikeState[s->source_nrn->nid])
					{
						b->depol += (s->iltp + s->eltp + s->weight) * ((s->source_nrn->type == 'I') ? -3.0 : 3.0); // XXX
					}
				}


				b->depol = branch_nonlinearity(b->depol);

				for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s = *si;
					if (!s->isPlastic) continue;
#if STDP
					if (spikeState[src->nid])
					{
						if (t > 50 && !spikeState[s->target_nrn->nid] && s->isPlastic)
						{
							stdp -= 0.9*exp(-(float(t - lastSpikeT[s->target_nrn->nid])) / 30.0);
						}
					}

					if (spikeState[s->target_nrn->nid] && t > 10 && s->isPlastic)
					{
						stdp += exp(-float(t-lastSpikeT[src->nid])/30.0);
					}

					stdp *= 0.1;
					if (stdp != 0.0  && this->enablePlasticity && s->isPlastic)
					{
						adjust_with_bounds(s->stag, stdp, 1.0, -1.0);
					}

					//if (stdp != 0.0) if (s->sid == 603) printf("%s %f  %f %d lpost=%d lpre=%d\n", spikeState[s->target_nrn->nid] ? "post" : "pre", s->stag, stdp, t, lastSpikeT[s->target_nrn->nid] , lastSpikeT[src->nid]);



#else
					if (spikeState[s->source_nrn->nid])
					{
						float ff = (1.0-s->calcium)*0.1*calcium_step(b->depol);
						s->calcium +=  ff;
						//if (s->sid == 603) printf("tag=%f calc=%f depol=%f ff=%f \n",  s->stag, s->calcium, b->depol , ff);
					}

					if (s->calcium > 0)
						s->calcium -= s->calcium / param_tau_calcium_syn;

					if (t > 50 &&  t%100 ==0)
					{
						float dir  = caDP(s->calcium)/10.0;
						adjust_with_bounds(s->stag, dir, 1.0, -1.0);
					}
#endif

					if (s->source_nrn->nid == 60) fprintf(spfile, " %f", s->stag);
				}

				sv += b->depol*(1.0 / float(N_BRANCHES_PER_NEURON));

				if (b->depol != 0) b->depol -= b->depol/param_tau_branch;
			}

			if (n->type == 'S')
			{
				LAInput* in = (LAInput*)n;
				if (t >= in->nextSpikeT && in->nextSpikeT >0)
				{
					if (in->curSpike < in->totalSpikes)
						in->nextSpikeT = in->spikeTimes[in->curSpike++];
					else 
						in->nextSpikeT = -1;
					spikeState[in->nid] = 1;
					lastSpikeT[in->nid] = t;
					in->total_spikes++;
					in->V = 0.0;

					//if (in->nid == 60) printf ("T=%d, next=%d\n", t, in->nextSpikeT);
				}
				else
				{
					spikeState[in->nid] = 0;
					in->V = param_E_L;
				}
			}
			else
			{
				if (spikeState[n->nid])
				{
					n->wadapt += param_beta;
					spikeState[n->nid] =0;
					n->V = param_E_L;
				}

				n->wadapt += (param_alpha*(n->V-param_E_L) - n->wadapt) / (n->type == 'I' ? param_tau_wadapt_in : (param_tau_wadapt_pyr*(1.0 - n->crebLevel*0.3)));


				float ww = ((param_g_L * (param_E_L - n->V))+(param_g_L*param_delta_t*exp((n->V-param_thresh)/param_delta_t)) - n->wadapt + 14.0*sv )/param_C;

				//if (n->nid == 60) printf("SV=%f ww=%f\n", sv, ww);

				n->V += ww; // + rgen()*2.0;

				//if (n->V > (param_thresh - n->crebLevel*10.0))
				if (n->V > (param_thresh))
				{
					spikeState[n->nid] = 1;
					lastSpikeT[n->nid] = t;
					n->total_spikes++;
					n->V = 0.0;
				}
				else
				{
					spikeState[n->nid] = 0;
				}
			}

			fprintf(outfile, " %f ", n->V);
#if CIMG
			vimg.atXY(t, n->nid) =  n->V;
#endif
		}

		fprintf(outfile, "\n");
		fprintf(spfile, "\n");

		if (this->wx_panel && t%200==0)
		{
			this->wx_panel->Refresh();
			this->wx_panel->Update();
		}
	}

	
	for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ++ni)
	{
		LANeuron* n = *ni;

		float ltp[2] = {0};
		float oltp[2] = {0};

		for (syn_iter si = n->incoming.begin(); si != n->incoming.end(); ++si)
		{
			LASynapse* so = *si; 
			ltp[so->stag>0 ? 0 : 1] += fabs(so->stag);
		}

		for (syn_iter si = n->outgoing.begin(); si != n->outgoing.end(); ++si)
		{
			LASynapse* so = *si; 
			oltp[so->stag>0 ? 0 : 1] += fabs(so->stag);
		}


		printf("#%d (%c) Spikes : %d  LTPIN: +%f -%f,  LTPOUT: +%f -%f \n", (*ni)->nid, n->type, (*ni)->total_spikes, ltp[0],ltp[1], oltp[0], oltp[1]);

	}


	//gplot << LINES << vsurf;

	fclose(outfile);
	fclose(spfile);

#if CIMG
	//synimg.save("synplot.png");
	CImgDisplay window(vimg);
	while (!window.is_closed()) window.wait(10);
#endif

}



void LANetwork::Interstim(int duration)
{
	int tstop = T + duration;
	
	printf("Interstim %d seconds (T=%d) ... \n", duration, T);
	float tstep = 50.0;
	for (; T < tstop; T+= tstep)
	{
		for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ni++)
		{
			LANeuron*  nrn = *ni;
			float totalSynapseWeight =0.0;
			int totalSynapses =0;

			for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s =*si;
					if (s->stag !=0.0)
					{
						if (s->stag >0.4)
							s->eltp += tstep*(s->stag)*(2.0 - (s->weight + s->iltp + s->eltp)) / param_tau_eltprise;
						else if (s->stag <-0.4) s->eltp -= tstep*fabs(s->stag)*((s->iltp + s->weight + s->eltp) - 0.05) / param_tau_eltprise;
						if (s->stag >= 0.1 || s->stag <= -0.1)
						{
							float val = fabs(s->stag);
							val = (val / (0.5+val));
							b->proteinRate += 0.7*tstep*  val /1000.0;
							nrn->proteinRate += 0.25*tstep* val/1000.0;
						}
						s->stag -= tstep* s->stag / param_tau_stag;
					}

					if (s->trigger > 0 )
						s->trigger -= tstep*s->trigger / param_tau_trigger;


					float tprotein = b->protein+nrn->protein;

					if (tprotein > 0.1 && s->eltp != 0)
					{
						float q = tprotein*tstep* s->eltp / 1000.0;
						s->weight  += q;
						s->eltp  -= q;
						//s->eltp   -= q;
					}

					//if (s->eltp != 0.0) s->eltp -= tstep* (s->eltp) / param_tau_eltp;
					if (s->eltp != 0.0) s->eltp -= tstep* s->eltp/param_tau_eltp;

					///if (T%10==0 && s->sid == 530) printf("%d %d (%d->%d) stag=%f eltp=%f iltp=%f weight=%f prateBranch=%f pbranch=%f prateSoma=%f pSoma=%f \n", T, b->bid, s->source_nrn->nid, s->target_nrn->nid, s->stag, s->eltp, s->iltp, s->weight, b->proteinRate, b->protein, nrn->proteinRate, nrn->protein);
					//if (s->target_nrn->nid == 49) printf("SID=%d\n", s->sid);

					if (s->sid == 236) fprintf(eltpfile, "%f %f %f %f %f %f %f %f %f\n", s->stag, s->eltp, b->proteinRate, b->protein, nrn->proteinRate, nrn->protein, s->iltp, s->weight, s->trigger );

					totalSynapseWeight += s->weight;
					totalSynapses++;
					s->weight += tstep*s->weight*(1.0 - nrn->synScaling) / param_tau_scaling;
				}

				if (b->proteinRate > 0.1)
				{
					float v = b->proteinRate*b->proteinRate;
					b->protein += tstep* (1.0-b->protein)*(v/(.2+v))/ (900.0);
					//b->protein += tstep* (1.0 - b->protein) * (b->protein / (1.0 + b->protein));
				}

				if (b->proteinRate > 0.0) b->proteinRate -= tstep*b->proteinRate / 3600.0;
				//param_tau_proteinrate;

				if (b->protein > 0.0) 
				{
					b->protein -= tstep*b->protein/param_tau_protein_branch;
					//param_tau_protein_branch;
					//if (b->protein <0.0) b->protein = 0.0;
				}
			}

			if (nrn->proteinRate > 0.1)
			{
				//nrn->crebLevel = 0.9;
				float v = nrn->proteinRate*nrn->proteinRate;
				nrn->protein += tstep* (1.0-nrn->protein)* (v/(0.2+v))/ (900.0);
			}
			else 
				nrn->crebLevel =0.0;

			if (nrn->proteinRate > 0.0) nrn->proteinRate -= tstep*nrn->proteinRate / (param_tau_proteinrate);

			if (nrn->protein > 0.0) 
				nrn->protein -= tstep*nrn->protein / param_tau_protein_neuron;

			if (totalSynapses != 0.0)
				nrn->synScaling = totalSynapseWeight/float(totalSynapses);
			else
				nrn->synScaling = 1.0;
			
		}

		if (this->wx_panel)
		{
			this->wx_panel->Refresh();
			this->wx_panel->Update();
		}
	}
	fflush(eltpfile);
}



void LANetwork::Stimulate2(int nInput, int mode)
{
	this->Stimulate(nInput, 1);
	if (mode)
	{
		this->Interstim(8*60);
		this->Stimulate(nInput, 1);
	}
}


void LANetwork::Stimulate(int nInput, int mode)
{
	printf("Stimulating input #%d (%s)... \n", nInput, (mode ? "CS+US" : "CS"));

	int dur = 2000;
	for (nrn_iter n = this->inputs_cs[nInput].begin(); n != this->inputs_cs[nInput].end(); ++n)
	{
		LAInput* in = (LAInput*)*n;
		in->reset();
		in->program(200, dur, 60, 0.1);
	}

	for (nrn_iter n = this->inputs_us[nInput].begin(); n != this->inputs_us[nInput].end(); ++n)
	{
		LAInput* in = (LAInput*)*n;
		in->reset();
		if (mode ==1)
			in->program(200, dur, 60, 0.4);
	}

	this->StimDynamics(dur + 600);
	printf("Done.\n");
}


void LANetwork::Begin()
{
	
}


void LANetwork::Cleanup(void)
{
	fclose(eltpfile);
}


