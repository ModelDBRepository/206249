//
//  Version: $Id: wxmodel.cpp 162 2014-01-08 15:09:14Z gk $
//  Implementation of GUI interface
//
/* 
    Version: $Id: constructs.h 170 2014-01-29 14:00:04Z gk $
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
#include "wxmodel.h"


enum {
	MENU_RUN=1,
	MENU_1Hour, 
	MENU_REFRESH,
	MENU_QUIT,
	MENU_CS = 100,
	MENU_CSUS = 300,
	MENU_STRONGTRAIN = 400,
	MENU_INTERSTIM = 500,
	MENU_PROTO = 600,
};



LAWindow::LAWindow(const wxString& title, LANetwork* net) 
	:wxFrame(NULL, wxID_ANY, title, wxDefaultPosition,wxSize(1000, 600))
{
	CreateStatusBar();
	SetStatusText(wxT("Ready"));
	Centre(); 

	this->network = new LANetwork();
	//this->network_panel = new LAPanel(this);
	this->network_panel->network = this->network;
	this->network->wx_window = this;

	//this->network_panel->SetFocus();

	wxMenu* m = new wxMenu;
	m->Append(MENU_RUN, _("&Run"));
	m->Append(MENU_REFRESH, _("Re&fresh"));
	m->Append(MENU_QUIT, _("&Quit"));

	wxMenuBar* mb = new wxMenuBar;
	mb->Append(m,  _("&Simulation"));

	m = new wxMenu;
	for (int i =0; i < network->n_inputs; i++)
		m->Append(MENU_CS+i, wxString::Format(wxT("CS %d"), i));
	mb->Append(m,  _("Run CS"));

	m = new wxMenu;
	for (int i =0; i < network->n_inputs; i++)
		m->Append(MENU_CSUS+i, wxString::Format(wxT("CS+US %d"), i));

	mb->Append(m,  _("Train"));

	m = new wxMenu;
	for (int i =1; i < 5; i++)
		m->Append(MENU_INTERSTIM+i, wxString::Format(wxT("Interstim %d hours"), i));
	mb->Append(m,  _("Interstim"));

	wxString  mprotos[] = {
		wxT("Coallocation Protocol"),
		wxT("Capacity Protocol"), 
	};

	m = new wxMenu;
	for (size_t i =0; i < 2; i++)
	{
		m->Append(MENU_PROTO+i, mprotos[i]);
	}
	mb->Append(m,  _("Protocols"));

	SetMenuBar(mb);
	//Refresh();
}



void LAWindow::OnQuit(wxCommandEvent& ev)
{
	printf("cleaning up... \n");
	this->network->Cleanup();
	Close(true);
	printf("Done\n");
}



/* Menu handlers */
void LAWindow::OnRun(wxCommandEvent& ev)
{
	int id = ev.GetId();
	printf("Running %d!\n", id);

	SetStatusText(wxString::Format(wxT("Wait, running %d..."), id));

	if (id >= MENU_PROTO) // Protocols
	{
		int n = id - MENU_PROTO;
		switch (n)
		{
			case 0:
			{
				this->network->CreateFearNet(100,20, 10, 4);
				for (int i =0;   i < network->n_inputs; i++)
				{
					SetStatusText(wxString::Format(wxT("Training # %d..."), i));
					network->Stimulate(i, 2);
					int mins = 50;
					SetStatusText(wxString::Format(wxT("Waiting %d minutes ..."), mins) );
					network->Interstim(mins*60);
				}
				SetStatusText(wxString::Format(wxT("Waiting 8 Hours...")));
				network->Interstim((60*8)*60);

				FILE* f = fopen("./data/0/memorytest.dat", "w");
				network->enablePlasticity = false;
				network->ResetCrebLevels();

				for (int i =0; i < network->n_inputs; i++)
				{
					SetStatusText(wxString::Format(wxT("Testing # %d..."), i));
					network->Stimulate(i, 1);
					for (nrn_iter ni = network->pyr_list.begin(); ni != network->pyr_list.end(); ++ni)
					{
						fprintf(f ,  "%d  ", (*ni)->total_spikes);
					}
					fprintf(f, "\n");
				}
				

				network->enablePlasticity = true;
				fclose(f);

				f = fopen("./data/0/synweights.txt", "w") ;
				for (syn_iter si = network->synapses.begin(); si != network->synapses.end(); ++si)
				{
					LASynapse* syn = *si;
					if (syn->isPlastic != -1)
						fprintf(f, "%d %d %d %d %f\n", syn->sid, syn->target_branch->bid, syn->target_nrn->nid, syn->source_nrn->input_id, syn->weight); 
				}
				fclose(f);
			}
			break;

			case 1:
			{
				this->network->CreateFearNet(100, 20, 10, 6);
				this->network->RunStoreTest(10, 1, 60, 0);
				printf("Done\n");
				//this->network->StoreDataFiles("./data/0/");
			}		
			break;

			case 2:

			break;

			case 3:
			break;

			case 4:
			break;
		}
	}
	else if (id >= MENU_INTERSTIM)
		network->Interstim((id - MENU_INTERSTIM)*3600);

	else if (id >= MENU_CSUS)
		network->Stimulate( id - MENU_CSUS, 2);

	else if (id >= MENU_CS)
		network->Stimulate( id - MENU_CS, 1);

	SetStatusText(wxT("Done!"));
	printf("Done running\n");
}



void LAPanel::OnTimer(wxCommandEvent& ev)
{

}


void LAPanel::OnClick(wxMouseEvent& e)
{

	/* Find out the nearest clicked neuron, and print properties in the console  */
	for (pt_iterator pi = nrnpoints.begin(); pi != nrnpoints.end(); ++pi)
	{

		wxPoint p = pi->second;

		if ( e.m_x >=p.x && e.m_x <= p.x + 80 && e.m_y <= p.y && e.m_y >= p.y-80)
		{
			int nid = pi->first;
			LANeuron* nrn = network->neurons[nid];
			printf("NID %d (%s) p=%f pr=%f\n", nid, nrn->type == 'I'? "in" :"pyr" , nrn->protein, nrn->proteinRate);

			for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				printf(" branch %d str=%f (tag %f) pb=%f pbr=%f depol=%f\n", b->bid, b->strength, b->strengthTag, b->protein, b->proteinRate, b->depol);
				for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse*  s = *si;
					printf("  syn %d from=%d (%d) calc=%f stag=%f eltp=%f w=%f\n", s->sid, s->source_nrn->input_id, s->source_nrn->nid, s->calcium, s->stag, s->eltp, s->weight);
				}

			}


		}

	}
	fflush(stdout);
}



LAPanel::LAPanel(wxFrame* parent)
	: wxPanel(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxBORDER_NONE | wxFULL_REPAINT_ON_RESIZE)
{


	Connect(wxEVT_PAINT, wxPaintEventHandler(LAPanel::OnPaint));
	Connect(wxEVT_LEFT_UP, wxMouseEventHandler(LAPanel::OnClick));


	hasStarted =0;
}



inline static void colormap1(wxBrush& pen, float value )
{
	pen.SetColour(250, 250-value*200, value*200); 
}


inline static void colormap1(wxPen& pen, float value )
{
	pen.SetColour(value*240.0, 0, 250 - value*200.0);
}


/*
static int palette[][3] = {
	{4, 3,200},
	{4, 200,3},
	{220, 200,1},
	{200, 190,3},
	{255,255,0},
	{255,0,255 },
	{255,0,0 },
	{128,0,0 },
	{192,192,192 },
	{128,128,128 },
	{250, 3,3},
	{0,255,0},
	{128,128,0},
	{200, 190,3},
	{0,128,0},
	{128,0,128},
	{0,255,255},
	{0,128,128},
	{ 50,50,255 },
	{4, 200,3},
};
*/


static wxPoint somap[7] = {
	wxPoint(30, 10),
	wxPoint(30, 10),
	wxPoint(40, 10),
	wxPoint(40, 10),
	wxPoint(40, 10),
	wxPoint(35, 1),
	wxPoint(30, 10)
};



static wxColour pyrcolor(150, 200, 150);
static wxColour incolor(200, 100, 100);
static wxColour syncolor(30, 200, 20);
static wxColour eltpcolor(200, 30, 20);
static wxColour calccolor(20, 30, 200);

static wxColour firecolor(255, 30,30);
static wxColour branchfirecolor(250, 100,30);
static wxPen pen;

void LAPanel::drawPyrCell(wxPaintDC& dc, LANeuron* n, int x , int y)
{
	int nrnheight = this->network->n_branches_per_neuron* 16;
	dc.SetFont(*wxSMALL_FONT);

	pen.SetColour(0,0,0);
	pen.SetStyle(wxCAP_PROJECTING);

	float v = (n->V - (-param_E_L)) / 20.0;
	if (v > 1.0) v = 1.0;
	else if (v  < 0.0) v = 0.0;

	wxColor* nrncolor = &pyrcolor;
	if (n->V > 35.0)
		nrncolor = &firecolor;

	//pen.SetColour(0,0,0);
	dc.SetBrush(wxBrush(*nrncolor, wxSOLID));
	pen.SetColour(*nrncolor);

	pen.SetWidth(3);
	dc.SetPen(pen);

	wxPoint px1(x+60, y);
	wxPoint px2(x+60, y-nrnheight-10);
	dc.DrawLine(px1, px2);

	pen.SetWidth(2);

	dc.SetPen(pen);
	dc.DrawPolygon(7, somap, x+25, y);
	dc.DrawText(wxString::Format(wxT("%f %d  %f"), n->exc_cur, n->total_spikes, n->totcalc), x+50, y+14);

	float yy = y-10;
	int nbr =0;

	for (branch_iter bi = n->branches.begin(); bi != n->branches.end(); ++bi)
	{
		//nbr = 0; // XXX
		LABranch* b = *bi;
		yy -= 15;
		wxPoint pt1(x+60, yy);

		int dx = x+60+56;
		wxPoint pt2(dx, yy);

		if (nbr%2)
		{
			pt2.x = x;
		}

		float hh = b->strength*3.0;
		pen.SetWidth(hh);
		if (b->dspike > 10.0)
			pen.SetColour( branchfirecolor);
		else 
			pen.SetColour( pyrcolor);

		dc.SetPen(pen);
		dc.DrawLine(pt2, pt1);

		dc.DrawText(wxString::Format(wxT("%d  %.3f "),  b->branch_spikes, b->totcalc ), x+110, yy);

		int xxx = x+66;
		if (nbr%2) xxx = x+54;
		int lastinput = -1;

		for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
		{
			LASynapse*  s = *si;

			if (s->source_nrn->isSpiking)
			{
				pen.SetColour(firecolor);
			}
			else if (s->source_nrn->input_type == TYPE_CS && s->source_nrn->input_id >=0)
			{
				pen.SetColour(syncolor);
				//pen.SetColour( palette[s->source_nrn->input_id][0], palette[s->source_nrn->input_id][1], palette[s->source_nrn->input_id][2]);
			}
			else
			{
				//continue;
				pen.SetColour(220, 220, 220);
				//pen.SetColour(255.*(s->stag+1.0) / 2.0, 255.*(s->eltp+1.0)/2.0, 200. );
			}


			int space =1;
			if (lastinput == -1) lastinput = s->source_nrn->input_id;
			else if ( s->source_nrn->input_id != lastinput)
			{
				space = 4;
				lastinput = s->source_nrn->input_id;
			}

			if ( nbr%2)
				xxx -= space;
			else
				xxx += space;


			wxPoint sp1(xxx, yy-1);

			int ddd = yy -1 - ( s->weight)*10.0;
			int ddd2 = ddd - (s->stag)*10.0;

			wxPoint sp2(xxx, ddd-1);
			wxPoint sp3(xxx, ddd2-1);
			wxPoint sp4(xxx, ddd2-1 - s->calcium*10.);

			pen.SetWidth(1);
			dc.SetPen(pen);
			dc.DrawLine(sp1, sp2);

			pen.SetColour(eltpcolor);
			dc.SetPen(pen);
			dc.DrawLine(sp2, sp3);

			//pen.SetColour(calccolor);
			//dc.SetPen(pen);
			//dc.DrawLine(sp3, sp4);

		}
		nbr++;
	}
}


void LAPanel::drawInCell(wxPaintDC& dc, LANeuron* n, int x , int y)
{
	int nrnheight = this->network->n_branches_per_neuron* 16;

	pen.SetColour(0,0,0);
	pen.SetStyle(wxCAP_PROJECTING);

	//pen.SetColour(0,0,0);
	wxColor* nrncolor = &incolor;

	if (n->V > 18.)
		nrncolor = &firecolor;

	pen.SetColour( *nrncolor);
	pen.SetWidth(3);
	dc.SetBrush(wxBrush(*nrncolor, wxSOLID));
	dc.SetPen(pen);

	dc.DrawCircle(x+60, y+8, 5);

	wxPoint px1(x+60, y);
	wxPoint px2(x+60, y-nrnheight-10);
	dc.DrawLine(px1, px2);

/*
		if (n->V ==0.0)
			pen.SetColour(240, 100, 20);
		else
			pen.SetColour(20, 180, 20);
*/


	pen.SetWidth(2);
	dc.SetPen(pen);
	//dc.DrawRectangle(x-4, y-4, 8,8);


	dc.DrawText(wxString::Format(wxT("%d"), n->total_spikes), x, y+4);

	//float xx = x;
	float yy = y-10;
	int nbr =0;

	for (branch_iter bi = n->branches.begin(); bi != n->branches.end(); ++bi)
	{
		LABranch* b = *bi;
		yy -= 15;

		//float v = b->depol / 20.0;
		//if (v > 1.0) v = 1.0;
		wxPoint pt1(x+60, yy);

		int dx = x+60+56;
		wxPoint pt2(dx, yy);

		if (nbr%2)
		{
			pt2.x = x;
		}

		//dc.DrawText(wxString::Format(wxT("%d"), b->branch_spikes), x+70, yy-6);
		if (b->strength > 1.0)
			pen.SetWidth( b->strength*1.0);
		else
			pen.SetWidth(1);

		pen.SetColour( *nrncolor);
		dc.SetPen(pen);
		dc.DrawLine(pt2, pt1);

		int xxx = x+66;
		if (nbr%2) xxx = x+54;
		int lastinput = -1;

		for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
		{
			LASynapse*  s = *si;

			int space =1;
			if (lastinput == -1) lastinput = s->source_nrn->input_id;
			else if ( s->source_nrn->input_id != lastinput)
			{
				space = 3;
				lastinput = s->source_nrn->input_id;
			}

			if ( nbr%2)
				xxx -= space;
			else
				xxx+= space;


			wxPoint sp1(xxx, yy);

			int ddd = yy - ( s->weight )*7.0;
			wxPoint sp2(xxx, ddd);
			//if (s->source_nrn->V > -10)
			//	pen.SetColour(200.0, 180.0, 0);

			pen.SetColour(200, 200, 200);

			pen.SetWidth(1);
			dc.SetPen(pen);
			dc.DrawLine(sp1, sp2);
		}
		//nbr++;
	}
}


void LAPanel::OnPaint(wxPaintEvent& event)
{
	wxPaintDC dc(this);
	wxSize size = GetClientSize();
	int width = size.GetWidth();
	dc.Clear();
	wxPen pen ;
	wxBrush brush;


	if (!this->network || !this->network->neurons.size()) return;

	int nrnheight = this->network->n_branches_per_neuron* 20;
	int nrnwidth = 160;
	int x = 20;
	int y = nrnheight+40;

	pen.SetColour(0,0,0);
	pen.SetStyle(wxCAP_PROJECTING);

	nrn_iter pc  = network->pyr_list.begin();
	nrn_iter ic = network->in_list.begin();
	int nni =0;

	int totshown =0;
	while(totshown++ < 10)
	{
		if (x+100 >= width)
		{
			x = 20;
			y +=nrnheight+80;
		}

		LANeuron* n;

		if (nni++%5 !=0)
		{
			n = *pc;
			drawPyrCell(dc, n, x, y);
			pc++;
		}
		else
		{

			n = *ic;
			drawInCell(dc, n, x, y);
			ic++;
		}

		nrnpoints[n->nid] = wxPoint(x, y);
		x += nrnwidth;

		if (pc == network->pyr_list.end() || ic == network->in_list.end())
			break;
	}


}





bool LAApp::OnInit()
{
	LAWindow* window = new LAWindow(wxT("LA Model"), NULL);
	//window->network->window = window;
	window->Show(true);
	return true;
}




BEGIN_EVENT_TABLE(LAWindow, wxFrame)

	EVT_MENU(MENU_RUN,  LAWindow::OnRun)
	EVT_MENU_RANGE(MENU_CS, MENU_CS+20,  LAWindow::OnRun)
	//EVT_MENU_RANGE(MENU_CSUS, MENU_CSUS+20,  LAWindow::OnRun)
	EVT_MENU_RANGE(MENU_INTERSTIM, MENU_INTERSTIM+5,  LAWindow::OnRun)
	//EVT_MENU_RANGE(MENU_STRONGTRAIN, MENU_STRONGTRAIN,  LAWindow::OnRun)
	EVT_MENU_RANGE(MENU_PROTO, MENU_PROTO+4,  LAWindow::OnRun)
	EVT_MENU(MENU_QUIT, LAWindow::OnQuit)

END_EVENT_TABLE()




IMPLEMENT_APP(LAApp)


