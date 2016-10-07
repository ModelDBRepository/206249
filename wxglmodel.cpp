/* 
 
    Version: $Id: wxmodel.cpp 134 2013-11-22 14:30:58Z gk $
    Implementation of GUI interface

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
#include "wxglmodel.h"
#include "wx/sizer.h"
#include "wx/glcanvas.h"
#include <GL/glu.h>
#include <GL/gl.h>
#include <GL/glut.h>


// The next global variable controls the animation's state and speed.
float RotateAngle = 0.0f;		// Angle in degrees of rotation around y-axis
float Azimuth = 20.0;			// Rotated up or down by this amount

float AngleStepSize = 3.0f;		// Step three degrees at a time
const float AngleStepMax = 10.0f;
const float AngleStepMin = 0.1f;

int WireFrameOn = 1;			// == 1 for wire frame mode

GLfloat light_position2[] = { 0, 1,1,0};

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




/*
wxDECLARE_EVENT(wxEVT_COMMAND_MYTHREAD_COMPLETED, wxThreadEvent);
wxDECLARE_EVENT(wxEVT_COMMAND_MYTHREAD_UPDATE, wxThreadEvent);
*/


GLuint raw_texture_load(const char *filename, int width, int height)
{
     GLuint texture =0;
     unsigned char *data;
     FILE *file;
 
     // open texture data
     file = fopen(filename, "rb");
     if (file == NULL) return 0;
 
     // allocate buffer
     int sq;
     data = (unsigned char*) malloc(sq = width * height *4  );
     cout << "sq=" << sq<<endl;
 
     // read texture data
     sq = fread(data,1,  width * height *3 ,  file);

     cout << "read "<<sq << endl;
     fclose(file);
 
     // allocate a texture name
     glEnable(GL_TEXTURE_2D);
     glGenTextures(1, &texture);
 
     // select our current texture
     glBindTexture(GL_TEXTURE_2D, texture);
 
     // select modulate to mix texture with color for shading
     glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
 
 
     // when texture area is small, bilinear filter the closest mipmap
     //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
     // when texture area is large, bilinear filter the first mipmap
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
 
     // texture should tile
     glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
     glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
 
     // build our texture mipmaps
       //float pixels[] = { 0.0f, 0.0f, 0.0f,   1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,   0.0f, 0.0f, 0.0f };
    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 2, 2, 0, GL_RGB, GL_FLOAT, pixels);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);

    // gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, width, height, GL_RGB, GL_UNSIGNED_BYTE, data);
    glDisable(GL_TEXTURE_2D);
 
     free(data);
 
     return texture;
 }


class MyThread : public wxThread
{
	 public:
	 LANetwork* m_network;
	 MyThread(LAWindow *handler) : wxThread(wxTHREAD_JOINABLE)
	 { 
		m_pHandler = handler ;
		m_network = handler->network;
	 }
	 ~MyThread() {}

	 protected:
	 virtual ExitCode Entry();
	 LAWindow *m_pHandler;
};



wxThread::ExitCode MyThread::Entry()
{
	this->m_network->RunStore2(3, 1, 180, 0);
	exit(0);
	return (wxThread::ExitCode)0; // success
}


LAWindow::LAWindow(const wxString& title, LANetwork* netw) 
	:wxFrame(NULL, wxID_ANY, title, wxDefaultPosition,wxSize(1024, 768))
{
	Centre(); 
	timersRun =0;

	int args[] = {WX_GL_RGBA, WX_GL_DOUBLEBUFFER, WX_GL_DEPTH_SIZE, 16, 0};
	this->network = netw;
	this->network_glpanel = new LAGLPane(this, this->network, args);

	//this->network_panel->network = this->network;
	//this->network->wx_window = this;
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

	//SetMenuBar(mb);
	//
	Connect(wxEVT_TIMER, wxCommandEventHandler(LAWindow::OnTimer));
}



void LAWindow::OnQuit(wxCommandEvent& ev)
{
	printf("cleaning up... \n");
	this->network->Cleanup();
	Close(true);
	printf("Done\n");
}


void LAWindow::OnTimer(wxCommandEvent& ev)
{
	//printf("Window Timer!\n");
	timersRun++;

	float angle = this->network_glpanel->angle;
	angle += 0.2;
	if (angle>360.0) angle =0.;
	this->network_glpanel->angle =angle;
	this->network_glpanel->Refresh(false);

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
				this->network->CreateFearNet(30, 100, 10, 4);
				for (int i =0;   i < network->n_inputs; i++)
				{
					SetStatusText(wxString::Format(wxT("Training # %d..."), i));
					//network->Stimulate(i, 2);
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
					//network->Stimulate(i, 1);
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


	SetStatusText(wxT("Done!"));
	printf("Done running\n");
}






inline static void colormap1(wxBrush& pen, float value )
{
	pen.SetColour(250, 250-value*200, value*200); 
}


inline static void colormap1(wxPen& pen, float value )
{
	pen.SetColour(value*240.0, 0, 250 - value*200.0);
}



static int palette[][3] = {
	{200, 50,50},
	{50, 200,50},
	{250 , 250,100},
	{200, 190,3},
	{255,0,255 },
	{4, 3,200},
	{255,0,0 },
	{128,0,0 },
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


bool LAApp::OnInit()
{

	this->network = new LANetwork();
	LANetwork::SetRandomSeed(1980);
	this->network->CreateFearNet(20, 20, 10, 6);

	LAWindow* window = new LAWindow(wxT("LA Model"), this->network);
	window->network = network;

	window->timer = new wxTimer(window, 13);
	window->timer->Start(30);

	

	MyThread* m_pThread = new MyThread(window);
	cout << "Starting thread!!\n"<<endl;

	int err;
	err = m_pThread->Create();
	if (err != wxTHREAD_NO_ERROR)
	{
		perror("TA");
		cout << ("Can't create the thread!") << endl;
		printf("Error %d\n",  err);
		delete m_pThread;
		m_pThread = NULL;
	}
	else 
	{

		err =  m_pThread->Run() ;
		if ( err != wxTHREAD_NO_ERROR )
		{
			perror("TA");
			cout << ("Can't create the thread!") << endl;
			printf("Error %d\n",  err);
			delete m_pThread;
			m_pThread = NULL;
		}
	}


	window->Show(true);


	return true;
}



IMPLEMENT_APP(LAApp)



BEGIN_EVENT_TABLE(LAWindow, wxFrame)

	EVT_MENU(MENU_RUN,  LAWindow::OnRun)
	EVT_MENU_RANGE(MENU_CS, MENU_CS+20,  LAWindow::OnRun)
	//EVT_MENU_RANGE(MENU_CSUS, MENU_CSUS+20,  LAWindow::OnRun)
	EVT_MENU_RANGE(MENU_INTERSTIM, MENU_INTERSTIM+5,  LAWindow::OnRun)
	//EVT_MENU_RANGE(MENU_STRONGTRAIN, MENU_STRONGTRAIN,  LAWindow::OnRun)
	EVT_MENU_RANGE(MENU_PROTO, MENU_PROTO+4,  LAWindow::OnRun)
	EVT_MENU(MENU_QUIT, LAWindow::OnQuit)

END_EVENT_TABLE()




BEGIN_EVENT_TABLE(LAGLPane, wxGLCanvas)
	//EVT_MOTION(BasicGLPane::mouseMoved)
	//EVT_LEFT_DOWN(BasicGLPane::mouseDown)
	//EVT_LEFT_UP(BasicGLPane::mouseReleased)
	//EVT_RIGHT_DOWN(BasicGLPane::rightClick)
	//EVT_LEAVE_WINDOW(BasicGLPane::mouseLeftWindow)
	//EVT_SIZE(BasicGLPane::resized)
	//EVT_KEY_DOWN(BasicGLPane::keyPressed)
	//EVT_KEY_UP(BasicGLPane::keyReleased)
	//EVT_MOUSEWHEEL(BasicGLPane::mouseWheelMoved)
	EVT_PAINT(LAGLPane::render)
END_EVENT_TABLE()







// Vertices and faces of a simple cube to demonstrate 3D render
// source: http://www.opengl.org/resources/code/samples/glut_examples/examples/cube.c
GLfloat v[8][3];
GLint faces[6][4] = {  /* Vertex indices for the 6 faces of a cube. */
    {0, 1, 2, 3}, {3, 2, 6, 7}, {7, 6, 5, 4},
    {4, 5, 1, 0}, {5, 6, 2, 1}, {7, 4, 0, 3} };
 
 
GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 2.0};  /* Red diffuse light. */

LAGLPane::LAGLPane(wxFrame* parent, LANetwork* net, int* args) :
    wxGLCanvas(parent, wxID_ANY, args, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE)
{
	this->m_context = new wxGLContext(this);
	SetBackgroundStyle(wxBG_STYLE_CUSTOM);
	m_sphereTexture = 0;

	//glClearColor (0.0, 0.0, 0.0, 0.0);
	this->m_quadric =  gluNewQuadric();
	m_network = net;
	angle =0;

	// To avoid flashing on MSW

	this->m_sphereq = gluNewQuadric();

}


 
LAGLPane::~LAGLPane()
{
	delete m_context;
}
 
void LAGLPane::resized(wxSizeEvent& evt)
{
//	wxGLCanvas::OnSize(evt);
 
    Refresh();
}
 


GLfloat light_position[] = { 2, 2,2,0};
GLfloat light_amb[] = { 0.4,0.4, 0.4 } ;

GLfloat myglfog[] = {0.8, 0.8,1.0 , 1.};


/** Inits the OpenGL viewport for drawing in 3D. */
void LAGLPane::prepare3DViewport(int topleft_x, int topleft_y, int bottomrigth_x, int bottomrigth_y)
{

    //glClearColor(0.0, 0.0f, 0.0f, 1.0f); // Black Background
    
}


 
/** Inits the OpenGL viewport for drawing in 2D. */
void LAGLPane::prepare2DViewport(int topleft_x, int topleft_y, int bottomrigth_x, int bottomrigth_y)
{	
}


 
int LAGLPane::getWidth()
{
    return GetSize().x;
}
 


int LAGLPane::getHeight()
{
    return GetSize().y;
}
 


void LAdrawNeuron()
{
}


GLUquadricObj* myReusableQuadric = 0;

void LAGLPane::drawPyrNeuron( LANeuron* n)
{
	glPushMatrix();
	glRotatef(-90, 1,0,0);

	float b =n->V/90;
	if (b < -0.1) b = -0.1;

	if (n->isSpiking)
		glColor3f(1.0, 1.0, 1.0);
	else
		glColor3f(0.5,  0.5+b,  0.5-b);



	gluCylinder(this->m_quadric, 0.7, 0.02, 0.9, 16, 16);
	gluCylinder(this->m_quadric, 0.05, 0.01, 8.2, 8, 8);

	glRotatef(90, 1,0,0);
	glTranslatef(0, 1.0, 0);

	float ang = 0.;
	pthread_mutex_lock(&n->network->synapses_mutex);

	for (branch_iter bi = n->branches.begin(); bi != n->branches.end(); ++bi)
	{
		LABranch* b = *bi;
		glTranslatef(0, 0.3, 0);

		glPushMatrix();

		glRotatef(ang, 0,1,0);
		ang+= 55;

		glColor3f(0.6, 0.6, 0.6);

		glBegin(GL_LINES);
			glVertex3f(0, 0, 0);
			glVertex3f(0, 0, 2);
		glEnd();

		for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
		{
			LASynapse* s = *si;
			if (s->source_nrn->input_id >=0 && s->source_nrn->input_id < 3 )
			{
			
				int iid = s->source_nrn->input_id;

				if (s->source_nrn->isSpiking)
					glColor3f(255, 255, 255);
				else if (s->source_nrn->input_id <0)
				{
					glColor3f(50, 50, 50);
				}
				else
				{
					int* p = palette[iid];
					if (s->source_nrn->nid%2)
						glColor3f(p[0]/500., p[1]/500., p[2]/500.);
					else
						glColor3f(p[0]/300., p[1]/300., p[2]/300.);
				}

				glTranslatef(0.0, 0.0, 0.2);
				//glRotatef(-90, 1,0,0);

				float height = 0.05 + (0.12*(s->weight/3.0));
				gluDisk(this->m_quadric, 0, height, 9, 9);
				//gluCylinder(this->m_quadric, 0.03, 0.03, height, 8, 8);
				//glRotatef(90, 1,0,0);
			}
		}

		glTranslatef(0, -0.08, 0);
		glPopMatrix();
	}

	pthread_mutex_unlock(&n->network->synapses_mutex);
	glPopMatrix();
	glFlush();
}




void LAGLPane::render( wxPaintEvent& evt)
{
	if(!IsShown()) return;
	cout<< "Rendering" <<endl;

	if (m_sphereTexture ==0)
	{
		glEnable(GL_TEXTURE_2D);
		m_sphereTexture = raw_texture_load("./texture.raw", 256, 256 );
		printf("Loaded texture id %d \n", m_sphereTexture);
	}

	wxGLCanvas::SetCurrent(*m_context);
 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	//prepare3DViewport(0 ,0,getWidth(), getHeight());
	int topleft_x= 0, topleft_y =0;
	int bottomrigth_x= getWidth(), bottomrigth_y = getHeight();

	glClearColor(0.0, 0.0, 0.0, 1.0f); // Black Background
	glClearDepth(1.0f);	// Depth Buffer Setup
	glEnable(GL_DEPTH_TEST); // Enables Depth Testing
	glDepthFunc(GL_LEQUAL); // The Type Of Depth Testing To Do
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);


	glViewport(topleft_x, topleft_y, bottomrigth_x-topleft_x, bottomrigth_y-topleft_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	float ratio_w_h = (float)(bottomrigth_x-topleft_x)/(float)(bottomrigth_y-topleft_y);
	gluPerspective(35 /*view angle*/, ratio_w_h, 0.05 /*clip close*/, 200 /*clip far*/);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glLightfv(GL_LIGHT0, GL_SPECULAR, light_diffuse);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);


	glEnable( GL_POLYGON_SMOOTH );
	glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	glLineWidth(2.0);
	//glLightfv(GL_LIGHT0, GL_SPECULAR, light_amb);




	glTranslatef(1, -4, -18);
	glRotatef( this->angle, 0.0, 1, .05);


	glPushMatrix();
	glRotatef( 90, 1, 0, 0);
	gluQuadricDrawStyle(m_sphereq, GLU_LINE);
	glColor3f(0.3, 0.4, 0.5);
	gluSphere(m_sphereq, 20.2, 32, 16);
	glPopMatrix();

/*
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, m_sphereTexture);
	gluSphere(m_sphereq, 1.2, 32, 32);
	glDisable(GL_TEXTURE_2D);
	*/




	if (this->m_network && this->m_network->pyr_list.size())
	{
		int tot =0;
		for (nrn_iter ni = this->m_network->pyr_list.begin(); ni != this->m_network->pyr_list.end(); ++ni)
		{
		    LANeuron* n= *ni;

		    glPushMatrix();
		    glTranslatef((tot%4)*6.0 - 9., 0., (tot/4)*6.0 - 9.);
		    drawPyrNeuron(n);
		    glPopMatrix();

		    if (tot++>10) break;
		}
	}
 
    //glFlush();
    SwapBuffers();


}




