//
// Version: $Id: wxmodel.h 162 2014-01-08 15:09:14Z gk $
//
// WxWidgets GUI for visualization and demonstration purposes
//
/* 
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





#include <wx/wx.h>
#include "wx/glcanvas.h"
#include <wx/thread.h>
#include <wx/event.h>

// extern
class LANetwork;

class LAApp: public wxApp
{
	public:
	LAWindow* window;
	LANetwork* network;
	virtual bool OnInit();

};


typedef map<int, wxPoint>::iterator pt_iterator;

class LAPanel: public wxPanel
{
	public:
	LAPanel(wxFrame*);
	LANetwork* network;
	bool hasStarted ;
	map<int, wxPoint> nrnpoints;


	protected:

	void drawPyrCell(wxPaintDC& dc, LANeuron* n, int x , int y);
	void drawInCell(wxPaintDC& dc, LANeuron* n, int x , int y);
	void OnPaint(wxPaintEvent& event);
	void OnTimer(wxCommandEvent& );
	void OnClick(wxMouseEvent& );
};


class LAGLPane;


class LAWindow: public wxFrame
{
	public:
	LAPanel* network_panel;
	LAGLPane* network_glpanel;
	LANetwork* network;
	wxTimer* timer;
	int timersRun;

	LAWindow(const wxString& title, LANetwork* net);
	void OnRun(wxCommandEvent& ev);
	void OnRefresh(wxCommandEvent& ev);
	void OnQuit(wxCommandEvent& ev);

	void UpdatePanel()
	{
		this->network_panel->Refresh();
		this->SetStatusText(wxString::Format(wxT("T=%d"), this->network->T));
		this->Update();
	}


	protected:

	void OnTimer(wxCommandEvent& );


	DECLARE_EVENT_TABLE()
};



class LAGLPane : public wxGLCanvas
{
    wxGLContext*	m_context;
    GLUquadric* m_quadric;
 
	public:
	LANetwork* m_network;

	LAGLPane(wxFrame* parent, LANetwork* net, int* args);
	virtual ~LAGLPane();
 
	void resized(wxSizeEvent& evt);
 
	int getWidth();
	int getHeight();
 
	void render(wxPaintEvent& evt);
	void prepare3DViewport(int topleft_x, int topleft_y, int bottomrigth_x, int bottomrigth_y);
	void prepare2DViewport(int topleft_x, int topleft_y, int bottomrigth_x, int bottomrigth_y);
 
	void drawPyrNeuron( LANeuron* n);
	// events
	//void mouseMoved(wxMouseEvent& event);
	//void mouseDown(wxMouseEvent& event);
	//void mouseWheelMoved(wxMouseEvent& event);
	//void mouseReleased(wxMouseEvent& event);
	//void rightClick(wxMouseEvent& event);
	//void mouseLeftWindow(wxMouseEvent& event);
	//void keyPressed(wxKeyEvent& event);
	//void keyReleased(wxKeyEvent& event);
	//
	DECLARE_EVENT_TABLE()
};



