/*
 * This is basicaly a colormap, just it works...
 * Only handles diffuse colors for now (rest is constant)
 * Peter-Pike Sloan
 */

#include <Classlib/NotFinished.h>
#include <Dataflow/Module.h>
#include <Datatypes/ColormapPort.h>
#include <Datatypes/Colormap.h>
#include <Datatypes/GeometryPort.h>

#include <Math/CatmullRomSpline.h>
#include <Malloc/Allocator.h>

#include <TCL/TCLvar.h>
#include <TCL/TCLTask.h>

#include <Geom/Color.h>

#define Colormap XColormap

#include <Geom/GeomOpenGL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

#include <tcl/tcl/tcl.h>
#include <tcl/tk/tk.h>

#undef Colormap

#include <stdio.h>

// tcl interpreter corresponding to this module

extern Tcl_Interp* the_interp;

// the OpenGL context structure

extern "C" GLXContext OpenGLGetContext(Tcl_Interp*, char*);


struct ColorPoint {
  float  _t;  // time - 
  Color _rgb;
  HSVColor _hsv;  // you have both reps...
};

class GenTransferFunc : public Module {

  ColormapOPort         *outport;  // outputs a colormap
  GeometryOPort         *ogeom;  
  TCLint		RGBorHSV; // which mode
  TCLint		lineVSspline; // linear vs. spline interpolate

  Array1< ColorPoint > points;   // actual line(s)

  Array1< float >	alphas;   // alpha polyline
  Array1< float >       aTimes;   // alpha times

  Array1< Point >       otherGraph; // xyz to RGB or HSV

  int			activeLine; // active RGB,HSV,alpha,-1
  int			selNode;    // selected node
  int                   graphStat;  // status of otherGraph

  int			bdown;      // wether mouse button is down
  int 			whichWin;   // which window is selected

  int			winX[3],winY[3];  // integer size
  float			winDX[2],winDY[2]; // scales (so it's square)

  GLXContext		ctxs[3];    // OpenGL Contexts
  Display		*dpy[3];
  Window		win0;
  Window		win1;
  Window		win2;

  ColormapHandle	cmap;  // created once, first execute...

public:
  GenTransferFunc( const clString& id);
  GenTransferFunc( const GenTransferFunc&, int deep);

  void DrawGraphs(void); // this function just blasts away...

  void DoMotion(int win, int x, int y);

  void DoDown(int win, int x, int y, int button);
  void DoRelease(int win, int x, int y, int button);

  float GetTime(int x, int y) 
    { 
      float v = x/(1.0*winX[whichWin]);
      if (v > 1.0) v = 1.0;	
      if (v < 0.0) v = 0.0;
      return v;
    }
  float GetVal(int x, int y)  
    { 
      float v= (winY[whichWin]-y)/(1.0*winY[whichWin]); 
      if (v > 1.0) v = 1.0;	
      if (v < 0.0) v = 0.0;
      return v;
    }

  void GetClosest(float time, float val, int& cline, int& cpoint);

  void Resize(int win);

  virtual ~GenTransferFunc();
  virtual Module* clone( int deep );
  virtual void execute();
  void tcl_command( TCLArgs&, void* );

  int makeCurrent(void);

  void BuildOther(void); // builds array in other space (every 4 pixels)

};

extern "C" {
  Module* make_GenTransferFunc( const clString& id)
    {
      return scinew GenTransferFunc(id);
    }

};

GenTransferFunc::GenTransferFunc( const clString& id)
:Module("GenTransferFunc",id,Source),
 RGBorHSV("rgbhsv",id,this),lineVSspline("linespline",id,this),
 activeLine(-1),selNode(-1),graphStat(-1),bdown(-1),whichWin(-1)
{

  cmap = 0; // start as nothing...

  // Create the output port
  outport = scinew ColormapOPort(this,"Colormap",ColormapIPort::Atomic);
  add_oport(outport);

  ogeom=scinew GeometryOPort(this, "Geometry", GeometryIPort::Atomic);
  add_oport(ogeom);
 
  // initialize the transfer function

  ColorPoint p;

  p._t = 0.0;
  p._rgb = Color(1,0,0);
  p._hsv = HSVColor(p._rgb);

  points.add(p);

  p._t = 1.0;
  p._rgb = Color(0,0,1);
  p._hsv = HSVColor(p._rgb);

  points.add(p);

  alphas.add(1);
  alphas.add(1);
  aTimes.add(0);
  aTimes.add(1);

  ctxs[0] = ctxs[1] = ctxs[2] = 0;
}

GenTransferFunc::GenTransferFunc(const GenTransferFunc& copy, int deep)
: Module(copy,deep),
  RGBorHSV("rgbhsv",id,this),lineVSspline("linespline",id,this),
  activeLine(-1),selNode(-1),graphStat(-1)
{
  NOT_FINISHED("GenTransferFunc copy constructor type thing");
}

GenTransferFunc::~GenTransferFunc()
{

}

Module* GenTransferFunc::clone(int deep)
{
  return scinew GenTransferFunc(*this,deep);
}

void GenTransferFunc::BuildOther(void)
{
//  if (RGBorHSV) {
    
//  }
}

// draws/updates the graphs...

void drawBox(float x, float y, float dx, float dy)
{
  glBegin(GL_LINE_LOOP);
  glVertex2f(x-dx,y-dy);
  glVertex2f(x-dx,y+dy);
  glVertex2f(x+dx,y+dy);
  glVertex2f(x+dx,y-dy);
  glEnd();
}

void GenTransferFunc::Resize(int win)
{

  // do a make current...

  if ( ! makeCurrent() )
    return;

  switch (win) {
  case 0:
    glXMakeCurrent(dpy[win],win0,ctxs[win]);
    break;
  case 1:
    glXMakeCurrent(dpy[win],win1,ctxs[win]);
    break;
  case 2:
    glXMakeCurrent(dpy[win],win2,ctxs[win]);
    break;
  }

  if (win < 2) {  // DX DY only for first 2...
    int sx,sy;
    sx = winX[win];
    sy = winY[win];

    winDX[win] = 1.0/100.0;  // X is normaly bigger...

    winDY[win] = sx/sy*winDX[win];
  }

  // make current...

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(-1,-1,0);
  glScalef(2,2,2);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

//  gluOrtho2D(0,winX[win],0,winY[win]);
  

  glEnable(GL_COLOR_MATERIAL);  // more state?
  glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);
  glDisable(GL_TEXTURE_2D);

  glViewport(0,0,winX[win],winY[win]);
  glClear(GL_COLOR_BUFFER_BIT);

  TCLTask::unlock();

}

void GenTransferFunc::DrawGraphs()
{

  static double  colors[] = { 1.,0.,0.,   0.,1.,0.,   0.,0.,1. };

  if (!makeCurrent())
    return; // don't care if you can't make current...

  glXMakeCurrent(dpy[0],win0,ctxs[0]);
  // other state should be fine...

  glDrawBuffer(GL_BACK);
  glClear(GL_COLOR_BUFFER_BIT);
  for(int i=0;i<3;i++) {

    float mul=0.5;
    if (activeLine == i*2)
      mul = 1.0;

    glColor3f(colors[i*3 + 0]*mul,colors[i*3 + 1]*mul,colors[i*3 + 2]*mul);
    glBegin(GL_LINE_STRIP);
    for(int j=0;j<points.size();j++) {
      glVertex2f(points[j]._t,points[j]._rgb[i]);
    }
    glEnd();
  }

  // now draw the alpha curve...

  float mul = 0.5;
  if (activeLine == 7)
    mul = 0.7;

  glColor3f(mul,mul,mul);
  glBegin(GL_LINE_STRIP);
  
  for(i=0;i<alphas.size();i++) {
    glVertex2f(aTimes[i],alphas[i]);
  }
  glEnd();

  for(i=0;i<alphas.size();i++) {
    mul = 0.5;
    if (i == selNode)
      mul = 0.7;
    glColor3f(mul,mul,mul);
    drawBox(aTimes[i],alphas[i],winDX[0],winDY[0]);
  }

  // now the polylines have been draw, draw the points...

  for(i=0;i<3;i++) {
    for(int j=0;j<points.size();j++) {
      float mul = 0.6;
      if (selNode == j*2)
	mul = 1.0;

      glColor3f(colors[i*3 + 0]*mul,colors[i*3 + 1]*mul,colors[i*3 + 2]*mul);
      drawBox(points[j]._t,points[j]._rgb[i],winDX[0],winDY[0]);
    }
  }

  glXSwapBuffers(dpy[0],win0);

  // go to the last window...

  glXMakeCurrent(dpy[2],win2,ctxs[2]);
  glDrawBuffer(GL_BACK);
  glClear(GL_COLOR_BUFFER_BIT);

  // 1st lay down the "pattern"

  glDisable(GL_BLEND);

  glBegin(GL_QUAD_STRIP);
  for(i=0;i<10;i++) {
    glColor3f(i&1,i&1,i&1);
    glVertex2f(0,i/9.0);
    glVertex2f(1,i/9.0);
  }
  glEnd();
  glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_TRUE);
  
  glBegin(GL_QUAD_STRIP);
  for(i=0;i<alphas.size();i++) {
    glColor4f(0,0,0,alphas[i]);
    glVertex2f(aTimes[i],0);
    glVertex2f(aTimes[i],1);
  }
  glEnd();

  glBlendFunc(GL_DST_ALPHA,GL_ONE_MINUS_DST_ALPHA);
  glEnable(GL_BLEND);
  glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);

  glBegin(GL_QUAD_STRIP);
  for(int j=0;j<points.size();j++) {
    glColor3f(points[j]._rgb[0],points[j]._rgb[1],points[j]._rgb[2]);
    glVertex2f(points[j]._t,0);
    glVertex2f(points[j]._t,1);
  }
  glEnd();
  glXSwapBuffers(dpy[2],win2);

  
  // here you update this colormap thing...
  
  if (cmap.get_rep()) {
    Array1<Color> ncolors(points.size());
    Array1<float> times(points.size());

    for(int i=0;i<points.size();i++) {
      ncolors[i] = points[i]._rgb;
      times[i] = points[i]._t;
    }
    
    cmap.get_rep()->SetRaw(ncolors,times,alphas,aTimes);

    glXMakeCurrent(dpy[0],None,NULL);

    TCLTask::unlock();

    ogeom->flushViews();

    return;

  }

  glXMakeCurrent(dpy[0],None,NULL);

  TCLTask::unlock();

}

// this is how the gui talks with this thing

void GenTransferFunc::tcl_command( TCLArgs& args, void* userdata)
{
  if (args.count() < 2) {
    args.error("No command for GenTransferFunc");
    return;
  }

  // mouse commands are motion, down, release

  if (args[1] == "mouse") {
    int x,y,whichw,whichb;
    
    args[4].get_int(x);
    args[5].get_int(y);

    args[2].get_int(whichw); // which window it was in

    if (args[3] == "motion") {
      if (bdown == -1) // not buttons down!
	return; 
      DoMotion(whichw,x,y);
    } else {
      args[6].get_int(whichb); // which button it was
      if (args[3] == "down") {
	DoDown(whichw,x,y,whichb);
      } else { // must be release
	DoRelease(whichw,x,y,whichb);
      }
    }
  } else  if (args[1] == "resize") { // resize event!
    int whichwin;
    args[4].get_int(whichwin);
    args[2].get_int(winX[whichwin]);
    args[3].get_int(winY[whichwin]);
    Resize(whichwin);  // kick of the resize function...
  } else if (args[1] == "expose") {
    int whichwin;
    args[2].get_int(whichwin);
    Resize(whichwin); // just sets up OGL stuff...
    
  }else {
    Module::tcl_command(args, userdata);
  }
}

// you just have to find the time interval first
// then look at the values...

void GenTransferFunc::
GetClosest(float time, float val, int& cline, int& cpoint)
{
  int i;
  
  float minP=10000,minA=10000;

  for(i=0;i<points.size() && (points[i]._t <= time); i++) {
    ; // do nothing - just getting index...
  }

  // now i is the value 1 greater...

  if (i == points.size()) {
    i = points.size() -1;
  } 

  // now find the closest point

  float d1,d2;

  d1 = time - points[i-1]._t;  d1 = d1*d1;
  d2 = points[i]._t - time;  d2 = d2*d2;

  if (!whichWin) { // RGB space is 0...
    
    float dy1,dy2;
    for(int j=0;j<3;j++) {
      dy1 = points[i-1]._rgb[j] - val;
      dy2 = points[i]._rgb[j]-val;

      float min1 = d1 + dy1*dy1;
      float min2 = d2 + dy2*dy2;

      if (min1 < minP) {
	minP = min1;
	cpoint = i-1;
	cline = j*2;
      }
      if (min2 < minP) {
	minP = min2;
	cpoint = i;
	cline = j*2;
      }
    }
  } else { // HSV space...
    printf("Fuck!\n");
  }

  // now check alpha...

  for(i=0;i<alphas.size() && (aTimes[i] <= time); i++) {
    ; // do nothing - just getting index...
  }

  // now i is the value 1 greater...

  if (i == alphas.size()) {
    i = alphas.size() -1;
  } 


  d1 = time - aTimes[i-1];  d1 = d1*d1;
  d2 = aTimes[i] - time;  d2 = d2*d2;

  d1 += (alphas[i-1]-val)*(alphas[i-1]-val);
  d2 += (alphas[i]-val)*(alphas[i]-val);
  
  if (d1 < d2) {
    if (d1 < minP) {
      cline = 7;
      cpoint = i-1;
    }
  } else {
    if (d2 < minP) {
      cline = 7;
      cpoint = i;
    }
  }

  // now it should be ok...

}	

// button must be down

void GenTransferFunc::DoMotion(int win, int x, int y)
{
  if ((selNode == -1) || (activeLine == -1))  // this shouldn't happen!
    return;

  // this means you have to move the selected node...

  float time = GetTime(x,y); // get the valid time value...

  float val = GetVal(x,y);   // this remaps the Y coord!

  // end conditions are special cases - can't change x!

  if (activeLine == 7) {
    if (selNode && selNode < alphas.size()-1) {
      if (aTimes[selNode-1] > time)
	time = aTimes[selNode-1];
      if (aTimes[selNode+1] < time)
	time = aTimes[selNode+1];
    } else {
      time = aTimes[selNode];
    }
    aTimes[selNode] = time;
  } else {
    if (selNode && selNode < points.size()-1) {
      if (points[selNode-1]._t > time)
	time = points[selNode-1]._t;
      if (points[selNode+1]._t < time)
	time = points[selNode+1]._t;
    } else {
      time = points[selNode]._t;
    }
    points[selNode]._t = time;
  }

  // now update this guy

  if (activeLine < 6) {
    if (activeLine&1) {
      if (activeLine == 1) { // this is hue - 0-360...
	val *= 360;
      }
      points[selNode]._hsv[activeLine>>1] = val;
      points[selNode]._rgb = Color(points[selNode]._hsv);

    } else {
      points[selNode]._rgb[activeLine>>1] = val;
      points[selNode]._hsv = HSVColor(points[selNode]._rgb);
    }
  } else { // it is the alpha value...
    alphas[selNode] = val;
  }

  // now you just have to redraw this

  DrawGraphs();

}

// Button 1 is modify, 2 is insert, 3 is delete

void GenTransferFunc::DoDown(int win, int x, int y, int button)
{
  // you have to find the point to select...

  float time,val;
  int cline,cpoint; // closest info...

  time = GetTime(x,y);
  val  = GetVal(x,y);

  whichWin = win;  // this has to be pre closest selection...
  GetClosest(time,val,cline,cpoint);
  
  ColorPoint p;

  if (button == 2) {
    // you have to find what side of the closest point
    // you are on...

    activeLine = cline;
    
    p._t = time;

    p._rgb[cline>>1] = val;
    p._rgb[((cline>>1)+1)%3] = 0;
    p._rgb[((cline>>1)+2)%3] = 1;
    
    p._hsv = HSVColor(p._rgb);

    if (cline != 7) {
      if (cpoint && cpoint != points.size()-1) {
	if (points[cpoint]._t <= time)
	  cpoint = cpoint+1;
	
      } else {
	if (!cpoint) {
	  cpoint = 1;
	}
      }

      points.insert(cpoint,p);

    } else { // this is the alpha curve...
      if (cpoint && cpoint != alphas.size()-1) {
	if (aTimes[cpoint] <= time)
	  cpoint = cpoint+1;
      }	
      else {
	if (!cpoint) {
	  cpoint = 1;
	}
      }
  
      alphas.insert(cpoint,val);
      aTimes.insert(cpoint,time);
    }

  }

  selNode = cpoint;
  activeLine = cline;
  
  bdown = button;

  DrawGraphs(); // just update the display (diferent box for selected)

}

void GenTransferFunc::DoRelease(int win, int x, int y, int button)
{
  // just fake a move for the final posistion...

  DoMotion(win,x,y);

  printf("%d Button!\n",button);

  switch(button) {
  case 1:
  case 2:
    // do nothing - just move/add...
    break;
  case 3:
    if (activeLine == 7) {
      if (selNode && (selNode != alphas.size()-1)) {
	alphas.remove(selNode);
	aTimes.remove(selNode);
      }

    }  else {
      if (selNode && (selNode != points.size()-1))
	points.remove(selNode);
    }
    break;
      
  }

  selNode = -1; // leave activeLine...
  bdown = -1;

  DrawGraphs();

}

void GenTransferFunc::execute(void)
{
  if (!cmap.get_rep()) {
    Array1<Color> ncolors(points.size());
    Array1<float> times(points.size());

    for(int i=0;i<points.size();i++) {
      ncolors[i] = points[i]._rgb;
      times[i] = points[i]._t;
    }
    cmap = scinew ColorMap(ncolors,times,alphas,aTimes);
  }
  outport->send(cmap);
}

int GenTransferFunc::makeCurrent(void)
{
  Tk_Window tkwin;
  Window win;

  // lock a mutex
  TCLTask::lock();

  if (!ctxs[0]) {
    clString myname(clString(".ui")+id+".f.gl1.gl");
    tkwin = Tk_NameToWindow(the_interp, myname(),Tk_MainWindow(the_interp));

    if (!tkwin) {
      cerr << "Unable to locate window!\n";

      // unlock mutex
      TCLTask::unlock();
      return 0;
    }
    winX[0] = Tk_Width(tkwin);
    winY[0] = Tk_Height(tkwin);

    dpy[0] = Tk_Display(tkwin);
    win0 = Tk_WindowId(tkwin);

    ctxs[0] = OpenGLGetContext(the_interp,myname());

    // check if it was created
    if(!ctxs[0])
      {
	cerr << "Unable to create OpenGL Context!\n";
	TCLTask::unlock();
	return 0;
      }
  }	

  if (!ctxs[1]) {
    clString myname(clString(".ui")+id+".f.gl2.gl");
    tkwin = Tk_NameToWindow(the_interp, myname(),Tk_MainWindow(the_interp));

    if (!tkwin) {
      cerr << "Unable to locate window!\n";

      // unlock mutex
      TCLTask::unlock();
      return 0;
    }
    
    winX[1] = Tk_Width(tkwin);
    winY[1] = Tk_Height(tkwin);

    dpy[1] = Tk_Display(tkwin);
    win1 = Tk_WindowId(tkwin);

    ctxs[1] = OpenGLGetContext(the_interp,myname());

    // check if it was created
    if(!ctxs[1])
      {
	cerr << "Unable to create OpenGL Context!\n";
	TCLTask::unlock();
	return 0;
      }
  }	

  if (!ctxs[2]) {
    clString myname(clString(".ui")+id+".f.gl3.gl");
    tkwin = Tk_NameToWindow(the_interp, myname(),Tk_MainWindow(the_interp));

    if (!tkwin) {
      cerr << "Unable to locate window!\n";

      // unlock mutex
      TCLTask::unlock();
      return 0;
    }
    winX[2] = Tk_Width(tkwin);
    winY[2] = Tk_Height(tkwin);

    dpy[2] = Tk_Display(tkwin);
    win2 = Tk_WindowId(tkwin);

    ctxs[2] = OpenGLGetContext(the_interp,myname());

    // check if it was created
    if(!ctxs[2])
      {
	cerr << "Unable to create OpenGL Context!\n";
	TCLTask::unlock();
	return 0;
      }
  }	

  // ok, 3 contexts are created - now you only need to
  // do the other crap...

  return 1;

}
