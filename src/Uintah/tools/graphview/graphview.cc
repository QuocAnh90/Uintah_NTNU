#include <errno.h>
#include <fstream>
#include <iostream>
#include <string>
#include "DaVinci.h"
#include "graphview.h"
#include "TaskGraph.h"

using namespace std;

// global variable instatiations
bool gQuit = false;
queue<Event> gEventQueue;

static TaskGraph* gGraph = 0;
static DaVinci* gDavinci = 0;

static void handle_event(const Event& event);
static void handle_console_input();

int
main(int argc, char* argv[])
{
  if (argc < 2) {
    cerr << "usage: " << argv[0] << " <uda directory>" << endl;
    return 1;
  }
  string infile = string(argv[1]) + "/taskgraph.xml";
  
  gGraph = TaskGraph::inflate(infile);
  if (!gGraph) {
    cerr << "Failed reading task graph, quitting" << endl;
    return 1;
  }
  
  gDavinci = DaVinci::run();
  gDavinci->setOrientation(DaVinci::BOTTOM_UP);
  gDavinci->setGraph(gGraph);
  
  while (!gQuit) {
    while (!gEventQueue.empty()) {
      Event event = gEventQueue.front();
      gEventQueue.pop();
      
      cout << "Handling event (type=" << event.type() << ")\n";
      cout << endl;
      handle_event(event);
    }
    
    /* Now that we've handled all pending events, we wait for input either
     * from the user at the console or from daVinci (presumably also
     * initiated by the user)
     */
    // select() changes the FD sets, so we have to reset them every time
    fd_set inputs, errors;
    int dv_fd = gDavinci->getOutput();
    FD_ZERO(&inputs);
    FD_ZERO(&errors);
    FD_SET(0 /* stdin */, &inputs);
    FD_SET(0 /* stdin */, &errors);
    FD_SET(dv_fd, &inputs);
    FD_SET(dv_fd, &errors);
	
    cout << endl << "? " << flush;
    
    while ((select(dv_fd + 1, &inputs, 0, &errors, 0) == -1) &&
	   (errno == EINTR))
      ;
    if ((FD_ISSET(0, &errors)) || FD_ISSET(dv_fd, &errors))
      gQuit = true;
    else if (FD_ISSET(0, &inputs)) {
      handle_console_input();
    } else if (FD_ISSET(dv_fd, &inputs))
      gDavinci->handleInput();
  }
  
  delete gDavinci;
  delete gGraph;
  
  return 0;
}

static
void
handle_event(const Event& event)
{ 
  switch (event.type()) {
  case DaVinci::EVT_DV_QUIT:
    gQuit = true;
    break;

  case DaVinci::EVT_DV_SELECT_NODE:
    {
      list<string> selected_nodes = gDavinci->getSelectedNodes();
      if (selected_nodes.size() == 1) {
	cout << selected_nodes.front() << endl;
	Task* pTask = gGraph->findTask(selected_nodes.front());
	if (pTask != NULL) {
	  cout << "\tCost (duration): " << pTask->getDuration() << endl;
	  cout << "\tMax Path Cost: " << pTask->getMaxInclusivePathCost()
	       << endl;
	  cout << "\tMax Path Percent: " << pTask->getMaxPathPercent()
	       << endl;
	}
	else {
	  cout << "\tError, task not found\n"; 
	}
      }
      else if (selected_nodes.size() > 1) {
	double total_cost = 0;
	for (list<string>::iterator iter = selected_nodes.begin();
	     iter != selected_nodes.end(); iter++) {
	  cout << *iter;
	  Task* pTask = gGraph->findTask(*iter);
	  if (pTask == NULL) {
	    cout << "\n\tError, task not found\n";
	    return;
	  }
	  cout << "\t(" << pTask->getDuration() << ")\n";
	  total_cost += pTask->getDuration();
	}
	cout << "\nTotal cost (duration): " << total_cost << endl;

      }
    }
    break;

  case DaVinci::EVT_DV_SELECT_EDGE:
    cout << gDavinci->getSelectedEdge() << endl;
    Edge* pEdge = gGraph->findEdge(gDavinci->getSelectedEdge());
    if (pEdge != NULL) {
      cout << "\tMax Path: " << pEdge->getMaxInclusivePathCost() << endl;
      cout << "\tMax Path Percent: " << pEdge->getMaxPathPercent()
	   << endl;
    }
    else {
      cout << "\tError, edge not found\n"; 
    }
    break;
  }
}

static void handle_console_input()
{
  string cmd;
  cin >> cmd;

  if ((cmd == "prune") || (cmd == "p") || (cmd == "P")) {
    float percent;
    cin >> percent;
    if (percent < 0) percent = 0;
    if (percent > 1) percent = 1;
    if (gGraph != NULL && gDavinci != NULL) {
      cout << "\nSetting threshold... " << percent << endl << endl;
      gGraph->setThresholdPercent(percent);
      gDavinci->setGraph(gGraph); // refresh graph
    }
  }

}
