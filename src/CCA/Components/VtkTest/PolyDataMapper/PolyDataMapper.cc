/*
  The contents of this file are subject to the University of Utah Public
  License (the "License"); you may not use this file except in compliance
  with the License.
  
  Software distributed under the License is distributed on an "AS IS"
  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
  License for the specific language governing rights and limitations under
  the License.
  
  The Original Source Code is SCIRun, released March 12, 2001.
  
  The Original Source Code was developed by the University of Utah.
  Portions created by UNIVERSITY are Copyright (C) 2001, 1994 
  University of Utah. All Rights Reserved.
*/

/*
 *  PolyDataMapper.cc:
 *
 *  Written by:
 *   Keming Zhang
 *   Department of Computer Science
 *   University of Utah
 *   January 2004
 *
 */

#include <iostream>
#include <SCIRun/Vtk/Port.h>
#include <CCA/Components/VtkTest/PolyDataMapper/PolyDataMapper.h>

#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkContourFilter.h"

using namespace std;
using namespace SCIRun;
using namespace vtk;

extern "C" vtk::Component* make_Vtk_PolyDataMapper()
{
  return new PolyDataMapper;
}

//Input Port
IPort::IPort(vtkPolyDataMapper *mapper){
  this->mapper=mapper;  
}

IPort::~IPort(){

}

bool
IPort::isInput(){
  return true;
}

std::string
IPort::getName(){
  return "PolyDataMapper::input";
}

bool 
IPort::accept(Port* port){
  return dynamic_cast<vtkPolyData*>(port->getObj())!=0;
}

void
IPort::connect(Port* port){
  mapper->SetInput(dynamic_cast<vtkPolyData*>(port->getObj()));
  mapper->ScalarVisibilityOff();
}


//Output Port

OPort::OPort(vtkPolyDataMapper *mapper){
  this->mapper=mapper;
}

OPort::~OPort(){

}

bool
OPort::isInput(){
  return false;
}

std::string
OPort::getName(){
  return "PolyDataMapper::output";
}

vtkObject *
OPort::getObj(){
  return mapper;
}


PolyDataMapper::PolyDataMapper(){

  mapper=vtkPolyDataMapper::New();
  iports.push_back(new IPort(mapper));
  oports.push_back(new OPort(mapper));
}

PolyDataMapper::~PolyDataMapper(){
  for(unsigned int i=0; i<iports.size(); i++){
    delete iports[i];
  }
  for(unsigned int i=0; i<oports.size(); i++){
    delete oports[i];
  }
  mapper->Delete();
}
