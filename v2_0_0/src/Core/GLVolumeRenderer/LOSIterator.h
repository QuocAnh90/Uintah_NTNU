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

#ifndef LOSITERATOR_H
#define LOSITERATOR_H

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Ray.h>
#include <Core/GLVolumeRenderer/GLTextureIterator.h>
#include <Core/GLVolumeRenderer/GLTexture3D.h>

namespace SCIRun {


/**************************************

CLASS
   LOSIterator
   
   Line Of Site IteratorClass.

GENERAL INFORMATION

   LOSIterator.h

   Kurt Zimmerman
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   Copyright (C) 2000 SCI Group

KEYWORDS
   Texture

DESCRIPTION
   Line Of Site Iterator  class.  If we have a multiresolution texture,
   we may only want to resolve it along some path.  This iterator will
   step through the Bricks and return those along the path at full
   resolution.  Textures further away from the line of sight will get
   rendered at lower resolution.
  
WARNING
  
****************************************/

class LOSIterator : public GLTextureIterator {
public:
  // GROUP: Constructors:
  //////////
  // Constructor
  LOSIterator(const GLTexture3D* tex,  Ray view,
	      Point  control);
  // GROUP: Destructors
  //////////
  // Destructor
  ~LOSIterator(){}
 
  // GROUP: Access
  //////////
  // get first brick
  virtual Brick* Start();
  //////////
  // get next brick
  virtual Brick* Next();
  // GROUP: Query
  //////////
  // are we finished?
  virtual bool isDone();

protected:
  void SetNext();
private:

};

} // End namespace SCIRun
#endif
