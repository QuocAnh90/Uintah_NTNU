//static char *id="@(#) $Id"

/*
 *  ImageToGeom.cc:  Unfinished modules
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <Util/NotFinished.h>
#include <Dataflow/Module.h>
#include <CommonDatatypes/GeometryPort.h>
#include <Datatypes/Image/ImagePort.h>
#include <Geom/GeomGrid.h>
#include <Malloc/Allocator.h>
#include <TclInterface/TCLvar.h>

namespace SCIRun {
namespace Modules {

using namespace SCICore::TclInterface;

using namespace PSECommon::Dataflow;
using namespace PSECommon::CommonDatatypes;

using namespace SCIRun::Datatypes;

class ImageToGeom : public Module {
    ImageIPort* iport;
    GeometryOPort* ogeom;
    double scale;
    int have_scale;
public:
    TCLint format;
    TCLdouble heightscale;
    TCLint downsample;
    ImageToGeom(const clString& id);
    ImageToGeom(const ImageToGeom& copy, int deep);
    virtual ~ImageToGeom();
    virtual Module* clone(int deep);
    virtual void execute();

    int oldid;
};

extern "C" {
  Module* make_ImageToGeom(const clString& id)
    {
      return scinew ImageToGeom(id);
    }
}

ImageToGeom::ImageToGeom(const clString& id)
: Module("ImageToGeom", id, Filter), format("format", id, this),
  heightscale("heightscale", id, this), downsample("downsample", id, this)
{
    iport=scinew ImageIPort(this, "Image", ImageIPort::Atomic);
    add_iport(iport);
    // Create the output port
    ogeom=scinew GeometryOPort(this, "Geom object", ImageIPort::Atomic);
    add_oport(ogeom);
    oldid=0;
    have_scale=0;
}

ImageToGeom::ImageToGeom(const ImageToGeom& copy, int deep)
: Module(copy, deep), format("format", id, this),
  heightscale("heightscale", id, this), downsample("downsample", id, this)
{
    NOT_FINISHED("ImageToGeom::ImageToGeom");
}

ImageToGeom::~ImageToGeom()
{
}

Module* ImageToGeom::clone(int deep)
{
    return scinew ImageToGeom(*this, deep);
}

void ImageToGeom::execute()
{
    ImageHandle image;
    if(!iport->get(image))
	return;
    //cerr << "generation=" << image->generation << endl;
    int form=format.get();
    GeomObj* obj=0;
    double hs(heightscale.get());
    switch(form){
    case 0:
	{
	    int xres=image->xres();
	    int yres=image->yres();
#if 0
	    GeomGrid* grid=new GeomGrid(image->xres(), image->yres(),
					Point(0,0,0),
					Vector(xres-1,0,0),
					Vector(0,yres-1,0),
					GeomGrid::WithNormals);
	    for(int j=0;j<yres;j++){
		for(int i=0;i<xres;i++){
		    float value=image->getr(i,j);
		    grid->setn(i,j,hs*value);
		}
	    }
	    grid->compute_normals();
	    obj=grid;
#endif
	}
	break;
    default:
	error("Unknown format");
	break;
    }

    if(oldid != 0)
	ogeom->delObj(oldid);
    if(obj != 0){	
	oldid=ogeom->addObj(obj, "Image");
    } else {
	oldid=0;
    }
    ogeom->flushViews();
}

} // End namespace Modules
} // End namespace SCIRun

//
// $Log$
// Revision 1.1  1999/07/27 16:58:53  mcq
// Initial commit
//
// Revision 1.2  1999/04/30 01:11:53  dav
// moved TiffReader to SCIRun from PSECommon
//
// Revision 1.1  1999/04/29 22:26:33  dav
// Added image files to SCIRun
//
//
