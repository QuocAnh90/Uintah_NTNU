//  ContourGeom.h - A base class for regular geometries with alligned axes
//
//  Written by:
//   Eric Kuehne
//   Department of Computer Science
//   University of Utah
//   April 2000
//
//  Copyright (C) 2000 SCI Institute


#ifndef SCI_project_ContourGeom_h
#define SCI_project_ContourGeom_h 1

#include <SCICore/Geometry/Vector.h>
#include <SCICore/Geometry/Point.h>
#include <SCICore/Geom/GeomTriangles.h>
#include <SCICore/Geom/GeomPolyline.h>
#include <SCICore/Datatypes/PointCloudGeom.h>
#include <SCICore/Containers/LockingHandle.h>
#include <SCICore/Math/MiscMath.h>
#include <SCICore/Util/DebugStream.h>
//#include <sstream>
#include <vector>
#include <string>
//#include <set>


namespace SCICore {
namespace Datatypes {

using SCICore::Geometry::Vector;
using SCICore::Geometry::Point;
using SCICore::Geometry::Min;
using SCICore::Geometry::Max;
using SCICore::GeomSpace::GeomTrianglesP;
using SCICore::GeomSpace::GeomPolyline;
using std::vector;
using std::string;
using SCICore::PersistentSpace::Piostream;
using SCICore::PersistentSpace::PersistentTypeID;
using SCICore::Math::Interpolate;
using SCICore::Util::DebugStream;


class ContourGeom : public PointCloudGeom
{
public:

  ContourGeom();
  ContourGeom(const vector<NodeSimp>&, const vector<EdgeSimp>&);
  ~ContourGeom();

  virtual string getInfo();
  
  //////////
  // Interpolate
  template <class A>
  int slinterpolate(A* att, elem_t, const Point& p, double& outval,
		    double eps=1.0e-6);

  //////////
  // Deletes these pointers if they are already set.
  void setEdges(const vector<EdgeSimp>&);

  ///////////
  // Persistent representation...
  virtual void io(Piostream&);
  static PersistentTypeID type_id;

protected:

  bool d_has_neighbors;
 
  vector<EdgeSimp> d_edge;

private:
  static DebugStream dbg;
};

template <class A>
int ContourGeom::slinterpolate(A* att, elem_t elem_type, const Point& p,
			       double& outval, double eps)
{
}

} // end Datatypes
} // end SCICore


#endif
