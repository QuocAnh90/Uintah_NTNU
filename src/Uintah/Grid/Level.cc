/* REFERENCED */
static char *id="@(#) $Id$";

#include "Level.h"
#include "Handle.h"
#include "Region.h"

#include <iostream>
using SCICore::Geometry::Point;
using std::cerr;

namespace Uintah {
namespace Grid {

Level::Level()
{
}

Level::~Level()
{
  // Delete all of the regions managed by this level
  for(regionIterator iter=d_regions.begin(); iter != d_regions.end(); iter++)
    delete *iter;
}

Level::const_regionIterator Level::regionsBegin() const
{
    return d_regions.begin();
}

Level::const_regionIterator Level::regionsEnd() const
{
    return d_regions.end();
}

Level::regionIterator Level::regionsBegin()
{
    return d_regions.begin();
}

Level::regionIterator Level::regionsEnd()
{
    return d_regions.end();
}

Region* Level::addRegion(const Point& lower, const Point& upper,
			 const IntVector& res)
{
    Region* r = new Region(lower, upper, res);
    d_regions.push_back(r);
    return r;
}

void Level::performConsistencyCheck() const
{
}

} // end namespace Grid
} // end namespace Uintah

//
// $Log$
// Revision 1.3  2000/04/12 23:00:47  sparker
// Starting problem setup code
// Other compilation fixes
//
// Revision 1.2  2000/03/16 22:07:59  dav
// Added the beginnings of cocoon docs.  Added namespaces.  Did a few other coding standards updates too
//
//
