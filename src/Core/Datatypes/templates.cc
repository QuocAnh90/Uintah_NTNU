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
 * Manual template instantiations
 */


/*
 * These aren't used by Datatypes directly, but since they are used in
 * a lot of different modules, we instantiate them here to avoid bloat
 *
 * Find the bloaters with:
find . -name "*.ii" -print | xargs cat | sort | uniq -c | sort -nr | more
 */

#include <Core/Containers/LockingHandle.h>
#include <Core/Malloc/Allocator.h>

using namespace SCIRun;
#ifdef __sgi
#pragma set woff 1468
#endif


#include <Core/Datatypes/ColumnMatrix.h>
template class LockingHandle<ColumnMatrix>;

#include <Core/Datatypes/Matrix.h>
template class LockingHandle<Matrix>;

#include <Core/Datatypes/TetVol.h>
// linux needs these explicit declarations so that type_id is initialized.
template class TetVol<double>;
template class GenericField<TetVolMesh, vector<double> >;

#include <Core/Geometry/Tensor.h>
template class TetVol<Tensor>;
template class GenericField<TetVolMesh, vector<Tensor> >;

#include <Core/Datatypes/LatticeVol.h>
template class LatticeVol<double>;
template class GenericField<LatVolMesh, FData3d<double> >;

#include <Core/Geometry/Vector.h>
template class LatticeVol<Vector>;
template class GenericField<LatVolMesh, FData3d<Vector> >;

#include <Core/Datatypes/ContourField.h>
template class ContourField<double>;
template class GenericField<ContourMesh, vector<double> >;

#include <Core/Datatypes/GenericField.h>
#include <Core/Persistent/PersistentSTL.h>

#include <Core/Datatypes/TriSurfMesh.h>
template class GenericField<TriSurfMesh, vector<double> >;

#include <Core/Datatypes/PropertyManager.h>
template class Property<string>;
template class Property<pair<double,double> >;
template class Property<Array1<double> >;
template class Property<Array1<Tensor> >;
template class Property<pair<int,double> >;


//! Compute the gradient g in cell ci.
template <>
Vector TetVol<Vector>::cell_gradient(TetVolMesh::cell_index ci)
{
  ASSERT(type_name(1) != "Vector");  // redundant, useful error message
  return Vector(0, 0, 0);
}


template <>
Vector TetVol<Tensor>::cell_gradient(TetVolMesh::cell_index ci)
{
  ASSERT(type_name(1) != "Tensor");  // redundant, useful error message
  return Vector(0, 0, 0);
}

#ifdef __sgi
#pragma reset woff 1468
#endif










