//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2004 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : HexTriquadraticLgn.h
//    Author : Martin Cole, Frank B. Sachse
//    Date   : Dec 3 2004

#if !defined(HexTriquadraticLgn_h)
#define HexTriquadraticLgn_h

#include <Core/Basis/HexTrilinearLgn.h>

namespace SCIRun {

//! Class for describing unit geometry of HexTriquadraticLgn 
class HexTriquadraticLgnUnitElement {
public:
  static double unit_vertices[20][3]; //!< Parametric coordinates of vertices of unit edge 
  static int unit_edges[12][2];  //!< References to vertices of unit edge  
  static int unit_faces[6][4];  //!< References to vertices of unit face 
 
  HexTriquadraticLgnUnitElement() {};
  virtual ~HexTriquadraticLgnUnitElement() {};
  
  static int domain_dimension() { return 3; }; //! return dimension of domain 
  
  static int number_of_vertices() { return 8; }; //! return number of vertices
  static int number_of_edges() { return 12; }; //! return number of edges
  
  static int vertices_of_face() { return 4; }; //! return number of vertices per face 

  static int faces_of_cell() { return 6; }; //! return number of faces per cell 
};


//! Class for handling of element of type hexahedron with 
//! triquadratic lagrangian interpolation
template <class T>
class HexTriquadraticLgn : public HexApprox, 
			   public HexGaussian3<double>, 
			   public HexTriquadraticLgnUnitElement
{
public:
  typedef T value_type;

  HexTriquadraticLgn() {}
  virtual ~HexTriquadraticLgn() {}

  int polynomial_order() const { return 2; }

  //! get weight factors at parametric coordinate 
  inline
  int get_weights(const vector<double> &coords, double *w) const
  {
    const double x=coords[0], y=coords[1], z=coords[2];
    w[0]  = (-1 + x)*(-1 + y)*(-1 + z)*(-1 + 2*x + 2*y + 2*z);
    w[1]  = +x*(-1 + y)*(-1 + 2*x - 2*y - 2*z)*(-1 + z);
    w[2]  = -(x*y*(-3 + 2*x + 2*y - 2*z)*(-1 + z));
    w[3]  = -((-1 + x)*y*(-1 + z)*(1 + 2*x - 2*y + 2*z));
    w[4]  = -((-1 + x)*(-1 + y)*(1 + 2*x + 2*y - 2*z)*z);
    w[5]  = +x*(-1 + y)*(3 - 2*x + 2*y - 2*z)*z;
    w[6]  = +x*y*z*(-5 + 2*x + 2*y + 2*z);
    w[7]  = +(-1 + x)*y*(3 + 2*x - 2*y - 2*z)*z;
    w[8]  = -4*(-1 + x)*x*(-1 + y)*(-1 + z);
    w[9]  = +4*x*(-1 +y)*y*(-1 + z);
    w[10] = +4*(-1 + x)*x*y*(-1 + z);
    w[11] = -4*(-1 + x)*(-1 + y)*y*(-1 + z);
    w[12] = -4*(-1 + x)*(-1 + y)*(-1 + z)*z;
    w[13] = +4*x*(-1 + y)*(-1 + z)*z;
    w[14] = +4*(-1 + x)*y*(-1 + z)*z;
    w[15] = -4*x*y*(-1 + z)*z;
    w[16] = +4*(-1 + x)*x*(-1 + y)*z;
    w[17] = -4*x*(-1 + y)*y*z;
    w[18] = -4*(-1 + x)*x*y*z;
    w[19] = +4*(-1 + x)*(-1 + y)*y*z;
    return 20;
  }
  //! get value at parametric coordinate 
  template <class ElemData>
  T interpolate(const vector<double> &coords, const ElemData &cd) const
  {
    double w[20];
    get_weights(coords, w); 
    return (T)(w[0]  * cd.node0() +
	       w[1]  * cd.node1() +
	       w[2]  * cd.node2() +
	       w[3]  * cd.node3() +
	       w[4]  * cd.node4() +
	       w[5]  * cd.node5() +
	       w[6]  * cd.node6() +
	       w[7]  * cd.node7() +
	       w[8]  * nodes_[cd.edge0_index()] +
	       w[9]  * nodes_[cd.edge1_index()] +
	       w[10] * nodes_[cd.edge2_index()] +
	       w[11] * nodes_[cd.edge3_index()] +
	       w[12] * nodes_[cd.edge4_index()] +
	       w[13] * nodes_[cd.edge5_index()] +
	       w[14] * nodes_[cd.edge6_index()] +
	       w[15] * nodes_[cd.edge7_index()] +
	       w[16] * nodes_[cd.edge8_index()] +
	       w[17] * nodes_[cd.edge9_index()] +
	       w[18] * nodes_[cd.edge10_index()] +
	       w[19] * nodes_[cd.edge11_index()]);
  }
  
  //! get first derivative at parametric coordinate
  template <class ElemData>
  void derivate(const vector<double> &coords, const ElemData &cd, 
		vector<T> &derivs) const
  {
    const double x=coords[0], y=coords[1], z=coords[2];

    derivs.resize(3);

    derivs[0]=T((-1 + y)*(-1 + z)*(-3 + 4*x + 2*y + 2*z)*cd.node0()
		-((-1 + y)*(-1 + z)*(1 - 4*x + 2*y + 2*z))*cd.node1()
		-(y*(-3 + 4*x + 2*y - 2*z)*(-1 + z))*cd.node2()
		+y*(1 - 4*x + 2*y - 2*z)*(-1 + z)*cd.node3()
		-((-1 + y)*(-1 + 4*x + 2*y - 2*z)*z)*cd.node4()
		+(-1 + y)*(3 - 4*x + 2*y - 2*z)*z*cd.node5()
		+y*z*(-5 + 4*x + 2*y + 2*z)*cd.node6()
		+y*(1 + 4*x - 2*y - 2*z)*z*cd.node7()
		-4*(-1 + 2*x)*(-1 + y)*(-1 + z)*nodes_[cd.edge0_index()]
		+4*(-1 + y)*y*(-1 + z)*nodes_[cd.edge1_index()]
		+4*(-1 + 2*x)*y*(-1 + z)*nodes_[cd.edge2_index()]
		+4*y*(-1 + y + z - y*z)*nodes_[cd.edge3_index()]
		+4*z*(-1 + y + z - y*z)*nodes_[cd.edge4_index()]
		+4*(-1 + y)*(-1 + z)*z*nodes_[cd.edge5_index()]
		+4*y*(-1 + z)*z*nodes_[cd.edge6_index()]
		-4*y*(-1 + z)*z*nodes_[cd.edge7_index()]
		+4*(-1 + 2*x)*(-1 + y)*z*nodes_[cd.edge8_index()]
		-4*(-1 + y)*y*z*nodes_[cd.edge9_index()]
		+4*(1 - 2*x)*y*z*nodes_[cd.edge10_index()]
		+4*(-1 + y)*y*z*nodes_[cd.edge11_index()]);
    
    derivs[1]=T((-1 + x)*(-1 + z)*(-3 + 2*x + 4*y + 2*z)*cd.node0()
		+x*(1 + 2*x - 4*y - 2*z)*(-1 + z)*cd.node1()
		-(x*(-3 + 2*x + 4*y - 2*z)*(-1 + z))*cd.node2()
		-((-1 + x)*(-1 + z)*(1 + 2*x - 4*y + 2*z))*cd.node3()
		-((-1 + x)*(-1 + 2*x + 4*y - 2*z)*z)*cd.node4()
		+x*(1 - 2*x + 4*y - 2*z)*z*cd.node5()
		+x*z*(-5 + 2*x + 4*y + 2*z)*cd.node6()
		+(-1 + x)*(3 + 2*x - 4*y - 2*z)*z*cd.node7()
		+4*x*(-1 + x + z - x*z)*nodes_[cd.edge0_index()]
		+4*x*(-1 + 2*y)*(-1 + z)*nodes_[cd.edge1_index()]
		+4*(-1 + x)*x*(-1 + z)*nodes_[cd.edge2_index()]
		-4*(-1 + x)*(-1 + 2*y)*(-1 + z)*nodes_[cd.edge3_index()]
		+4*z*(-1 + x + z - x*z)*nodes_[cd.edge4_index()]
		+4*x*(-1 + z)*z*nodes_[cd.edge5_index()]
		+4*(-1 + x)*(-1 + z)*z*nodes_[cd.edge6_index()]
		-4*x*(-1 + z)*z*nodes_[cd.edge7_index()]
		+4*(-1 + x)*x*z*nodes_[cd.edge8_index()]
		+4*x*(1 - 2*y)*z*nodes_[cd.edge9_index()]
		-4*(-1 + x)*x*z*nodes_[cd.edge10_index()]
		+4*(-1 + x)*(-1 + 2*y)*z*nodes_[cd.edge11_index()]);
    
    derivs[2]=T((-1 + x)*(-1 + y)*(-3 + 2*x + 2*y + 4*z)*cd.node0()
		+x*(-1 + y)*(1 + 2*x - 2*y - 4*z)*cd.node1()
		+x*y*(1 - 2*x - 2*y + 4*z)*cd.node2()
		-((-1 + x)*y*(-1 + 2*x - 2*y + 4*z))*cd.node3()
		-((-1 + x)*(-1 + y)*(1 + 2*x + 2*y - 4*z))*cd.node4()
		+x*(-1 + y)*(3 - 2*x + 2*y - 4*z)*cd.node5()
		+x*y*(-5 + 2*x + 2*y + 4*z)*cd.node6()
 		+(-1 + x)*y*(3 + 2*x - 2*y - 4*z)*cd.node7()
		+4*x*(-1 + x + y - x*y)*nodes_[cd.edge0_index()]
		+4*x*(-1 + y)*y*nodes_[cd.edge1_index()]
		+4*(-1 + x)*x*y*nodes_[cd.edge2_index()]
		+4*y*(-1 + x + y - x*y)*nodes_[cd.edge3_index()]
		-4*(-1 + x)*(-1 + y)*(-1 + 2*z)*nodes_[cd.edge4_index()]
		+4*x*(-1 + y)*(-1 + 2*z)*nodes_[cd.edge5_index()]
		+4*(-1 + x)*y*(-1 + 2*z)*nodes_[cd.edge6_index()]
		+4*x*y*(1 - 2*z)*nodes_[cd.edge7_index()]
		+4*(-1 + x)*x*(-1 + y)*nodes_[cd.edge8_index()]
		-4*x*(-1 + y)*y*nodes_[cd.edge9_index()]
		-4*(-1 + x)*x*y*nodes_[cd.edge10_index()]
		+4*(-1 + x)*(-1 + y)*y*nodes_[cd.edge11_index()]); 
    
  }

  //! get parametric coordinate for value within the element
  template <class ElemData>
  bool get_coords(vector<double> &coords, const T& value, 
		  const ElemData &cd) const  
  {
    HexLocate< HexTriquadraticLgn<T> > CL;
    return CL.get_coords(this, coords, value, cd);
  };
    
  //! add a node value corresponding to edge
  void add_node_value(const T &p) { nodes_.push_back(p); }
 
  static  const string type_name(int n = -1);

  virtual void io (Piostream& str);

protected:
  //! Additional support values.
  //! Triquadratic Lagrangian only needs additional nodes stored for each edge
  //! in the topology.
  vector<T>          nodes_; 
};


template <class T>
const TypeDescription* get_type_description(HexTriquadraticLgn<T> *)
{
  static TypeDescription* td = 0;
  if(!td){
    const TypeDescription *sub = get_type_description((T*)0);
    TypeDescription::td_vec *subs = scinew TypeDescription::td_vec(1);
    (*subs)[0] = sub;
    td = scinew TypeDescription(HexTriquadraticLgn<T>::type_name(0), subs, 
				string(__FILE__),
				"SCIRun", 
				TypeDescription::BASIS_E);
  }
  return td;
}

template <class T>
const string
HexTriquadraticLgn<T>::type_name(int n)
{
  ASSERT((n >= -1) && n <= 1);
  if (n == -1)
  {
    static const string name = type_name(0) + FTNS + type_name(1) + FTNE;
    return name;
  }
  else if (n == 0)
  {
    static const string nm("HexTriquadraticLgn");
    return nm;
  } else {
    return find_type_name((T *)0);
  }
}


const int HEXTRIQUADRATICLGN_VERSION = 1;
template <class T>
void
HexTriquadraticLgn<T>::io(Piostream &stream)
{
  stream.begin_class(type_name(-1), HEXTRIQUADRATICLGN_VERSION);
  Pio(stream, nodes_);
  stream.end_class();
}

} //namespace SCIRun

#endif // HexTriquadraticLgn_h
