
#include <Packages/Uintah/Core/Datatypes/VectorParticlesPort.h>
#include <Packages/Uintah/share/share.h>
#include <Core/Malloc/Allocator.h>

using namespace SCIRun;
using namespace Uintah;

extern "C" {
  UINTAHSHARE IPort* make_VectorParticlesIPort(Module* module,
					       const string& name) {
    return scinew SimpleIPort<VectorParticlesHandle>(module,name);
  }
  UINTAHSHARE OPort* make_VectorParticlesOPort(Module* module,
					       const string& name) {
    return scinew SimpleOPort<VectorParticlesHandle>(module,name);
  }
}

template<> string SimpleIPort<VectorParticlesHandle>::port_type("VectorParticles");
template<> string SimpleIPort<VectorParticlesHandle>::port_color("chartreuse3");


