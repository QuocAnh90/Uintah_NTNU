
/*
 *  SoundMixer.h: The 
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#ifndef SCI_project_module_SoundMixer_h
#define SCI_project_module_SoundMixer_h

#include <UserModule.h>
class MUI_slider_real;
class SoundIPort;
class SoundOPort;

class SoundMixer : public UserModule {
    struct PortInfo {
	double gain;
	MUI_slider_real* interface;
	SoundIPort* isound;
    };
    Array1<PortInfo*> portinfo;
    double overall_gain;
    SoundOPort* osound;
public:
    SoundMixer();
    SoundMixer(const SoundMixer&, int deep);
    virtual ~SoundMixer();
    virtual Module* clone(int deep);
    virtual void connection(Module::ConnectionMode, int, int);
    virtual void execute();
};

#endif
