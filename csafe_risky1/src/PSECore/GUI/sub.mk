#
# Makefile fragment for this subdirectory
# $Id$
#

SRCDIR := PSECore/GUI

ALLTARGETS := $(ALLTARGETS) $(SRCDIR)/tclIndex

$(SRCDIR)/tclIndex: $(SRCDIR)/ArrowWidget.tcl $(SRCDIR)/BaseWidget.tcl \
	$(SRCDIR)/BoxWidget.tcl $(SRCDIR)/CriticalPointWidget.tcl \
	$(SRCDIR)/CrosshairWidget.tcl $(SRCDIR)/FrameWidget.tcl \
	$(SRCDIR)/GaugeWidget.tcl $(SRCDIR)/LightWidget.tcl \
	$(SRCDIR)/Module.tcl $(SRCDIR)/NetworkEditor.tcl \
	$(SRCDIR)/PathWidget.tcl $(SRCDIR)/PointWidget.tcl \
	$(SRCDIR)/RingWidget.tcl $(SRCDIR)/ScaledBoxWidget.tcl \
	$(SRCDIR)/ScaledFrameWidget.tcl $(SRCDIR)/ViewWidget.tcl \
	$(SRCDIR)/defaults.tcl $(SRCDIR)/devices.tcl \
	$(SRCDIR)/platformSpecific.tcl
	$(SRCTOP)/scripts/createTclIndex $(SRCTOP)/PSECore/GUI

CLEANPROGS := $(CLEANPROGS) $(SRCDIR)/tclIndex

#
# $Log$
# Revision 1.2  2000/03/20 19:37:23  sparker
# Added VPATH support
#
# Revision 1.1  2000/03/17 09:28:02  sparker
# New makefile scheme: sub.mk instead of Makefile.in
# Use XML-based files for module repository
# Plus many other changes to make these two things work
#
#
