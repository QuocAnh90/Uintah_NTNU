#Makefile fragment for the Packages/McQ/Core directory

include $(SRCTOP_ABS)/scripts/largeso_prologue.mk

SRCDIR := Packages/McQ/Core
SUBDIRS := \
	$(SRCDIR)/Algorithms \
	$(SRCDIR)/Datatypes \
	$(SRCDIR)/ThirdParty \

include $(SRCTOP_ABS)/scripts/recurse.mk

PSELIBS := 
LIBS := $(TK_LIBRARY) $(GL_LIBS) $(M_LIBRARY)

include $(SRCTOP_ABS)/scripts/largeso_epilogue.mk
