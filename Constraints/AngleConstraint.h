
/*
 *  AngleConstraint.h
 *
 *  Written by:
 *   James Purciful
 *   Department of Computer Science
 *   University of Utah
 *   Aug. 1994
 *
 *  Copyright (C) 1994 SCI Group
 */


#ifndef SCI_project_Angle_Constraint_h
#define SCI_project_Angle_Constraint_h 1

#include <Constraints/BaseConstraint.h>


class AngleConstraint : public BaseConstraint {
public:
   AngleConstraint( const clString& name,
		    const Index numSchemes,
		    Variable* center, Variable* end1,
		    Variable* end2, Variable* p,
		    Variable* angleInX );
   virtual ~AngleConstraint();

protected:
   virtual void Satisfy( const Index index, const Scheme scheme );
};

#endif
