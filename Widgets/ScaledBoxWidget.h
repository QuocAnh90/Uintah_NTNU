
/*
 *  ScaledBoxWidget.h
 *
 *  Written by:
 *   James Purciful
 *   Department of Computer Science
 *   University of Utah
 *   Jan. 1995
 *
 *  Copyright (C) 1995 SCI Group
 */


#ifndef SCI_project_ScaledBox_Widget_h
#define SCI_project_ScaledBox_Widget_h 1

#include <Widgets/BaseWidget.h>


class ScaledBoxWidget : public BaseWidget {
public:
   ScaledBoxWidget( Module* module, CrowdMonitor* lock, double widget_scale );
   ScaledBoxWidget( const ScaledBoxWidget& );
   virtual ~ScaledBoxWidget();

   virtual void widget_execute();
   virtual void geom_moved(int, double, const Vector&, int);

   virtual void MoveDelta( const Vector& delta );
   virtual Point ReferencePoint() const;

   Real GetRatioR() const;
   Real GetRatioD() const;
   Real GetRatioI() const;

   const Vector& GetRightAxis();
   const Vector& GetDownAxis();
   const Vector& GetInAxis();

   // Variable indexs
   enum { CenterVar, PointRVar, PointDVar, PointIVar,
	  DistRVar, DistDVar, DistIVar, HypoRDVar, HypoDIVar, HypoIRVar,
	  SDistRVar, RatioRVar, SDistDVar, RatioDVar, SDistIVar, RatioIVar };

   // Material indexs
   enum { PointMatl, EdgeMatl, ResizeMatl, SliderMatl };

private:
   Vector oldrightaxis, olddownaxis, oldinaxis;
};


#endif
