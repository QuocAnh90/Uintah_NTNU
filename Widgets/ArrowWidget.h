
/*
 *  ArrowWidget.h
 *
 *  Written by:
 *   James Purciful
 *   Department of Computer Science
 *   University of Utah
 *   Aug. 1994
 *
 *  Copyright (C) 1994 SCI Group
 */


#ifndef SCI_project_Arrow_Widget_h
#define SCI_project_Arrow_Widget_h 1

#include <Widgets/BaseWidget.h>


// Variable indexs
enum { ArrowW_Point };
// Material indexs
enum { ArrowW_PointMatl, ArrowW_EdgeMatl, ArrowW_HighMatl };


class ArrowWidget : public BaseWidget {
public:
   ArrowWidget( Module* module, CrowdMonitor* lock, double widget_scale );
   ArrowWidget( const ArrowWidget& );
   ~ArrowWidget();

   virtual void widget_execute();
   virtual void geom_moved(int, double, const Vector&, void*);

   void SetPosition( const Point& );
   const Point& GetPosition() const;
   
   void SetDirect( const Vector& v );
   const Vector& GetDirect() const;

private:
   Vector direction;
};


#endif
