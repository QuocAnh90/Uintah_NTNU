
/*
 *  ScaledFrameWidget.h
 *
 *  Written by:
 *   James Purciful
 *   Department of Computer Science
 *   University of Utah
 *   Jan. 1995
 *
 *  Copyright (C) 1995 SCI Group
 */


#ifndef SCI_project_ScaledFrame_Widget_h
#define SCI_project_ScaledFrame_Widget_h 1

#include <Widgets/BaseWidget.h>


// Variable indexs
enum { SFrameW_PointUL, SFrameW_PointUR, SFrameW_PointDR, SFrameW_PointDL,
       SFrameW_Dist1, SFrameW_Dist2, SFrameW_Hypo,
       SFrameW_Slider1, SFrameW_SDist1, SFrameW_Ratio1,
       SFrameW_Slider2, SFrameW_SDist2, SFrameW_Ratio2 };
// Material indexs
enum { SFrameW_PointMatl, SFrameW_EdgeMatl, SFrameW_SliderMatl, SFrameW_HighMatl };


class ScaledFrameWidget : public BaseWidget {
public:
   ScaledFrameWidget( Module* module, CrowdMonitor* lock, double widget_scale );
   ScaledFrameWidget( const ScaledFrameWidget& );
   ~ScaledFrameWidget();

   virtual void widget_execute();
   virtual void geom_moved(int, double, const Vector&, void*);

   inline Real GetRatio1() const;
   inline Real GetRatio2() const;
   inline Vector GetAxis1();
   inline Vector GetAxis2();

private:
   Vector oldaxis1, oldaxis2;
};


inline Real
ScaledFrameWidget::GetRatio1() const
{
   return (variables[SFrameW_Ratio1]->Get().x());
}


inline Real
ScaledFrameWidget::GetRatio2() const
{
   return (variables[SFrameW_Ratio2]->Get().x());
}


inline Vector
ScaledFrameWidget::GetAxis1()
{
   Vector axis(variables[SFrameW_PointUR]->Get() - variables[SFrameW_PointUL]->Get());
   if (axis.length2() <= 1e-6)
      return oldaxis1;
   else
      return (oldaxis1 = axis.normal());
}


inline Vector
ScaledFrameWidget::GetAxis2()
{
   Vector axis(variables[SFrameW_PointDL]->Get() - variables[SFrameW_PointUL]->Get());
   if (axis.length2() <= 1e-6)
      return oldaxis2;
   else
      return (oldaxis2 = axis.normal());
}


#endif
