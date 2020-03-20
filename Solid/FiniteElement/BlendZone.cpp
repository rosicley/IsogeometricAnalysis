#include "BlendZone.h"

CP_BlendingZone::CP_BlendingZone() {}

CP_BlendingZone::CP_BlendingZone(ControlPoint *cp, Element *el, const bounded_vector<double, 2> &xsi)
{
    controlPoint_ = cp;
    element_ = el;
    xsi_ = xsi;
}

CP_BlendingZone::~CP_BlendingZone() {}

ControlPoint *CP_BlendingZone::getControlPoint()
{
    return controlPoint_;
}

Element *CP_BlendingZone::getElement()
{
   return element_; 
}

bounded_vector<double, 2> CP_BlendingZone::getDimensionlessCoordinates()
{
    return xsi_;
}

void CP_BlendingZone::interpolateGlobalCoordinate()
{
    bounded_vector<double, 2> coord;
    coord = element_ ->calculateGlobalCoordinate(xsi_);
    controlPoint_->setCurrentCoordinate(coord);
}