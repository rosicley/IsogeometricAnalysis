#pragma once
#include "Element.h"

class CP_BlendingZone //control point inside blending zone
{
public:
    CP_BlendingZone();

    CP_BlendingZone(ControlPoint *cp, Element* el, const bounded_vector<double, 2> &xsi);

    ~CP_BlendingZone();

    ControlPoint *getControlPoint();

    Element *getElement();

    bounded_vector<double, 2> getDimensionlessCoordinates();

    void interpolateGlobalCoordinate();

private:
    ControlPoint* controlPoint_;

    Element* element_;

    bounded_vector<double, 2> xsi_;

};

