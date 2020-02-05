#pragma once

#include "../Isogeometric/ControlPoint.h"

class NeumannCondition
{
public:
    NeumannCondition();

    NeumannCondition(ControlPoint *point,
                       const int &direction,
                       const double &value);

    ~NeumannCondition();

    ControlPoint *getControlPoint();

    int getDirection();

    double getValue();

private:
    ControlPoint *point_;

    int direction_;

    double value_;
};