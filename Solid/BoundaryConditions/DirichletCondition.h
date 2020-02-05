#pragma once

#include "../Isogeometric/ControlPoint.h"

class DirichletCondition
{
public:
    DirichletCondition();

    DirichletCondition(ControlPoint *point,
                       const int &direction,
                       const double &value);

    ~DirichletCondition();

    ControlPoint *getControlPoint();

    int getDirection();

    double getValue();

private:
    ControlPoint *point_;

    int direction_;

    double value_;
};