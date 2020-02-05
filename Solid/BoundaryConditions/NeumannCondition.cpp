#include "NeumannCondition.h"

NeumannCondition::NeumannCondition() {}

NeumannCondition::NeumannCondition(ControlPoint *point,
                                   const int &direction,
                                   const double &value)
{
    point_ = point;
    direction_ = direction;
    value_ = value;
}

NeumannCondition::~NeumannCondition() {}

ControlPoint *NeumannCondition::getControlPoint()
{
    return point_;
}

int NeumannCondition::getDirection()
{
    return direction_;
}

double NeumannCondition::getValue()
{
    return value_;
}