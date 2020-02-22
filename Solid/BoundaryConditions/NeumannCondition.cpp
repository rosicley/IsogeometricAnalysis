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


NeumannConditionFE::NeumannConditionFE() {}

NeumannConditionFE::NeumannConditionFE(Node *node,
                                   const int &direction,
                                   const double &value)
{
    node_ = node;
    direction_ = direction;
    value_ = value;
}

NeumannConditionFE::~NeumannConditionFE() {}

Node *NeumannConditionFE::getNode()
{
    return node_;
}

int NeumannConditionFE::getDirection()
{
    return direction_;
}

double NeumannConditionFE::getValue()
{
    return value_;
}