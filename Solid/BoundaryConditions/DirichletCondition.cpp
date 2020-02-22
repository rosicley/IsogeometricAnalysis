#include "DirichletCondition.h"

DirichletCondition::DirichletCondition() {}

DirichletCondition::DirichletCondition(ControlPoint *point,
                                       const int &direction,
                                       const double &value)
{
    point_ = point;
    direction_ = direction;
    value_ = value;
}

DirichletCondition::~DirichletCondition() {}

ControlPoint *DirichletCondition::getControlPoint()
{
    return point_;
}

int DirichletCondition::getDirection()
{
    return direction_;
}

double DirichletCondition::getValue()
{
    return value_;
}


DirichletConditionFE::DirichletConditionFE() {}

DirichletConditionFE::DirichletConditionFE(Node *node,
                                       const int &direction,
                                       const double &value)
{
    node_ = node;
    direction_ = direction;
    value_ = value;
}

DirichletConditionFE::~DirichletConditionFE() {}

Node *DirichletConditionFE::getNode()
{
    return node_;
}

int DirichletConditionFE::getDirection()
{
    return direction_;
}

double DirichletConditionFE::getValue()
{
    return value_;
}