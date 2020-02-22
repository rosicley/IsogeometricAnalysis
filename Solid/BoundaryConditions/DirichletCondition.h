#pragma once

#include "../Isogeometric/ControlPoint.h"
#include "../FiniteElement/Node.h"


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

class DirichletConditionFE
{
public:
    DirichletConditionFE();

    DirichletConditionFE(Node *node,
                       const int &direction,
                       const double &value);

    ~DirichletConditionFE();

    Node *getNode();

    int getDirection();

    double getValue();

private:
    Node *node_;

    int direction_;

    double value_;
};