#pragma once

#include "../Isogeometric/ControlPoint.h"
#include "../FiniteElement/Node.h"

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

class NeumannConditionFE
{
public:
    NeumannConditionFE();

    NeumannConditionFE(Node *node,
                       const int &direction,
                       const double &value);

    ~NeumannConditionFE();

    Node *getNode();

    int getDirection();

    double getValue();

private:
    Node *node_;

    int direction_;

    double value_;
};