#pragma once
#include "Element.h"

class InactiveCP //control point inside blending zone
{
public:
    InactiveCP();

    InactiveCP(ControlPoint *cp, Element* el, const bounded_vector<double, 2> &xsi);

    ~InactiveCP();

    ControlPoint *getControlPoint();

    Element *getElement();

    bounded_vector<double, 2> getDimensionlessCoordinates();

    void interpolateGlobalCoordinate();

private:
    ControlPoint* controlPoint_;

    Element* element_;

    bounded_vector<double, 2> xsi_;

};

class InactiveNode //control point inside blending zone
{
public:
    InactiveNode();

    InactiveNode(Node *no, Cell* cell, const bounded_vector<double, 2> &xsi);

    ~InactiveNode();

    Node *getNode();

    Cell *getCell();

    bounded_vector<double, 2> getDimensionlessCoordinates();

    void interpolateGlobalCoordinate();

private:
    Node* node_;

    Cell* cell_;

    bounded_vector<double, 2> xsi_;

};

