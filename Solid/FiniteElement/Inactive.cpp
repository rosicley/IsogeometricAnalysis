#include "Inactive.h"

InactiveCP::InactiveCP() {}

InactiveCP::InactiveCP(ControlPoint *cp, Element *el, const bounded_vector<double, 2> &xsi)
{
    controlPoint_ = cp;
    element_ = el;
    xsi_ = xsi;
}

InactiveCP::~InactiveCP() {}

void InactiveCP::interpolateGlobalCoordinate()
{
    bounded_vector<double, 2> coord;
    coord = element_->calculateGlobalCoordinate(xsi_);
    controlPoint_->setCurrentCoordinate(coord);
}

InactiveNode::InactiveNode() {}

InactiveNode::InactiveNode(Node *no, Cell *cell, const bounded_vector<double, 2> &xsi)
{
    node_ = no;
    cell_ = cell;
    xsi_ = xsi;
}

InactiveNode::~InactiveNode() {}

void InactiveNode::interpolateGlobalCoordinate()
{
    bounded_vector<double, 2> coord;
    coord = cell_->calculateGlobalCoordinate(xsi_);
    node_->setCurrentCoordinate(coord);
}