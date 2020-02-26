#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include "Node.h"
#include "Mesh.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

class BoundaryElement
{
public:
    BoundaryElement();
    BoundaryElement(const int &boundaryIndex, std::vector<Node *> connection); //, Mesh *mesh
    ~BoundaryElement();

    int getBoundaryIndex();

    std::vector<Node *> getNodes();

    Node *getNode(const int &index);

    matrix<double> shapeFunctionsAndDerivates(const double &xsi);

    vector<double> computeDistribuitedLoads(const bounded_vector<double, 2> &value, const int &quadraturePoints);

    matrix<double> boundaryIsoQuadrature(const int &points);

    //matrix<double> derivatesShapeFunctions(const double &xsi);

private:
    int boundaryIndex_; //

    //int lineNumber_;

    std::vector<Node *> connection_; //OK

    //Patch *patch_; //OK
};