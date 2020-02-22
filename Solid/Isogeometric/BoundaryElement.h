#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include "ControlPoint.h"
#include "../Material.h"
#include "Patch.h"
// #include "ShapeFunction.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

class BoundaryElement
{
public:
    BoundaryElement();
    BoundaryElement(const int &index, Patch *patch, std::vector<ControlPoint *> controlPoints, const int &curveNumber);
    ~BoundaryElement();

    std::vector<ControlPoint *> getControlPoints();

    ControlPoint *getControlPoint(const int &index);

    matrix<double> isoQuadrature(const int &points);

    std::pair<vector<double>, vector<double>> shapeFunctionAndDerivates(const double &xsi,
                                                                        const vector<double> &wpc,
                                                                        const bounded_vector<double, 2> inc);

    vector<double> computeDistribuitedLoads(const bounded_vector<double, 2> &value, const int &quadraturePoints);

private:
    int index_; //

    int cunverNumber_;

    std::vector<ControlPoint *> controlPoints_; //OK

    Patch *patch_; //OK
};