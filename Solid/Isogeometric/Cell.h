#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include "ControlPoint.h"
#include "../Material.h"
#include "Patch.h"
#include "../FiniteElement/BoundaryElement.h"
// #include "ShapeFunction.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

class Cell
{
public:
    Cell();

    Cell(const int &index,
         Patch *patch,
         const std::vector<ControlPoint *> &controlPoints);

    ~Cell();

    int getIndex();

    std::vector<ControlPoint *> getControlPoints();

    ControlPoint *getControlPoint(const int &index);

    bounded_matrix<double, 2, 2> referenceJacobianMatrix(const matrix<double> &dphi_dxsi, const vector<double> &wpc);

    bounded_matrix<double, 2, 2> currentJacobianMatrix(const matrix<double> &dphi_dxsi, const vector<double> &wpc);

    double jacobianDeterminant(const bounded_matrix<double, 2, 2> &jacobianMatrix);

    std::pair<vector<double>, matrix<double>> cellContributions(const std::string &ep, const std::string &typeAnalyze,
                                                                const int &step, const int &numberOfStep, const double &deltat,
                                                                const double &beta, const double &gama,
                                                                const int &pointsQuadrature);

    void setShapeForce(const bounded_vector<double, 2> &shapeForce);

    void StressCalculate(const std::string &ep);

    matrix<double> massMatrix(const int &pointsQuadrature);

    Patch *getPatch();

    matrix<double> isoQuadrature(const int &points);

    bounded_vector<double, 4> getCauchStress(const bounded_vector<double, 2> &qxsi, const std::string &ep);

    bounded_vector<double, 4> getGreen(const bounded_vector<double, 2> &qxsi, const std::string &ep);

    vector<double> shapeFunction(const bounded_vector<double, 2> &qxsi,
                                 const vector<double> &wpc,
                                 const bounded_vector<double, 2> inc);

    std::pair<vector<double>, matrix<double>> shapeFunctionAndDerivates(const bounded_vector<double, 2> &qxsi,
                                                                        const vector<double> &wpc,
                                                                        const bounded_vector<double, 2> inc);

    matrix<double> boundaryIsoQuadrature(const int &points);

    std::pair<vector<double>, vector<double>> boundaryShapeFunctionAndDerivates(const double &xsi,
                                                                        const vector<double> &wpc,
                                                                        const bounded_vector<double, 2> inc, const int &curveNumber);

    vector<double> computeDistribuitedLoads(const bounded_vector<double, 2> &value, const int &quadraturePoints, const int &curveNumber);

    std::vector<ControlPoint * > getControlPointsOnSide(const int &side);

    void computeDistanceFromFEBoundary(const int &pointsQuadrature, std::vector<BoundaryElement *> boundaryFE);

private:
    int index_; //cell number //OK

    bounded_vector<double, 2> shapeForce_;

    std::vector<ControlPoint *> controlPoints_; //OK

    Patch *patch_; //OK

    std::vector<double> distanceFE_;
};