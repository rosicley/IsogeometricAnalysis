#pragma once
#include <iostream>
#include <vector>
#include <string>
#include "Node.h"
#include "Mesh.h"
#include "BoundaryElement.h"
#include "../Material.h"
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "../Isogeometric/Cell.h"
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>

using namespace boost::numeric::ublas;

class Element
{
public:
    Element();

    Element(const int &index,
            Mesh *mesh,
            const std::vector<Node *> &connection);

    ~Element();

    int getIndex();

    std::vector<Node *> getConnection();

    Node *getNode(int index);

    Mesh *getMesh();

    std::vector<Node *> getSideNodes(const int &side);

    //Material *getMaterial();

    vector<double> domainShapeFunction(const double &xsi1, const double &xsi2);

    matrix<double> domainDerivativeShapeFunction(const double &xsi1, const double &xsi2);

    matrix<double> domainSecondDerivativeShapeFunction(const double &xsi1, const double &xsi2);

    matrix<double> hammerQuadrature(const int &nH);

    bounded_matrix<double, 2, 2> referenceJacobianMatrix(const double &xsi1, const double &xsi2);

    bounded_matrix<double, 2, 2> currentJacobianMatrix(const double &xsi1, const double &xsi2);

    double jacobianDeterminant(const bounded_matrix<double, 2, 2> &jacobianMatrix);

    std::pair<vector<double>, matrix<double>> elementContributions(const std::string &ep, const std::string &typeAnalyze,
                                                                   const int &step, const int &numberOfStep, const double &deltat,
                                                                   const double &beta, const double &gama,
                                                                   const int &nHammer, const int &nHammerBlendZone);

    void setShapeForce(const bounded_vector<double, 2> &shapeForce);

    void setAnalysisParameters(const double &numberOfDomainIntegrationPoints, const double &deltat, const double &beta, const double &gamma);

    matrix<double> massMatrix(const int &nHammer, const int &nHammerBlendZone);

    bounded_vector<double, 2> calculateGlobalCoordinate(const bounded_vector<double, 2> &qxsi);

    void setIncidenceOnGlobal(const std::vector<bool> &inside, std::vector<Cell *> cells, const std::vector<bounded_vector<double, 2>> xsis,
                              const std::vector<double> &bValue, const std::vector<bounded_vector<double, 2>> &db_dxsiValue);

    //void setDistanceFromFEBoundary(const std::vector<double> &distance);

    //std::vector<double> getDistanceFromFEBoundary();

    //double blendFunction(const double &distance, const double &blendZoneThickness);

    bounded_matrix<double, 2, 2> referenceJacobianMatrixBlendZone(const matrix<double> &dphi_dxsiBlended, std::vector<ControlPoint *> cpsGLobal);

    bounded_matrix<double, 2, 2> currentJacobianMatrixBlendZone(const matrix<double> &dphi_dxsiBlended, std::vector<ControlPoint *> cpsGLobal);

    bounded_matrix<double, 2, 2> inverseMatrix(const bounded_matrix<double, 2, 2> &matrix);

    std::vector<int> getFreedomDegree();

    std::vector<int> getInactiveControlPoints();

    std::vector<bool> ipInsideBlendZone();

    std::vector<bounded_vector<double, 2>> xsiIncidenceCell();

    std::vector<Cell *> incidenceCell();

    void setBlendFunctionValues(const std::vector<double> &b, const std::vector<bounded_vector<double, 2>> &db_dxsi);

    std::vector<double> getBlendValue();

    bounded_vector<double, 4> getCauchyStress(const bounded_vector<double, 2> &qxsi, const std::string &ep);

    bounded_vector<double, 4> getCauchyStressBlendZone(const bounded_vector<double, 2> &qxsiLocal, const std::string &ep, Cell *cell, const bounded_vector<double, 2> &xsiGlobal, const bounded_vector<double, 2> &blendFunctions);

    vector<double> diagonalMass(const int &nHammer, const int &nHammerBlendZone);

    bounded_vector<double, 2> contributionJ_Integral4(const int &quadraturePoints, const int &side, const std::string &ep, const double &rotation, Node *tipNode, const bool &crack);

    bounded_vector<double, 2> contributionJ_IntegralInitial(const int &quadraturePoints, const int &side, const std::string &ep, const double &rotation, const bounded_vector<double, 2> &tipCoordinates);

    bounded_vector<double, 2> contributionJ_Integral5(const int &quadraturePoints, const int &side, const std::string &ep, const double &rotation, Node *tipNode, const bool &crack);

    bounded_vector<double, 2> contributionJ_IntegralFromRice(const int &quadraturePoints, const int &side, const std::string &ep, const double &rotation, Node *tipNode);

    // void setCauchyStressToNode(const std::string &planeState);

private:
    int index_;

    int order_;

    Mesh *mesh_;

    bounded_vector<double, 2> shapeForce_;

    std::vector<Node *> connection_;

    std::vector<Cell *> incidenceCell_;

    std::vector<bounded_vector<double, 2>> xsiIncidenceCell_;

    std::vector<bool> insideBlendZone_;

    std::vector<double> bValue_;

    std::vector<bounded_vector<double, 2>> db_dxsiValue_;
};