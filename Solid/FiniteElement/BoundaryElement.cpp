#include "BoundaryElement.h"

BoundaryElement::BoundaryElement() {}

BoundaryElement::BoundaryElement(const int &boundaryIndex, std::vector<Node *> connection)
{
    boundaryIndex_ = boundaryIndex;
    connection_ = connection;
    //mesh_ = mesh;
    //std::string BoundaryElementType = mesh->getBoundaryElementType();
    // if (BoundaryElementType == "T3")
    // {
    //     order_ = 1;
    // }
    // else if (BoundaryElementType == "T6")
    // {
    //     order_ = 2;
    // }
    // else if (BoundaryElementType == "T10")
    // {
    //     order_ = 3;
    // }

    //int n = (order_ + 1) * (order_ + 2) / 2.0;

    //connection_.reserve(n);
}

BoundaryElement::~BoundaryElement() {}

int BoundaryElement::getBoundaryIndex()
{
    return boundaryIndex_;
}

std::vector<Node *> BoundaryElement::getNodes()
{
    return connection_;
}

Node *BoundaryElement::getNode(const int &index)
{
    return connection_[index];
}

matrix<double> BoundaryElement::shapeFunctionsAndDerivates(const double &xsi)
{
    int nconec = connection_.size();
    matrix<double> phi(nconec, 3, 0.0); //phi, phi', phi''

    if (nconec == 2) //linear
    {
        phi(0, 0) = 0.5 * (1.0 - xsi);
        phi(1, 0) = 0.5 * (xsi + 1.0);

        phi(0, 1) = -0.5;
        phi(1, 1) = 0.5;

        phi(0, 2) = 0.0;
        phi(1, 2) = 0.0;
    }
    else if (nconec == 3) //quadratic
    {
        phi(0, 0) = 0.5 * xsi * (xsi - 1.0);
        phi(2, 0) = (1.0 + xsi) * (1.0 - xsi);
        phi(1, 0) = 0.5 * (xsi + 1.0) * xsi;

        phi(0, 1) = xsi - 0.5;
        phi(2, 1) = -2.0 * xsi;
        phi(1, 1) = xsi + 0.5;

        phi(0, 2) = 1.0;
        phi(2, 2) = -2.0;
        phi(1, 2) = 1.0;
    }
    else if (nconec == 4) //cubic
    {
        phi(0, 0) = -(9.0 * xsi * xsi * xsi - 9.0 * xsi * xsi - xsi + 1.0) / 16.0;
        phi(2, 0) = (27.0 * xsi * xsi * xsi - 9.0 * xsi * xsi - 27.0 * xsi + 9.0) / 16.0;
        phi(3, 0) = -(27.0 * xsi * xsi * xsi + 9.0 * xsi * xsi - 27.0 * xsi - 9.0) / 16.0;
        phi(1, 0) = (9.0 * xsi * xsi * xsi + 9.0 * xsi * xsi - xsi - 1.0) / 16.0;

        phi(0, 1) = -(27.0 * xsi * xsi - 18.0 * xsi - 1.0) / 16.0;
        phi(2, 1) = (81.0 * xsi * xsi - 18.0 * xsi - 27.0) / 16.0;
        phi(3, 1) = -(81.0 * xsi * xsi + 18.0 * xsi - 27.0) / 16.0;
        phi(1, 1) = (27.0 * xsi * xsi + 18.0 * xsi - 1.0) / 16.0;

        phi(0, 2) = -(27.0 * xsi - 9.0) / 8.0;
        phi(2, 2) = (81.0 * xsi - 9.0) / 8.0;
        phi(3, 2) = -(81.0 * xsi + 9.0) / 8.0;
        phi(1, 2) = (27.0 * xsi + 9.0) / 8.0;
    }
    return phi;
}

matrix<double> BoundaryElement::boundaryIsoQuadrature(const int &points)
{
    matrix<double> pointCoordIso(points, 2); //xsi1, xsi2, weight

    if (points == 3)
    {
        pointCoordIso(0, 0) = -0.774596669241483;
        pointCoordIso(0, 1) = 0.555555555555556;

        pointCoordIso(1, 0) = 0.774596669241483;
        pointCoordIso(1, 1) = 0.555555555555556;

        pointCoordIso(2, 0) = 0.0;
        pointCoordIso(2, 1) = 0.888888888888889;
    }
    else if (points == 4)
    {
        pointCoordIso(0, 0) = -0.861136311594053;
        pointCoordIso(0, 1) = 0.347854845137454;

        pointCoordIso(1, 0) = 0.861136311594053;
        pointCoordIso(1, 1) = 0.347854845137454;

        pointCoordIso(2, 0) = -0.339981043584856;
        pointCoordIso(2, 1) = 0.652145154862546;

        pointCoordIso(3, 0) = 0.339981043584856;
        pointCoordIso(3, 1) = 0.652145154862546;
    }

    if (points == 5)
    {
        pointCoordIso(0, 0) = -0.906179845938664;
        pointCoordIso(0, 1) = 0.236926885056189;

        pointCoordIso(1, 0) = 0.906179845938664;
        pointCoordIso(1, 1) = 0.236926885056189;

        pointCoordIso(2, 0) = -0.538469310195683;
        pointCoordIso(2, 1) = 0.478628670499366;

        pointCoordIso(3, 0) = 0.538469310195683;
        pointCoordIso(3, 1) = 0.478628670499366;

        pointCoordIso(4, 0) = 0.0;
        pointCoordIso(4, 1) = 0.568888888888889;
    }

    return pointCoordIso;
}

vector<double> BoundaryElement::computeDistribuitedLoads(const bounded_vector<double, 2> &value, const int &quadraturePoints)
{
    int nnode = connection_.size();
    vector<double> distribuitedLoad(2 * nnode, 0.0);
    matrix<double> integrationPoints(quadraturePoints, 2);
    integrationPoints = boundaryIsoQuadrature(quadraturePoints);
    //double thickness = patch_->getThickness();
    for (int iq = 0; iq < quadraturePoints; iq++)
    {
        double xsi = integrationPoints(iq, 0);
        double weight = integrationPoints(iq, 1);

        matrix<double> functions = shapeFunctionsAndDerivates(xsi); //phi, phi', phi''
        bounded_vector<double, 2> tangent;
        tangent(0) = 0.0;
        tangent(1) = 0.0;
        for (int ih = 0; ih < nnode; ih++)
        {
            tangent(0) += functions(ih, 1) * connection_[ih]->getCurrentCoordinate()(0);
            tangent(1) += functions(ih, 1) * connection_[ih]->getCurrentCoordinate()(1);
        }

        double jacobian = sqrt(pow(tangent(0), 2) + pow(tangent(1), 2));

        // double aux = 0.0;
        // for (int m = 0; m < npc; m++)
        // {
        //     aux += functions.first(m);
        // }

        for (int ih = 0; ih < nnode; ih++)
        {
            distribuitedLoad(2 * ih) += value(0) * functions(ih, 0) * weight * jacobian;
            distribuitedLoad(2 * ih + 1) += value(1) * functions(ih, 0) * weight * jacobian;
        }
    }
    return distribuitedLoad;
}
