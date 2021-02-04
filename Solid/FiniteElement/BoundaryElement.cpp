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
        phi(1, 1) = (27.0 * xsi * xsi + 18.0 * xsi - 1.0) / 16.0;
        phi(2, 1) = (81.0 * xsi * xsi - 18.0 * xsi - 27.0) / 16.0;
        phi(3, 1) = -(81.0 * xsi * xsi + 18.0 * xsi - 27.0) / 16.0;

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
    else if (points == 5)
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
    else if (points == 6)
    {
        pointCoordIso(0, 0) = -0.9324695142;
        pointCoordIso(0, 1) = 0.1713244923;

        pointCoordIso(1, 0) = -0.6612093864;
        pointCoordIso(1, 1) = 0.360761573;

        pointCoordIso(2, 0) = -0.238619186;
        pointCoordIso(2, 1) = 0.4679139345;

        pointCoordIso(3, 0) = 0.238619186;
        pointCoordIso(3, 1) = 0.4679139345;

        pointCoordIso(4, 0) = 0.6612093864;
        pointCoordIso(4, 1) = 0.360761573;

        pointCoordIso(5, 0) = 0.9324695142;
        pointCoordIso(5, 1) = 0.1713244923;
    }
    else if (points == 7)
    {
        pointCoordIso(0, 0) = -0.9491079123;
        pointCoordIso(0, 1) = 0.1294849661;

        pointCoordIso(1, 0) = -0.7415311855;
        pointCoordIso(1, 1) = 0.2797053914;

        pointCoordIso(2, 0) = -0.4058451513;
        pointCoordIso(2, 1) = 0.3818300505;

        pointCoordIso(3, 0) = 0.0;
        pointCoordIso(3, 1) = 0.4179591836;

        pointCoordIso(4, 0) = 0.4058451513;
        pointCoordIso(4, 1) = 0.3818300505;

        pointCoordIso(5, 0) = 0.7415311855;
        pointCoordIso(5, 1) = 0.2797053914;

        pointCoordIso(6, 0) = 0.9491079123;
        pointCoordIso(6, 1) = 0.1294849661;
    }
    else if (points == 8)
    {
        pointCoordIso(0, 0) = -0.9602898564;
        pointCoordIso(0, 1) = 0.1012285362;

        pointCoordIso(1, 0) = -0.7966664774;
        pointCoordIso(1, 1) = 0.2223810344;

        pointCoordIso(2, 0) = -0.5255324099;
        pointCoordIso(2, 1) = 0.3137066458;

        pointCoordIso(3, 0) = -0.1834346424;
        pointCoordIso(3, 1) = 0.3626837833;

        pointCoordIso(4, 0) = 0.1834346424;
        pointCoordIso(4, 1) = 0.3626837833;

        pointCoordIso(5, 0) = 0.5255324099;
        pointCoordIso(5, 1) = 0.3137066458;

        pointCoordIso(6, 0) = 0.7966664774;
        pointCoordIso(6, 1) = 0.2223810344;

        pointCoordIso(7, 0) = 0.9602898564;
        pointCoordIso(7, 1) = 0.1012285362;
    }
    else if (points == 10)
    {
        pointCoordIso(0, 0) = 0.148874338981631;
        pointCoordIso(0, 1) = 0.295524224714753;

        pointCoordIso(1, 0) = 0.433395394129247;
        pointCoordIso(1, 1) = 0.269266719309996;

        pointCoordIso(2, 0) = 0.679409568299024;
        pointCoordIso(2, 1) = 0.219086362515982;

        pointCoordIso(3, 0) = 0.865063366688985;
        pointCoordIso(3, 1) = 0.149451349150581;

        pointCoordIso(4, 0) = 0.973906528517172;
        pointCoordIso(4, 1) = 0.0666713443086883;

        pointCoordIso(5, 0) = -0.148874338981631;
        pointCoordIso(5, 1) = 0.295524224714753;

        pointCoordIso(6, 0) = -0.433395394129247;
        pointCoordIso(6, 1) = 0.269266719309996;

        pointCoordIso(7, 0) = -0.679409568299024;
        pointCoordIso(7, 1) = 0.219086362515982;

        pointCoordIso(8, 0) = -0.865063366688985;
        pointCoordIso(8, 1) = 0.149451349150581;

        pointCoordIso(9, 0) = -0.973906528517172;
        pointCoordIso(9, 1) = 0.066671344308688;
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

        // bounded_vector<double, 2> normal;
        // normal(0)=tangent(1)/jacobian;
        // normal(1)=-tangent(0)/jacobian;

        // double aux = 0.0;
        // for (int m = 0; m < npc; m++)
        // {
        //     aux += functions.first(m);
        // }

        for (int ih = 0; ih < nnode; ih++)
        {
            distribuitedLoad(2 * ih) += value(0) * functions(ih, 0) * weight * jacobian;     // * normal(0);
            distribuitedLoad(2 * ih + 1) += value(1) * functions(ih, 0) * weight * jacobian; // * normal(1);
        }
    }
    return distribuitedLoad;
}

bounded_vector<double, 2> BoundaryElement::getVectorTangente(const double &xsi, const std::string &type)
{

    matrix<double> functions = shapeFunctionsAndDerivates(xsi); //phi, phi', phi''
    bounded_vector<double, 2> tangent;
    tangent(0) = 0.0;
    tangent(1) = 0.0;

    if (type == "current")
    {
        for (int ih = 0; ih < connection_.size(); ih++)
        {
            tangent(0) += functions(ih, 1) * connection_[ih]->getCurrentCoordinate()(0);
            tangent(1) += functions(ih, 1) * connection_[ih]->getCurrentCoordinate()(1);
        }
    }
    else if (type == "initial")
    {
        for (int ih = 0; ih < connection_.size(); ih++)
        {
            tangent(0) += functions(ih, 1) * connection_[ih]->getInitialCoordinate()(0);
            tangent(1) += functions(ih, 1) * connection_[ih]->getInitialCoordinate()(1);
        }
    }
    else
    {
        std::cout << "Something is wrong..." << std::endl;
    }

    return (1.0 / norm_2(tangent)) * tangent;
}
