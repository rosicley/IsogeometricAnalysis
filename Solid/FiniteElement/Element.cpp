#include "Element.h"

Element::Element() {}

Element::Element(const int &index,
                 Mesh *mesh,
                 const std::vector<Node *> &connection)
{
    index_ = index;
    connection_ = connection;
    shapeForce_(0) = 0.0;
    shapeForce_(1) = 0.0;
    mesh_ = mesh;
    order_ = mesh->getOrder();
    //int n = (order_ + 1) * (order_ + 2) / 2.0;

    //connection_.reserve(n);
}

Element::~Element() {}

int Element::getIndex()
{
    return index_;
}

std::vector<Node *> Element::getConnection()
{
    return connection_;
}

Node *Element::getNode(int index)
{
    return connection_[index];
}

Mesh *Element::getMesh()
{
    return mesh_;
}

void Element::setShapeForce(const bounded_vector<double, 2> &shapeForce)
{
    shapeForce_ = shapeForce;
}

vector<double> Element::domainShapeFunction(const double &xsi1, const double &xsi2)
{
    int n = (order_ + 1) * (order_ + 2) / 2.0; //number of nodes per element
    vector<double> phi(n, 0.0);

    if (order_ == 1)
    {
        phi(0) = xsi1;
        phi(1) = xsi2;
        phi(2) = 1.0 - xsi1 - xsi2;
    }
    else if (order_ == 2)
    {
        phi(0) = xsi1 * (2.0 * xsi1 - 1.0);
        phi(1) = xsi2 * (2.0 * xsi2 - 1.0);
        phi(2) = (xsi2 + xsi1 - 1.0) * (2.0 * xsi2 + 2.0 * xsi1 - 1.0);
        phi(3) = 4.0 * xsi1 * xsi2;
        phi(4) = -4.0 * xsi2 * (xsi2 + xsi1 - 1.0);
        phi(5) = -4.0 * xsi1 * (xsi2 + xsi1 - 1.0);
    }
    else if (order_ == 3)
    {
        phi(0) = (xsi1 * (3.0 * xsi1 - 2.0) * (3.0 * xsi1 - 1.0)) / 2.0;
        phi(1) = (xsi2 * (3.0 * xsi2 - 2.0) * (3.0 * xsi2 - 1.0)) / 2.0;
        phi(2) = -((xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0) * (3.0 * xsi2 + 3.0 * xsi1 - 1.0)) / 2.0;
        phi(3) = (9.0 * xsi1 * xsi2 * (3.0 * xsi1 - 1.0)) / 2.0;
        phi(4) = (9.0 * xsi1 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        phi(5) = -(9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 - 1.0)) / 2.0;
        phi(6) = (9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
        phi(7) = (9.0 * xsi1 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
        phi(8) = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0) * (xsi2 + xsi1 - 1.0)) / 2.0;
        phi(9) = -27.0 * xsi1 * xsi2 * (xsi2 + xsi1 - 1.0);
    }

    return phi;
}

matrix<double> Element::domainDerivativeShapeFunction(const double &xsi1, const double &xsi2)
{
    int n = (order_ + 1) * (order_ + 2) / 2.0; //number of nodes per element
    matrix<double> dphi_dxsi(2, n, 0.0);

    if (order_ == 1)
    {
        dphi_dxsi(0, 0) = 1.0;
        dphi_dxsi(0, 1) = 0.0;
        dphi_dxsi(0, 2) = -1.0;

        dphi_dxsi(1, 0) = 0.0;
        dphi_dxsi(1, 1) = 1.0;
        dphi_dxsi(1, 2) = -1.0;
    }
    else if (order_ == 2)
    {
        dphi_dxsi(0, 0) = 4.0 * xsi1 - 1.0;
        dphi_dxsi(0, 1) = 0.0;
        dphi_dxsi(0, 2) = 4.0 * xsi2 + 4.0 * xsi1 - 3.0;
        dphi_dxsi(0, 3) = 4.0 * xsi2;
        dphi_dxsi(0, 4) = -4.0 * xsi2;
        dphi_dxsi(0, 5) = -4.0 * (xsi2 + 2.0 * xsi1 - 1.0);

        dphi_dxsi(1, 0) = 0.0;
        dphi_dxsi(1, 1) = 4.0 * xsi2 - 1.0;
        dphi_dxsi(1, 2) = 4.0 * xsi2 + 4.0 * xsi1 - 3.0;
        dphi_dxsi(1, 3) = 4.0 * xsi1;
        dphi_dxsi(1, 4) = -4.0 * (2.0 * xsi2 + xsi1 - 1.0);
        dphi_dxsi(1, 5) = -4.0 * xsi1;
    }
    else if (order_ == 3)
    {
        dphi_dxsi(0, 0) = (27.0 * xsi1 * xsi1 - 18.0 * xsi1 + 2.0) / 2.0;
        dphi_dxsi(0, 1) = 0.0;
        dphi_dxsi(0, 2) = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
        dphi_dxsi(0, 3) = (9.0 * xsi2 * (6.0 * xsi1 - 1.0)) / 2.0;
        dphi_dxsi(0, 4) = (9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        dphi_dxsi(0, 5) = -(9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        dphi_dxsi(0, 6) = (9.0 * xsi2 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
        dphi_dxsi(0, 7) = (9.0 * (3.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 5.0 * xsi2 + 9.0 * xsi1 * xsi1 - 10.0 * xsi1 + 2.0)) / 2.0;
        dphi_dxsi(0, 8) = -(9.0 * (6.0 * xsi1 * xsi2 - xsi2 + 9.0 * xsi1 * xsi1 - 8.0 * xsi1 + 1.0)) / 2.0;
        dphi_dxsi(0, 9) = -27.0 * xsi2 * (xsi2 + 2.0 * xsi1 - 1.0);

        dphi_dxsi(1, 0) = 0.0;
        dphi_dxsi(1, 1) = (27.0 * xsi2 * xsi2 - 18.0 * xsi2 + 2) / 2.0;
        dphi_dxsi(1, 2) = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
        dphi_dxsi(1, 3) = (9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
        dphi_dxsi(1, 4) = (9.0 * xsi1 * (6.0 * xsi2 - 1.0)) / 2.0;
        dphi_dxsi(1, 5) = -(9.0 * (9.0 * xsi2 * xsi2 + 6.0 * xsi1 * xsi2 - 8.0 * xsi2 - xsi1 + 1.0)) / 2.0;
        dphi_dxsi(1, 6) = (9.0 * (9.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 10.0 * xsi2 + 3.0 * xsi1 * xsi1 - 5.0 * xsi1 + 2.0)) / 2.0;
        dphi_dxsi(1, 7) = (9.0 * xsi1 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
        dphi_dxsi(1, 8) = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
        dphi_dxsi(1, 9) = -27.0 * xsi1 * (2.0 * xsi2 + xsi1 - 1.0);
    }

    return dphi_dxsi;
}

matrix<double> Element::domainSecondDerivativeShapeFunction(const double &xsi1, const double &xsi2)
{
    int n = (order_ + 1) * (order_ + 2) / 2.0; //number of nodes per element
    matrix<double> d2phi_dxsi2(3, n, 0.0);

    if (order_ == 1)
    {
    }
    else if (order_ == 2)
    {
        // d2Phi/dsi1 dxsi1
        d2phi_dxsi2(0, 0) = 4.0;
        d2phi_dxsi2(0, 1) = 0.0;
        d2phi_dxsi2(0, 2) = 4.0;
        d2phi_dxsi2(0, 3) = 0.0;
        d2phi_dxsi2(0, 4) = 0.0;
        d2phi_dxsi2(0, 5) = -8.0;

        // d2Phi/dsi1 dxsi2
        d2phi_dxsi2(1, 0) = 0.0;
        d2phi_dxsi2(1, 1) = 0.0;
        d2phi_dxsi2(1, 2) = 4.0;
        d2phi_dxsi2(1, 3) = 4.0;
        d2phi_dxsi2(1, 4) = -4.0;
        d2phi_dxsi2(1, 5) = -4.0;

        // d2Phi/dsi2 dxsi2
        d2phi_dxsi2(2, 0) = 0.0;
        d2phi_dxsi2(2, 1) = 4.0;
        d2phi_dxsi2(2, 2) = 4.0;
        d2phi_dxsi2(2, 3) = 0.0;
        d2phi_dxsi2(2, 4) = -8.0;
        d2phi_dxsi2(2, 5) = 0.0;
    }
    else if (order_ == 3)
    {
        // d2Phi/dsi1 dxsi1
        d2phi_dxsi2(0, 0) = 0.5 * (54.0 * xsi1 - 18.0);
        d2phi_dxsi2(0, 1) = 0.0;
        d2phi_dxsi2(0, 2) = -0.5 * (54.0 * xsi2 + 54.0 * xsi1 - 36.0);
        d2phi_dxsi2(0, 3) = 0.5 * (54.0 * xsi2);
        d2phi_dxsi2(0, 4) = 0.0;
        d2phi_dxsi2(0, 5) = 0.0;
        d2phi_dxsi2(0, 6) = 0.5 * (54.0 * xsi2);
        d2phi_dxsi2(0, 7) = 0.5 * (108.0 * xsi2 + 162.0 * xsi1 - 90.0);
        d2phi_dxsi2(0, 8) = -0.5 * (54.0 * xsi2 + 162.0 * xsi1 - 72.0);
        d2phi_dxsi2(0, 9) = -54.0 * xsi2;

        // d2Phi/dsi1 dxsi2
        d2phi_dxsi2(1, 0) = 0.0;
        d2phi_dxsi2(1, 1) = 0.0;
        d2phi_dxsi2(1, 2) = -0.5 * (54.0 * xsi2 + 54.0 * xsi1 - 36.0);
        d2phi_dxsi2(1, 3) = 0.5 * (54.0 * xsi1 - 9.0);
        d2phi_dxsi2(1, 4) = 0.5 * (54.0 * xsi2 - 9.0);
        d2phi_dxsi2(1, 5) = -0.5 * (54.0 * xsi2 - 9.0);
        d2phi_dxsi2(1, 6) = 0.5 * (108.0 * xsi2 + 54.0 * xsi1 - 45.0);
        d2phi_dxsi2(1, 7) = 0.5 * (54.0 * xsi2 + 108.0 * xsi1 - 45.0);
        d2phi_dxsi2(1, 8) = -0.5 * (54.0 * xsi1 - 9.0);
        d2phi_dxsi2(1, 9) = -54.0 * xsi2 - 54.0 * xsi1 + 27.0;

        // d2Phi/dsi2 dxsi2
        d2phi_dxsi2(2, 0) = 0.0;
        d2phi_dxsi2(2, 1) = 0.5 * (54.0 * xsi2 - 18.0);
        d2phi_dxsi2(2, 2) = -0.5 * (54.0 * xsi2 + 54.0 * xsi1 - 36.0);
        d2phi_dxsi2(2, 3) = 0.0;
        d2phi_dxsi2(2, 4) = 0.5 * (54.0 * xsi1);
        d2phi_dxsi2(2, 5) = -0.5 * (162.0 * xsi2 + 54.0 * xsi1 - 72.0);
        d2phi_dxsi2(2, 6) = 0.5 * (162.0 * xsi2 + 108.0 * xsi1 - 90.0);
        d2phi_dxsi2(2, 7) = 0.5 * (54.0 * xsi1);
        d2phi_dxsi2(2, 8) = 0.0;
        d2phi_dxsi2(2, 9) = -54.0 * xsi1;
    }

    return d2phi_dxsi2;
}

matrix<double> Element::hammerQuadrature(const int &nH)
{
    matrix<double> pointCoord(nH, 3);

    if (nH == 7)
    {
        pointCoord(0, 0) = 0.333333333333333;
        pointCoord(0, 1) = 0.333333333333333;
        pointCoord(0, 2) = 0.112500000000000;
        pointCoord(1, 0) = 0.797426985353087;
        pointCoord(1, 1) = 0.101286507323456;
        pointCoord(1, 2) = 0.062969590272414;
        pointCoord(2, 0) = 0.101286507323456;
        pointCoord(2, 1) = 0.797426985353087;
        pointCoord(2, 2) = 0.062969590272414;
        pointCoord(3, 0) = 0.101286507323456;
        pointCoord(3, 1) = 0.101286507323456;
        pointCoord(3, 2) = 0.062969590272414;
        pointCoord(4, 0) = 0.470142064105115;
        pointCoord(4, 1) = 0.470142064105115;
        pointCoord(4, 2) = 0.066197076394253;
        pointCoord(5, 0) = 0.059715871789770;
        pointCoord(5, 1) = 0.470142064105115;
        pointCoord(5, 2) = 0.066197076394253;
        pointCoord(6, 0) = 0.470142064105115;
        pointCoord(6, 1) = 0.059715871789770;
        pointCoord(6, 2) = 0.066197076394253;
    }
    else if (nH == 12)
    {
        pointCoord(0, 0) = 0.501426509658179;
        pointCoord(0, 1) = 0.249286745170910;
        pointCoord(0, 2) = 0.05839313786319;
        pointCoord(1, 0) = 0.249286745170910;
        pointCoord(1, 1) = 0.249286745170910;
        pointCoord(1, 2) = 0.05839313786319;
        pointCoord(2, 0) = 0.249286745170910;
        pointCoord(2, 1) = 0.501426509658179;
        pointCoord(2, 2) = 0.05839313786319;
        pointCoord(3, 0) = 0.873821971016996;
        pointCoord(3, 1) = 0.063089014491502;
        pointCoord(3, 2) = 0.025422453185104;
        pointCoord(4, 0) = 0.063089014491502;
        pointCoord(4, 1) = 0.063089014491502;
        pointCoord(4, 2) = 0.025422453185104;
        pointCoord(5, 0) = 0.063089014491502;
        pointCoord(5, 1) = 0.873821971016996;
        pointCoord(5, 2) = 0.025422453185104;
        pointCoord(6, 0) = 0.053145049844816;
        pointCoord(6, 1) = 0.310352451033785;
        pointCoord(6, 2) = 0.041425537809187;
        pointCoord(7, 0) = 0.310352451033785;
        pointCoord(7, 1) = 0.636502499121399;
        pointCoord(7, 2) = 0.041425537809187;
        pointCoord(8, 0) = 0.636502499121399;
        pointCoord(8, 1) = 0.053145049844816;
        pointCoord(8, 2) = 0.041425537809187;
        pointCoord(9, 0) = 0.310352451033785;
        pointCoord(9, 1) = 0.053145049844816;
        pointCoord(9, 2) = 0.041425537809187;
        pointCoord(10, 0) = 0.636502499121399;
        pointCoord(10, 1) = 0.310352451033785;
        pointCoord(10, 2) = 0.041425537809187;
        pointCoord(11, 0) = 0.053145049844816;
        pointCoord(11, 1) = 0.636502499121399;
        pointCoord(11, 2) = 0.041425537809187;
    }
    else if (nH == 28)
    {
        pointCoord(0, 0) = 0.166666666666667;
        pointCoord(1, 0) = 0.050643253661728;
        pointCoord(2, 0) = 0.398713492676543;
        pointCoord(3, 0) = 0.050643253661728;
        pointCoord(4, 0) = 0.235071032052557;
        pointCoord(5, 0) = 0.235071032052557;
        pointCoord(6, 0) = 0.029857935894885;
        pointCoord(7, 0) = 0.333333333333333;
        pointCoord(8, 0) = 0.101286507323457;
        pointCoord(9, 0) = 0.449356746338272;
        pointCoord(10, 0) = 0.449356746338272;
        pointCoord(11, 0) = 0.264928967947442;
        pointCoord(12, 0) = 0.470142064105115;
        pointCoord(13, 0) = 0.264928967947442;
        pointCoord(14, 0) = 0.666666666666667;
        pointCoord(15, 0) = 0.550643253661728;
        pointCoord(16, 0) = 0.898713492676543;
        pointCoord(17, 0) = 0.550643253661728;
        pointCoord(18, 0) = 0.735071032052557;
        pointCoord(19, 0) = 0.735071032052558;
        pointCoord(20, 0) = 0.529857935894885;
        pointCoord(21, 0) = 0.166666666666667;
        pointCoord(22, 0) = 0.050643253661728;
        pointCoord(23, 0) = 0.398713492676543;
        pointCoord(24, 0) = 0.050643253661728;
        pointCoord(25, 0) = 0.235071032052557;
        pointCoord(26, 0) = 0.235071032052557;
        pointCoord(27, 0) = 0.029857935894885;
        pointCoord(0, 1) = 0.166666666666667;
        pointCoord(1, 1) = 0.050643253661729;
        pointCoord(2, 1) = 0.050643253661729;
        pointCoord(3, 1) = 0.398713492676544;
        pointCoord(4, 1) = 0.029857935894885;
        pointCoord(5, 1) = 0.235071032052557;
        pointCoord(6, 1) = 0.235071032052557;
        pointCoord(7, 1) = 0.333333333333333;
        pointCoord(8, 1) = 0.449356746338272;
        pointCoord(9, 1) = 0.101286507323457;
        pointCoord(10, 1) = 0.449356746338272;
        pointCoord(11, 1) = 0.264928967947442;
        pointCoord(12, 1) = 0.264928967947442;
        pointCoord(13, 1) = 0.470142064105115;
        pointCoord(14, 1) = 0.166666666666667;
        pointCoord(15, 1) = 0.050643253661729;
        pointCoord(16, 1) = 0.050643253661729;
        pointCoord(17, 1) = 0.398713492676544;
        pointCoord(18, 1) = 0.029857935894885;
        pointCoord(19, 1) = 0.235071032052557;
        pointCoord(20, 1) = 0.235071032052557;
        pointCoord(21, 1) = 0.666666666666667;
        pointCoord(22, 1) = 0.550643253661729;
        pointCoord(23, 1) = 0.550643253661729;
        pointCoord(24, 1) = 0.898713492676544;
        pointCoord(25, 1) = 0.529857935894885;
        pointCoord(26, 1) = 0.735071032052558;
        pointCoord(27, 1) = 0.735071032052557;
        pointCoord(0, 2) = 0.0281250000000000;
        pointCoord(1, 2) = 0.0157423975681034;
        pointCoord(2, 2) = 0.0157423975681034;
        pointCoord(3, 2) = 0.0157423975681034;
        pointCoord(4, 2) = 0.0165492690985632;
        pointCoord(5, 2) = 0.0165492690985632;
        pointCoord(6, 2) = 0.0165492690985632;
        pointCoord(7, 2) = 0.0281250000000000;
        pointCoord(8, 2) = 0.0157423975681034;
        pointCoord(9, 2) = 0.0157423975681034;
        pointCoord(10, 2) = 0.0157423975681034;
        pointCoord(11, 2) = 0.0165492690985632;
        pointCoord(12, 2) = 0.0165492690985632;
        pointCoord(13, 2) = 0.0165492690985632;
        pointCoord(14, 2) = 0.0281250000000000;
        pointCoord(15, 2) = 0.0157423975681034;
        pointCoord(16, 2) = 0.0157423975681034;
        pointCoord(17, 2) = 0.0157423975681034;
        pointCoord(18, 2) = 0.0165492690985632;
        pointCoord(19, 2) = 0.0165492690985632;
        pointCoord(20, 2) = 0.0165492690985632;
        pointCoord(21, 2) = 0.0281250000000000;
        pointCoord(22, 2) = 0.0157423975681034;
        pointCoord(23, 2) = 0.0157423975681034;
        pointCoord(24, 2) = 0.0157423975681034;
        pointCoord(25, 2) = 0.0165492690985632;
        pointCoord(26, 2) = 0.0165492690985632;
        pointCoord(27, 2) = 0.0165492690985632;
    }
    else if (nH == 112)
    {
        pointCoord(0, 0) = 0.083333333333333;
        pointCoord(1, 0) = 0.025321626830864;
        pointCoord(2, 0) = 0.199356746338272;
        pointCoord(3, 0) = 0.025321626830864;
        pointCoord(4, 0) = 0.117535516026279;
        pointCoord(5, 0) = 0.117535516026279;
        pointCoord(6, 0) = 0.014928967947443;
        pointCoord(7, 0) = 0.166666666666667;
        pointCoord(8, 0) = 0.050643253661728;
        pointCoord(9, 0) = 0.224678373169136;
        pointCoord(10, 0) = 0.224678373169136;
        pointCoord(11, 0) = 0.132464483973721;
        pointCoord(12, 0) = 0.235071032052557;
        pointCoord(13, 0) = 0.132464483973721;
        pointCoord(14, 0) = 0.333333333333333;
        pointCoord(15, 0) = 0.275321626830864;
        pointCoord(16, 0) = 0.449356746338272;
        pointCoord(17, 0) = 0.275321626830864;
        pointCoord(18, 0) = 0.367535516026279;
        pointCoord(19, 0) = 0.367535516026279;
        pointCoord(20, 0) = 0.264928967947442;
        pointCoord(21, 0) = 0.416666666666667;
        pointCoord(22, 0) = 0.300643253661728;
        pointCoord(23, 0) = 0.474678373169136;
        pointCoord(24, 0) = 0.474678373169136;
        pointCoord(25, 0) = 0.382464483973721;
        pointCoord(26, 0) = 0.485071032052558;
        pointCoord(27, 0) = 0.382464483973721;
        pointCoord(28, 0) = 0.583333333333333;
        pointCoord(29, 0) = 0.525321626830864;
        pointCoord(30, 0) = 0.699356746338272;
        pointCoord(31, 0) = 0.525321626830864;
        pointCoord(32, 0) = 0.617535516026279;
        pointCoord(33, 0) = 0.617535516026279;
        pointCoord(34, 0) = 0.514928967947442;
        pointCoord(35, 0) = 0.666666666666667;
        pointCoord(36, 0) = 0.550643253661728;
        pointCoord(37, 0) = 0.724678373169136;
        pointCoord(38, 0) = 0.724678373169136;
        pointCoord(39, 0) = 0.632464483973721;
        pointCoord(40, 0) = 0.735071032052558;
        pointCoord(41, 0) = 0.632464483973721;
        pointCoord(42, 0) = 0.833333333333333;
        pointCoord(43, 0) = 0.775321626830864;
        pointCoord(44, 0) = 0.949356746338272;
        pointCoord(45, 0) = 0.775321626830864;
        pointCoord(46, 0) = 0.867535516026279;
        pointCoord(47, 0) = 0.867535516026279;
        pointCoord(48, 0) = 0.764928967947442;
        pointCoord(49, 0) = 0.083333333333333;
        pointCoord(50, 0) = 0.025321626830864;
        pointCoord(51, 0) = 0.199356746338272;
        pointCoord(52, 0) = 0.025321626830864;
        pointCoord(53, 0) = 0.117535516026279;
        pointCoord(54, 0) = 0.117535516026279;
        pointCoord(55, 0) = 0.014928967947443;
        pointCoord(56, 0) = 0.166666666666667;
        pointCoord(57, 0) = 0.050643253661728;
        pointCoord(58, 0) = 0.224678373169136;
        pointCoord(59, 0) = 0.224678373169136;
        pointCoord(60, 0) = 0.132464483973721;
        pointCoord(61, 0) = 0.235071032052557;
        pointCoord(62, 0) = 0.132464483973721;
        pointCoord(63, 0) = 0.333333333333333;
        pointCoord(64, 0) = 0.275321626830864;
        pointCoord(65, 0) = 0.449356746338272;
        pointCoord(66, 0) = 0.275321626830864;
        pointCoord(67, 0) = 0.367535516026279;
        pointCoord(68, 0) = 0.367535516026279;
        pointCoord(69, 0) = 0.264928967947442;
        pointCoord(70, 0) = 0.416666666666667;
        pointCoord(71, 0) = 0.300643253661728;
        pointCoord(72, 0) = 0.474678373169136;
        pointCoord(73, 0) = 0.474678373169136;
        pointCoord(74, 0) = 0.382464483973721;
        pointCoord(75, 0) = 0.485071032052558;
        pointCoord(76, 0) = 0.382464483973721;
        pointCoord(77, 0) = 0.583333333333333;
        pointCoord(78, 0) = 0.525321626830864;
        pointCoord(79, 0) = 0.699356746338272;
        pointCoord(80, 0) = 0.525321626830864;
        pointCoord(81, 0) = 0.617535516026279;
        pointCoord(82, 0) = 0.617535516026279;
        pointCoord(83, 0) = 0.514928967947442;
        pointCoord(84, 0) = 0.083333333333333;
        pointCoord(85, 0) = 0.025321626830864;
        pointCoord(86, 0) = 0.199356746338272;
        pointCoord(87, 0) = 0.025321626830864;
        pointCoord(88, 0) = 0.117535516026279;
        pointCoord(89, 0) = 0.117535516026279;
        pointCoord(90, 0) = 0.014928967947443;
        pointCoord(91, 0) = 0.166666666666667;
        pointCoord(92, 0) = 0.050643253661728;
        pointCoord(93, 0) = 0.224678373169136;
        pointCoord(94, 0) = 0.224678373169136;
        pointCoord(95, 0) = 0.132464483973721;
        pointCoord(96, 0) = 0.235071032052557;
        pointCoord(97, 0) = 0.132464483973721;
        pointCoord(98, 0) = 0.333333333333333;
        pointCoord(99, 0) = 0.275321626830864;
        pointCoord(100, 0) = 0.449356746338272;
        pointCoord(101, 0) = 0.275321626830864;
        pointCoord(102, 0) = 0.367535516026279;
        pointCoord(103, 0) = 0.367535516026279;
        pointCoord(104, 0) = 0.264928967947442;
        pointCoord(105, 0) = 0.083333333333333;
        pointCoord(106, 0) = 0.025321626830864;
        pointCoord(107, 0) = 0.199356746338272;
        pointCoord(108, 0) = 0.025321626830864;
        pointCoord(109, 0) = 0.117535516026279;
        pointCoord(110, 0) = 0.117535516026279;
        pointCoord(111, 0) = 0.014928967947443;
        pointCoord(0, 1) = 0.166666666666667;
        pointCoord(1, 1) = 0.050643253661728;
        pointCoord(2, 1) = 0.224678373169136;
        pointCoord(3, 1) = 0.224678373169136;
        pointCoord(4, 1) = 0.132464483973721;
        pointCoord(5, 1) = 0.235071032052557;
        pointCoord(6, 1) = 0.132464483973721;
        pointCoord(7, 1) = 0.083333333333333;
        pointCoord(8, 1) = 0.025321626830864;
        pointCoord(9, 1) = 0.025321626830864;
        pointCoord(10, 1) = 0.199356746338272;
        pointCoord(11, 1) = 0.014928967947443;
        pointCoord(12, 1) = 0.117535516026279;
        pointCoord(13, 1) = 0.117535516026279;
        pointCoord(14, 1) = 0.083333333333334;
        pointCoord(15, 1) = 0.025321626830864;
        pointCoord(16, 1) = 0.025321626830864;
        pointCoord(17, 1) = 0.199356746338272;
        pointCoord(18, 1) = 0.014928967947443;
        pointCoord(19, 1) = 0.117535516026279;
        pointCoord(20, 1) = 0.117535516026279;
        pointCoord(21, 1) = 0.166666666666667;
        pointCoord(22, 1) = 0.224678373169136;
        pointCoord(23, 1) = 0.050643253661728;
        pointCoord(24, 1) = 0.224678373169136;
        pointCoord(25, 1) = 0.132464483973721;
        pointCoord(26, 1) = 0.132464483973721;
        pointCoord(27, 1) = 0.235071032052557;
        pointCoord(28, 1) = 0.166666666666667;
        pointCoord(29, 1) = 0.050643253661728;
        pointCoord(30, 1) = 0.224678373169136;
        pointCoord(31, 1) = 0.224678373169136;
        pointCoord(32, 1) = 0.132464483973721;
        pointCoord(33, 1) = 0.235071032052557;
        pointCoord(34, 1) = 0.132464483973721;
        pointCoord(35, 1) = 0.083333333333334;
        pointCoord(36, 1) = 0.025321626830864;
        pointCoord(37, 1) = 0.025321626830864;
        pointCoord(38, 1) = 0.199356746338272;
        pointCoord(39, 1) = 0.014928967947443;
        pointCoord(40, 1) = 0.117535516026279;
        pointCoord(41, 1) = 0.117535516026279;
        pointCoord(42, 1) = 0.083333333333334;
        pointCoord(43, 1) = 0.025321626830864;
        pointCoord(44, 1) = 0.025321626830864;
        pointCoord(45, 1) = 0.199356746338272;
        pointCoord(46, 1) = 0.014928967947443;
        pointCoord(47, 1) = 0.117535516026279;
        pointCoord(48, 1) = 0.117535516026279;
        pointCoord(49, 1) = 0.333333333333333;
        pointCoord(50, 1) = 0.275321626830864;
        pointCoord(51, 1) = 0.275321626830864;
        pointCoord(52, 1) = 0.449356746338272;
        pointCoord(53, 1) = 0.264928967947442;
        pointCoord(54, 1) = 0.367535516026279;
        pointCoord(55, 1) = 0.367535516026279;
        pointCoord(56, 1) = 0.416666666666667;
        pointCoord(57, 1) = 0.474678373169136;
        pointCoord(58, 1) = 0.300643253661728;
        pointCoord(59, 1) = 0.474678373169136;
        pointCoord(60, 1) = 0.382464483973721;
        pointCoord(61, 1) = 0.382464483973721;
        pointCoord(62, 1) = 0.485071032052557;
        pointCoord(63, 1) = 0.416666666666667;
        pointCoord(64, 1) = 0.300643253661728;
        pointCoord(65, 1) = 0.474678373169136;
        pointCoord(66, 1) = 0.474678373169136;
        pointCoord(67, 1) = 0.382464483973721;
        pointCoord(68, 1) = 0.485071032052558;
        pointCoord(69, 1) = 0.382464483973721;
        pointCoord(70, 1) = 0.333333333333333;
        pointCoord(71, 1) = 0.275321626830864;
        pointCoord(72, 1) = 0.275321626830864;
        pointCoord(73, 1) = 0.449356746338272;
        pointCoord(74, 1) = 0.264928967947442;
        pointCoord(75, 1) = 0.367535516026279;
        pointCoord(76, 1) = 0.367535516026279;
        pointCoord(77, 1) = 0.333333333333333;
        pointCoord(78, 1) = 0.275321626830864;
        pointCoord(79, 1) = 0.275321626830864;
        pointCoord(80, 1) = 0.449356746338272;
        pointCoord(81, 1) = 0.264928967947442;
        pointCoord(82, 1) = 0.367535516026279;
        pointCoord(83, 1) = 0.367535516026279;
        pointCoord(84, 1) = 0.666666666666667;
        pointCoord(85, 1) = 0.550643253661728;
        pointCoord(86, 1) = 0.724678373169136;
        pointCoord(87, 1) = 0.724678373169136;
        pointCoord(88, 1) = 0.632464483973721;
        pointCoord(89, 1) = 0.735071032052558;
        pointCoord(90, 1) = 0.632464483973721;
        pointCoord(91, 1) = 0.583333333333333;
        pointCoord(92, 1) = 0.525321626830864;
        pointCoord(93, 1) = 0.525321626830864;
        pointCoord(94, 1) = 0.699356746338272;
        pointCoord(95, 1) = 0.514928967947442;
        pointCoord(96, 1) = 0.617535516026279;
        pointCoord(97, 1) = 0.617535516026279;
        pointCoord(98, 1) = 0.583333333333333;
        pointCoord(99, 1) = 0.525321626830864;
        pointCoord(100, 1) = 0.525321626830864;
        pointCoord(101, 1) = 0.699356746338272;
        pointCoord(102, 1) = 0.514928967947442;
        pointCoord(103, 1) = 0.617535516026279;
        pointCoord(104, 1) = 0.617535516026279;
        pointCoord(105, 1) = 0.833333333333333;
        pointCoord(106, 1) = 0.775321626830864;
        pointCoord(107, 1) = 0.775321626830864;
        pointCoord(108, 1) = 0.949356746338272;
        pointCoord(109, 1) = 0.764928967947442;
        pointCoord(110, 1) = 0.867535516026279;
        pointCoord(111, 1) = 0.867535516026279;
        pointCoord(0, 2) = 0.00703125000000000;
        pointCoord(1, 2) = 0.00393559939202584;
        pointCoord(2, 2) = 0.00393559939202584;
        pointCoord(3, 2) = 0.00393559939202584;
        pointCoord(4, 2) = 0.00413731727464081;
        pointCoord(5, 2) = 0.00413731727464081;
        pointCoord(6, 2) = 0.00413731727464081;
        pointCoord(7, 2) = 0.00703125000000000;
        pointCoord(8, 2) = 0.00393559939202584;
        pointCoord(9, 2) = 0.00393559939202584;
        pointCoord(10, 2) = 0.00393559939202584;
        pointCoord(11, 2) = 0.00413731727464081;
        pointCoord(12, 2) = 0.00413731727464081;
        pointCoord(13, 2) = 0.00413731727464081;
        pointCoord(14, 2) = 0.00703125000000000;
        pointCoord(15, 2) = 0.00393559939202584;
        pointCoord(16, 2) = 0.00393559939202584;
        pointCoord(17, 2) = 0.00393559939202584;
        pointCoord(18, 2) = 0.00413731727464081;
        pointCoord(19, 2) = 0.00413731727464081;
        pointCoord(20, 2) = 0.00413731727464081;
        pointCoord(21, 2) = 0.00703125000000000;
        pointCoord(22, 2) = 0.00393559939202584;
        pointCoord(23, 2) = 0.00393559939202584;
        pointCoord(24, 2) = 0.00393559939202584;
        pointCoord(25, 2) = 0.00413731727464081;
        pointCoord(26, 2) = 0.00413731727464081;
        pointCoord(27, 2) = 0.00413731727464081;
        pointCoord(28, 2) = 0.00703125000000000;
        pointCoord(29, 2) = 0.00393559939202584;
        pointCoord(30, 2) = 0.00393559939202584;
        pointCoord(31, 2) = 0.00393559939202584;
        pointCoord(32, 2) = 0.00413731727464081;
        pointCoord(33, 2) = 0.00413731727464081;
        pointCoord(34, 2) = 0.00413731727464081;
        pointCoord(35, 2) = 0.00703125000000000;
        pointCoord(36, 2) = 0.00393559939202584;
        pointCoord(37, 2) = 0.00393559939202584;
        pointCoord(38, 2) = 0.00393559939202584;
        pointCoord(39, 2) = 0.00413731727464081;
        pointCoord(40, 2) = 0.00413731727464081;
        pointCoord(41, 2) = 0.00413731727464081;
        pointCoord(42, 2) = 0.00703125000000000;
        pointCoord(43, 2) = 0.00393559939202584;
        pointCoord(44, 2) = 0.00393559939202584;
        pointCoord(45, 2) = 0.00393559939202584;
        pointCoord(46, 2) = 0.00413731727464081;
        pointCoord(47, 2) = 0.00413731727464081;
        pointCoord(48, 2) = 0.00413731727464081;
        pointCoord(49, 2) = 0.00703125000000000;
        pointCoord(50, 2) = 0.00393559939202584;
        pointCoord(51, 2) = 0.00393559939202584;
        pointCoord(52, 2) = 0.00393559939202584;
        pointCoord(53, 2) = 0.00413731727464081;
        pointCoord(54, 2) = 0.00413731727464081;
        pointCoord(55, 2) = 0.00413731727464081;
        pointCoord(56, 2) = 0.00703125000000000;
        pointCoord(57, 2) = 0.00393559939202584;
        pointCoord(58, 2) = 0.00393559939202584;
        pointCoord(59, 2) = 0.00393559939202584;
        pointCoord(60, 2) = 0.00413731727464081;
        pointCoord(61, 2) = 0.00413731727464081;
        pointCoord(62, 2) = 0.00413731727464081;
        pointCoord(63, 2) = 0.00703125000000000;
        pointCoord(64, 2) = 0.00393559939202584;
        pointCoord(65, 2) = 0.00393559939202584;
        pointCoord(66, 2) = 0.00393559939202584;
        pointCoord(67, 2) = 0.00413731727464081;
        pointCoord(68, 2) = 0.00413731727464081;
        pointCoord(69, 2) = 0.00413731727464081;
        pointCoord(70, 2) = 0.00703125000000000;
        pointCoord(71, 2) = 0.00393559939202584;
        pointCoord(72, 2) = 0.00393559939202584;
        pointCoord(73, 2) = 0.00393559939202584;
        pointCoord(74, 2) = 0.00413731727464081;
        pointCoord(75, 2) = 0.00413731727464081;
        pointCoord(76, 2) = 0.00413731727464081;
        pointCoord(77, 2) = 0.00703125000000000;
        pointCoord(78, 2) = 0.00393559939202584;
        pointCoord(79, 2) = 0.00393559939202584;
        pointCoord(80, 2) = 0.00393559939202584;
        pointCoord(81, 2) = 0.00413731727464081;
        pointCoord(82, 2) = 0.00413731727464081;
        pointCoord(83, 2) = 0.00413731727464081;
        pointCoord(84, 2) = 0.00703125000000000;
        pointCoord(85, 2) = 0.00393559939202584;
        pointCoord(86, 2) = 0.00393559939202584;
        pointCoord(87, 2) = 0.00393559939202584;
        pointCoord(88, 2) = 0.00413731727464081;
        pointCoord(89, 2) = 0.00413731727464081;
        pointCoord(90, 2) = 0.00413731727464081;
        pointCoord(91, 2) = 0.00703125000000000;
        pointCoord(92, 2) = 0.00393559939202584;
        pointCoord(93, 2) = 0.00393559939202584;
        pointCoord(94, 2) = 0.00393559939202584;
        pointCoord(95, 2) = 0.00413731727464081;
        pointCoord(96, 2) = 0.00413731727464081;
        pointCoord(97, 2) = 0.00413731727464081;
        pointCoord(98, 2) = 0.00703125000000000;
        pointCoord(99, 2) = 0.00393559939202584;
        pointCoord(100, 2) = 0.00393559939202584;
        pointCoord(101, 2) = 0.00393559939202584;
        pointCoord(102, 2) = 0.00413731727464081;
        pointCoord(103, 2) = 0.00413731727464081;
        pointCoord(104, 2) = 0.00413731727464081;
        pointCoord(105, 2) = 0.00703125000000000;
        pointCoord(106, 2) = 0.00393559939202584;
        pointCoord(107, 2) = 0.00393559939202584;
        pointCoord(108, 2) = 0.00393559939202584;
        pointCoord(109, 2) = 0.00413731727464081;
        pointCoord(110, 2) = 0.00413731727464081;
        pointCoord(111, 2) = 0.00413731727464081;
    }
    else if (nH == 14)
    {
        pointCoord(0, 0) = 0.500000000000000;
        pointCoord(1, 0) = 0.151929760985185;
        pointCoord(2, 0) = 0.848070239014815;
        pointCoord(3, 0) = 0.500000000000000;
        pointCoord(4, 0) = 0.500000000000000;
        pointCoord(5, 0) = 0.705213096157672;
        pointCoord(6, 0) = 0.294786903842327;
        pointCoord(7, 0) = 0.166666666666667;
        pointCoord(8, 0) = 0.050643253661728;
        pointCoord(9, 0) = 0.398713492676543;
        pointCoord(10, 0) = 0.050643253661728;
        pointCoord(11, 0) = 0.235071032052557;
        pointCoord(12, 0) = 0.235071032052557;
        pointCoord(13, 0) = 0.029857935894885;
        pointCoord(0, 1) = 0.166666666666667;
        pointCoord(1, 1) = 0.050643253661729;
        pointCoord(2, 1) = 0.050643253661729;
        pointCoord(3, 1) = 0.398713492676544;
        pointCoord(4, 1) = 0.029857935894885;
        pointCoord(5, 1) = 0.235071032052557;
        pointCoord(6, 1) = 0.235071032052557;
        pointCoord(7, 1) = 0.500000000000000;
        pointCoord(8, 1) = 0.151929760985185;
        pointCoord(9, 1) = 0.500000000000000;
        pointCoord(10, 1) = 0.848070239014816;
        pointCoord(11, 1) = 0.294786903842327;
        pointCoord(12, 1) = 0.705213096157672;
        pointCoord(13, 1) = 0.500000000000000;
        pointCoord(0, 2) = 0.0562500000000000;
        pointCoord(1, 2) = 0.0314847951362068;
        pointCoord(2, 2) = 0.0314847951362068;
        pointCoord(3, 2) = 0.0314847951362068;
        pointCoord(4, 2) = 0.0330985381971264;
        pointCoord(5, 2) = 0.0330985381971264;
        pointCoord(6, 2) = 0.0330985381971264;
        pointCoord(7, 2) = 0.0562500000000000;
        pointCoord(8, 2) = 0.0314847951362068;
        pointCoord(9, 2) = 0.0314847951362068;
        pointCoord(10, 2) = 0.0314847951362068;
        pointCoord(11, 2) = 0.0330985381971264;
        pointCoord(12, 2) = 0.0330985381971264;
        pointCoord(13, 2) = 0.0330985381971264;
    }
    else if (nH == 224)
    {
        pointCoord(0, 0) = 0.125000000000000;
        pointCoord(1, 0) = 0.037982440246296;
        pointCoord(2, 0) = 0.212017559753704;
        pointCoord(3, 0) = 0.125000000000000;
        pointCoord(4, 0) = 0.125000000000000;
        pointCoord(5, 0) = 0.176303274039418;
        pointCoord(6, 0) = 0.073696725960582;
        pointCoord(7, 0) = 0.375000000000000;
        pointCoord(8, 0) = 0.287982440246296;
        pointCoord(9, 0) = 0.462017559753704;
        pointCoord(10, 0) = 0.375000000000000;
        pointCoord(11, 0) = 0.375000000000000;
        pointCoord(12, 0) = 0.426303274039418;
        pointCoord(13, 0) = 0.323696725960582;
        pointCoord(14, 0) = 0.625000000000000;
        pointCoord(15, 0) = 0.537982440246296;
        pointCoord(16, 0) = 0.712017559753704;
        pointCoord(17, 0) = 0.625000000000000;
        pointCoord(18, 0) = 0.625000000000000;
        pointCoord(19, 0) = 0.676303274039418;
        pointCoord(20, 0) = 0.573696725960582;
        pointCoord(21, 0) = 0.875000000000000;
        pointCoord(22, 0) = 0.787982440246296;
        pointCoord(23, 0) = 0.962017559753704;
        pointCoord(24, 0) = 0.875000000000000;
        pointCoord(25, 0) = 0.875000000000000;
        pointCoord(26, 0) = 0.926303274039418;
        pointCoord(27, 0) = 0.823696725960582;
        pointCoord(28, 0) = 0.041666666666667;
        pointCoord(29, 0) = 0.012660813415432;
        pointCoord(30, 0) = 0.099678373169136;
        pointCoord(31, 0) = 0.012660813415432;
        pointCoord(32, 0) = 0.058767758013139;
        pointCoord(33, 0) = 0.058767758013139;
        pointCoord(34, 0) = 0.007464483973721;
        pointCoord(35, 0) = 0.208333333333333;
        pointCoord(36, 0) = 0.237339186584568;
        pointCoord(37, 0) = 0.237339186584568;
        pointCoord(38, 0) = 0.150321626830864;
        pointCoord(39, 0) = 0.242535516026279;
        pointCoord(40, 0) = 0.191232241986861;
        pointCoord(41, 0) = 0.191232241986861;
        pointCoord(42, 0) = 0.291666666666667;
        pointCoord(43, 0) = 0.262660813415432;
        pointCoord(44, 0) = 0.349678373169136;
        pointCoord(45, 0) = 0.262660813415432;
        pointCoord(46, 0) = 0.308767758013139;
        pointCoord(47, 0) = 0.308767758013139;
        pointCoord(48, 0) = 0.257464483973721;
        pointCoord(49, 0) = 0.458333333333333;
        pointCoord(50, 0) = 0.487339186584568;
        pointCoord(51, 0) = 0.487339186584568;
        pointCoord(52, 0) = 0.400321626830864;
        pointCoord(53, 0) = 0.492535516026279;
        pointCoord(54, 0) = 0.441232241986861;
        pointCoord(55, 0) = 0.441232241986861;
        pointCoord(56, 0) = 0.541666666666667;
        pointCoord(57, 0) = 0.512660813415432;
        pointCoord(58, 0) = 0.599678373169136;
        pointCoord(59, 0) = 0.512660813415432;
        pointCoord(60, 0) = 0.558767758013139;
        pointCoord(61, 0) = 0.558767758013139;
        pointCoord(62, 0) = 0.507464483973721;
        pointCoord(63, 0) = 0.708333333333333;
        pointCoord(64, 0) = 0.737339186584568;
        pointCoord(65, 0) = 0.737339186584568;
        pointCoord(66, 0) = 0.650321626830864;
        pointCoord(67, 0) = 0.742535516026279;
        pointCoord(68, 0) = 0.691232241986861;
        pointCoord(69, 0) = 0.691232241986861;
        pointCoord(70, 0) = 0.791666666666667;
        pointCoord(71, 0) = 0.762660813415432;
        pointCoord(72, 0) = 0.849678373169136;
        pointCoord(73, 0) = 0.762660813415432;
        pointCoord(74, 0) = 0.808767758013139;
        pointCoord(75, 0) = 0.808767758013139;
        pointCoord(76, 0) = 0.757464483973721;
        pointCoord(77, 0) = 0.125000000000000;
        pointCoord(78, 0) = 0.037982440246296;
        pointCoord(79, 0) = 0.125000000000000;
        pointCoord(80, 0) = 0.212017559753704;
        pointCoord(81, 0) = 0.073696725960582;
        pointCoord(82, 0) = 0.176303274039418;
        pointCoord(83, 0) = 0.125000000000000;
        pointCoord(84, 0) = 0.375000000000000;
        pointCoord(85, 0) = 0.287982440246296;
        pointCoord(86, 0) = 0.375000000000000;
        pointCoord(87, 0) = 0.462017559753704;
        pointCoord(88, 0) = 0.323696725960582;
        pointCoord(89, 0) = 0.426303274039418;
        pointCoord(90, 0) = 0.375000000000000;
        pointCoord(91, 0) = 0.541666666666667;
        pointCoord(92, 0) = 0.338625693908024;
        pointCoord(93, 0) = 0.599678373169136;
        pointCoord(94, 0) = 0.686695932922840;
        pointCoord(95, 0) = 0.456161209934303;
        pointCoord(96, 0) = 0.661374306091976;
        pointCoord(97, 0) = 0.507464483973721;
        pointCoord(98, 0) = 0.125000000000000;
        pointCoord(99, 0) = 0.037982440246296;
        pointCoord(100, 0) = 0.212017559753704;
        pointCoord(101, 0) = 0.125000000000000;
        pointCoord(102, 0) = 0.125000000000000;
        pointCoord(103, 0) = 0.176303274039418;
        pointCoord(104, 0) = 0.073696725960582;
        pointCoord(105, 0) = 0.375000000000000;
        pointCoord(106, 0) = 0.287982440246296;
        pointCoord(107, 0) = 0.462017559753704;
        pointCoord(108, 0) = 0.375000000000000;
        pointCoord(109, 0) = 0.375000000000000;
        pointCoord(110, 0) = 0.426303274039418;
        pointCoord(111, 0) = 0.323696725960582;
        pointCoord(112, 0) = 0.625000000000000;
        pointCoord(113, 0) = 0.537982440246296;
        pointCoord(114, 0) = 0.712017559753704;
        pointCoord(115, 0) = 0.625000000000000;
        pointCoord(116, 0) = 0.625000000000000;
        pointCoord(117, 0) = 0.676303274039418;
        pointCoord(118, 0) = 0.573696725960582;
        pointCoord(119, 0) = 0.041666666666667;
        pointCoord(120, 0) = 0.012660813415432;
        pointCoord(121, 0) = 0.099678373169136;
        pointCoord(122, 0) = 0.012660813415432;
        pointCoord(123, 0) = 0.058767758013139;
        pointCoord(124, 0) = 0.058767758013139;
        pointCoord(125, 0) = 0.007464483973721;
        pointCoord(126, 0) = 0.208333333333333;
        pointCoord(127, 0) = 0.150321626830864;
        pointCoord(128, 0) = 0.237339186584568;
        pointCoord(129, 0) = 0.237339186584568;
        pointCoord(130, 0) = 0.191232241986861;
        pointCoord(131, 0) = 0.242535516026279;
        pointCoord(132, 0) = 0.191232241986861;
        pointCoord(133, 0) = 0.291666666666667;
        pointCoord(134, 0) = 0.262660813415432;
        pointCoord(135, 0) = 0.349678373169136;
        pointCoord(136, 0) = 0.262660813415432;
        pointCoord(137, 0) = 0.308767758013139;
        pointCoord(138, 0) = 0.308767758013139;
        pointCoord(139, 0) = 0.257464483973721;
        pointCoord(140, 0) = 0.458333333333333;
        pointCoord(141, 0) = 0.400321626830864;
        pointCoord(142, 0) = 0.487339186584568;
        pointCoord(143, 0) = 0.487339186584568;
        pointCoord(144, 0) = 0.441232241986861;
        pointCoord(145, 0) = 0.492535516026279;
        pointCoord(146, 0) = 0.441232241986861;
        pointCoord(147, 0) = 0.541666666666667;
        pointCoord(148, 0) = 0.512660813415432;
        pointCoord(149, 0) = 0.599678373169136;
        pointCoord(150, 0) = 0.512660813415432;
        pointCoord(151, 0) = 0.558767758013139;
        pointCoord(152, 0) = 0.558767758013139;
        pointCoord(153, 0) = 0.507464483973721;
        pointCoord(154, 0) = 0.125000000000000;
        pointCoord(155, 0) = 0.037982440246296;
        pointCoord(156, 0) = 0.125000000000000;
        pointCoord(157, 0) = 0.212017559753704;
        pointCoord(158, 0) = 0.073696725960582;
        pointCoord(159, 0) = 0.176303274039418;
        pointCoord(160, 0) = 0.125000000000000;
        pointCoord(161, 0) = 0.375000000000000;
        pointCoord(162, 0) = 0.287982440246296;
        pointCoord(163, 0) = 0.375000000000000;
        pointCoord(164, 0) = 0.462017559753704;
        pointCoord(165, 0) = 0.323696725960582;
        pointCoord(166, 0) = 0.426303274039418;
        pointCoord(167, 0) = 0.375000000000000;
        pointCoord(168, 0) = 0.125000000000000;
        pointCoord(169, 0) = 0.037982440246296;
        pointCoord(170, 0) = 0.212017559753704;
        pointCoord(171, 0) = 0.125000000000000;
        pointCoord(172, 0) = 0.125000000000000;
        pointCoord(173, 0) = 0.176303274039418;
        pointCoord(174, 0) = 0.073696725960582;
        pointCoord(175, 0) = 0.375000000000000;
        pointCoord(176, 0) = 0.287982440246296;
        pointCoord(177, 0) = 0.462017559753704;
        pointCoord(178, 0) = 0.375000000000000;
        pointCoord(179, 0) = 0.375000000000000;
        pointCoord(180, 0) = 0.426303274039418;
        pointCoord(181, 0) = 0.323696725960582;
        pointCoord(182, 0) = 0.041666666666667;
        pointCoord(183, 0) = 0.012660813415432;
        pointCoord(184, 0) = 0.099678373169136;
        pointCoord(185, 0) = 0.012660813415432;
        pointCoord(186, 0) = 0.058767758013139;
        pointCoord(187, 0) = 0.058767758013139;
        pointCoord(188, 0) = 0.007464483973721;
        pointCoord(189, 0) = 0.208333333333333;
        pointCoord(190, 0) = 0.150321626830864;
        pointCoord(191, 0) = 0.237339186584568;
        pointCoord(192, 0) = 0.237339186584568;
        pointCoord(193, 0) = 0.191232241986861;
        pointCoord(194, 0) = 0.242535516026279;
        pointCoord(195, 0) = 0.191232241986861;
        pointCoord(196, 0) = 0.291666666666667;
        pointCoord(197, 0) = 0.262660813415432;
        pointCoord(198, 0) = 0.349678373169136;
        pointCoord(199, 0) = 0.262660813415432;
        pointCoord(200, 0) = 0.308767758013139;
        pointCoord(201, 0) = 0.308767758013139;
        pointCoord(202, 0) = 0.257464483973721;
        pointCoord(203, 0) = 0.125000000000000;
        pointCoord(204, 0) = 0.037982440246296;
        pointCoord(205, 0) = 0.125000000000000;
        pointCoord(206, 0) = 0.212017559753704;
        pointCoord(207, 0) = 0.073696725960582;
        pointCoord(208, 0) = 0.176303274039418;
        pointCoord(209, 0) = 0.125000000000000;
        pointCoord(210, 0) = 0.125000000000000;
        pointCoord(211, 0) = 0.037982440246296;
        pointCoord(212, 0) = 0.212017559753704;
        pointCoord(213, 0) = 0.125000000000000;
        pointCoord(214, 0) = 0.125000000000000;
        pointCoord(215, 0) = 0.176303274039418;
        pointCoord(216, 0) = 0.073696725960582;
        pointCoord(217, 0) = 0.041666666666667;
        pointCoord(218, 0) = 0.012660813415432;
        pointCoord(219, 0) = 0.099678373169136;
        pointCoord(220, 0) = 0.012660813415432;
        pointCoord(221, 0) = 0.058767758013139;
        pointCoord(222, 0) = 0.058767758013139;
        pointCoord(223, 0) = 0.007464483973721;
        pointCoord(0, 1) = 0.041666666666667;
        pointCoord(1, 1) = 0.012660813415432;
        pointCoord(2, 1) = 0.012660813415432;
        pointCoord(3, 1) = 0.099678373169136;
        pointCoord(4, 1) = 0.007464483973721;
        pointCoord(5, 1) = 0.058767758013139;
        pointCoord(6, 1) = 0.058767758013139;
        pointCoord(7, 1) = 0.041666666666667;
        pointCoord(8, 1) = 0.012660813415432;
        pointCoord(9, 1) = 0.012660813415432;
        pointCoord(10, 1) = 0.099678373169136;
        pointCoord(11, 1) = 0.007464483973721;
        pointCoord(12, 1) = 0.058767758013139;
        pointCoord(13, 1) = 0.058767758013139;
        pointCoord(14, 1) = 0.041666666666667;
        pointCoord(15, 1) = 0.012660813415432;
        pointCoord(16, 1) = 0.012660813415432;
        pointCoord(17, 1) = 0.099678373169136;
        pointCoord(18, 1) = 0.007464483973721;
        pointCoord(19, 1) = 0.058767758013139;
        pointCoord(20, 1) = 0.058767758013139;
        pointCoord(21, 1) = 0.041666666666667;
        pointCoord(22, 1) = 0.012660813415432;
        pointCoord(23, 1) = 0.012660813415432;
        pointCoord(24, 1) = 0.099678373169136;
        pointCoord(25, 1) = 0.007464483973721;
        pointCoord(26, 1) = 0.058767758013139;
        pointCoord(27, 1) = 0.058767758013139;
        pointCoord(28, 1) = 0.125000000000000;
        pointCoord(29, 1) = 0.037982440246296;
        pointCoord(30, 1) = 0.125000000000000;
        pointCoord(31, 1) = 0.212017559753704;
        pointCoord(32, 1) = 0.073696725960582;
        pointCoord(33, 1) = 0.176303274039418;
        pointCoord(34, 1) = 0.125000000000000;
        pointCoord(35, 1) = 0.125000000000000;
        pointCoord(36, 1) = 0.037982440246296;
        pointCoord(37, 1) = 0.212017559753704;
        pointCoord(38, 1) = 0.125000000000000;
        pointCoord(39, 1) = 0.125000000000000;
        pointCoord(40, 1) = 0.176303274039418;
        pointCoord(41, 1) = 0.073696725960582;
        pointCoord(42, 1) = 0.125000000000000;
        pointCoord(43, 1) = 0.037982440246296;
        pointCoord(44, 1) = 0.125000000000000;
        pointCoord(45, 1) = 0.212017559753704;
        pointCoord(46, 1) = 0.073696725960582;
        pointCoord(47, 1) = 0.176303274039418;
        pointCoord(48, 1) = 0.125000000000000;
        pointCoord(49, 1) = 0.125000000000000;
        pointCoord(50, 1) = 0.037982440246296;
        pointCoord(51, 1) = 0.212017559753704;
        pointCoord(52, 1) = 0.125000000000000;
        pointCoord(53, 1) = 0.125000000000000;
        pointCoord(54, 1) = 0.176303274039418;
        pointCoord(55, 1) = 0.073696725960582;
        pointCoord(56, 1) = 0.125000000000000;
        pointCoord(57, 1) = 0.037982440246296;
        pointCoord(58, 1) = 0.125000000000000;
        pointCoord(59, 1) = 0.212017559753704;
        pointCoord(60, 1) = 0.073696725960582;
        pointCoord(61, 1) = 0.176303274039418;
        pointCoord(62, 1) = 0.125000000000000;
        pointCoord(63, 1) = 0.125000000000000;
        pointCoord(64, 1) = 0.037982440246296;
        pointCoord(65, 1) = 0.212017559753704;
        pointCoord(66, 1) = 0.125000000000000;
        pointCoord(67, 1) = 0.125000000000000;
        pointCoord(68, 1) = 0.176303274039418;
        pointCoord(69, 1) = 0.073696725960582;
        pointCoord(70, 1) = 0.125000000000000;
        pointCoord(71, 1) = 0.037982440246296;
        pointCoord(72, 1) = 0.125000000000000;
        pointCoord(73, 1) = 0.212017559753704;
        pointCoord(74, 1) = 0.073696725960582;
        pointCoord(75, 1) = 0.176303274039418;
        pointCoord(76, 1) = 0.125000000000000;
        pointCoord(77, 1) = 0.208333333333333;
        pointCoord(78, 1) = 0.237339186584568;
        pointCoord(79, 1) = 0.150321626830864;
        pointCoord(80, 1) = 0.237339186584568;
        pointCoord(81, 1) = 0.191232241986861;
        pointCoord(82, 1) = 0.191232241986861;
        pointCoord(83, 1) = 0.242535516026279;
        pointCoord(84, 1) = 0.208333333333333;
        pointCoord(85, 1) = 0.237339186584568;
        pointCoord(86, 1) = 0.150321626830864;
        pointCoord(87, 1) = 0.237339186584568;
        pointCoord(88, 1) = 0.191232241986861;
        pointCoord(89, 1) = 0.191232241986861;
        pointCoord(90, 1) = 0.242535516026279;
        pointCoord(91, 1) = 0.291666666666667;
        pointCoord(92, 1) = 0.436695932922840;
        pointCoord(93, 1) = 0.175643253661728;
        pointCoord(94, 1) = 0.262660813415432;
        pointCoord(95, 1) = 0.308767758013139;
        pointCoord(96, 1) = 0.206161209934303;
        pointCoord(97, 1) = 0.360071032052558;
        pointCoord(98, 1) = 0.291666666666667;
        pointCoord(99, 1) = 0.262660813415432;
        pointCoord(100, 1) = 0.262660813415432;
        pointCoord(101, 1) = 0.349678373169136;
        pointCoord(102, 1) = 0.257464483973721;
        pointCoord(103, 1) = 0.308767758013139;
        pointCoord(104, 1) = 0.308767758013139;
        pointCoord(105, 1) = 0.291666666666667;
        pointCoord(106, 1) = 0.262660813415432;
        pointCoord(107, 1) = 0.262660813415432;
        pointCoord(108, 1) = 0.349678373169136;
        pointCoord(109, 1) = 0.257464483973721;
        pointCoord(110, 1) = 0.308767758013139;
        pointCoord(111, 1) = 0.308767758013139;
        pointCoord(112, 1) = 0.291666666666667;
        pointCoord(113, 1) = 0.262660813415432;
        pointCoord(114, 1) = 0.262660813415432;
        pointCoord(115, 1) = 0.349678373169136;
        pointCoord(116, 1) = 0.257464483973721;
        pointCoord(117, 1) = 0.308767758013139;
        pointCoord(118, 1) = 0.308767758013139;
        pointCoord(119, 1) = 0.375000000000000;
        pointCoord(120, 1) = 0.287982440246296;
        pointCoord(121, 1) = 0.375000000000000;
        pointCoord(122, 1) = 0.462017559753704;
        pointCoord(123, 1) = 0.323696725960582;
        pointCoord(124, 1) = 0.426303274039418;
        pointCoord(125, 1) = 0.375000000000000;
        pointCoord(126, 1) = 0.375000000000000;
        pointCoord(127, 1) = 0.375000000000000;
        pointCoord(128, 1) = 0.287982440246296;
        pointCoord(129, 1) = 0.462017559753704;
        pointCoord(130, 1) = 0.323696725960582;
        pointCoord(131, 1) = 0.375000000000000;
        pointCoord(132, 1) = 0.426303274039418;
        pointCoord(133, 1) = 0.375000000000000;
        pointCoord(134, 1) = 0.287982440246296;
        pointCoord(135, 1) = 0.375000000000000;
        pointCoord(136, 1) = 0.462017559753704;
        pointCoord(137, 1) = 0.323696725960582;
        pointCoord(138, 1) = 0.426303274039418;
        pointCoord(139, 1) = 0.375000000000000;
        pointCoord(140, 1) = 0.375000000000000;
        pointCoord(141, 1) = 0.375000000000000;
        pointCoord(142, 1) = 0.287982440246296;
        pointCoord(143, 1) = 0.462017559753704;
        pointCoord(144, 1) = 0.323696725960582;
        pointCoord(145, 1) = 0.375000000000000;
        pointCoord(146, 1) = 0.426303274039418;
        pointCoord(147, 1) = 0.375000000000000;
        pointCoord(148, 1) = 0.287982440246296;
        pointCoord(149, 1) = 0.375000000000000;
        pointCoord(150, 1) = 0.462017559753704;
        pointCoord(151, 1) = 0.323696725960582;
        pointCoord(152, 1) = 0.426303274039418;
        pointCoord(153, 1) = 0.375000000000000;
        pointCoord(154, 1) = 0.458333333333333;
        pointCoord(155, 1) = 0.487339186584568;
        pointCoord(156, 1) = 0.400321626830864;
        pointCoord(157, 1) = 0.487339186584568;
        pointCoord(158, 1) = 0.441232241986861;
        pointCoord(159, 1) = 0.441232241986861;
        pointCoord(160, 1) = 0.492535516026279;
        pointCoord(161, 1) = 0.458333333333333;
        pointCoord(162, 1) = 0.487339186584568;
        pointCoord(163, 1) = 0.400321626830864;
        pointCoord(164, 1) = 0.487339186584568;
        pointCoord(165, 1) = 0.441232241986861;
        pointCoord(166, 1) = 0.441232241986861;
        pointCoord(167, 1) = 0.492535516026279;
        pointCoord(168, 1) = 0.541666666666667;
        pointCoord(169, 1) = 0.512660813415432;
        pointCoord(170, 1) = 0.512660813415432;
        pointCoord(171, 1) = 0.599678373169136;
        pointCoord(172, 1) = 0.507464483973721;
        pointCoord(173, 1) = 0.558767758013139;
        pointCoord(174, 1) = 0.558767758013139;
        pointCoord(175, 1) = 0.541666666666667;
        pointCoord(176, 1) = 0.512660813415432;
        pointCoord(177, 1) = 0.512660813415432;
        pointCoord(178, 1) = 0.599678373169136;
        pointCoord(179, 1) = 0.507464483973721;
        pointCoord(180, 1) = 0.558767758013139;
        pointCoord(181, 1) = 0.558767758013139;
        pointCoord(182, 1) = 0.625000000000000;
        pointCoord(183, 1) = 0.537982440246296;
        pointCoord(184, 1) = 0.625000000000000;
        pointCoord(185, 1) = 0.712017559753704;
        pointCoord(186, 1) = 0.573696725960582;
        pointCoord(187, 1) = 0.676303274039418;
        pointCoord(188, 1) = 0.625000000000000;
        pointCoord(189, 1) = 0.625000000000000;
        pointCoord(190, 1) = 0.625000000000000;
        pointCoord(191, 1) = 0.537982440246296;
        pointCoord(192, 1) = 0.712017559753704;
        pointCoord(193, 1) = 0.573696725960582;
        pointCoord(194, 1) = 0.625000000000000;
        pointCoord(195, 1) = 0.676303274039418;
        pointCoord(196, 1) = 0.625000000000000;
        pointCoord(197, 1) = 0.537982440246296;
        pointCoord(198, 1) = 0.625000000000000;
        pointCoord(199, 1) = 0.712017559753704;
        pointCoord(200, 1) = 0.573696725960582;
        pointCoord(201, 1) = 0.676303274039418;
        pointCoord(202, 1) = 0.625000000000000;
        pointCoord(203, 1) = 0.708333333333333;
        pointCoord(204, 1) = 0.737339186584568;
        pointCoord(205, 1) = 0.650321626830864;
        pointCoord(206, 1) = 0.737339186584568;
        pointCoord(207, 1) = 0.691232241986861;
        pointCoord(208, 1) = 0.691232241986861;
        pointCoord(209, 1) = 0.742535516026279;
        pointCoord(210, 1) = 0.791666666666667;
        pointCoord(211, 1) = 0.762660813415432;
        pointCoord(212, 1) = 0.762660813415432;
        pointCoord(213, 1) = 0.849678373169136;
        pointCoord(214, 1) = 0.757464483973721;
        pointCoord(215, 1) = 0.808767758013139;
        pointCoord(216, 1) = 0.808767758013139;
        pointCoord(217, 1) = 0.875000000000000;
        pointCoord(218, 1) = 0.787982440246296;
        pointCoord(219, 1) = 0.875000000000000;
        pointCoord(220, 1) = 0.962017559753704;
        pointCoord(221, 1) = 0.823696725960582;
        pointCoord(222, 1) = 0.926303274039418;
        pointCoord(223, 1) = 0.875000000000000;
        pointCoord(0, 2) = 0.00351562500000000;
        pointCoord(1, 2) = 0.00196779969601292;
        pointCoord(2, 2) = 0.00196779969601292;
        pointCoord(3, 2) = 0.00196779969601292;
        pointCoord(4, 2) = 0.00206865863732041;
        pointCoord(5, 2) = 0.00206865863732041;
        pointCoord(6, 2) = 0.00206865863732041;
        pointCoord(7, 2) = 0.00351562500000000;
        pointCoord(8, 2) = 0.00196779969601292;
        pointCoord(9, 2) = 0.00196779969601292;
        pointCoord(10, 2) = 0.00196779969601292;
        pointCoord(11, 2) = 0.00206865863732041;
        pointCoord(12, 2) = 0.00206865863732041;
        pointCoord(13, 2) = 0.00206865863732041;
        pointCoord(14, 2) = 0.00351562500000000;
        pointCoord(15, 2) = 0.00196779969601292;
        pointCoord(16, 2) = 0.00196779969601292;
        pointCoord(17, 2) = 0.00196779969601292;
        pointCoord(18, 2) = 0.00206865863732041;
        pointCoord(19, 2) = 0.00206865863732041;
        pointCoord(20, 2) = 0.00206865863732041;
        pointCoord(21, 2) = 0.00351562500000000;
        pointCoord(22, 2) = 0.00196779969601292;
        pointCoord(23, 2) = 0.00196779969601292;
        pointCoord(24, 2) = 0.00196779969601292;
        pointCoord(25, 2) = 0.00206865863732041;
        pointCoord(26, 2) = 0.00206865863732041;
        pointCoord(27, 2) = 0.00206865863732041;
        pointCoord(28, 2) = 0.00351562500000000;
        pointCoord(29, 2) = 0.00196779969601292;
        pointCoord(30, 2) = 0.00196779969601292;
        pointCoord(31, 2) = 0.00196779969601292;
        pointCoord(32, 2) = 0.00206865863732041;
        pointCoord(33, 2) = 0.00206865863732041;
        pointCoord(34, 2) = 0.00206865863732041;
        pointCoord(35, 2) = 0.00351562500000000;
        pointCoord(36, 2) = 0.00196779969601292;
        pointCoord(37, 2) = 0.00196779969601292;
        pointCoord(38, 2) = 0.00196779969601292;
        pointCoord(39, 2) = 0.00206865863732041;
        pointCoord(40, 2) = 0.00206865863732041;
        pointCoord(41, 2) = 0.00206865863732041;
        pointCoord(42, 2) = 0.00351562500000000;
        pointCoord(43, 2) = 0.00196779969601292;
        pointCoord(44, 2) = 0.00196779969601292;
        pointCoord(45, 2) = 0.00196779969601292;
        pointCoord(46, 2) = 0.00206865863732041;
        pointCoord(47, 2) = 0.00206865863732041;
        pointCoord(48, 2) = 0.00206865863732041;
        pointCoord(49, 2) = 0.00351562500000000;
        pointCoord(50, 2) = 0.00196779969601292;
        pointCoord(51, 2) = 0.00196779969601292;
        pointCoord(52, 2) = 0.00196779969601292;
        pointCoord(53, 2) = 0.00206865863732041;
        pointCoord(54, 2) = 0.00206865863732041;
        pointCoord(55, 2) = 0.00206865863732041;
        pointCoord(56, 2) = 0.00351562500000000;
        pointCoord(57, 2) = 0.00196779969601292;
        pointCoord(58, 2) = 0.00196779969601292;
        pointCoord(59, 2) = 0.00196779969601292;
        pointCoord(60, 2) = 0.00206865863732041;
        pointCoord(61, 2) = 0.00206865863732041;
        pointCoord(62, 2) = 0.00206865863732041;
        pointCoord(63, 2) = 0.00351562500000000;
        pointCoord(64, 2) = 0.00196779969601292;
        pointCoord(65, 2) = 0.00196779969601292;
        pointCoord(66, 2) = 0.00196779969601292;
        pointCoord(67, 2) = 0.00206865863732041;
        pointCoord(68, 2) = 0.00206865863732041;
        pointCoord(69, 2) = 0.00206865863732041;
        pointCoord(70, 2) = 0.00351562500000000;
        pointCoord(71, 2) = 0.00196779969601292;
        pointCoord(72, 2) = 0.00196779969601292;
        pointCoord(73, 2) = 0.00196779969601292;
        pointCoord(74, 2) = 0.00206865863732041;
        pointCoord(75, 2) = 0.00206865863732041;
        pointCoord(76, 2) = 0.00206865863732041;
        pointCoord(77, 2) = 0.00351562500000000;
        pointCoord(78, 2) = 0.00196779969601292;
        pointCoord(79, 2) = 0.00196779969601292;
        pointCoord(80, 2) = 0.00196779969601292;
        pointCoord(81, 2) = 0.00206865863732041;
        pointCoord(82, 2) = 0.00206865863732041;
        pointCoord(83, 2) = 0.00206865863732041;
        pointCoord(84, 2) = 0.00351562500000000;
        pointCoord(85, 2) = 0.00196779969601292;
        pointCoord(86, 2) = 0.00196779969601292;
        pointCoord(87, 2) = 0.00196779969601292;
        pointCoord(88, 2) = 0.00206865863732041;
        pointCoord(89, 2) = 0.00206865863732041;
        pointCoord(90, 2) = 0.00206865863732041;
        pointCoord(91, 2) = 0.00351562500000000;
        pointCoord(92, 2) = 0.00196779969601292;
        pointCoord(93, 2) = 0.00196779969601292;
        pointCoord(94, 2) = 0.00196779969601292;
        pointCoord(95, 2) = 0.00206865863732041;
        pointCoord(96, 2) = 0.00206865863732041;
        pointCoord(97, 2) = 0.00206865863732041;
        pointCoord(98, 2) = 0.00351562500000000;
        pointCoord(99, 2) = 0.00196779969601292;
        pointCoord(100, 2) = 0.00196779969601292;
        pointCoord(101, 2) = 0.00196779969601292;
        pointCoord(102, 2) = 0.00206865863732041;
        pointCoord(103, 2) = 0.00206865863732041;
        pointCoord(104, 2) = 0.00206865863732041;
        pointCoord(105, 2) = 0.00351562500000000;
        pointCoord(106, 2) = 0.00196779969601292;
        pointCoord(107, 2) = 0.00196779969601292;
        pointCoord(108, 2) = 0.00196779969601292;
        pointCoord(109, 2) = 0.00206865863732041;
        pointCoord(110, 2) = 0.00206865863732041;
        pointCoord(111, 2) = 0.00206865863732041;
        pointCoord(112, 2) = 0.00351562500000000;
        pointCoord(113, 2) = 0.00196779969601292;
        pointCoord(114, 2) = 0.00196779969601292;
        pointCoord(115, 2) = 0.00196779969601292;
        pointCoord(116, 2) = 0.00206865863732041;
        pointCoord(117, 2) = 0.00206865863732041;
        pointCoord(118, 2) = 0.00206865863732041;
        pointCoord(119, 2) = 0.00351562500000000;
        pointCoord(120, 2) = 0.00196779969601292;
        pointCoord(121, 2) = 0.00196779969601292;
        pointCoord(122, 2) = 0.00196779969601292;
        pointCoord(123, 2) = 0.00206865863732041;
        pointCoord(124, 2) = 0.00206865863732041;
        pointCoord(125, 2) = 0.00206865863732041;
        pointCoord(126, 2) = 0.00351562500000000;
        pointCoord(127, 2) = 0.00196779969601292;
        pointCoord(128, 2) = 0.00196779969601292;
        pointCoord(129, 2) = 0.00196779969601292;
        pointCoord(130, 2) = 0.00206865863732041;
        pointCoord(131, 2) = 0.00206865863732041;
        pointCoord(132, 2) = 0.00206865863732041;
        pointCoord(133, 2) = 0.00351562500000000;
        pointCoord(134, 2) = 0.00196779969601292;
        pointCoord(135, 2) = 0.00196779969601292;
        pointCoord(136, 2) = 0.00196779969601292;
        pointCoord(137, 2) = 0.00206865863732041;
        pointCoord(138, 2) = 0.00206865863732041;
        pointCoord(139, 2) = 0.00206865863732041;
        pointCoord(140, 2) = 0.00351562500000000;
        pointCoord(141, 2) = 0.00196779969601292;
        pointCoord(142, 2) = 0.00196779969601292;
        pointCoord(143, 2) = 0.00196779969601292;
        pointCoord(144, 2) = 0.00206865863732041;
        pointCoord(145, 2) = 0.00206865863732041;
        pointCoord(146, 2) = 0.00206865863732041;
        pointCoord(147, 2) = 0.00351562500000000;
        pointCoord(148, 2) = 0.00196779969601292;
        pointCoord(149, 2) = 0.00196779969601292;
        pointCoord(150, 2) = 0.00196779969601292;
        pointCoord(151, 2) = 0.00206865863732041;
        pointCoord(152, 2) = 0.00206865863732041;
        pointCoord(153, 2) = 0.00206865863732041;
        pointCoord(154, 2) = 0.00351562500000000;
        pointCoord(155, 2) = 0.00196779969601292;
        pointCoord(156, 2) = 0.00196779969601292;
        pointCoord(157, 2) = 0.00196779969601292;
        pointCoord(158, 2) = 0.00206865863732041;
        pointCoord(159, 2) = 0.00206865863732041;
        pointCoord(160, 2) = 0.00206865863732041;
        pointCoord(161, 2) = 0.00351562500000000;
        pointCoord(162, 2) = 0.00196779969601292;
        pointCoord(163, 2) = 0.00196779969601292;
        pointCoord(164, 2) = 0.00196779969601292;
        pointCoord(165, 2) = 0.00206865863732041;
        pointCoord(166, 2) = 0.00206865863732041;
        pointCoord(167, 2) = 0.00206865863732041;
        pointCoord(168, 2) = 0.00351562500000000;
        pointCoord(169, 2) = 0.00196779969601292;
        pointCoord(170, 2) = 0.00196779969601292;
        pointCoord(171, 2) = 0.00196779969601292;
        pointCoord(172, 2) = 0.00206865863732041;
        pointCoord(173, 2) = 0.00206865863732041;
        pointCoord(174, 2) = 0.00206865863732041;
        pointCoord(175, 2) = 0.00351562500000000;
        pointCoord(176, 2) = 0.00196779969601292;
        pointCoord(177, 2) = 0.00196779969601292;
        pointCoord(178, 2) = 0.00196779969601292;
        pointCoord(179, 2) = 0.00206865863732041;
        pointCoord(180, 2) = 0.00206865863732041;
        pointCoord(181, 2) = 0.00206865863732041;
        pointCoord(182, 2) = 0.00351562500000000;
        pointCoord(183, 2) = 0.00196779969601292;
        pointCoord(184, 2) = 0.00196779969601292;
        pointCoord(185, 2) = 0.00196779969601292;
        pointCoord(186, 2) = 0.00206865863732041;
        pointCoord(187, 2) = 0.00206865863732041;
        pointCoord(188, 2) = 0.00206865863732041;
        pointCoord(189, 2) = 0.00351562500000000;
        pointCoord(190, 2) = 0.00196779969601292;
        pointCoord(191, 2) = 0.00196779969601292;
        pointCoord(192, 2) = 0.00196779969601292;
        pointCoord(193, 2) = 0.00206865863732041;
        pointCoord(194, 2) = 0.00206865863732041;
        pointCoord(195, 2) = 0.00206865863732041;
        pointCoord(196, 2) = 0.00351562500000000;
        pointCoord(197, 2) = 0.00196779969601292;
        pointCoord(198, 2) = 0.00196779969601292;
        pointCoord(199, 2) = 0.00196779969601292;
        pointCoord(200, 2) = 0.00206865863732041;
        pointCoord(201, 2) = 0.00206865863732041;
        pointCoord(202, 2) = 0.00206865863732041;
        pointCoord(203, 2) = 0.00351562500000000;
        pointCoord(204, 2) = 0.00196779969601292;
        pointCoord(205, 2) = 0.00196779969601292;
        pointCoord(206, 2) = 0.00196779969601292;
        pointCoord(207, 2) = 0.00206865863732041;
        pointCoord(208, 2) = 0.00206865863732041;
        pointCoord(209, 2) = 0.00206865863732041;
        pointCoord(210, 2) = 0.00351562500000000;
        pointCoord(211, 2) = 0.00196779969601292;
        pointCoord(212, 2) = 0.00196779969601292;
        pointCoord(213, 2) = 0.00196779969601292;
        pointCoord(214, 2) = 0.00206865863732041;
        pointCoord(215, 2) = 0.00206865863732041;
        pointCoord(216, 2) = 0.00206865863732041;
        pointCoord(217, 2) = 0.00351562500000000;
        pointCoord(218, 2) = 0.00196779969601292;
        pointCoord(219, 2) = 0.00196779969601292;
        pointCoord(220, 2) = 0.00196779969601292;
        pointCoord(221, 2) = 0.00206865863732041;
        pointCoord(222, 2) = 0.00206865863732041;
        pointCoord(223, 2) = 0.00206865863732041;
    }
    else if (nH == 4)
    {
        pointCoord(0, 0) = 1.0 / 3.0;
        pointCoord(0, 1) = 1.0 / 3.0;
        pointCoord(0, 2) = -27.0 / 96.0;
        pointCoord(1, 0) = 0.6;
        pointCoord(1, 1) = 0.2;
        pointCoord(1, 2) = 25.0 / 96.0;
        pointCoord(2, 0) = 0.2;
        pointCoord(2, 1) = 0.6;
        pointCoord(2, 2) = 25.0 / 96.0;
        pointCoord(3, 0) = 0.2;
        pointCoord(3, 1) = 0.2;
        pointCoord(3, 2) = 25.0 / 96.0;
    }
    else if (nH == 25)
    {
        pointCoord(0, 0) = 2.200555327023207e-03;
        pointCoord(0, 1) = 4.470952170364481e-02;
        pointCoord(0, 2) = 6.583166573007299e-04;
        pointCoord(1, 0) = 1.082522010747988e-02;
        pointCoord(1, 1) = 3.608485692318814e-02;
        pointCoord(1, 2) = 1.329900683819439e-03;
        pointCoord(2, 0) = 2.345503851533401e-02;
        pointCoord(2, 1) = 2.345503851533401e-02;
        pointCoord(2, 2) = 1.580694532070693e-03;
        pointCoord(3, 0) = 3.608485692318814e-02;
        pointCoord(3, 1) = 1.082522010747988e-02;
        pointCoord(3, 2) = 1.329900683819439e-03;
        pointCoord(4, 0) = 4.470952170364481e-02;
        pointCoord(4, 1) = 2.200555327023207e-03;
        pointCoord(4, 2) = 6.583166573007299e-04;
        pointCoord(5, 0) = 1.082522010747988e-02;
        pointCoord(5, 1) = 2.199401248396786e-01;
        pointCoord(5, 2) = 6.542197529251943e-03;
        pointCoord(6, 0) = 5.325264442858103e-02;
        pointCoord(6, 1) = 1.775127005185774e-01;
        pointCoord(6, 2) = 1.321624308202714e-02;
        pointCoord(7, 0) = 1.153826724735792e-01;
        pointCoord(7, 1) = 1.153826724735792e-01;
        pointCoord(7, 2) = 1.570857390213492e-02;
        pointCoord(8, 0) = 1.775127005185774e-01;
        pointCoord(8, 1) = 5.325264442858103e-02;
        pointCoord(8, 2) = 1.321624308202714e-02;
        pointCoord(9, 0) = 2.199401248396786e-01;
        pointCoord(9, 1) = 1.082522010747988e-02;
        pointCoord(9, 2) = 6.542197529251943e-03;
        pointCoord(10, 0) = 2.345503851533401e-02;
        pointCoord(10, 1) = 4.765449614846660e-01;
        pointCoord(10, 2) = 1.684813404844011e-02;
        pointCoord(11, 0) = 1.153826724735792e-01;
        pointCoord(11, 1) = 3.846173275264207e-01;
        pointCoord(11, 2) = 3.403581656884384e-02;
        pointCoord(12, 0) = 2.500000000000000e-01;
        pointCoord(12, 1) = 2.500000000000000e-01;
        pointCoord(12, 2) = 4.045432098765432e-02;
        pointCoord(13, 0) = 3.846173275264207e-01;
        pointCoord(13, 1) = 1.153826724735792e-01;
        pointCoord(13, 2) = 3.403581656884384e-02;
        pointCoord(14, 0) = 4.765449614846660e-01;
        pointCoord(14, 1) = 2.345503851533401e-02;
        pointCoord(14, 2) = 1.684813404844011e-02;
        pointCoord(15, 0) = 3.608485692318814e-02;
        pointCoord(15, 1) = 7.331497981296533e-01;
        pointCoord(15, 2) = 2.180780247074806e-02;
        pointCoord(16, 0) = 1.775127005185774e-01;
        pointCoord(16, 1) = 5.917219545342640e-01;
        pointCoord(16, 2) = 4.405510797397064e-02;
        pointCoord(17, 0) = 3.846173275264207e-01;
        pointCoord(17, 1) = 3.846173275264207e-01;
        pointCoord(17, 2) = 5.236305923555275e-02;
        pointCoord(18, 0) = 5.917219545342640e-01;
        pointCoord(18, 1) = 1.775127005185774e-01;
        pointCoord(18, 2) = 4.405510797397064e-02;
        pointCoord(19, 0) = 7.331497981296533e-01;
        pointCoord(19, 1) = 3.608485692318814e-02;
        pointCoord(19, 2) = 2.180780247074806e-02;
        pointCoord(20, 0) = 4.470952170364481e-02;
        pointCoord(20, 1) = 9.083804012656871e-01;
        pointCoord(20, 2) = 1.337527055830643e-02;
        pointCoord(21, 0) = 2.199401248396786e-01;
        pointCoord(21, 1) = 7.331497981296533e-01;
        pointCoord(21, 2) = 2.702009931618056e-02;
        pointCoord(22, 0) = 4.765449614846660e-01;
        pointCoord(22, 1) = 4.765449614846660e-01;
        pointCoord(22, 2) = 3.211557356480953e-02;
        pointCoord(23, 0) = 7.331497981296533e-01;
        pointCoord(23, 1) = 2.199401248396786e-01;
        pointCoord(23, 2) = 2.702009931618056e-02;
        pointCoord(24, 0) = 9.083804012656871e-01;
        pointCoord(24, 1) = 4.470952170364481e-02;
        pointCoord(24, 2) = 1.337527055830643e-02;
    }
    else if (nH == 16)
    {
        pointCoord(0, 0) = 4.820780989426014e-03;
        pointCoord(0, 1) = 6.461106321354769e-02;
        pointCoord(0, 2) = 2.100365244474847e-03;
        pointCoord(1, 0) = 2.291316667641278e-02;
        pointCoord(1, 1) = 4.651867752656094e-02;
        pointCoord(1, 2) = 3.937685608733463e-03;
        pointCoord(2, 0) = 4.651867752656094e-02;
        pointCoord(2, 1) = 2.291316667641278e-02;
        pointCoord(2, 2) = 3.937685608733463e-03;
        pointCoord(3, 0) = 6.461106321354769e-02;
        pointCoord(3, 1) = 4.820780989426014e-03;
        pointCoord(3, 2) = 2.100365244474847e-03;
        pointCoord(4, 0) = 2.291316667641278e-02;
        pointCoord(4, 1) = 3.070963115311591e-01;
        pointCoord(4, 2) = 1.871581531501276e-02;
        pointCoord(5, 0) = 1.089062557068338e-01;
        pointCoord(5, 1) = 2.211032225007380e-01;
        pointCoord(5, 2) = 3.508770525293353e-02;
        pointCoord(6, 0) = 2.211032225007380e-01;
        pointCoord(6, 1) = 1.089062557068338e-01;
        pointCoord(6, 2) = 3.508770525293353e-02;
        pointCoord(7, 0) = 3.070963115311591e-01;
        pointCoord(7, 1) = 2.291316667641278e-02;
        pointCoord(7, 2) = 1.871581531501276e-02;
        pointCoord(8, 0) = 4.651867752656094e-02;
        pointCoord(8, 1) = 6.234718442658671e-01;
        pointCoord(8, 2) = 3.799714764795021e-02;
        pointCoord(9, 0) = 2.211032225007380e-01;
        pointCoord(9, 1) = 4.488872992916901e-01;
        pointCoord(9, 2) = 7.123562049974014e-02;
        pointCoord(10, 0) = 4.488872992916901e-01;
        pointCoord(10, 1) = 2.211032225007380e-01;
        pointCoord(10, 2) = 7.123562049974014e-02;
        pointCoord(11, 0) = 6.234718442658671e-01;
        pointCoord(11, 1) = 4.651867752656094e-02;
        pointCoord(11, 2) = 3.799714764795021e-02;
        pointCoord(12, 0) = 6.461106321354769e-02;
        pointCoord(12, 1) = 8.659570925834785e-01;
        pointCoord(12, 2) = 2.815038307692563e-02;
        pointCoord(13, 0) = 3.070963115311591e-01;
        pointCoord(13, 1) = 6.234718442658671e-01;
        pointCoord(13, 2) = 5.277527735422950e-02;
        pointCoord(14, 0) = 6.234718442658671e-01;
        pointCoord(14, 1) = 3.070963115311591e-01;
        pointCoord(14, 2) = 5.277527735422950e-02;
        pointCoord(15, 0) = 8.659570925834785e-01;
        pointCoord(15, 1) = 6.461106321354769e-02;
        pointCoord(15, 2) = 2.815038307692563e-02;
    }
    else if (nH == 9)
    {
        pointCoord(0, 0) = 1.270166537925831e-02;
        pointCoord(0, 1) = 9.999999999999999e-02;
        pointCoord(0, 2) = 8.696116155806878e-03;
        pointCoord(1, 0) = 5.635083268962915e-02;
        pointCoord(1, 1) = 5.635083268962915e-02;
        pointCoord(1, 2) = 1.391378584929107e-02;
        pointCoord(2, 0) = 9.999999999999999e-02;
        pointCoord(2, 1) = 1.270166537925831e-02;
        pointCoord(2, 2) = 8.696116155806878e-03;
        pointCoord(3, 0) = 5.635083268962915e-02;
        pointCoord(3, 1) = 4.436491673103709e-01;
        pointCoord(3, 2) = 6.172839506172807e-02;
        pointCoord(4, 0) = 2.500000000000000e-01;
        pointCoord(4, 1) = 2.500000000000000e-01;
        pointCoord(4, 2) = 9.876543209876543e-02;
        pointCoord(5, 0) = 4.436491673103709e-01;
        pointCoord(5, 1) = 5.635083268962915e-02;
        pointCoord(5, 2) = 6.172839506172807e-02;
        pointCoord(6, 0) = 9.999999999999999e-02;
        pointCoord(6, 1) = 7.872983346207417e-01;
        pointCoord(6, 2) = 6.846437767135283e-02;
        pointCoord(7, 0) = 4.436491673103709e-01;
        pointCoord(7, 1) = 4.436491673103709e-01;
        pointCoord(7, 2) = 1.095430042741651e-01;
        pointCoord(8, 0) = 7.872983346207417e-01;
        pointCoord(8, 1) = 9.999999999999999e-02;
        pointCoord(8, 2) = 6.846437767135283e-02;
    }
    else if (nH == 192)
    {
        pointCoord(0, 0) = 1.157197403381969e-02;
        pointCoord(0, 1) = 2.234048456942410e-02;
        pointCoord(0, 2) = 1.622253649952248e-03;
        pointCoord(1, 0) = 1.157197403381968e-02;
        pointCoord(1, 1) = 1.061842982897120e-01;
        pointCoord(1, 2) = 3.041340008804794e-03;
        pointCoord(2, 0) = 1.157197403381968e-02;
        pointCoord(2, 1) = 2.155770610095427e-01;
        pointCoord(2, 2) = 3.041340008809772e-03;
        pointCoord(3, 0) = 1.157197403381968e-02;
        pointCoord(3, 1) = 2.994208747300729e-01;
        pointCoord(3, 2) = 1.622253649958975e-03;
        pointCoord(4, 0) = 5.500157970124948e-02;
        pointCoord(4, 1) = 1.932508695492488e-02;
        pointCoord(4, 2) = 2.630836405859374e-03;
        pointCoord(5, 0) = 5.500157970124944e-02;
        pointCoord(5, 1) = 9.185211678468351e-02;
        pointCoord(5, 2) = 4.932192951453974e-03;
        pointCoord(6, 0) = 5.500157970124945e-02;
        pointCoord(6, 1) = 1.864796368472894e-01;
        pointCoord(6, 2) = 4.932192951460071e-03;
        pointCoord(7, 0) = 5.500157970124946e-02;
        pointCoord(7, 1) = 2.590066666772063e-01;
        pointCoord(7, 2) = 2.630836405867609e-03;
        pointCoord(8, 0) = 1.116650869654173e-01;
        pointCoord(8, 1) = 1.539083514655823e-02;
        pointCoord(8, 2) = 2.095243841054768e-03;
        pointCoord(9, 0) = 1.116650869654172e-01;
        pointCoord(9, 1) = 7.315262231904138e-02;
        pointCoord(9, 2) = 3.928084194600115e-03;
        pointCoord(10, 0) = 1.116650869654172e-01;
        pointCoord(10, 1) = 1.485156240488007e-01;
        pointCoord(10, 2) = 3.928084194599062e-03;
        pointCoord(11, 0) = 1.116650869654172e-01;
        pointCoord(11, 1) = 2.062774112212564e-01;
        pointCoord(11, 2) = 2.095243841053343e-03;
        pointCoord(12, 0) = 1.550946926328470e-01;
        pointCoord(12, 1) = 1.237543753205690e-02;
        pointCoord(12, 2) = 8.986420434954288e-04;
        pointCoord(13, 0) = 1.550946926328470e-01;
        pointCoord(13, 1) = 5.882044081397868e-02;
        pointCoord(13, 2) = 1.684740238101779e-03;
        pointCoord(14, 0) = 1.550946926328470e-01;
        pointCoord(14, 1) = 1.194181998864780e-01;
        pointCoord(14, 2) = 1.684740238100998e-03;
        pointCoord(15, 0) = 1.550946926328470e-01;
        pointCoord(15, 1) = 1.658632031683617e-01;
        pointCoord(15, 2) = 8.986420434943734e-04;
        pointCoord(16, 0) = 1.157197403381969e-02;
        pointCoord(16, 1) = 3.341367968315672e-01;
        pointCoord(16, 2) = 8.986420434905156e-04;
        pointCoord(17, 0) = 1.157197403381968e-02;
        pointCoord(17, 1) = 3.805818001134422e-01;
        pointCoord(17, 2) = 1.684740238105775e-03;
        pointCoord(18, 0) = 1.157197403381968e-02;
        pointCoord(18, 1) = 4.411795591861550e-01;
        pointCoord(18, 2) = 1.684740238105773e-03;
        pointCoord(19, 0) = 1.157197403381968e-02;
        pointCoord(19, 1) = 4.876245624680299e-01;
        pointCoord(19, 2) = 8.986420434905170e-04;
        pointCoord(20, 0) = 5.500157970124948e-02;
        pointCoord(20, 1) = 2.937225887787070e-01;
        pointCoord(20, 2) = 2.095243841049653e-03;
        pointCoord(21, 0) = 5.500157970124944e-02;
        pointCoord(21, 1) = 3.514843759511614e-01;
        pointCoord(21, 2) = 3.928084194603993e-03;
        pointCoord(22, 0) = 5.500157970124945e-02;
        pointCoord(22, 1) = 4.268473776810351e-01;
        pointCoord(22, 2) = 3.928084194603994e-03;
        pointCoord(23, 0) = 5.500157970124946e-02;
        pointCoord(23, 1) = 4.846091648534898e-01;
        pointCoord(23, 2) = 2.095243841049653e-03;
        pointCoord(24, 0) = 1.116650869654173e-01;
        pointCoord(24, 1) = 2.409933333227800e-01;
        pointCoord(24, 2) = 2.630836405865475e-03;
        pointCoord(25, 0) = 1.116650869654172e-01;
        pointCoord(25, 1) = 3.135203631526244e-01;
        pointCoord(25, 2) = 4.932192951455038e-03;
        pointCoord(26, 0) = 1.116650869654172e-01;
        pointCoord(26, 1) = 4.081478832151789e-01;
        pointCoord(26, 2) = 4.932192951455039e-03;
        pointCoord(27, 0) = 1.116650869654172e-01;
        pointCoord(27, 1) = 4.806749130450234e-01;
        pointCoord(27, 2) = 2.630836405865472e-03;
        pointCoord(28, 0) = 1.550946926328470e-01;
        pointCoord(28, 1) = 2.005791252699198e-01;
        pointCoord(28, 2) = 1.622253649957981e-03;
        pointCoord(29, 0) = 1.550946926328470e-01;
        pointCoord(29, 1) = 2.844229389903439e-01;
        pointCoord(29, 2) = 3.041340008804913e-03;
        pointCoord(30, 0) = 1.550946926328470e-01;
        pointCoord(30, 1) = 3.938157017100590e-01;
        pointCoord(30, 2) = 3.041340008804911e-03;
        pointCoord(31, 0) = 1.550946926328470e-01;
        pointCoord(31, 1) = 4.776595154304831e-01;
        pointCoord(31, 2) = 1.622253649957981e-03;
        pointCoord(32, 0) = 1.157197403381969e-02;
        pointCoord(32, 1) = 5.223404845695169e-01;
        pointCoord(32, 2) = 1.622253649957981e-03;
        pointCoord(33, 0) = 1.157197403381968e-02;
        pointCoord(33, 1) = 6.061842982899409e-01;
        pointCoord(33, 2) = 3.041340008804914e-03;
        pointCoord(34, 0) = 1.157197403381968e-02;
        pointCoord(34, 1) = 7.155770610096562e-01;
        pointCoord(34, 2) = 3.041340008804914e-03;
        pointCoord(35, 0) = 1.157197403381968e-02;
        pointCoord(35, 1) = 7.994208747300804e-01;
        pointCoord(35, 2) = 1.622253649957981e-03;
        pointCoord(36, 0) = 5.500157970124948e-02;
        pointCoord(36, 1) = 5.193250869549767e-01;
        pointCoord(36, 2) = 2.630836405865474e-03;
        pointCoord(37, 0) = 5.500157970124944e-02;
        pointCoord(37, 1) = 5.918521167848209e-01;
        pointCoord(37, 2) = 4.932192951455039e-03;
        pointCoord(38, 0) = 5.500157970124945e-02;
        pointCoord(38, 1) = 6.864796368473756e-01;
        pointCoord(38, 2) = 4.932192951455043e-03;
        pointCoord(39, 0) = 5.500157970124946e-02;
        pointCoord(39, 1) = 7.590066666772201e-01;
        pointCoord(39, 2) = 2.630836405865474e-03;
        pointCoord(40, 0) = 1.116650869654173e-01;
        pointCoord(40, 1) = 5.153908351465103e-01;
        pointCoord(40, 2) = 2.095243841049657e-03;
        pointCoord(41, 0) = 1.116650869654172e-01;
        pointCoord(41, 1) = 5.731526223189646e-01;
        pointCoord(41, 2) = 3.928084194603990e-03;
        pointCoord(42, 0) = 1.116650869654172e-01;
        pointCoord(42, 1) = 6.485156240488386e-01;
        pointCoord(42, 2) = 3.928084194603990e-03;
        pointCoord(43, 0) = 1.116650869654172e-01;
        pointCoord(43, 1) = 7.062774112212932e-01;
        pointCoord(43, 2) = 2.095243841049653e-03;
        pointCoord(44, 0) = 1.550946926328470e-01;
        pointCoord(44, 1) = 5.123754375319701e-01;
        pointCoord(44, 2) = 8.986420434905160e-04;
        pointCoord(45, 0) = 1.550946926328470e-01;
        pointCoord(45, 1) = 5.588204408138451e-01;
        pointCoord(45, 2) = 1.684740238105771e-03;
        pointCoord(46, 0) = 1.550946926328470e-01;
        pointCoord(46, 1) = 6.194181998865580e-01;
        pointCoord(46, 2) = 1.684740238105775e-03;
        pointCoord(47, 0) = 1.550946926328470e-01;
        pointCoord(47, 1) = 6.658632031684328e-01;
        pointCoord(47, 2) = 8.986420434905156e-04;
        pointCoord(48, 0) = 1.237543753205690e-02;
        pointCoord(48, 1) = 8.325298698350996e-01;
        pointCoord(48, 2) = 8.986420434901241e-04;
        pointCoord(49, 0) = 1.539083514655823e-02;
        pointCoord(49, 1) = 8.729440778879596e-01;
        pointCoord(49, 2) = 2.095243841057306e-03;
        pointCoord(50, 0) = 1.932508695492488e-02;
        pointCoord(50, 1) = 9.256733333438868e-01;
        pointCoord(50, 2) = 2.630836405863360e-03;
        pointCoord(51, 0) = 2.234048456942410e-02;
        pointCoord(51, 1) = 9.660875413967470e-01;
        pointCoord(51, 2) = 1.622253649943573e-03;
        pointCoord(52, 0) = 5.882044081397868e-02;
        pointCoord(52, 1) = 7.860848665532245e-01;
        pointCoord(52, 2) = 1.684740238095509e-03;
        pointCoord(53, 0) = 7.315262231904136e-02;
        pointCoord(53, 1) = 8.151822907155051e-01;
        pointCoord(53, 2) = 3.928084194600797e-03;
        pointCoord(54, 0) = 9.185211678468352e-02;
        pointCoord(54, 1) = 8.531463035140421e-01;
        pointCoord(54, 2) = 4.932192951451356e-03;
        pointCoord(55, 0) = 1.061842982897120e-01;
        pointCoord(55, 1) = 8.822437276763228e-01;
        pointCoord(55, 2) = 3.041340008787695e-03;
        pointCoord(56, 0) = 1.194181998864780e-01;
        pointCoord(56, 1) = 7.254871074805118e-01;
        pointCoord(56, 2) = 1.684740238098923e-03;
        pointCoord(57, 0) = 1.485156240488006e-01;
        pointCoord(57, 1) = 7.398192889856315e-01;
        pointCoord(57, 2) = 3.928084194613241e-03;
        pointCoord(58, 0) = 1.864796368472894e-01;
        pointCoord(58, 1) = 7.585187834514878e-01;
        pointCoord(58, 2) = 4.932192951465895e-03;
        pointCoord(59, 0) = 2.155770610095426e-01;
        pointCoord(59, 1) = 7.728509649566077e-01;
        pointCoord(59, 2) = 3.041340008789994e-03;
        pointCoord(60, 0) = 1.658632031683617e-01;
        pointCoord(60, 1) = 6.790421041986368e-01;
        pointCoord(60, 2) = 8.986420434867307e-04;
        pointCoord(61, 0) = 2.062774112212564e-01;
        pointCoord(61, 1) = 6.820575018131768e-01;
        pointCoord(61, 2) = 2.095243841057333e-03;
        pointCoord(62, 0) = 2.590066666772065e-01;
        pointCoord(62, 1) = 6.859917536216433e-01;
        pointCoord(62, 2) = 2.630836405873973e-03;
        pointCoord(63, 0) = 2.994208747300729e-01;
        pointCoord(63, 1) = 6.890071512361834e-01;
        pointCoord(63, 2) = 1.622253649949223e-03;
        pointCoord(64, 0) = 1.890071512361835e-01;
        pointCoord(64, 1) = 1.157197403381968e-02;
        pointCoord(64, 2) = 1.622253649957981e-03;
        pointCoord(65, 0) = 1.859917536216433e-01;
        pointCoord(65, 1) = 5.500157970124947e-02;
        pointCoord(65, 2) = 2.630836405865474e-03;
        pointCoord(66, 0) = 1.820575018131769e-01;
        pointCoord(66, 1) = 1.116650869654172e-01;
        pointCoord(66, 2) = 2.095243841049651e-03;
        pointCoord(67, 0) = 1.790421041986367e-01;
        pointCoord(67, 1) = 1.550946926328470e-01;
        pointCoord(67, 2) = 8.986420434905154e-04;
        pointCoord(68, 0) = 2.728509649566077e-01;
        pointCoord(68, 1) = 1.157197403381968e-02;
        pointCoord(68, 2) = 3.041340008804913e-03;
        pointCoord(69, 0) = 2.585187834514877e-01;
        pointCoord(69, 1) = 5.500157970124944e-02;
        pointCoord(69, 2) = 4.932192951455040e-03;
        pointCoord(70, 0) = 2.398192889856315e-01;
        pointCoord(70, 1) = 1.116650869654172e-01;
        pointCoord(70, 2) = 3.928084194603993e-03;
        pointCoord(71, 0) = 2.254871074805117e-01;
        pointCoord(71, 1) = 1.550946926328470e-01;
        pointCoord(71, 2) = 1.684740238105774e-03;
        pointCoord(72, 0) = 3.822437276763229e-01;
        pointCoord(72, 1) = 1.157197403381968e-02;
        pointCoord(72, 2) = 3.041340008804913e-03;
        pointCoord(73, 0) = 3.531463035140420e-01;
        pointCoord(73, 1) = 5.500157970124944e-02;
        pointCoord(73, 2) = 4.932192951455039e-03;
        pointCoord(74, 0) = 3.151822907155051e-01;
        pointCoord(74, 1) = 1.116650869654172e-01;
        pointCoord(74, 2) = 3.928084194603992e-03;
        pointCoord(75, 0) = 2.860848665532246e-01;
        pointCoord(75, 1) = 1.550946926328470e-01;
        pointCoord(75, 2) = 1.684740238105774e-03;
        pointCoord(76, 0) = 4.660875413967471e-01;
        pointCoord(76, 1) = 1.157197403381968e-02;
        pointCoord(76, 2) = 1.622253649957980e-03;
        pointCoord(77, 0) = 4.256733333438866e-01;
        pointCoord(77, 1) = 5.500157970124945e-02;
        pointCoord(77, 2) = 2.630836405865472e-03;
        pointCoord(78, 0) = 3.729440778879598e-01;
        pointCoord(78, 1) = 1.116650869654172e-01;
        pointCoord(78, 2) = 2.095243841049652e-03;
        pointCoord(79, 0) = 3.325298698350995e-01;
        pointCoord(79, 1) = 1.550946926328470e-01;
        pointCoord(79, 2) = 8.986420434905160e-04;
        pointCoord(80, 0) = 1.782386407004029e-01;
        pointCoord(80, 1) = 1.890071512361835e-01;
        pointCoord(80, 2) = 1.622253649949041e-03;
        pointCoord(81, 0) = 1.782386407004029e-01;
        pointCoord(81, 1) = 2.728509649566077e-01;
        pointCoord(81, 2) = 3.041340008788154e-03;
        pointCoord(82, 0) = 1.782386407004029e-01;
        pointCoord(82, 1) = 3.822437276763229e-01;
        pointCoord(82, 2) = 3.041340008788151e-03;
        pointCoord(83, 0) = 1.782386407004029e-01;
        pointCoord(83, 1) = 4.660875413967469e-01;
        pointCoord(83, 2) = 1.622253649949040e-03;
        pointCoord(84, 0) = 2.216682463678034e-01;
        pointCoord(84, 1) = 1.859917536216433e-01;
        pointCoord(84, 2) = 2.630836405873208e-03;
        pointCoord(85, 0) = 2.216682463678033e-01;
        pointCoord(85, 1) = 2.585187834514877e-01;
        pointCoord(85, 2) = 4.932192951469537e-03;
        pointCoord(86, 0) = 2.216682463678033e-01;
        pointCoord(86, 1) = 3.531463035140421e-01;
        pointCoord(86, 2) = 4.932192951469540e-03;
        pointCoord(87, 0) = 2.216682463678034e-01;
        pointCoord(87, 1) = 4.256733333438866e-01;
        pointCoord(87, 2) = 2.630836405873206e-03;
        pointCoord(88, 0) = 2.783317536321966e-01;
        pointCoord(88, 1) = 1.820575018131769e-01;
        pointCoord(88, 2) = 2.095243841055812e-03;
        pointCoord(89, 0) = 2.783317536321966e-01;
        pointCoord(89, 1) = 2.398192889856314e-01;
        pointCoord(89, 2) = 3.928084194615539e-03;
        pointCoord(90, 0) = 2.783317536321966e-01;
        pointCoord(90, 1) = 3.151822907155051e-01;
        pointCoord(90, 2) = 3.928084194615539e-03;
        pointCoord(91, 0) = 2.783317536321967e-01;
        pointCoord(91, 1) = 3.729440778879598e-01;
        pointCoord(91, 2) = 2.095243841055812e-03;
        pointCoord(92, 0) = 3.217613592995972e-01;
        pointCoord(92, 1) = 1.790421041986368e-01;
        pointCoord(92, 2) = 8.986420434855628e-04;
        pointCoord(93, 0) = 3.217613592995971e-01;
        pointCoord(93, 1) = 2.254871074805117e-01;
        pointCoord(93, 2) = 1.684740238096489e-03;
        pointCoord(94, 0) = 3.217613592995971e-01;
        pointCoord(94, 1) = 2.860848665532246e-01;
        pointCoord(94, 2) = 1.684740238096490e-03;
        pointCoord(95, 0) = 3.217613592995971e-01;
        pointCoord(95, 1) = 3.325298698350995e-01;
        pointCoord(95, 2) = 8.986420434855651e-04;
        pointCoord(96, 0) = 1.790421041986367e-01;
        pointCoord(96, 1) = 4.991965365017663e-01;
        pointCoord(96, 2) = 8.986420434855649e-04;
        pointCoord(97, 0) = 1.820575018131770e-01;
        pointCoord(97, 1) = 5.396107445546265e-01;
        pointCoord(97, 2) = 2.095243841055814e-03;
        pointCoord(98, 0) = 1.859917536216432e-01;
        pointCoord(98, 1) = 5.923400000105535e-01;
        pointCoord(98, 2) = 2.630836405873210e-03;
        pointCoord(99, 0) = 1.890071512361835e-01;
        pointCoord(99, 1) = 6.327542080634135e-01;
        pointCoord(99, 2) = 1.622253649949043e-03;
        pointCoord(100, 0) = 2.254871074805117e-01;
        pointCoord(100, 1) = 4.527515332198913e-01;
        pointCoord(100, 2) = 1.684740238096492e-03;
        pointCoord(101, 0) = 2.398192889856314e-01;
        pointCoord(101, 1) = 4.818489573821716e-01;
        pointCoord(101, 2) = 3.928084194615543e-03;
        pointCoord(102, 0) = 2.585187834514877e-01;
        pointCoord(102, 1) = 5.198129701807088e-01;
        pointCoord(102, 2) = 4.932192951469535e-03;
        pointCoord(103, 0) = 2.728509649566076e-01;
        pointCoord(103, 1) = 5.489103943429896e-01;
        pointCoord(103, 2) = 3.041340008788149e-03;
        pointCoord(104, 0) = 2.860848665532246e-01;
        pointCoord(104, 1) = 3.921537741471786e-01;
        pointCoord(104, 2) = 1.684740238096498e-03;
        pointCoord(105, 0) = 3.151822907155051e-01;
        pointCoord(105, 1) = 4.064859556522981e-01;
        pointCoord(105, 2) = 3.928084194615536e-03;
        pointCoord(106, 0) = 3.531463035140421e-01;
        pointCoord(106, 1) = 4.251854501181545e-01;
        pointCoord(106, 2) = 4.932192951469534e-03;
        pointCoord(107, 0) = 3.822437276763229e-01;
        pointCoord(107, 1) = 4.395176316232743e-01;
        pointCoord(107, 2) = 3.041340008788155e-03;
        pointCoord(108, 0) = 3.325298698350996e-01;
        pointCoord(108, 1) = 3.457087708653034e-01;
        pointCoord(108, 2) = 8.986420434855642e-04;
        pointCoord(109, 0) = 3.729440778879599e-01;
        pointCoord(109, 1) = 3.487241684798436e-01;
        pointCoord(109, 2) = 2.095243841055809e-03;
        pointCoord(110, 0) = 4.256733333438867e-01;
        pointCoord(110, 1) = 3.526584202883100e-01;
        pointCoord(110, 2) = 2.630836405873202e-03;
        pointCoord(111, 0) = 4.660875413967469e-01;
        pointCoord(111, 1) = 3.556738179028502e-01;
        pointCoord(111, 2) = 1.622253649949044e-03;
        pointCoord(112, 0) = 2.005791252699198e-01;
        pointCoord(112, 1) = 6.443261820971499e-01;
        pointCoord(112, 2) = 1.622253649949039e-03;
        pointCoord(113, 0) = 2.409933333227800e-01;
        pointCoord(113, 1) = 6.473415797116902e-01;
        pointCoord(113, 2) = 2.630836405873203e-03;
        pointCoord(114, 0) = 2.937225887787069e-01;
        pointCoord(114, 1) = 6.512758315201566e-01;
        pointCoord(114, 2) = 2.095243841055811e-03;
        pointCoord(115, 0) = 3.341367968315672e-01;
        pointCoord(115, 1) = 6.542912291346966e-01;
        pointCoord(115, 2) = 8.986420434855642e-04;
        pointCoord(116, 0) = 2.844229389903439e-01;
        pointCoord(116, 1) = 5.604823683767259e-01;
        pointCoord(116, 2) = 3.041340008788156e-03;
        pointCoord(117, 0) = 3.135203631526244e-01;
        pointCoord(117, 1) = 5.748145498818454e-01;
        pointCoord(117, 2) = 4.932192951469537e-03;
        pointCoord(118, 0) = 3.514843759511614e-01;
        pointCoord(118, 1) = 5.935140443477018e-01;
        pointCoord(118, 2) = 3.928084194615538e-03;
        pointCoord(119, 0) = 3.805818001134422e-01;
        pointCoord(119, 1) = 6.078462258528217e-01;
        pointCoord(119, 2) = 1.684740238096484e-03;
        pointCoord(120, 0) = 3.938157017100593e-01;
        pointCoord(120, 1) = 4.510896056570105e-01;
        pointCoord(120, 2) = 3.041340008788162e-03;
        pointCoord(121, 0) = 4.081478832151789e-01;
        pointCoord(121, 1) = 4.801870298192910e-01;
        pointCoord(121, 2) = 4.932192951469535e-03;
        pointCoord(122, 0) = 4.268473776810351e-01;
        pointCoord(122, 1) = 5.181510426178281e-01;
        pointCoord(122, 2) = 3.928084194615541e-03;
        pointCoord(123, 0) = 4.411795591861550e-01;
        pointCoord(123, 1) = 5.472484667801089e-01;
        pointCoord(123, 2) = 1.684740238096494e-03;
        pointCoord(124, 0) = 4.776595154304833e-01;
        pointCoord(124, 1) = 3.672457919365865e-01;
        pointCoord(124, 2) = 1.622253649949044e-03;
        pointCoord(125, 0) = 4.806749130450235e-01;
        pointCoord(125, 1) = 4.076599999894468e-01;
        pointCoord(125, 2) = 2.630836405873204e-03;
        pointCoord(126, 0) = 4.846091648534898e-01;
        pointCoord(126, 1) = 4.603892554453735e-01;
        pointCoord(126, 2) = 2.095243841055809e-03;
        pointCoord(127, 0) = 4.876245624680299e-01;
        pointCoord(127, 1) = 5.008034634982338e-01;
        pointCoord(127, 2) = 8.986420434855657e-04;
        pointCoord(128, 0) = 3.449053073670695e-01;
        pointCoord(128, 1) = 1.674701301648327e-01;
        pointCoord(128, 2) = 8.986420434866303e-04;
        pointCoord(129, 0) = 3.449053073670696e-01;
        pointCoord(129, 1) = 2.139151334467493e-01;
        pointCoord(129, 2) = 1.684740238097554e-03;
        pointCoord(130, 0) = 3.449053073670696e-01;
        pointCoord(130, 1) = 2.745128925194848e-01;
        pointCoord(130, 2) = 1.684740238096771e-03;
        pointCoord(131, 0) = 3.449053073670696e-01;
        pointCoord(131, 1) = 3.209578958013633e-01;
        pointCoord(131, 2) = 8.986420434855741e-04;
        pointCoord(132, 0) = 3.883349130344700e-01;
        pointCoord(132, 1) = 1.270559221119396e-01;
        pointCoord(132, 2) = 2.095243841057444e-03;
        pointCoord(133, 0) = 3.883349130344700e-01;
        pointCoord(133, 1) = 1.848177092844367e-01;
        pointCoord(133, 2) = 3.928084194618207e-03;
        pointCoord(134, 0) = 3.883349130344700e-01;
        pointCoord(134, 1) = 2.601807110143525e-01;
        pointCoord(134, 2) = 3.928084194617158e-03;
        pointCoord(135, 0) = 3.883349130344700e-01;
        pointCoord(135, 1) = 3.179424981868223e-01;
        pointCoord(135, 2) = 2.095243841056018e-03;
        pointCoord(136, 0) = 4.449984202988633e-01;
        pointCoord(136, 1) = 7.432666665616368e-02;
        pointCoord(136, 2) = 2.630836405865751e-03;
        pointCoord(137, 0) = 4.449984202988632e-01;
        pointCoord(137, 1) = 1.468536964858919e-01;
        pointCoord(137, 2) = 4.932192951466768e-03;
        pointCoord(138, 0) = 4.449984202988632e-01;
        pointCoord(138, 1) = 2.414812165484672e-01;
        pointCoord(138, 2) = 4.932192951472861e-03;
        pointCoord(139, 0) = 4.449984202988633e-01;
        pointCoord(139, 1) = 3.140082463783538e-01;
        pointCoord(139, 2) = 2.630836405873983e-03;
        pointCoord(140, 0) = 4.884280259662638e-01;
        pointCoord(140, 1) = 3.391245860324224e-02;
        pointCoord(140, 2) = 1.622253649943122e-03;
        pointCoord(141, 0) = 4.884280259662638e-01;
        pointCoord(141, 1) = 1.177562723235098e-01;
        pointCoord(141, 2) = 3.041340008787009e-03;
        pointCoord(142, 0) = 4.884280259662639e-01;
        pointCoord(142, 1) = 2.271490350433006e-01;
        pointCoord(142, 2) = 3.041340008791988e-03;
        pointCoord(143, 0) = 4.884280259662637e-01;
        pointCoord(143, 1) = 3.109928487638107e-01;
        pointCoord(143, 2) = 1.622253649949850e-03;
        pointCoord(144, 0) = 5.115719740337362e-01;
        pointCoord(144, 1) = 2.234048456942410e-02;
        pointCoord(144, 2) = 1.622253649943309e-03;
        pointCoord(145, 0) = 5.115719740337363e-01;
        pointCoord(145, 1) = 1.061842982897120e-01;
        pointCoord(145, 2) = 3.041340008788037e-03;
        pointCoord(146, 0) = 5.115719740337363e-01;
        pointCoord(146, 1) = 2.155770610095427e-01;
        pointCoord(146, 2) = 3.041340008793014e-03;
        pointCoord(147, 0) = 5.115719740337362e-01;
        pointCoord(147, 1) = 2.994208747300729e-01;
        pointCoord(147, 2) = 1.622253649950036e-03;
        pointCoord(148, 0) = 5.550015797011367e-01;
        pointCoord(148, 1) = 1.932508695492488e-02;
        pointCoord(148, 2) = 2.630836405867108e-03;
        pointCoord(149, 0) = 5.550015797011366e-01;
        pointCoord(149, 1) = 9.185211678468351e-02;
        pointCoord(149, 2) = 4.932192951468473e-03;
        pointCoord(150, 0) = 5.550015797011366e-01;
        pointCoord(150, 1) = 1.864796368472894e-01;
        pointCoord(150, 2) = 4.932192951474568e-03;
        pointCoord(151, 0) = 5.550015797011367e-01;
        pointCoord(151, 1) = 2.590066666772063e-01;
        pointCoord(151, 2) = 2.630836405875342e-03;
        pointCoord(152, 0) = 6.116650869655300e-01;
        pointCoord(152, 1) = 1.539083514655823e-02;
        pointCoord(152, 2) = 2.095243841060927e-03;
        pointCoord(153, 0) = 6.116650869655298e-01;
        pointCoord(153, 1) = 7.315262231904138e-02;
        pointCoord(153, 2) = 3.928084194611662e-03;
        pointCoord(154, 0) = 6.116650869655300e-01;
        pointCoord(154, 1) = 1.485156240488007e-01;
        pointCoord(154, 2) = 3.928084194610608e-03;
        pointCoord(155, 0) = 6.116650869655301e-01;
        pointCoord(155, 1) = 2.062774112212564e-01;
        pointCoord(155, 2) = 2.095243841059504e-03;
        pointCoord(156, 0) = 6.550946926329305e-01;
        pointCoord(156, 1) = 1.237543753205690e-02;
        pointCoord(156, 2) = 8.986420434904767e-04;
        pointCoord(157, 0) = 6.550946926329305e-01;
        pointCoord(157, 1) = 5.882044081397868e-02;
        pointCoord(157, 2) = 1.684740238092495e-03;
        pointCoord(158, 0) = 6.550946926329305e-01;
        pointCoord(158, 1) = 1.194181998864780e-01;
        pointCoord(158, 2) = 1.684740238091714e-03;
        pointCoord(159, 0) = 6.550946926329304e-01;
        pointCoord(159, 1) = 1.658632031683617e-01;
        pointCoord(159, 2) = 8.986420434894220e-04;
        pointCoord(160, 0) = 5.123754375319700e-01;
        pointCoord(160, 1) = 3.325298698350996e-01;
        pointCoord(160, 2) = 8.986420434855647e-04;
        pointCoord(161, 0) = 5.153908351465103e-01;
        pointCoord(161, 1) = 3.729440778879600e-01;
        pointCoord(161, 2) = 2.095243841055814e-03;
        pointCoord(162, 0) = 5.193250869549766e-01;
        pointCoord(162, 1) = 4.256733333438867e-01;
        pointCoord(162, 2) = 2.630836405873206e-03;
        pointCoord(163, 0) = 5.223404845695169e-01;
        pointCoord(163, 1) = 4.660875413967469e-01;
        pointCoord(163, 2) = 1.622253649949039e-03;
        pointCoord(164, 0) = 5.588204408138450e-01;
        pointCoord(164, 1) = 2.860848665532246e-01;
        pointCoord(164, 2) = 1.684740238096494e-03;
        pointCoord(165, 0) = 5.731526223189647e-01;
        pointCoord(165, 1) = 3.151822907155051e-01;
        pointCoord(165, 2) = 3.928084194615540e-03;
        pointCoord(166, 0) = 5.918521167848211e-01;
        pointCoord(166, 1) = 3.531463035140422e-01;
        pointCoord(166, 2) = 4.932192951469534e-03;
        pointCoord(167, 0) = 6.061842982899409e-01;
        pointCoord(167, 1) = 3.822437276763229e-01;
        pointCoord(167, 2) = 3.041340008788153e-03;
        pointCoord(168, 0) = 6.194181998865580e-01;
        pointCoord(168, 1) = 2.254871074805117e-01;
        pointCoord(168, 2) = 1.684740238096494e-03;
        pointCoord(169, 0) = 6.485156240488382e-01;
        pointCoord(169, 1) = 2.398192889856314e-01;
        pointCoord(169, 2) = 3.928084194615536e-03;
        pointCoord(170, 0) = 6.864796368473755e-01;
        pointCoord(170, 1) = 2.585187834514878e-01;
        pointCoord(170, 2) = 4.932192951469540e-03;
        pointCoord(171, 0) = 7.155770610096562e-01;
        pointCoord(171, 1) = 2.728509649566076e-01;
        pointCoord(171, 2) = 3.041340008788154e-03;
        pointCoord(172, 0) = 6.658632031684329e-01;
        pointCoord(172, 1) = 1.790421041986367e-01;
        pointCoord(172, 2) = 8.986420434855635e-04;
        pointCoord(173, 0) = 7.062774112212931e-01;
        pointCoord(173, 1) = 1.820575018131769e-01;
        pointCoord(173, 2) = 2.095243841055808e-03;
        pointCoord(174, 0) = 7.590066666772201e-01;
        pointCoord(174, 1) = 1.859917536216433e-01;
        pointCoord(174, 2) = 2.630836405873206e-03;
        pointCoord(175, 0) = 7.994208747300803e-01;
        pointCoord(175, 1) = 1.890071512361835e-01;
        pointCoord(175, 2) = 1.622253649949043e-03;
        pointCoord(176, 0) = 6.890071512361834e-01;
        pointCoord(176, 1) = 1.157197403381968e-02;
        pointCoord(176, 2) = 1.622253649957982e-03;
        pointCoord(177, 0) = 6.859917536216433e-01;
        pointCoord(177, 1) = 5.500157970124947e-02;
        pointCoord(177, 2) = 2.630836405865475e-03;
        pointCoord(178, 0) = 6.820575018131770e-01;
        pointCoord(178, 1) = 1.116650869654172e-01;
        pointCoord(178, 2) = 2.095243841049653e-03;
        pointCoord(179, 0) = 6.790421041986366e-01;
        pointCoord(179, 1) = 1.550946926328470e-01;
        pointCoord(179, 2) = 8.986420434905149e-04;
        pointCoord(180, 0) = 7.728509649566077e-01;
        pointCoord(180, 1) = 1.157197403381968e-02;
        pointCoord(180, 2) = 3.041340008804914e-03;
        pointCoord(181, 0) = 7.585187834514876e-01;
        pointCoord(181, 1) = 5.500157970124944e-02;
        pointCoord(181, 2) = 4.932192951455040e-03;
        pointCoord(182, 0) = 7.398192889856315e-01;
        pointCoord(182, 1) = 1.116650869654172e-01;
        pointCoord(182, 2) = 3.928084194603994e-03;
        pointCoord(183, 0) = 7.254871074805117e-01;
        pointCoord(183, 1) = 1.550946926328470e-01;
        pointCoord(183, 2) = 1.684740238105776e-03;
        pointCoord(184, 0) = 8.822437276763228e-01;
        pointCoord(184, 1) = 1.157197403381968e-02;
        pointCoord(184, 2) = 3.041340008804913e-03;
        pointCoord(185, 0) = 8.531463035140420e-01;
        pointCoord(185, 1) = 5.500157970124944e-02;
        pointCoord(185, 2) = 4.932192951455038e-03;
        pointCoord(186, 0) = 8.151822907155051e-01;
        pointCoord(186, 1) = 1.116650869654172e-01;
        pointCoord(186, 2) = 3.928084194603993e-03;
        pointCoord(187, 0) = 7.860848665532245e-01;
        pointCoord(187, 1) = 1.550946926328470e-01;
        pointCoord(187, 2) = 1.684740238105773e-03;
        pointCoord(188, 0) = 9.660875413967469e-01;
        pointCoord(188, 1) = 1.157197403381968e-02;
        pointCoord(188, 2) = 1.622253649957980e-03;
        pointCoord(189, 0) = 9.256733333438866e-01;
        pointCoord(189, 1) = 5.500157970124945e-02;
        pointCoord(189, 2) = 2.630836405865475e-03;
        pointCoord(190, 0) = 8.729440778879599e-01;
        pointCoord(190, 1) = 1.116650869654172e-01;
        pointCoord(190, 2) = 2.095243841049651e-03;
        pointCoord(191, 0) = 8.325298698350996e-01;
        pointCoord(191, 1) = 1.550946926328470e-01;
        pointCoord(191, 2) = 8.986420434905168e-04;
    }
    else if (nH == 300)
    {
        pointCoord(0, 0) = 7.818346171771254e-03;
        pointCoord(0, 1) = 1.526993312238961e-02;
        pointCoord(0, 2) = 7.613571603859526e-04;
        pointCoord(1, 0) = 7.818346171771256e-03;
        pointCoord(1, 1) = 7.511757829776841e-02;
        pointCoord(1, 2) = 1.538058314335309e-03;
        pointCoord(2, 0) = 7.818346171771252e-03;
        pointCoord(2, 1) = 1.627574935806183e-01;
        pointCoord(2, 2) = 1.828106712821296e-03;
        pointCoord(3, 0) = 7.818346171771252e-03;
        pointCoord(3, 1) = 2.503974088635801e-01;
        pointCoord(3, 2) = 1.538058314339238e-03;
        pointCoord(4, 0) = 7.818346171771252e-03;
        pointCoord(4, 1) = 3.102450540391642e-01;
        pointCoord(4, 2) = 7.613571603892265e-04;
        pointCoord(5, 0) = 3.846089082451049e-02;
        pointCoord(5, 1) = 1.383248899231319e-02;
        pointCoord(5, 2) = 1.393272290854339e-03;
        pointCoord(6, 0) = 3.846089082451049e-02;
        pointCoord(6, 1) = 6.804634091093671e-02;
        pointCoord(6, 2) = 2.814623861941490e-03;
        pointCoord(7, 0) = 3.846089082451048e-02;
        pointCoord(7, 1) = 1.474362212543184e-01;
        pointCoord(7, 2) = 3.345408121477312e-03;
        pointCoord(8, 0) = 3.846089082451048e-02;
        pointCoord(8, 1) = 2.268261015977956e-01;
        pointCoord(8, 2) = 2.814623861948262e-03;
        pointCoord(9, 0) = 3.846089082451048e-02;
        pointCoord(9, 1) = 2.810399535165941e-01;
        pointCoord(9, 2) = 1.393272290859980e-03;
        pointCoord(10, 0) = 8.333333333333333e-02;
        pointCoord(10, 1) = 1.172751925766622e-02;
        pointCoord(10, 2) = 1.404011170703860e-03;
        pointCoord(11, 0) = 8.333333333333333e-02;
        pointCoord(11, 1) = 5.769133623677425e-02;
        pointCoord(11, 2) = 2.836318047403788e-03;
        pointCoord(12, 0) = 8.333333333333333e-02;
        pointCoord(12, 1) = 1.249999999999531e-01;
        pointCoord(12, 2) = 3.371193415638282e-03;
        pointCoord(13, 0) = 8.333333333333333e-02;
        pointCoord(13, 1) = 1.923086637631592e-01;
        pointCoord(13, 2) = 2.836318047406078e-03;
        pointCoord(14, 0) = 8.333333333333333e-02;
        pointCoord(14, 1) = 2.382724807423170e-01;
        pointCoord(14, 2) = 1.404011170705770e-03;
        pointCoord(15, 0) = 1.282057758421562e-01;
        pointCoord(15, 1) = 9.622549523018796e-03;
        pointCoord(15, 2) = 9.692277091454175e-04;
        pointCoord(16, 0) = 1.282057758421562e-01;
        pointCoord(16, 1) = 4.733633156260292e-02;
        pointCoord(16, 2) = 1.957988726055851e-03;
        pointCoord(17, 0) = 1.282057758421562e-01;
        pointCoord(17, 1) = 1.025637787455607e-01;
        pointCoord(17, 2) = 2.327227973327619e-03;
        pointCoord(18, 0) = 1.282057758421562e-01;
        pointCoord(18, 1) = 1.577912259284930e-01;
        pointCoord(18, 2) = 1.957988726054050e-03;
        pointCoord(19, 0) = 1.282057758421562e-01;
        pointCoord(19, 1) = 1.955050079680306e-01;
        pointCoord(19, 2) = 9.692277091439168e-04;
        pointCoord(20, 0) = 1.588483204948954e-01;
        pointCoord(20, 1) = 8.185105392941547e-03;
        pointCoord(20, 2) = 4.081084409135109e-04;
        pointCoord(21, 0) = 1.588483204948954e-01;
        pointCoord(21, 1) = 4.026509417575489e-02;
        pointCoord(21, 2) = 8.244416856610159e-04;
        pointCoord(22, 0) = 1.588483204948954e-01;
        pointCoord(22, 1) = 8.724250641921094e-02;
        pointCoord(22, 2) = 9.799156285814153e-04;
        pointCoord(23, 0) = 1.588483204948954e-01;
        pointCoord(23, 1) = 1.342199186626540e-01;
        pointCoord(23, 2) = 8.244416856605599e-04;
        pointCoord(24, 0) = 1.588483204948954e-01;
        pointCoord(24, 1) = 1.662999074454436e-01;
        pointCoord(24, 2) = 4.081084409131312e-04;
        pointCoord(25, 0) = 7.818346171771254e-03;
        pointCoord(25, 1) = 3.337000925545017e-01;
        pointCoord(25, 2) = 4.081084409108359e-04;
        pointCoord(26, 0) = 7.818346171771256e-03;
        pointCoord(26, 1) = 3.657800813372232e-01;
        pointCoord(26, 2) = 8.244416856612353e-04;
        pointCoord(27, 0) = 7.818346171771252e-03;
        pointCoord(27, 1) = 4.127574935808113e-01;
        pointCoord(27, 2) = 9.799156285854949e-04;
        pointCoord(28, 0) = 7.818346171771252e-03;
        pointCoord(28, 1) = 4.597349058243998e-01;
        pointCoord(28, 2) = 8.244416856612347e-04;
        pointCoord(29, 0) = 7.818346171771252e-03;
        pointCoord(29, 1) = 4.918148946071214e-01;
        pointCoord(29, 2) = 4.081084409108353e-04;
        pointCoord(30, 0) = 3.846089082451049e-02;
        pointCoord(30, 1) = 3.044949920319341e-01;
        pointCoord(30, 2) = 9.692277091409942e-04;
        pointCoord(31, 0) = 3.846089082451049e-02;
        pointCoord(31, 1) = 3.422087740714365e-01;
        pointCoord(31, 2) = 1.957988726055611e-03;
        pointCoord(32, 0) = 3.846089082451048e-02;
        pointCoord(32, 1) = 3.974362212544831e-01;
        pointCoord(32, 2) = 2.327227973333649e-03;
        pointCoord(33, 0) = 3.846089082451048e-02;
        pointCoord(33, 1) = 4.526636684375298e-01;
        pointCoord(33, 2) = 1.957988726055609e-03;
        pointCoord(34, 0) = 3.846089082451048e-02;
        pointCoord(34, 1) = 4.903774504770321e-01;
        pointCoord(34, 2) = 9.692277091409940e-04;
        pointCoord(35, 0) = 8.333333333333333e-02;
        pointCoord(35, 1) = 2.617275192576670e-01;
        pointCoord(35, 2) = 1.404011170704044e-03;
        pointCoord(36, 0) = 8.333333333333333e-02;
        pointCoord(36, 1) = 3.076913362367896e-01;
        pointCoord(36, 2) = 2.836318047405072e-03;
        pointCoord(37, 0) = 8.333333333333333e-02;
        pointCoord(37, 1) = 3.750000000000000e-01;
        pointCoord(37, 2) = 3.371193415639545e-03;
        pointCoord(38, 0) = 8.333333333333333e-02;
        pointCoord(38, 1) = 4.423086637632104e-01;
        pointCoord(38, 2) = 2.836318047405072e-03;
        pointCoord(39, 0) = 8.333333333333333e-02;
        pointCoord(39, 1) = 4.882724807423330e-01;
        pointCoord(39, 2) = 1.404011170704044e-03;
        pointCoord(40, 0) = 1.282057758421562e-01;
        pointCoord(40, 1) = 2.189600464833999e-01;
        pointCoord(40, 2) = 1.393272290859161e-03;
        pointCoord(41, 0) = 1.282057758421562e-01;
        pointCoord(41, 1) = 2.731738984021428e-01;
        pointCoord(41, 2) = 2.814623861944518e-03;
        pointCoord(42, 0) = 1.282057758421562e-01;
        pointCoord(42, 1) = 3.525637787455169e-01;
        pointCoord(42, 2) = 3.345408121474027e-03;
        pointCoord(43, 0) = 1.282057758421562e-01;
        pointCoord(43, 1) = 4.319536590888911e-01;
        pointCoord(43, 2) = 2.814623861944516e-03;
        pointCoord(44, 0) = 1.282057758421562e-01;
        pointCoord(44, 1) = 4.861675110076338e-01;
        pointCoord(44, 2) = 1.393272290859160e-03;
        pointCoord(45, 0) = 1.588483204948954e-01;
        pointCoord(45, 1) = 1.897549459608325e-01;
        pointCoord(45, 2) = 7.613571603889063e-04;
        pointCoord(46, 0) = 1.588483204948954e-01;
        pointCoord(46, 1) = 2.496025911363562e-01;
        pointCoord(46, 2) = 1.538058314337037e-03;
        pointCoord(47, 0) = 1.588483204948954e-01;
        pointCoord(47, 1) = 3.372425064191886e-01;
        pointCoord(47, 2) = 1.828106712819135e-03;
        pointCoord(48, 0) = 1.588483204948954e-01;
        pointCoord(48, 1) = 4.248824217020209e-01;
        pointCoord(48, 2) = 1.538058314337035e-03;
        pointCoord(49, 0) = 1.588483204948954e-01;
        pointCoord(49, 1) = 4.847300668775447e-01;
        pointCoord(49, 2) = 7.613571603889057e-04;
        pointCoord(50, 0) = 7.818346171771254e-03;
        pointCoord(50, 1) = 5.152699331224554e-01;
        pointCoord(50, 2) = 7.613571603889057e-04;
        pointCoord(51, 0) = 7.818346171771256e-03;
        pointCoord(51, 1) = 5.751175782979793e-01;
        pointCoord(51, 2) = 1.538058314337038e-03;
        pointCoord(52, 0) = 7.818346171771252e-03;
        pointCoord(52, 1) = 6.627574935808115e-01;
        pointCoord(52, 2) = 1.828106712819137e-03;
        pointCoord(53, 0) = 7.818346171771252e-03;
        pointCoord(53, 1) = 7.503974088636437e-01;
        pointCoord(53, 2) = 1.538058314337037e-03;
        pointCoord(54, 0) = 7.818346171771252e-03;
        pointCoord(54, 1) = 8.102450540391676e-01;
        pointCoord(54, 2) = 7.613571603889051e-04;
        pointCoord(55, 0) = 3.846089082451049e-02;
        pointCoord(55, 1) = 5.138324889923662e-01;
        pointCoord(55, 2) = 1.393272290859161e-03;
        pointCoord(56, 0) = 3.846089082451049e-02;
        pointCoord(56, 1) = 5.680463409111093e-01;
        pointCoord(56, 2) = 2.814623861944515e-03;
        pointCoord(57, 0) = 3.846089082451048e-02;
        pointCoord(57, 1) = 6.474362212544831e-01;
        pointCoord(57, 2) = 3.345408121474031e-03;
        pointCoord(58, 0) = 3.846089082451048e-02;
        pointCoord(58, 1) = 7.268261015978573e-01;
        pointCoord(58, 2) = 2.814623861944519e-03;
        pointCoord(59, 0) = 3.846089082451048e-02;
        pointCoord(59, 1) = 7.810399535166002e-01;
        pointCoord(59, 2) = 1.393272290859163e-03;
        pointCoord(60, 0) = 8.333333333333333e-02;
        pointCoord(60, 1) = 5.117275192576670e-01;
        pointCoord(60, 2) = 1.404011170704043e-03;
        pointCoord(61, 0) = 8.333333333333333e-02;
        pointCoord(61, 1) = 5.576913362367897e-01;
        pointCoord(61, 2) = 2.836318047405074e-03;
        pointCoord(62, 0) = 8.333333333333333e-02;
        pointCoord(62, 1) = 6.250000000000000e-01;
        pointCoord(62, 2) = 3.371193415639546e-03;
        pointCoord(63, 0) = 8.333333333333333e-02;
        pointCoord(63, 1) = 6.923086637632103e-01;
        pointCoord(63, 2) = 2.836318047405073e-03;
        pointCoord(64, 0) = 8.333333333333333e-02;
        pointCoord(64, 1) = 7.382724807423330e-01;
        pointCoord(64, 2) = 1.404011170704046e-03;
        pointCoord(65, 0) = 1.282057758421562e-01;
        pointCoord(65, 1) = 5.096225495229679e-01;
        pointCoord(65, 2) = 9.692277091409942e-04;
        pointCoord(66, 0) = 1.282057758421562e-01;
        pointCoord(66, 1) = 5.473363315624705e-01;
        pointCoord(66, 2) = 1.957988726055611e-03;
        pointCoord(67, 0) = 1.282057758421562e-01;
        pointCoord(67, 1) = 6.025637787455169e-01;
        pointCoord(67, 2) = 2.327227973333649e-03;
        pointCoord(68, 0) = 1.282057758421562e-01;
        pointCoord(68, 1) = 6.577912259285635e-01;
        pointCoord(68, 2) = 1.957988726055610e-03;
        pointCoord(69, 0) = 1.282057758421562e-01;
        pointCoord(69, 1) = 6.955050079680659e-01;
        pointCoord(69, 2) = 9.692277091409909e-04;
        pointCoord(70, 0) = 1.588483204948954e-01;
        pointCoord(70, 1) = 5.081851053928789e-01;
        pointCoord(70, 2) = 4.081084409108353e-04;
        pointCoord(71, 0) = 1.588483204948954e-01;
        pointCoord(71, 1) = 5.402650941756003e-01;
        pointCoord(71, 2) = 8.244416856612347e-04;
        pointCoord(72, 0) = 1.588483204948954e-01;
        pointCoord(72, 1) = 5.872425064191885e-01;
        pointCoord(72, 2) = 9.799156285854923e-04;
        pointCoord(73, 0) = 1.588483204948954e-01;
        pointCoord(73, 1) = 6.342199186627768e-01;
        pointCoord(73, 2) = 8.244416856612338e-04;
        pointCoord(74, 0) = 1.588483204948954e-01;
        pointCoord(74, 1) = 6.662999074454984e-01;
        pointCoord(74, 2) = 4.081084409108362e-04;
        pointCoord(75, 0) = 8.185105392941545e-03;
        pointCoord(75, 1) = 8.329665741121656e-01;
        pointCoord(75, 2) = 4.081084409107431e-04;
        pointCoord(76, 0) = 9.622549523018798e-03;
        pointCoord(76, 1) = 8.621716746347330e-01;
        pointCoord(76, 2) = 9.692277091449929e-04;
        pointCoord(77, 0) = 1.172751925766622e-02;
        pointCoord(77, 1) = 9.049391474089994e-01;
        pointCoord(77, 2) = 1.404011170708332e-03;
        pointCoord(78, 0) = 1.383248899231319e-02;
        pointCoord(78, 1) = 9.477066201832667e-01;
        pointCoord(78, 2) = 1.393272290854221e-03;
        pointCoord(79, 0) = 1.526993312238960e-02;
        pointCoord(79, 1) = 9.769117207058344e-01;
        pointCoord(79, 2) = 7.613571603811338e-04;
        pointCoord(80, 0) = 4.026509417575489e-02;
        pointCoord(80, 1) = 8.008865853294439e-01;
        pointCoord(80, 2) = 8.244416856564127e-04;
        pointCoord(81, 0) = 4.733633156260292e-02;
        pointCoord(81, 1) = 8.244578925952301e-01;
        pointCoord(81, 2) = 1.957988726052771e-03;
        pointCoord(82, 0) = 5.769133623677425e-02;
        pointCoord(82, 1) = 8.589753304298769e-01;
        pointCoord(82, 2) = 2.836318047405692e-03;
        pointCoord(83, 0) = 6.804634091093668e-02;
        pointCoord(83, 1) = 8.934927682645241e-01;
        pointCoord(83, 2) = 2.814623861936720e-03;
        pointCoord(84, 0) = 7.511757829776837e-02;
        pointCoord(84, 1) = 9.170640755303107e-01;
        pointCoord(84, 2) = 1.538058314326002e-03;
        pointCoord(85, 0) = 8.724250641921094e-02;
        pointCoord(85, 1) = 7.539091730858550e-01;
        pointCoord(85, 2) = 9.799156285793162e-04;
        pointCoord(86, 0) = 1.025637787455607e-01;
        pointCoord(86, 1) = 7.692304454121838e-01;
        pointCoord(86, 2) = 2.327227973330276e-03;
        pointCoord(87, 0) = 1.249999999999531e-01;
        pointCoord(87, 1) = 7.916666666666666e-01;
        pointCoord(87, 2) = 3.371193415644602e-03;
        pointCoord(88, 0) = 1.474362212543184e-01;
        pointCoord(88, 1) = 8.141028879211499e-01;
        pointCoord(88, 2) = 3.345408121471103e-03;
        pointCoord(89, 0) = 1.627574935806183e-01;
        pointCoord(89, 1) = 8.294241602474780e-01;
        pointCoord(89, 2) = 1.828106712808768e-03;
        pointCoord(90, 0) = 1.342199186626541e-01;
        pointCoord(90, 1) = 7.069317608422671e-01;
        pointCoord(90, 2) = 8.244416856576413e-04;
        pointCoord(91, 0) = 1.577912259284931e-01;
        pointCoord(91, 1) = 7.140029982291373e-01;
        pointCoord(91, 2) = 1.957988726058514e-03;
        pointCoord(92, 0) = 1.923086637631592e-01;
        pointCoord(92, 1) = 7.243580029034563e-01;
        pointCoord(92, 2) = 2.836318047417890e-03;
        pointCoord(93, 0) = 2.268261015977956e-01;
        pointCoord(93, 1) = 7.347130075777756e-01;
        pointCoord(93, 2) = 2.814623861946352e-03;
        pointCoord(94, 0) = 2.503974088635803e-01;
        pointCoord(94, 1) = 7.417842449646458e-01;
        pointCoord(94, 2) = 1.538058314327700e-03;
        pointCoord(95, 0) = 1.662999074454437e-01;
        pointCoord(95, 1) = 6.748517720595456e-01;
        pointCoord(95, 2) = 4.081084409085538e-04;
        pointCoord(96, 0) = 1.955050079680307e-01;
        pointCoord(96, 1) = 6.762892161896349e-01;
        pointCoord(96, 2) = 9.692277091422312e-04;
        pointCoord(97, 0) = 2.382724807423170e-01;
        pointCoord(97, 1) = 6.783941859243336e-01;
        pointCoord(97, 2) = 1.404011170710890e-03;
        pointCoord(98, 0) = 2.810399535165941e-01;
        pointCoord(98, 1) = 6.804991556590326e-01;
        pointCoord(98, 2) = 1.393272290860203e-03;
        pointCoord(99, 0) = 3.102450540391642e-01;
        pointCoord(99, 1) = 6.819365997891221e-01;
        pointCoord(99, 2) = 7.613571603839364e-04;
        pointCoord(100, 0) = 1.819365997891220e-01;
        pointCoord(100, 1) = 7.818346171771254e-03;
        pointCoord(100, 2) = 7.613571603889064e-04;
        pointCoord(101, 0) = 1.804991556590328e-01;
        pointCoord(101, 1) = 3.846089082451050e-02;
        pointCoord(101, 2) = 1.393272290859161e-03;
        pointCoord(102, 0) = 1.783941859243336e-01;
        pointCoord(102, 1) = 8.333333333333333e-02;
        pointCoord(102, 2) = 1.404011170704044e-03;
        pointCoord(103, 0) = 1.762892161896346e-01;
        pointCoord(103, 1) = 1.282057758421562e-01;
        pointCoord(103, 2) = 9.692277091409926e-04;
        pointCoord(104, 0) = 1.748517720595454e-01;
        pointCoord(104, 1) = 1.588483204948954e-01;
        pointCoord(104, 2) = 4.081084409108355e-04;
        pointCoord(105, 0) = 2.417842449646457e-01;
        pointCoord(105, 1) = 7.818346171771254e-03;
        pointCoord(105, 2) = 1.538058314337037e-03;
        pointCoord(106, 0) = 2.347130075777757e-01;
        pointCoord(106, 1) = 3.846089082451048e-02;
        pointCoord(106, 2) = 2.814623861944517e-03;
        pointCoord(107, 0) = 2.243580029034563e-01;
        pointCoord(107, 1) = 8.333333333333333e-02;
        pointCoord(107, 2) = 2.836318047405072e-03;
        pointCoord(108, 0) = 2.140029982291369e-01;
        pointCoord(108, 1) = 1.282057758421562e-01;
        pointCoord(108, 2) = 1.957988726055610e-03;
        pointCoord(109, 0) = 2.069317608422668e-01;
        pointCoord(109, 1) = 1.588483204948954e-01;
        pointCoord(109, 2) = 8.244416856612345e-04;
        pointCoord(110, 0) = 3.294241602474781e-01;
        pointCoord(110, 1) = 7.818346171771249e-03;
        pointCoord(110, 2) = 1.828106712819136e-03;
        pointCoord(111, 0) = 3.141028879211498e-01;
        pointCoord(111, 1) = 3.846089082451047e-02;
        pointCoord(111, 2) = 3.345408121474029e-03;
        pointCoord(112, 0) = 2.916666666666667e-01;
        pointCoord(112, 1) = 8.333333333333333e-02;
        pointCoord(112, 2) = 3.371193415639546e-03;
        pointCoord(113, 0) = 2.692304454121835e-01;
        pointCoord(113, 1) = 1.282057758421562e-01;
        pointCoord(113, 2) = 2.327227973333647e-03;
        pointCoord(114, 0) = 2.539091730858553e-01;
        pointCoord(114, 1) = 1.588483204948954e-01;
        pointCoord(114, 2) = 9.799156285854936e-04;
        pointCoord(115, 0) = 4.170640755303106e-01;
        pointCoord(115, 1) = 7.818346171771254e-03;
        pointCoord(115, 2) = 1.538058314337037e-03;
        pointCoord(116, 0) = 3.934927682645241e-01;
        pointCoord(116, 1) = 3.846089082451049e-02;
        pointCoord(116, 2) = 2.814623861944518e-03;
        pointCoord(117, 0) = 3.589753304298771e-01;
        pointCoord(117, 1) = 8.333333333333333e-02;
        pointCoord(117, 2) = 2.836318047405072e-03;
        pointCoord(118, 0) = 3.244578925952302e-01;
        pointCoord(118, 1) = 1.282057758421562e-01;
        pointCoord(118, 2) = 1.957988726055611e-03;
        pointCoord(119, 0) = 3.008865853294436e-01;
        pointCoord(119, 1) = 1.588483204948954e-01;
        pointCoord(119, 2) = 8.244416856612344e-04;
        pointCoord(120, 0) = 4.769117207058344e-01;
        pointCoord(120, 1) = 7.818346171771256e-03;
        pointCoord(120, 2) = 7.613571603889065e-04;
        pointCoord(121, 0) = 4.477066201832667e-01;
        pointCoord(121, 1) = 3.846089082451048e-02;
        pointCoord(121, 2) = 1.393272290859160e-03;
        pointCoord(122, 0) = 4.049391474089996e-01;
        pointCoord(122, 1) = 8.333333333333333e-02;
        pointCoord(122, 2) = 1.404011170704044e-03;
        pointCoord(123, 0) = 3.621716746347325e-01;
        pointCoord(123, 1) = 1.282057758421562e-01;
        pointCoord(123, 2) = 9.692277091409923e-04;
        pointCoord(124, 0) = 3.329665741121651e-01;
        pointCoord(124, 1) = 1.588483204948954e-01;
        pointCoord(124, 2) = 4.081084409108357e-04;
        pointCoord(125, 0) = 1.744850128383772e-01;
        pointCoord(125, 1) = 1.819365997891220e-01;
        pointCoord(125, 2) = 7.613571603838925e-04;
        pointCoord(126, 0) = 1.744850128383772e-01;
        pointCoord(126, 1) = 2.417842449646458e-01;
        pointCoord(126, 2) = 1.538058314326908e-03;
        pointCoord(127, 0) = 1.744850128383771e-01;
        pointCoord(127, 1) = 3.294241602474781e-01;
        pointCoord(127, 2) = 1.828106712807097e-03;
        pointCoord(128, 0) = 1.744850128383772e-01;
        pointCoord(128, 1) = 4.170640755303105e-01;
        pointCoord(128, 2) = 1.538058314326907e-03;
        pointCoord(129, 0) = 1.744850128383772e-01;
        pointCoord(129, 1) = 4.769117207058342e-01;
        pointCoord(129, 2) = 7.613571603838910e-04;
        pointCoord(130, 0) = 2.051275574910337e-01;
        pointCoord(130, 1) = 1.804991556590328e-01;
        pointCoord(130, 2) = 1.393272290859977e-03;
        pointCoord(131, 0) = 2.051275574910338e-01;
        pointCoord(131, 1) = 2.347130075777757e-01;
        pointCoord(131, 2) = 2.814623861946165e-03;
        pointCoord(132, 0) = 2.051275574910338e-01;
        pointCoord(132, 1) = 3.141028879211498e-01;
        pointCoord(132, 2) = 3.345408121475988e-03;
        pointCoord(133, 0) = 2.051275574910338e-01;
        pointCoord(133, 1) = 3.934927682645240e-01;
        pointCoord(133, 2) = 2.814623861946165e-03;
        pointCoord(134, 0) = 2.051275574910337e-01;
        pointCoord(134, 1) = 4.477066201832667e-01;
        pointCoord(134, 2) = 1.393272290859977e-03;
        pointCoord(135, 0) = 2.500000000000000e-01;
        pointCoord(135, 1) = 1.783941859243336e-01;
        pointCoord(135, 2) = 1.404011170710362e-03;
        pointCoord(136, 0) = 2.500000000000000e-01;
        pointCoord(136, 1) = 2.243580029034562e-01;
        pointCoord(136, 2) = 2.836318047417836e-03;
        pointCoord(137, 0) = 2.500000000000000e-01;
        pointCoord(137, 1) = 2.916666666666667e-01;
        pointCoord(137, 2) = 3.371193415654717e-03;
        pointCoord(138, 0) = 2.500000000000000e-01;
        pointCoord(138, 1) = 3.589753304298771e-01;
        pointCoord(138, 2) = 2.836318047417834e-03;
        pointCoord(139, 0) = 2.500000000000000e-01;
        pointCoord(139, 1) = 4.049391474089997e-01;
        pointCoord(139, 2) = 1.404011170710363e-03;
        pointCoord(140, 0) = 2.948724425089663e-01;
        pointCoord(140, 1) = 1.762892161896346e-01;
        pointCoord(140, 2) = 9.692277091415605e-04;
        pointCoord(141, 0) = 2.948724425089662e-01;
        pointCoord(141, 1) = 2.140029982291370e-01;
        pointCoord(141, 2) = 1.957988726056759e-03;
        pointCoord(142, 0) = 2.948724425089662e-01;
        pointCoord(142, 1) = 2.692304454121835e-01;
        pointCoord(142, 2) = 2.327227973335011e-03;
        pointCoord(143, 0) = 2.948724425089662e-01;
        pointCoord(143, 1) = 3.244578925952302e-01;
        pointCoord(143, 2) = 1.957988726056759e-03;
        pointCoord(144, 0) = 2.948724425089662e-01;
        pointCoord(144, 1) = 3.621716746347325e-01;
        pointCoord(144, 2) = 9.692277091415591e-04;
        pointCoord(145, 0) = 3.255149871616230e-01;
        pointCoord(145, 1) = 1.748517720595454e-01;
        pointCoord(145, 2) = 4.081084409081483e-04;
        pointCoord(146, 0) = 3.255149871616229e-01;
        pointCoord(146, 1) = 2.069317608422669e-01;
        pointCoord(146, 2) = 8.244416856558047e-04;
        pointCoord(147, 0) = 3.255149871616228e-01;
        pointCoord(147, 1) = 2.539091730858553e-01;
        pointCoord(147, 2) = 9.799156285790393e-04;
        pointCoord(148, 0) = 3.255149871616228e-01;
        pointCoord(148, 1) = 3.008865853294436e-01;
        pointCoord(148, 2) = 8.244416856558036e-04;
        pointCoord(149, 0) = 3.255149871616229e-01;
        pointCoord(149, 1) = 3.329665741121651e-01;
        pointCoord(149, 2) = 4.081084409081480e-04;
        pointCoord(150, 0) = 1.748517720595455e-01;
        pointCoord(150, 1) = 4.996332407788319e-01;
        pointCoord(150, 2) = 4.081084409081488e-04;
        pointCoord(151, 0) = 1.762892161896346e-01;
        pointCoord(151, 1) = 5.288383413013994e-01;
        pointCoord(151, 2) = 9.692277091415623e-04;
        pointCoord(152, 0) = 1.783941859243336e-01;
        pointCoord(152, 1) = 5.716058140756662e-01;
        pointCoord(152, 2) = 1.404011170710364e-03;
        pointCoord(153, 0) = 1.804991556590328e-01;
        pointCoord(153, 1) = 6.143732868499333e-01;
        pointCoord(153, 2) = 1.393272290859978e-03;
        pointCoord(154, 0) = 1.819365997891220e-01;
        pointCoord(154, 1) = 6.435783873725009e-01;
        pointCoord(154, 2) = 7.613571603838906e-04;
        pointCoord(155, 0) = 2.069317608422669e-01;
        pointCoord(155, 1) = 4.675532519961104e-01;
        pointCoord(155, 2) = 8.244416856558047e-04;
        pointCoord(156, 0) = 2.140029982291370e-01;
        pointCoord(156, 1) = 4.911245592618969e-01;
        pointCoord(156, 2) = 1.957988726056756e-03;
        pointCoord(157, 0) = 2.243580029034563e-01;
        pointCoord(157, 1) = 5.256419970965437e-01;
        pointCoord(157, 2) = 2.836318047417835e-03;
        pointCoord(158, 0) = 2.347130075777757e-01;
        pointCoord(158, 1) = 5.601594349311906e-01;
        pointCoord(158, 2) = 2.814623861946166e-03;
        pointCoord(159, 0) = 2.417842449646457e-01;
        pointCoord(159, 1) = 5.837307421969772e-01;
        pointCoord(159, 2) = 1.538058314326912e-03;
        pointCoord(160, 0) = 2.539091730858552e-01;
        pointCoord(160, 1) = 4.205758397525219e-01;
        pointCoord(160, 2) = 9.799156285790406e-04;
        pointCoord(161, 0) = 2.692304454121835e-01;
        pointCoord(161, 1) = 4.358971120788502e-01;
        pointCoord(161, 2) = 2.327227973335012e-03;
        pointCoord(162, 0) = 2.916666666666667e-01;
        pointCoord(162, 1) = 4.583333333333333e-01;
        pointCoord(162, 2) = 3.371193415654718e-03;
        pointCoord(163, 0) = 3.141028879211498e-01;
        pointCoord(163, 1) = 4.807695545878165e-01;
        pointCoord(163, 2) = 3.345408121475987e-03;
        pointCoord(164, 0) = 3.294241602474781e-01;
        pointCoord(164, 1) = 4.960908269141448e-01;
        pointCoord(164, 2) = 1.828106712807097e-03;
        pointCoord(165, 0) = 3.008865853294437e-01;
        pointCoord(165, 1) = 3.735984275089335e-01;
        pointCoord(165, 2) = 8.244416856558042e-04;
        pointCoord(166, 0) = 3.244578925952302e-01;
        pointCoord(166, 1) = 3.806696648958037e-01;
        pointCoord(166, 2) = 1.957988726056757e-03;
        pointCoord(167, 0) = 3.589753304298770e-01;
        pointCoord(167, 1) = 3.910246695701230e-01;
        pointCoord(167, 2) = 2.836318047417832e-03;
        pointCoord(168, 0) = 3.934927682645240e-01;
        pointCoord(168, 1) = 4.013796742444424e-01;
        pointCoord(168, 2) = 2.814623861946166e-03;
        pointCoord(169, 0) = 4.170640755303105e-01;
        pointCoord(169, 1) = 4.084509116313124e-01;
        pointCoord(169, 2) = 1.538058314326906e-03;
        pointCoord(170, 0) = 3.329665741121652e-01;
        pointCoord(170, 1) = 3.415184387262122e-01;
        pointCoord(170, 2) = 4.081084409081477e-04;
        pointCoord(171, 0) = 3.621716746347325e-01;
        pointCoord(171, 1) = 3.429558828563013e-01;
        pointCoord(171, 2) = 9.692277091415618e-04;
        pointCoord(172, 0) = 4.049391474089996e-01;
        pointCoord(172, 1) = 3.450608525910003e-01;
        pointCoord(172, 2) = 1.404011170710361e-03;
        pointCoord(173, 0) = 4.477066201832667e-01;
        pointCoord(173, 1) = 3.471658223256994e-01;
        pointCoord(173, 2) = 1.393272290859977e-03;
        pointCoord(174, 0) = 4.769117207058342e-01;
        pointCoord(174, 1) = 3.486032664557886e-01;
        pointCoord(174, 2) = 7.613571603838919e-04;
        pointCoord(175, 0) = 1.897549459608325e-01;
        pointCoord(175, 1) = 6.513967335442116e-01;
        pointCoord(175, 2) = 7.613571603838933e-04;
        pointCoord(176, 0) = 2.189600464833999e-01;
        pointCoord(176, 1) = 6.528341776743006e-01;
        pointCoord(176, 2) = 1.393272290859980e-03;
        pointCoord(177, 0) = 2.617275192576670e-01;
        pointCoord(177, 1) = 6.549391474089995e-01;
        pointCoord(177, 2) = 1.404011170710364e-03;
        pointCoord(178, 0) = 3.044949920319341e-01;
        pointCoord(178, 1) = 6.570441171436988e-01;
        pointCoord(178, 2) = 9.692277091415618e-04;
        pointCoord(179, 0) = 3.337000925545016e-01;
        pointCoord(179, 1) = 6.584815612737880e-01;
        pointCoord(179, 2) = 4.081084409081461e-04;
        pointCoord(180, 0) = 2.496025911363562e-01;
        pointCoord(180, 1) = 5.915490883686877e-01;
        pointCoord(180, 2) = 1.538058314326907e-03;
        pointCoord(181, 0) = 2.731738984021428e-01;
        pointCoord(181, 1) = 5.986203257555579e-01;
        pointCoord(181, 2) = 2.814623861946164e-03;
        pointCoord(182, 0) = 3.076913362367896e-01;
        pointCoord(182, 1) = 6.089753304298771e-01;
        pointCoord(182, 2) = 2.836318047417837e-03;
        pointCoord(183, 0) = 3.422087740714365e-01;
        pointCoord(183, 1) = 6.193303351041964e-01;
        pointCoord(183, 2) = 1.957988726056759e-03;
        pointCoord(184, 0) = 3.657800813372230e-01;
        pointCoord(184, 1) = 6.264015724910664e-01;
        pointCoord(184, 2) = 8.244416856558064e-04;
        pointCoord(185, 0) = 3.372425064191886e-01;
        pointCoord(185, 1) = 5.039091730858553e-01;
        pointCoord(185, 2) = 1.828106712807096e-03;
        pointCoord(186, 0) = 3.525637787455169e-01;
        pointCoord(186, 1) = 5.192304454121835e-01;
        pointCoord(186, 2) = 3.345408121475984e-03;
        pointCoord(187, 0) = 3.750000000000000e-01;
        pointCoord(187, 1) = 5.416666666666666e-01;
        pointCoord(187, 2) = 3.371193415654719e-03;
        pointCoord(188, 0) = 3.974362212544831e-01;
        pointCoord(188, 1) = 5.641028879211498e-01;
        pointCoord(188, 2) = 2.327227973335008e-03;
        pointCoord(189, 0) = 4.127574935808115e-01;
        pointCoord(189, 1) = 5.794241602474780e-01;
        pointCoord(189, 2) = 9.799156285790400e-04;
        pointCoord(190, 0) = 4.248824217020210e-01;
        pointCoord(190, 1) = 4.162692578030228e-01;
        pointCoord(190, 2) = 1.538058314326908e-03;
        pointCoord(191, 0) = 4.319536590888911e-01;
        pointCoord(191, 1) = 4.398405650688095e-01;
        pointCoord(191, 2) = 2.814623861946166e-03;
        pointCoord(192, 0) = 4.423086637632105e-01;
        pointCoord(192, 1) = 4.743580029034563e-01;
        pointCoord(192, 2) = 2.836318047417831e-03;
        pointCoord(193, 0) = 4.526636684375298e-01;
        pointCoord(193, 1) = 5.088754407381032e-01;
        pointCoord(193, 2) = 1.957988726056758e-03;
        pointCoord(194, 0) = 4.597349058243999e-01;
        pointCoord(194, 1) = 5.324467480038898e-01;
        pointCoord(194, 2) = 8.244416856558055e-04;
        pointCoord(195, 0) = 4.847300668775448e-01;
        pointCoord(195, 1) = 3.564216126274991e-01;
        pointCoord(195, 2) = 7.613571603838913e-04;
        pointCoord(196, 0) = 4.861675110076339e-01;
        pointCoord(196, 1) = 3.856267131500666e-01;
        pointCoord(196, 2) = 1.393272290859975e-03;
        pointCoord(197, 0) = 4.882724807423328e-01;
        pointCoord(197, 1) = 4.283941859243336e-01;
        pointCoord(197, 2) = 1.404011170710359e-03;
        pointCoord(198, 0) = 4.903774504770320e-01;
        pointCoord(198, 1) = 4.711616586986007e-01;
        pointCoord(198, 2) = 9.692277091415589e-04;
        pointCoord(199, 0) = 4.918148946071213e-01;
        pointCoord(199, 1) = 5.003667592211682e-01;
        pointCoord(199, 2) = 4.081084409081484e-04;
        pointCoord(200, 0) = 3.411516795050438e-01;
        pointCoord(200, 1) = 1.670334258877822e-01;
        pointCoord(200, 2) = 4.081084409085300e-04;
        pointCoord(201, 0) = 3.411516795050439e-01;
        pointCoord(201, 1) = 1.991134146705283e-01;
        pointCoord(201, 2) = 8.244416856563124e-04;
        pointCoord(202, 0) = 3.411516795050438e-01;
        pointCoord(202, 1) = 2.460908269141367e-01;
        pointCoord(202, 2) = 9.799156285793017e-04;
        pointCoord(203, 0) = 3.411516795050438e-01;
        pointCoord(203, 1) = 2.930682391577322e-01;
        pointCoord(203, 2) = 8.244416856558556e-04;
        pointCoord(204, 0) = 3.411516795050438e-01;
        pointCoord(204, 1) = 3.251482279404546e-01;
        pointCoord(204, 2) = 4.081084409081496e-04;
        pointCoord(205, 0) = 3.717942241577004e-01;
        pointCoord(205, 1) = 1.378283253651400e-01;
        pointCoord(205, 2) = 9.692277091430978e-04;
        pointCoord(206, 0) = 3.717942241577006e-01;
        pointCoord(206, 1) = 1.755421074046936e-01;
        pointCoord(206, 2) = 1.957988726059005e-03;
        pointCoord(207, 0) = 3.717942241577005e-01;
        pointCoord(207, 1) = 2.307695545877885e-01;
        pointCoord(207, 2) = 2.327227973336458e-03;
        pointCoord(208, 0) = 3.717942241577005e-01;
        pointCoord(208, 1) = 2.859970017708580e-01;
        pointCoord(208, 2) = 1.957988726057206e-03;
        pointCoord(209, 0) = 3.717942241577004e-01;
        pointCoord(209, 1) = 3.237107838103652e-01;
        pointCoord(209, 2) = 9.692277091415980e-04;
        pointCoord(210, 0) = 4.166666666666667e-01;
        pointCoord(210, 1) = 9.506085259098436e-02;
        pointCoord(210, 2) = 1.404011170708637e-03;
        pointCoord(211, 0) = 4.166666666666666e-01;
        pointCoord(211, 1) = 1.410246695700717e-01;
        pointCoord(211, 2) = 2.836318047416827e-03;
        pointCoord(212, 0) = 4.166666666666667e-01;
        pointCoord(212, 1) = 2.083333333332865e-01;
        pointCoord(212, 2) = 3.371193415655979e-03;
        pointCoord(213, 0) = 4.166666666666667e-01;
        pointCoord(213, 1) = 2.756419970965283e-01;
        pointCoord(213, 2) = 2.836318047419117e-03;
        pointCoord(214, 0) = 4.166666666666667e-01;
        pointCoord(214, 1) = 3.216058140756655e-01;
        pointCoord(214, 2) = 1.404011170710547e-03;
        pointCoord(215, 0) = 4.615391091756330e-01;
        pointCoord(215, 1) = 5.229337981681954e-02;
        pointCoord(215, 2) = 1.393272290854696e-03;
        pointCoord(216, 0) = 4.615391091756330e-01;
        pointCoord(216, 1) = 1.065072317354202e-01;
        pointCoord(216, 2) = 2.814623861941732e-03;
        pointCoord(217, 0) = 4.615391091756330e-01;
        pointCoord(217, 1) = 1.858971120787572e-01;
        pointCoord(217, 2) = 3.345408121477283e-03;
        pointCoord(218, 0) = 4.615391091756329e-01;
        pointCoord(218, 1) = 2.652869924221897e-01;
        pointCoord(218, 2) = 2.814623861948503e-03;
        pointCoord(219, 0) = 4.615391091756330e-01;
        pointCoord(219, 1) = 3.195008443409654e-01;
        pointCoord(219, 2) = 1.393272290860338e-03;
        pointCoord(220, 0) = 4.921816538282898e-01;
        pointCoord(220, 1) = 2.308827929416035e-02;
        pointCoord(220, 2) = 7.613571603808961e-04;
        pointCoord(221, 0) = 4.921816538282895e-01;
        pointCoord(221, 1) = 8.293592446953112e-02;
        pointCoord(221, 2) = 1.538058314324874e-03;
        pointCoord(222, 0) = 4.921816538282893e-01;
        pointCoord(222, 1) = 1.705758397523592e-01;
        pointCoord(222, 2) = 1.828106712808753e-03;
        pointCoord(223, 0) = 4.921816538282894e-01;
        pointCoord(223, 1) = 2.582157550352992e-01;
        pointCoord(223, 2) = 1.538058314328803e-03;
        pointCoord(224, 0) = 4.921816538282895e-01;
        pointCoord(224, 1) = 3.180634002108751e-01;
        pointCoord(224, 2) = 7.613571603841691e-04;
        pointCoord(225, 0) = 5.078183461717106e-01;
        pointCoord(225, 1) = 1.526993312238961e-02;
        pointCoord(225, 2) = 7.613571603809396e-04;
        pointCoord(226, 0) = 5.078183461717106e-01;
        pointCoord(226, 1) = 7.511757829776841e-02;
        pointCoord(226, 2) = 1.538058314325179e-03;
        pointCoord(227, 0) = 5.078183461717104e-01;
        pointCoord(227, 1) = 1.627574935806183e-01;
        pointCoord(227, 2) = 1.828106712809257e-03;
        pointCoord(228, 0) = 5.078183461717104e-01;
        pointCoord(228, 1) = 2.503974088635801e-01;
        pointCoord(228, 2) = 1.538058314329110e-03;
        pointCoord(229, 0) = 5.078183461717106e-01;
        pointCoord(229, 1) = 3.102450540391642e-01;
        pointCoord(229, 2) = 7.613571603842124e-04;
        pointCoord(230, 0) = 5.384608908243670e-01;
        pointCoord(230, 1) = 1.383248899231319e-02;
        pointCoord(230, 2) = 1.393272290855154e-03;
        pointCoord(231, 0) = 5.384608908243671e-01;
        pointCoord(231, 1) = 6.804634091093671e-02;
        pointCoord(231, 2) = 2.814623861943139e-03;
        pointCoord(232, 0) = 5.384608908243671e-01;
        pointCoord(232, 1) = 1.474362212543184e-01;
        pointCoord(232, 2) = 3.345408121479272e-03;
        pointCoord(233, 0) = 5.384608908243672e-01;
        pointCoord(233, 1) = 2.268261015977956e-01;
        pointCoord(233, 2) = 2.814623861949911e-03;
        pointCoord(234, 0) = 5.384608908243671e-01;
        pointCoord(234, 1) = 2.810399535165941e-01;
        pointCoord(234, 2) = 1.393272290860796e-03;
        pointCoord(235, 0) = 5.833333333333334e-01;
        pointCoord(235, 1) = 1.172751925766622e-02;
        pointCoord(235, 2) = 1.404011170710178e-03;
        pointCoord(236, 0) = 5.833333333333334e-01;
        pointCoord(236, 1) = 5.769133623677425e-02;
        pointCoord(236, 2) = 2.836318047416551e-03;
        pointCoord(237, 0) = 5.833333333333334e-01;
        pointCoord(237, 1) = 1.249999999999531e-01;
        pointCoord(237, 2) = 3.371193415653452e-03;
        pointCoord(238, 0) = 5.833333333333333e-01;
        pointCoord(238, 1) = 1.923086637631592e-01;
        pointCoord(238, 2) = 2.836318047418840e-03;
        pointCoord(239, 0) = 5.833333333333334e-01;
        pointCoord(239, 1) = 2.382724807423170e-01;
        pointCoord(239, 2) = 1.404011170712085e-03;
        pointCoord(240, 0) = 6.282057758422998e-01;
        pointCoord(240, 1) = 9.622549523018796e-03;
        pointCoord(240, 2) = 9.692277091459844e-04;
        pointCoord(241, 0) = 6.282057758422998e-01;
        pointCoord(241, 1) = 4.733633156260292e-02;
        pointCoord(241, 2) = 1.957988726056997e-03;
        pointCoord(242, 0) = 6.282057758422995e-01;
        pointCoord(242, 1) = 1.025637787455607e-01;
        pointCoord(242, 2) = 2.327227973328982e-03;
        pointCoord(243, 0) = 6.282057758422995e-01;
        pointCoord(243, 1) = 1.577912259284930e-01;
        pointCoord(243, 2) = 1.957988726055198e-03;
        pointCoord(244, 0) = 6.282057758422996e-01;
        pointCoord(244, 1) = 1.955050079680306e-01;
        pointCoord(244, 2) = 9.692277091444852e-04;
        pointCoord(245, 0) = 6.588483204949563e-01;
        pointCoord(245, 1) = 8.185105392941547e-03;
        pointCoord(245, 2) = 4.081084409108241e-04;
        pointCoord(246, 0) = 6.588483204949563e-01;
        pointCoord(246, 1) = 4.026509417575489e-02;
        pointCoord(246, 2) = 8.244416856555859e-04;
        pointCoord(247, 0) = 6.588483204949560e-01;
        pointCoord(247, 1) = 8.724250641921094e-02;
        pointCoord(247, 2) = 9.799156285749623e-04;
        pointCoord(248, 0) = 6.588483204949561e-01;
        pointCoord(248, 1) = 1.342199186626540e-01;
        pointCoord(248, 2) = 8.244416856551300e-04;
        pointCoord(249, 0) = 6.588483204949562e-01;
        pointCoord(249, 1) = 1.662999074454436e-01;
        pointCoord(249, 2) = 4.081084409104437e-04;
        pointCoord(250, 0) = 5.081851053928789e-01;
        pointCoord(250, 1) = 3.329665741121652e-01;
        pointCoord(250, 2) = 4.081084409081491e-04;
        pointCoord(251, 0) = 5.096225495229681e-01;
        pointCoord(251, 1) = 3.621716746347326e-01;
        pointCoord(251, 2) = 9.692277091415620e-04;
        pointCoord(252, 0) = 5.117275192576670e-01;
        pointCoord(252, 1) = 4.049391474089996e-01;
        pointCoord(252, 2) = 1.404011170710363e-03;
        pointCoord(253, 0) = 5.138324889923661e-01;
        pointCoord(253, 1) = 4.477066201832668e-01;
        pointCoord(253, 2) = 1.393272290859976e-03;
        pointCoord(254, 0) = 5.152699331224554e-01;
        pointCoord(254, 1) = 4.769117207058343e-01;
        pointCoord(254, 2) = 7.613571603838885e-04;
        pointCoord(255, 0) = 5.402650941756003e-01;
        pointCoord(255, 1) = 3.008865853294437e-01;
        pointCoord(255, 2) = 8.244416856558055e-04;
        pointCoord(256, 0) = 5.473363315624704e-01;
        pointCoord(256, 1) = 3.244578925952302e-01;
        pointCoord(256, 2) = 1.957988726056759e-03;
        pointCoord(257, 0) = 5.576913362367896e-01;
        pointCoord(257, 1) = 3.589753304298771e-01;
        pointCoord(257, 2) = 2.836318047417838e-03;
        pointCoord(258, 0) = 5.680463409111091e-01;
        pointCoord(258, 1) = 3.934927682645240e-01;
        pointCoord(258, 2) = 2.814623861946163e-03;
        pointCoord(259, 0) = 5.751175782979790e-01;
        pointCoord(259, 1) = 4.170640755303105e-01;
        pointCoord(259, 2) = 1.538058314326909e-03;
        pointCoord(260, 0) = 5.872425064191885e-01;
        pointCoord(260, 1) = 2.539091730858552e-01;
        pointCoord(260, 2) = 9.799156285790383e-04;
        pointCoord(261, 0) = 6.025637787455170e-01;
        pointCoord(261, 1) = 2.692304454121835e-01;
        pointCoord(261, 2) = 2.327227973335010e-03;
        pointCoord(262, 0) = 6.250000000000000e-01;
        pointCoord(262, 1) = 2.916666666666667e-01;
        pointCoord(262, 2) = 3.371193415654718e-03;
        pointCoord(263, 0) = 6.474362212544831e-01;
        pointCoord(263, 1) = 3.141028879211498e-01;
        pointCoord(263, 2) = 3.345408121475985e-03;
        pointCoord(264, 0) = 6.627574935808115e-01;
        pointCoord(264, 1) = 3.294241602474781e-01;
        pointCoord(264, 2) = 1.828106712807099e-03;
        pointCoord(265, 0) = 6.342199186627769e-01;
        pointCoord(265, 1) = 2.069317608422669e-01;
        pointCoord(265, 2) = 8.244416856558036e-04;
        pointCoord(266, 0) = 6.577912259285638e-01;
        pointCoord(266, 1) = 2.140029982291370e-01;
        pointCoord(266, 2) = 1.957988726056757e-03;
        pointCoord(267, 0) = 6.923086637632103e-01;
        pointCoord(267, 1) = 2.243580029034563e-01;
        pointCoord(267, 2) = 2.836318047417835e-03;
        pointCoord(268, 0) = 7.268261015978572e-01;
        pointCoord(268, 1) = 2.347130075777757e-01;
        pointCoord(268, 2) = 2.814623861946169e-03;
        pointCoord(269, 0) = 7.503974088636438e-01;
        pointCoord(269, 1) = 2.417842449646457e-01;
        pointCoord(269, 2) = 1.538058314326903e-03;
        pointCoord(270, 0) = 6.662999074454985e-01;
        pointCoord(270, 1) = 1.748517720595455e-01;
        pointCoord(270, 2) = 4.081084409081470e-04;
        pointCoord(271, 0) = 6.955050079680659e-01;
        pointCoord(271, 1) = 1.762892161896346e-01;
        pointCoord(271, 2) = 9.692277091415605e-04;
        pointCoord(272, 0) = 7.382724807423330e-01;
        pointCoord(272, 1) = 1.783941859243336e-01;
        pointCoord(272, 2) = 1.404011170710362e-03;
        pointCoord(273, 0) = 7.810399535166000e-01;
        pointCoord(273, 1) = 1.804991556590328e-01;
        pointCoord(273, 2) = 1.393272290859974e-03;
        pointCoord(274, 0) = 8.102450540391676e-01;
        pointCoord(274, 1) = 1.819365997891220e-01;
        pointCoord(274, 2) = 7.613571603838920e-04;
        pointCoord(275, 0) = 6.819365997891222e-01;
        pointCoord(275, 1) = 7.818346171771254e-03;
        pointCoord(275, 2) = 7.613571603889068e-04;
        pointCoord(276, 0) = 6.804991556590330e-01;
        pointCoord(276, 1) = 3.846089082451050e-02;
        pointCoord(276, 2) = 1.393272290859161e-03;
        pointCoord(277, 0) = 6.783941859243336e-01;
        pointCoord(277, 1) = 8.333333333333333e-02;
        pointCoord(277, 2) = 1.404011170704043e-03;
        pointCoord(278, 0) = 6.762892161896344e-01;
        pointCoord(278, 1) = 1.282057758421562e-01;
        pointCoord(278, 2) = 9.692277091409926e-04;
        pointCoord(279, 0) = 6.748517720595455e-01;
        pointCoord(279, 1) = 1.588483204948954e-01;
        pointCoord(279, 2) = 4.081084409108353e-04;
        pointCoord(280, 0) = 7.417842449646458e-01;
        pointCoord(280, 1) = 7.818346171771254e-03;
        pointCoord(280, 2) = 1.538058314337037e-03;
        pointCoord(281, 0) = 7.347130075777758e-01;
        pointCoord(281, 1) = 3.846089082451048e-02;
        pointCoord(281, 2) = 2.814623861944516e-03;
        pointCoord(282, 0) = 7.243580029034562e-01;
        pointCoord(282, 1) = 8.333333333333333e-02;
        pointCoord(282, 2) = 2.836318047405072e-03;
        pointCoord(283, 0) = 7.140029982291369e-01;
        pointCoord(283, 1) = 1.282057758421562e-01;
        pointCoord(283, 2) = 1.957988726055609e-03;
        pointCoord(284, 0) = 7.069317608422668e-01;
        pointCoord(284, 1) = 1.588483204948954e-01;
        pointCoord(284, 2) = 8.244416856612340e-04;
        pointCoord(285, 0) = 8.294241602474780e-01;
        pointCoord(285, 1) = 7.818346171771249e-03;
        pointCoord(285, 2) = 1.828106712819136e-03;
        pointCoord(286, 0) = 8.141028879211499e-01;
        pointCoord(286, 1) = 3.846089082451047e-02;
        pointCoord(286, 2) = 3.345408121474028e-03;
        pointCoord(287, 0) = 7.916666666666666e-01;
        pointCoord(287, 1) = 8.333333333333333e-02;
        pointCoord(287, 2) = 3.371193415639546e-03;
        pointCoord(288, 0) = 7.692304454121835e-01;
        pointCoord(288, 1) = 1.282057758421562e-01;
        pointCoord(288, 2) = 2.327227973333647e-03;
        pointCoord(289, 0) = 7.539091730858553e-01;
        pointCoord(289, 1) = 1.588483204948954e-01;
        pointCoord(289, 2) = 9.799156285854940e-04;
        pointCoord(290, 0) = 9.170640755303108e-01;
        pointCoord(290, 1) = 7.818346171771254e-03;
        pointCoord(290, 2) = 1.538058314337038e-03;
        pointCoord(291, 0) = 8.934927682645238e-01;
        pointCoord(291, 1) = 3.846089082451049e-02;
        pointCoord(291, 2) = 2.814623861944519e-03;
        pointCoord(292, 0) = 8.589753304298772e-01;
        pointCoord(292, 1) = 8.333333333333333e-02;
        pointCoord(292, 2) = 2.836318047405070e-03;
        pointCoord(293, 0) = 8.244578925952304e-01;
        pointCoord(293, 1) = 1.282057758421562e-01;
        pointCoord(293, 2) = 1.957988726055611e-03;
        pointCoord(294, 0) = 8.008865853294437e-01;
        pointCoord(294, 1) = 1.588483204948954e-01;
        pointCoord(294, 2) = 8.244416856612343e-04;
        pointCoord(295, 0) = 9.769117207058344e-01;
        pointCoord(295, 1) = 7.818346171771256e-03;
        pointCoord(295, 2) = 7.613571603889076e-04;
        pointCoord(296, 0) = 9.477066201832668e-01;
        pointCoord(296, 1) = 3.846089082451048e-02;
        pointCoord(296, 2) = 1.393272290859161e-03;
        pointCoord(297, 0) = 9.049391474089995e-01;
        pointCoord(297, 1) = 8.333333333333333e-02;
        pointCoord(297, 2) = 1.404011170704043e-03;
        pointCoord(298, 0) = 8.621716746347324e-01;
        pointCoord(298, 1) = 1.282057758421562e-01;
        pointCoord(298, 2) = 9.692277091409905e-04;
        pointCoord(299, 0) = 8.329665741121651e-01;
        pointCoord(299, 1) = 1.588483204948954e-01;
        pointCoord(299, 2) = 4.081084409108356e-04;
    }
    else if (nH == 100)
    {
        pointCoord(0, 0) = 1.702173135062929e-04;
        pointCoord(0, 1) = 1.287651842790784e-02;
        pointCoord(0, 2) = 1.449840738269199e-05;
        pointCoord(1, 0) = 8.802412983224589e-04;
        pointCoord(1, 1) = 1.216649444309167e-02;
        pointCoord(1, 2) = 3.249981782047017e-05;
        pointCoord(2, 0) = 2.091329321814251e-03;
        pointCoord(2, 1) = 1.095540641959988e-02;
        pointCoord(2, 2) = 4.764270720329757e-05;
        pointCoord(3, 0) = 3.696170281331908e-03;
        pointCoord(3, 1) = 9.350565460082221e-03;
        pointCoord(3, 2) = 5.855497037951345e-05;
        pointCoord(4, 0) = 5.552205791021539e-03;
        pointCoord(4, 1) = 7.494529950392589e-03;
        pointCoord(4, 2) = 6.426494989408299e-05;
        pointCoord(5, 0) = 7.494529950392589e-03;
        pointCoord(5, 1) = 5.552205791021539e-03;
        pointCoord(5, 2) = 6.426494989408299e-05;
        pointCoord(6, 0) = 9.350565460082221e-03;
        pointCoord(6, 1) = 3.696170281331908e-03;
        pointCoord(6, 2) = 5.855497037951345e-05;
        pointCoord(7, 0) = 1.095540641959988e-02;
        pointCoord(7, 1) = 2.091329321814251e-03;
        pointCoord(7, 2) = 4.764270720329757e-05;
        pointCoord(8, 0) = 1.216649444309167e-02;
        pointCoord(8, 1) = 8.802412983224589e-04;
        pointCoord(8, 2) = 3.249981782047017e-05;
        pointCoord(9, 0) = 1.287651842790784e-02;
        pointCoord(9, 1) = 1.702173135062929e-04;
        pointCoord(9, 2) = 1.449840738269199e-05;
        pointCoord(10, 0) = 8.802412983224589e-04;
        pointCoord(10, 1) = 6.658807535718528e-02;
        pointCoord(10, 2) = 1.680656405875918e-04;
        pointCoord(11, 0) = 4.551973752327862e-03;
        pointCoord(11, 1) = 6.291634290317988e-02;
        pointCoord(11, 2) = 3.767381172843810e-04;
        pointCoord(12, 0) = 1.081484838136367e-02;
        pointCoord(12, 1) = 5.665346827414406e-02;
        pointCoord(12, 2) = 5.522745977608584e-04;
        pointCoord(13, 0) = 1.911392948367855e-02;
        pointCoord(13, 1) = 4.835438717182918e-02;
        pointCoord(13, 2) = 6.787696294262320e-04;
        pointCoord(14, 0) = 2.871200780560782e-02;
        pointCoord(14, 1) = 3.875630884989991e-02;
        pointCoord(14, 2) = 7.449597522119785e-04;
        pointCoord(15, 0) = 3.875630884989991e-02;
        pointCoord(15, 1) = 2.871200780560782e-02;
        pointCoord(15, 2) = 7.449597522119785e-04;
        pointCoord(16, 0) = 4.835438717182918e-02;
        pointCoord(16, 1) = 1.911392948367855e-02;
        pointCoord(16, 2) = 6.787696294262320e-04;
        pointCoord(17, 0) = 5.665346827414406e-02;
        pointCoord(17, 1) = 1.081484838136367e-02;
        pointCoord(17, 2) = 5.522745977608584e-04;
        pointCoord(18, 0) = 6.291634290317988e-02;
        pointCoord(18, 1) = 4.551973752327862e-03;
        pointCoord(18, 2) = 3.767381172843810e-04;
        pointCoord(19, 0) = 6.658807535718528e-02;
        pointCoord(19, 1) = 8.802412983224589e-04;
        pointCoord(19, 2) = 1.680656405875918e-04;
        pointCoord(20, 0) = 2.091329321814251e-03;
        pointCoord(20, 1) = 1.582038865286735e-01;
        pointCoord(20, 2) = 5.853493307611375e-04;
        pointCoord(21, 0) = 1.081484838136367e-02;
        pointCoord(21, 1) = 1.494803674691241e-01;
        pointCoord(21, 2) = 1.312126643218855e-03;
        pointCoord(22, 0) = 2.569455622455447e-02;
        pointCoord(22, 1) = 1.346006596259333e-01;
        pointCoord(22, 2) = 1.923495873787554e-03;
        pointCoord(23, 0) = 4.541200379996644e-02;
        pointCoord(23, 1) = 1.148832120505213e-01;
        pointCoord(23, 2) = 2.364060535731918e-03;
        pointCoord(24, 0) = 6.821568577441427e-02;
        pointCoord(24, 1) = 9.207953007607352e-02;
        pointCoord(24, 2) = 2.594591558849886e-03;
        pointCoord(25, 0) = 9.207953007607352e-02;
        pointCoord(25, 1) = 6.821568577441427e-02;
        pointCoord(25, 2) = 2.594591558849886e-03;
        pointCoord(26, 0) = 1.148832120505213e-01;
        pointCoord(26, 1) = 4.541200379996644e-02;
        pointCoord(26, 2) = 2.364060535731918e-03;
        pointCoord(27, 0) = 1.346006596259333e-01;
        pointCoord(27, 1) = 2.569455622455447e-02;
        pointCoord(27, 2) = 1.923495873787554e-03;
        pointCoord(28, 0) = 1.494803674691241e-01;
        pointCoord(28, 1) = 1.081484838136367e-02;
        pointCoord(28, 2) = 1.312126643218855e-03;
        pointCoord(29, 0) = 1.582038865286735e-01;
        pointCoord(29, 1) = 2.091329321814251e-03;
        pointCoord(29, 2) = 5.853493307611375e-04;
        pointCoord(30, 0) = 3.696170281331908e-03;
        pointCoord(30, 1) = 2.796061326540445e-01;
        pointCoord(30, 2) = 1.271487235245470e-03;
        pointCoord(31, 0) = 1.911392948367855e-02;
        pointCoord(31, 1) = 2.641883734516979e-01;
        pointCoord(31, 2) = 2.850182259043298e-03;
        pointCoord(32, 0) = 4.541200379996644e-02;
        pointCoord(32, 1) = 2.378902991354100e-01;
        pointCoord(32, 2) = 4.178189539207349e-03;
        pointCoord(33, 0) = 8.026019484848777e-02;
        pointCoord(33, 1) = 2.030421080868886e-01;
        pointCoord(33, 2) = 5.135177639345937e-03;
        pointCoord(34, 0) = 1.205629299269492e-01;
        pointCoord(34, 1) = 1.627393730084272e-01;
        pointCoord(34, 2) = 5.635933748251760e-03;
        pointCoord(35, 0) = 1.627393730084272e-01;
        pointCoord(35, 1) = 1.205629299269492e-01;
        pointCoord(35, 2) = 5.635933748251760e-03;
        pointCoord(36, 0) = 2.030421080868886e-01;
        pointCoord(36, 1) = 8.026019484848777e-02;
        pointCoord(36, 2) = 5.135177639345937e-03;
        pointCoord(37, 0) = 2.378902991354100e-01;
        pointCoord(37, 1) = 4.541200379996644e-02;
        pointCoord(37, 2) = 4.178189539207349e-03;
        pointCoord(38, 0) = 2.641883734516979e-01;
        pointCoord(38, 1) = 1.911392948367855e-02;
        pointCoord(38, 2) = 2.850182259043298e-03;
        pointCoord(39, 0) = 2.796061326540445e-01;
        pointCoord(39, 1) = 3.696170281331908e-03;
        pointCoord(39, 2) = 1.271487235245470e-03;
        pointCoord(40, 0) = 5.552205791021539e-03;
        pointCoord(40, 1) = 4.200106247181629e-01;
        pointCoord(40, 2) = 2.096215829116851e-03;
        pointCoord(41, 0) = 2.871200780560782e-02;
        pointCoord(41, 1) = 3.968508227035766e-01;
        pointCoord(41, 2) = 4.698904559683718e-03;
        pointCoord(42, 0) = 6.821568577441427e-02;
        pointCoord(42, 1) = 3.573471447347701e-01;
        pointCoord(42, 2) = 6.888301200637704e-03;
        pointCoord(43, 0) = 1.205629299269492e-01;
        pointCoord(43, 1) = 3.049999005822352e-01;
        pointCoord(43, 2) = 8.466023373679956e-03;
        pointCoord(44, 0) = 1.811037227109888e-01;
        pointCoord(44, 1) = 2.444591077981956e-01;
        pointCoord(44, 2) = 9.291586425292262e-03;
        pointCoord(45, 0) = 2.444591077981956e-01;
        pointCoord(45, 1) = 1.811037227109888e-01;
        pointCoord(45, 2) = 9.291586425292262e-03;
        pointCoord(46, 0) = 3.049999005822352e-01;
        pointCoord(46, 1) = 1.205629299269492e-01;
        pointCoord(46, 2) = 8.466023373679956e-03;
        pointCoord(47, 0) = 3.573471447347701e-01;
        pointCoord(47, 1) = 6.821568577441427e-02;
        pointCoord(47, 2) = 6.888301200637704e-03;
        pointCoord(48, 0) = 3.968508227035766e-01;
        pointCoord(48, 1) = 2.871200780560782e-02;
        pointCoord(48, 2) = 4.698904559683718e-03;
        pointCoord(49, 0) = 4.200106247181629e-01;
        pointCoord(49, 1) = 5.552205791021539e-03;
        pointCoord(49, 2) = 2.096215829116851e-03;
        pointCoord(50, 0) = 7.494529950392589e-03;
        pointCoord(50, 1) = 5.669426395404230e-01;
        pointCoord(50, 2) = 2.829533505261662e-03;
        pointCoord(51, 0) = 3.875630884989991e-02;
        pointCoord(51, 1) = 5.356808606409157e-01;
        pointCoord(51, 2) = 6.342718962891070e-03;
        pointCoord(52, 0) = 9.207953007607352e-02;
        pointCoord(52, 1) = 4.823576394147421e-01;
        pointCoord(52, 2) = 9.298030656390024e-03;
        pointCoord(53, 0) = 1.627393730084272e-01;
        pointCoord(53, 1) = 4.116977964823884e-01;
        pointCoord(53, 2) = 1.142768624271297e-02;
        pointCoord(54, 0) = 2.444591077981956e-01;
        pointCoord(54, 1) = 3.299780616926200e-01;
        pointCoord(54, 2) = 1.254205542302168e-02;
        pointCoord(55, 0) = 3.299780616926200e-01;
        pointCoord(55, 1) = 2.444591077981956e-01;
        pointCoord(55, 2) = 1.254205542302168e-02;
        pointCoord(56, 0) = 4.116977964823884e-01;
        pointCoord(56, 1) = 1.627393730084272e-01;
        pointCoord(56, 2) = 1.142768624271297e-02;
        pointCoord(57, 0) = 4.823576394147421e-01;
        pointCoord(57, 1) = 9.207953007607352e-02;
        pointCoord(57, 2) = 9.298030656390024e-03;
        pointCoord(58, 0) = 5.356808606409157e-01;
        pointCoord(58, 1) = 3.875630884989991e-02;
        pointCoord(58, 2) = 6.342718962891070e-03;
        pointCoord(59, 0) = 5.669426395404230e-01;
        pointCoord(59, 1) = 7.494529950392589e-03;
        pointCoord(59, 2) = 2.829533505261662e-03;
        pointCoord(60, 0) = 9.350565460082221e-03;
        pointCoord(60, 1) = 7.073471316045414e-01;
        pointCoord(60, 2) = 3.216606303251134e-03;
        pointCoord(61, 0) = 4.835438717182918e-02;
        pointCoord(61, 1) = 6.683433098927944e-01;
        pointCoord(61, 2) = 7.210386361514114e-03;
        pointCoord(62, 0) = 1.148832120505213e-01;
        pointCoord(62, 1) = 6.018144850141023e-01;
        pointCoord(62, 2) = 1.056997698085242e-02;
        pointCoord(63, 0) = 2.030421080868886e-01;
        pointCoord(63, 1) = 5.136555889777350e-01;
        pointCoord(63, 2) = 1.299096389264618e-02;
        pointCoord(64, 0) = 3.049999005822352e-01;
        pointCoord(64, 1) = 4.116977964823884e-01;
        pointCoord(64, 2) = 1.425777586814117e-02;
        pointCoord(65, 0) = 4.116977964823884e-01;
        pointCoord(65, 1) = 3.049999005822352e-01;
        pointCoord(65, 2) = 1.425777586814117e-02;
        pointCoord(66, 0) = 5.136555889777350e-01;
        pointCoord(66, 1) = 2.030421080868886e-01;
        pointCoord(66, 2) = 1.299096389264618e-02;
        pointCoord(67, 0) = 6.018144850141023e-01;
        pointCoord(67, 1) = 1.148832120505213e-01;
        pointCoord(67, 2) = 1.056997698085242e-02;
        pointCoord(68, 0) = 6.683433098927944e-01;
        pointCoord(68, 1) = 4.835438717182918e-02;
        pointCoord(68, 2) = 7.210386361514114e-03;
        pointCoord(69, 0) = 7.073471316045414e-01;
        pointCoord(69, 1) = 9.350565460082221e-03;
        pointCoord(69, 2) = 3.216606303251134e-03;
        pointCoord(70, 0) = 1.095540641959988e-02;
        pointCoord(70, 1) = 8.287493777299123e-01;
        pointCoord(70, 2) = 3.066346246398886e-03;
        pointCoord(71, 0) = 5.665346827414406e-02;
        pointCoord(71, 1) = 7.830513158753682e-01;
        pointCoord(71, 2) = 6.873561471407818e-03;
        pointCoord(72, 0) = 1.346006596259333e-01;
        pointCoord(72, 1) = 7.051041245235790e-01;
        pointCoord(72, 2) = 1.007621268633352e-02;
        pointCoord(73, 0) = 2.378902991354100e-01;
        pointCoord(73, 1) = 6.018144850141023e-01;
        pointCoord(73, 2) = 1.238410598432785e-02;
        pointCoord(74, 0) = 3.573471447347701e-01;
        pointCoord(74, 1) = 4.823576394147421e-01;
        pointCoord(74, 2) = 1.359174029817784e-02;
        pointCoord(75, 0) = 4.823576394147421e-01;
        pointCoord(75, 1) = 3.573471447347701e-01;
        pointCoord(75, 2) = 1.359174029817784e-02;
        pointCoord(76, 0) = 6.018144850141023e-01;
        pointCoord(76, 1) = 2.378902991354100e-01;
        pointCoord(76, 2) = 1.238410598432785e-02;
        pointCoord(77, 0) = 7.051041245235790e-01;
        pointCoord(77, 1) = 1.346006596259333e-01;
        pointCoord(77, 2) = 1.007621268633352e-02;
        pointCoord(78, 0) = 7.830513158753682e-01;
        pointCoord(78, 1) = 5.665346827414406e-02;
        pointCoord(78, 2) = 6.873561471407818e-03;
        pointCoord(79, 0) = 8.287493777299123e-01;
        pointCoord(79, 1) = 1.095540641959988e-02;
        pointCoord(79, 2) = 3.066346246398886e-03;
        pointCoord(80, 0) = 1.216649444309167e-02;
        pointCoord(80, 1) = 9.203651889014006e-01;
        pointCoord(80, 2) = 2.322964948566316e-03;
        pointCoord(81, 0) = 6.291634290317988e-02;
        pointCoord(81, 1) = 8.696153404413125e-01;
        pointCoord(81, 2) = 5.207188323447800e-03;
        pointCoord(82, 0) = 1.494803674691241e-01;
        pointCoord(82, 1) = 7.830513158753682e-01;
        pointCoord(82, 2) = 7.633413516865814e-03;
        pointCoord(83, 0) = 2.641883734516979e-01;
        pointCoord(83, 1) = 6.683433098927944e-01;
        pointCoord(83, 2) = 9.381798991131181e-03;
        pointCoord(84, 0) = 3.968508227035766e-01;
        pointCoord(84, 1) = 5.356808606409157e-01;
        pointCoord(84, 2) = 1.029666377036281e-02;
        pointCoord(85, 0) = 5.356808606409157e-01;
        pointCoord(85, 1) = 3.968508227035766e-01;
        pointCoord(85, 2) = 1.029666377036281e-02;
        pointCoord(86, 0) = 6.683433098927944e-01;
        pointCoord(86, 1) = 2.641883734516979e-01;
        pointCoord(86, 2) = 9.381798991131181e-03;
        pointCoord(87, 0) = 7.830513158753682e-01;
        pointCoord(87, 1) = 1.494803674691241e-01;
        pointCoord(87, 2) = 7.633413516865814e-03;
        pointCoord(88, 0) = 8.696153404413125e-01;
        pointCoord(88, 1) = 6.291634290317988e-02;
        pointCoord(88, 2) = 5.207188323447800e-03;
        pointCoord(89, 0) = 9.203651889014006e-01;
        pointCoord(89, 1) = 1.216649444309167e-02;
        pointCoord(89, 2) = 2.322964948566316e-03;
        pointCoord(90, 0) = 1.287651842790784e-02;
        pointCoord(90, 1) = 9.740767458306780e-01;
        pointCoord(90, 2) = 1.096768630599064e-03;
        pointCoord(91, 0) = 6.658807535718528e-02;
        pointCoord(91, 1) = 9.203651889014006e-01;
        pointCoord(91, 2) = 2.458530771333437e-03;
        pointCoord(92, 0) = 1.582038865286735e-01;
        pointCoord(92, 1) = 8.287493777299123e-01;
        pointCoord(92, 2) = 3.604052869956726e-03;
        pointCoord(93, 0) = 2.796061326540445e-01;
        pointCoord(93, 1) = 7.073471316045414e-01;
        pointCoord(93, 2) = 4.429538568117091e-03;
        pointCoord(94, 0) = 4.200106247181629e-01;
        pointCoord(94, 1) = 5.669426395404230e-01;
        pointCoord(94, 2) = 4.861484384484430e-03;
        pointCoord(95, 0) = 5.669426395404230e-01;
        pointCoord(95, 1) = 4.200106247181629e-01;
        pointCoord(95, 2) = 4.861484384484430e-03;
        pointCoord(96, 0) = 7.073471316045414e-01;
        pointCoord(96, 1) = 2.796061326540445e-01;
        pointCoord(96, 2) = 4.429538568117091e-03;
        pointCoord(97, 0) = 8.287493777299123e-01;
        pointCoord(97, 1) = 1.582038865286735e-01;
        pointCoord(97, 2) = 3.604052869956726e-03;
        pointCoord(98, 0) = 9.203651889014006e-01;
        pointCoord(98, 1) = 6.658807535718528e-02;
        pointCoord(98, 2) = 2.458530771333437e-03;
        pointCoord(99, 0) = 9.740767458306780e-01;
        pointCoord(99, 1) = 1.287651842790784e-02;
        pointCoord(99, 2) = 1.096768630599064e-03;
    }
    else if (nH == 225)
    {
        pointCoord(0, 0) = 5.985718536821210e-03;
        pointCoord(0, 1) = 5.985718536821210e-03;
        pointCoord(0, 2) = 2.350209459566036e-04;
        pointCoord(1, 0) = 5.909592413459181e-03;
        pointCoord(1, 1) = 3.126915522334901e-02;
        pointCoord(1, 2) = 5.308883103915073e-04;
        pointCoord(2, 0) = 5.775908900468726e-03;
        pointCoord(2, 1) = 7.566887620549786e-02;
        pointCoord(2, 2) = 7.901355606254923e-04;
        pointCoord(3, 0) = 5.590109849186445e-03;
        pointCoord(3, 1) = 1.373775031793442e-01;
        pointCoord(3, 2) = 9.959122492752457e-04;
        pointCoord(4, 0) = 5.359798001493101e-03;
        pointCoord(4, 1) = 2.138699707074664e-01;
        pointCoord(4, 2) = 1.137382193294796e-03;
        pointCoord(5, 0) = 5.094401391972344e-03;
        pointCoord(5, 1) = 3.020149868634334e-01;
        pointCoord(5, 2) = 1.210184826439239e-03;
        pointCoord(6, 0) = 4.804785049575310e-03;
        pointCoord(6, 1) = 3.982039970611008e-01;
        pointCoord(6, 2) = 1.216358063642314e-03;
        pointCoord(7, 0) = 4.502805742317942e-03;
        pointCoord(7, 1) = 4.984990647525607e-01;
        pointCoord(7, 2) = 1.163437951663850e-03;
        pointCoord(8, 0) = 4.200826435060573e-03;
        pointCoord(8, 1) = 5.987941324440206e-01;
        pointCoord(8, 2) = 1.062886919828535e-03;
        pointCoord(9, 0) = 3.911210092663540e-03;
        pointCoord(9, 1) = 6.949831426416880e-01;
        pointCoord(9, 2) = 9.281175941868199e-04;
        pointCoord(10, 0) = 3.645813483142782e-03;
        pointCoord(10, 1) = 7.831281587976550e-01;
        pointCoord(10, 2) = 7.724369699259093e-04;
        pointCoord(11, 0) = 3.415501635449438e-03;
        pointCoord(11, 1) = 8.596206263257772e-01;
        pointCoord(11, 2) = 6.072394228417237e-04;
        pointCoord(12, 0) = 3.229702584167157e-03;
        pointCoord(12, 1) = 9.213292532996236e-01;
        pointCoord(12, 2) = 4.407281567922316e-04;
        pointCoord(13, 0) = 3.096019071176703e-03;
        pointCoord(13, 1) = 9.657289742817724e-01;
        pointCoord(13, 2) = 2.773577221348077e-04;
        pointCoord(14, 0) = 3.019892947814674e-03;
        pointCoord(14, 1) = 9.910124109683001e-01;
        pointCoord(14, 2) = 1.182202366589500e-04;
        pointCoord(15, 0) = 3.126915522334901e-02;
        pointCoord(15, 1) = 5.909592413459181e-03;
        pointCoord(15, 2) = 5.308883103915073e-04;
        pointCoord(16, 0) = 3.087147538703260e-02;
        pointCoord(16, 1) = 3.087147538703260e-02;
        pointCoord(16, 2) = 1.199022245971729e-03;
        pointCoord(17, 0) = 3.017311803982579e-02;
        pointCoord(17, 1) = 7.470652253496511e-02;
        pointCoord(17, 2) = 1.783995166196046e-03;
        pointCoord(18, 0) = 2.920251119636035e-02;
        pointCoord(18, 1) = 1.356303417166282e-01;
        pointCoord(18, 2) = 2.247600239327195e-03;
        pointCoord(19, 0) = 2.799937127740184e-02;
        pointCoord(19, 1) = 2.111499811734854e-01;
        pointCoord(19, 2) = 2.565340287901570e-03;
        pointCoord(20, 0) = 2.661294996009375e-02;
        pointCoord(20, 1) = 2.981739726216650e-01;
        pointCoord(20, 2) = 2.727481943612962e-03;
        pointCoord(21, 0) = 2.510000572291938e-02;
        pointCoord(21, 1) = 3.931396549245550e-01;
        pointCoord(21, 2) = 2.738869826410303e-03;
        pointCoord(22, 0) = 2.352247784973531e-02;
        pointCoord(22, 1) = 4.921591740500882e-01;
        pointCoord(22, 2) = 2.616859123730370e-03;
        pointCoord(23, 0) = 2.194494997655123e-02;
        pointCoord(23, 1) = 5.911786931756214e-01;
        pointCoord(23, 2) = 2.387714725591376e-03;
        pointCoord(24, 0) = 2.043200573937687e-02;
        pointCoord(24, 1) = 6.861443754785116e-01;
        pointCoord(24, 2) = 2.082087995813234e-03;
        pointCoord(25, 0) = 1.904558442206878e-02;
        pointCoord(25, 1) = 7.731683669266911e-01;
        pointCoord(25, 2) = 1.730314409284368e-03;
        pointCoord(26, 0) = 1.784244450311027e-02;
        pointCoord(26, 1) = 8.486880063835482e-01;
        pointCoord(26, 2) = 1.358283575880347e-03;
        pointCoord(27, 0) = 1.687183765964483e-02;
        pointCoord(27, 1) = 9.096118255652115e-01;
        pointCoord(27, 2) = 9.845211473614196e-04;
        pointCoord(28, 0) = 1.617348031243802e-02;
        pointCoord(28, 1) = 9.534468727131440e-01;
        pointCoord(28, 2) = 6.189225798873323e-04;
        pointCoord(29, 0) = 1.577580047612162e-02;
        pointCoord(29, 1) = 9.784087556867173e-01;
        pointCoord(29, 2) = 2.636382995431647e-04;
        pointCoord(30, 0) = 7.566887620549786e-02;
        pointCoord(30, 1) = 5.775908900468726e-03;
        pointCoord(30, 2) = 7.901355606254923e-04;
        pointCoord(31, 0) = 7.470652253496511e-02;
        pointCoord(31, 1) = 3.017311803982579e-02;
        pointCoord(31, 2) = 1.783995166196046e-03;
        pointCoord(32, 0) = 7.301655312979445e-02;
        pointCoord(32, 1) = 7.301655312979445e-02;
        pointCoord(32, 2) = 2.652892287898599e-03;
        pointCoord(33, 0) = 7.066776153124324e-02;
        pointCoord(33, 1) = 1.325621875563718e-01;
        pointCoord(33, 2) = 3.339574226016945e-03;
        pointCoord(34, 0) = 6.775625832831747e-02;
        pointCoord(34, 1) = 2.063734637292617e-01;
        pointCoord(34, 2) = 3.807528754956306e-03;
        pointCoord(35, 0) = 6.440122867437553e-02;
        pointCoord(35, 1) = 2.914288468408074e-01;
        pointCoord(35, 2) = 4.042585575420045e-03;
        pointCoord(36, 0) = 6.074002358677907e-02;
        pointCoord(36, 1) = 3.842462682932755e-01;
        pointCoord(36, 2) = 4.052608385115568e-03;
        pointCoord(37, 0) = 5.692253122108980e-02;
        pointCoord(37, 1) = 4.810258229263034e-01;
        pointCoord(37, 2) = 3.864326798602145e-03;
        pointCoord(38, 0) = 5.310503885540051e-02;
        pointCoord(38, 1) = 5.778053775593315e-01;
        pointCoord(38, 2) = 3.517840436508153e-03;
        pointCoord(39, 0) = 4.944383376780406e-02;
        pointCoord(39, 1) = 6.706227990117994e-01;
        pointCoord(39, 2) = 3.059726443647184e-03;
        pointCoord(40, 0) = 4.608880411386213e-02;
        pointCoord(40, 1) = 7.556781821233451e-01;
        pointCoord(40, 2) = 2.535882479276867e-03;
        pointCoord(41, 0) = 4.317730091093636e-02;
        pointCoord(41, 1) = 8.294894582962350e-01;
        pointCoord(41, 2) = 1.985249465633941e-03;
        pointCoord(42, 0) = 4.082850931238514e-02;
        pointCoord(42, 1) = 8.890350927228124e-01;
        pointCoord(42, 2) = 1.435387316391502e-03;
        pointCoord(43, 0) = 3.913853990721448e-02;
        pointCoord(43, 1) = 9.318785278127811e-01;
        pointCoord(43, 2) = 9.005715516839962e-04;
        pointCoord(44, 0) = 3.817618623668173e-02;
        pointCoord(44, 1) = 9.562757369521381e-01;
        pointCoord(44, 2) = 3.831452029932598e-04;
        pointCoord(45, 0) = 1.373775031793442e-01;
        pointCoord(45, 1) = 5.590109849186445e-03;
        pointCoord(45, 2) = 9.959122492752457e-04;
        pointCoord(46, 0) = 1.356303417166282e-01;
        pointCoord(46, 1) = 2.920251119636035e-02;
        pointCoord(46, 2) = 2.247600239327195e-03;
        pointCoord(47, 0) = 1.325621875563718e-01;
        pointCoord(47, 1) = 7.066776153124324e-02;
        pointCoord(47, 2) = 3.339574226016945e-03;
        pointCoord(48, 0) = 1.282979359713305e-01;
        pointCoord(48, 1) = 1.282979359713305e-01;
        pointCoord(48, 2) = 4.198951600984662e-03;
        pointCoord(49, 0) = 1.230120765721454e-01;
        pointCoord(49, 1) = 1.997348559479611e-01;
        pointCoord(49, 2) = 4.779614109778093e-03;
        pointCoord(50, 0) = 1.169209910418212e-01;
        pointCoord(50, 1) = 2.820541831831245e-01;
        pointCoord(50, 2) = 5.064286617308996e-03;
        pointCoord(51, 0) = 1.102740413475297e-01;
        pointCoord(51, 1) = 3.718858600288976e-01;
        pointCoord(51, 2) = 5.064090884582665e-03;
        pointCoord(52, 0) = 1.033433507399362e-01;
        pointCoord(52, 1) = 4.655522164200213e-01;
        pointCoord(52, 2) = 4.814383713682779e-03;
        pointCoord(53, 0) = 9.641266013234269e-02;
        pointCoord(53, 1) = 5.592185728111451e-01;
        pointCoord(53, 2) = 4.367576624821128e-03;
        pointCoord(54, 0) = 8.976571043805123e-02;
        pointCoord(54, 1) = 6.490502496569180e-01;
        pointCoord(54, 2) = 3.784151182118204e-03;
        pointCoord(55, 0) = 8.367462490772701e-02;
        pointCoord(55, 1) = 7.313695768920815e-01;
        pointCoord(55, 2) = 3.123344756519591e-03;
        pointCoord(56, 0) = 7.838876550854190e-02;
        pointCoord(56, 1) = 8.028064968687120e-01;
        pointCoord(56, 2) = 2.434996767095785e-03;
        pointCoord(57, 0) = 7.412451392350064e-02;
        pointCoord(57, 1) = 8.604366713087993e-01;
        pointCoord(57, 2) = 1.753821796026420e-03;
        pointCoord(58, 0) = 7.105635976324422e-02;
        pointCoord(58, 1) = 9.019019216436822e-01;
        pointCoord(58, 2) = 1.096975661844450e-03;
        pointCoord(59, 0) = 6.930919830052830e-02;
        pointCoord(59, 1) = 9.255143229908560e-01;
        pointCoord(59, 2) = 4.658232856144815e-04;
        pointCoord(60, 0) = 2.138699707074664e-01;
        pointCoord(60, 1) = 5.359798001493101e-03;
        pointCoord(60, 2) = 1.137382193294796e-03;
        pointCoord(61, 0) = 2.111499811734854e-01;
        pointCoord(61, 1) = 2.799937127740184e-02;
        pointCoord(61, 2) = 2.565340287901570e-03;
        pointCoord(62, 0) = 2.063734637292617e-01;
        pointCoord(62, 1) = 6.775625832831747e-02;
        pointCoord(62, 2) = 3.807528754956306e-03;
        pointCoord(63, 0) = 1.997348559479611e-01;
        pointCoord(63, 1) = 1.230120765721454e-01;
        pointCoord(63, 2) = 4.779614109778093e-03;
        pointCoord(64, 0) = 1.915058041112009e-01;
        pointCoord(64, 1) = 1.915058041112009e-01;
        pointCoord(64, 2) = 5.428778846101111e-03;
        pointCoord(65, 0) = 1.820231722843110e-01;
        pointCoord(65, 1) = 2.704335850497988e-01;
        pointCoord(65, 2) = 5.736187419275157e-03;
        pointCoord(66, 0) = 1.716751684007620e-01;
        pointCoord(66, 1) = 3.565642077063142e-01;
        pointCoord(66, 2) = 5.716387981049921e-03;
        pointCoord(67, 0) = 1.608854352717979e-01;
        pointCoord(67, 1) = 4.463715215760674e-01;
        pointCoord(67, 2) = 5.412301897912824e-03;
        pointCoord(68, 0) = 1.500957021428338e-01;
        pointCoord(68, 1) = 5.361788354458206e-01;
        pointCoord(68, 2) = 4.886637244494526e-03;
        pointCoord(69, 0) = 1.397476982592848e-01;
        pointCoord(69, 1) = 6.223094581023360e-01;
        pointCoord(69, 2) = 4.211174391395614e-03;
        pointCoord(70, 0) = 1.302650664323950e-01;
        pointCoord(70, 1) = 7.012372390409338e-01;
        pointCoord(70, 2) = 3.455681100376738e-03;
        pointCoord(71, 0) = 1.220360145956348e-01;
        pointCoord(71, 1) = 7.697309665799893e-01;
        pointCoord(71, 2) = 2.678231687010688e-03;
        pointCoord(72, 0) = 1.153974068143342e-01;
        pointCoord(72, 1) = 8.249867848238173e-01;
        pointCoord(72, 2) = 1.918437141484353e-03;
        pointCoord(73, 0) = 1.106208893701105e-01;
        pointCoord(73, 1) = 8.647436718747329e-01;
        pointCoord(73, 2) = 1.194612298797779e-03;
        pointCoord(74, 0) = 1.079008998361295e-01;
        pointCoord(74, 1) = 8.873832451506416e-01;
        pointCoord(74, 2) = 5.058923108220498e-04;
        pointCoord(75, 0) = 3.020149868634334e-01;
        pointCoord(75, 1) = 5.094401391972344e-03;
        pointCoord(75, 2) = 1.210184826439239e-03;
        pointCoord(76, 0) = 2.981739726216650e-01;
        pointCoord(76, 1) = 2.661294996009375e-02;
        pointCoord(76, 2) = 2.727481943612962e-03;
        pointCoord(77, 0) = 2.914288468408074e-01;
        pointCoord(77, 1) = 6.440122867437553e-02;
        pointCoord(77, 2) = 4.042585575420045e-03;
        pointCoord(78, 0) = 2.820541831831245e-01;
        pointCoord(78, 1) = 1.169209910418212e-01;
        pointCoord(78, 2) = 5.064286617308996e-03;
        pointCoord(79, 0) = 2.704335850497988e-01;
        pointCoord(79, 1) = 1.820231722843110e-01;
        pointCoord(79, 2) = 5.736187419275157e-03;
        pointCoord(80, 0) = 2.570427526802269e-01;
        pointCoord(80, 1) = 2.570427526802269e-01;
        pointCoord(80, 2) = 6.039449332731192e-03;
        pointCoord(81, 0) = 2.424298911989507e-01;
        pointCoord(81, 1) = 3.389085177390151e-01;
        pointCoord(81, 2) = 5.992036836059638e-03;
        pointCoord(82, 0) = 2.271932448459137e-01;
        pointCoord(82, 1) = 4.242689183846954e-01;
        pointCoord(82, 2) = 5.643039878353591e-03;
        pointCoord(83, 0) = 2.119565984928767e-01;
        pointCoord(83, 1) = 5.096293190303758e-01;
        pointCoord(83, 2) = 5.063017983038727e-03;
        pointCoord(84, 0) = 1.973437370116005e-01;
        pointCoord(84, 1) = 5.914950840891640e-01;
        pointCoord(84, 2) = 4.331989740849268e-03;
        pointCoord(85, 0) = 1.839529046420287e-01;
        pointCoord(85, 1) = 6.665146644850798e-01;
        pointCoord(85, 2) = 3.527036015275618e-03;
        pointCoord(86, 0) = 1.723323065087029e-01;
        pointCoord(86, 1) = 7.316168457275697e-01;
        pointCoord(86, 2) = 2.711503061777502e-03;
        pointCoord(87, 0) = 1.629576428510200e-01;
        pointCoord(87, 1) = 7.841366080950154e-01;
        pointCoord(87, 2) = 1.927490467117033e-03;
        pointCoord(88, 0) = 1.562125170701625e-01;
        pointCoord(88, 1) = 8.219248868092972e-01;
        pointCoord(88, 2) = 1.192765446068943e-03;
        pointCoord(89, 0) = 1.523715028283941e-01;
        pointCoord(89, 1) = 8.434434353774185e-01;
        pointCoord(89, 2) = 5.031459767425094e-04;
        pointCoord(90, 0) = 3.982039970611008e-01;
        pointCoord(90, 1) = 4.804785049575310e-03;
        pointCoord(90, 2) = 1.216358063642314e-03;
        pointCoord(91, 0) = 3.931396549245550e-01;
        pointCoord(91, 1) = 2.510000572291938e-02;
        pointCoord(91, 2) = 2.738869826410303e-03;
        pointCoord(92, 0) = 3.842462682932755e-01;
        pointCoord(92, 1) = 6.074002358677907e-02;
        pointCoord(92, 2) = 4.052608385115568e-03;
        pointCoord(93, 0) = 3.718858600288976e-01;
        pointCoord(93, 1) = 1.102740413475297e-01;
        pointCoord(93, 2) = 5.064090884582665e-03;
        pointCoord(94, 0) = 3.565642077063142e-01;
        pointCoord(94, 1) = 1.716751684007620e-01;
        pointCoord(94, 2) = 5.716387981049921e-03;
        pointCoord(95, 0) = 3.389085177390151e-01;
        pointCoord(95, 1) = 2.424298911989507e-01;
        pointCoord(95, 2) = 5.992036836059638e-03;
        pointCoord(96, 0) = 3.196415935682103e-01;
        pointCoord(96, 1) = 3.196415935682103e-01;
        pointCoord(96, 2) = 5.912135344877402e-03;
        pointCoord(97, 0) = 2.995522147509621e-01;
        pointCoord(97, 1) = 4.001492617496793e-01;
        pointCoord(97, 2) = 5.530211452584761e-03;
        pointCoord(98, 0) = 2.794628359337138e-01;
        pointCoord(98, 1) = 4.806569299311483e-01;
        pointCoord(98, 2) = 4.921881796140457e-03;
        pointCoord(99, 0) = 2.601959117629090e-01;
        pointCoord(99, 1) = 5.578686323004080e-01;
        pointCoord(99, 2) = 4.172032952728380e-03;
        pointCoord(100, 0) = 2.425402217956099e-01;
        pointCoord(100, 1) = 6.286233550985966e-01;
        pointCoord(100, 2) = 3.361624124111803e-03;
        pointCoord(101, 0) = 2.272185694730266e-01;
        pointCoord(101, 1) = 6.900244821518289e-01;
        pointCoord(101, 2) = 2.556227607428541e-03;
        pointCoord(102, 0) = 2.148581612086487e-01;
        pointCoord(102, 1) = 7.395584999125796e-01;
        pointCoord(102, 2) = 1.798100384440942e-03;
        pointCoord(103, 0) = 2.059647745773691e-01;
        pointCoord(103, 1) = 7.751985177764392e-01;
        pointCoord(103, 2) = 1.102995104324453e-03;
        pointCoord(104, 0) = 2.009004324408233e-01;
        pointCoord(104, 1) = 7.954937384497832e-01;
        pointCoord(104, 2) = 4.627159521497556e-04;
        pointCoord(105, 0) = 4.984990647525607e-01;
        pointCoord(105, 1) = 4.502805742317942e-03;
        pointCoord(105, 2) = 1.163437951663850e-03;
        pointCoord(106, 0) = 4.921591740500882e-01;
        pointCoord(106, 1) = 2.352247784973531e-02;
        pointCoord(106, 2) = 2.616859123730370e-03;
        pointCoord(107, 0) = 4.810258229263034e-01;
        pointCoord(107, 1) = 5.692253122108980e-02;
        pointCoord(107, 2) = 3.864326798602145e-03;
        pointCoord(108, 0) = 4.655522164200213e-01;
        pointCoord(108, 1) = 1.033433507399362e-01;
        pointCoord(108, 2) = 4.814383713682779e-03;
        pointCoord(109, 0) = 4.463715215760674e-01;
        pointCoord(109, 1) = 1.608854352717979e-01;
        pointCoord(109, 2) = 5.412301897912824e-03;
        pointCoord(110, 0) = 4.242689183846954e-01;
        pointCoord(110, 1) = 2.271932448459137e-01;
        pointCoord(110, 2) = 5.643039878353591e-03;
        pointCoord(111, 0) = 4.001492617496793e-01;
        pointCoord(111, 1) = 2.995522147509621e-01;
        pointCoord(111, 2) = 5.530211452584761e-03;
        pointCoord(112, 0) = 3.750000000000000e-01;
        pointCoord(112, 1) = 3.750000000000000e-01;
        pointCoord(112, 2) = 5.129743012706404e-03;
        pointCoord(113, 0) = 3.498507382503206e-01;
        pointCoord(113, 1) = 4.504477852490379e-01;
        pointCoord(113, 2) = 4.519263907476252e-03;
        pointCoord(114, 0) = 3.257310816153046e-01;
        pointCoord(114, 1) = 5.228067551540863e-01;
        pointCoord(114, 2) = 3.785002146210653e-03;
        pointCoord(115, 0) = 3.036284784239326e-01;
        pointCoord(115, 1) = 5.891145647282021e-01;
        pointCoord(115, 2) = 3.008328952278663e-03;
        pointCoord(116, 0) = 2.844477835799787e-01;
        pointCoord(116, 1) = 6.466566492600638e-01;
        pointCoord(116, 2) = 2.254111925976988e-03;
        pointCoord(117, 0) = 2.689741770736966e-01;
        pointCoord(117, 1) = 6.930774687789103e-01;
        pointCoord(117, 2) = 1.562704823486195e-03;
        pointCoord(118, 0) = 2.578408259499118e-01;
        pointCoord(118, 1) = 7.264775221502647e-01;
        pointCoord(118, 2) = 9.467984241174950e-04;
        pointCoord(119, 0) = 2.515009352474393e-01;
        pointCoord(119, 1) = 7.454971942576820e-01;
        pointCoord(119, 2) = 3.940464726072197e-04;
        pointCoord(120, 0) = 5.987941324440206e-01;
        pointCoord(120, 1) = 4.200826435060573e-03;
        pointCoord(120, 2) = 1.062886919828536e-03;
        pointCoord(121, 0) = 5.911786931756214e-01;
        pointCoord(121, 1) = 2.194494997655123e-02;
        pointCoord(121, 2) = 2.387714725591376e-03;
        pointCoord(122, 0) = 5.778053775593315e-01;
        pointCoord(122, 1) = 5.310503885540051e-02;
        pointCoord(122, 2) = 3.517840436508153e-03;
        pointCoord(123, 0) = 5.592185728111451e-01;
        pointCoord(123, 1) = 9.641266013234269e-02;
        pointCoord(123, 2) = 4.367576624821128e-03;
        pointCoord(124, 0) = 5.361788354458206e-01;
        pointCoord(124, 1) = 1.500957021428338e-01;
        pointCoord(124, 2) = 4.886637244494526e-03;
        pointCoord(125, 0) = 5.096293190303758e-01;
        pointCoord(125, 1) = 2.119565984928767e-01;
        pointCoord(125, 2) = 5.063017983038727e-03;
        pointCoord(126, 0) = 4.806569299311483e-01;
        pointCoord(126, 1) = 2.794628359337138e-01;
        pointCoord(126, 2) = 4.921881796140457e-03;
        pointCoord(127, 0) = 4.504477852490379e-01;
        pointCoord(127, 1) = 3.498507382503206e-01;
        pointCoord(127, 2) = 4.519263907476252e-03;
        pointCoord(128, 0) = 4.202386405669276e-01;
        pointCoord(128, 1) = 4.202386405669276e-01;
        pointCoord(128, 2) = 3.931628247403513e-03;
        pointCoord(129, 0) = 3.912662514677002e-01;
        pointCoord(129, 1) = 4.877448780077646e-01;
        pointCoord(129, 2) = 3.243014099707468e-03;
        pointCoord(130, 0) = 3.647167350522554e-01;
        pointCoord(130, 1) = 5.496057743578074e-01;
        pointCoord(130, 2) = 2.531873387556408e-03;
        pointCoord(131, 0) = 3.416769976869309e-01;
        pointCoord(131, 1) = 6.032888163682986e-01;
        pointCoord(131, 2) = 1.859713347667004e-03;
        pointCoord(132, 0) = 3.230901929387445e-01;
        pointCoord(132, 1) = 6.465964376452408e-01;
        pointCoord(132, 2) = 1.263332435833526e-03;
        pointCoord(133, 0) = 3.097168773224545e-01;
        pointCoord(133, 1) = 6.777565265240901e-01;
        pointCoord(133, 2) = 7.518400035055260e-04;
        pointCoord(134, 0) = 3.021014380540554e-01;
        pointCoord(134, 1) = 6.955006500655807e-01;
        pointCoord(134, 2) = 3.092448083359775e-04;
        pointCoord(135, 0) = 6.949831426416880e-01;
        pointCoord(135, 1) = 3.911210092663540e-03;
        pointCoord(135, 2) = 9.281175941868199e-04;
        pointCoord(136, 0) = 6.861443754785116e-01;
        pointCoord(136, 1) = 2.043200573937687e-02;
        pointCoord(136, 2) = 2.082087995813234e-03;
        pointCoord(137, 0) = 6.706227990117994e-01;
        pointCoord(137, 1) = 4.944383376780406e-02;
        pointCoord(137, 2) = 3.059726443647184e-03;
        pointCoord(138, 0) = 6.490502496569180e-01;
        pointCoord(138, 1) = 8.976571043805123e-02;
        pointCoord(138, 2) = 3.784151182118204e-03;
        pointCoord(139, 0) = 6.223094581023360e-01;
        pointCoord(139, 1) = 1.397476982592848e-01;
        pointCoord(139, 2) = 4.211174391395614e-03;
        pointCoord(140, 0) = 5.914950840891640e-01;
        pointCoord(140, 1) = 1.973437370116005e-01;
        pointCoord(140, 2) = 4.331989740849268e-03;
        pointCoord(141, 0) = 5.578686323004080e-01;
        pointCoord(141, 1) = 2.601959117629090e-01;
        pointCoord(141, 2) = 4.172032952728379e-03;
        pointCoord(142, 0) = 5.228067551540863e-01;
        pointCoord(142, 1) = 3.257310816153046e-01;
        pointCoord(142, 2) = 3.785002146210653e-03;
        pointCoord(143, 0) = 4.877448780077646e-01;
        pointCoord(143, 1) = 3.912662514677002e-01;
        pointCoord(143, 2) = 3.243014099707468e-03;
        pointCoord(144, 0) = 4.541184262190086e-01;
        pointCoord(144, 1) = 4.541184262190086e-01;
        pointCoord(144, 2) = 2.624530148967344e-03;
        pointCoord(145, 0) = 4.233040522058366e-01;
        pointCoord(145, 1) = 5.117144649713243e-01;
        pointCoord(145, 2) = 2.002022987396075e-03;
        pointCoord(146, 0) = 3.965632606512546e-01;
        pointCoord(146, 1) = 5.616964527925579e-01;
        pointCoord(146, 2) = 1.431367626586710e-03;
        pointCoord(147, 0) = 3.749907112963732e-01;
        pointCoord(147, 1) = 6.020183294628051e-01;
        pointCoord(147, 2) = 9.446313353441727e-04;
        pointCoord(148, 0) = 3.594691348296611e-01;
        pointCoord(148, 1) = 6.310301574912323e-01;
        pointCoord(148, 2) = 5.473714982692157e-04;
        pointCoord(149, 0) = 3.506303676664846e-01;
        pointCoord(149, 1) = 6.475509531379455e-01;
        pointCoord(149, 2) = 2.210787444900900e-04;
        pointCoord(150, 0) = 7.831281587976550e-01;
        pointCoord(150, 1) = 3.645813483142782e-03;
        pointCoord(150, 2) = 7.724369699259093e-04;
        pointCoord(151, 0) = 7.731683669266911e-01;
        pointCoord(151, 1) = 1.904558442206878e-02;
        pointCoord(151, 2) = 1.730314409284368e-03;
        pointCoord(152, 0) = 7.556781821233451e-01;
        pointCoord(152, 1) = 4.608880411386213e-02;
        pointCoord(152, 2) = 2.535882479276867e-03;
        pointCoord(153, 0) = 7.313695768920815e-01;
        pointCoord(153, 1) = 8.367462490772701e-02;
        pointCoord(153, 2) = 3.123344756519591e-03;
        pointCoord(154, 0) = 7.012372390409338e-01;
        pointCoord(154, 1) = 1.302650664323950e-01;
        pointCoord(154, 2) = 3.455681100376738e-03;
        pointCoord(155, 0) = 6.665146644850798e-01;
        pointCoord(155, 1) = 1.839529046420287e-01;
        pointCoord(155, 2) = 3.527036015275618e-03;
        pointCoord(156, 0) = 6.286233550985966e-01;
        pointCoord(156, 1) = 2.425402217956099e-01;
        pointCoord(156, 2) = 3.361624124111803e-03;
        pointCoord(157, 0) = 5.891145647282021e-01;
        pointCoord(157, 1) = 3.036284784239326e-01;
        pointCoord(157, 2) = 3.008328952278663e-03;
        pointCoord(158, 0) = 5.496057743578074e-01;
        pointCoord(158, 1) = 3.647167350522554e-01;
        pointCoord(158, 2) = 2.531873387556408e-03;
        pointCoord(159, 0) = 5.117144649713243e-01;
        pointCoord(159, 1) = 4.233040522058366e-01;
        pointCoord(159, 2) = 2.002022987396075e-03;
        pointCoord(160, 0) = 4.769918904154703e-01;
        pointCoord(160, 1) = 4.769918904154703e-01;
        pointCoord(160, 2) = 1.482583354652366e-03;
        pointCoord(161, 0) = 4.468595525643226e-01;
        pointCoord(161, 1) = 5.235823319401383e-01;
        pointCoord(161, 2) = 1.021962333752187e-03;
        pointCoord(162, 0) = 4.225509473330590e-01;
        pointCoord(162, 1) = 5.611681527340032e-01;
        pointCoord(162, 2) = 6.467908658049148e-04;
        pointCoord(163, 0) = 4.050607625297130e-01;
        pointCoord(163, 1) = 5.882113724257966e-01;
        pointCoord(163, 2) = 3.595864201805764e-04;
        pointCoord(164, 0) = 3.951009706587492e-01;
        pointCoord(164, 1) = 6.036111433647224e-01;
        pointCoord(164, 2) = 1.409470874531633e-04;
        pointCoord(165, 0) = 8.596206263257772e-01;
        pointCoord(165, 1) = 3.415501635449438e-03;
        pointCoord(165, 2) = 6.072394228417237e-04;
        pointCoord(166, 0) = 8.486880063835482e-01;
        pointCoord(166, 1) = 1.784244450311027e-02;
        pointCoord(166, 2) = 1.358283575880347e-03;
        pointCoord(167, 0) = 8.294894582962350e-01;
        pointCoord(167, 1) = 4.317730091093636e-02;
        pointCoord(167, 2) = 1.985249465633941e-03;
        pointCoord(168, 0) = 8.028064968687120e-01;
        pointCoord(168, 1) = 7.838876550854190e-02;
        pointCoord(168, 2) = 2.434996767095785e-03;
        pointCoord(169, 0) = 7.697309665799893e-01;
        pointCoord(169, 1) = 1.220360145956348e-01;
        pointCoord(169, 2) = 2.678231687010688e-03;
        pointCoord(170, 0) = 7.316168457275697e-01;
        pointCoord(170, 1) = 1.723323065087029e-01;
        pointCoord(170, 2) = 2.711503061777502e-03;
        pointCoord(171, 0) = 6.900244821518289e-01;
        pointCoord(171, 1) = 2.272185694730266e-01;
        pointCoord(171, 2) = 2.556227607428541e-03;
        pointCoord(172, 0) = 6.466566492600638e-01;
        pointCoord(172, 1) = 2.844477835799787e-01;
        pointCoord(172, 2) = 2.254111925976988e-03;
        pointCoord(173, 0) = 6.032888163682986e-01;
        pointCoord(173, 1) = 3.416769976869309e-01;
        pointCoord(173, 2) = 1.859713347667004e-03;
        pointCoord(174, 0) = 5.616964527925579e-01;
        pointCoord(174, 1) = 3.965632606512546e-01;
        pointCoord(174, 2) = 1.431367626586710e-03;
        pointCoord(175, 0) = 5.235823319401383e-01;
        pointCoord(175, 1) = 4.468595525643226e-01;
        pointCoord(175, 2) = 1.021962333752187e-03;
        pointCoord(176, 0) = 4.905068016514156e-01;
        pointCoord(176, 1) = 4.905068016514156e-01;
        pointCoord(176, 2) = 6.710419332069080e-04;
        pointCoord(177, 0) = 4.638238402238926e-01;
        pointCoord(177, 1) = 5.257182662490212e-01;
        pointCoord(177, 2) = 3.994970356434159e-04;
        pointCoord(178, 0) = 4.446252921365793e-01;
        pointCoord(178, 1) = 5.510531226568472e-01;
        pointCoord(178, 2) = 2.076589983976010e-04;
        pointCoord(179, 0) = 4.336926721943504e-01;
        pointCoord(179, 1) = 5.654800655245080e-01;
        pointCoord(179, 2) = 7.715045918095962e-05;
        pointCoord(180, 0) = 9.213292532996236e-01;
        pointCoord(180, 1) = 3.229702584167157e-03;
        pointCoord(180, 2) = 4.407281567922316e-04;
        pointCoord(181, 0) = 9.096118255652115e-01;
        pointCoord(181, 1) = 1.687183765964483e-02;
        pointCoord(181, 2) = 9.845211473614196e-04;
        pointCoord(182, 0) = 8.890350927228124e-01;
        pointCoord(182, 1) = 4.082850931238514e-02;
        pointCoord(182, 2) = 1.435387316391502e-03;
        pointCoord(183, 0) = 8.604366713087993e-01;
        pointCoord(183, 1) = 7.412451392350064e-02;
        pointCoord(183, 2) = 1.753821796026420e-03;
        pointCoord(184, 0) = 8.249867848238173e-01;
        pointCoord(184, 1) = 1.153974068143342e-01;
        pointCoord(184, 2) = 1.918437141484353e-03;
        pointCoord(185, 0) = 7.841366080950154e-01;
        pointCoord(185, 1) = 1.629576428510200e-01;
        pointCoord(185, 2) = 1.927490467117033e-03;
        pointCoord(186, 0) = 7.395584999125796e-01;
        pointCoord(186, 1) = 2.148581612086487e-01;
        pointCoord(186, 2) = 1.798100384440941e-03;
        pointCoord(187, 0) = 6.930774687789103e-01;
        pointCoord(187, 1) = 2.689741770736966e-01;
        pointCoord(187, 2) = 1.562704823486195e-03;
        pointCoord(188, 0) = 6.465964376452408e-01;
        pointCoord(188, 1) = 3.230901929387445e-01;
        pointCoord(188, 2) = 1.263332435833526e-03;
        pointCoord(189, 0) = 6.020183294628051e-01;
        pointCoord(189, 1) = 3.749907112963732e-01;
        pointCoord(189, 2) = 9.446313353441727e-04;
        pointCoord(190, 0) = 5.611681527340032e-01;
        pointCoord(190, 1) = 4.225509473330590e-01;
        pointCoord(190, 2) = 6.467908658049148e-04;
        pointCoord(191, 0) = 5.257182662490212e-01;
        pointCoord(191, 1) = 4.638238402238926e-01;
        pointCoord(191, 2) = 3.994970356434159e-04;
        pointCoord(192, 0) = 4.971198448350080e-01;
        pointCoord(192, 1) = 4.971198448350080e-01;
        pointCoord(192, 2) = 2.178823448844040e-04;
        pointCoord(193, 0) = 4.765431119926091e-01;
        pointCoord(193, 1) = 5.210765164877483e-01;
        pointCoord(193, 2) = 1.010975328493697e-04;
        pointCoord(194, 0) = 4.648256842581970e-01;
        pointCoord(194, 1) = 5.347186515632260e-01;
        pointCoord(194, 2) = 3.373779915999911e-05;
        pointCoord(195, 0) = 9.657289742817724e-01;
        pointCoord(195, 1) = 3.096019071176703e-03;
        pointCoord(195, 2) = 2.773577221348077e-04;
        pointCoord(196, 0) = 9.534468727131440e-01;
        pointCoord(196, 1) = 1.617348031243802e-02;
        pointCoord(196, 2) = 6.189225798873323e-04;
        pointCoord(197, 0) = 9.318785278127811e-01;
        pointCoord(197, 1) = 3.913853990721448e-02;
        pointCoord(197, 2) = 9.005715516839962e-04;
        pointCoord(198, 0) = 9.019019216436822e-01;
        pointCoord(198, 1) = 7.105635976324422e-02;
        pointCoord(198, 2) = 1.096975661844450e-03;
        pointCoord(199, 0) = 8.647436718747329e-01;
        pointCoord(199, 1) = 1.106208893701105e-01;
        pointCoord(199, 2) = 1.194612298797779e-03;
        pointCoord(200, 0) = 8.219248868092972e-01;
        pointCoord(200, 1) = 1.562125170701625e-01;
        pointCoord(200, 2) = 1.192765446068943e-03;
        pointCoord(201, 0) = 7.751985177764392e-01;
        pointCoord(201, 1) = 2.059647745773691e-01;
        pointCoord(201, 2) = 1.102995104324453e-03;
        pointCoord(202, 0) = 7.264775221502647e-01;
        pointCoord(202, 1) = 2.578408259499118e-01;
        pointCoord(202, 2) = 9.467984241174950e-04;
        pointCoord(203, 0) = 6.777565265240901e-01;
        pointCoord(203, 1) = 3.097168773224545e-01;
        pointCoord(203, 2) = 7.518400035055260e-04;
        pointCoord(204, 0) = 6.310301574912323e-01;
        pointCoord(204, 1) = 3.594691348296611e-01;
        pointCoord(204, 2) = 5.473714982692157e-04;
        pointCoord(205, 0) = 5.882113724257966e-01;
        pointCoord(205, 1) = 4.050607625297130e-01;
        pointCoord(205, 2) = 3.595864201805764e-04;
        pointCoord(206, 0) = 5.510531226568472e-01;
        pointCoord(206, 1) = 4.446252921365793e-01;
        pointCoord(206, 2) = 2.076589983976010e-04;
        pointCoord(207, 0) = 5.210765164877483e-01;
        pointCoord(207, 1) = 4.765431119926091e-01;
        pointCoord(207, 2) = 1.010975328493697e-04;
        pointCoord(208, 0) = 4.995081715873855e-01;
        pointCoord(208, 1) = 4.995081715873855e-01;
        pointCoord(208, 2) = 3.882291380293542e-05;
        pointCoord(209, 0) = 4.872260700187571e-01;
        pointCoord(209, 1) = 5.125856328286468e-01;
        pointCoord(209, 2) = 1.010771128646518e-05;
        pointCoord(210, 0) = 9.910124109683001e-01;
        pointCoord(210, 1) = 3.019892947814674e-03;
        pointCoord(210, 2) = 1.182202366589500e-04;
        pointCoord(211, 0) = 9.784087556867173e-01;
        pointCoord(211, 1) = 1.577580047612162e-02;
        pointCoord(211, 2) = 2.636382995431647e-04;
        pointCoord(212, 0) = 9.562757369521381e-01;
        pointCoord(212, 1) = 3.817618623668173e-02;
        pointCoord(212, 2) = 3.831452029932598e-04;
        pointCoord(213, 0) = 9.255143229908560e-01;
        pointCoord(213, 1) = 6.930919830052830e-02;
        pointCoord(213, 2) = 4.658232856144815e-04;
        pointCoord(214, 0) = 8.873832451506416e-01;
        pointCoord(214, 1) = 1.079008998361295e-01;
        pointCoord(214, 2) = 5.058923108220498e-04;
        pointCoord(215, 0) = 8.434434353774185e-01;
        pointCoord(215, 1) = 1.523715028283941e-01;
        pointCoord(215, 2) = 5.031459767425094e-04;
        pointCoord(216, 0) = 7.954937384497832e-01;
        pointCoord(216, 1) = 2.009004324408233e-01;
        pointCoord(216, 2) = 4.627159521497556e-04;
        pointCoord(217, 0) = 7.454971942576820e-01;
        pointCoord(217, 1) = 2.515009352474393e-01;
        pointCoord(217, 2) = 3.940464726072197e-04;
        pointCoord(218, 0) = 6.955006500655807e-01;
        pointCoord(218, 1) = 3.021014380540554e-01;
        pointCoord(218, 2) = 3.092448083359775e-04;
        pointCoord(219, 0) = 6.475509531379455e-01;
        pointCoord(219, 1) = 3.506303676664846e-01;
        pointCoord(219, 2) = 2.210787444900900e-04;
        pointCoord(220, 0) = 6.036111433647224e-01;
        pointCoord(220, 1) = 3.951009706587492e-01;
        pointCoord(220, 2) = 1.409470874531633e-04;
        pointCoord(221, 0) = 5.654800655245080e-01;
        pointCoord(221, 1) = 4.336926721943504e-01;
        pointCoord(221, 2) = 7.715045918095962e-05;
        pointCoord(222, 0) = 5.347186515632260e-01;
        pointCoord(222, 1) = 4.648256842581970e-01;
        pointCoord(222, 2) = 3.373779915999911e-05;
        pointCoord(223, 0) = 5.125856328286468e-01;
        pointCoord(223, 1) = 4.872260700187571e-01;
        pointCoord(223, 2) = 1.010771128646518e-05;
        pointCoord(224, 0) = 4.999819775470640e-01;
        pointCoord(224, 1) = 4.999819775470640e-01;
        pointCoord(224, 2) = 1.419527361296296e-06;
    }
    else if (nH == 400)
    {
      pointCoord(0, 0) = 5.501388317558e-04; pointCoord(0, 1) =1.117738042591e-02; pointCoord(0, 2) =4.114479108130e-05;
pointCoord(1, 0) = 2.505501388318e-01; pointCoord(1, 1) =1.117738042591e-02; pointCoord(1, 2) =4.114479108130e-05;
pointCoord(2, 0) = 5.005501388318e-01; pointCoord(2, 1) =1.117738042591e-02; pointCoord(2, 2) =4.114479108130e-05;
pointCoord(3, 0) = 7.505501388318e-01; pointCoord(3, 1) =1.117738042591e-02; pointCoord(3, 2) =4.114479108130e-05;
pointCoord(4, 0) = 5.501388317558e-04; pointCoord(4, 1) =2.611773804259e-01; pointCoord(4, 2) =4.114479108130e-05;
pointCoord(5, 0) = 2.505501388318e-01; pointCoord(5, 1) =2.611773804259e-01; pointCoord(5, 2) =4.114479108130e-05;
pointCoord(6, 0) = 5.005501388318e-01; pointCoord(6, 1) =2.611773804259e-01; pointCoord(6, 2) =4.114479108130e-05;
pointCoord(7, 0) = 5.501388317558e-04; pointCoord(7, 1) =5.111773804259e-01; pointCoord(7, 2) =4.114479108130e-05;
pointCoord(8, 0) = 2.505501388318e-01; pointCoord(8, 1) =5.111773804259e-01; pointCoord(8, 2) =4.114479108130e-05;
pointCoord(9, 0) = 5.501388317558e-04; pointCoord(9, 1) =7.611773804259e-01; pointCoord(9, 2) =4.114479108130e-05;
pointCoord(10, 0) = 2.706305026870e-03; pointCoord(10, 1) =9.021214230797e-03; pointCoord(10, 2) =8.311879273871e-05;
pointCoord(11, 0) = 2.527063050269e-01; pointCoord(11, 1) =9.021214230797e-03; pointCoord(11, 2) =8.311879273871e-05;
pointCoord(12, 0) = 5.027063050269e-01; pointCoord(12, 1) =9.021214230797e-03; pointCoord(12, 2) =8.311879273871e-05;
pointCoord(13, 0) = 7.527063050269e-01; pointCoord(13, 1) =9.021214230797e-03; pointCoord(13, 2) =8.311879273871e-05;
pointCoord(14, 0) = 2.706305026870e-03; pointCoord(14, 1) =2.590212142308e-01; pointCoord(14, 2) =8.311879273871e-05;
pointCoord(15, 0) = 2.527063050269e-01; pointCoord(15, 1) =2.590212142308e-01; pointCoord(15, 2) =8.311879273871e-05;
pointCoord(16, 0) = 5.027063050269e-01; pointCoord(16, 1) =2.590212142308e-01; pointCoord(16, 2) =8.311879273871e-05;
pointCoord(17, 0) = 2.706305026870e-03; pointCoord(17, 1) =5.090212142308e-01; pointCoord(17, 2) =8.311879273871e-05;
pointCoord(18, 0) = 2.527063050269e-01; pointCoord(18, 1) =5.090212142308e-01; pointCoord(18, 2) =8.311879273871e-05;
pointCoord(19, 0) = 2.706305026870e-03; pointCoord(19, 1) =7.590212142308e-01; pointCoord(19, 2) =8.311879273871e-05;
pointCoord(20, 0) = 5.863759628834e-03; pointCoord(20, 1) =5.863759628834e-03; pointCoord(20, 2) =9.879340825442e-05;
pointCoord(21, 0) = 2.558637596288e-01; pointCoord(21, 1) =5.863759628834e-03; pointCoord(21, 2) =9.879340825442e-05;
pointCoord(22, 0) = 5.058637596288e-01; pointCoord(22, 1) =5.863759628834e-03; pointCoord(22, 2) =9.879340825442e-05;
pointCoord(23, 0) = 7.558637596288e-01; pointCoord(23, 1) =5.863759628834e-03; pointCoord(23, 2) =9.879340825442e-05;
pointCoord(24, 0) = 5.863759628834e-03; pointCoord(24, 1) =2.558637596288e-01; pointCoord(24, 2) =9.879340825442e-05;
pointCoord(25, 0) = 2.558637596288e-01; pointCoord(25, 1) =2.558637596288e-01; pointCoord(25, 2) =9.879340825442e-05;
pointCoord(26, 0) = 5.058637596288e-01; pointCoord(26, 1) =2.558637596288e-01; pointCoord(26, 2) =9.879340825442e-05;
pointCoord(27, 0) = 5.863759628834e-03; pointCoord(27, 1) =5.058637596288e-01; pointCoord(27, 2) =9.879340825442e-05;
pointCoord(28, 0) = 2.558637596288e-01; pointCoord(28, 1) =5.058637596288e-01; pointCoord(28, 2) =9.879340825442e-05;
pointCoord(29, 0) = 5.863759628834e-03; pointCoord(29, 1) =7.558637596288e-01; pointCoord(29, 2) =9.879340825442e-05;
pointCoord(30, 0) = 9.021214230797e-03; pointCoord(30, 1) =2.706305026870e-03; pointCoord(30, 2) =8.311879273871e-05;
pointCoord(31, 0) = 2.590212142308e-01; pointCoord(31, 1) =2.706305026870e-03; pointCoord(31, 2) =8.311879273871e-05;
pointCoord(32, 0) = 5.090212142308e-01; pointCoord(32, 1) =2.706305026870e-03; pointCoord(32, 2) =8.311879273871e-05;
pointCoord(33, 0) = 7.590212142308e-01; pointCoord(33, 1) =2.706305026870e-03; pointCoord(33, 2) =8.311879273871e-05;
pointCoord(34, 0) = 9.021214230797e-03; pointCoord(34, 1) =2.527063050269e-01; pointCoord(34, 2) =8.311879273871e-05;
pointCoord(35, 0) = 2.590212142308e-01; pointCoord(35, 1) =2.527063050269e-01; pointCoord(35, 2) =8.311879273871e-05;
pointCoord(36, 0) = 5.090212142308e-01; pointCoord(36, 1) =2.527063050269e-01; pointCoord(36, 2) =8.311879273871e-05;
pointCoord(37, 0) = 9.021214230797e-03; pointCoord(37, 1) =5.027063050269e-01; pointCoord(37, 2) =8.311879273871e-05;
pointCoord(38, 0) = 2.590212142308e-01; pointCoord(38, 1) =5.027063050269e-01; pointCoord(38, 2) =8.311879273871e-05;
pointCoord(39, 0) = 9.021214230797e-03; pointCoord(39, 1) =7.527063050269e-01; pointCoord(39, 2) =8.311879273871e-05;
pointCoord(40, 0) = 1.117738042591e-02; pointCoord(40, 1) =5.501388317558e-04; pointCoord(40, 2) =4.114479108130e-05;
pointCoord(41, 0) = 2.611773804259e-01; pointCoord(41, 1) =5.501388317558e-04; pointCoord(41, 2) =4.114479108130e-05;
pointCoord(42, 0) = 5.111773804259e-01; pointCoord(42, 1) =5.501388317558e-04; pointCoord(42, 2) =4.114479108130e-05;
pointCoord(43, 0) = 7.611773804259e-01; pointCoord(43, 1) =5.501388317558e-04; pointCoord(43, 2) =4.114479108130e-05;
pointCoord(44, 0) = 1.117738042591e-02; pointCoord(44, 1) =2.505501388318e-01; pointCoord(44, 2) =4.114479108130e-05;
pointCoord(45, 0) = 2.611773804259e-01; pointCoord(45, 1) =2.505501388318e-01; pointCoord(45, 2) =4.114479108130e-05;
pointCoord(46, 0) = 5.111773804259e-01; pointCoord(46, 1) =2.505501388318e-01; pointCoord(46, 2) =4.114479108130e-05;
pointCoord(47, 0) = 1.117738042591e-02; pointCoord(47, 1) =5.005501388318e-01; pointCoord(47, 2) =4.114479108130e-05;
pointCoord(48, 0) = 2.611773804259e-01; pointCoord(48, 1) =5.005501388318e-01; pointCoord(48, 2) =4.114479108130e-05;
pointCoord(49, 0) = 1.117738042591e-02; pointCoord(49, 1) =7.505501388318e-01; pointCoord(49, 2) =4.114479108130e-05;
pointCoord(50, 0) = 2.706305026870e-03; pointCoord(50, 1) =5.498503120992e-02; pointCoord(50, 2) =4.088873455782e-04;
pointCoord(51, 0) = 2.527063050269e-01; pointCoord(51, 1) =5.498503120992e-02; pointCoord(51, 2) =4.088873455782e-04;
pointCoord(52, 0) = 5.027063050269e-01; pointCoord(52, 1) =5.498503120992e-02; pointCoord(52, 2) =4.088873455782e-04;
pointCoord(53, 0) = 7.527063050269e-01; pointCoord(53, 1) =5.498503120992e-02; pointCoord(53, 2) =4.088873455782e-04;
pointCoord(54, 0) = 2.706305026870e-03; pointCoord(54, 1) =3.049850312099e-01; pointCoord(54, 2) =4.088873455782e-04;
pointCoord(55, 0) = 2.527063050269e-01; pointCoord(55, 1) =3.049850312099e-01; pointCoord(55, 2) =4.088873455782e-04;
pointCoord(56, 0) = 5.027063050269e-01; pointCoord(56, 1) =3.049850312099e-01; pointCoord(56, 2) =4.088873455782e-04;
pointCoord(57, 0) = 2.706305026870e-03; pointCoord(57, 1) =5.549850312099e-01; pointCoord(57, 2) =4.088873455782e-04;
pointCoord(58, 0) = 2.527063050269e-01; pointCoord(58, 1) =5.549850312099e-01; pointCoord(58, 2) =4.088873455782e-04;
pointCoord(59, 0) = 2.706305026870e-03; pointCoord(59, 1) =8.049850312099e-01; pointCoord(59, 2) =4.088873455782e-04;
pointCoord(60, 0) = 1.331316110715e-02; pointCoord(60, 1) =4.437817512964e-02; pointCoord(60, 2) =8.260151926267e-04;
pointCoord(61, 0) = 2.633131611071e-01; pointCoord(61, 1) =4.437817512964e-02; pointCoord(61, 2) =8.260151926267e-04;
pointCoord(62, 0) = 5.133131611071e-01; pointCoord(62, 1) =4.437817512964e-02; pointCoord(62, 2) =8.260151926267e-04;
pointCoord(63, 0) = 7.633131611071e-01; pointCoord(63, 1) =4.437817512964e-02; pointCoord(63, 2) =8.260151926267e-04;
pointCoord(64, 0) = 1.331316110715e-02; pointCoord(64, 1) =2.943781751296e-01; pointCoord(64, 2) =8.260151926267e-04;
pointCoord(65, 0) = 2.633131611071e-01; pointCoord(65, 1) =2.943781751296e-01; pointCoord(65, 2) =8.260151926267e-04;
pointCoord(66, 0) = 5.133131611071e-01; pointCoord(66, 1) =2.943781751296e-01; pointCoord(66, 2) =8.260151926267e-04;
pointCoord(67, 0) = 1.331316110715e-02; pointCoord(67, 1) =5.443781751296e-01; pointCoord(67, 2) =8.260151926267e-04;
pointCoord(68, 0) = 2.633131611071e-01; pointCoord(68, 1) =5.443781751296e-01; pointCoord(68, 2) =8.260151926267e-04;
pointCoord(69, 0) = 1.331316110715e-02; pointCoord(69, 1) =7.943781751296e-01; pointCoord(69, 2) =8.260151926267e-04;
pointCoord(70, 0) = 2.884566811839e-02; pointCoord(70, 1) =2.884566811839e-02; pointCoord(70, 2) =9.817858688834e-04;
pointCoord(71, 0) = 2.788456681184e-01; pointCoord(71, 1) =2.884566811839e-02; pointCoord(71, 2) =9.817858688834e-04;
pointCoord(72, 0) = 5.288456681184e-01; pointCoord(72, 1) =2.884566811839e-02; pointCoord(72, 2) =9.817858688834e-04;
pointCoord(73, 0) = 7.788456681184e-01; pointCoord(73, 1) =2.884566811839e-02; pointCoord(73, 2) =9.817858688834e-04;
pointCoord(74, 0) = 2.884566811839e-02; pointCoord(74, 1) =2.788456681184e-01; pointCoord(74, 2) =9.817858688834e-04;
pointCoord(75, 0) = 2.788456681184e-01; pointCoord(75, 1) =2.788456681184e-01; pointCoord(75, 2) =9.817858688834e-04;
pointCoord(76, 0) = 5.288456681184e-01; pointCoord(76, 1) =2.788456681184e-01; pointCoord(76, 2) =9.817858688834e-04;
pointCoord(77, 0) = 2.884566811839e-02; pointCoord(77, 1) =5.288456681184e-01; pointCoord(77, 2) =9.817858688834e-04;
pointCoord(78, 0) = 2.788456681184e-01; pointCoord(78, 1) =5.288456681184e-01; pointCoord(78, 2) =9.817858688834e-04;
pointCoord(79, 0) = 2.884566811839e-02; pointCoord(79, 1) =7.788456681184e-01; pointCoord(79, 2) =9.817858688834e-04;
pointCoord(80, 0) = 4.437817512964e-02; pointCoord(80, 1) =1.331316110715e-02; pointCoord(80, 2) =8.260151926267e-04;
pointCoord(81, 0) = 2.943781751296e-01; pointCoord(81, 1) =1.331316110715e-02; pointCoord(81, 2) =8.260151926267e-04;
pointCoord(82, 0) = 5.443781751296e-01; pointCoord(82, 1) =1.331316110715e-02; pointCoord(82, 2) =8.260151926267e-04;
pointCoord(83, 0) = 7.943781751296e-01; pointCoord(83, 1) =1.331316110715e-02; pointCoord(83, 2) =8.260151926267e-04;
pointCoord(84, 0) = 4.437817512964e-02; pointCoord(84, 1) =2.633131611071e-01; pointCoord(84, 2) =8.260151926267e-04;
pointCoord(85, 0) = 2.943781751296e-01; pointCoord(85, 1) =2.633131611071e-01; pointCoord(85, 2) =8.260151926267e-04;
pointCoord(86, 0) = 5.443781751296e-01; pointCoord(86, 1) =2.633131611071e-01; pointCoord(86, 2) =8.260151926267e-04;
pointCoord(87, 0) = 4.437817512964e-02; pointCoord(87, 1) =5.133131611071e-01; pointCoord(87, 2) =8.260151926267e-04;
pointCoord(88, 0) = 2.943781751296e-01; pointCoord(88, 1) =5.133131611071e-01; pointCoord(88, 2) =8.260151926267e-04;
pointCoord(89, 0) = 4.437817512964e-02; pointCoord(89, 1) =7.633131611071e-01; pointCoord(89, 2) =8.260151926267e-04;
pointCoord(90, 0) = 5.498503120992e-02; pointCoord(90, 1) =2.706305026870e-03; pointCoord(90, 2) =4.088873455782e-04;
pointCoord(91, 0) = 3.049850312099e-01; pointCoord(91, 1) =2.706305026870e-03; pointCoord(91, 2) =4.088873455782e-04;
pointCoord(92, 0) = 5.549850312099e-01; pointCoord(92, 1) =2.706305026870e-03; pointCoord(92, 2) =4.088873455782e-04;
pointCoord(93, 0) = 8.049850312099e-01; pointCoord(93, 1) =2.706305026870e-03; pointCoord(93, 2) =4.088873455782e-04;
pointCoord(94, 0) = 5.498503120992e-02; pointCoord(94, 1) =2.527063050269e-01; pointCoord(94, 2) =4.088873455782e-04;
pointCoord(95, 0) = 3.049850312099e-01; pointCoord(95, 1) =2.527063050269e-01; pointCoord(95, 2) =4.088873455782e-04;
pointCoord(96, 0) = 5.549850312099e-01; pointCoord(96, 1) =2.527063050269e-01; pointCoord(96, 2) =4.088873455782e-04;
pointCoord(97, 0) = 5.498503120992e-02; pointCoord(97, 1) =5.027063050269e-01; pointCoord(97, 2) =4.088873455782e-04;
pointCoord(98, 0) = 3.049850312099e-01; pointCoord(98, 1) =5.027063050269e-01; pointCoord(98, 2) =4.088873455782e-04;
pointCoord(99, 0) = 5.498503120992e-02; pointCoord(99, 1) =7.527063050269e-01; pointCoord(99, 2) =4.088873455782e-04;
pointCoord(100, 0) = 5.863759628834e-03; pointCoord(100, 1) =1.191362403712e-01; pointCoord(100, 2) =1.053008378028e-03;
pointCoord(101, 0) = 2.558637596288e-01; pointCoord(101, 1) =1.191362403712e-01; pointCoord(101, 2) =1.053008378028e-03;
pointCoord(102, 0) = 5.058637596288e-01; pointCoord(102, 1) =1.191362403712e-01; pointCoord(102, 2) =1.053008378028e-03;
pointCoord(103, 0) = 7.558637596288e-01; pointCoord(103, 1) =1.191362403712e-01; pointCoord(103, 2) =1.053008378028e-03;
pointCoord(104, 0) = 5.863759628834e-03; pointCoord(104, 1) =3.691362403712e-01; pointCoord(104, 2) =1.053008378028e-03;
pointCoord(105, 0) = 2.558637596288e-01; pointCoord(105, 1) =3.691362403712e-01; pointCoord(105, 2) =1.053008378028e-03;
pointCoord(106, 0) = 5.058637596288e-01; pointCoord(106, 1) =3.691362403712e-01; pointCoord(106, 2) =1.053008378028e-03;
pointCoord(107, 0) = 5.863759628834e-03; pointCoord(107, 1) =6.191362403712e-01; pointCoord(107, 2) =1.053008378028e-03;
pointCoord(108, 0) = 2.558637596288e-01; pointCoord(108, 1) =6.191362403712e-01; pointCoord(108, 2) =1.053008378028e-03;
pointCoord(109, 0) = 5.863759628834e-03; pointCoord(109, 1) =8.691362403712e-01; pointCoord(109, 2) =1.053008378028e-03;
pointCoord(110, 0) = 2.884566811839e-02; pointCoord(110, 1) =9.615433188161e-02; pointCoord(110, 2) =2.127238535553e-03;
pointCoord(111, 0) = 2.788456681184e-01; pointCoord(111, 1) =9.615433188161e-02; pointCoord(111, 2) =2.127238535553e-03;
pointCoord(112, 0) = 5.288456681184e-01; pointCoord(112, 1) =9.615433188161e-02; pointCoord(112, 2) =2.127238535553e-03;
pointCoord(113, 0) = 7.788456681184e-01; pointCoord(113, 1) =9.615433188161e-02; pointCoord(113, 2) =2.127238535553e-03;
pointCoord(114, 0) = 2.884566811839e-02; pointCoord(114, 1) =3.461543318816e-01; pointCoord(114, 2) =2.127238535553e-03;
pointCoord(115, 0) = 2.788456681184e-01; pointCoord(115, 1) =3.461543318816e-01; pointCoord(115, 2) =2.127238535553e-03;
pointCoord(116, 0) = 5.288456681184e-01; pointCoord(116, 1) =3.461543318816e-01; pointCoord(116, 2) =2.127238535553e-03;
pointCoord(117, 0) = 2.884566811839e-02; pointCoord(117, 1) =5.961543318816e-01; pointCoord(117, 2) =2.127238535553e-03;
pointCoord(118, 0) = 2.788456681184e-01; pointCoord(118, 1) =5.961543318816e-01; pointCoord(118, 2) =2.127238535553e-03;
pointCoord(119, 0) = 2.884566811839e-02; pointCoord(119, 1) =8.461543318816e-01; pointCoord(119, 2) =2.127238535553e-03;
pointCoord(120, 0) = 6.250000000000e-02; pointCoord(120, 1) =6.250000000000e-02; pointCoord(120, 2) =2.528395061728e-03;
pointCoord(121, 0) = 3.125000000000e-01; pointCoord(121, 1) =6.250000000000e-02; pointCoord(121, 2) =2.528395061728e-03;
pointCoord(122, 0) = 5.625000000000e-01; pointCoord(122, 1) =6.250000000000e-02; pointCoord(122, 2) =2.528395061728e-03;
pointCoord(123, 0) = 8.125000000000e-01; pointCoord(123, 1) =6.250000000000e-02; pointCoord(123, 2) =2.528395061728e-03;
pointCoord(124, 0) = 6.250000000000e-02; pointCoord(124, 1) =3.125000000000e-01; pointCoord(124, 2) =2.528395061728e-03;
pointCoord(125, 0) = 3.125000000000e-01; pointCoord(125, 1) =3.125000000000e-01; pointCoord(125, 2) =2.528395061728e-03;
pointCoord(126, 0) = 5.625000000000e-01; pointCoord(126, 1) =3.125000000000e-01; pointCoord(126, 2) =2.528395061728e-03;
pointCoord(127, 0) = 6.250000000000e-02; pointCoord(127, 1) =5.625000000000e-01; pointCoord(127, 2) =2.528395061728e-03;
pointCoord(128, 0) = 3.125000000000e-01; pointCoord(128, 1) =5.625000000000e-01; pointCoord(128, 2) =2.528395061728e-03;
pointCoord(129, 0) = 6.250000000000e-02; pointCoord(129, 1) =8.125000000000e-01; pointCoord(129, 2) =2.528395061728e-03;
pointCoord(130, 0) = 9.615433188161e-02; pointCoord(130, 1) =2.884566811839e-02; pointCoord(130, 2) =2.127238535553e-03;
pointCoord(131, 0) = 3.461543318816e-01; pointCoord(131, 1) =2.884566811839e-02; pointCoord(131, 2) =2.127238535553e-03;
pointCoord(132, 0) = 5.961543318816e-01; pointCoord(132, 1) =2.884566811839e-02; pointCoord(132, 2) =2.127238535553e-03;
pointCoord(133, 0) = 8.461543318816e-01; pointCoord(133, 1) =2.884566811839e-02; pointCoord(133, 2) =2.127238535553e-03;
pointCoord(134, 0) = 9.615433188161e-02; pointCoord(134, 1) =2.788456681184e-01; pointCoord(134, 2) =2.127238535553e-03;
pointCoord(135, 0) = 3.461543318816e-01; pointCoord(135, 1) =2.788456681184e-01; pointCoord(135, 2) =2.127238535553e-03;
pointCoord(136, 0) = 5.961543318816e-01; pointCoord(136, 1) =2.788456681184e-01; pointCoord(136, 2) =2.127238535553e-03;
pointCoord(137, 0) = 9.615433188161e-02; pointCoord(137, 1) =5.288456681184e-01; pointCoord(137, 2) =2.127238535553e-03;
pointCoord(138, 0) = 3.461543318816e-01; pointCoord(138, 1) =5.288456681184e-01; pointCoord(138, 2) =2.127238535553e-03;
pointCoord(139, 0) = 9.615433188161e-02; pointCoord(139, 1) =7.788456681184e-01; pointCoord(139, 2) =2.127238535553e-03;
pointCoord(140, 0) = 1.191362403712e-01; pointCoord(140, 1) =5.863759628834e-03; pointCoord(140, 2) =1.053008378028e-03;
pointCoord(141, 0) = 3.691362403712e-01; pointCoord(141, 1) =5.863759628834e-03; pointCoord(141, 2) =1.053008378028e-03;
pointCoord(142, 0) = 6.191362403712e-01; pointCoord(142, 1) =5.863759628834e-03; pointCoord(142, 2) =1.053008378028e-03;
pointCoord(143, 0) = 8.691362403712e-01; pointCoord(143, 1) =5.863759628834e-03; pointCoord(143, 2) =1.053008378028e-03;
pointCoord(144, 0) = 1.191362403712e-01; pointCoord(144, 1) =2.558637596288e-01; pointCoord(144, 2) =1.053008378028e-03;
pointCoord(145, 0) = 3.691362403712e-01; pointCoord(145, 1) =2.558637596288e-01; pointCoord(145, 2) =1.053008378028e-03;
pointCoord(146, 0) = 6.191362403712e-01; pointCoord(146, 1) =2.558637596288e-01; pointCoord(146, 2) =1.053008378028e-03;
pointCoord(147, 0) = 1.191362403712e-01; pointCoord(147, 1) =5.058637596288e-01; pointCoord(147, 2) =1.053008378028e-03;
pointCoord(148, 0) = 3.691362403712e-01; pointCoord(148, 1) =5.058637596288e-01; pointCoord(148, 2) =1.053008378028e-03;
pointCoord(149, 0) = 1.191362403712e-01; pointCoord(149, 1) =7.558637596288e-01; pointCoord(149, 2) =1.053008378028e-03;
pointCoord(150, 0) = 9.021214230797e-03; pointCoord(150, 1) =1.832874495324e-01; pointCoord(150, 2) =1.362987654422e-03;
pointCoord(151, 0) = 2.590212142308e-01; pointCoord(151, 1) =1.832874495324e-01; pointCoord(151, 2) =1.362987654422e-03;
pointCoord(152, 0) = 5.090212142308e-01; pointCoord(152, 1) =1.832874495324e-01; pointCoord(152, 2) =1.362987654422e-03;
pointCoord(153, 0) = 7.590212142308e-01; pointCoord(153, 1) =1.832874495324e-01; pointCoord(153, 2) =1.362987654422e-03;
pointCoord(154, 0) = 9.021214230797e-03; pointCoord(154, 1) =4.332874495324e-01; pointCoord(154, 2) =1.362987654422e-03;
pointCoord(155, 0) = 2.590212142308e-01; pointCoord(155, 1) =4.332874495324e-01; pointCoord(155, 2) =1.362987654422e-03;
pointCoord(156, 0) = 5.090212142308e-01; pointCoord(156, 1) =4.332874495324e-01; pointCoord(156, 2) =1.362987654422e-03;
pointCoord(157, 0) = 9.021214230797e-03; pointCoord(157, 1) =6.832874495324e-01; pointCoord(157, 2) =1.362987654422e-03;
pointCoord(158, 0) = 2.590212142308e-01; pointCoord(158, 1) =6.832874495324e-01; pointCoord(158, 2) =1.362987654422e-03;
pointCoord(159, 0) = 9.021214230797e-03; pointCoord(159, 1) =9.332874495324e-01; pointCoord(159, 2) =1.362987654422e-03;
pointCoord(160, 0) = 4.437817512964e-02; pointCoord(160, 1) =1.479304886336e-01; pointCoord(160, 2) =2.753444248373e-03;
pointCoord(161, 0) = 2.943781751296e-01; pointCoord(161, 1) =1.479304886336e-01; pointCoord(161, 2) =2.753444248373e-03;
pointCoord(162, 0) = 5.443781751296e-01; pointCoord(162, 1) =1.479304886336e-01; pointCoord(162, 2) =2.753444248373e-03;
pointCoord(163, 0) = 7.943781751296e-01; pointCoord(163, 1) =1.479304886336e-01; pointCoord(163, 2) =2.753444248373e-03;
pointCoord(164, 0) = 4.437817512964e-02; pointCoord(164, 1) =3.979304886336e-01; pointCoord(164, 2) =2.753444248373e-03;
pointCoord(165, 0) = 2.943781751296e-01; pointCoord(165, 1) =3.979304886336e-01; pointCoord(165, 2) =2.753444248373e-03;
pointCoord(166, 0) = 5.443781751296e-01; pointCoord(166, 1) =3.979304886336e-01; pointCoord(166, 2) =2.753444248373e-03;
pointCoord(167, 0) = 4.437817512964e-02; pointCoord(167, 1) =6.479304886336e-01; pointCoord(167, 2) =2.753444248373e-03;
pointCoord(168, 0) = 2.943781751296e-01; pointCoord(168, 1) =6.479304886336e-01; pointCoord(168, 2) =2.753444248373e-03;
pointCoord(169, 0) = 4.437817512964e-02; pointCoord(169, 1) =8.979304886336e-01; pointCoord(169, 2) =2.753444248373e-03;
pointCoord(170, 0) = 9.615433188161e-02; pointCoord(170, 1) =9.615433188161e-02; pointCoord(170, 2) =3.272691202222e-03;
pointCoord(171, 0) = 3.461543318816e-01; pointCoord(171, 1) =9.615433188161e-02; pointCoord(171, 2) =3.272691202222e-03;
pointCoord(172, 0) = 5.961543318816e-01; pointCoord(172, 1) =9.615433188161e-02; pointCoord(172, 2) =3.272691202222e-03;
pointCoord(173, 0) = 8.461543318816e-01; pointCoord(173, 1) =9.615433188161e-02; pointCoord(173, 2) =3.272691202222e-03;
pointCoord(174, 0) = 9.615433188161e-02; pointCoord(174, 1) =3.461543318816e-01; pointCoord(174, 2) =3.272691202222e-03;
pointCoord(175, 0) = 3.461543318816e-01; pointCoord(175, 1) =3.461543318816e-01; pointCoord(175, 2) =3.272691202222e-03;
pointCoord(176, 0) = 5.961543318816e-01; pointCoord(176, 1) =3.461543318816e-01; pointCoord(176, 2) =3.272691202222e-03;
pointCoord(177, 0) = 9.615433188161e-02; pointCoord(177, 1) =5.961543318816e-01; pointCoord(177, 2) =3.272691202222e-03;
pointCoord(178, 0) = 3.461543318816e-01; pointCoord(178, 1) =5.961543318816e-01; pointCoord(178, 2) =3.272691202222e-03;
pointCoord(179, 0) = 9.615433188161e-02; pointCoord(179, 1) =8.461543318816e-01; pointCoord(179, 2) =3.272691202222e-03;
pointCoord(180, 0) = 1.479304886336e-01; pointCoord(180, 1) =4.437817512964e-02; pointCoord(180, 2) =2.753444248373e-03;
pointCoord(181, 0) = 3.979304886336e-01; pointCoord(181, 1) =4.437817512964e-02; pointCoord(181, 2) =2.753444248373e-03;
pointCoord(182, 0) = 6.479304886336e-01; pointCoord(182, 1) =4.437817512964e-02; pointCoord(182, 2) =2.753444248373e-03;
pointCoord(183, 0) = 8.979304886336e-01; pointCoord(183, 1) =4.437817512964e-02; pointCoord(183, 2) =2.753444248373e-03;
pointCoord(184, 0) = 1.479304886336e-01; pointCoord(184, 1) =2.943781751296e-01; pointCoord(184, 2) =2.753444248373e-03;
pointCoord(185, 0) = 3.979304886336e-01; pointCoord(185, 1) =2.943781751296e-01; pointCoord(185, 2) =2.753444248373e-03;
pointCoord(186, 0) = 6.479304886336e-01; pointCoord(186, 1) =2.943781751296e-01; pointCoord(186, 2) =2.753444248373e-03;
pointCoord(187, 0) = 1.479304886336e-01; pointCoord(187, 1) =5.443781751296e-01; pointCoord(187, 2) =2.753444248373e-03;
pointCoord(188, 0) = 3.979304886336e-01; pointCoord(188, 1) =5.443781751296e-01; pointCoord(188, 2) =2.753444248373e-03;
pointCoord(189, 0) = 1.479304886336e-01; pointCoord(189, 1) =7.943781751296e-01; pointCoord(189, 2) =2.753444248373e-03;
pointCoord(190, 0) = 1.832874495324e-01; pointCoord(190, 1) =9.021214230797e-03; pointCoord(190, 2) =1.362987654422e-03;
pointCoord(191, 0) = 4.332874495324e-01; pointCoord(191, 1) =9.021214230797e-03; pointCoord(191, 2) =1.362987654422e-03;
pointCoord(192, 0) = 6.832874495324e-01; pointCoord(192, 1) =9.021214230797e-03; pointCoord(192, 2) =1.362987654422e-03;
pointCoord(193, 0) = 9.332874495324e-01; pointCoord(193, 1) =9.021214230797e-03; pointCoord(193, 2) =1.362987654422e-03;
pointCoord(194, 0) = 1.832874495324e-01; pointCoord(194, 1) =2.590212142308e-01; pointCoord(194, 2) =1.362987654422e-03;
pointCoord(195, 0) = 4.332874495324e-01; pointCoord(195, 1) =2.590212142308e-01; pointCoord(195, 2) =1.362987654422e-03;
pointCoord(196, 0) = 6.832874495324e-01; pointCoord(196, 1) =2.590212142308e-01; pointCoord(196, 2) =1.362987654422e-03;
pointCoord(197, 0) = 1.832874495324e-01; pointCoord(197, 1) =5.090212142308e-01; pointCoord(197, 2) =1.362987654422e-03;
pointCoord(198, 0) = 4.332874495324e-01; pointCoord(198, 1) =5.090212142308e-01; pointCoord(198, 2) =1.362987654422e-03;
pointCoord(199, 0) = 1.832874495324e-01; pointCoord(199, 1) =7.590212142308e-01; pointCoord(199, 2) =1.362987654422e-03;
pointCoord(200, 0) = 1.117738042591e-02; pointCoord(200, 1) =2.270951003164e-01; pointCoord(200, 2) =8.359544098942e-04;
pointCoord(201, 0) = 2.611773804259e-01; pointCoord(201, 1) =2.270951003164e-01; pointCoord(201, 2) =8.359544098942e-04;
pointCoord(202, 0) = 5.111773804259e-01; pointCoord(202, 1) =2.270951003164e-01; pointCoord(202, 2) =8.359544098942e-04;
pointCoord(203, 0) = 7.611773804259e-01; pointCoord(203, 1) =2.270951003164e-01; pointCoord(203, 2) =8.359544098942e-04;
pointCoord(204, 0) = 1.117738042591e-02; pointCoord(204, 1) =4.770951003164e-01; pointCoord(204, 2) =8.359544098942e-04;
pointCoord(205, 0) = 2.611773804259e-01; pointCoord(205, 1) =4.770951003164e-01; pointCoord(205, 2) =8.359544098942e-04;
pointCoord(206, 0) = 5.111773804259e-01; pointCoord(206, 1) =4.770951003164e-01; pointCoord(206, 2) =8.359544098942e-04;
pointCoord(207, 0) = 1.117738042591e-02; pointCoord(207, 1) =7.270951003164e-01; pointCoord(207, 2) =8.359544098942e-04;
pointCoord(208, 0) = 2.611773804259e-01; pointCoord(208, 1) =7.270951003164e-01; pointCoord(208, 2) =8.359544098942e-04;
pointCoord(209, 0) = 1.117738042591e-02; pointCoord(209, 1) =9.770951003164e-01; pointCoord(209, 2) =8.359544098942e-04;
pointCoord(210, 0) = 5.498503120992e-02; pointCoord(210, 1) =1.832874495324e-01; pointCoord(210, 2) =1.688756207261e-03;
pointCoord(211, 0) = 3.049850312099e-01; pointCoord(211, 1) =1.832874495324e-01; pointCoord(211, 2) =1.688756207261e-03;
pointCoord(212, 0) = 5.549850312099e-01; pointCoord(212, 1) =1.832874495324e-01; pointCoord(212, 2) =1.688756207261e-03;
pointCoord(213, 0) = 8.049850312099e-01; pointCoord(213, 1) =1.832874495324e-01; pointCoord(213, 2) =1.688756207261e-03;
pointCoord(214, 0) = 5.498503120992e-02; pointCoord(214, 1) =4.332874495324e-01; pointCoord(214, 2) =1.688756207261e-03;
pointCoord(215, 0) = 3.049850312099e-01; pointCoord(215, 1) =4.332874495324e-01; pointCoord(215, 2) =1.688756207261e-03;
pointCoord(216, 0) = 5.549850312099e-01; pointCoord(216, 1) =4.332874495324e-01; pointCoord(216, 2) =1.688756207261e-03;
pointCoord(217, 0) = 5.498503120992e-02; pointCoord(217, 1) =6.832874495324e-01; pointCoord(217, 2) =1.688756207261e-03;
pointCoord(218, 0) = 3.049850312099e-01; pointCoord(218, 1) =6.832874495324e-01; pointCoord(218, 2) =1.688756207261e-03;
pointCoord(219, 0) = 5.498503120992e-02; pointCoord(219, 1) =9.332874495324e-01; pointCoord(219, 2) =1.688756207261e-03;
pointCoord(220, 0) = 1.191362403712e-01; pointCoord(220, 1) =1.191362403712e-01; pointCoord(220, 2) =2.007223347801e-03;
pointCoord(221, 0) = 3.691362403712e-01; pointCoord(221, 1) =1.191362403712e-01; pointCoord(221, 2) =2.007223347801e-03;
pointCoord(222, 0) = 6.191362403712e-01; pointCoord(222, 1) =1.191362403712e-01; pointCoord(222, 2) =2.007223347801e-03;
pointCoord(223, 0) = 8.691362403712e-01; pointCoord(223, 1) =1.191362403712e-01; pointCoord(223, 2) =2.007223347801e-03;
pointCoord(224, 0) = 1.191362403712e-01; pointCoord(224, 1) =3.691362403712e-01; pointCoord(224, 2) =2.007223347801e-03;
pointCoord(225, 0) = 3.691362403712e-01; pointCoord(225, 1) =3.691362403712e-01; pointCoord(225, 2) =2.007223347801e-03;
pointCoord(226, 0) = 6.191362403712e-01; pointCoord(226, 1) =3.691362403712e-01; pointCoord(226, 2) =2.007223347801e-03;
pointCoord(227, 0) = 1.191362403712e-01; pointCoord(227, 1) =6.191362403712e-01; pointCoord(227, 2) =2.007223347801e-03;
pointCoord(228, 0) = 3.691362403712e-01; pointCoord(228, 1) =6.191362403712e-01; pointCoord(228, 2) =2.007223347801e-03;
pointCoord(229, 0) = 1.191362403712e-01; pointCoord(229, 1) =8.691362403712e-01; pointCoord(229, 2) =2.007223347801e-03;
pointCoord(230, 0) = 1.832874495324e-01; pointCoord(230, 1) =5.498503120992e-02; pointCoord(230, 2) =1.688756207261e-03;
pointCoord(231, 0) = 4.332874495324e-01; pointCoord(231, 1) =5.498503120992e-02; pointCoord(231, 2) =1.688756207261e-03;
pointCoord(232, 0) = 6.832874495324e-01; pointCoord(232, 1) =5.498503120992e-02; pointCoord(232, 2) =1.688756207261e-03;
pointCoord(233, 0) = 9.332874495324e-01; pointCoord(233, 1) =5.498503120992e-02; pointCoord(233, 2) =1.688756207261e-03;
pointCoord(234, 0) = 1.832874495324e-01; pointCoord(234, 1) =3.049850312099e-01; pointCoord(234, 2) =1.688756207261e-03;
pointCoord(235, 0) = 4.332874495324e-01; pointCoord(235, 1) =3.049850312099e-01; pointCoord(235, 2) =1.688756207261e-03;
pointCoord(236, 0) = 6.832874495324e-01; pointCoord(236, 1) =3.049850312099e-01; pointCoord(236, 2) =1.688756207261e-03;
pointCoord(237, 0) = 1.832874495324e-01; pointCoord(237, 1) =5.549850312099e-01; pointCoord(237, 2) =1.688756207261e-03;
pointCoord(238, 0) = 4.332874495324e-01; pointCoord(238, 1) =5.549850312099e-01; pointCoord(238, 2) =1.688756207261e-03;
pointCoord(239, 0) = 1.832874495324e-01; pointCoord(239, 1) =8.049850312099e-01; pointCoord(239, 2) =1.688756207261e-03;
pointCoord(240, 0) = 2.270951003164e-01; pointCoord(240, 1) =1.117738042591e-02; pointCoord(240, 2) =8.359544098942e-04;
pointCoord(241, 0) = 4.770951003164e-01; pointCoord(241, 1) =1.117738042591e-02; pointCoord(241, 2) =8.359544098942e-04;
pointCoord(242, 0) = 7.270951003164e-01; pointCoord(242, 1) =1.117738042591e-02; pointCoord(242, 2) =8.359544098942e-04;
pointCoord(243, 0) = 9.770951003164e-01; pointCoord(243, 1) =1.117738042591e-02; pointCoord(243, 2) =8.359544098942e-04;
pointCoord(244, 0) = 2.270951003164e-01; pointCoord(244, 1) =2.611773804259e-01; pointCoord(244, 2) =8.359544098942e-04;
pointCoord(245, 0) = 4.770951003164e-01; pointCoord(245, 1) =2.611773804259e-01; pointCoord(245, 2) =8.359544098942e-04;
pointCoord(246, 0) = 7.270951003164e-01; pointCoord(246, 1) =2.611773804259e-01; pointCoord(246, 2) =8.359544098942e-04;
pointCoord(247, 0) = 2.270951003164e-01; pointCoord(247, 1) =5.111773804259e-01; pointCoord(247, 2) =8.359544098942e-04;
pointCoord(248, 0) = 4.770951003164e-01; pointCoord(248, 1) =5.111773804259e-01; pointCoord(248, 2) =8.359544098942e-04;
pointCoord(249, 0) = 2.270951003164e-01; pointCoord(249, 1) =7.611773804259e-01; pointCoord(249, 2) =8.359544098942e-04;
pointCoord(250, 0) = 2.494498611682e-01; pointCoord(250, 1) =2.388226195741e-01; pointCoord(250, 2) =4.114479108130e-05;
pointCoord(251, 0) = 4.994498611682e-01; pointCoord(251, 1) =2.388226195741e-01; pointCoord(251, 2) =4.114479108130e-05;
pointCoord(252, 0) = 7.494498611682e-01; pointCoord(252, 1) =2.388226195741e-01; pointCoord(252, 2) =4.114479108130e-05;
pointCoord(253, 0) = 2.494498611682e-01; pointCoord(253, 1) =4.888226195741e-01; pointCoord(253, 2) =4.114479108130e-05;
pointCoord(254, 0) = 4.994498611682e-01; pointCoord(254, 1) =4.888226195741e-01; pointCoord(254, 2) =4.114479108130e-05;
pointCoord(255, 0) = 2.494498611682e-01; pointCoord(255, 1) =7.388226195741e-01; pointCoord(255, 2) =4.114479108130e-05;
pointCoord(256, 0) = 2.472936949731e-01; pointCoord(256, 1) =2.409787857692e-01; pointCoord(256, 2) =8.311879273871e-05;
pointCoord(257, 0) = 4.972936949731e-01; pointCoord(257, 1) =2.409787857692e-01; pointCoord(257, 2) =8.311879273871e-05;
pointCoord(258, 0) = 7.472936949731e-01; pointCoord(258, 1) =2.409787857692e-01; pointCoord(258, 2) =8.311879273871e-05;
pointCoord(259, 0) = 2.472936949731e-01; pointCoord(259, 1) =4.909787857692e-01; pointCoord(259, 2) =8.311879273871e-05;
pointCoord(260, 0) = 4.972936949731e-01; pointCoord(260, 1) =4.909787857692e-01; pointCoord(260, 2) =8.311879273871e-05;
pointCoord(261, 0) = 2.472936949731e-01; pointCoord(261, 1) =7.409787857692e-01; pointCoord(261, 2) =8.311879273871e-05;
pointCoord(262, 0) = 2.441362403712e-01; pointCoord(262, 1) =2.441362403712e-01; pointCoord(262, 2) =9.879340825442e-05;
pointCoord(263, 0) = 4.941362403712e-01; pointCoord(263, 1) =2.441362403712e-01; pointCoord(263, 2) =9.879340825442e-05;
pointCoord(264, 0) = 7.441362403712e-01; pointCoord(264, 1) =2.441362403712e-01; pointCoord(264, 2) =9.879340825442e-05;
pointCoord(265, 0) = 2.441362403712e-01; pointCoord(265, 1) =4.941362403712e-01; pointCoord(265, 2) =9.879340825442e-05;
pointCoord(266, 0) = 4.941362403712e-01; pointCoord(266, 1) =4.941362403712e-01; pointCoord(266, 2) =9.879340825442e-05;
pointCoord(267, 0) = 2.441362403712e-01; pointCoord(267, 1) =7.441362403712e-01; pointCoord(267, 2) =9.879340825442e-05;
pointCoord(268, 0) = 2.409787857692e-01; pointCoord(268, 1) =2.472936949731e-01; pointCoord(268, 2) =8.311879273871e-05;
pointCoord(269, 0) = 4.909787857692e-01; pointCoord(269, 1) =2.472936949731e-01; pointCoord(269, 2) =8.311879273871e-05;
pointCoord(270, 0) = 7.409787857692e-01; pointCoord(270, 1) =2.472936949731e-01; pointCoord(270, 2) =8.311879273871e-05;
pointCoord(271, 0) = 2.409787857692e-01; pointCoord(271, 1) =4.972936949731e-01; pointCoord(271, 2) =8.311879273871e-05;
pointCoord(272, 0) = 4.909787857692e-01; pointCoord(272, 1) =4.972936949731e-01; pointCoord(272, 2) =8.311879273871e-05;
pointCoord(273, 0) = 2.409787857692e-01; pointCoord(273, 1) =7.472936949731e-01; pointCoord(273, 2) =8.311879273871e-05;
pointCoord(274, 0) = 2.388226195741e-01; pointCoord(274, 1) =2.494498611682e-01; pointCoord(274, 2) =4.114479108130e-05;
pointCoord(275, 0) = 4.888226195741e-01; pointCoord(275, 1) =2.494498611682e-01; pointCoord(275, 2) =4.114479108130e-05;
pointCoord(276, 0) = 7.388226195741e-01; pointCoord(276, 1) =2.494498611682e-01; pointCoord(276, 2) =4.114479108130e-05;
pointCoord(277, 0) = 2.388226195741e-01; pointCoord(277, 1) =4.994498611682e-01; pointCoord(277, 2) =4.114479108130e-05;
pointCoord(278, 0) = 4.888226195741e-01; pointCoord(278, 1) =4.994498611682e-01; pointCoord(278, 2) =4.114479108130e-05;
pointCoord(279, 0) = 2.388226195741e-01; pointCoord(279, 1) =7.494498611682e-01; pointCoord(279, 2) =4.114479108130e-05;
pointCoord(280, 0) = 2.472936949731e-01; pointCoord(280, 1) =1.950149687901e-01; pointCoord(280, 2) =4.088873455782e-04;
pointCoord(281, 0) = 4.972936949731e-01; pointCoord(281, 1) =1.950149687901e-01; pointCoord(281, 2) =4.088873455782e-04;
pointCoord(282, 0) = 7.472936949731e-01; pointCoord(282, 1) =1.950149687901e-01; pointCoord(282, 2) =4.088873455782e-04;
pointCoord(283, 0) = 2.472936949731e-01; pointCoord(283, 1) =4.450149687901e-01; pointCoord(283, 2) =4.088873455782e-04;
pointCoord(284, 0) = 4.972936949731e-01; pointCoord(284, 1) =4.450149687901e-01; pointCoord(284, 2) =4.088873455782e-04;
pointCoord(285, 0) = 2.472936949731e-01; pointCoord(285, 1) =6.950149687901e-01; pointCoord(285, 2) =4.088873455782e-04;
pointCoord(286, 0) = 2.366868388929e-01; pointCoord(286, 1) =2.056218248704e-01; pointCoord(286, 2) =8.260151926267e-04;
pointCoord(287, 0) = 4.866868388929e-01; pointCoord(287, 1) =2.056218248704e-01; pointCoord(287, 2) =8.260151926267e-04;
pointCoord(288, 0) = 7.366868388929e-01; pointCoord(288, 1) =2.056218248704e-01; pointCoord(288, 2) =8.260151926267e-04;
pointCoord(289, 0) = 2.366868388929e-01; pointCoord(289, 1) =4.556218248704e-01; pointCoord(289, 2) =8.260151926267e-04;
pointCoord(290, 0) = 4.866868388929e-01; pointCoord(290, 1) =4.556218248704e-01; pointCoord(290, 2) =8.260151926267e-04;
pointCoord(291, 0) = 2.366868388929e-01; pointCoord(291, 1) =7.056218248704e-01; pointCoord(291, 2) =8.260151926267e-04;
pointCoord(292, 0) = 2.211543318816e-01; pointCoord(292, 1) =2.211543318816e-01; pointCoord(292, 2) =9.817858688834e-04;
pointCoord(293, 0) = 4.711543318816e-01; pointCoord(293, 1) =2.211543318816e-01; pointCoord(293, 2) =9.817858688834e-04;
pointCoord(294, 0) = 7.211543318816e-01; pointCoord(294, 1) =2.211543318816e-01; pointCoord(294, 2) =9.817858688834e-04;
pointCoord(295, 0) = 2.211543318816e-01; pointCoord(295, 1) =4.711543318816e-01; pointCoord(295, 2) =9.817858688834e-04;
pointCoord(296, 0) = 4.711543318816e-01; pointCoord(296, 1) =4.711543318816e-01; pointCoord(296, 2) =9.817858688834e-04;
pointCoord(297, 0) = 2.211543318816e-01; pointCoord(297, 1) =7.211543318816e-01; pointCoord(297, 2) =9.817858688834e-04;
pointCoord(298, 0) = 2.056218248704e-01; pointCoord(298, 1) =2.366868388929e-01; pointCoord(298, 2) =8.260151926267e-04;
pointCoord(299, 0) = 4.556218248704e-01; pointCoord(299, 1) =2.366868388929e-01; pointCoord(299, 2) =8.260151926267e-04;
pointCoord(300, 0) = 7.056218248704e-01; pointCoord(300, 1) =2.366868388929e-01; pointCoord(300, 2) =8.260151926267e-04;
pointCoord(301, 0) = 2.056218248704e-01; pointCoord(301, 1) =4.866868388929e-01; pointCoord(301, 2) =8.260151926267e-04;
pointCoord(302, 0) = 4.556218248704e-01; pointCoord(302, 1) =4.866868388929e-01; pointCoord(302, 2) =8.260151926267e-04;
pointCoord(303, 0) = 2.056218248704e-01; pointCoord(303, 1) =7.366868388929e-01; pointCoord(303, 2) =8.260151926267e-04;
pointCoord(304, 0) = 1.950149687901e-01; pointCoord(304, 1) =2.472936949731e-01; pointCoord(304, 2) =4.088873455782e-04;
pointCoord(305, 0) = 4.450149687901e-01; pointCoord(305, 1) =2.472936949731e-01; pointCoord(305, 2) =4.088873455782e-04;
pointCoord(306, 0) = 6.950149687901e-01; pointCoord(306, 1) =2.472936949731e-01; pointCoord(306, 2) =4.088873455782e-04;
pointCoord(307, 0) = 1.950149687901e-01; pointCoord(307, 1) =4.972936949731e-01; pointCoord(307, 2) =4.088873455782e-04;
pointCoord(308, 0) = 4.450149687901e-01; pointCoord(308, 1) =4.972936949731e-01; pointCoord(308, 2) =4.088873455782e-04;
pointCoord(309, 0) = 1.950149687901e-01; pointCoord(309, 1) =7.472936949731e-01; pointCoord(309, 2) =4.088873455782e-04;
pointCoord(310, 0) = 2.441362403712e-01; pointCoord(310, 1) =1.308637596288e-01; pointCoord(310, 2) =1.053008378028e-03;
pointCoord(311, 0) = 4.941362403712e-01; pointCoord(311, 1) =1.308637596288e-01; pointCoord(311, 2) =1.053008378028e-03;
pointCoord(312, 0) = 7.441362403712e-01; pointCoord(312, 1) =1.308637596288e-01; pointCoord(312, 2) =1.053008378028e-03;
pointCoord(313, 0) = 2.441362403712e-01; pointCoord(313, 1) =3.808637596288e-01; pointCoord(313, 2) =1.053008378028e-03;
pointCoord(314, 0) = 4.941362403712e-01; pointCoord(314, 1) =3.808637596288e-01; pointCoord(314, 2) =1.053008378028e-03;
pointCoord(315, 0) = 2.441362403712e-01; pointCoord(315, 1) =6.308637596288e-01; pointCoord(315, 2) =1.053008378028e-03;
pointCoord(316, 0) = 2.211543318816e-01; pointCoord(316, 1) =1.538456681184e-01; pointCoord(316, 2) =2.127238535553e-03;
pointCoord(317, 0) = 4.711543318816e-01; pointCoord(317, 1) =1.538456681184e-01; pointCoord(317, 2) =2.127238535553e-03;
pointCoord(318, 0) = 7.211543318816e-01; pointCoord(318, 1) =1.538456681184e-01; pointCoord(318, 2) =2.127238535553e-03;
pointCoord(319, 0) = 2.211543318816e-01; pointCoord(319, 1) =4.038456681184e-01; pointCoord(319, 2) =2.127238535553e-03;
pointCoord(320, 0) = 4.711543318816e-01; pointCoord(320, 1) =4.038456681184e-01; pointCoord(320, 2) =2.127238535553e-03;
pointCoord(321, 0) = 2.211543318816e-01; pointCoord(321, 1) =6.538456681184e-01; pointCoord(321, 2) =2.127238535553e-03;
pointCoord(322, 0) = 1.875000000000e-01; pointCoord(322, 1) =1.875000000000e-01; pointCoord(322, 2) =2.528395061728e-03;
pointCoord(323, 0) = 4.375000000000e-01; pointCoord(323, 1) =1.875000000000e-01; pointCoord(323, 2) =2.528395061728e-03;
pointCoord(324, 0) = 6.875000000000e-01; pointCoord(324, 1) =1.875000000000e-01; pointCoord(324, 2) =2.528395061728e-03;
pointCoord(325, 0) = 1.875000000000e-01; pointCoord(325, 1) =4.375000000000e-01; pointCoord(325, 2) =2.528395061728e-03;
pointCoord(326, 0) = 4.375000000000e-01; pointCoord(326, 1) =4.375000000000e-01; pointCoord(326, 2) =2.528395061728e-03;
pointCoord(327, 0) = 1.875000000000e-01; pointCoord(327, 1) =6.875000000000e-01; pointCoord(327, 2) =2.528395061728e-03;
pointCoord(328, 0) = 1.538456681184e-01; pointCoord(328, 1) =2.211543318816e-01; pointCoord(328, 2) =2.127238535553e-03;
pointCoord(329, 0) = 4.038456681184e-01; pointCoord(329, 1) =2.211543318816e-01; pointCoord(329, 2) =2.127238535553e-03;
pointCoord(330, 0) = 6.538456681184e-01; pointCoord(330, 1) =2.211543318816e-01; pointCoord(330, 2) =2.127238535553e-03;
pointCoord(331, 0) = 1.538456681184e-01; pointCoord(331, 1) =4.711543318816e-01; pointCoord(331, 2) =2.127238535553e-03;
pointCoord(332, 0) = 4.038456681184e-01; pointCoord(332, 1) =4.711543318816e-01; pointCoord(332, 2) =2.127238535553e-03;
pointCoord(333, 0) = 1.538456681184e-01; pointCoord(333, 1) =7.211543318816e-01; pointCoord(333, 2) =2.127238535553e-03;
pointCoord(334, 0) = 1.308637596288e-01; pointCoord(334, 1) =2.441362403712e-01; pointCoord(334, 2) =1.053008378028e-03;
pointCoord(335, 0) = 3.808637596288e-01; pointCoord(335, 1) =2.441362403712e-01; pointCoord(335, 2) =1.053008378028e-03;
pointCoord(336, 0) = 6.308637596288e-01; pointCoord(336, 1) =2.441362403712e-01; pointCoord(336, 2) =1.053008378028e-03;
pointCoord(337, 0) = 1.308637596288e-01; pointCoord(337, 1) =4.941362403712e-01; pointCoord(337, 2) =1.053008378028e-03;
pointCoord(338, 0) = 3.808637596288e-01; pointCoord(338, 1) =4.941362403712e-01; pointCoord(338, 2) =1.053008378028e-03;
pointCoord(339, 0) = 1.308637596288e-01; pointCoord(339, 1) =7.441362403712e-01; pointCoord(339, 2) =1.053008378028e-03;
pointCoord(340, 0) = 2.409787857692e-01; pointCoord(340, 1) =6.671255046759e-02; pointCoord(340, 2) =1.362987654422e-03;
pointCoord(341, 0) = 4.909787857692e-01; pointCoord(341, 1) =6.671255046759e-02; pointCoord(341, 2) =1.362987654422e-03;
pointCoord(342, 0) = 7.409787857692e-01; pointCoord(342, 1) =6.671255046759e-02; pointCoord(342, 2) =1.362987654422e-03;
pointCoord(343, 0) = 2.409787857692e-01; pointCoord(343, 1) =3.167125504676e-01; pointCoord(343, 2) =1.362987654422e-03;
pointCoord(344, 0) = 4.909787857692e-01; pointCoord(344, 1) =3.167125504676e-01; pointCoord(344, 2) =1.362987654422e-03;
pointCoord(345, 0) = 2.409787857692e-01; pointCoord(345, 1) =5.667125504676e-01; pointCoord(345, 2) =1.362987654422e-03;
pointCoord(346, 0) = 2.056218248704e-01; pointCoord(346, 1) =1.020695113664e-01; pointCoord(346, 2) =2.753444248373e-03;
pointCoord(347, 0) = 4.556218248704e-01; pointCoord(347, 1) =1.020695113664e-01; pointCoord(347, 2) =2.753444248373e-03;
pointCoord(348, 0) = 7.056218248704e-01; pointCoord(348, 1) =1.020695113664e-01; pointCoord(348, 2) =2.753444248373e-03;
pointCoord(349, 0) = 2.056218248704e-01; pointCoord(349, 1) =3.520695113664e-01; pointCoord(349, 2) =2.753444248373e-03;
pointCoord(350, 0) = 4.556218248704e-01; pointCoord(350, 1) =3.520695113664e-01; pointCoord(350, 2) =2.753444248373e-03;
pointCoord(351, 0) = 2.056218248704e-01; pointCoord(351, 1) =6.020695113664e-01; pointCoord(351, 2) =2.753444248373e-03;
pointCoord(352, 0) = 1.538456681184e-01; pointCoord(352, 1) =1.538456681184e-01; pointCoord(352, 2) =3.272691202222e-03;
pointCoord(353, 0) = 4.038456681184e-01; pointCoord(353, 1) =1.538456681184e-01; pointCoord(353, 2) =3.272691202222e-03;
pointCoord(354, 0) = 6.538456681184e-01; pointCoord(354, 1) =1.538456681184e-01; pointCoord(354, 2) =3.272691202222e-03;
pointCoord(355, 0) = 1.538456681184e-01; pointCoord(355, 1) =4.038456681184e-01; pointCoord(355, 2) =3.272691202222e-03;
pointCoord(356, 0) = 4.038456681184e-01; pointCoord(356, 1) =4.038456681184e-01; pointCoord(356, 2) =3.272691202222e-03;
pointCoord(357, 0) = 1.538456681184e-01; pointCoord(357, 1) =6.538456681184e-01; pointCoord(357, 2) =3.272691202222e-03;
pointCoord(358, 0) = 1.020695113664e-01; pointCoord(358, 1) =2.056218248704e-01; pointCoord(358, 2) =2.753444248373e-03;
pointCoord(359, 0) = 3.520695113664e-01; pointCoord(359, 1) =2.056218248704e-01; pointCoord(359, 2) =2.753444248373e-03;
pointCoord(360, 0) = 6.020695113664e-01; pointCoord(360, 1) =2.056218248704e-01; pointCoord(360, 2) =2.753444248373e-03;
pointCoord(361, 0) = 1.020695113664e-01; pointCoord(361, 1) =4.556218248704e-01; pointCoord(361, 2) =2.753444248373e-03;
pointCoord(362, 0) = 3.520695113664e-01; pointCoord(362, 1) =4.556218248704e-01; pointCoord(362, 2) =2.753444248373e-03;
pointCoord(363, 0) = 1.020695113664e-01; pointCoord(363, 1) =7.056218248704e-01; pointCoord(363, 2) =2.753444248373e-03;
pointCoord(364, 0) = 6.671255046759e-02; pointCoord(364, 1) =2.409787857692e-01; pointCoord(364, 2) =1.362987654422e-03;
pointCoord(365, 0) = 3.167125504676e-01; pointCoord(365, 1) =2.409787857692e-01; pointCoord(365, 2) =1.362987654422e-03;
pointCoord(366, 0) = 5.667125504676e-01; pointCoord(366, 1) =2.409787857692e-01; pointCoord(366, 2) =1.362987654422e-03;
pointCoord(367, 0) = 6.671255046759e-02; pointCoord(367, 1) =4.909787857692e-01; pointCoord(367, 2) =1.362987654422e-03;
pointCoord(368, 0) = 3.167125504676e-01; pointCoord(368, 1) =4.909787857692e-01; pointCoord(368, 2) =1.362987654422e-03;
pointCoord(369, 0) = 6.671255046759e-02; pointCoord(369, 1) =7.409787857692e-01; pointCoord(369, 2) =1.362987654422e-03;
pointCoord(370, 0) = 2.388226195741e-01; pointCoord(370, 1) =2.290489968358e-02; pointCoord(370, 2) =8.359544098942e-04;
pointCoord(371, 0) = 4.888226195741e-01; pointCoord(371, 1) =2.290489968358e-02; pointCoord(371, 2) =8.359544098942e-04;
pointCoord(372, 0) = 7.388226195741e-01; pointCoord(372, 1) =2.290489968358e-02; pointCoord(372, 2) =8.359544098942e-04;
pointCoord(373, 0) = 2.388226195741e-01; pointCoord(373, 1) =2.729048996836e-01; pointCoord(373, 2) =8.359544098942e-04;
pointCoord(374, 0) = 4.888226195741e-01; pointCoord(374, 1) =2.729048996836e-01; pointCoord(374, 2) =8.359544098942e-04;
pointCoord(375, 0) = 2.388226195741e-01; pointCoord(375, 1) =5.229048996836e-01; pointCoord(375, 2) =8.359544098942e-04;
pointCoord(376, 0) = 1.950149687901e-01; pointCoord(376, 1) =6.671255046759e-02; pointCoord(376, 2) =1.688756207261e-03;
pointCoord(377, 0) = 4.450149687901e-01; pointCoord(377, 1) =6.671255046759e-02; pointCoord(377, 2) =1.688756207261e-03;
pointCoord(378, 0) = 6.950149687901e-01; pointCoord(378, 1) =6.671255046759e-02; pointCoord(378, 2) =1.688756207261e-03;
pointCoord(379, 0) = 1.950149687901e-01; pointCoord(379, 1) =3.167125504676e-01; pointCoord(379, 2) =1.688756207261e-03;
pointCoord(380, 0) = 4.450149687901e-01; pointCoord(380, 1) =3.167125504676e-01; pointCoord(380, 2) =1.688756207261e-03;
pointCoord(381, 0) = 1.950149687901e-01; pointCoord(381, 1) =5.667125504676e-01; pointCoord(381, 2) =1.688756207261e-03;
pointCoord(382, 0) = 1.308637596288e-01; pointCoord(382, 1) =1.308637596288e-01; pointCoord(382, 2) =2.007223347801e-03;
pointCoord(383, 0) = 3.808637596288e-01; pointCoord(383, 1) =1.308637596288e-01; pointCoord(383, 2) =2.007223347801e-03;
pointCoord(384, 0) = 6.308637596288e-01; pointCoord(384, 1) =1.308637596288e-01; pointCoord(384, 2) =2.007223347801e-03;
pointCoord(385, 0) = 1.308637596288e-01; pointCoord(385, 1) =3.808637596288e-01; pointCoord(385, 2) =2.007223347801e-03;
pointCoord(386, 0) = 3.808637596288e-01; pointCoord(386, 1) =3.808637596288e-01; pointCoord(386, 2) =2.007223347801e-03;
pointCoord(387, 0) = 1.308637596288e-01; pointCoord(387, 1) =6.308637596288e-01; pointCoord(387, 2) =2.007223347801e-03;
pointCoord(388, 0) = 6.671255046759e-02; pointCoord(388, 1) =1.950149687901e-01; pointCoord(388, 2) =1.688756207261e-03;
pointCoord(389, 0) = 3.167125504676e-01; pointCoord(389, 1) =1.950149687901e-01; pointCoord(389, 2) =1.688756207261e-03;
pointCoord(390, 0) = 5.667125504676e-01; pointCoord(390, 1) =1.950149687901e-01; pointCoord(390, 2) =1.688756207261e-03;
pointCoord(391, 0) = 6.671255046759e-02; pointCoord(391, 1) =4.450149687901e-01; pointCoord(391, 2) =1.688756207261e-03;
pointCoord(392, 0) = 3.167125504676e-01; pointCoord(392, 1) =4.450149687901e-01; pointCoord(392, 2) =1.688756207261e-03;
pointCoord(393, 0) = 6.671255046759e-02; pointCoord(393, 1) =6.950149687901e-01; pointCoord(393, 2) =1.688756207261e-03;
pointCoord(394, 0) = 2.290489968358e-02; pointCoord(394, 1) =2.388226195741e-01; pointCoord(394, 2) =8.359544098942e-04;
pointCoord(395, 0) = 2.729048996836e-01; pointCoord(395, 1) =2.388226195741e-01; pointCoord(395, 2) =8.359544098942e-04;
pointCoord(396, 0) = 5.229048996836e-01; pointCoord(396, 1) =2.388226195741e-01; pointCoord(396, 2) =8.359544098942e-04;
pointCoord(397, 0) = 2.290489968358e-02; pointCoord(397, 1) =4.888226195741e-01; pointCoord(397, 2) =8.359544098942e-04;
pointCoord(398, 0) = 2.729048996836e-01; pointCoord(398, 1) =4.888226195741e-01; pointCoord(398, 2) =8.359544098942e-04;
pointCoord(399, 0) = 2.290489968358e-02; pointCoord(399, 1) =7.388226195741e-01; pointCoord(399, 2) =8.359544098942e-04;  }
    else if (nH == 256)
    {
        pointCoord(0, 0) = 1.971119371232826e-04;
        pointCoord(0, 1) = 9.730423938492647e-03;
        pointCoord(0, 2) = 1.271620125244547e-05;
        pointCoord(1, 0) = 5.001971119371232e-01;
        pointCoord(1, 1) = 9.730423938492647e-03;
        pointCoord(1, 2) = 1.271620125244547e-05;
        pointCoord(2, 0) = 1.971119371232826e-04;
        pointCoord(2, 1) = 5.097304239384927e-01;
        pointCoord(2, 2) = 1.271620125244547e-05;
        pointCoord(3, 0) = 1.009300420095791e-03;
        pointCoord(3, 1) = 8.918235455520137e-03;
        pointCoord(3, 2) = 2.793522550523083e-05;
        pointCoord(4, 0) = 5.010093004200958e-01;
        pointCoord(4, 1) = 8.918235455520137e-03;
        pointCoord(4, 2) = 2.793522550523083e-05;
        pointCoord(5, 0) = 1.009300420095791e-03;
        pointCoord(5, 1) = 5.089182354555202e-01;
        pointCoord(5, 2) = 2.793522550523083e-05;
        pointCoord(6, 0) = 2.355147011186338e-03;
        pointCoord(6, 1) = 7.572388864429590e-03;
        pointCoord(6, 2) = 3.940743380670695e-05;
        pointCoord(7, 0) = 5.023551470111863e-01;
        pointCoord(7, 1) = 7.572388864429590e-03;
        pointCoord(7, 2) = 3.940743380670695e-05;
        pointCoord(8, 0) = 2.355147011186338e-03;
        pointCoord(8, 1) = 5.075723888644296e-01;
        pointCoord(8, 2) = 3.940743380670695e-05;
        pointCoord(9, 0) = 4.053240940704791e-03;
        pointCoord(9, 1) = 5.874294934911136e-03;
        pointCoord(9, 2) = 4.555988014296735e-05;
        pointCoord(10, 0) = 5.040532409407048e-01;
        pointCoord(10, 1) = 5.874294934911136e-03;
        pointCoord(10, 2) = 4.555988014296735e-05;
        pointCoord(11, 0) = 4.053240940704791e-03;
        pointCoord(11, 1) = 5.058742949349111e-01;
        pointCoord(11, 2) = 4.555988014296735e-05;
        pointCoord(12, 0) = 5.874294934911136e-03;
        pointCoord(12, 1) = 4.053240940704791e-03;
        pointCoord(12, 2) = 4.555988014296735e-05;
        pointCoord(13, 0) = 5.058742949349111e-01;
        pointCoord(13, 1) = 4.053240940704791e-03;
        pointCoord(13, 2) = 4.555988014296735e-05;
        pointCoord(14, 0) = 5.874294934911136e-03;
        pointCoord(14, 1) = 5.040532409407048e-01;
        pointCoord(14, 2) = 4.555988014296735e-05;
        pointCoord(15, 0) = 7.572388864429590e-03;
        pointCoord(15, 1) = 2.355147011186338e-03;
        pointCoord(15, 2) = 3.940743380670695e-05;
        pointCoord(16, 0) = 5.075723888644296e-01;
        pointCoord(16, 1) = 2.355147011186338e-03;
        pointCoord(16, 2) = 3.940743380670695e-05;
        pointCoord(17, 0) = 7.572388864429590e-03;
        pointCoord(17, 1) = 5.023551470111863e-01;
        pointCoord(17, 2) = 3.940743380670695e-05;
        pointCoord(18, 0) = 8.918235455520137e-03;
        pointCoord(18, 1) = 1.009300420095791e-03;
        pointCoord(18, 2) = 2.793522550523083e-05;
        pointCoord(19, 0) = 5.089182354555202e-01;
        pointCoord(19, 1) = 1.009300420095791e-03;
        pointCoord(19, 2) = 2.793522550523083e-05;
        pointCoord(20, 0) = 8.918235455520137e-03;
        pointCoord(20, 1) = 5.010093004200958e-01;
        pointCoord(20, 2) = 2.793522550523083e-05;
        pointCoord(21, 0) = 9.730423938492647e-03;
        pointCoord(21, 1) = 1.971119371232826e-04;
        pointCoord(21, 2) = 1.271620125244547e-05;
        pointCoord(22, 0) = 5.097304239384927e-01;
        pointCoord(22, 1) = 1.971119371232826e-04;
        pointCoord(22, 2) = 1.271620125244547e-05;
        pointCoord(23, 0) = 9.730423938492647e-03;
        pointCoord(23, 1) = 5.001971119371232e-01;
        pointCoord(23, 2) = 1.271620125244547e-05;
        pointCoord(24, 0) = 1.009300420095791e-03;
        pointCoord(24, 1) = 4.982408022649753e-02;
        pointCoord(24, 2) = 1.430407272608036e-04;
        pointCoord(25, 0) = 5.010093004200958e-01;
        pointCoord(25, 1) = 4.982408022649753e-02;
        pointCoord(25, 2) = 1.430407272608036e-04;
        pointCoord(26, 0) = 1.009300420095791e-03;
        pointCoord(26, 1) = 5.498240802264975e-01;
        pointCoord(26, 2) = 1.430407272608036e-04;
        pointCoord(27, 0) = 5.168065175922896e-03;
        pointCoord(27, 1) = 4.566531547067042e-02;
        pointCoord(27, 2) = 3.142349584703465e-04;
        pointCoord(28, 0) = 5.051680651759229e-01;
        pointCoord(28, 1) = 4.566531547067042e-02;
        pointCoord(28, 2) = 3.142349584703465e-04;
        pointCoord(29, 0) = 5.168065175922896e-03;
        pointCoord(29, 1) = 5.456653154706704e-01;
        pointCoord(29, 2) = 3.142349584703465e-04;
        pointCoord(30, 0) = 1.205939580559753e-02;
        pointCoord(30, 1) = 3.877398484099579e-02;
        pointCoord(30, 2) = 4.432823827878085e-04;
        pointCoord(31, 0) = 5.120593958055976e-01;
        pointCoord(31, 1) = 3.877398484099579e-02;
        pointCoord(31, 2) = 4.432823827878085e-04;
        pointCoord(32, 0) = 1.205939580559753e-02;
        pointCoord(32, 1) = 5.387739848409958e-01;
        pointCoord(32, 2) = 4.432823827878085e-04;
        pointCoord(33, 0) = 2.075438882042010e-02;
        pointCoord(33, 1) = 3.007899182617322e-02;
        pointCoord(33, 2) = 5.124894030999879e-04;
        pointCoord(34, 0) = 5.207543888204201e-01;
        pointCoord(34, 1) = 3.007899182617322e-02;
        pointCoord(34, 2) = 5.124894030999879e-04;
        pointCoord(35, 0) = 2.075438882042010e-02;
        pointCoord(35, 1) = 5.300789918261732e-01;
        pointCoord(35, 2) = 5.124894030999879e-04;
        pointCoord(36, 0) = 3.007899182617322e-02;
        pointCoord(36, 1) = 2.075438882042010e-02;
        pointCoord(36, 2) = 5.124894030999879e-04;
        pointCoord(37, 0) = 5.300789918261732e-01;
        pointCoord(37, 1) = 2.075438882042010e-02;
        pointCoord(37, 2) = 5.124894030999879e-04;
        pointCoord(38, 0) = 3.007899182617322e-02;
        pointCoord(38, 1) = 5.207543888204201e-01;
        pointCoord(38, 2) = 5.124894030999879e-04;
        pointCoord(39, 0) = 3.877398484099579e-02;
        pointCoord(39, 1) = 1.205939580559753e-02;
        pointCoord(39, 2) = 4.432823827878085e-04;
        pointCoord(40, 0) = 5.387739848409958e-01;
        pointCoord(40, 1) = 1.205939580559753e-02;
        pointCoord(40, 2) = 4.432823827878085e-04;
        pointCoord(41, 0) = 3.877398484099579e-02;
        pointCoord(41, 1) = 5.120593958055976e-01;
        pointCoord(41, 2) = 4.432823827878085e-04;
        pointCoord(42, 0) = 4.566531547067042e-02;
        pointCoord(42, 1) = 5.168065175922896e-03;
        pointCoord(42, 2) = 3.142349584703465e-04;
        pointCoord(43, 0) = 5.456653154706704e-01;
        pointCoord(43, 1) = 5.168065175922896e-03;
        pointCoord(43, 2) = 3.142349584703465e-04;
        pointCoord(44, 0) = 4.566531547067042e-02;
        pointCoord(44, 1) = 5.051680651759229e-01;
        pointCoord(44, 2) = 3.142349584703465e-04;
        pointCoord(45, 0) = 4.982408022649753e-02;
        pointCoord(45, 1) = 1.009300420095791e-03;
        pointCoord(45, 2) = 1.430407272608036e-04;
        pointCoord(46, 0) = 5.498240802264975e-01;
        pointCoord(46, 1) = 1.009300420095791e-03;
        pointCoord(46, 2) = 1.430407272608036e-04;
        pointCoord(47, 0) = 4.982408022649753e-02;
        pointCoord(47, 1) = 5.010093004200958e-01;
        pointCoord(47, 2) = 1.430407272608036e-04;
        pointCoord(48, 0) = 2.355147011186338e-03;
        pointCoord(48, 1) = 1.162617505097314e-01;
        pointCoord(48, 2) = 4.708507323447469e-04;
        pointCoord(49, 0) = 5.023551470111863e-01;
        pointCoord(49, 1) = 1.162617505097314e-01;
        pointCoord(49, 2) = 4.708507323447469e-04;
        pointCoord(50, 0) = 2.355147011186338e-03;
        pointCoord(50, 1) = 6.162617505097314e-01;
        pointCoord(50, 2) = 4.708507323447469e-04;
        pointCoord(51, 0) = 1.205939580559753e-02;
        pointCoord(51, 1) = 1.065575017153202e-01;
        pointCoord(51, 2) = 1.034375056373385e-03;
        pointCoord(52, 0) = 5.120593958055976e-01;
        pointCoord(52, 1) = 1.065575017153202e-01;
        pointCoord(52, 2) = 1.034375056373385e-03;
        pointCoord(53, 0) = 1.205939580559753e-02;
        pointCoord(53, 1) = 6.065575017153202e-01;
        pointCoord(53, 2) = 1.034375056373385e-03;
        pointCoord(54, 0) = 2.813993675497581e-02;
        pointCoord(54, 1) = 9.047696076594194e-02;
        pointCoord(54, 2) = 1.459163684134582e-03;
        pointCoord(55, 0) = 5.281399367549758e-01;
        pointCoord(55, 1) = 9.047696076594194e-02;
        pointCoord(55, 2) = 1.459163684134582e-03;
        pointCoord(56, 0) = 2.813993675497581e-02;
        pointCoord(56, 1) = 5.904769607659419e-01;
        pointCoord(56, 2) = 1.459163684134582e-03;
        pointCoord(57, 0) = 4.842922466511254e-02;
        pointCoord(57, 1) = 7.018767285580521e-02;
        pointCoord(57, 2) = 1.686974160363311e-03;
        pointCoord(58, 0) = 5.484292246651126e-01;
        pointCoord(58, 1) = 7.018767285580521e-02;
        pointCoord(58, 2) = 1.686974160363311e-03;
        pointCoord(59, 0) = 4.842922466511254e-02;
        pointCoord(59, 1) = 5.701876728558052e-01;
        pointCoord(59, 2) = 1.686974160363311e-03;
        pointCoord(60, 0) = 7.018767285580521e-02;
        pointCoord(60, 1) = 4.842922466511254e-02;
        pointCoord(60, 2) = 1.686974160363311e-03;
        pointCoord(61, 0) = 5.701876728558052e-01;
        pointCoord(61, 1) = 4.842922466511254e-02;
        pointCoord(61, 2) = 1.686974160363311e-03;
        pointCoord(62, 0) = 7.018767285580521e-02;
        pointCoord(62, 1) = 5.484292246651126e-01;
        pointCoord(62, 2) = 1.686974160363311e-03;
        pointCoord(63, 0) = 9.047696076594194e-02;
        pointCoord(63, 1) = 2.813993675497581e-02;
        pointCoord(63, 2) = 1.459163684134582e-03;
        pointCoord(64, 0) = 5.904769607659419e-01;
        pointCoord(64, 1) = 2.813993675497581e-02;
        pointCoord(64, 2) = 1.459163684134582e-03;
        pointCoord(65, 0) = 9.047696076594194e-02;
        pointCoord(65, 1) = 5.281399367549758e-01;
        pointCoord(65, 2) = 1.459163684134582e-03;
        pointCoord(66, 0) = 1.065575017153202e-01;
        pointCoord(66, 1) = 1.205939580559753e-02;
        pointCoord(66, 2) = 1.034375056373385e-03;
        pointCoord(67, 0) = 6.065575017153202e-01;
        pointCoord(67, 1) = 1.205939580559753e-02;
        pointCoord(67, 2) = 1.034375056373385e-03;
        pointCoord(68, 0) = 1.065575017153202e-01;
        pointCoord(68, 1) = 5.120593958055976e-01;
        pointCoord(68, 2) = 1.034375056373385e-03;
        pointCoord(69, 0) = 1.162617505097314e-01;
        pointCoord(69, 1) = 2.355147011186338e-03;
        pointCoord(69, 2) = 4.708507323447469e-04;
        pointCoord(70, 0) = 6.162617505097314e-01;
        pointCoord(70, 1) = 2.355147011186338e-03;
        pointCoord(70, 2) = 4.708507323447469e-04;
        pointCoord(71, 0) = 1.162617505097314e-01;
        pointCoord(71, 1) = 5.023551470111863e-01;
        pointCoord(71, 2) = 4.708507323447469e-04;
        pointCoord(72, 0) = 4.053240940704791e-03;
        pointCoord(72, 1) = 2.000880984353828e-01;
        pointCoord(72, 2) = 9.368543282773415e-04;
        pointCoord(73, 0) = 5.040532409407048e-01;
        pointCoord(73, 1) = 2.000880984353828e-01;
        pointCoord(73, 2) = 9.368543282773415e-04;
        pointCoord(74, 0) = 4.053240940704791e-03;
        pointCoord(74, 1) = 7.000880984353828e-01;
        pointCoord(74, 2) = 9.368543282773415e-04;
        pointCoord(75, 0) = 2.075438882042010e-02;
        pointCoord(75, 1) = 1.833869505556674e-01;
        pointCoord(75, 2) = 2.058101818807411e-03;
        pointCoord(76, 0) = 5.207543888204201e-01;
        pointCoord(76, 1) = 1.833869505556674e-01;
        pointCoord(76, 2) = 2.058101818807411e-03;
        pointCoord(77, 0) = 2.075438882042010e-02;
        pointCoord(77, 1) = 6.833869505556674e-01;
        pointCoord(77, 2) = 2.058101818807411e-03;
        pointCoord(78, 0) = 4.842922466511254e-02;
        pointCoord(78, 1) = 1.557121147109750e-01;
        pointCoord(78, 2) = 2.903306120687285e-03;
        pointCoord(79, 0) = 5.484292246651126e-01;
        pointCoord(79, 1) = 1.557121147109750e-01;
        pointCoord(79, 2) = 2.903306120687285e-03;
        pointCoord(80, 0) = 4.842922466511254e-02;
        pointCoord(80, 1) = 6.557121147109750e-01;
        pointCoord(80, 2) = 2.903306120687285e-03;
        pointCoord(81, 0) = 8.334737288452591e-02;
        pointCoord(81, 1) = 1.207939664915616e-01;
        pointCoord(81, 2) = 3.356581895833668e-03;
        pointCoord(82, 0) = 5.833473728845259e-01;
        pointCoord(82, 1) = 1.207939664915616e-01;
        pointCoord(82, 2) = 3.356581895833668e-03;
        pointCoord(83, 0) = 8.334737288452591e-02;
        pointCoord(83, 1) = 6.207939664915616e-01;
        pointCoord(83, 2) = 3.356581895833668e-03;
        pointCoord(84, 0) = 1.207939664915616e-01;
        pointCoord(84, 1) = 8.334737288452591e-02;
        pointCoord(84, 2) = 3.356581895833668e-03;
        pointCoord(85, 0) = 6.207939664915616e-01;
        pointCoord(85, 1) = 8.334737288452591e-02;
        pointCoord(85, 2) = 3.356581895833668e-03;
        pointCoord(86, 0) = 1.207939664915616e-01;
        pointCoord(86, 1) = 5.833473728845259e-01;
        pointCoord(86, 2) = 3.356581895833668e-03;
        pointCoord(87, 0) = 1.557121147109750e-01;
        pointCoord(87, 1) = 4.842922466511254e-02;
        pointCoord(87, 2) = 2.903306120687285e-03;
        pointCoord(88, 0) = 6.557121147109750e-01;
        pointCoord(88, 1) = 4.842922466511254e-02;
        pointCoord(88, 2) = 2.903306120687285e-03;
        pointCoord(89, 0) = 1.557121147109750e-01;
        pointCoord(89, 1) = 5.484292246651126e-01;
        pointCoord(89, 2) = 2.903306120687285e-03;
        pointCoord(90, 0) = 1.833869505556674e-01;
        pointCoord(90, 1) = 2.075438882042010e-02;
        pointCoord(90, 2) = 2.058101818807411e-03;
        pointCoord(91, 0) = 6.833869505556674e-01;
        pointCoord(91, 1) = 2.075438882042010e-02;
        pointCoord(91, 2) = 2.058101818807411e-03;
        pointCoord(92, 0) = 1.833869505556674e-01;
        pointCoord(92, 1) = 5.207543888204201e-01;
        pointCoord(92, 2) = 2.058101818807411e-03;
        pointCoord(93, 0) = 2.000880984353828e-01;
        pointCoord(93, 1) = 4.053240940704791e-03;
        pointCoord(93, 2) = 9.368543282773415e-04;
        pointCoord(94, 0) = 7.000880984353828e-01;
        pointCoord(94, 1) = 4.053240940704791e-03;
        pointCoord(94, 2) = 9.368543282773415e-04;
        pointCoord(95, 0) = 2.000880984353828e-01;
        pointCoord(95, 1) = 5.040532409407048e-01;
        pointCoord(95, 2) = 9.368543282773415e-04;
        pointCoord(96, 0) = 5.874294934911136e-03;
        pointCoord(96, 1) = 2.899843656890013e-01;
        pointCoord(96, 2) = 1.357767454700637e-03;
        pointCoord(97, 0) = 5.058742949349111e-01;
        pointCoord(97, 1) = 2.899843656890013e-01;
        pointCoord(97, 2) = 1.357767454700637e-03;
        pointCoord(98, 0) = 5.874294934911136e-03;
        pointCoord(98, 1) = 7.899843656890013e-01;
        pointCoord(98, 2) = 1.357767454700637e-03;
        pointCoord(99, 0) = 3.007899182617322e-02;
        pointCoord(99, 1) = 2.657796687977392e-01;
        pointCoord(99, 2) = 2.982772864139071e-03;
        pointCoord(100, 0) = 5.300789918261732e-01;
        pointCoord(100, 1) = 2.657796687977392e-01;
        pointCoord(100, 2) = 2.982772864139071e-03;
        pointCoord(101, 0) = 3.007899182617322e-02;
        pointCoord(101, 1) = 7.657796687977392e-01;
        pointCoord(101, 2) = 2.982772864139071e-03;
        pointCoord(102, 0) = 7.018767285580521e-02;
        pointCoord(102, 1) = 2.256709877681072e-01;
        pointCoord(102, 2) = 4.207713454183224e-03;
        pointCoord(103, 0) = 5.701876728558052e-01;
        pointCoord(103, 1) = 2.256709877681072e-01;
        pointCoord(103, 2) = 4.207713454183224e-03;
        pointCoord(104, 0) = 7.018767285580521e-02;
        pointCoord(104, 1) = 7.256709877681072e-01;
        pointCoord(104, 2) = 4.207713454183224e-03;
        pointCoord(105, 0) = 1.207939664915616e-01;
        pointCoord(105, 1) = 1.750646941323508e-01;
        pointCoord(105, 2) = 4.864638524518993e-03;
        pointCoord(106, 0) = 6.207939664915616e-01;
        pointCoord(106, 1) = 1.750646941323508e-01;
        pointCoord(106, 2) = 4.864638524518993e-03;
        pointCoord(107, 0) = 1.207939664915616e-01;
        pointCoord(107, 1) = 6.750646941323508e-01;
        pointCoord(107, 2) = 4.864638524518993e-03;
        pointCoord(108, 0) = 1.750646941323508e-01;
        pointCoord(108, 1) = 1.207939664915616e-01;
        pointCoord(108, 2) = 4.864638524518993e-03;
        pointCoord(109, 0) = 6.750646941323508e-01;
        pointCoord(109, 1) = 1.207939664915616e-01;
        pointCoord(109, 2) = 4.864638524518993e-03;
        pointCoord(110, 0) = 1.750646941323508e-01;
        pointCoord(110, 1) = 6.207939664915616e-01;
        pointCoord(110, 2) = 4.864638524518993e-03;
        pointCoord(111, 0) = 2.256709877681072e-01;
        pointCoord(111, 1) = 7.018767285580521e-02;
        pointCoord(111, 2) = 4.207713454183224e-03;
        pointCoord(112, 0) = 7.256709877681072e-01;
        pointCoord(112, 1) = 7.018767285580521e-02;
        pointCoord(112, 2) = 4.207713454183224e-03;
        pointCoord(113, 0) = 2.256709877681072e-01;
        pointCoord(113, 1) = 5.701876728558052e-01;
        pointCoord(113, 2) = 4.207713454183224e-03;
        pointCoord(114, 0) = 2.657796687977392e-01;
        pointCoord(114, 1) = 3.007899182617322e-02;
        pointCoord(114, 2) = 2.982772864139071e-03;
        pointCoord(115, 0) = 7.657796687977392e-01;
        pointCoord(115, 1) = 3.007899182617322e-02;
        pointCoord(115, 2) = 2.982772864139071e-03;
        pointCoord(116, 0) = 2.657796687977392e-01;
        pointCoord(116, 1) = 5.300789918261732e-01;
        pointCoord(116, 2) = 2.982772864139071e-03;
        pointCoord(117, 0) = 2.899843656890013e-01;
        pointCoord(117, 1) = 5.874294934911136e-03;
        pointCoord(117, 2) = 1.357767454700637e-03;
        pointCoord(118, 0) = 7.899843656890013e-01;
        pointCoord(118, 1) = 5.874294934911136e-03;
        pointCoord(118, 2) = 1.357767454700637e-03;
        pointCoord(119, 0) = 2.899843656890013e-01;
        pointCoord(119, 1) = 5.058742949349111e-01;
        pointCoord(119, 2) = 1.357767454700637e-03;
        pointCoord(120, 0) = 7.572388864429590e-03;
        pointCoord(120, 1) = 3.738107136146527e-01;
        pointCoord(120, 2) = 1.513903304329133e-03;
        pointCoord(121, 0) = 5.075723888644296e-01;
        pointCoord(121, 1) = 3.738107136146527e-01;
        pointCoord(121, 2) = 1.513903304329133e-03;
        pointCoord(122, 0) = 7.572388864429590e-03;
        pointCoord(122, 1) = 8.738107136146527e-01;
        pointCoord(122, 2) = 1.513903304329133e-03;
        pointCoord(123, 0) = 3.877398484099579e-02;
        pointCoord(123, 1) = 3.426091176380864e-01;
        pointCoord(123, 2) = 3.325775470203051e-03;
        pointCoord(124, 0) = 5.387739848409958e-01;
        pointCoord(124, 1) = 3.426091176380864e-01;
        pointCoord(124, 2) = 3.325775470203051e-03;
        pointCoord(125, 0) = 3.877398484099579e-02;
        pointCoord(125, 1) = 8.426091176380864e-01;
        pointCoord(125, 2) = 3.325775470203051e-03;
        pointCoord(126, 0) = 9.047696076594194e-02;
        pointCoord(126, 1) = 2.909061417131403e-01;
        pointCoord(126, 2) = 4.691577545112551e-03;
        pointCoord(127, 0) = 5.904769607659419e-01;
        pointCoord(127, 1) = 2.909061417131403e-01;
        pointCoord(127, 2) = 4.691577545112551e-03;
        pointCoord(128, 0) = 9.047696076594194e-02;
        pointCoord(128, 1) = 7.909061417131402e-01;
        pointCoord(128, 2) = 4.691577545112551e-03;
        pointCoord(129, 0) = 1.557121147109750e-01;
        pointCoord(129, 1) = 2.256709877681072e-01;
        pointCoord(129, 2) = 5.424045414507199e-03;
        pointCoord(130, 0) = 6.557121147109750e-01;
        pointCoord(130, 1) = 2.256709877681072e-01;
        pointCoord(130, 2) = 5.424045414507199e-03;
        pointCoord(131, 0) = 1.557121147109750e-01;
        pointCoord(131, 1) = 7.256709877681072e-01;
        pointCoord(131, 2) = 5.424045414507199e-03;
        pointCoord(132, 0) = 2.256709877681072e-01;
        pointCoord(132, 1) = 1.557121147109750e-01;
        pointCoord(132, 2) = 5.424045414507199e-03;
        pointCoord(133, 0) = 7.256709877681072e-01;
        pointCoord(133, 1) = 1.557121147109750e-01;
        pointCoord(133, 2) = 5.424045414507199e-03;
        pointCoord(134, 0) = 2.256709877681072e-01;
        pointCoord(134, 1) = 6.557121147109750e-01;
        pointCoord(134, 2) = 5.424045414507199e-03;
        pointCoord(135, 0) = 2.909061417131403e-01;
        pointCoord(135, 1) = 9.047696076594194e-02;
        pointCoord(135, 2) = 4.691577545112551e-03;
        pointCoord(136, 0) = 7.909061417131402e-01;
        pointCoord(136, 1) = 9.047696076594194e-02;
        pointCoord(136, 2) = 4.691577545112551e-03;
        pointCoord(137, 0) = 2.909061417131403e-01;
        pointCoord(137, 1) = 5.904769607659419e-01;
        pointCoord(137, 2) = 4.691577545112551e-03;
        pointCoord(138, 0) = 3.426091176380864e-01;
        pointCoord(138, 1) = 3.877398484099579e-02;
        pointCoord(138, 2) = 3.325775470203051e-03;
        pointCoord(139, 0) = 8.426091176380864e-01;
        pointCoord(139, 1) = 3.877398484099579e-02;
        pointCoord(139, 2) = 3.325775470203051e-03;
        pointCoord(140, 0) = 3.426091176380864e-01;
        pointCoord(140, 1) = 5.387739848409958e-01;
        pointCoord(140, 2) = 3.325775470203051e-03;
        pointCoord(141, 0) = 3.738107136146527e-01;
        pointCoord(141, 1) = 7.572388864429590e-03;
        pointCoord(141, 2) = 1.513903304329133e-03;
        pointCoord(142, 0) = 8.738107136146527e-01;
        pointCoord(142, 1) = 7.572388864429590e-03;
        pointCoord(142, 2) = 1.513903304329133e-03;
        pointCoord(143, 0) = 3.738107136146527e-01;
        pointCoord(143, 1) = 5.075723888644296e-01;
        pointCoord(143, 2) = 1.513903304329133e-03;
        pointCoord(144, 0) = 8.918235455520137e-03;
        pointCoord(144, 1) = 4.402483838978866e-01;
        pointCoord(144, 2) = 1.263915936267630e-03;
        pointCoord(145, 0) = 5.089182354555202e-01;
        pointCoord(145, 1) = 4.402483838978866e-01;
        pointCoord(145, 2) = 1.263915936267630e-03;
        pointCoord(146, 0) = 8.918235455520137e-03;
        pointCoord(146, 1) = 9.402483838978866e-01;
        pointCoord(146, 2) = 1.263915936267630e-03;
        pointCoord(147, 0) = 4.566531547067042e-02;
        pointCoord(147, 1) = 4.035013038827363e-01;
        pointCoord(147, 2) = 2.776597821814211e-03;
        pointCoord(148, 0) = 5.456653154706704e-01;
        pointCoord(148, 1) = 4.035013038827363e-01;
        pointCoord(148, 2) = 2.776597821814211e-03;
        pointCoord(149, 0) = 4.566531547067042e-02;
        pointCoord(149, 1) = 9.035013038827362e-01;
        pointCoord(149, 2) = 2.776597821814211e-03;
        pointCoord(150, 0) = 1.065575017153202e-01;
        pointCoord(150, 1) = 3.426091176380864e-01;
        pointCoord(150, 2) = 3.916868143788627e-03;
        pointCoord(151, 0) = 6.065575017153202e-01;
        pointCoord(151, 1) = 3.426091176380864e-01;
        pointCoord(151, 2) = 3.916868143788627e-03;
        pointCoord(152, 0) = 1.065575017153202e-01;
        pointCoord(152, 1) = 8.426091176380864e-01;
        pointCoord(152, 2) = 3.916868143788627e-03;
        pointCoord(153, 0) = 1.833869505556674e-01;
        pointCoord(153, 1) = 2.657796687977392e-01;
        pointCoord(153, 2) = 4.528385279846495e-03;
        pointCoord(154, 0) = 6.833869505556674e-01;
        pointCoord(154, 1) = 2.657796687977392e-01;
        pointCoord(154, 2) = 4.528385279846495e-03;
        pointCoord(155, 0) = 1.833869505556674e-01;
        pointCoord(155, 1) = 7.657796687977392e-01;
        pointCoord(155, 2) = 4.528385279846495e-03;
        pointCoord(156, 0) = 2.657796687977392e-01;
        pointCoord(156, 1) = 1.833869505556674e-01;
        pointCoord(156, 2) = 4.528385279846495e-03;
        pointCoord(157, 0) = 7.657796687977392e-01;
        pointCoord(157, 1) = 1.833869505556674e-01;
        pointCoord(157, 2) = 4.528385279846495e-03;
        pointCoord(158, 0) = 2.657796687977392e-01;
        pointCoord(158, 1) = 6.833869505556674e-01;
        pointCoord(158, 2) = 4.528385279846495e-03;
        pointCoord(159, 0) = 3.426091176380864e-01;
        pointCoord(159, 1) = 1.065575017153202e-01;
        pointCoord(159, 2) = 3.916868143788627e-03;
        pointCoord(160, 0) = 8.426091176380864e-01;
        pointCoord(160, 1) = 1.065575017153202e-01;
        pointCoord(160, 2) = 3.916868143788627e-03;
        pointCoord(161, 0) = 3.426091176380864e-01;
        pointCoord(161, 1) = 6.065575017153202e-01;
        pointCoord(161, 2) = 3.916868143788627e-03;
        pointCoord(162, 0) = 4.035013038827363e-01;
        pointCoord(162, 1) = 4.566531547067042e-02;
        pointCoord(162, 2) = 2.776597821814211e-03;
        pointCoord(163, 0) = 9.035013038827362e-01;
        pointCoord(163, 1) = 4.566531547067042e-02;
        pointCoord(163, 2) = 2.776597821814211e-03;
        pointCoord(164, 0) = 4.035013038827363e-01;
        pointCoord(164, 1) = 5.456653154706704e-01;
        pointCoord(164, 2) = 2.776597821814211e-03;
        pointCoord(165, 0) = 4.402483838978866e-01;
        pointCoord(165, 1) = 8.918235455520137e-03;
        pointCoord(165, 2) = 1.263915936267630e-03;
        pointCoord(166, 0) = 9.402483838978866e-01;
        pointCoord(166, 1) = 8.918235455520137e-03;
        pointCoord(166, 2) = 1.263915936267630e-03;
        pointCoord(167, 0) = 4.402483838978866e-01;
        pointCoord(167, 1) = 5.089182354555202e-01;
        pointCoord(167, 2) = 1.263915936267630e-03;
        pointCoord(168, 0) = 9.730423938492647e-03;
        pointCoord(168, 1) = 4.803420401858915e-01;
        pointCoord(168, 2) = 6.277348337158126e-04;
        pointCoord(169, 0) = 5.097304239384927e-01;
        pointCoord(169, 1) = 4.803420401858915e-01;
        pointCoord(169, 2) = 6.277348337158126e-04;
        pointCoord(170, 0) = 9.730423938492647e-03;
        pointCoord(170, 1) = 9.803420401858916e-01;
        pointCoord(170, 2) = 6.277348337158126e-04;
        pointCoord(171, 0) = 4.982408022649753e-02;
        pointCoord(171, 1) = 4.402483838978866e-01;
        pointCoord(171, 2) = 1.379021438023203e-03;
        pointCoord(172, 0) = 5.498240802264975e-01;
        pointCoord(172, 1) = 4.402483838978866e-01;
        pointCoord(172, 2) = 1.379021438023203e-03;
        pointCoord(173, 0) = 4.982408022649753e-02;
        pointCoord(173, 1) = 9.402483838978866e-01;
        pointCoord(173, 2) = 1.379021438023203e-03;
        pointCoord(174, 0) = 1.162617505097314e-01;
        pointCoord(174, 1) = 3.738107136146527e-01;
        pointCoord(174, 2) = 1.945346602867173e-03;
        pointCoord(175, 0) = 6.162617505097314e-01;
        pointCoord(175, 1) = 3.738107136146527e-01;
        pointCoord(175, 2) = 1.945346602867173e-03;
        pointCoord(176, 0) = 1.162617505097314e-01;
        pointCoord(176, 1) = 8.738107136146527e-01;
        pointCoord(176, 2) = 1.945346602867173e-03;
        pointCoord(177, 0) = 2.000880984353828e-01;
        pointCoord(177, 1) = 2.899843656890013e-01;
        pointCoord(177, 2) = 2.249061902835012e-03;
        pointCoord(178, 0) = 7.000880984353828e-01;
        pointCoord(178, 1) = 2.899843656890013e-01;
        pointCoord(178, 2) = 2.249061902835012e-03;
        pointCoord(179, 0) = 2.000880984353828e-01;
        pointCoord(179, 1) = 7.899843656890013e-01;
        pointCoord(179, 2) = 2.249061902835012e-03;
        pointCoord(180, 0) = 2.899843656890013e-01;
        pointCoord(180, 1) = 2.000880984353828e-01;
        pointCoord(180, 2) = 2.249061902835012e-03;
        pointCoord(181, 0) = 7.899843656890013e-01;
        pointCoord(181, 1) = 2.000880984353828e-01;
        pointCoord(181, 2) = 2.249061902835012e-03;
        pointCoord(182, 0) = 2.899843656890013e-01;
        pointCoord(182, 1) = 7.000880984353828e-01;
        pointCoord(182, 2) = 2.249061902835012e-03;
        pointCoord(183, 0) = 3.738107136146527e-01;
        pointCoord(183, 1) = 1.162617505097314e-01;
        pointCoord(183, 2) = 1.945346602867173e-03;
        pointCoord(184, 0) = 8.738107136146527e-01;
        pointCoord(184, 1) = 1.162617505097314e-01;
        pointCoord(184, 2) = 1.945346602867173e-03;
        pointCoord(185, 0) = 3.738107136146527e-01;
        pointCoord(185, 1) = 6.162617505097314e-01;
        pointCoord(185, 2) = 1.945346602867173e-03;
        pointCoord(186, 0) = 4.402483838978866e-01;
        pointCoord(186, 1) = 4.982408022649753e-02;
        pointCoord(186, 2) = 1.379021438023203e-03;
        pointCoord(187, 0) = 9.402483838978866e-01;
        pointCoord(187, 1) = 4.982408022649753e-02;
        pointCoord(187, 2) = 1.379021438023203e-03;
        pointCoord(188, 0) = 4.402483838978866e-01;
        pointCoord(188, 1) = 5.498240802264975e-01;
        pointCoord(188, 2) = 1.379021438023203e-03;
        pointCoord(189, 0) = 4.803420401858915e-01;
        pointCoord(189, 1) = 9.730423938492647e-03;
        pointCoord(189, 2) = 6.277348337158126e-04;
        pointCoord(190, 0) = 9.803420401858916e-01;
        pointCoord(190, 1) = 9.730423938492647e-03;
        pointCoord(190, 2) = 6.277348337158126e-04;
        pointCoord(191, 0) = 4.803420401858915e-01;
        pointCoord(191, 1) = 5.097304239384927e-01;
        pointCoord(191, 2) = 6.277348337158126e-04;
        pointCoord(192, 0) = 4.998028880628767e-01;
        pointCoord(192, 1) = 4.902695760615073e-01;
        pointCoord(192, 2) = 1.271620125244547e-05;
        pointCoord(193, 0) = 4.989906995799042e-01;
        pointCoord(193, 1) = 4.910817645444799e-01;
        pointCoord(193, 2) = 2.793522550523083e-05;
        pointCoord(194, 0) = 4.976448529888137e-01;
        pointCoord(194, 1) = 4.924276111355704e-01;
        pointCoord(194, 2) = 3.940743380670695e-05;
        pointCoord(195, 0) = 4.959467590592952e-01;
        pointCoord(195, 1) = 4.941257050650888e-01;
        pointCoord(195, 2) = 4.555988014296735e-05;
        pointCoord(196, 0) = 4.941257050650888e-01;
        pointCoord(196, 1) = 4.959467590592952e-01;
        pointCoord(196, 2) = 4.555988014296735e-05;
        pointCoord(197, 0) = 4.924276111355704e-01;
        pointCoord(197, 1) = 4.976448529888137e-01;
        pointCoord(197, 2) = 3.940743380670695e-05;
        pointCoord(198, 0) = 4.910817645444799e-01;
        pointCoord(198, 1) = 4.989906995799042e-01;
        pointCoord(198, 2) = 2.793522550523083e-05;
        pointCoord(199, 0) = 4.902695760615073e-01;
        pointCoord(199, 1) = 4.998028880628767e-01;
        pointCoord(199, 2) = 1.271620125244547e-05;
        pointCoord(200, 0) = 4.989906995799042e-01;
        pointCoord(200, 1) = 4.501759197735025e-01;
        pointCoord(200, 2) = 1.430407272608036e-04;
        pointCoord(201, 0) = 4.948319348240771e-01;
        pointCoord(201, 1) = 4.543346845293296e-01;
        pointCoord(201, 2) = 3.142349584703465e-04;
        pointCoord(202, 0) = 4.879406041944025e-01;
        pointCoord(202, 1) = 4.612260151590042e-01;
        pointCoord(202, 2) = 4.432823827878085e-04;
        pointCoord(203, 0) = 4.792456111795799e-01;
        pointCoord(203, 1) = 4.699210081738268e-01;
        pointCoord(203, 2) = 5.124894030999879e-04;
        pointCoord(204, 0) = 4.699210081738268e-01;
        pointCoord(204, 1) = 4.792456111795799e-01;
        pointCoord(204, 2) = 5.124894030999879e-04;
        pointCoord(205, 0) = 4.612260151590042e-01;
        pointCoord(205, 1) = 4.879406041944025e-01;
        pointCoord(205, 2) = 4.432823827878085e-04;
        pointCoord(206, 0) = 4.543346845293296e-01;
        pointCoord(206, 1) = 4.948319348240771e-01;
        pointCoord(206, 2) = 3.142349584703465e-04;
        pointCoord(207, 0) = 4.501759197735025e-01;
        pointCoord(207, 1) = 4.989906995799042e-01;
        pointCoord(207, 2) = 1.430407272608036e-04;
        pointCoord(208, 0) = 4.976448529888137e-01;
        pointCoord(208, 1) = 3.837382494902686e-01;
        pointCoord(208, 2) = 4.708507323447469e-04;
        pointCoord(209, 0) = 4.879406041944025e-01;
        pointCoord(209, 1) = 3.934424982846798e-01;
        pointCoord(209, 2) = 1.034375056373385e-03;
        pointCoord(210, 0) = 4.718600632450242e-01;
        pointCoord(210, 1) = 4.095230392340581e-01;
        pointCoord(210, 2) = 1.459163684134582e-03;
        pointCoord(211, 0) = 4.515707753348874e-01;
        pointCoord(211, 1) = 4.298123271441948e-01;
        pointCoord(211, 2) = 1.686974160363311e-03;
        pointCoord(212, 0) = 4.298123271441948e-01;
        pointCoord(212, 1) = 4.515707753348874e-01;
        pointCoord(212, 2) = 1.686974160363311e-03;
        pointCoord(213, 0) = 4.095230392340581e-01;
        pointCoord(213, 1) = 4.718600632450242e-01;
        pointCoord(213, 2) = 1.459163684134582e-03;
        pointCoord(214, 0) = 3.934424982846798e-01;
        pointCoord(214, 1) = 4.879406041944025e-01;
        pointCoord(214, 2) = 1.034375056373385e-03;
        pointCoord(215, 0) = 3.837382494902686e-01;
        pointCoord(215, 1) = 4.976448529888137e-01;
        pointCoord(215, 2) = 4.708507323447469e-04;
        pointCoord(216, 0) = 4.959467590592952e-01;
        pointCoord(216, 1) = 2.999119015646172e-01;
        pointCoord(216, 2) = 9.368543282773415e-04;
        pointCoord(217, 0) = 4.792456111795799e-01;
        pointCoord(217, 1) = 3.166130494443326e-01;
        pointCoord(217, 2) = 2.058101818807411e-03;
        pointCoord(218, 0) = 4.515707753348874e-01;
        pointCoord(218, 1) = 3.442878852890250e-01;
        pointCoord(218, 2) = 2.903306120687285e-03;
        pointCoord(219, 0) = 4.166526271154741e-01;
        pointCoord(219, 1) = 3.792060335084384e-01;
        pointCoord(219, 2) = 3.356581895833668e-03;
        pointCoord(220, 0) = 3.792060335084384e-01;
        pointCoord(220, 1) = 4.166526271154741e-01;
        pointCoord(220, 2) = 3.356581895833668e-03;
        pointCoord(221, 0) = 3.442878852890250e-01;
        pointCoord(221, 1) = 4.515707753348874e-01;
        pointCoord(221, 2) = 2.903306120687285e-03;
        pointCoord(222, 0) = 3.166130494443326e-01;
        pointCoord(222, 1) = 4.792456111795799e-01;
        pointCoord(222, 2) = 2.058101818807411e-03;
        pointCoord(223, 0) = 2.999119015646172e-01;
        pointCoord(223, 1) = 4.959467590592952e-01;
        pointCoord(223, 2) = 9.368543282773415e-04;
        pointCoord(224, 0) = 4.941257050650888e-01;
        pointCoord(224, 1) = 2.100156343109987e-01;
        pointCoord(224, 2) = 1.357767454700637e-03;
        pointCoord(225, 0) = 4.699210081738268e-01;
        pointCoord(225, 1) = 2.342203312022608e-01;
        pointCoord(225, 2) = 2.982772864139071e-03;
        pointCoord(226, 0) = 4.298123271441948e-01;
        pointCoord(226, 1) = 2.743290122318928e-01;
        pointCoord(226, 2) = 4.207713454183224e-03;
        pointCoord(227, 0) = 3.792060335084384e-01;
        pointCoord(227, 1) = 3.249353058676492e-01;
        pointCoord(227, 2) = 4.864638524518993e-03;
        pointCoord(228, 0) = 3.249353058676492e-01;
        pointCoord(228, 1) = 3.792060335084384e-01;
        pointCoord(228, 2) = 4.864638524518993e-03;
        pointCoord(229, 0) = 2.743290122318928e-01;
        pointCoord(229, 1) = 4.298123271441948e-01;
        pointCoord(229, 2) = 4.207713454183224e-03;
        pointCoord(230, 0) = 2.342203312022608e-01;
        pointCoord(230, 1) = 4.699210081738268e-01;
        pointCoord(230, 2) = 2.982772864139071e-03;
        pointCoord(231, 0) = 2.100156343109987e-01;
        pointCoord(231, 1) = 4.941257050650888e-01;
        pointCoord(231, 2) = 1.357767454700637e-03;
        pointCoord(232, 0) = 4.924276111355704e-01;
        pointCoord(232, 1) = 1.261892863853473e-01;
        pointCoord(232, 2) = 1.513903304329133e-03;
        pointCoord(233, 0) = 4.612260151590042e-01;
        pointCoord(233, 1) = 1.573908823619136e-01;
        pointCoord(233, 2) = 3.325775470203051e-03;
        pointCoord(234, 0) = 4.095230392340581e-01;
        pointCoord(234, 1) = 2.090938582868597e-01;
        pointCoord(234, 2) = 4.691577545112551e-03;
        pointCoord(235, 0) = 3.442878852890250e-01;
        pointCoord(235, 1) = 2.743290122318928e-01;
        pointCoord(235, 2) = 5.424045414507199e-03;
        pointCoord(236, 0) = 2.743290122318928e-01;
        pointCoord(236, 1) = 3.442878852890250e-01;
        pointCoord(236, 2) = 5.424045414507199e-03;
        pointCoord(237, 0) = 2.090938582868597e-01;
        pointCoord(237, 1) = 4.095230392340581e-01;
        pointCoord(237, 2) = 4.691577545112551e-03;
        pointCoord(238, 0) = 1.573908823619136e-01;
        pointCoord(238, 1) = 4.612260151590042e-01;
        pointCoord(238, 2) = 3.325775470203051e-03;
        pointCoord(239, 0) = 1.261892863853473e-01;
        pointCoord(239, 1) = 4.924276111355704e-01;
        pointCoord(239, 2) = 1.513903304329133e-03;
        pointCoord(240, 0) = 4.910817645444799e-01;
        pointCoord(240, 1) = 5.975161610211344e-02;
        pointCoord(240, 2) = 1.263915936267630e-03;
        pointCoord(241, 0) = 4.543346845293296e-01;
        pointCoord(241, 1) = 9.649869611726375e-02;
        pointCoord(241, 2) = 2.776597821814211e-03;
        pointCoord(242, 0) = 3.934424982846798e-01;
        pointCoord(242, 1) = 1.573908823619136e-01;
        pointCoord(242, 2) = 3.916868143788627e-03;
        pointCoord(243, 0) = 3.166130494443326e-01;
        pointCoord(243, 1) = 2.342203312022608e-01;
        pointCoord(243, 2) = 4.528385279846495e-03;
        pointCoord(244, 0) = 2.342203312022608e-01;
        pointCoord(244, 1) = 3.166130494443326e-01;
        pointCoord(244, 2) = 4.528385279846495e-03;
        pointCoord(245, 0) = 1.573908823619136e-01;
        pointCoord(245, 1) = 3.934424982846798e-01;
        pointCoord(245, 2) = 3.916868143788627e-03;
        pointCoord(246, 0) = 9.649869611726375e-02;
        pointCoord(246, 1) = 4.543346845293296e-01;
        pointCoord(246, 2) = 2.776597821814211e-03;
        pointCoord(247, 0) = 5.975161610211344e-02;
        pointCoord(247, 1) = 4.910817645444799e-01;
        pointCoord(247, 2) = 1.263915936267630e-03;
        pointCoord(248, 0) = 4.902695760615073e-01;
        pointCoord(248, 1) = 1.965795981410851e-02;
        pointCoord(248, 2) = 6.277348337158126e-04;
        pointCoord(249, 0) = 4.501759197735025e-01;
        pointCoord(249, 1) = 5.975161610211344e-02;
        pointCoord(249, 2) = 1.379021438023203e-03;
        pointCoord(250, 0) = 3.837382494902686e-01;
        pointCoord(250, 1) = 1.261892863853473e-01;
        pointCoord(250, 2) = 1.945346602867173e-03;
        pointCoord(251, 0) = 2.999119015646172e-01;
        pointCoord(251, 1) = 2.100156343109987e-01;
        pointCoord(251, 2) = 2.249061902835012e-03;
        pointCoord(252, 0) = 2.100156343109987e-01;
        pointCoord(252, 1) = 2.999119015646172e-01;
        pointCoord(252, 2) = 2.249061902835012e-03;
        pointCoord(253, 0) = 1.261892863853473e-01;
        pointCoord(253, 1) = 3.837382494902686e-01;
        pointCoord(253, 2) = 1.945346602867173e-03;
        pointCoord(254, 0) = 5.975161610211344e-02;
        pointCoord(254, 1) = 4.501759197735025e-01;
        pointCoord(254, 2) = 1.379021438023203e-03;
        pointCoord(255, 0) = 1.965795981410851e-02;
        pointCoord(255, 1) = 4.902695760615073e-01;
        pointCoord(255, 2) = 6.277348337158126e-04;
    }
    else if (nH == 1024)
    {
        pointCoord(0, 0) = 9.855596856164131e-05; pointCoord(0, 1) =4.865211969246323e-03; pointCoord(0, 2) =3.179050313111366e-06;
pointCoord(1, 0) = 2.500985559685616e-01; pointCoord(1, 1) =4.865211969246323e-03; pointCoord(1, 2) =3.179050313111366e-06;
pointCoord(2, 0) = 5.000985559685617e-01; pointCoord(2, 1) =4.865211969246323e-03; pointCoord(2, 2) =3.179050313111366e-06;
pointCoord(3, 0) = 7.500985559685617e-01; pointCoord(3, 1) =4.865211969246323e-03; pointCoord(3, 2) =3.179050313111366e-06;
pointCoord(4, 0) = 9.855596856164131e-05; pointCoord(4, 1) =2.548652119692463e-01; pointCoord(4, 2) =3.179050313111366e-06;
pointCoord(5, 0) = 2.500985559685616e-01; pointCoord(5, 1) =2.548652119692463e-01; pointCoord(5, 2) =3.179050313111366e-06;
pointCoord(6, 0) = 5.000985559685617e-01; pointCoord(6, 1) =2.548652119692463e-01; pointCoord(6, 2) =3.179050313111366e-06;
pointCoord(7, 0) = 9.855596856164131e-05; pointCoord(7, 1) =5.048652119692463e-01; pointCoord(7, 2) =3.179050313111366e-06;
pointCoord(8, 0) = 2.500985559685616e-01; pointCoord(8, 1) =5.048652119692463e-01; pointCoord(8, 2) =3.179050313111366e-06;
pointCoord(9, 0) = 9.855596856164131e-05; pointCoord(9, 1) =7.548652119692463e-01; pointCoord(9, 2) =3.179050313111366e-06;
pointCoord(10, 0) = 5.046502100478956e-04; pointCoord(10, 1) =4.459117727760069e-03; pointCoord(10, 2) =6.983806376307709e-06;
pointCoord(11, 0) = 2.505046502100479e-01; pointCoord(11, 1) =4.459117727760069e-03; pointCoord(11, 2) =6.983806376307709e-06;
pointCoord(12, 0) = 5.005046502100479e-01; pointCoord(12, 1) =4.459117727760069e-03; pointCoord(12, 2) =6.983806376307709e-06;
pointCoord(13, 0) = 7.505046502100479e-01; pointCoord(13, 1) =4.459117727760069e-03; pointCoord(13, 2) =6.983806376307709e-06;
pointCoord(14, 0) = 5.046502100478956e-04; pointCoord(14, 1) =2.544591177277601e-01; pointCoord(14, 2) =6.983806376307709e-06;
pointCoord(15, 0) = 2.505046502100479e-01; pointCoord(15, 1) =2.544591177277601e-01; pointCoord(15, 2) =6.983806376307709e-06;
pointCoord(16, 0) = 5.005046502100479e-01; pointCoord(16, 1) =2.544591177277601e-01; pointCoord(16, 2) =6.983806376307709e-06;
pointCoord(17, 0) = 5.046502100478956e-04; pointCoord(17, 1) =5.044591177277601e-01; pointCoord(17, 2) =6.983806376307709e-06;
pointCoord(18, 0) = 2.505046502100479e-01; pointCoord(18, 1) =5.044591177277601e-01; pointCoord(18, 2) =6.983806376307709e-06;
pointCoord(19, 0) = 5.046502100478956e-04; pointCoord(19, 1) =7.544591177277601e-01; pointCoord(19, 2) =6.983806376307709e-06;
pointCoord(20, 0) = 1.177573505593169e-03; pointCoord(20, 1) =3.786194432214795e-03; pointCoord(20, 2) =9.851858451676738e-06;
pointCoord(21, 0) = 2.511775735055932e-01; pointCoord(21, 1) =3.786194432214795e-03; pointCoord(21, 2) =9.851858451676738e-06;
pointCoord(22, 0) = 5.011775735055932e-01; pointCoord(22, 1) =3.786194432214795e-03; pointCoord(22, 2) =9.851858451676738e-06;
pointCoord(23, 0) = 7.511775735055932e-01; pointCoord(23, 1) =3.786194432214795e-03; pointCoord(23, 2) =9.851858451676738e-06;
pointCoord(24, 0) = 1.177573505593169e-03; pointCoord(24, 1) =2.537861944322148e-01; pointCoord(24, 2) =9.851858451676738e-06;
pointCoord(25, 0) = 2.511775735055932e-01; pointCoord(25, 1) =2.537861944322148e-01; pointCoord(25, 2) =9.851858451676738e-06;
pointCoord(26, 0) = 5.011775735055932e-01; pointCoord(26, 1) =2.537861944322148e-01; pointCoord(26, 2) =9.851858451676738e-06;
pointCoord(27, 0) = 1.177573505593169e-03; pointCoord(27, 1) =5.037861944322148e-01; pointCoord(27, 2) =9.851858451676738e-06;
pointCoord(28, 0) = 2.511775735055932e-01; pointCoord(28, 1) =5.037861944322148e-01; pointCoord(28, 2) =9.851858451676738e-06;
pointCoord(29, 0) = 1.177573505593169e-03; pointCoord(29, 1) =7.537861944322148e-01; pointCoord(29, 2) =9.851858451676738e-06;
pointCoord(30, 0) = 2.026620470352396e-03; pointCoord(30, 1) =2.937147467455568e-03; pointCoord(30, 2) =1.138997003574184e-05;
pointCoord(31, 0) = 2.520266204703524e-01; pointCoord(31, 1) =2.937147467455568e-03; pointCoord(31, 2) =1.138997003574184e-05;
pointCoord(32, 0) = 5.020266204703524e-01; pointCoord(32, 1) =2.937147467455568e-03; pointCoord(32, 2) =1.138997003574184e-05;
pointCoord(33, 0) = 7.520266204703524e-01; pointCoord(33, 1) =2.937147467455568e-03; pointCoord(33, 2) =1.138997003574184e-05;
pointCoord(34, 0) = 2.026620470352396e-03; pointCoord(34, 1) =2.529371474674555e-01; pointCoord(34, 2) =1.138997003574184e-05;
pointCoord(35, 0) = 2.520266204703524e-01; pointCoord(35, 1) =2.529371474674555e-01; pointCoord(35, 2) =1.138997003574184e-05;
pointCoord(36, 0) = 5.020266204703524e-01; pointCoord(36, 1) =2.529371474674555e-01; pointCoord(36, 2) =1.138997003574184e-05;
pointCoord(37, 0) = 2.026620470352396e-03; pointCoord(37, 1) =5.029371474674555e-01; pointCoord(37, 2) =1.138997003574184e-05;
pointCoord(38, 0) = 2.520266204703524e-01; pointCoord(38, 1) =5.029371474674555e-01; pointCoord(38, 2) =1.138997003574184e-05;
pointCoord(39, 0) = 2.026620470352396e-03; pointCoord(39, 1) =7.529371474674555e-01; pointCoord(39, 2) =1.138997003574184e-05;
pointCoord(40, 0) = 2.937147467455568e-03; pointCoord(40, 1) =2.026620470352396e-03; pointCoord(40, 2) =1.138997003574184e-05;
pointCoord(41, 0) = 2.529371474674555e-01; pointCoord(41, 1) =2.026620470352396e-03; pointCoord(41, 2) =1.138997003574184e-05;
pointCoord(42, 0) = 5.029371474674555e-01; pointCoord(42, 1) =2.026620470352396e-03; pointCoord(42, 2) =1.138997003574184e-05;
pointCoord(43, 0) = 7.529371474674555e-01; pointCoord(43, 1) =2.026620470352396e-03; pointCoord(43, 2) =1.138997003574184e-05;
pointCoord(44, 0) = 2.937147467455568e-03; pointCoord(44, 1) =2.520266204703524e-01; pointCoord(44, 2) =1.138997003574184e-05;
pointCoord(45, 0) = 2.529371474674555e-01; pointCoord(45, 1) =2.520266204703524e-01; pointCoord(45, 2) =1.138997003574184e-05;
pointCoord(46, 0) = 5.029371474674555e-01; pointCoord(46, 1) =2.520266204703524e-01; pointCoord(46, 2) =1.138997003574184e-05;
pointCoord(47, 0) = 2.937147467455568e-03; pointCoord(47, 1) =5.020266204703524e-01; pointCoord(47, 2) =1.138997003574184e-05;
pointCoord(48, 0) = 2.529371474674555e-01; pointCoord(48, 1) =5.020266204703524e-01; pointCoord(48, 2) =1.138997003574184e-05;
pointCoord(49, 0) = 2.937147467455568e-03; pointCoord(49, 1) =7.520266204703524e-01; pointCoord(49, 2) =1.138997003574184e-05;
pointCoord(50, 0) = 3.786194432214795e-03; pointCoord(50, 1) =1.177573505593169e-03; pointCoord(50, 2) =9.851858451676738e-06;
pointCoord(51, 0) = 2.537861944322148e-01; pointCoord(51, 1) =1.177573505593169e-03; pointCoord(51, 2) =9.851858451676738e-06;
pointCoord(52, 0) = 5.037861944322148e-01; pointCoord(52, 1) =1.177573505593169e-03; pointCoord(52, 2) =9.851858451676738e-06;
pointCoord(53, 0) = 7.537861944322148e-01; pointCoord(53, 1) =1.177573505593169e-03; pointCoord(53, 2) =9.851858451676738e-06;
pointCoord(54, 0) = 3.786194432214795e-03; pointCoord(54, 1) =2.511775735055932e-01; pointCoord(54, 2) =9.851858451676738e-06;
pointCoord(55, 0) = 2.537861944322148e-01; pointCoord(55, 1) =2.511775735055932e-01; pointCoord(55, 2) =9.851858451676738e-06;
pointCoord(56, 0) = 5.037861944322148e-01; pointCoord(56, 1) =2.511775735055932e-01; pointCoord(56, 2) =9.851858451676738e-06;
pointCoord(57, 0) = 3.786194432214795e-03; pointCoord(57, 1) =5.011775735055932e-01; pointCoord(57, 2) =9.851858451676738e-06;
pointCoord(58, 0) = 2.537861944322148e-01; pointCoord(58, 1) =5.011775735055932e-01; pointCoord(58, 2) =9.851858451676738e-06;
pointCoord(59, 0) = 3.786194432214795e-03; pointCoord(59, 1) =7.511775735055932e-01; pointCoord(59, 2) =9.851858451676738e-06;
pointCoord(60, 0) = 4.459117727760069e-03; pointCoord(60, 1) =5.046502100478956e-04; pointCoord(60, 2) =6.983806376307709e-06;
pointCoord(61, 0) = 2.544591177277601e-01; pointCoord(61, 1) =5.046502100478956e-04; pointCoord(61, 2) =6.983806376307709e-06;
pointCoord(62, 0) = 5.044591177277601e-01; pointCoord(62, 1) =5.046502100478956e-04; pointCoord(62, 2) =6.983806376307709e-06;
pointCoord(63, 0) = 7.544591177277601e-01; pointCoord(63, 1) =5.046502100478956e-04; pointCoord(63, 2) =6.983806376307709e-06;
pointCoord(64, 0) = 4.459117727760069e-03; pointCoord(64, 1) =2.505046502100479e-01; pointCoord(64, 2) =6.983806376307709e-06;
pointCoord(65, 0) = 2.544591177277601e-01; pointCoord(65, 1) =2.505046502100479e-01; pointCoord(65, 2) =6.983806376307709e-06;
pointCoord(66, 0) = 5.044591177277601e-01; pointCoord(66, 1) =2.505046502100479e-01; pointCoord(66, 2) =6.983806376307709e-06;
pointCoord(67, 0) = 4.459117727760069e-03; pointCoord(67, 1) =5.005046502100479e-01; pointCoord(67, 2) =6.983806376307709e-06;
pointCoord(68, 0) = 2.544591177277601e-01; pointCoord(68, 1) =5.005046502100479e-01; pointCoord(68, 2) =6.983806376307709e-06;
pointCoord(69, 0) = 4.459117727760069e-03; pointCoord(69, 1) =7.505046502100479e-01; pointCoord(69, 2) =6.983806376307709e-06;
pointCoord(70, 0) = 4.865211969246323e-03; pointCoord(70, 1) =9.855596856164131e-05; pointCoord(70, 2) =3.179050313111366e-06;
pointCoord(71, 0) = 2.548652119692463e-01; pointCoord(71, 1) =9.855596856164131e-05; pointCoord(71, 2) =3.179050313111366e-06;
pointCoord(72, 0) = 5.048652119692463e-01; pointCoord(72, 1) =9.855596856164131e-05; pointCoord(72, 2) =3.179050313111366e-06;
pointCoord(73, 0) = 7.548652119692463e-01; pointCoord(73, 1) =9.855596856164131e-05; pointCoord(73, 2) =3.179050313111366e-06;
pointCoord(74, 0) = 4.865211969246323e-03; pointCoord(74, 1) =2.500985559685616e-01; pointCoord(74, 2) =3.179050313111366e-06;
pointCoord(75, 0) = 2.548652119692463e-01; pointCoord(75, 1) =2.500985559685616e-01; pointCoord(75, 2) =3.179050313111366e-06;
pointCoord(76, 0) = 5.048652119692463e-01; pointCoord(76, 1) =2.500985559685616e-01; pointCoord(76, 2) =3.179050313111366e-06;
pointCoord(77, 0) = 4.865211969246323e-03; pointCoord(77, 1) =5.000985559685617e-01; pointCoord(77, 2) =3.179050313111366e-06;
pointCoord(78, 0) = 2.548652119692463e-01; pointCoord(78, 1) =5.000985559685617e-01; pointCoord(78, 2) =3.179050313111366e-06;
pointCoord(79, 0) = 4.865211969246323e-03; pointCoord(79, 1) =7.500985559685617e-01; pointCoord(79, 2) =3.179050313111366e-06;
pointCoord(80, 0) = 5.046502100478956e-04; pointCoord(80, 1) =2.491204011324876e-02; pointCoord(80, 2) =3.576018181520091e-05;
pointCoord(81, 0) = 2.505046502100479e-01; pointCoord(81, 1) =2.491204011324876e-02; pointCoord(81, 2) =3.576018181520091e-05;
pointCoord(82, 0) = 5.005046502100479e-01; pointCoord(82, 1) =2.491204011324876e-02; pointCoord(82, 2) =3.576018181520091e-05;
pointCoord(83, 0) = 7.505046502100479e-01; pointCoord(83, 1) =2.491204011324876e-02; pointCoord(83, 2) =3.576018181520091e-05;
pointCoord(84, 0) = 5.046502100478956e-04; pointCoord(84, 1) =2.749120401132488e-01; pointCoord(84, 2) =3.576018181520091e-05;
pointCoord(85, 0) = 2.505046502100479e-01; pointCoord(85, 1) =2.749120401132488e-01; pointCoord(85, 2) =3.576018181520091e-05;
pointCoord(86, 0) = 5.005046502100479e-01; pointCoord(86, 1) =2.749120401132488e-01; pointCoord(86, 2) =3.576018181520091e-05;
pointCoord(87, 0) = 5.046502100478956e-04; pointCoord(87, 1) =5.249120401132488e-01; pointCoord(87, 2) =3.576018181520091e-05;
pointCoord(88, 0) = 2.505046502100479e-01; pointCoord(88, 1) =5.249120401132488e-01; pointCoord(88, 2) =3.576018181520091e-05;
pointCoord(89, 0) = 5.046502100478956e-04; pointCoord(89, 1) =7.749120401132488e-01; pointCoord(89, 2) =3.576018181520091e-05;
pointCoord(90, 0) = 2.584032587961448e-03; pointCoord(90, 1) =2.283265773533521e-02; pointCoord(90, 2) =7.855873961758663e-05;
pointCoord(91, 0) = 2.525840325879615e-01; pointCoord(91, 1) =2.283265773533521e-02; pointCoord(91, 2) =7.855873961758663e-05;
pointCoord(92, 0) = 5.025840325879615e-01; pointCoord(92, 1) =2.283265773533521e-02; pointCoord(92, 2) =7.855873961758663e-05;
pointCoord(93, 0) = 7.525840325879615e-01; pointCoord(93, 1) =2.283265773533521e-02; pointCoord(93, 2) =7.855873961758663e-05;
pointCoord(94, 0) = 2.584032587961448e-03; pointCoord(94, 1) =2.728326577353352e-01; pointCoord(94, 2) =7.855873961758663e-05;
pointCoord(95, 0) = 2.525840325879615e-01; pointCoord(95, 1) =2.728326577353352e-01; pointCoord(95, 2) =7.855873961758663e-05;
pointCoord(96, 0) = 5.025840325879615e-01; pointCoord(96, 1) =2.728326577353352e-01; pointCoord(96, 2) =7.855873961758663e-05;
pointCoord(97, 0) = 2.584032587961448e-03; pointCoord(97, 1) =5.228326577353352e-01; pointCoord(97, 2) =7.855873961758663e-05;
pointCoord(98, 0) = 2.525840325879615e-01; pointCoord(98, 1) =5.228326577353352e-01; pointCoord(98, 2) =7.855873961758663e-05;
pointCoord(99, 0) = 2.584032587961448e-03; pointCoord(99, 1) =7.728326577353352e-01; pointCoord(99, 2) =7.855873961758663e-05;
pointCoord(100, 0) = 6.029697902798763e-03; pointCoord(100, 1) =1.938699242049790e-02; pointCoord(100, 2) =1.108205956969521e-04;
pointCoord(101, 0) = 2.560296979027988e-01; pointCoord(101, 1) =1.938699242049790e-02; pointCoord(101, 2) =1.108205956969521e-04;
pointCoord(102, 0) = 5.060296979027987e-01; pointCoord(102, 1) =1.938699242049790e-02; pointCoord(102, 2) =1.108205956969521e-04;
pointCoord(103, 0) = 7.560296979027987e-01; pointCoord(103, 1) =1.938699242049790e-02; pointCoord(103, 2) =1.108205956969521e-04;
pointCoord(104, 0) = 6.029697902798763e-03; pointCoord(104, 1) =2.693869924204979e-01; pointCoord(104, 2) =1.108205956969521e-04;
pointCoord(105, 0) = 2.560296979027988e-01; pointCoord(105, 1) =2.693869924204979e-01; pointCoord(105, 2) =1.108205956969521e-04;
pointCoord(106, 0) = 5.060296979027987e-01; pointCoord(106, 1) =2.693869924204979e-01; pointCoord(106, 2) =1.108205956969521e-04;
pointCoord(107, 0) = 6.029697902798763e-03; pointCoord(107, 1) =5.193869924204979e-01; pointCoord(107, 2) =1.108205956969521e-04;
pointCoord(108, 0) = 2.560296979027988e-01; pointCoord(108, 1) =5.193869924204979e-01; pointCoord(108, 2) =1.108205956969521e-04;
pointCoord(109, 0) = 6.029697902798763e-03; pointCoord(109, 1) =7.693869924204979e-01; pointCoord(109, 2) =1.108205956969521e-04;
pointCoord(110, 0) = 1.037719441021005e-02; pointCoord(110, 1) =1.503949591308661e-02; pointCoord(110, 2) =1.281223507749970e-04;
pointCoord(111, 0) = 2.603771944102101e-01; pointCoord(111, 1) =1.503949591308661e-02; pointCoord(111, 2) =1.281223507749970e-04;
pointCoord(112, 0) = 5.103771944102100e-01; pointCoord(112, 1) =1.503949591308661e-02; pointCoord(112, 2) =1.281223507749970e-04;
pointCoord(113, 0) = 7.603771944102100e-01; pointCoord(113, 1) =1.503949591308661e-02; pointCoord(113, 2) =1.281223507749970e-04;
pointCoord(114, 0) = 1.037719441021005e-02; pointCoord(114, 1) =2.650394959130866e-01; pointCoord(114, 2) =1.281223507749970e-04;
pointCoord(115, 0) = 2.603771944102101e-01; pointCoord(115, 1) =2.650394959130866e-01; pointCoord(115, 2) =1.281223507749970e-04;
pointCoord(116, 0) = 5.103771944102100e-01; pointCoord(116, 1) =2.650394959130866e-01; pointCoord(116, 2) =1.281223507749970e-04;
pointCoord(117, 0) = 1.037719441021005e-02; pointCoord(117, 1) =5.150394959130866e-01; pointCoord(117, 2) =1.281223507749970e-04;
pointCoord(118, 0) = 2.603771944102101e-01; pointCoord(118, 1) =5.150394959130866e-01; pointCoord(118, 2) =1.281223507749970e-04;
pointCoord(119, 0) = 1.037719441021005e-02; pointCoord(119, 1) =7.650394959130866e-01; pointCoord(119, 2) =1.281223507749970e-04;
pointCoord(120, 0) = 1.503949591308661e-02; pointCoord(120, 1) =1.037719441021005e-02; pointCoord(120, 2) =1.281223507749970e-04;
pointCoord(121, 0) = 2.650394959130866e-01; pointCoord(121, 1) =1.037719441021005e-02; pointCoord(121, 2) =1.281223507749970e-04;
pointCoord(122, 0) = 5.150394959130866e-01; pointCoord(122, 1) =1.037719441021005e-02; pointCoord(122, 2) =1.281223507749970e-04;
pointCoord(123, 0) = 7.650394959130866e-01; pointCoord(123, 1) =1.037719441021005e-02; pointCoord(123, 2) =1.281223507749970e-04;
pointCoord(124, 0) = 1.503949591308661e-02; pointCoord(124, 1) =2.603771944102101e-01; pointCoord(124, 2) =1.281223507749970e-04;
pointCoord(125, 0) = 2.650394959130866e-01; pointCoord(125, 1) =2.603771944102101e-01; pointCoord(125, 2) =1.281223507749970e-04;
pointCoord(126, 0) = 5.150394959130866e-01; pointCoord(126, 1) =2.603771944102101e-01; pointCoord(126, 2) =1.281223507749970e-04;
pointCoord(127, 0) = 1.503949591308661e-02; pointCoord(127, 1) =5.103771944102100e-01; pointCoord(127, 2) =1.281223507749970e-04;
pointCoord(128, 0) = 2.650394959130866e-01; pointCoord(128, 1) =5.103771944102100e-01; pointCoord(128, 2) =1.281223507749970e-04;
pointCoord(129, 0) = 1.503949591308661e-02; pointCoord(129, 1) =7.603771944102100e-01; pointCoord(129, 2) =1.281223507749970e-04;
pointCoord(130, 0) = 1.938699242049790e-02; pointCoord(130, 1) =6.029697902798763e-03; pointCoord(130, 2) =1.108205956969521e-04;
pointCoord(131, 0) = 2.693869924204979e-01; pointCoord(131, 1) =6.029697902798763e-03; pointCoord(131, 2) =1.108205956969521e-04;
pointCoord(132, 0) = 5.193869924204979e-01; pointCoord(132, 1) =6.029697902798763e-03; pointCoord(132, 2) =1.108205956969521e-04;
pointCoord(133, 0) = 7.693869924204979e-01; pointCoord(133, 1) =6.029697902798763e-03; pointCoord(133, 2) =1.108205956969521e-04;
pointCoord(134, 0) = 1.938699242049790e-02; pointCoord(134, 1) =2.560296979027988e-01; pointCoord(134, 2) =1.108205956969521e-04;
pointCoord(135, 0) = 2.693869924204979e-01; pointCoord(135, 1) =2.560296979027988e-01; pointCoord(135, 2) =1.108205956969521e-04;
pointCoord(136, 0) = 5.193869924204979e-01; pointCoord(136, 1) =2.560296979027988e-01; pointCoord(136, 2) =1.108205956969521e-04;
pointCoord(137, 0) = 1.938699242049790e-02; pointCoord(137, 1) =5.060296979027987e-01; pointCoord(137, 2) =1.108205956969521e-04;
pointCoord(138, 0) = 2.693869924204979e-01; pointCoord(138, 1) =5.060296979027987e-01; pointCoord(138, 2) =1.108205956969521e-04;
pointCoord(139, 0) = 1.938699242049790e-02; pointCoord(139, 1) =7.560296979027987e-01; pointCoord(139, 2) =1.108205956969521e-04;
pointCoord(140, 0) = 2.283265773533521e-02; pointCoord(140, 1) =2.584032587961448e-03; pointCoord(140, 2) =7.855873961758663e-05;
pointCoord(141, 0) = 2.728326577353352e-01; pointCoord(141, 1) =2.584032587961448e-03; pointCoord(141, 2) =7.855873961758663e-05;
pointCoord(142, 0) = 5.228326577353352e-01; pointCoord(142, 1) =2.584032587961448e-03; pointCoord(142, 2) =7.855873961758663e-05;
pointCoord(143, 0) = 7.728326577353352e-01; pointCoord(143, 1) =2.584032587961448e-03; pointCoord(143, 2) =7.855873961758663e-05;
pointCoord(144, 0) = 2.283265773533521e-02; pointCoord(144, 1) =2.525840325879615e-01; pointCoord(144, 2) =7.855873961758663e-05;
pointCoord(145, 0) = 2.728326577353352e-01; pointCoord(145, 1) =2.525840325879615e-01; pointCoord(145, 2) =7.855873961758663e-05;
pointCoord(146, 0) = 5.228326577353352e-01; pointCoord(146, 1) =2.525840325879615e-01; pointCoord(146, 2) =7.855873961758663e-05;
pointCoord(147, 0) = 2.283265773533521e-02; pointCoord(147, 1) =5.025840325879615e-01; pointCoord(147, 2) =7.855873961758663e-05;
pointCoord(148, 0) = 2.728326577353352e-01; pointCoord(148, 1) =5.025840325879615e-01; pointCoord(148, 2) =7.855873961758663e-05;
pointCoord(149, 0) = 2.283265773533521e-02; pointCoord(149, 1) =7.525840325879615e-01; pointCoord(149, 2) =7.855873961758663e-05;
pointCoord(150, 0) = 2.491204011324876e-02; pointCoord(150, 1) =5.046502100478956e-04; pointCoord(150, 2) =3.576018181520091e-05;
pointCoord(151, 0) = 2.749120401132488e-01; pointCoord(151, 1) =5.046502100478956e-04; pointCoord(151, 2) =3.576018181520091e-05;
pointCoord(152, 0) = 5.249120401132488e-01; pointCoord(152, 1) =5.046502100478956e-04; pointCoord(152, 2) =3.576018181520091e-05;
pointCoord(153, 0) = 7.749120401132488e-01; pointCoord(153, 1) =5.046502100478956e-04; pointCoord(153, 2) =3.576018181520091e-05;
pointCoord(154, 0) = 2.491204011324876e-02; pointCoord(154, 1) =2.505046502100479e-01; pointCoord(154, 2) =3.576018181520091e-05;
pointCoord(155, 0) = 2.749120401132488e-01; pointCoord(155, 1) =2.505046502100479e-01; pointCoord(155, 2) =3.576018181520091e-05;
pointCoord(156, 0) = 5.249120401132488e-01; pointCoord(156, 1) =2.505046502100479e-01; pointCoord(156, 2) =3.576018181520091e-05;
pointCoord(157, 0) = 2.491204011324876e-02; pointCoord(157, 1) =5.005046502100479e-01; pointCoord(157, 2) =3.576018181520091e-05;
pointCoord(158, 0) = 2.749120401132488e-01; pointCoord(158, 1) =5.005046502100479e-01; pointCoord(158, 2) =3.576018181520091e-05;
pointCoord(159, 0) = 2.491204011324876e-02; pointCoord(159, 1) =7.505046502100479e-01; pointCoord(159, 2) =3.576018181520091e-05;
pointCoord(160, 0) = 1.177573505593169e-03; pointCoord(160, 1) =5.813087525486571e-02; pointCoord(160, 2) =1.177126830861867e-04;
pointCoord(161, 0) = 2.511775735055932e-01; pointCoord(161, 1) =5.813087525486571e-02; pointCoord(161, 2) =1.177126830861867e-04;
pointCoord(162, 0) = 5.011775735055932e-01; pointCoord(162, 1) =5.813087525486571e-02; pointCoord(162, 2) =1.177126830861867e-04;
pointCoord(163, 0) = 7.511775735055932e-01; pointCoord(163, 1) =5.813087525486571e-02; pointCoord(163, 2) =1.177126830861867e-04;
pointCoord(164, 0) = 1.177573505593169e-03; pointCoord(164, 1) =3.081308752548657e-01; pointCoord(164, 2) =1.177126830861867e-04;
pointCoord(165, 0) = 2.511775735055932e-01; pointCoord(165, 1) =3.081308752548657e-01; pointCoord(165, 2) =1.177126830861867e-04;
pointCoord(166, 0) = 5.011775735055932e-01; pointCoord(166, 1) =3.081308752548657e-01; pointCoord(166, 2) =1.177126830861867e-04;
pointCoord(167, 0) = 1.177573505593169e-03; pointCoord(167, 1) =5.581308752548657e-01; pointCoord(167, 2) =1.177126830861867e-04;
pointCoord(168, 0) = 2.511775735055932e-01; pointCoord(168, 1) =5.581308752548657e-01; pointCoord(168, 2) =1.177126830861867e-04;
pointCoord(169, 0) = 1.177573505593169e-03; pointCoord(169, 1) =8.081308752548657e-01; pointCoord(169, 2) =1.177126830861867e-04;
pointCoord(170, 0) = 6.029697902798763e-03; pointCoord(170, 1) =5.327875085766011e-02; pointCoord(170, 2) =2.585937640933464e-04;
pointCoord(171, 0) = 2.560296979027988e-01; pointCoord(171, 1) =5.327875085766011e-02; pointCoord(171, 2) =2.585937640933464e-04;
pointCoord(172, 0) = 5.060296979027987e-01; pointCoord(172, 1) =5.327875085766011e-02; pointCoord(172, 2) =2.585937640933464e-04;
pointCoord(173, 0) = 7.560296979027987e-01; pointCoord(173, 1) =5.327875085766011e-02; pointCoord(173, 2) =2.585937640933464e-04;
pointCoord(174, 0) = 6.029697902798763e-03; pointCoord(174, 1) =3.032787508576601e-01; pointCoord(174, 2) =2.585937640933464e-04;
pointCoord(175, 0) = 2.560296979027988e-01; pointCoord(175, 1) =3.032787508576601e-01; pointCoord(175, 2) =2.585937640933464e-04;
pointCoord(176, 0) = 5.060296979027987e-01; pointCoord(176, 1) =3.032787508576601e-01; pointCoord(176, 2) =2.585937640933464e-04;
pointCoord(177, 0) = 6.029697902798763e-03; pointCoord(177, 1) =5.532787508576601e-01; pointCoord(177, 2) =2.585937640933464e-04;
pointCoord(178, 0) = 2.560296979027988e-01; pointCoord(178, 1) =5.532787508576601e-01; pointCoord(178, 2) =2.585937640933464e-04;
pointCoord(179, 0) = 6.029697902798763e-03; pointCoord(179, 1) =8.032787508576601e-01; pointCoord(179, 2) =2.585937640933464e-04;
pointCoord(180, 0) = 1.406996837748790e-02; pointCoord(180, 1) =4.523848038297097e-02; pointCoord(180, 2) =3.647909210336454e-04;
pointCoord(181, 0) = 2.640699683774879e-01; pointCoord(181, 1) =4.523848038297097e-02; pointCoord(181, 2) =3.647909210336454e-04;
pointCoord(182, 0) = 5.140699683774879e-01; pointCoord(182, 1) =4.523848038297097e-02; pointCoord(182, 2) =3.647909210336454e-04;
pointCoord(183, 0) = 7.640699683774879e-01; pointCoord(183, 1) =4.523848038297097e-02; pointCoord(183, 2) =3.647909210336454e-04;
pointCoord(184, 0) = 1.406996837748790e-02; pointCoord(184, 1) =2.952384803829710e-01; pointCoord(184, 2) =3.647909210336454e-04;
pointCoord(185, 0) = 2.640699683774879e-01; pointCoord(185, 1) =2.952384803829710e-01; pointCoord(185, 2) =3.647909210336454e-04;
pointCoord(186, 0) = 5.140699683774879e-01; pointCoord(186, 1) =2.952384803829710e-01; pointCoord(186, 2) =3.647909210336454e-04;
pointCoord(187, 0) = 1.406996837748790e-02; pointCoord(187, 1) =5.452384803829710e-01; pointCoord(187, 2) =3.647909210336454e-04;
pointCoord(188, 0) = 2.640699683774879e-01; pointCoord(188, 1) =5.452384803829710e-01; pointCoord(188, 2) =3.647909210336454e-04;
pointCoord(189, 0) = 1.406996837748790e-02; pointCoord(189, 1) =7.952384803829710e-01; pointCoord(189, 2) =3.647909210336454e-04;
pointCoord(190, 0) = 2.421461233255627e-02; pointCoord(190, 1) =3.509383642790261e-02; pointCoord(190, 2) =4.217435400908277e-04;
pointCoord(191, 0) = 2.742146123325563e-01; pointCoord(191, 1) =3.509383642790261e-02; pointCoord(191, 2) =4.217435400908277e-04;
pointCoord(192, 0) = 5.242146123325563e-01; pointCoord(192, 1) =3.509383642790261e-02; pointCoord(192, 2) =4.217435400908277e-04;
pointCoord(193, 0) = 7.742146123325563e-01; pointCoord(193, 1) =3.509383642790261e-02; pointCoord(193, 2) =4.217435400908277e-04;
pointCoord(194, 0) = 2.421461233255627e-02; pointCoord(194, 1) =2.850938364279026e-01; pointCoord(194, 2) =4.217435400908277e-04;
pointCoord(195, 0) = 2.742146123325563e-01; pointCoord(195, 1) =2.850938364279026e-01; pointCoord(195, 2) =4.217435400908277e-04;
pointCoord(196, 0) = 5.242146123325563e-01; pointCoord(196, 1) =2.850938364279026e-01; pointCoord(196, 2) =4.217435400908277e-04;
pointCoord(197, 0) = 2.421461233255627e-02; pointCoord(197, 1) =5.350938364279026e-01; pointCoord(197, 2) =4.217435400908277e-04;
pointCoord(198, 0) = 2.742146123325563e-01; pointCoord(198, 1) =5.350938364279026e-01; pointCoord(198, 2) =4.217435400908277e-04;
pointCoord(199, 0) = 2.421461233255627e-02; pointCoord(199, 1) =7.850938364279026e-01; pointCoord(199, 2) =4.217435400908277e-04;
pointCoord(200, 0) = 3.509383642790261e-02; pointCoord(200, 1) =2.421461233255627e-02; pointCoord(200, 2) =4.217435400908277e-04;
pointCoord(201, 0) = 2.850938364279026e-01; pointCoord(201, 1) =2.421461233255627e-02; pointCoord(201, 2) =4.217435400908277e-04;
pointCoord(202, 0) = 5.350938364279026e-01; pointCoord(202, 1) =2.421461233255627e-02; pointCoord(202, 2) =4.217435400908277e-04;
pointCoord(203, 0) = 7.850938364279026e-01; pointCoord(203, 1) =2.421461233255627e-02; pointCoord(203, 2) =4.217435400908277e-04;
pointCoord(204, 0) = 3.509383642790261e-02; pointCoord(204, 1) =2.742146123325563e-01; pointCoord(204, 2) =4.217435400908277e-04;
pointCoord(205, 0) = 2.850938364279026e-01; pointCoord(205, 1) =2.742146123325563e-01; pointCoord(205, 2) =4.217435400908277e-04;
pointCoord(206, 0) = 5.350938364279026e-01; pointCoord(206, 1) =2.742146123325563e-01; pointCoord(206, 2) =4.217435400908277e-04;
pointCoord(207, 0) = 3.509383642790261e-02; pointCoord(207, 1) =5.242146123325563e-01; pointCoord(207, 2) =4.217435400908277e-04;
pointCoord(208, 0) = 2.850938364279026e-01; pointCoord(208, 1) =5.242146123325563e-01; pointCoord(208, 2) =4.217435400908277e-04;
pointCoord(209, 0) = 3.509383642790261e-02; pointCoord(209, 1) =7.742146123325563e-01; pointCoord(209, 2) =4.217435400908277e-04;
pointCoord(210, 0) = 4.523848038297097e-02; pointCoord(210, 1) =1.406996837748790e-02; pointCoord(210, 2) =3.647909210336454e-04;
pointCoord(211, 0) = 2.952384803829710e-01; pointCoord(211, 1) =1.406996837748790e-02; pointCoord(211, 2) =3.647909210336454e-04;
pointCoord(212, 0) = 5.452384803829710e-01; pointCoord(212, 1) =1.406996837748790e-02; pointCoord(212, 2) =3.647909210336454e-04;
pointCoord(213, 0) = 7.952384803829710e-01; pointCoord(213, 1) =1.406996837748790e-02; pointCoord(213, 2) =3.647909210336454e-04;
pointCoord(214, 0) = 4.523848038297097e-02; pointCoord(214, 1) =2.640699683774879e-01; pointCoord(214, 2) =3.647909210336454e-04;
pointCoord(215, 0) = 2.952384803829710e-01; pointCoord(215, 1) =2.640699683774879e-01; pointCoord(215, 2) =3.647909210336454e-04;
pointCoord(216, 0) = 5.452384803829710e-01; pointCoord(216, 1) =2.640699683774879e-01; pointCoord(216, 2) =3.647909210336454e-04;
pointCoord(217, 0) = 4.523848038297097e-02; pointCoord(217, 1) =5.140699683774879e-01; pointCoord(217, 2) =3.647909210336454e-04;
pointCoord(218, 0) = 2.952384803829710e-01; pointCoord(218, 1) =5.140699683774879e-01; pointCoord(218, 2) =3.647909210336454e-04;
pointCoord(219, 0) = 4.523848038297097e-02; pointCoord(219, 1) =7.640699683774879e-01; pointCoord(219, 2) =3.647909210336454e-04;
pointCoord(220, 0) = 5.327875085766011e-02; pointCoord(220, 1) =6.029697902798763e-03; pointCoord(220, 2) =2.585937640933464e-04;
pointCoord(221, 0) = 3.032787508576601e-01; pointCoord(221, 1) =6.029697902798763e-03; pointCoord(221, 2) =2.585937640933464e-04;
pointCoord(222, 0) = 5.532787508576601e-01; pointCoord(222, 1) =6.029697902798763e-03; pointCoord(222, 2) =2.585937640933464e-04;
pointCoord(223, 0) = 8.032787508576601e-01; pointCoord(223, 1) =6.029697902798763e-03; pointCoord(223, 2) =2.585937640933464e-04;
pointCoord(224, 0) = 5.327875085766011e-02; pointCoord(224, 1) =2.560296979027988e-01; pointCoord(224, 2) =2.585937640933464e-04;
pointCoord(225, 0) = 3.032787508576601e-01; pointCoord(225, 1) =2.560296979027988e-01; pointCoord(225, 2) =2.585937640933464e-04;
pointCoord(226, 0) = 5.532787508576601e-01; pointCoord(226, 1) =2.560296979027988e-01; pointCoord(226, 2) =2.585937640933464e-04;
pointCoord(227, 0) = 5.327875085766011e-02; pointCoord(227, 1) =5.060296979027987e-01; pointCoord(227, 2) =2.585937640933464e-04;
pointCoord(228, 0) = 3.032787508576601e-01; pointCoord(228, 1) =5.060296979027987e-01; pointCoord(228, 2) =2.585937640933464e-04;
pointCoord(229, 0) = 5.327875085766011e-02; pointCoord(229, 1) =7.560296979027987e-01; pointCoord(229, 2) =2.585937640933464e-04;
pointCoord(230, 0) = 5.813087525486571e-02; pointCoord(230, 1) =1.177573505593169e-03; pointCoord(230, 2) =1.177126830861867e-04;
pointCoord(231, 0) = 3.081308752548657e-01; pointCoord(231, 1) =1.177573505593169e-03; pointCoord(231, 2) =1.177126830861867e-04;
pointCoord(232, 0) = 5.581308752548657e-01; pointCoord(232, 1) =1.177573505593169e-03; pointCoord(232, 2) =1.177126830861867e-04;
pointCoord(233, 0) = 8.081308752548657e-01; pointCoord(233, 1) =1.177573505593169e-03; pointCoord(233, 2) =1.177126830861867e-04;
pointCoord(234, 0) = 5.813087525486571e-02; pointCoord(234, 1) =2.511775735055932e-01; pointCoord(234, 2) =1.177126830861867e-04;
pointCoord(235, 0) = 3.081308752548657e-01; pointCoord(235, 1) =2.511775735055932e-01; pointCoord(235, 2) =1.177126830861867e-04;
pointCoord(236, 0) = 5.581308752548657e-01; pointCoord(236, 1) =2.511775735055932e-01; pointCoord(236, 2) =1.177126830861867e-04;
pointCoord(237, 0) = 5.813087525486571e-02; pointCoord(237, 1) =5.011775735055932e-01; pointCoord(237, 2) =1.177126830861867e-04;
pointCoord(238, 0) = 3.081308752548657e-01; pointCoord(238, 1) =5.011775735055932e-01; pointCoord(238, 2) =1.177126830861867e-04;
pointCoord(239, 0) = 5.813087525486571e-02; pointCoord(239, 1) =7.511775735055932e-01; pointCoord(239, 2) =1.177126830861867e-04;
pointCoord(240, 0) = 2.026620470352396e-03; pointCoord(240, 1) =1.000440492176914e-01; pointCoord(240, 2) =2.342135820693354e-04;
pointCoord(241, 0) = 2.520266204703524e-01; pointCoord(241, 1) =1.000440492176914e-01; pointCoord(241, 2) =2.342135820693354e-04;
pointCoord(242, 0) = 5.020266204703524e-01; pointCoord(242, 1) =1.000440492176914e-01; pointCoord(242, 2) =2.342135820693354e-04;
pointCoord(243, 0) = 7.520266204703524e-01; pointCoord(243, 1) =1.000440492176914e-01; pointCoord(243, 2) =2.342135820693354e-04;
pointCoord(244, 0) = 2.026620470352396e-03; pointCoord(244, 1) =3.500440492176914e-01; pointCoord(244, 2) =2.342135820693354e-04;
pointCoord(245, 0) = 2.520266204703524e-01; pointCoord(245, 1) =3.500440492176914e-01; pointCoord(245, 2) =2.342135820693354e-04;
pointCoord(246, 0) = 5.020266204703524e-01; pointCoord(246, 1) =3.500440492176914e-01; pointCoord(246, 2) =2.342135820693354e-04;
pointCoord(247, 0) = 2.026620470352396e-03; pointCoord(247, 1) =6.000440492176914e-01; pointCoord(247, 2) =2.342135820693354e-04;
pointCoord(248, 0) = 2.520266204703524e-01; pointCoord(248, 1) =6.000440492176914e-01; pointCoord(248, 2) =2.342135820693354e-04;
pointCoord(249, 0) = 2.026620470352396e-03; pointCoord(249, 1) =8.500440492176914e-01; pointCoord(249, 2) =2.342135820693354e-04;
pointCoord(250, 0) = 1.037719441021005e-02; pointCoord(250, 1) =9.169347527783372e-02; pointCoord(250, 2) =5.145254547018528e-04;
pointCoord(251, 0) = 2.603771944102101e-01; pointCoord(251, 1) =9.169347527783372e-02; pointCoord(251, 2) =5.145254547018528e-04;
pointCoord(252, 0) = 5.103771944102100e-01; pointCoord(252, 1) =9.169347527783372e-02; pointCoord(252, 2) =5.145254547018528e-04;
pointCoord(253, 0) = 7.603771944102100e-01; pointCoord(253, 1) =9.169347527783372e-02; pointCoord(253, 2) =5.145254547018528e-04;
pointCoord(254, 0) = 1.037719441021005e-02; pointCoord(254, 1) =3.416934752778337e-01; pointCoord(254, 2) =5.145254547018528e-04;
pointCoord(255, 0) = 2.603771944102101e-01; pointCoord(255, 1) =3.416934752778337e-01; pointCoord(255, 2) =5.145254547018528e-04;
pointCoord(256, 0) = 5.103771944102100e-01; pointCoord(256, 1) =3.416934752778337e-01; pointCoord(256, 2) =5.145254547018528e-04;
pointCoord(257, 0) = 1.037719441021005e-02; pointCoord(257, 1) =5.916934752778338e-01; pointCoord(257, 2) =5.145254547018528e-04;
pointCoord(258, 0) = 2.603771944102101e-01; pointCoord(258, 1) =5.916934752778338e-01; pointCoord(258, 2) =5.145254547018528e-04;
pointCoord(259, 0) = 1.037719441021005e-02; pointCoord(259, 1) =8.416934752778338e-01; pointCoord(259, 2) =5.145254547018528e-04;
pointCoord(260, 0) = 2.421461233255627e-02; pointCoord(260, 1) =7.785605735548751e-02; pointCoord(260, 2) =7.258265301718213e-04;
pointCoord(261, 0) = 2.742146123325563e-01; pointCoord(261, 1) =7.785605735548751e-02; pointCoord(261, 2) =7.258265301718213e-04;
pointCoord(262, 0) = 5.242146123325563e-01; pointCoord(262, 1) =7.785605735548751e-02; pointCoord(262, 2) =7.258265301718213e-04;
pointCoord(263, 0) = 7.742146123325563e-01; pointCoord(263, 1) =7.785605735548751e-02; pointCoord(263, 2) =7.258265301718213e-04;
pointCoord(264, 0) = 2.421461233255627e-02; pointCoord(264, 1) =3.278560573554875e-01; pointCoord(264, 2) =7.258265301718213e-04;
pointCoord(265, 0) = 2.742146123325563e-01; pointCoord(265, 1) =3.278560573554875e-01; pointCoord(265, 2) =7.258265301718213e-04;
pointCoord(266, 0) = 5.242146123325563e-01; pointCoord(266, 1) =3.278560573554875e-01; pointCoord(266, 2) =7.258265301718213e-04;
pointCoord(267, 0) = 2.421461233255627e-02; pointCoord(267, 1) =5.778560573554875e-01; pointCoord(267, 2) =7.258265301718213e-04;
pointCoord(268, 0) = 2.742146123325563e-01; pointCoord(268, 1) =5.778560573554875e-01; pointCoord(268, 2) =7.258265301718213e-04;
pointCoord(269, 0) = 2.421461233255627e-02; pointCoord(269, 1) =8.278560573554875e-01; pointCoord(269, 2) =7.258265301718213e-04;
pointCoord(270, 0) = 4.167368644226296e-02; pointCoord(270, 1) =6.039698324578081e-02; pointCoord(270, 2) =8.391454739584170e-04;
pointCoord(271, 0) = 2.916736864422629e-01; pointCoord(271, 1) =6.039698324578081e-02; pointCoord(271, 2) =8.391454739584170e-04;
pointCoord(272, 0) = 5.416736864422630e-01; pointCoord(272, 1) =6.039698324578081e-02; pointCoord(272, 2) =8.391454739584170e-04;
pointCoord(273, 0) = 7.916736864422630e-01; pointCoord(273, 1) =6.039698324578081e-02; pointCoord(273, 2) =8.391454739584170e-04;
pointCoord(274, 0) = 4.167368644226296e-02; pointCoord(274, 1) =3.103969832457808e-01; pointCoord(274, 2) =8.391454739584170e-04;
pointCoord(275, 0) = 2.916736864422629e-01; pointCoord(275, 1) =3.103969832457808e-01; pointCoord(275, 2) =8.391454739584170e-04;
pointCoord(276, 0) = 5.416736864422630e-01; pointCoord(276, 1) =3.103969832457808e-01; pointCoord(276, 2) =8.391454739584170e-04;
pointCoord(277, 0) = 4.167368644226296e-02; pointCoord(277, 1) =5.603969832457808e-01; pointCoord(277, 2) =8.391454739584170e-04;
pointCoord(278, 0) = 2.916736864422629e-01; pointCoord(278, 1) =5.603969832457808e-01; pointCoord(278, 2) =8.391454739584170e-04;
pointCoord(279, 0) = 4.167368644226296e-02; pointCoord(279, 1) =8.103969832457808e-01; pointCoord(279, 2) =8.391454739584170e-04;
pointCoord(280, 0) = 6.039698324578081e-02; pointCoord(280, 1) =4.167368644226296e-02; pointCoord(280, 2) =8.391454739584170e-04;
pointCoord(281, 0) = 3.103969832457808e-01; pointCoord(281, 1) =4.167368644226296e-02; pointCoord(281, 2) =8.391454739584170e-04;
pointCoord(282, 0) = 5.603969832457808e-01; pointCoord(282, 1) =4.167368644226296e-02; pointCoord(282, 2) =8.391454739584170e-04;
pointCoord(283, 0) = 8.103969832457808e-01; pointCoord(283, 1) =4.167368644226296e-02; pointCoord(283, 2) =8.391454739584170e-04;
pointCoord(284, 0) = 6.039698324578081e-02; pointCoord(284, 1) =2.916736864422629e-01; pointCoord(284, 2) =8.391454739584170e-04;
pointCoord(285, 0) = 3.103969832457808e-01; pointCoord(285, 1) =2.916736864422629e-01; pointCoord(285, 2) =8.391454739584170e-04;
pointCoord(286, 0) = 5.603969832457808e-01; pointCoord(286, 1) =2.916736864422629e-01; pointCoord(286, 2) =8.391454739584170e-04;
pointCoord(287, 0) = 6.039698324578081e-02; pointCoord(287, 1) =5.416736864422630e-01; pointCoord(287, 2) =8.391454739584170e-04;
pointCoord(288, 0) = 3.103969832457808e-01; pointCoord(288, 1) =5.416736864422630e-01; pointCoord(288, 2) =8.391454739584170e-04;
pointCoord(289, 0) = 6.039698324578081e-02; pointCoord(289, 1) =7.916736864422630e-01; pointCoord(289, 2) =8.391454739584170e-04;
pointCoord(290, 0) = 7.785605735548751e-02; pointCoord(290, 1) =2.421461233255627e-02; pointCoord(290, 2) =7.258265301718213e-04;
pointCoord(291, 0) = 3.278560573554875e-01; pointCoord(291, 1) =2.421461233255627e-02; pointCoord(291, 2) =7.258265301718213e-04;
pointCoord(292, 0) = 5.778560573554875e-01; pointCoord(292, 1) =2.421461233255627e-02; pointCoord(292, 2) =7.258265301718213e-04;
pointCoord(293, 0) = 8.278560573554875e-01; pointCoord(293, 1) =2.421461233255627e-02; pointCoord(293, 2) =7.258265301718213e-04;
pointCoord(294, 0) = 7.785605735548751e-02; pointCoord(294, 1) =2.742146123325563e-01; pointCoord(294, 2) =7.258265301718213e-04;
pointCoord(295, 0) = 3.278560573554875e-01; pointCoord(295, 1) =2.742146123325563e-01; pointCoord(295, 2) =7.258265301718213e-04;
pointCoord(296, 0) = 5.778560573554875e-01; pointCoord(296, 1) =2.742146123325563e-01; pointCoord(296, 2) =7.258265301718213e-04;
pointCoord(297, 0) = 7.785605735548751e-02; pointCoord(297, 1) =5.242146123325563e-01; pointCoord(297, 2) =7.258265301718213e-04;
pointCoord(298, 0) = 3.278560573554875e-01; pointCoord(298, 1) =5.242146123325563e-01; pointCoord(298, 2) =7.258265301718213e-04;
pointCoord(299, 0) = 7.785605735548751e-02; pointCoord(299, 1) =7.742146123325563e-01; pointCoord(299, 2) =7.258265301718213e-04;
pointCoord(300, 0) = 9.169347527783372e-02; pointCoord(300, 1) =1.037719441021005e-02; pointCoord(300, 2) =5.145254547018528e-04;
pointCoord(301, 0) = 3.416934752778337e-01; pointCoord(301, 1) =1.037719441021005e-02; pointCoord(301, 2) =5.145254547018528e-04;
pointCoord(302, 0) = 5.916934752778338e-01; pointCoord(302, 1) =1.037719441021005e-02; pointCoord(302, 2) =5.145254547018528e-04;
pointCoord(303, 0) = 8.416934752778338e-01; pointCoord(303, 1) =1.037719441021005e-02; pointCoord(303, 2) =5.145254547018528e-04;
pointCoord(304, 0) = 9.169347527783372e-02; pointCoord(304, 1) =2.603771944102101e-01; pointCoord(304, 2) =5.145254547018528e-04;
pointCoord(305, 0) = 3.416934752778337e-01; pointCoord(305, 1) =2.603771944102101e-01; pointCoord(305, 2) =5.145254547018528e-04;
pointCoord(306, 0) = 5.916934752778338e-01; pointCoord(306, 1) =2.603771944102101e-01; pointCoord(306, 2) =5.145254547018528e-04;
pointCoord(307, 0) = 9.169347527783372e-02; pointCoord(307, 1) =5.103771944102100e-01; pointCoord(307, 2) =5.145254547018528e-04;
pointCoord(308, 0) = 3.416934752778337e-01; pointCoord(308, 1) =5.103771944102100e-01; pointCoord(308, 2) =5.145254547018528e-04;
pointCoord(309, 0) = 9.169347527783372e-02; pointCoord(309, 1) =7.603771944102100e-01; pointCoord(309, 2) =5.145254547018528e-04;
pointCoord(310, 0) = 1.000440492176914e-01; pointCoord(310, 1) =2.026620470352396e-03; pointCoord(310, 2) =2.342135820693354e-04;
pointCoord(311, 0) = 3.500440492176914e-01; pointCoord(311, 1) =2.026620470352396e-03; pointCoord(311, 2) =2.342135820693354e-04;
pointCoord(312, 0) = 6.000440492176914e-01; pointCoord(312, 1) =2.026620470352396e-03; pointCoord(312, 2) =2.342135820693354e-04;
pointCoord(313, 0) = 8.500440492176914e-01; pointCoord(313, 1) =2.026620470352396e-03; pointCoord(313, 2) =2.342135820693354e-04;
pointCoord(314, 0) = 1.000440492176914e-01; pointCoord(314, 1) =2.520266204703524e-01; pointCoord(314, 2) =2.342135820693354e-04;
pointCoord(315, 0) = 3.500440492176914e-01; pointCoord(315, 1) =2.520266204703524e-01; pointCoord(315, 2) =2.342135820693354e-04;
pointCoord(316, 0) = 6.000440492176914e-01; pointCoord(316, 1) =2.520266204703524e-01; pointCoord(316, 2) =2.342135820693354e-04;
pointCoord(317, 0) = 1.000440492176914e-01; pointCoord(317, 1) =5.020266204703524e-01; pointCoord(317, 2) =2.342135820693354e-04;
pointCoord(318, 0) = 3.500440492176914e-01; pointCoord(318, 1) =5.020266204703524e-01; pointCoord(318, 2) =2.342135820693354e-04;
pointCoord(319, 0) = 1.000440492176914e-01; pointCoord(319, 1) =7.520266204703524e-01; pointCoord(319, 2) =2.342135820693354e-04;
pointCoord(320, 0) = 2.937147467455568e-03; pointCoord(320, 1) =1.449921828445007e-01; pointCoord(320, 2) =3.394418636751593e-04;
pointCoord(321, 0) = 2.529371474674555e-01; pointCoord(321, 1) =1.449921828445007e-01; pointCoord(321, 2) =3.394418636751593e-04;
pointCoord(322, 0) = 5.029371474674555e-01; pointCoord(322, 1) =1.449921828445007e-01; pointCoord(322, 2) =3.394418636751593e-04;
pointCoord(323, 0) = 7.529371474674555e-01; pointCoord(323, 1) =1.449921828445007e-01; pointCoord(323, 2) =3.394418636751593e-04;
pointCoord(324, 0) = 2.937147467455568e-03; pointCoord(324, 1) =3.949921828445007e-01; pointCoord(324, 2) =3.394418636751593e-04;
pointCoord(325, 0) = 2.529371474674555e-01; pointCoord(325, 1) =3.949921828445007e-01; pointCoord(325, 2) =3.394418636751593e-04;
pointCoord(326, 0) = 5.029371474674555e-01; pointCoord(326, 1) =3.949921828445007e-01; pointCoord(326, 2) =3.394418636751593e-04;
pointCoord(327, 0) = 2.937147467455568e-03; pointCoord(327, 1) =6.449921828445007e-01; pointCoord(327, 2) =3.394418636751593e-04;
pointCoord(328, 0) = 2.529371474674555e-01; pointCoord(328, 1) =6.449921828445007e-01; pointCoord(328, 2) =3.394418636751593e-04;
pointCoord(329, 0) = 2.937147467455568e-03; pointCoord(329, 1) =8.949921828445007e-01; pointCoord(329, 2) =3.394418636751593e-04;
pointCoord(330, 0) = 1.503949591308661e-02; pointCoord(330, 1) =1.328898343988696e-01; pointCoord(330, 2) =7.456932160347678e-04;
pointCoord(331, 0) = 2.650394959130866e-01; pointCoord(331, 1) =1.328898343988696e-01; pointCoord(331, 2) =7.456932160347678e-04;
pointCoord(332, 0) = 5.150394959130866e-01; pointCoord(332, 1) =1.328898343988696e-01; pointCoord(332, 2) =7.456932160347678e-04;
pointCoord(333, 0) = 7.650394959130866e-01; pointCoord(333, 1) =1.328898343988696e-01; pointCoord(333, 2) =7.456932160347678e-04;
pointCoord(334, 0) = 1.503949591308661e-02; pointCoord(334, 1) =3.828898343988696e-01; pointCoord(334, 2) =7.456932160347678e-04;
pointCoord(335, 0) = 2.650394959130866e-01; pointCoord(335, 1) =3.828898343988696e-01; pointCoord(335, 2) =7.456932160347678e-04;
pointCoord(336, 0) = 5.150394959130866e-01; pointCoord(336, 1) =3.828898343988696e-01; pointCoord(336, 2) =7.456932160347678e-04;
pointCoord(337, 0) = 1.503949591308661e-02; pointCoord(337, 1) =6.328898343988696e-01; pointCoord(337, 2) =7.456932160347678e-04;
pointCoord(338, 0) = 2.650394959130866e-01; pointCoord(338, 1) =6.328898343988696e-01; pointCoord(338, 2) =7.456932160347678e-04;
pointCoord(339, 0) = 1.503949591308661e-02; pointCoord(339, 1) =8.828898343988696e-01; pointCoord(339, 2) =7.456932160347678e-04;
pointCoord(340, 0) = 3.509383642790261e-02; pointCoord(340, 1) =1.128354938840536e-01; pointCoord(340, 2) =1.051928363545806e-03;
pointCoord(341, 0) = 2.850938364279026e-01; pointCoord(341, 1) =1.128354938840536e-01; pointCoord(341, 2) =1.051928363545806e-03;
pointCoord(342, 0) = 5.350938364279026e-01; pointCoord(342, 1) =1.128354938840536e-01; pointCoord(342, 2) =1.051928363545806e-03;
pointCoord(343, 0) = 7.850938364279026e-01; pointCoord(343, 1) =1.128354938840536e-01; pointCoord(343, 2) =1.051928363545806e-03;
pointCoord(344, 0) = 3.509383642790261e-02; pointCoord(344, 1) =3.628354938840536e-01; pointCoord(344, 2) =1.051928363545806e-03;
pointCoord(345, 0) = 2.850938364279026e-01; pointCoord(345, 1) =3.628354938840536e-01; pointCoord(345, 2) =1.051928363545806e-03;
pointCoord(346, 0) = 5.350938364279026e-01; pointCoord(346, 1) =3.628354938840536e-01; pointCoord(346, 2) =1.051928363545806e-03;
pointCoord(347, 0) = 3.509383642790261e-02; pointCoord(347, 1) =6.128354938840536e-01; pointCoord(347, 2) =1.051928363545806e-03;
pointCoord(348, 0) = 2.850938364279026e-01; pointCoord(348, 1) =6.128354938840536e-01; pointCoord(348, 2) =1.051928363545806e-03;
pointCoord(349, 0) = 3.509383642790261e-02; pointCoord(349, 1) =8.628354938840536e-01; pointCoord(349, 2) =1.051928363545806e-03;
pointCoord(350, 0) = 6.039698324578081e-02; pointCoord(350, 1) =8.753234706617538e-02; pointCoord(350, 2) =1.216159631129748e-03;
pointCoord(351, 0) = 3.103969832457808e-01; pointCoord(351, 1) =8.753234706617538e-02; pointCoord(351, 2) =1.216159631129748e-03;
pointCoord(352, 0) = 5.603969832457808e-01; pointCoord(352, 1) =8.753234706617538e-02; pointCoord(352, 2) =1.216159631129748e-03;
pointCoord(353, 0) = 8.103969832457808e-01; pointCoord(353, 1) =8.753234706617538e-02; pointCoord(353, 2) =1.216159631129748e-03;
pointCoord(354, 0) = 6.039698324578081e-02; pointCoord(354, 1) =3.375323470661754e-01; pointCoord(354, 2) =1.216159631129748e-03;
pointCoord(355, 0) = 3.103969832457808e-01; pointCoord(355, 1) =3.375323470661754e-01; pointCoord(355, 2) =1.216159631129748e-03;
pointCoord(356, 0) = 5.603969832457808e-01; pointCoord(356, 1) =3.375323470661754e-01; pointCoord(356, 2) =1.216159631129748e-03;
pointCoord(357, 0) = 6.039698324578081e-02; pointCoord(357, 1) =5.875323470661754e-01; pointCoord(357, 2) =1.216159631129748e-03;
pointCoord(358, 0) = 3.103969832457808e-01; pointCoord(358, 1) =5.875323470661754e-01; pointCoord(358, 2) =1.216159631129748e-03;
pointCoord(359, 0) = 6.039698324578081e-02; pointCoord(359, 1) =8.375323470661754e-01; pointCoord(359, 2) =1.216159631129748e-03;
pointCoord(360, 0) = 8.753234706617538e-02; pointCoord(360, 1) =6.039698324578081e-02; pointCoord(360, 2) =1.216159631129748e-03;
pointCoord(361, 0) = 3.375323470661754e-01; pointCoord(361, 1) =6.039698324578081e-02; pointCoord(361, 2) =1.216159631129748e-03;
pointCoord(362, 0) = 5.875323470661754e-01; pointCoord(362, 1) =6.039698324578081e-02; pointCoord(362, 2) =1.216159631129748e-03;
pointCoord(363, 0) = 8.375323470661754e-01; pointCoord(363, 1) =6.039698324578081e-02; pointCoord(363, 2) =1.216159631129748e-03;
pointCoord(364, 0) = 8.753234706617538e-02; pointCoord(364, 1) =3.103969832457808e-01; pointCoord(364, 2) =1.216159631129748e-03;
pointCoord(365, 0) = 3.375323470661754e-01; pointCoord(365, 1) =3.103969832457808e-01; pointCoord(365, 2) =1.216159631129748e-03;
pointCoord(366, 0) = 5.875323470661754e-01; pointCoord(366, 1) =3.103969832457808e-01; pointCoord(366, 2) =1.216159631129748e-03;
pointCoord(367, 0) = 8.753234706617538e-02; pointCoord(367, 1) =5.603969832457808e-01; pointCoord(367, 2) =1.216159631129748e-03;
pointCoord(368, 0) = 3.375323470661754e-01; pointCoord(368, 1) =5.603969832457808e-01; pointCoord(368, 2) =1.216159631129748e-03;
pointCoord(369, 0) = 8.753234706617538e-02; pointCoord(369, 1) =8.103969832457808e-01; pointCoord(369, 2) =1.216159631129748e-03;
pointCoord(370, 0) = 1.128354938840536e-01; pointCoord(370, 1) =3.509383642790261e-02; pointCoord(370, 2) =1.051928363545806e-03;
pointCoord(371, 0) = 3.628354938840536e-01; pointCoord(371, 1) =3.509383642790261e-02; pointCoord(371, 2) =1.051928363545806e-03;
pointCoord(372, 0) = 6.128354938840536e-01; pointCoord(372, 1) =3.509383642790261e-02; pointCoord(372, 2) =1.051928363545806e-03;
pointCoord(373, 0) = 8.628354938840536e-01; pointCoord(373, 1) =3.509383642790261e-02; pointCoord(373, 2) =1.051928363545806e-03;
pointCoord(374, 0) = 1.128354938840536e-01; pointCoord(374, 1) =2.850938364279026e-01; pointCoord(374, 2) =1.051928363545806e-03;
pointCoord(375, 0) = 3.628354938840536e-01; pointCoord(375, 1) =2.850938364279026e-01; pointCoord(375, 2) =1.051928363545806e-03;
pointCoord(376, 0) = 6.128354938840536e-01; pointCoord(376, 1) =2.850938364279026e-01; pointCoord(376, 2) =1.051928363545806e-03;
pointCoord(377, 0) = 1.128354938840536e-01; pointCoord(377, 1) =5.350938364279026e-01; pointCoord(377, 2) =1.051928363545806e-03;
pointCoord(378, 0) = 3.628354938840536e-01; pointCoord(378, 1) =5.350938364279026e-01; pointCoord(378, 2) =1.051928363545806e-03;
pointCoord(379, 0) = 1.128354938840536e-01; pointCoord(379, 1) =7.850938364279026e-01; pointCoord(379, 2) =1.051928363545806e-03;
pointCoord(380, 0) = 1.328898343988696e-01; pointCoord(380, 1) =1.503949591308661e-02; pointCoord(380, 2) =7.456932160347678e-04;
pointCoord(381, 0) = 3.828898343988696e-01; pointCoord(381, 1) =1.503949591308661e-02; pointCoord(381, 2) =7.456932160347678e-04;
pointCoord(382, 0) = 6.328898343988696e-01; pointCoord(382, 1) =1.503949591308661e-02; pointCoord(382, 2) =7.456932160347678e-04;
pointCoord(383, 0) = 8.828898343988696e-01; pointCoord(383, 1) =1.503949591308661e-02; pointCoord(383, 2) =7.456932160347678e-04;
pointCoord(384, 0) = 1.328898343988696e-01; pointCoord(384, 1) =2.650394959130866e-01; pointCoord(384, 2) =7.456932160347678e-04;
pointCoord(385, 0) = 3.828898343988696e-01; pointCoord(385, 1) =2.650394959130866e-01; pointCoord(385, 2) =7.456932160347678e-04;
pointCoord(386, 0) = 6.328898343988696e-01; pointCoord(386, 1) =2.650394959130866e-01; pointCoord(386, 2) =7.456932160347678e-04;
pointCoord(387, 0) = 1.328898343988696e-01; pointCoord(387, 1) =5.150394959130866e-01; pointCoord(387, 2) =7.456932160347678e-04;
pointCoord(388, 0) = 3.828898343988696e-01; pointCoord(388, 1) =5.150394959130866e-01; pointCoord(388, 2) =7.456932160347678e-04;
pointCoord(389, 0) = 1.328898343988696e-01; pointCoord(389, 1) =7.650394959130866e-01; pointCoord(389, 2) =7.456932160347678e-04;
pointCoord(390, 0) = 1.449921828445007e-01; pointCoord(390, 1) =2.937147467455568e-03; pointCoord(390, 2) =3.394418636751593e-04;
pointCoord(391, 0) = 3.949921828445007e-01; pointCoord(391, 1) =2.937147467455568e-03; pointCoord(391, 2) =3.394418636751593e-04;
pointCoord(392, 0) = 6.449921828445007e-01; pointCoord(392, 1) =2.937147467455568e-03; pointCoord(392, 2) =3.394418636751593e-04;
pointCoord(393, 0) = 8.949921828445007e-01; pointCoord(393, 1) =2.937147467455568e-03; pointCoord(393, 2) =3.394418636751593e-04;
pointCoord(394, 0) = 1.449921828445007e-01; pointCoord(394, 1) =2.529371474674555e-01; pointCoord(394, 2) =3.394418636751593e-04;
pointCoord(395, 0) = 3.949921828445007e-01; pointCoord(395, 1) =2.529371474674555e-01; pointCoord(395, 2) =3.394418636751593e-04;
pointCoord(396, 0) = 6.449921828445007e-01; pointCoord(396, 1) =2.529371474674555e-01; pointCoord(396, 2) =3.394418636751593e-04;
pointCoord(397, 0) = 1.449921828445007e-01; pointCoord(397, 1) =5.029371474674555e-01; pointCoord(397, 2) =3.394418636751593e-04;
pointCoord(398, 0) = 3.949921828445007e-01; pointCoord(398, 1) =5.029371474674555e-01; pointCoord(398, 2) =3.394418636751593e-04;
pointCoord(399, 0) = 1.449921828445007e-01; pointCoord(399, 1) =7.529371474674555e-01; pointCoord(399, 2) =3.394418636751593e-04;
pointCoord(400, 0) = 3.786194432214795e-03; pointCoord(400, 1) =1.869053568073263e-01; pointCoord(400, 2) =3.784758260822833e-04;
pointCoord(401, 0) = 2.537861944322148e-01; pointCoord(401, 1) =1.869053568073263e-01; pointCoord(401, 2) =3.784758260822833e-04;
pointCoord(402, 0) = 5.037861944322148e-01; pointCoord(402, 1) =1.869053568073263e-01; pointCoord(402, 2) =3.784758260822833e-04;
pointCoord(403, 0) = 7.537861944322148e-01; pointCoord(403, 1) =1.869053568073263e-01; pointCoord(403, 2) =3.784758260822833e-04;
pointCoord(404, 0) = 3.786194432214795e-03; pointCoord(404, 1) =4.369053568073263e-01; pointCoord(404, 2) =3.784758260822833e-04;
pointCoord(405, 0) = 2.537861944322148e-01; pointCoord(405, 1) =4.369053568073263e-01; pointCoord(405, 2) =3.784758260822833e-04;
pointCoord(406, 0) = 5.037861944322148e-01; pointCoord(406, 1) =4.369053568073263e-01; pointCoord(406, 2) =3.784758260822833e-04;
pointCoord(407, 0) = 3.786194432214795e-03; pointCoord(407, 1) =6.869053568073263e-01; pointCoord(407, 2) =3.784758260822833e-04;
pointCoord(408, 0) = 2.537861944322148e-01; pointCoord(408, 1) =6.869053568073263e-01; pointCoord(408, 2) =3.784758260822833e-04;
pointCoord(409, 0) = 3.786194432214795e-03; pointCoord(409, 1) =9.369053568073263e-01; pointCoord(409, 2) =3.784758260822833e-04;
pointCoord(410, 0) = 1.938699242049790e-02; pointCoord(410, 1) =1.713045588190432e-01; pointCoord(410, 2) =8.314438675507627e-04;
pointCoord(411, 0) = 2.693869924204979e-01; pointCoord(411, 1) =1.713045588190432e-01; pointCoord(411, 2) =8.314438675507627e-04;
pointCoord(412, 0) = 5.193869924204979e-01; pointCoord(412, 1) =1.713045588190432e-01; pointCoord(412, 2) =8.314438675507627e-04;
pointCoord(413, 0) = 7.693869924204979e-01; pointCoord(413, 1) =1.713045588190432e-01; pointCoord(413, 2) =8.314438675507627e-04;
pointCoord(414, 0) = 1.938699242049790e-02; pointCoord(414, 1) =4.213045588190432e-01; pointCoord(414, 2) =8.314438675507627e-04;
pointCoord(415, 0) = 2.693869924204979e-01; pointCoord(415, 1) =4.213045588190432e-01; pointCoord(415, 2) =8.314438675507627e-04;
pointCoord(416, 0) = 5.193869924204979e-01; pointCoord(416, 1) =4.213045588190432e-01; pointCoord(416, 2) =8.314438675507627e-04;
pointCoord(417, 0) = 1.938699242049790e-02; pointCoord(417, 1) =6.713045588190432e-01; pointCoord(417, 2) =8.314438675507627e-04;
pointCoord(418, 0) = 2.693869924204979e-01; pointCoord(418, 1) =6.713045588190432e-01; pointCoord(418, 2) =8.314438675507627e-04;
pointCoord(419, 0) = 1.938699242049790e-02; pointCoord(419, 1) =9.213045588190432e-01; pointCoord(419, 2) =8.314438675507627e-04;
pointCoord(420, 0) = 4.523848038297097e-02; pointCoord(420, 1) =1.454530708565701e-01; pointCoord(420, 2) =1.172894386278138e-03;
pointCoord(421, 0) = 2.952384803829710e-01; pointCoord(421, 1) =1.454530708565701e-01; pointCoord(421, 2) =1.172894386278138e-03;
pointCoord(422, 0) = 5.452384803829710e-01; pointCoord(422, 1) =1.454530708565701e-01; pointCoord(422, 2) =1.172894386278138e-03;
pointCoord(423, 0) = 7.952384803829710e-01; pointCoord(423, 1) =1.454530708565701e-01; pointCoord(423, 2) =1.172894386278138e-03;
pointCoord(424, 0) = 4.523848038297097e-02; pointCoord(424, 1) =3.954530708565701e-01; pointCoord(424, 2) =1.172894386278138e-03;
pointCoord(425, 0) = 2.952384803829710e-01; pointCoord(425, 1) =3.954530708565701e-01; pointCoord(425, 2) =1.172894386278138e-03;
pointCoord(426, 0) = 5.452384803829710e-01; pointCoord(426, 1) =3.954530708565701e-01; pointCoord(426, 2) =1.172894386278138e-03;
pointCoord(427, 0) = 4.523848038297097e-02; pointCoord(427, 1) =6.454530708565701e-01; pointCoord(427, 2) =1.172894386278138e-03;
pointCoord(428, 0) = 2.952384803829710e-01; pointCoord(428, 1) =6.454530708565701e-01; pointCoord(428, 2) =1.172894386278138e-03;
pointCoord(429, 0) = 4.523848038297097e-02; pointCoord(429, 1) =8.954530708565701e-01; pointCoord(429, 2) =1.172894386278138e-03;
pointCoord(430, 0) = 7.785605735548751e-02; pointCoord(430, 1) =1.128354938840536e-01; pointCoord(430, 2) =1.356011353626800e-03;
pointCoord(431, 0) = 3.278560573554875e-01; pointCoord(431, 1) =1.128354938840536e-01; pointCoord(431, 2) =1.356011353626800e-03;
pointCoord(432, 0) = 5.778560573554875e-01; pointCoord(432, 1) =1.128354938840536e-01; pointCoord(432, 2) =1.356011353626800e-03;
pointCoord(433, 0) = 8.278560573554875e-01; pointCoord(433, 1) =1.128354938840536e-01; pointCoord(433, 2) =1.356011353626800e-03;
pointCoord(434, 0) = 7.785605735548751e-02; pointCoord(434, 1) =3.628354938840536e-01; pointCoord(434, 2) =1.356011353626800e-03;
pointCoord(435, 0) = 3.278560573554875e-01; pointCoord(435, 1) =3.628354938840536e-01; pointCoord(435, 2) =1.356011353626800e-03;
pointCoord(436, 0) = 5.778560573554875e-01; pointCoord(436, 1) =3.628354938840536e-01; pointCoord(436, 2) =1.356011353626800e-03;
pointCoord(437, 0) = 7.785605735548751e-02; pointCoord(437, 1) =6.128354938840536e-01; pointCoord(437, 2) =1.356011353626800e-03;
pointCoord(438, 0) = 3.278560573554875e-01; pointCoord(438, 1) =6.128354938840536e-01; pointCoord(438, 2) =1.356011353626800e-03;
pointCoord(439, 0) = 7.785605735548751e-02; pointCoord(439, 1) =8.628354938840536e-01; pointCoord(439, 2) =1.356011353626800e-03;
pointCoord(440, 0) = 1.128354938840536e-01; pointCoord(440, 1) =7.785605735548751e-02; pointCoord(440, 2) =1.356011353626800e-03;
pointCoord(441, 0) = 3.628354938840536e-01; pointCoord(441, 1) =7.785605735548751e-02; pointCoord(441, 2) =1.356011353626800e-03;
pointCoord(442, 0) = 6.128354938840536e-01; pointCoord(442, 1) =7.785605735548751e-02; pointCoord(442, 2) =1.356011353626800e-03;
pointCoord(443, 0) = 8.628354938840536e-01; pointCoord(443, 1) =7.785605735548751e-02; pointCoord(443, 2) =1.356011353626800e-03;
pointCoord(444, 0) = 1.128354938840536e-01; pointCoord(444, 1) =3.278560573554875e-01; pointCoord(444, 2) =1.356011353626800e-03;
pointCoord(445, 0) = 3.628354938840536e-01; pointCoord(445, 1) =3.278560573554875e-01; pointCoord(445, 2) =1.356011353626800e-03;
pointCoord(446, 0) = 6.128354938840536e-01; pointCoord(446, 1) =3.278560573554875e-01; pointCoord(446, 2) =1.356011353626800e-03;
pointCoord(447, 0) = 1.128354938840536e-01; pointCoord(447, 1) =5.778560573554875e-01; pointCoord(447, 2) =1.356011353626800e-03;
pointCoord(448, 0) = 3.628354938840536e-01; pointCoord(448, 1) =5.778560573554875e-01; pointCoord(448, 2) =1.356011353626800e-03;
pointCoord(449, 0) = 1.128354938840536e-01; pointCoord(449, 1) =8.278560573554875e-01; pointCoord(449, 2) =1.356011353626800e-03;
pointCoord(450, 0) = 1.454530708565701e-01; pointCoord(450, 1) =4.523848038297097e-02; pointCoord(450, 2) =1.172894386278138e-03;
pointCoord(451, 0) = 3.954530708565701e-01; pointCoord(451, 1) =4.523848038297097e-02; pointCoord(451, 2) =1.172894386278138e-03;
pointCoord(452, 0) = 6.454530708565701e-01; pointCoord(452, 1) =4.523848038297097e-02; pointCoord(452, 2) =1.172894386278138e-03;
pointCoord(453, 0) = 8.954530708565701e-01; pointCoord(453, 1) =4.523848038297097e-02; pointCoord(453, 2) =1.172894386278138e-03;
pointCoord(454, 0) = 1.454530708565701e-01; pointCoord(454, 1) =2.952384803829710e-01; pointCoord(454, 2) =1.172894386278138e-03;
pointCoord(455, 0) = 3.954530708565701e-01; pointCoord(455, 1) =2.952384803829710e-01; pointCoord(455, 2) =1.172894386278138e-03;
pointCoord(456, 0) = 6.454530708565701e-01; pointCoord(456, 1) =2.952384803829710e-01; pointCoord(456, 2) =1.172894386278138e-03;
pointCoord(457, 0) = 1.454530708565701e-01; pointCoord(457, 1) =5.452384803829710e-01; pointCoord(457, 2) =1.172894386278138e-03;
pointCoord(458, 0) = 3.954530708565701e-01; pointCoord(458, 1) =5.452384803829710e-01; pointCoord(458, 2) =1.172894386278138e-03;
pointCoord(459, 0) = 1.454530708565701e-01; pointCoord(459, 1) =7.952384803829710e-01; pointCoord(459, 2) =1.172894386278138e-03;
pointCoord(460, 0) = 1.713045588190432e-01; pointCoord(460, 1) =1.938699242049790e-02; pointCoord(460, 2) =8.314438675507627e-04;
pointCoord(461, 0) = 4.213045588190432e-01; pointCoord(461, 1) =1.938699242049790e-02; pointCoord(461, 2) =8.314438675507627e-04;
pointCoord(462, 0) = 6.713045588190432e-01; pointCoord(462, 1) =1.938699242049790e-02; pointCoord(462, 2) =8.314438675507627e-04;
pointCoord(463, 0) = 9.213045588190432e-01; pointCoord(463, 1) =1.938699242049790e-02; pointCoord(463, 2) =8.314438675507627e-04;
pointCoord(464, 0) = 1.713045588190432e-01; pointCoord(464, 1) =2.693869924204979e-01; pointCoord(464, 2) =8.314438675507627e-04;
pointCoord(465, 0) = 4.213045588190432e-01; pointCoord(465, 1) =2.693869924204979e-01; pointCoord(465, 2) =8.314438675507627e-04;
pointCoord(466, 0) = 6.713045588190432e-01; pointCoord(466, 1) =2.693869924204979e-01; pointCoord(466, 2) =8.314438675507627e-04;
pointCoord(467, 0) = 1.713045588190432e-01; pointCoord(467, 1) =5.193869924204979e-01; pointCoord(467, 2) =8.314438675507627e-04;
pointCoord(468, 0) = 4.213045588190432e-01; pointCoord(468, 1) =5.193869924204979e-01; pointCoord(468, 2) =8.314438675507627e-04;
pointCoord(469, 0) = 1.713045588190432e-01; pointCoord(469, 1) =7.693869924204979e-01; pointCoord(469, 2) =8.314438675507627e-04;
pointCoord(470, 0) = 1.869053568073263e-01; pointCoord(470, 1) =3.786194432214795e-03; pointCoord(470, 2) =3.784758260822833e-04;
pointCoord(471, 0) = 4.369053568073263e-01; pointCoord(471, 1) =3.786194432214795e-03; pointCoord(471, 2) =3.784758260822833e-04;
pointCoord(472, 0) = 6.869053568073263e-01; pointCoord(472, 1) =3.786194432214795e-03; pointCoord(472, 2) =3.784758260822833e-04;
pointCoord(473, 0) = 9.369053568073263e-01; pointCoord(473, 1) =3.786194432214795e-03; pointCoord(473, 2) =3.784758260822833e-04;
pointCoord(474, 0) = 1.869053568073263e-01; pointCoord(474, 1) =2.537861944322148e-01; pointCoord(474, 2) =3.784758260822833e-04;
pointCoord(475, 0) = 4.369053568073263e-01; pointCoord(475, 1) =2.537861944322148e-01; pointCoord(475, 2) =3.784758260822833e-04;
pointCoord(476, 0) = 6.869053568073263e-01; pointCoord(476, 1) =2.537861944322148e-01; pointCoord(476, 2) =3.784758260822833e-04;
pointCoord(477, 0) = 1.869053568073263e-01; pointCoord(477, 1) =5.037861944322148e-01; pointCoord(477, 2) =3.784758260822833e-04;
pointCoord(478, 0) = 4.369053568073263e-01; pointCoord(478, 1) =5.037861944322148e-01; pointCoord(478, 2) =3.784758260822833e-04;
pointCoord(479, 0) = 1.869053568073263e-01; pointCoord(479, 1) =7.537861944322148e-01; pointCoord(479, 2) =3.784758260822833e-04;
pointCoord(480, 0) = 4.459117727760069e-03; pointCoord(480, 1) =2.201241919489433e-01; pointCoord(480, 2) =3.159789840669076e-04;
pointCoord(481, 0) = 2.544591177277601e-01; pointCoord(481, 1) =2.201241919489433e-01; pointCoord(481, 2) =3.159789840669076e-04;
pointCoord(482, 0) = 5.044591177277601e-01; pointCoord(482, 1) =2.201241919489433e-01; pointCoord(482, 2) =3.159789840669076e-04;
pointCoord(483, 0) = 7.544591177277601e-01; pointCoord(483, 1) =2.201241919489433e-01; pointCoord(483, 2) =3.159789840669076e-04;
pointCoord(484, 0) = 4.459117727760069e-03; pointCoord(484, 1) =4.701241919489433e-01; pointCoord(484, 2) =3.159789840669076e-04;
pointCoord(485, 0) = 2.544591177277601e-01; pointCoord(485, 1) =4.701241919489433e-01; pointCoord(485, 2) =3.159789840669076e-04;
pointCoord(486, 0) = 5.044591177277601e-01; pointCoord(486, 1) =4.701241919489433e-01; pointCoord(486, 2) =3.159789840669076e-04;
pointCoord(487, 0) = 4.459117727760069e-03; pointCoord(487, 1) =7.201241919489433e-01; pointCoord(487, 2) =3.159789840669076e-04;
pointCoord(488, 0) = 2.544591177277601e-01; pointCoord(488, 1) =7.201241919489433e-01; pointCoord(488, 2) =3.159789840669076e-04;
pointCoord(489, 0) = 4.459117727760069e-03; pointCoord(489, 1) =9.701241919489433e-01; pointCoord(489, 2) =3.159789840669076e-04;
pointCoord(490, 0) = 2.283265773533521e-02; pointCoord(490, 1) =2.017506519413681e-01; pointCoord(490, 2) =6.941494554535528e-04;
pointCoord(491, 0) = 2.728326577353352e-01; pointCoord(491, 1) =2.017506519413681e-01; pointCoord(491, 2) =6.941494554535528e-04;
pointCoord(492, 0) = 5.228326577353352e-01; pointCoord(492, 1) =2.017506519413681e-01; pointCoord(492, 2) =6.941494554535528e-04;
pointCoord(493, 0) = 7.728326577353352e-01; pointCoord(493, 1) =2.017506519413681e-01; pointCoord(493, 2) =6.941494554535528e-04;
pointCoord(494, 0) = 2.283265773533521e-02; pointCoord(494, 1) =4.517506519413681e-01; pointCoord(494, 2) =6.941494554535528e-04;
pointCoord(495, 0) = 2.728326577353352e-01; pointCoord(495, 1) =4.517506519413681e-01; pointCoord(495, 2) =6.941494554535528e-04;
pointCoord(496, 0) = 5.228326577353352e-01; pointCoord(496, 1) =4.517506519413681e-01; pointCoord(496, 2) =6.941494554535528e-04;
pointCoord(497, 0) = 2.283265773533521e-02; pointCoord(497, 1) =7.017506519413681e-01; pointCoord(497, 2) =6.941494554535528e-04;
pointCoord(498, 0) = 2.728326577353352e-01; pointCoord(498, 1) =7.017506519413681e-01; pointCoord(498, 2) =6.941494554535528e-04;
pointCoord(499, 0) = 2.283265773533521e-02; pointCoord(499, 1) =9.517506519413681e-01; pointCoord(499, 2) =6.941494554535528e-04;
pointCoord(500, 0) = 5.327875085766011e-02; pointCoord(500, 1) =1.713045588190432e-01; pointCoord(500, 2) =9.792170359471568e-04;
pointCoord(501, 0) = 3.032787508576601e-01; pointCoord(501, 1) =1.713045588190432e-01; pointCoord(501, 2) =9.792170359471568e-04;
pointCoord(502, 0) = 5.532787508576601e-01; pointCoord(502, 1) =1.713045588190432e-01; pointCoord(502, 2) =9.792170359471568e-04;
pointCoord(503, 0) = 8.032787508576601e-01; pointCoord(503, 1) =1.713045588190432e-01; pointCoord(503, 2) =9.792170359471568e-04;
pointCoord(504, 0) = 5.327875085766011e-02; pointCoord(504, 1) =4.213045588190432e-01; pointCoord(504, 2) =9.792170359471568e-04;
pointCoord(505, 0) = 3.032787508576601e-01; pointCoord(505, 1) =4.213045588190432e-01; pointCoord(505, 2) =9.792170359471568e-04;
pointCoord(506, 0) = 5.532787508576601e-01; pointCoord(506, 1) =4.213045588190432e-01; pointCoord(506, 2) =9.792170359471568e-04;
pointCoord(507, 0) = 5.327875085766011e-02; pointCoord(507, 1) =6.713045588190432e-01; pointCoord(507, 2) =9.792170359471568e-04;
pointCoord(508, 0) = 3.032787508576601e-01; pointCoord(508, 1) =6.713045588190432e-01; pointCoord(508, 2) =9.792170359471568e-04;
pointCoord(509, 0) = 5.327875085766011e-02; pointCoord(509, 1) =9.213045588190432e-01; pointCoord(509, 2) =9.792170359471568e-04;
pointCoord(510, 0) = 9.169347527783372e-02; pointCoord(510, 1) =1.328898343988696e-01; pointCoord(510, 2) =1.132096319961624e-03;
pointCoord(511, 0) = 3.416934752778337e-01; pointCoord(511, 1) =1.328898343988696e-01; pointCoord(511, 2) =1.132096319961624e-03;
pointCoord(512, 0) = 5.916934752778338e-01; pointCoord(512, 1) =1.328898343988696e-01; pointCoord(512, 2) =1.132096319961624e-03;
pointCoord(513, 0) = 8.416934752778338e-01; pointCoord(513, 1) =1.328898343988696e-01; pointCoord(513, 2) =1.132096319961624e-03;
pointCoord(514, 0) = 9.169347527783372e-02; pointCoord(514, 1) =3.828898343988696e-01; pointCoord(514, 2) =1.132096319961624e-03;
pointCoord(515, 0) = 3.416934752778337e-01; pointCoord(515, 1) =3.828898343988696e-01; pointCoord(515, 2) =1.132096319961624e-03;
pointCoord(516, 0) = 5.916934752778338e-01; pointCoord(516, 1) =3.828898343988696e-01; pointCoord(516, 2) =1.132096319961624e-03;
pointCoord(517, 0) = 9.169347527783372e-02; pointCoord(517, 1) =6.328898343988696e-01; pointCoord(517, 2) =1.132096319961624e-03;
pointCoord(518, 0) = 3.416934752778337e-01; pointCoord(518, 1) =6.328898343988696e-01; pointCoord(518, 2) =1.132096319961624e-03;
pointCoord(519, 0) = 9.169347527783372e-02; pointCoord(519, 1) =8.828898343988696e-01; pointCoord(519, 2) =1.132096319961624e-03;
pointCoord(520, 0) = 1.328898343988696e-01; pointCoord(520, 1) =9.169347527783372e-02; pointCoord(520, 2) =1.132096319961624e-03;
pointCoord(521, 0) = 3.828898343988696e-01; pointCoord(521, 1) =9.169347527783372e-02; pointCoord(521, 2) =1.132096319961624e-03;
pointCoord(522, 0) = 6.328898343988696e-01; pointCoord(522, 1) =9.169347527783372e-02; pointCoord(522, 2) =1.132096319961624e-03;
pointCoord(523, 0) = 8.828898343988696e-01; pointCoord(523, 1) =9.169347527783372e-02; pointCoord(523, 2) =1.132096319961624e-03;
pointCoord(524, 0) = 1.328898343988696e-01; pointCoord(524, 1) =3.416934752778337e-01; pointCoord(524, 2) =1.132096319961624e-03;
pointCoord(525, 0) = 3.828898343988696e-01; pointCoord(525, 1) =3.416934752778337e-01; pointCoord(525, 2) =1.132096319961624e-03;
pointCoord(526, 0) = 6.328898343988696e-01; pointCoord(526, 1) =3.416934752778337e-01; pointCoord(526, 2) =1.132096319961624e-03;
pointCoord(527, 0) = 1.328898343988696e-01; pointCoord(527, 1) =5.916934752778338e-01; pointCoord(527, 2) =1.132096319961624e-03;
pointCoord(528, 0) = 3.828898343988696e-01; pointCoord(528, 1) =5.916934752778338e-01; pointCoord(528, 2) =1.132096319961624e-03;
pointCoord(529, 0) = 1.328898343988696e-01; pointCoord(529, 1) =8.416934752778338e-01; pointCoord(529, 2) =1.132096319961624e-03;
pointCoord(530, 0) = 1.713045588190432e-01; pointCoord(530, 1) =5.327875085766011e-02; pointCoord(530, 2) =9.792170359471568e-04;
pointCoord(531, 0) = 4.213045588190432e-01; pointCoord(531, 1) =5.327875085766011e-02; pointCoord(531, 2) =9.792170359471568e-04;
pointCoord(532, 0) = 6.713045588190432e-01; pointCoord(532, 1) =5.327875085766011e-02; pointCoord(532, 2) =9.792170359471568e-04;
pointCoord(533, 0) = 9.213045588190432e-01; pointCoord(533, 1) =5.327875085766011e-02; pointCoord(533, 2) =9.792170359471568e-04;
pointCoord(534, 0) = 1.713045588190432e-01; pointCoord(534, 1) =3.032787508576601e-01; pointCoord(534, 2) =9.792170359471568e-04;
pointCoord(535, 0) = 4.213045588190432e-01; pointCoord(535, 1) =3.032787508576601e-01; pointCoord(535, 2) =9.792170359471568e-04;
pointCoord(536, 0) = 6.713045588190432e-01; pointCoord(536, 1) =3.032787508576601e-01; pointCoord(536, 2) =9.792170359471568e-04;
pointCoord(537, 0) = 1.713045588190432e-01; pointCoord(537, 1) =5.532787508576601e-01; pointCoord(537, 2) =9.792170359471568e-04;
pointCoord(538, 0) = 4.213045588190432e-01; pointCoord(538, 1) =5.532787508576601e-01; pointCoord(538, 2) =9.792170359471568e-04;
pointCoord(539, 0) = 1.713045588190432e-01; pointCoord(539, 1) =8.032787508576601e-01; pointCoord(539, 2) =9.792170359471568e-04;
pointCoord(540, 0) = 2.017506519413681e-01; pointCoord(540, 1) =2.283265773533521e-02; pointCoord(540, 2) =6.941494554535528e-04;
pointCoord(541, 0) = 4.517506519413681e-01; pointCoord(541, 1) =2.283265773533521e-02; pointCoord(541, 2) =6.941494554535528e-04;
pointCoord(542, 0) = 7.017506519413681e-01; pointCoord(542, 1) =2.283265773533521e-02; pointCoord(542, 2) =6.941494554535528e-04;
pointCoord(543, 0) = 9.517506519413681e-01; pointCoord(543, 1) =2.283265773533521e-02; pointCoord(543, 2) =6.941494554535528e-04;
pointCoord(544, 0) = 2.017506519413681e-01; pointCoord(544, 1) =2.728326577353352e-01; pointCoord(544, 2) =6.941494554535528e-04;
pointCoord(545, 0) = 4.517506519413681e-01; pointCoord(545, 1) =2.728326577353352e-01; pointCoord(545, 2) =6.941494554535528e-04;
pointCoord(546, 0) = 7.017506519413681e-01; pointCoord(546, 1) =2.728326577353352e-01; pointCoord(546, 2) =6.941494554535528e-04;
pointCoord(547, 0) = 2.017506519413681e-01; pointCoord(547, 1) =5.228326577353352e-01; pointCoord(547, 2) =6.941494554535528e-04;
pointCoord(548, 0) = 4.517506519413681e-01; pointCoord(548, 1) =5.228326577353352e-01; pointCoord(548, 2) =6.941494554535528e-04;
pointCoord(549, 0) = 2.017506519413681e-01; pointCoord(549, 1) =7.728326577353352e-01; pointCoord(549, 2) =6.941494554535528e-04;
pointCoord(550, 0) = 2.201241919489433e-01; pointCoord(550, 1) =4.459117727760069e-03; pointCoord(550, 2) =3.159789840669076e-04;
pointCoord(551, 0) = 4.701241919489433e-01; pointCoord(551, 1) =4.459117727760069e-03; pointCoord(551, 2) =3.159789840669076e-04;
pointCoord(552, 0) = 7.201241919489433e-01; pointCoord(552, 1) =4.459117727760069e-03; pointCoord(552, 2) =3.159789840669076e-04;
pointCoord(553, 0) = 9.701241919489433e-01; pointCoord(553, 1) =4.459117727760069e-03; pointCoord(553, 2) =3.159789840669076e-04;
pointCoord(554, 0) = 2.201241919489433e-01; pointCoord(554, 1) =2.544591177277601e-01; pointCoord(554, 2) =3.159789840669076e-04;
pointCoord(555, 0) = 4.701241919489433e-01; pointCoord(555, 1) =2.544591177277601e-01; pointCoord(555, 2) =3.159789840669076e-04;
pointCoord(556, 0) = 7.201241919489433e-01; pointCoord(556, 1) =2.544591177277601e-01; pointCoord(556, 2) =3.159789840669076e-04;
pointCoord(557, 0) = 2.201241919489433e-01; pointCoord(557, 1) =5.044591177277601e-01; pointCoord(557, 2) =3.159789840669076e-04;
pointCoord(558, 0) = 4.701241919489433e-01; pointCoord(558, 1) =5.044591177277601e-01; pointCoord(558, 2) =3.159789840669076e-04;
pointCoord(559, 0) = 2.201241919489433e-01; pointCoord(559, 1) =7.544591177277601e-01; pointCoord(559, 2) =3.159789840669076e-04;
pointCoord(560, 0) = 4.865211969246323e-03; pointCoord(560, 1) =2.401710200929457e-01; pointCoord(560, 2) =1.569337084289531e-04;
pointCoord(561, 0) = 2.548652119692463e-01; pointCoord(561, 1) =2.401710200929457e-01; pointCoord(561, 2) =1.569337084289531e-04;
pointCoord(562, 0) = 5.048652119692463e-01; pointCoord(562, 1) =2.401710200929457e-01; pointCoord(562, 2) =1.569337084289531e-04;
pointCoord(563, 0) = 7.548652119692463e-01; pointCoord(563, 1) =2.401710200929457e-01; pointCoord(563, 2) =1.569337084289531e-04;
pointCoord(564, 0) = 4.865211969246323e-03; pointCoord(564, 1) =4.901710200929458e-01; pointCoord(564, 2) =1.569337084289531e-04;
pointCoord(565, 0) = 2.548652119692463e-01; pointCoord(565, 1) =4.901710200929458e-01; pointCoord(565, 2) =1.569337084289531e-04;
pointCoord(566, 0) = 5.048652119692463e-01; pointCoord(566, 1) =4.901710200929458e-01; pointCoord(566, 2) =1.569337084289531e-04;
pointCoord(567, 0) = 4.865211969246323e-03; pointCoord(567, 1) =7.401710200929458e-01; pointCoord(567, 2) =1.569337084289531e-04;
pointCoord(568, 0) = 2.548652119692463e-01; pointCoord(568, 1) =7.401710200929458e-01; pointCoord(568, 2) =1.569337084289531e-04;
pointCoord(569, 0) = 4.865211969246323e-03; pointCoord(569, 1) =9.901710200929458e-01; pointCoord(569, 2) =1.569337084289531e-04;
pointCoord(570, 0) = 2.491204011324876e-02; pointCoord(570, 1) =2.201241919489433e-01; pointCoord(570, 2) =3.447553595058008e-04;
pointCoord(571, 0) = 2.749120401132488e-01; pointCoord(571, 1) =2.201241919489433e-01; pointCoord(571, 2) =3.447553595058008e-04;
pointCoord(572, 0) = 5.249120401132488e-01; pointCoord(572, 1) =2.201241919489433e-01; pointCoord(572, 2) =3.447553595058008e-04;
pointCoord(573, 0) = 7.749120401132488e-01; pointCoord(573, 1) =2.201241919489433e-01; pointCoord(573, 2) =3.447553595058008e-04;
pointCoord(574, 0) = 2.491204011324876e-02; pointCoord(574, 1) =4.701241919489433e-01; pointCoord(574, 2) =3.447553595058008e-04;
pointCoord(575, 0) = 2.749120401132488e-01; pointCoord(575, 1) =4.701241919489433e-01; pointCoord(575, 2) =3.447553595058008e-04;
pointCoord(576, 0) = 5.249120401132488e-01; pointCoord(576, 1) =4.701241919489433e-01; pointCoord(576, 2) =3.447553595058008e-04;
pointCoord(577, 0) = 2.491204011324876e-02; pointCoord(577, 1) =7.201241919489433e-01; pointCoord(577, 2) =3.447553595058008e-04;
pointCoord(578, 0) = 2.749120401132488e-01; pointCoord(578, 1) =7.201241919489433e-01; pointCoord(578, 2) =3.447553595058008e-04;
pointCoord(579, 0) = 2.491204011324876e-02; pointCoord(579, 1) =9.701241919489433e-01; pointCoord(579, 2) =3.447553595058008e-04;
pointCoord(580, 0) = 5.813087525486571e-02; pointCoord(580, 1) =1.869053568073263e-01; pointCoord(580, 2) =4.863366507167933e-04;
pointCoord(581, 0) = 3.081308752548657e-01; pointCoord(581, 1) =1.869053568073263e-01; pointCoord(581, 2) =4.863366507167933e-04;
pointCoord(582, 0) = 5.581308752548657e-01; pointCoord(582, 1) =1.869053568073263e-01; pointCoord(582, 2) =4.863366507167933e-04;
pointCoord(583, 0) = 8.081308752548657e-01; pointCoord(583, 1) =1.869053568073263e-01; pointCoord(583, 2) =4.863366507167933e-04;
pointCoord(584, 0) = 5.813087525486571e-02; pointCoord(584, 1) =4.369053568073263e-01; pointCoord(584, 2) =4.863366507167933e-04;
pointCoord(585, 0) = 3.081308752548657e-01; pointCoord(585, 1) =4.369053568073263e-01; pointCoord(585, 2) =4.863366507167933e-04;
pointCoord(586, 0) = 5.581308752548657e-01; pointCoord(586, 1) =4.369053568073263e-01; pointCoord(586, 2) =4.863366507167933e-04;
pointCoord(587, 0) = 5.813087525486571e-02; pointCoord(587, 1) =6.869053568073263e-01; pointCoord(587, 2) =4.863366507167933e-04;
pointCoord(588, 0) = 3.081308752548657e-01; pointCoord(588, 1) =6.869053568073263e-01; pointCoord(588, 2) =4.863366507167933e-04;
pointCoord(589, 0) = 5.813087525486571e-02; pointCoord(589, 1) =9.369053568073263e-01; pointCoord(589, 2) =4.863366507167933e-04;
pointCoord(590, 0) = 1.000440492176914e-01; pointCoord(590, 1) =1.449921828445007e-01; pointCoord(590, 2) =5.622654757087529e-04;
pointCoord(591, 0) = 3.500440492176914e-01; pointCoord(591, 1) =1.449921828445007e-01; pointCoord(591, 2) =5.622654757087529e-04;
pointCoord(592, 0) = 6.000440492176914e-01; pointCoord(592, 1) =1.449921828445007e-01; pointCoord(592, 2) =5.622654757087529e-04;
pointCoord(593, 0) = 8.500440492176914e-01; pointCoord(593, 1) =1.449921828445007e-01; pointCoord(593, 2) =5.622654757087529e-04;
pointCoord(594, 0) = 1.000440492176914e-01; pointCoord(594, 1) =3.949921828445007e-01; pointCoord(594, 2) =5.622654757087529e-04;
pointCoord(595, 0) = 3.500440492176914e-01; pointCoord(595, 1) =3.949921828445007e-01; pointCoord(595, 2) =5.622654757087529e-04;
pointCoord(596, 0) = 6.000440492176914e-01; pointCoord(596, 1) =3.949921828445007e-01; pointCoord(596, 2) =5.622654757087529e-04;
pointCoord(597, 0) = 1.000440492176914e-01; pointCoord(597, 1) =6.449921828445007e-01; pointCoord(597, 2) =5.622654757087529e-04;
pointCoord(598, 0) = 3.500440492176914e-01; pointCoord(598, 1) =6.449921828445007e-01; pointCoord(598, 2) =5.622654757087529e-04;
pointCoord(599, 0) = 1.000440492176914e-01; pointCoord(599, 1) =8.949921828445007e-01; pointCoord(599, 2) =5.622654757087529e-04;
pointCoord(600, 0) = 1.449921828445007e-01; pointCoord(600, 1) =1.000440492176914e-01; pointCoord(600, 2) =5.622654757087529e-04;
pointCoord(601, 0) = 3.949921828445007e-01; pointCoord(601, 1) =1.000440492176914e-01; pointCoord(601, 2) =5.622654757087529e-04;
pointCoord(602, 0) = 6.449921828445007e-01; pointCoord(602, 1) =1.000440492176914e-01; pointCoord(602, 2) =5.622654757087529e-04;
pointCoord(603, 0) = 8.949921828445007e-01; pointCoord(603, 1) =1.000440492176914e-01; pointCoord(603, 2) =5.622654757087529e-04;
pointCoord(604, 0) = 1.449921828445007e-01; pointCoord(604, 1) =3.500440492176914e-01; pointCoord(604, 2) =5.622654757087529e-04;
pointCoord(605, 0) = 3.949921828445007e-01; pointCoord(605, 1) =3.500440492176914e-01; pointCoord(605, 2) =5.622654757087529e-04;
pointCoord(606, 0) = 6.449921828445007e-01; pointCoord(606, 1) =3.500440492176914e-01; pointCoord(606, 2) =5.622654757087529e-04;
pointCoord(607, 0) = 1.449921828445007e-01; pointCoord(607, 1) =6.000440492176914e-01; pointCoord(607, 2) =5.622654757087529e-04;
pointCoord(608, 0) = 3.949921828445007e-01; pointCoord(608, 1) =6.000440492176914e-01; pointCoord(608, 2) =5.622654757087529e-04;
pointCoord(609, 0) = 1.449921828445007e-01; pointCoord(609, 1) =8.500440492176914e-01; pointCoord(609, 2) =5.622654757087529e-04;
pointCoord(610, 0) = 1.869053568073263e-01; pointCoord(610, 1) =5.813087525486571e-02; pointCoord(610, 2) =4.863366507167933e-04;
pointCoord(611, 0) = 4.369053568073263e-01; pointCoord(611, 1) =5.813087525486571e-02; pointCoord(611, 2) =4.863366507167933e-04;
pointCoord(612, 0) = 6.869053568073263e-01; pointCoord(612, 1) =5.813087525486571e-02; pointCoord(612, 2) =4.863366507167933e-04;
pointCoord(613, 0) = 9.369053568073263e-01; pointCoord(613, 1) =5.813087525486571e-02; pointCoord(613, 2) =4.863366507167933e-04;
pointCoord(614, 0) = 1.869053568073263e-01; pointCoord(614, 1) =3.081308752548657e-01; pointCoord(614, 2) =4.863366507167933e-04;
pointCoord(615, 0) = 4.369053568073263e-01; pointCoord(615, 1) =3.081308752548657e-01; pointCoord(615, 2) =4.863366507167933e-04;
pointCoord(616, 0) = 6.869053568073263e-01; pointCoord(616, 1) =3.081308752548657e-01; pointCoord(616, 2) =4.863366507167933e-04;
pointCoord(617, 0) = 1.869053568073263e-01; pointCoord(617, 1) =5.581308752548657e-01; pointCoord(617, 2) =4.863366507167933e-04;
pointCoord(618, 0) = 4.369053568073263e-01; pointCoord(618, 1) =5.581308752548657e-01; pointCoord(618, 2) =4.863366507167933e-04;
pointCoord(619, 0) = 1.869053568073263e-01; pointCoord(619, 1) =8.081308752548657e-01; pointCoord(619, 2) =4.863366507167933e-04;
pointCoord(620, 0) = 2.201241919489433e-01; pointCoord(620, 1) =2.491204011324876e-02; pointCoord(620, 2) =3.447553595058008e-04;
pointCoord(621, 0) = 4.701241919489433e-01; pointCoord(621, 1) =2.491204011324876e-02; pointCoord(621, 2) =3.447553595058008e-04;
pointCoord(622, 0) = 7.201241919489433e-01; pointCoord(622, 1) =2.491204011324876e-02; pointCoord(622, 2) =3.447553595058008e-04;
pointCoord(623, 0) = 9.701241919489433e-01; pointCoord(623, 1) =2.491204011324876e-02; pointCoord(623, 2) =3.447553595058008e-04;
pointCoord(624, 0) = 2.201241919489433e-01; pointCoord(624, 1) =2.749120401132488e-01; pointCoord(624, 2) =3.447553595058008e-04;
pointCoord(625, 0) = 4.701241919489433e-01; pointCoord(625, 1) =2.749120401132488e-01; pointCoord(625, 2) =3.447553595058008e-04;
pointCoord(626, 0) = 7.201241919489433e-01; pointCoord(626, 1) =2.749120401132488e-01; pointCoord(626, 2) =3.447553595058008e-04;
pointCoord(627, 0) = 2.201241919489433e-01; pointCoord(627, 1) =5.249120401132488e-01; pointCoord(627, 2) =3.447553595058008e-04;
pointCoord(628, 0) = 4.701241919489433e-01; pointCoord(628, 1) =5.249120401132488e-01; pointCoord(628, 2) =3.447553595058008e-04;
pointCoord(629, 0) = 2.201241919489433e-01; pointCoord(629, 1) =7.749120401132488e-01; pointCoord(629, 2) =3.447553595058008e-04;
pointCoord(630, 0) = 2.401710200929457e-01; pointCoord(630, 1) =4.865211969246323e-03; pointCoord(630, 2) =1.569337084289531e-04;
pointCoord(631, 0) = 4.901710200929458e-01; pointCoord(631, 1) =4.865211969246323e-03; pointCoord(631, 2) =1.569337084289531e-04;
pointCoord(632, 0) = 7.401710200929458e-01; pointCoord(632, 1) =4.865211969246323e-03; pointCoord(632, 2) =1.569337084289531e-04;
pointCoord(633, 0) = 9.901710200929458e-01; pointCoord(633, 1) =4.865211969246323e-03; pointCoord(633, 2) =1.569337084289531e-04;
pointCoord(634, 0) = 2.401710200929457e-01; pointCoord(634, 1) =2.548652119692463e-01; pointCoord(634, 2) =1.569337084289531e-04;
pointCoord(635, 0) = 4.901710200929458e-01; pointCoord(635, 1) =2.548652119692463e-01; pointCoord(635, 2) =1.569337084289531e-04;
pointCoord(636, 0) = 7.401710200929458e-01; pointCoord(636, 1) =2.548652119692463e-01; pointCoord(636, 2) =1.569337084289531e-04;
pointCoord(637, 0) = 2.401710200929457e-01; pointCoord(637, 1) =5.048652119692463e-01; pointCoord(637, 2) =1.569337084289531e-04;
pointCoord(638, 0) = 4.901710200929458e-01; pointCoord(638, 1) =5.048652119692463e-01; pointCoord(638, 2) =1.569337084289531e-04;
pointCoord(639, 0) = 2.401710200929457e-01; pointCoord(639, 1) =7.548652119692463e-01; pointCoord(639, 2) =1.569337084289531e-04;
pointCoord(640, 0) = 2.499014440314384e-01; pointCoord(640, 1) =2.451347880307537e-01; pointCoord(640, 2) =3.179050313111366e-06;
pointCoord(641, 0) = 4.999014440314384e-01; pointCoord(641, 1) =2.451347880307537e-01; pointCoord(641, 2) =3.179050313111366e-06;
pointCoord(642, 0) = 7.499014440314383e-01; pointCoord(642, 1) =2.451347880307537e-01; pointCoord(642, 2) =3.179050313111366e-06;
pointCoord(643, 0) = 2.499014440314384e-01; pointCoord(643, 1) =4.951347880307537e-01; pointCoord(643, 2) =3.179050313111366e-06;
pointCoord(644, 0) = 4.999014440314384e-01; pointCoord(644, 1) =4.951347880307537e-01; pointCoord(644, 2) =3.179050313111366e-06;
pointCoord(645, 0) = 2.499014440314384e-01; pointCoord(645, 1) =7.451347880307537e-01; pointCoord(645, 2) =3.179050313111366e-06;
pointCoord(646, 0) = 2.494953497899521e-01; pointCoord(646, 1) =2.455408822722399e-01; pointCoord(646, 2) =6.983806376307709e-06;
pointCoord(647, 0) = 4.994953497899521e-01; pointCoord(647, 1) =2.455408822722399e-01; pointCoord(647, 2) =6.983806376307709e-06;
pointCoord(648, 0) = 7.494953497899521e-01; pointCoord(648, 1) =2.455408822722399e-01; pointCoord(648, 2) =6.983806376307709e-06;
pointCoord(649, 0) = 2.494953497899521e-01; pointCoord(649, 1) =4.955408822722399e-01; pointCoord(649, 2) =6.983806376307709e-06;
pointCoord(650, 0) = 4.994953497899521e-01; pointCoord(650, 1) =4.955408822722399e-01; pointCoord(650, 2) =6.983806376307709e-06;
pointCoord(651, 0) = 2.494953497899521e-01; pointCoord(651, 1) =7.455408822722399e-01; pointCoord(651, 2) =6.983806376307709e-06;
pointCoord(652, 0) = 2.488224264944068e-01; pointCoord(652, 1) =2.462138055677852e-01; pointCoord(652, 2) =9.851858451676738e-06;
pointCoord(653, 0) = 4.988224264944068e-01; pointCoord(653, 1) =2.462138055677852e-01; pointCoord(653, 2) =9.851858451676738e-06;
pointCoord(654, 0) = 7.488224264944068e-01; pointCoord(654, 1) =2.462138055677852e-01; pointCoord(654, 2) =9.851858451676738e-06;
pointCoord(655, 0) = 2.488224264944068e-01; pointCoord(655, 1) =4.962138055677852e-01; pointCoord(655, 2) =9.851858451676738e-06;
pointCoord(656, 0) = 4.988224264944068e-01; pointCoord(656, 1) =4.962138055677852e-01; pointCoord(656, 2) =9.851858451676738e-06;
pointCoord(657, 0) = 2.488224264944068e-01; pointCoord(657, 1) =7.462138055677852e-01; pointCoord(657, 2) =9.851858451676738e-06;
pointCoord(658, 0) = 2.479733795296476e-01; pointCoord(658, 1) =2.470628525325444e-01; pointCoord(658, 2) =1.138997003574184e-05;
pointCoord(659, 0) = 4.979733795296476e-01; pointCoord(659, 1) =2.470628525325444e-01; pointCoord(659, 2) =1.138997003574184e-05;
pointCoord(660, 0) = 7.479733795296476e-01; pointCoord(660, 1) =2.470628525325444e-01; pointCoord(660, 2) =1.138997003574184e-05;
pointCoord(661, 0) = 2.479733795296476e-01; pointCoord(661, 1) =4.970628525325445e-01; pointCoord(661, 2) =1.138997003574184e-05;
pointCoord(662, 0) = 4.979733795296476e-01; pointCoord(662, 1) =4.970628525325445e-01; pointCoord(662, 2) =1.138997003574184e-05;
pointCoord(663, 0) = 2.479733795296476e-01; pointCoord(663, 1) =7.470628525325445e-01; pointCoord(663, 2) =1.138997003574184e-05;
pointCoord(664, 0) = 2.470628525325444e-01; pointCoord(664, 1) =2.479733795296476e-01; pointCoord(664, 2) =1.138997003574184e-05;
pointCoord(665, 0) = 4.970628525325445e-01; pointCoord(665, 1) =2.479733795296476e-01; pointCoord(665, 2) =1.138997003574184e-05;
pointCoord(666, 0) = 7.470628525325445e-01; pointCoord(666, 1) =2.479733795296476e-01; pointCoord(666, 2) =1.138997003574184e-05;
pointCoord(667, 0) = 2.470628525325444e-01; pointCoord(667, 1) =4.979733795296476e-01; pointCoord(667, 2) =1.138997003574184e-05;
pointCoord(668, 0) = 4.970628525325445e-01; pointCoord(668, 1) =4.979733795296476e-01; pointCoord(668, 2) =1.138997003574184e-05;
pointCoord(669, 0) = 2.470628525325444e-01; pointCoord(669, 1) =7.479733795296476e-01; pointCoord(669, 2) =1.138997003574184e-05;
pointCoord(670, 0) = 2.462138055677852e-01; pointCoord(670, 1) =2.488224264944068e-01; pointCoord(670, 2) =9.851858451676738e-06;
pointCoord(671, 0) = 4.962138055677852e-01; pointCoord(671, 1) =2.488224264944068e-01; pointCoord(671, 2) =9.851858451676738e-06;
pointCoord(672, 0) = 7.462138055677852e-01; pointCoord(672, 1) =2.488224264944068e-01; pointCoord(672, 2) =9.851858451676738e-06;
pointCoord(673, 0) = 2.462138055677852e-01; pointCoord(673, 1) =4.988224264944068e-01; pointCoord(673, 2) =9.851858451676738e-06;
pointCoord(674, 0) = 4.962138055677852e-01; pointCoord(674, 1) =4.988224264944068e-01; pointCoord(674, 2) =9.851858451676738e-06;
pointCoord(675, 0) = 2.462138055677852e-01; pointCoord(675, 1) =7.488224264944068e-01; pointCoord(675, 2) =9.851858451676738e-06;
pointCoord(676, 0) = 2.455408822722399e-01; pointCoord(676, 1) =2.494953497899521e-01; pointCoord(676, 2) =6.983806376307709e-06;
pointCoord(677, 0) = 4.955408822722399e-01; pointCoord(677, 1) =2.494953497899521e-01; pointCoord(677, 2) =6.983806376307709e-06;
pointCoord(678, 0) = 7.455408822722399e-01; pointCoord(678, 1) =2.494953497899521e-01; pointCoord(678, 2) =6.983806376307709e-06;
pointCoord(679, 0) = 2.455408822722399e-01; pointCoord(679, 1) =4.994953497899521e-01; pointCoord(679, 2) =6.983806376307709e-06;
pointCoord(680, 0) = 4.955408822722399e-01; pointCoord(680, 1) =4.994953497899521e-01; pointCoord(680, 2) =6.983806376307709e-06;
pointCoord(681, 0) = 2.455408822722399e-01; pointCoord(681, 1) =7.494953497899521e-01; pointCoord(681, 2) =6.983806376307709e-06;
pointCoord(682, 0) = 2.451347880307537e-01; pointCoord(682, 1) =2.499014440314384e-01; pointCoord(682, 2) =3.179050313111366e-06;
pointCoord(683, 0) = 4.951347880307537e-01; pointCoord(683, 1) =2.499014440314384e-01; pointCoord(683, 2) =3.179050313111366e-06;
pointCoord(684, 0) = 7.451347880307537e-01; pointCoord(684, 1) =2.499014440314384e-01; pointCoord(684, 2) =3.179050313111366e-06;
pointCoord(685, 0) = 2.451347880307537e-01; pointCoord(685, 1) =4.999014440314384e-01; pointCoord(685, 2) =3.179050313111366e-06;
pointCoord(686, 0) = 4.951347880307537e-01; pointCoord(686, 1) =4.999014440314384e-01; pointCoord(686, 2) =3.179050313111366e-06;
pointCoord(687, 0) = 2.451347880307537e-01; pointCoord(687, 1) =7.499014440314383e-01; pointCoord(687, 2) =3.179050313111366e-06;
pointCoord(688, 0) = 2.494953497899521e-01; pointCoord(688, 1) =2.250879598867512e-01; pointCoord(688, 2) =3.576018181520091e-05;
pointCoord(689, 0) = 4.994953497899521e-01; pointCoord(689, 1) =2.250879598867512e-01; pointCoord(689, 2) =3.576018181520091e-05;
pointCoord(690, 0) = 7.494953497899521e-01; pointCoord(690, 1) =2.250879598867512e-01; pointCoord(690, 2) =3.576018181520091e-05;
pointCoord(691, 0) = 2.494953497899521e-01; pointCoord(691, 1) =4.750879598867512e-01; pointCoord(691, 2) =3.576018181520091e-05;
pointCoord(692, 0) = 4.994953497899521e-01; pointCoord(692, 1) =4.750879598867512e-01; pointCoord(692, 2) =3.576018181520091e-05;
pointCoord(693, 0) = 2.494953497899521e-01; pointCoord(693, 1) =7.250879598867512e-01; pointCoord(693, 2) =3.576018181520091e-05;
pointCoord(694, 0) = 2.474159674120386e-01; pointCoord(694, 1) =2.271673422646648e-01; pointCoord(694, 2) =7.855873961758663e-05;
pointCoord(695, 0) = 4.974159674120385e-01; pointCoord(695, 1) =2.271673422646648e-01; pointCoord(695, 2) =7.855873961758663e-05;
pointCoord(696, 0) = 7.474159674120385e-01; pointCoord(696, 1) =2.271673422646648e-01; pointCoord(696, 2) =7.855873961758663e-05;
pointCoord(697, 0) = 2.474159674120386e-01; pointCoord(697, 1) =4.771673422646648e-01; pointCoord(697, 2) =7.855873961758663e-05;
pointCoord(698, 0) = 4.974159674120385e-01; pointCoord(698, 1) =4.771673422646648e-01; pointCoord(698, 2) =7.855873961758663e-05;
pointCoord(699, 0) = 2.474159674120386e-01; pointCoord(699, 1) =7.271673422646648e-01; pointCoord(699, 2) =7.855873961758663e-05;
pointCoord(700, 0) = 2.439703020972012e-01; pointCoord(700, 1) =2.306130075795021e-01; pointCoord(700, 2) =1.108205956969521e-04;
pointCoord(701, 0) = 4.939703020972012e-01; pointCoord(701, 1) =2.306130075795021e-01; pointCoord(701, 2) =1.108205956969521e-04;
pointCoord(702, 0) = 7.439703020972013e-01; pointCoord(702, 1) =2.306130075795021e-01; pointCoord(702, 2) =1.108205956969521e-04;
pointCoord(703, 0) = 2.439703020972012e-01; pointCoord(703, 1) =4.806130075795021e-01; pointCoord(703, 2) =1.108205956969521e-04;
pointCoord(704, 0) = 4.939703020972012e-01; pointCoord(704, 1) =4.806130075795021e-01; pointCoord(704, 2) =1.108205956969521e-04;
pointCoord(705, 0) = 2.439703020972012e-01; pointCoord(705, 1) =7.306130075795021e-01; pointCoord(705, 2) =1.108205956969521e-04;
pointCoord(706, 0) = 2.396228055897900e-01; pointCoord(706, 1) =2.349605040869134e-01; pointCoord(706, 2) =1.281223507749970e-04;
pointCoord(707, 0) = 4.896228055897899e-01; pointCoord(707, 1) =2.349605040869134e-01; pointCoord(707, 2) =1.281223507749970e-04;
pointCoord(708, 0) = 7.396228055897900e-01; pointCoord(708, 1) =2.349605040869134e-01; pointCoord(708, 2) =1.281223507749970e-04;
pointCoord(709, 0) = 2.396228055897900e-01; pointCoord(709, 1) =4.849605040869134e-01; pointCoord(709, 2) =1.281223507749970e-04;
pointCoord(710, 0) = 4.896228055897899e-01; pointCoord(710, 1) =4.849605040869134e-01; pointCoord(710, 2) =1.281223507749970e-04;
pointCoord(711, 0) = 2.396228055897900e-01; pointCoord(711, 1) =7.349605040869134e-01; pointCoord(711, 2) =1.281223507749970e-04;
pointCoord(712, 0) = 2.349605040869134e-01; pointCoord(712, 1) =2.396228055897900e-01; pointCoord(712, 2) =1.281223507749970e-04;
pointCoord(713, 0) = 4.849605040869134e-01; pointCoord(713, 1) =2.396228055897900e-01; pointCoord(713, 2) =1.281223507749970e-04;
pointCoord(714, 0) = 7.349605040869134e-01; pointCoord(714, 1) =2.396228055897900e-01; pointCoord(714, 2) =1.281223507749970e-04;
pointCoord(715, 0) = 2.349605040869134e-01; pointCoord(715, 1) =4.896228055897899e-01; pointCoord(715, 2) =1.281223507749970e-04;
pointCoord(716, 0) = 4.849605040869134e-01; pointCoord(716, 1) =4.896228055897899e-01; pointCoord(716, 2) =1.281223507749970e-04;
pointCoord(717, 0) = 2.349605040869134e-01; pointCoord(717, 1) =7.396228055897900e-01; pointCoord(717, 2) =1.281223507749970e-04;
pointCoord(718, 0) = 2.306130075795021e-01; pointCoord(718, 1) =2.439703020972012e-01; pointCoord(718, 2) =1.108205956969521e-04;
pointCoord(719, 0) = 4.806130075795021e-01; pointCoord(719, 1) =2.439703020972012e-01; pointCoord(719, 2) =1.108205956969521e-04;
pointCoord(720, 0) = 7.306130075795021e-01; pointCoord(720, 1) =2.439703020972012e-01; pointCoord(720, 2) =1.108205956969521e-04;
pointCoord(721, 0) = 2.306130075795021e-01; pointCoord(721, 1) =4.939703020972012e-01; pointCoord(721, 2) =1.108205956969521e-04;
pointCoord(722, 0) = 4.806130075795021e-01; pointCoord(722, 1) =4.939703020972012e-01; pointCoord(722, 2) =1.108205956969521e-04;
pointCoord(723, 0) = 2.306130075795021e-01; pointCoord(723, 1) =7.439703020972013e-01; pointCoord(723, 2) =1.108205956969521e-04;
pointCoord(724, 0) = 2.271673422646648e-01; pointCoord(724, 1) =2.474159674120386e-01; pointCoord(724, 2) =7.855873961758663e-05;
pointCoord(725, 0) = 4.771673422646648e-01; pointCoord(725, 1) =2.474159674120386e-01; pointCoord(725, 2) =7.855873961758663e-05;
pointCoord(726, 0) = 7.271673422646648e-01; pointCoord(726, 1) =2.474159674120386e-01; pointCoord(726, 2) =7.855873961758663e-05;
pointCoord(727, 0) = 2.271673422646648e-01; pointCoord(727, 1) =4.974159674120385e-01; pointCoord(727, 2) =7.855873961758663e-05;
pointCoord(728, 0) = 4.771673422646648e-01; pointCoord(728, 1) =4.974159674120385e-01; pointCoord(728, 2) =7.855873961758663e-05;
pointCoord(729, 0) = 2.271673422646648e-01; pointCoord(729, 1) =7.474159674120385e-01; pointCoord(729, 2) =7.855873961758663e-05;
pointCoord(730, 0) = 2.250879598867512e-01; pointCoord(730, 1) =2.494953497899521e-01; pointCoord(730, 2) =3.576018181520091e-05;
pointCoord(731, 0) = 4.750879598867512e-01; pointCoord(731, 1) =2.494953497899521e-01; pointCoord(731, 2) =3.576018181520091e-05;
pointCoord(732, 0) = 7.250879598867512e-01; pointCoord(732, 1) =2.494953497899521e-01; pointCoord(732, 2) =3.576018181520091e-05;
pointCoord(733, 0) = 2.250879598867512e-01; pointCoord(733, 1) =4.994953497899521e-01; pointCoord(733, 2) =3.576018181520091e-05;
pointCoord(734, 0) = 4.750879598867512e-01; pointCoord(734, 1) =4.994953497899521e-01; pointCoord(734, 2) =3.576018181520091e-05;
pointCoord(735, 0) = 2.250879598867512e-01; pointCoord(735, 1) =7.494953497899521e-01; pointCoord(735, 2) =3.576018181520091e-05;
pointCoord(736, 0) = 2.488224264944068e-01; pointCoord(736, 1) =1.918691247451343e-01; pointCoord(736, 2) =1.177126830861867e-04;
pointCoord(737, 0) = 4.988224264944068e-01; pointCoord(737, 1) =1.918691247451343e-01; pointCoord(737, 2) =1.177126830861867e-04;
pointCoord(738, 0) = 7.488224264944068e-01; pointCoord(738, 1) =1.918691247451343e-01; pointCoord(738, 2) =1.177126830861867e-04;
pointCoord(739, 0) = 2.488224264944068e-01; pointCoord(739, 1) =4.418691247451343e-01; pointCoord(739, 2) =1.177126830861867e-04;
pointCoord(740, 0) = 4.988224264944068e-01; pointCoord(740, 1) =4.418691247451343e-01; pointCoord(740, 2) =1.177126830861867e-04;
pointCoord(741, 0) = 2.488224264944068e-01; pointCoord(741, 1) =6.918691247451343e-01; pointCoord(741, 2) =1.177126830861867e-04;
pointCoord(742, 0) = 2.439703020972012e-01; pointCoord(742, 1) =1.967212491423399e-01; pointCoord(742, 2) =2.585937640933464e-04;
pointCoord(743, 0) = 4.939703020972012e-01; pointCoord(743, 1) =1.967212491423399e-01; pointCoord(743, 2) =2.585937640933464e-04;
pointCoord(744, 0) = 7.439703020972013e-01; pointCoord(744, 1) =1.967212491423399e-01; pointCoord(744, 2) =2.585937640933464e-04;
pointCoord(745, 0) = 2.439703020972012e-01; pointCoord(745, 1) =4.467212491423399e-01; pointCoord(745, 2) =2.585937640933464e-04;
pointCoord(746, 0) = 4.939703020972012e-01; pointCoord(746, 1) =4.467212491423399e-01; pointCoord(746, 2) =2.585937640933464e-04;
pointCoord(747, 0) = 2.439703020972012e-01; pointCoord(747, 1) =6.967212491423399e-01; pointCoord(747, 2) =2.585937640933464e-04;
pointCoord(748, 0) = 2.359300316225121e-01; pointCoord(748, 1) =2.047615196170290e-01; pointCoord(748, 2) =3.647909210336454e-04;
pointCoord(749, 0) = 4.859300316225121e-01; pointCoord(749, 1) =2.047615196170290e-01; pointCoord(749, 2) =3.647909210336454e-04;
pointCoord(750, 0) = 7.359300316225121e-01; pointCoord(750, 1) =2.047615196170290e-01; pointCoord(750, 2) =3.647909210336454e-04;
pointCoord(751, 0) = 2.359300316225121e-01; pointCoord(751, 1) =4.547615196170290e-01; pointCoord(751, 2) =3.647909210336454e-04;
pointCoord(752, 0) = 4.859300316225121e-01; pointCoord(752, 1) =4.547615196170290e-01; pointCoord(752, 2) =3.647909210336454e-04;
pointCoord(753, 0) = 2.359300316225121e-01; pointCoord(753, 1) =7.047615196170290e-01; pointCoord(753, 2) =3.647909210336454e-04;
pointCoord(754, 0) = 2.257853876674437e-01; pointCoord(754, 1) =2.149061635720974e-01; pointCoord(754, 2) =4.217435400908277e-04;
pointCoord(755, 0) = 4.757853876674437e-01; pointCoord(755, 1) =2.149061635720974e-01; pointCoord(755, 2) =4.217435400908277e-04;
pointCoord(756, 0) = 7.257853876674437e-01; pointCoord(756, 1) =2.149061635720974e-01; pointCoord(756, 2) =4.217435400908277e-04;
pointCoord(757, 0) = 2.257853876674437e-01; pointCoord(757, 1) =4.649061635720974e-01; pointCoord(757, 2) =4.217435400908277e-04;
pointCoord(758, 0) = 4.757853876674437e-01; pointCoord(758, 1) =4.649061635720974e-01; pointCoord(758, 2) =4.217435400908277e-04;
pointCoord(759, 0) = 2.257853876674437e-01; pointCoord(759, 1) =7.149061635720974e-01; pointCoord(759, 2) =4.217435400908277e-04;
pointCoord(760, 0) = 2.149061635720974e-01; pointCoord(760, 1) =2.257853876674437e-01; pointCoord(760, 2) =4.217435400908277e-04;
pointCoord(761, 0) = 4.649061635720974e-01; pointCoord(761, 1) =2.257853876674437e-01; pointCoord(761, 2) =4.217435400908277e-04;
pointCoord(762, 0) = 7.149061635720974e-01; pointCoord(762, 1) =2.257853876674437e-01; pointCoord(762, 2) =4.217435400908277e-04;
pointCoord(763, 0) = 2.149061635720974e-01; pointCoord(763, 1) =4.757853876674437e-01; pointCoord(763, 2) =4.217435400908277e-04;
pointCoord(764, 0) = 4.649061635720974e-01; pointCoord(764, 1) =4.757853876674437e-01; pointCoord(764, 2) =4.217435400908277e-04;
pointCoord(765, 0) = 2.149061635720974e-01; pointCoord(765, 1) =7.257853876674437e-01; pointCoord(765, 2) =4.217435400908277e-04;
pointCoord(766, 0) = 2.047615196170290e-01; pointCoord(766, 1) =2.359300316225121e-01; pointCoord(766, 2) =3.647909210336454e-04;
pointCoord(767, 0) = 4.547615196170290e-01; pointCoord(767, 1) =2.359300316225121e-01; pointCoord(767, 2) =3.647909210336454e-04;
pointCoord(768, 0) = 7.047615196170290e-01; pointCoord(768, 1) =2.359300316225121e-01; pointCoord(768, 2) =3.647909210336454e-04;
pointCoord(769, 0) = 2.047615196170290e-01; pointCoord(769, 1) =4.859300316225121e-01; pointCoord(769, 2) =3.647909210336454e-04;
pointCoord(770, 0) = 4.547615196170290e-01; pointCoord(770, 1) =4.859300316225121e-01; pointCoord(770, 2) =3.647909210336454e-04;
pointCoord(771, 0) = 2.047615196170290e-01; pointCoord(771, 1) =7.359300316225121e-01; pointCoord(771, 2) =3.647909210336454e-04;
pointCoord(772, 0) = 1.967212491423399e-01; pointCoord(772, 1) =2.439703020972012e-01; pointCoord(772, 2) =2.585937640933464e-04;
pointCoord(773, 0) = 4.467212491423399e-01; pointCoord(773, 1) =2.439703020972012e-01; pointCoord(773, 2) =2.585937640933464e-04;
pointCoord(774, 0) = 6.967212491423399e-01; pointCoord(774, 1) =2.439703020972012e-01; pointCoord(774, 2) =2.585937640933464e-04;
pointCoord(775, 0) = 1.967212491423399e-01; pointCoord(775, 1) =4.939703020972012e-01; pointCoord(775, 2) =2.585937640933464e-04;
pointCoord(776, 0) = 4.467212491423399e-01; pointCoord(776, 1) =4.939703020972012e-01; pointCoord(776, 2) =2.585937640933464e-04;
pointCoord(777, 0) = 1.967212491423399e-01; pointCoord(777, 1) =7.439703020972013e-01; pointCoord(777, 2) =2.585937640933464e-04;
pointCoord(778, 0) = 1.918691247451343e-01; pointCoord(778, 1) =2.488224264944068e-01; pointCoord(778, 2) =1.177126830861867e-04;
pointCoord(779, 0) = 4.418691247451343e-01; pointCoord(779, 1) =2.488224264944068e-01; pointCoord(779, 2) =1.177126830861867e-04;
pointCoord(780, 0) = 6.918691247451343e-01; pointCoord(780, 1) =2.488224264944068e-01; pointCoord(780, 2) =1.177126830861867e-04;
pointCoord(781, 0) = 1.918691247451343e-01; pointCoord(781, 1) =4.988224264944068e-01; pointCoord(781, 2) =1.177126830861867e-04;
pointCoord(782, 0) = 4.418691247451343e-01; pointCoord(782, 1) =4.988224264944068e-01; pointCoord(782, 2) =1.177126830861867e-04;
pointCoord(783, 0) = 1.918691247451343e-01; pointCoord(783, 1) =7.488224264944068e-01; pointCoord(783, 2) =1.177126830861867e-04;
pointCoord(784, 0) = 2.479733795296476e-01; pointCoord(784, 1) =1.499559507823086e-01; pointCoord(784, 2) =2.342135820693354e-04;
pointCoord(785, 0) = 4.979733795296476e-01; pointCoord(785, 1) =1.499559507823086e-01; pointCoord(785, 2) =2.342135820693354e-04;
pointCoord(786, 0) = 7.479733795296476e-01; pointCoord(786, 1) =1.499559507823086e-01; pointCoord(786, 2) =2.342135820693354e-04;
pointCoord(787, 0) = 2.479733795296476e-01; pointCoord(787, 1) =3.999559507823086e-01; pointCoord(787, 2) =2.342135820693354e-04;
pointCoord(788, 0) = 4.979733795296476e-01; pointCoord(788, 1) =3.999559507823086e-01; pointCoord(788, 2) =2.342135820693354e-04;
pointCoord(789, 0) = 2.479733795296476e-01; pointCoord(789, 1) =6.499559507823086e-01; pointCoord(789, 2) =2.342135820693354e-04;
pointCoord(790, 0) = 2.396228055897900e-01; pointCoord(790, 1) =1.583065247221663e-01; pointCoord(790, 2) =5.145254547018528e-04;
pointCoord(791, 0) = 4.896228055897899e-01; pointCoord(791, 1) =1.583065247221663e-01; pointCoord(791, 2) =5.145254547018528e-04;
pointCoord(792, 0) = 7.396228055897900e-01; pointCoord(792, 1) =1.583065247221663e-01; pointCoord(792, 2) =5.145254547018528e-04;
pointCoord(793, 0) = 2.396228055897900e-01; pointCoord(793, 1) =4.083065247221663e-01; pointCoord(793, 2) =5.145254547018528e-04;
pointCoord(794, 0) = 4.896228055897899e-01; pointCoord(794, 1) =4.083065247221663e-01; pointCoord(794, 2) =5.145254547018528e-04;
pointCoord(795, 0) = 2.396228055897900e-01; pointCoord(795, 1) =6.583065247221662e-01; pointCoord(795, 2) =5.145254547018528e-04;
pointCoord(796, 0) = 2.257853876674437e-01; pointCoord(796, 1) =1.721439426445125e-01; pointCoord(796, 2) =7.258265301718213e-04;
pointCoord(797, 0) = 4.757853876674437e-01; pointCoord(797, 1) =1.721439426445125e-01; pointCoord(797, 2) =7.258265301718213e-04;
pointCoord(798, 0) = 7.257853876674437e-01; pointCoord(798, 1) =1.721439426445125e-01; pointCoord(798, 2) =7.258265301718213e-04;
pointCoord(799, 0) = 2.257853876674437e-01; pointCoord(799, 1) =4.221439426445125e-01; pointCoord(799, 2) =7.258265301718213e-04;
pointCoord(800, 0) = 4.757853876674437e-01; pointCoord(800, 1) =4.221439426445125e-01; pointCoord(800, 2) =7.258265301718213e-04;
pointCoord(801, 0) = 2.257853876674437e-01; pointCoord(801, 1) =6.721439426445125e-01; pointCoord(801, 2) =7.258265301718213e-04;
pointCoord(802, 0) = 2.083263135577370e-01; pointCoord(802, 1) =1.896030167542192e-01; pointCoord(802, 2) =8.391454739584170e-04;
pointCoord(803, 0) = 4.583263135577371e-01; pointCoord(803, 1) =1.896030167542192e-01; pointCoord(803, 2) =8.391454739584170e-04;
pointCoord(804, 0) = 7.083263135577370e-01; pointCoord(804, 1) =1.896030167542192e-01; pointCoord(804, 2) =8.391454739584170e-04;
pointCoord(805, 0) = 2.083263135577370e-01; pointCoord(805, 1) =4.396030167542192e-01; pointCoord(805, 2) =8.391454739584170e-04;
pointCoord(806, 0) = 4.583263135577371e-01; pointCoord(806, 1) =4.396030167542192e-01; pointCoord(806, 2) =8.391454739584170e-04;
pointCoord(807, 0) = 2.083263135577370e-01; pointCoord(807, 1) =6.896030167542192e-01; pointCoord(807, 2) =8.391454739584170e-04;
pointCoord(808, 0) = 1.896030167542192e-01; pointCoord(808, 1) =2.083263135577370e-01; pointCoord(808, 2) =8.391454739584170e-04;
pointCoord(809, 0) = 4.396030167542192e-01; pointCoord(809, 1) =2.083263135577370e-01; pointCoord(809, 2) =8.391454739584170e-04;
pointCoord(810, 0) = 6.896030167542192e-01; pointCoord(810, 1) =2.083263135577370e-01; pointCoord(810, 2) =8.391454739584170e-04;
pointCoord(811, 0) = 1.896030167542192e-01; pointCoord(811, 1) =4.583263135577371e-01; pointCoord(811, 2) =8.391454739584170e-04;
pointCoord(812, 0) = 4.396030167542192e-01; pointCoord(812, 1) =4.583263135577371e-01; pointCoord(812, 2) =8.391454739584170e-04;
pointCoord(813, 0) = 1.896030167542192e-01; pointCoord(813, 1) =7.083263135577370e-01; pointCoord(813, 2) =8.391454739584170e-04;
pointCoord(814, 0) = 1.721439426445125e-01; pointCoord(814, 1) =2.257853876674437e-01; pointCoord(814, 2) =7.258265301718213e-04;
pointCoord(815, 0) = 4.221439426445125e-01; pointCoord(815, 1) =2.257853876674437e-01; pointCoord(815, 2) =7.258265301718213e-04;
pointCoord(816, 0) = 6.721439426445125e-01; pointCoord(816, 1) =2.257853876674437e-01; pointCoord(816, 2) =7.258265301718213e-04;
pointCoord(817, 0) = 1.721439426445125e-01; pointCoord(817, 1) =4.757853876674437e-01; pointCoord(817, 2) =7.258265301718213e-04;
pointCoord(818, 0) = 4.221439426445125e-01; pointCoord(818, 1) =4.757853876674437e-01; pointCoord(818, 2) =7.258265301718213e-04;
pointCoord(819, 0) = 1.721439426445125e-01; pointCoord(819, 1) =7.257853876674437e-01; pointCoord(819, 2) =7.258265301718213e-04;
pointCoord(820, 0) = 1.583065247221663e-01; pointCoord(820, 1) =2.396228055897900e-01; pointCoord(820, 2) =5.145254547018528e-04;
pointCoord(821, 0) = 4.083065247221663e-01; pointCoord(821, 1) =2.396228055897900e-01; pointCoord(821, 2) =5.145254547018528e-04;
pointCoord(822, 0) = 6.583065247221662e-01; pointCoord(822, 1) =2.396228055897900e-01; pointCoord(822, 2) =5.145254547018528e-04;
pointCoord(823, 0) = 1.583065247221663e-01; pointCoord(823, 1) =4.896228055897899e-01; pointCoord(823, 2) =5.145254547018528e-04;
pointCoord(824, 0) = 4.083065247221663e-01; pointCoord(824, 1) =4.896228055897899e-01; pointCoord(824, 2) =5.145254547018528e-04;
pointCoord(825, 0) = 1.583065247221663e-01; pointCoord(825, 1) =7.396228055897900e-01; pointCoord(825, 2) =5.145254547018528e-04;
pointCoord(826, 0) = 1.499559507823086e-01; pointCoord(826, 1) =2.479733795296476e-01; pointCoord(826, 2) =2.342135820693354e-04;
pointCoord(827, 0) = 3.999559507823086e-01; pointCoord(827, 1) =2.479733795296476e-01; pointCoord(827, 2) =2.342135820693354e-04;
pointCoord(828, 0) = 6.499559507823086e-01; pointCoord(828, 1) =2.479733795296476e-01; pointCoord(828, 2) =2.342135820693354e-04;
pointCoord(829, 0) = 1.499559507823086e-01; pointCoord(829, 1) =4.979733795296476e-01; pointCoord(829, 2) =2.342135820693354e-04;
pointCoord(830, 0) = 3.999559507823086e-01; pointCoord(830, 1) =4.979733795296476e-01; pointCoord(830, 2) =2.342135820693354e-04;
pointCoord(831, 0) = 1.499559507823086e-01; pointCoord(831, 1) =7.479733795296476e-01; pointCoord(831, 2) =2.342135820693354e-04;
pointCoord(832, 0) = 2.470628525325444e-01; pointCoord(832, 1) =1.050078171554993e-01; pointCoord(832, 2) =3.394418636751593e-04;
pointCoord(833, 0) = 4.970628525325445e-01; pointCoord(833, 1) =1.050078171554993e-01; pointCoord(833, 2) =3.394418636751593e-04;
pointCoord(834, 0) = 7.470628525325445e-01; pointCoord(834, 1) =1.050078171554993e-01; pointCoord(834, 2) =3.394418636751593e-04;
pointCoord(835, 0) = 2.470628525325444e-01; pointCoord(835, 1) =3.550078171554993e-01; pointCoord(835, 2) =3.394418636751593e-04;
pointCoord(836, 0) = 4.970628525325445e-01; pointCoord(836, 1) =3.550078171554993e-01; pointCoord(836, 2) =3.394418636751593e-04;
pointCoord(837, 0) = 2.470628525325444e-01; pointCoord(837, 1) =6.050078171554993e-01; pointCoord(837, 2) =3.394418636751593e-04;
pointCoord(838, 0) = 2.349605040869134e-01; pointCoord(838, 1) =1.171101656011304e-01; pointCoord(838, 2) =7.456932160347678e-04;
pointCoord(839, 0) = 4.849605040869134e-01; pointCoord(839, 1) =1.171101656011304e-01; pointCoord(839, 2) =7.456932160347678e-04;
pointCoord(840, 0) = 7.349605040869134e-01; pointCoord(840, 1) =1.171101656011304e-01; pointCoord(840, 2) =7.456932160347678e-04;
pointCoord(841, 0) = 2.349605040869134e-01; pointCoord(841, 1) =3.671101656011304e-01; pointCoord(841, 2) =7.456932160347678e-04;
pointCoord(842, 0) = 4.849605040869134e-01; pointCoord(842, 1) =3.671101656011304e-01; pointCoord(842, 2) =7.456932160347678e-04;
pointCoord(843, 0) = 2.349605040869134e-01; pointCoord(843, 1) =6.171101656011304e-01; pointCoord(843, 2) =7.456932160347678e-04;
pointCoord(844, 0) = 2.149061635720974e-01; pointCoord(844, 1) =1.371645061159464e-01; pointCoord(844, 2) =1.051928363545806e-03;
pointCoord(845, 0) = 4.649061635720974e-01; pointCoord(845, 1) =1.371645061159464e-01; pointCoord(845, 2) =1.051928363545806e-03;
pointCoord(846, 0) = 7.149061635720974e-01; pointCoord(846, 1) =1.371645061159464e-01; pointCoord(846, 2) =1.051928363545806e-03;
pointCoord(847, 0) = 2.149061635720974e-01; pointCoord(847, 1) =3.871645061159464e-01; pointCoord(847, 2) =1.051928363545806e-03;
pointCoord(848, 0) = 4.649061635720974e-01; pointCoord(848, 1) =3.871645061159464e-01; pointCoord(848, 2) =1.051928363545806e-03;
pointCoord(849, 0) = 2.149061635720974e-01; pointCoord(849, 1) =6.371645061159464e-01; pointCoord(849, 2) =1.051928363545806e-03;
pointCoord(850, 0) = 1.896030167542192e-01; pointCoord(850, 1) =1.624676529338246e-01; pointCoord(850, 2) =1.216159631129748e-03;
pointCoord(851, 0) = 4.396030167542192e-01; pointCoord(851, 1) =1.624676529338246e-01; pointCoord(851, 2) =1.216159631129748e-03;
pointCoord(852, 0) = 6.896030167542192e-01; pointCoord(852, 1) =1.624676529338246e-01; pointCoord(852, 2) =1.216159631129748e-03;
pointCoord(853, 0) = 1.896030167542192e-01; pointCoord(853, 1) =4.124676529338246e-01; pointCoord(853, 2) =1.216159631129748e-03;
pointCoord(854, 0) = 4.396030167542192e-01; pointCoord(854, 1) =4.124676529338246e-01; pointCoord(854, 2) =1.216159631129748e-03;
pointCoord(855, 0) = 1.896030167542192e-01; pointCoord(855, 1) =6.624676529338246e-01; pointCoord(855, 2) =1.216159631129748e-03;
pointCoord(856, 0) = 1.624676529338246e-01; pointCoord(856, 1) =1.896030167542192e-01; pointCoord(856, 2) =1.216159631129748e-03;
pointCoord(857, 0) = 4.124676529338246e-01; pointCoord(857, 1) =1.896030167542192e-01; pointCoord(857, 2) =1.216159631129748e-03;
pointCoord(858, 0) = 6.624676529338246e-01; pointCoord(858, 1) =1.896030167542192e-01; pointCoord(858, 2) =1.216159631129748e-03;
pointCoord(859, 0) = 1.624676529338246e-01; pointCoord(859, 1) =4.396030167542192e-01; pointCoord(859, 2) =1.216159631129748e-03;
pointCoord(860, 0) = 4.124676529338246e-01; pointCoord(860, 1) =4.396030167542192e-01; pointCoord(860, 2) =1.216159631129748e-03;
pointCoord(861, 0) = 1.624676529338246e-01; pointCoord(861, 1) =6.896030167542192e-01; pointCoord(861, 2) =1.216159631129748e-03;
pointCoord(862, 0) = 1.371645061159464e-01; pointCoord(862, 1) =2.149061635720974e-01; pointCoord(862, 2) =1.051928363545806e-03;
pointCoord(863, 0) = 3.871645061159464e-01; pointCoord(863, 1) =2.149061635720974e-01; pointCoord(863, 2) =1.051928363545806e-03;
pointCoord(864, 0) = 6.371645061159464e-01; pointCoord(864, 1) =2.149061635720974e-01; pointCoord(864, 2) =1.051928363545806e-03;
pointCoord(865, 0) = 1.371645061159464e-01; pointCoord(865, 1) =4.649061635720974e-01; pointCoord(865, 2) =1.051928363545806e-03;
pointCoord(866, 0) = 3.871645061159464e-01; pointCoord(866, 1) =4.649061635720974e-01; pointCoord(866, 2) =1.051928363545806e-03;
pointCoord(867, 0) = 1.371645061159464e-01; pointCoord(867, 1) =7.149061635720974e-01; pointCoord(867, 2) =1.051928363545806e-03;
pointCoord(868, 0) = 1.171101656011304e-01; pointCoord(868, 1) =2.349605040869134e-01; pointCoord(868, 2) =7.456932160347678e-04;
pointCoord(869, 0) = 3.671101656011304e-01; pointCoord(869, 1) =2.349605040869134e-01; pointCoord(869, 2) =7.456932160347678e-04;
pointCoord(870, 0) = 6.171101656011304e-01; pointCoord(870, 1) =2.349605040869134e-01; pointCoord(870, 2) =7.456932160347678e-04;
pointCoord(871, 0) = 1.171101656011304e-01; pointCoord(871, 1) =4.849605040869134e-01; pointCoord(871, 2) =7.456932160347678e-04;
pointCoord(872, 0) = 3.671101656011304e-01; pointCoord(872, 1) =4.849605040869134e-01; pointCoord(872, 2) =7.456932160347678e-04;
pointCoord(873, 0) = 1.171101656011304e-01; pointCoord(873, 1) =7.349605040869134e-01; pointCoord(873, 2) =7.456932160347678e-04;
pointCoord(874, 0) = 1.050078171554993e-01; pointCoord(874, 1) =2.470628525325444e-01; pointCoord(874, 2) =3.394418636751593e-04;
pointCoord(875, 0) = 3.550078171554993e-01; pointCoord(875, 1) =2.470628525325444e-01; pointCoord(875, 2) =3.394418636751593e-04;
pointCoord(876, 0) = 6.050078171554993e-01; pointCoord(876, 1) =2.470628525325444e-01; pointCoord(876, 2) =3.394418636751593e-04;
pointCoord(877, 0) = 1.050078171554993e-01; pointCoord(877, 1) =4.970628525325445e-01; pointCoord(877, 2) =3.394418636751593e-04;
pointCoord(878, 0) = 3.550078171554993e-01; pointCoord(878, 1) =4.970628525325445e-01; pointCoord(878, 2) =3.394418636751593e-04;
pointCoord(879, 0) = 1.050078171554993e-01; pointCoord(879, 1) =7.470628525325445e-01; pointCoord(879, 2) =3.394418636751593e-04;
pointCoord(880, 0) = 2.462138055677852e-01; pointCoord(880, 1) =6.309464319267366e-02; pointCoord(880, 2) =3.784758260822833e-04;
pointCoord(881, 0) = 4.962138055677852e-01; pointCoord(881, 1) =6.309464319267366e-02; pointCoord(881, 2) =3.784758260822833e-04;
pointCoord(882, 0) = 7.462138055677852e-01; pointCoord(882, 1) =6.309464319267366e-02; pointCoord(882, 2) =3.784758260822833e-04;
pointCoord(883, 0) = 2.462138055677852e-01; pointCoord(883, 1) =3.130946431926737e-01; pointCoord(883, 2) =3.784758260822833e-04;
pointCoord(884, 0) = 4.962138055677852e-01; pointCoord(884, 1) =3.130946431926737e-01; pointCoord(884, 2) =3.784758260822833e-04;
pointCoord(885, 0) = 2.462138055677852e-01; pointCoord(885, 1) =5.630946431926737e-01; pointCoord(885, 2) =3.784758260822833e-04;
pointCoord(886, 0) = 2.306130075795021e-01; pointCoord(886, 1) =7.869544118095678e-02; pointCoord(886, 2) =8.314438675507627e-04;
pointCoord(887, 0) = 4.806130075795021e-01; pointCoord(887, 1) =7.869544118095678e-02; pointCoord(887, 2) =8.314438675507627e-04;
pointCoord(888, 0) = 7.306130075795021e-01; pointCoord(888, 1) =7.869544118095678e-02; pointCoord(888, 2) =8.314438675507627e-04;
pointCoord(889, 0) = 2.306130075795021e-01; pointCoord(889, 1) =3.286954411809568e-01; pointCoord(889, 2) =8.314438675507627e-04;
pointCoord(890, 0) = 4.806130075795021e-01; pointCoord(890, 1) =3.286954411809568e-01; pointCoord(890, 2) =8.314438675507627e-04;
pointCoord(891, 0) = 2.306130075795021e-01; pointCoord(891, 1) =5.786954411809568e-01; pointCoord(891, 2) =8.314438675507627e-04;
pointCoord(892, 0) = 2.047615196170290e-01; pointCoord(892, 1) =1.045469291434299e-01; pointCoord(892, 2) =1.172894386278138e-03;
pointCoord(893, 0) = 4.547615196170290e-01; pointCoord(893, 1) =1.045469291434299e-01; pointCoord(893, 2) =1.172894386278138e-03;
pointCoord(894, 0) = 7.047615196170290e-01; pointCoord(894, 1) =1.045469291434299e-01; pointCoord(894, 2) =1.172894386278138e-03;
pointCoord(895, 0) = 2.047615196170290e-01; pointCoord(895, 1) =3.545469291434299e-01; pointCoord(895, 2) =1.172894386278138e-03;
pointCoord(896, 0) = 4.547615196170290e-01; pointCoord(896, 1) =3.545469291434299e-01; pointCoord(896, 2) =1.172894386278138e-03;
pointCoord(897, 0) = 2.047615196170290e-01; pointCoord(897, 1) =6.045469291434299e-01; pointCoord(897, 2) =1.172894386278138e-03;
pointCoord(898, 0) = 1.721439426445125e-01; pointCoord(898, 1) =1.371645061159464e-01; pointCoord(898, 2) =1.356011353626800e-03;
pointCoord(899, 0) = 4.221439426445125e-01; pointCoord(899, 1) =1.371645061159464e-01; pointCoord(899, 2) =1.356011353626800e-03;
pointCoord(900, 0) = 6.721439426445125e-01; pointCoord(900, 1) =1.371645061159464e-01; pointCoord(900, 2) =1.356011353626800e-03;
pointCoord(901, 0) = 1.721439426445125e-01; pointCoord(901, 1) =3.871645061159464e-01; pointCoord(901, 2) =1.356011353626800e-03;
pointCoord(902, 0) = 4.221439426445125e-01; pointCoord(902, 1) =3.871645061159464e-01; pointCoord(902, 2) =1.356011353626800e-03;
pointCoord(903, 0) = 1.721439426445125e-01; pointCoord(903, 1) =6.371645061159464e-01; pointCoord(903, 2) =1.356011353626800e-03;
pointCoord(904, 0) = 1.371645061159464e-01; pointCoord(904, 1) =1.721439426445125e-01; pointCoord(904, 2) =1.356011353626800e-03;
pointCoord(905, 0) = 3.871645061159464e-01; pointCoord(905, 1) =1.721439426445125e-01; pointCoord(905, 2) =1.356011353626800e-03;
pointCoord(906, 0) = 6.371645061159464e-01; pointCoord(906, 1) =1.721439426445125e-01; pointCoord(906, 2) =1.356011353626800e-03;
pointCoord(907, 0) = 1.371645061159464e-01; pointCoord(907, 1) =4.221439426445125e-01; pointCoord(907, 2) =1.356011353626800e-03;
pointCoord(908, 0) = 3.871645061159464e-01; pointCoord(908, 1) =4.221439426445125e-01; pointCoord(908, 2) =1.356011353626800e-03;
pointCoord(909, 0) = 1.371645061159464e-01; pointCoord(909, 1) =6.721439426445125e-01; pointCoord(909, 2) =1.356011353626800e-03;
pointCoord(910, 0) = 1.045469291434299e-01; pointCoord(910, 1) =2.047615196170290e-01; pointCoord(910, 2) =1.172894386278138e-03;
pointCoord(911, 0) = 3.545469291434299e-01; pointCoord(911, 1) =2.047615196170290e-01; pointCoord(911, 2) =1.172894386278138e-03;
pointCoord(912, 0) = 6.045469291434299e-01; pointCoord(912, 1) =2.047615196170290e-01; pointCoord(912, 2) =1.172894386278138e-03;
pointCoord(913, 0) = 1.045469291434299e-01; pointCoord(913, 1) =4.547615196170290e-01; pointCoord(913, 2) =1.172894386278138e-03;
pointCoord(914, 0) = 3.545469291434299e-01; pointCoord(914, 1) =4.547615196170290e-01; pointCoord(914, 2) =1.172894386278138e-03;
pointCoord(915, 0) = 1.045469291434299e-01; pointCoord(915, 1) =7.047615196170290e-01; pointCoord(915, 2) =1.172894386278138e-03;
pointCoord(916, 0) = 7.869544118095678e-02; pointCoord(916, 1) =2.306130075795021e-01; pointCoord(916, 2) =8.314438675507627e-04;
pointCoord(917, 0) = 3.286954411809568e-01; pointCoord(917, 1) =2.306130075795021e-01; pointCoord(917, 2) =8.314438675507627e-04;
pointCoord(918, 0) = 5.786954411809568e-01; pointCoord(918, 1) =2.306130075795021e-01; pointCoord(918, 2) =8.314438675507627e-04;
pointCoord(919, 0) = 7.869544118095678e-02; pointCoord(919, 1) =4.806130075795021e-01; pointCoord(919, 2) =8.314438675507627e-04;
pointCoord(920, 0) = 3.286954411809568e-01; pointCoord(920, 1) =4.806130075795021e-01; pointCoord(920, 2) =8.314438675507627e-04;
pointCoord(921, 0) = 7.869544118095678e-02; pointCoord(921, 1) =7.306130075795021e-01; pointCoord(921, 2) =8.314438675507627e-04;
pointCoord(922, 0) = 6.309464319267366e-02; pointCoord(922, 1) =2.462138055677852e-01; pointCoord(922, 2) =3.784758260822833e-04;
pointCoord(923, 0) = 3.130946431926737e-01; pointCoord(923, 1) =2.462138055677852e-01; pointCoord(923, 2) =3.784758260822833e-04;
pointCoord(924, 0) = 5.630946431926737e-01; pointCoord(924, 1) =2.462138055677852e-01; pointCoord(924, 2) =3.784758260822833e-04;
pointCoord(925, 0) = 6.309464319267366e-02; pointCoord(925, 1) =4.962138055677852e-01; pointCoord(925, 2) =3.784758260822833e-04;
pointCoord(926, 0) = 3.130946431926737e-01; pointCoord(926, 1) =4.962138055677852e-01; pointCoord(926, 2) =3.784758260822833e-04;
pointCoord(927, 0) = 6.309464319267366e-02; pointCoord(927, 1) =7.462138055677852e-01; pointCoord(927, 2) =3.784758260822833e-04;
pointCoord(928, 0) = 2.455408822722399e-01; pointCoord(928, 1) =2.987580805105672e-02; pointCoord(928, 2) =3.159789840669076e-04;
pointCoord(929, 0) = 4.955408822722399e-01; pointCoord(929, 1) =2.987580805105672e-02; pointCoord(929, 2) =3.159789840669076e-04;
pointCoord(930, 0) = 7.455408822722399e-01; pointCoord(930, 1) =2.987580805105672e-02; pointCoord(930, 2) =3.159789840669076e-04;
pointCoord(931, 0) = 2.455408822722399e-01; pointCoord(931, 1) =2.798758080510567e-01; pointCoord(931, 2) =3.159789840669076e-04;
pointCoord(932, 0) = 4.955408822722399e-01; pointCoord(932, 1) =2.798758080510567e-01; pointCoord(932, 2) =3.159789840669076e-04;
pointCoord(933, 0) = 2.455408822722399e-01; pointCoord(933, 1) =5.298758080510567e-01; pointCoord(933, 2) =3.159789840669076e-04;
pointCoord(934, 0) = 2.271673422646648e-01; pointCoord(934, 1) =4.824934805863187e-02; pointCoord(934, 2) =6.941494554535528e-04;
pointCoord(935, 0) = 4.771673422646648e-01; pointCoord(935, 1) =4.824934805863187e-02; pointCoord(935, 2) =6.941494554535528e-04;
pointCoord(936, 0) = 7.271673422646648e-01; pointCoord(936, 1) =4.824934805863187e-02; pointCoord(936, 2) =6.941494554535528e-04;
pointCoord(937, 0) = 2.271673422646648e-01; pointCoord(937, 1) =2.982493480586319e-01; pointCoord(937, 2) =6.941494554535528e-04;
pointCoord(938, 0) = 4.771673422646648e-01; pointCoord(938, 1) =2.982493480586319e-01; pointCoord(938, 2) =6.941494554535528e-04;
pointCoord(939, 0) = 2.271673422646648e-01; pointCoord(939, 1) =5.482493480586319e-01; pointCoord(939, 2) =6.941494554535528e-04;
pointCoord(940, 0) = 1.967212491423399e-01; pointCoord(940, 1) =7.869544118095678e-02; pointCoord(940, 2) =9.792170359471568e-04;
pointCoord(941, 0) = 4.467212491423399e-01; pointCoord(941, 1) =7.869544118095678e-02; pointCoord(941, 2) =9.792170359471568e-04;
pointCoord(942, 0) = 6.967212491423399e-01; pointCoord(942, 1) =7.869544118095678e-02; pointCoord(942, 2) =9.792170359471568e-04;
pointCoord(943, 0) = 1.967212491423399e-01; pointCoord(943, 1) =3.286954411809568e-01; pointCoord(943, 2) =9.792170359471568e-04;
pointCoord(944, 0) = 4.467212491423399e-01; pointCoord(944, 1) =3.286954411809568e-01; pointCoord(944, 2) =9.792170359471568e-04;
pointCoord(945, 0) = 1.967212491423399e-01; pointCoord(945, 1) =5.786954411809568e-01; pointCoord(945, 2) =9.792170359471568e-04;
pointCoord(946, 0) = 1.583065247221663e-01; pointCoord(946, 1) =1.171101656011304e-01; pointCoord(946, 2) =1.132096319961624e-03;
pointCoord(947, 0) = 4.083065247221663e-01; pointCoord(947, 1) =1.171101656011304e-01; pointCoord(947, 2) =1.132096319961624e-03;
pointCoord(948, 0) = 6.583065247221662e-01; pointCoord(948, 1) =1.171101656011304e-01; pointCoord(948, 2) =1.132096319961624e-03;
pointCoord(949, 0) = 1.583065247221663e-01; pointCoord(949, 1) =3.671101656011304e-01; pointCoord(949, 2) =1.132096319961624e-03;
pointCoord(950, 0) = 4.083065247221663e-01; pointCoord(950, 1) =3.671101656011304e-01; pointCoord(950, 2) =1.132096319961624e-03;
pointCoord(951, 0) = 1.583065247221663e-01; pointCoord(951, 1) =6.171101656011304e-01; pointCoord(951, 2) =1.132096319961624e-03;
pointCoord(952, 0) = 1.171101656011304e-01; pointCoord(952, 1) =1.583065247221663e-01; pointCoord(952, 2) =1.132096319961624e-03;
pointCoord(953, 0) = 3.671101656011304e-01; pointCoord(953, 1) =1.583065247221663e-01; pointCoord(953, 2) =1.132096319961624e-03;
pointCoord(954, 0) = 6.171101656011304e-01; pointCoord(954, 1) =1.583065247221663e-01; pointCoord(954, 2) =1.132096319961624e-03;
pointCoord(955, 0) = 1.171101656011304e-01; pointCoord(955, 1) =4.083065247221663e-01; pointCoord(955, 2) =1.132096319961624e-03;
pointCoord(956, 0) = 3.671101656011304e-01; pointCoord(956, 1) =4.083065247221663e-01; pointCoord(956, 2) =1.132096319961624e-03;
pointCoord(957, 0) = 1.171101656011304e-01; pointCoord(957, 1) =6.583065247221662e-01; pointCoord(957, 2) =1.132096319961624e-03;
pointCoord(958, 0) = 7.869544118095678e-02; pointCoord(958, 1) =1.967212491423399e-01; pointCoord(958, 2) =9.792170359471568e-04;
pointCoord(959, 0) = 3.286954411809568e-01; pointCoord(959, 1) =1.967212491423399e-01; pointCoord(959, 2) =9.792170359471568e-04;
pointCoord(960, 0) = 5.786954411809568e-01; pointCoord(960, 1) =1.967212491423399e-01; pointCoord(960, 2) =9.792170359471568e-04;
pointCoord(961, 0) = 7.869544118095678e-02; pointCoord(961, 1) =4.467212491423399e-01; pointCoord(961, 2) =9.792170359471568e-04;
pointCoord(962, 0) = 3.286954411809568e-01; pointCoord(962, 1) =4.467212491423399e-01; pointCoord(962, 2) =9.792170359471568e-04;
pointCoord(963, 0) = 7.869544118095678e-02; pointCoord(963, 1) =6.967212491423399e-01; pointCoord(963, 2) =9.792170359471568e-04;
pointCoord(964, 0) = 4.824934805863187e-02; pointCoord(964, 1) =2.271673422646648e-01; pointCoord(964, 2) =6.941494554535528e-04;
pointCoord(965, 0) = 2.982493480586319e-01; pointCoord(965, 1) =2.271673422646648e-01; pointCoord(965, 2) =6.941494554535528e-04;
pointCoord(966, 0) = 5.482493480586319e-01; pointCoord(966, 1) =2.271673422646648e-01; pointCoord(966, 2) =6.941494554535528e-04;
pointCoord(967, 0) = 4.824934805863187e-02; pointCoord(967, 1) =4.771673422646648e-01; pointCoord(967, 2) =6.941494554535528e-04;
pointCoord(968, 0) = 2.982493480586319e-01; pointCoord(968, 1) =4.771673422646648e-01; pointCoord(968, 2) =6.941494554535528e-04;
pointCoord(969, 0) = 4.824934805863187e-02; pointCoord(969, 1) =7.271673422646648e-01; pointCoord(969, 2) =6.941494554535528e-04;
pointCoord(970, 0) = 2.987580805105672e-02; pointCoord(970, 1) =2.455408822722399e-01; pointCoord(970, 2) =3.159789840669076e-04;
pointCoord(971, 0) = 2.798758080510567e-01; pointCoord(971, 1) =2.455408822722399e-01; pointCoord(971, 2) =3.159789840669076e-04;
pointCoord(972, 0) = 5.298758080510567e-01; pointCoord(972, 1) =2.455408822722399e-01; pointCoord(972, 2) =3.159789840669076e-04;
pointCoord(973, 0) = 2.987580805105672e-02; pointCoord(973, 1) =4.955408822722399e-01; pointCoord(973, 2) =3.159789840669076e-04;
pointCoord(974, 0) = 2.798758080510567e-01; pointCoord(974, 1) =4.955408822722399e-01; pointCoord(974, 2) =3.159789840669076e-04;
pointCoord(975, 0) = 2.987580805105672e-02; pointCoord(975, 1) =7.455408822722399e-01; pointCoord(975, 2) =3.159789840669076e-04;
pointCoord(976, 0) = 2.451347880307537e-01; pointCoord(976, 1) =9.828979907054253e-03; pointCoord(976, 2) =1.569337084289531e-04;
pointCoord(977, 0) = 4.951347880307537e-01; pointCoord(977, 1) =9.828979907054253e-03; pointCoord(977, 2) =1.569337084289531e-04;
pointCoord(978, 0) = 7.451347880307537e-01; pointCoord(978, 1) =9.828979907054253e-03; pointCoord(978, 2) =1.569337084289531e-04;
pointCoord(979, 0) = 2.451347880307537e-01; pointCoord(979, 1) =2.598289799070542e-01; pointCoord(979, 2) =1.569337084289531e-04;
pointCoord(980, 0) = 4.951347880307537e-01; pointCoord(980, 1) =2.598289799070542e-01; pointCoord(980, 2) =1.569337084289531e-04;
pointCoord(981, 0) = 2.451347880307537e-01; pointCoord(981, 1) =5.098289799070542e-01; pointCoord(981, 2) =1.569337084289531e-04;
pointCoord(982, 0) = 2.250879598867512e-01; pointCoord(982, 1) =2.987580805105672e-02; pointCoord(982, 2) =3.447553595058008e-04;
pointCoord(983, 0) = 4.750879598867512e-01; pointCoord(983, 1) =2.987580805105672e-02; pointCoord(983, 2) =3.447553595058008e-04;
pointCoord(984, 0) = 7.250879598867512e-01; pointCoord(984, 1) =2.987580805105672e-02; pointCoord(984, 2) =3.447553595058008e-04;
pointCoord(985, 0) = 2.250879598867512e-01; pointCoord(985, 1) =2.798758080510567e-01; pointCoord(985, 2) =3.447553595058008e-04;
pointCoord(986, 0) = 4.750879598867512e-01; pointCoord(986, 1) =2.798758080510567e-01; pointCoord(986, 2) =3.447553595058008e-04;
pointCoord(987, 0) = 2.250879598867512e-01; pointCoord(987, 1) =5.298758080510567e-01; pointCoord(987, 2) =3.447553595058008e-04;
pointCoord(988, 0) = 1.918691247451343e-01; pointCoord(988, 1) =6.309464319267366e-02; pointCoord(988, 2) =4.863366507167933e-04;
pointCoord(989, 0) = 4.418691247451343e-01; pointCoord(989, 1) =6.309464319267366e-02; pointCoord(989, 2) =4.863366507167933e-04;
pointCoord(990, 0) = 6.918691247451343e-01; pointCoord(990, 1) =6.309464319267366e-02; pointCoord(990, 2) =4.863366507167933e-04;
pointCoord(991, 0) = 1.918691247451343e-01; pointCoord(991, 1) =3.130946431926737e-01; pointCoord(991, 2) =4.863366507167933e-04;
pointCoord(992, 0) = 4.418691247451343e-01; pointCoord(992, 1) =3.130946431926737e-01; pointCoord(992, 2) =4.863366507167933e-04;
pointCoord(993, 0) = 1.918691247451343e-01; pointCoord(993, 1) =5.630946431926737e-01; pointCoord(993, 2) =4.863366507167933e-04;
pointCoord(994, 0) = 1.499559507823086e-01; pointCoord(994, 1) =1.050078171554993e-01; pointCoord(994, 2) =5.622654757087529e-04;
pointCoord(995, 0) = 3.999559507823086e-01; pointCoord(995, 1) =1.050078171554993e-01; pointCoord(995, 2) =5.622654757087529e-04;
pointCoord(996, 0) = 6.499559507823086e-01; pointCoord(996, 1) =1.050078171554993e-01; pointCoord(996, 2) =5.622654757087529e-04;
pointCoord(997, 0) = 1.499559507823086e-01; pointCoord(997, 1) =3.550078171554993e-01; pointCoord(997, 2) =5.622654757087529e-04;
pointCoord(998, 0) = 3.999559507823086e-01; pointCoord(998, 1) =3.550078171554993e-01; pointCoord(998, 2) =5.622654757087529e-04;
pointCoord(999, 0) = 1.499559507823086e-01; pointCoord(999, 1) =6.050078171554993e-01; pointCoord(999, 2) =5.622654757087529e-04;
pointCoord(1000, 0) = 1.050078171554993e-01; pointCoord(1000, 1) =1.499559507823086e-01; pointCoord(1000, 2) =5.622654757087529e-04;
pointCoord(1001, 0) = 3.550078171554993e-01; pointCoord(1001, 1) =1.499559507823086e-01; pointCoord(1001, 2) =5.622654757087529e-04;
pointCoord(1002, 0) = 6.050078171554993e-01; pointCoord(1002, 1) =1.499559507823086e-01; pointCoord(1002, 2) =5.622654757087529e-04;
pointCoord(1003, 0) = 1.050078171554993e-01; pointCoord(1003, 1) =3.999559507823086e-01; pointCoord(1003, 2) =5.622654757087529e-04;
pointCoord(1004, 0) = 3.550078171554993e-01; pointCoord(1004, 1) =3.999559507823086e-01; pointCoord(1004, 2) =5.622654757087529e-04;
pointCoord(1005, 0) = 1.050078171554993e-01; pointCoord(1005, 1) =6.499559507823086e-01; pointCoord(1005, 2) =5.622654757087529e-04;
pointCoord(1006, 0) = 6.309464319267366e-02; pointCoord(1006, 1) =1.918691247451343e-01; pointCoord(1006, 2) =4.863366507167933e-04;
pointCoord(1007, 0) = 3.130946431926737e-01; pointCoord(1007, 1) =1.918691247451343e-01; pointCoord(1007, 2) =4.863366507167933e-04;
pointCoord(1008, 0) = 5.630946431926737e-01; pointCoord(1008, 1) =1.918691247451343e-01; pointCoord(1008, 2) =4.863366507167933e-04;
pointCoord(1009, 0) = 6.309464319267366e-02; pointCoord(1009, 1) =4.418691247451343e-01; pointCoord(1009, 2) =4.863366507167933e-04;
pointCoord(1010, 0) = 3.130946431926737e-01; pointCoord(1010, 1) =4.418691247451343e-01; pointCoord(1010, 2) =4.863366507167933e-04;
pointCoord(1011, 0) = 6.309464319267366e-02; pointCoord(1011, 1) =6.918691247451343e-01; pointCoord(1011, 2) =4.863366507167933e-04;
pointCoord(1012, 0) = 2.987580805105672e-02; pointCoord(1012, 1) =2.250879598867512e-01; pointCoord(1012, 2) =3.447553595058008e-04;
pointCoord(1013, 0) = 2.798758080510567e-01; pointCoord(1013, 1) =2.250879598867512e-01; pointCoord(1013, 2) =3.447553595058008e-04;
pointCoord(1014, 0) = 5.298758080510567e-01; pointCoord(1014, 1) =2.250879598867512e-01; pointCoord(1014, 2) =3.447553595058008e-04;
pointCoord(1015, 0) = 2.987580805105672e-02; pointCoord(1015, 1) =4.750879598867512e-01; pointCoord(1015, 2) =3.447553595058008e-04;
pointCoord(1016, 0) = 2.798758080510567e-01; pointCoord(1016, 1) =4.750879598867512e-01; pointCoord(1016, 2) =3.447553595058008e-04;
pointCoord(1017, 0) = 2.987580805105672e-02; pointCoord(1017, 1) =7.250879598867512e-01; pointCoord(1017, 2) =3.447553595058008e-04;
pointCoord(1018, 0) = 9.828979907054253e-03; pointCoord(1018, 1) =2.451347880307537e-01; pointCoord(1018, 2) =1.569337084289531e-04;
pointCoord(1019, 0) = 2.598289799070542e-01; pointCoord(1019, 1) =2.451347880307537e-01; pointCoord(1019, 2) =1.569337084289531e-04;
pointCoord(1020, 0) = 5.098289799070542e-01; pointCoord(1020, 1) =2.451347880307537e-01; pointCoord(1020, 2) =1.569337084289531e-04;
pointCoord(1021, 0) = 9.828979907054253e-03; pointCoord(1021, 1) =4.951347880307537e-01; pointCoord(1021, 2) =1.569337084289531e-04;
pointCoord(1022, 0) = 2.598289799070542e-01; pointCoord(1022, 1) =4.951347880307537e-01; pointCoord(1022, 2) =1.569337084289531e-04;
pointCoord(1023, 0) = 9.828979907054253e-03; pointCoord(1023, 1) =7.451347880307537e-01; pointCoord(1023, 2) =1.569337084289531e-04;
    }
    else if(nH == 1600)
    {
        pointCoord(0, 0) = 7.884477484931305e-05; pointCoord(0, 1) =3.892169575397059e-03; pointCoord(0, 2) =2.034592200391275e-06;
pointCoord(1, 0) = 2.000788447748493e-01; pointCoord(1, 1) =3.892169575397059e-03; pointCoord(1, 2) =2.034592200391275e-06;
pointCoord(2, 0) = 4.000788447748493e-01; pointCoord(2, 1) =3.892169575397059e-03; pointCoord(2, 2) =2.034592200391275e-06;
pointCoord(3, 0) = 6.000788447748493e-01; pointCoord(3, 1) =3.892169575397059e-03; pointCoord(3, 2) =2.034592200391275e-06;
pointCoord(4, 0) = 8.000788447748494e-01; pointCoord(4, 1) =3.892169575397059e-03; pointCoord(4, 2) =2.034592200391275e-06;
pointCoord(5, 0) = 7.884477484931305e-05; pointCoord(5, 1) =2.038921695753971e-01; pointCoord(5, 2) =2.034592200391275e-06;
pointCoord(6, 0) = 2.000788447748493e-01; pointCoord(6, 1) =2.038921695753971e-01; pointCoord(6, 2) =2.034592200391275e-06;
pointCoord(7, 0) = 4.000788447748493e-01; pointCoord(7, 1) =2.038921695753971e-01; pointCoord(7, 2) =2.034592200391275e-06;
pointCoord(8, 0) = 6.000788447748493e-01; pointCoord(8, 1) =2.038921695753971e-01; pointCoord(8, 2) =2.034592200391275e-06;
pointCoord(9, 0) = 7.884477484931305e-05; pointCoord(9, 1) =4.038921695753971e-01; pointCoord(9, 2) =2.034592200391275e-06;
pointCoord(10, 0) = 2.000788447748493e-01; pointCoord(10, 1) =4.038921695753971e-01; pointCoord(10, 2) =2.034592200391275e-06;
pointCoord(11, 0) = 4.000788447748493e-01; pointCoord(11, 1) =4.038921695753971e-01; pointCoord(11, 2) =2.034592200391275e-06;
pointCoord(12, 0) = 7.884477484931305e-05; pointCoord(12, 1) =6.038921695753970e-01; pointCoord(12, 2) =2.034592200391275e-06;
pointCoord(13, 0) = 2.000788447748493e-01; pointCoord(13, 1) =6.038921695753970e-01; pointCoord(13, 2) =2.034592200391275e-06;
pointCoord(14, 0) = 7.884477484931305e-05; pointCoord(14, 1) =8.038921695753971e-01; pointCoord(14, 2) =2.034592200391275e-06;
pointCoord(15, 0) = 4.037201680383164e-04; pointCoord(15, 1) =3.567294182208055e-03; pointCoord(15, 2) =4.469636080836934e-06;
pointCoord(16, 0) = 2.004037201680383e-01; pointCoord(16, 1) =3.567294182208055e-03; pointCoord(16, 2) =4.469636080836934e-06;
pointCoord(17, 0) = 4.004037201680383e-01; pointCoord(17, 1) =3.567294182208055e-03; pointCoord(17, 2) =4.469636080836934e-06;
pointCoord(18, 0) = 6.004037201680383e-01; pointCoord(18, 1) =3.567294182208055e-03; pointCoord(18, 2) =4.469636080836934e-06;
pointCoord(19, 0) = 8.004037201680384e-01; pointCoord(19, 1) =3.567294182208055e-03; pointCoord(19, 2) =4.469636080836934e-06;
pointCoord(20, 0) = 4.037201680383164e-04; pointCoord(20, 1) =2.035672941822081e-01; pointCoord(20, 2) =4.469636080836934e-06;
pointCoord(21, 0) = 2.004037201680383e-01; pointCoord(21, 1) =2.035672941822081e-01; pointCoord(21, 2) =4.469636080836934e-06;
pointCoord(22, 0) = 4.004037201680383e-01; pointCoord(22, 1) =2.035672941822081e-01; pointCoord(22, 2) =4.469636080836934e-06;
pointCoord(23, 0) = 6.004037201680383e-01; pointCoord(23, 1) =2.035672941822081e-01; pointCoord(23, 2) =4.469636080836934e-06;
pointCoord(24, 0) = 4.037201680383164e-04; pointCoord(24, 1) =4.035672941822081e-01; pointCoord(24, 2) =4.469636080836934e-06;
pointCoord(25, 0) = 2.004037201680383e-01; pointCoord(25, 1) =4.035672941822081e-01; pointCoord(25, 2) =4.469636080836934e-06;
pointCoord(26, 0) = 4.004037201680383e-01; pointCoord(26, 1) =4.035672941822081e-01; pointCoord(26, 2) =4.469636080836934e-06;
pointCoord(27, 0) = 4.037201680383164e-04; pointCoord(27, 1) =6.035672941822080e-01; pointCoord(27, 2) =4.469636080836934e-06;
pointCoord(28, 0) = 2.004037201680383e-01; pointCoord(28, 1) =6.035672941822080e-01; pointCoord(28, 2) =4.469636080836934e-06;
pointCoord(29, 0) = 4.037201680383164e-04; pointCoord(29, 1) =8.035672941822081e-01; pointCoord(29, 2) =4.469636080836934e-06;
pointCoord(30, 0) = 9.420588044745352e-04; pointCoord(30, 1) =3.028955545771836e-03; pointCoord(30, 2) =6.305189409073112e-06;
pointCoord(31, 0) = 2.009420588044745e-01; pointCoord(31, 1) =3.028955545771836e-03; pointCoord(31, 2) =6.305189409073112e-06;
pointCoord(32, 0) = 4.009420588044745e-01; pointCoord(32, 1) =3.028955545771836e-03; pointCoord(32, 2) =6.305189409073112e-06;
pointCoord(33, 0) = 6.009420588044745e-01; pointCoord(33, 1) =3.028955545771836e-03; pointCoord(33, 2) =6.305189409073112e-06;
pointCoord(34, 0) = 8.009420588044746e-01; pointCoord(34, 1) =3.028955545771836e-03; pointCoord(34, 2) =6.305189409073112e-06;
pointCoord(35, 0) = 9.420588044745352e-04; pointCoord(35, 1) =2.030289555457719e-01; pointCoord(35, 2) =6.305189409073112e-06;
pointCoord(36, 0) = 2.009420588044745e-01; pointCoord(36, 1) =2.030289555457719e-01; pointCoord(36, 2) =6.305189409073112e-06;
pointCoord(37, 0) = 4.009420588044745e-01; pointCoord(37, 1) =2.030289555457719e-01; pointCoord(37, 2) =6.305189409073112e-06;
pointCoord(38, 0) = 6.009420588044745e-01; pointCoord(38, 1) =2.030289555457719e-01; pointCoord(38, 2) =6.305189409073112e-06;
pointCoord(39, 0) = 9.420588044745352e-04; pointCoord(39, 1) =4.030289555457718e-01; pointCoord(39, 2) =6.305189409073112e-06;
pointCoord(40, 0) = 2.009420588044745e-01; pointCoord(40, 1) =4.030289555457718e-01; pointCoord(40, 2) =6.305189409073112e-06;
pointCoord(41, 0) = 4.009420588044745e-01; pointCoord(41, 1) =4.030289555457718e-01; pointCoord(41, 2) =6.305189409073112e-06;
pointCoord(42, 0) = 9.420588044745352e-04; pointCoord(42, 1) =6.030289555457718e-01; pointCoord(42, 2) =6.305189409073112e-06;
pointCoord(43, 0) = 2.009420588044745e-01; pointCoord(43, 1) =6.030289555457718e-01; pointCoord(43, 2) =6.305189409073112e-06;
pointCoord(44, 0) = 9.420588044745352e-04; pointCoord(44, 1) =8.030289555457719e-01; pointCoord(44, 2) =6.305189409073112e-06;
pointCoord(45, 0) = 1.621296376281917e-03; pointCoord(45, 1) =2.349717973964455e-03; pointCoord(45, 2) =7.289580822874776e-06;
pointCoord(46, 0) = 2.016212963762819e-01; pointCoord(46, 1) =2.349717973964455e-03; pointCoord(46, 2) =7.289580822874776e-06;
pointCoord(47, 0) = 4.016212963762820e-01; pointCoord(47, 1) =2.349717973964455e-03; pointCoord(47, 2) =7.289580822874776e-06;
pointCoord(48, 0) = 6.016212963762819e-01; pointCoord(48, 1) =2.349717973964455e-03; pointCoord(48, 2) =7.289580822874776e-06;
pointCoord(49, 0) = 8.016212963762820e-01; pointCoord(49, 1) =2.349717973964455e-03; pointCoord(49, 2) =7.289580822874776e-06;
pointCoord(50, 0) = 1.621296376281917e-03; pointCoord(50, 1) =2.023497179739645e-01; pointCoord(50, 2) =7.289580822874776e-06;
pointCoord(51, 0) = 2.016212963762819e-01; pointCoord(51, 1) =2.023497179739645e-01; pointCoord(51, 2) =7.289580822874776e-06;
pointCoord(52, 0) = 4.016212963762820e-01; pointCoord(52, 1) =2.023497179739645e-01; pointCoord(52, 2) =7.289580822874776e-06;
pointCoord(53, 0) = 6.016212963762819e-01; pointCoord(53, 1) =2.023497179739645e-01; pointCoord(53, 2) =7.289580822874776e-06;
pointCoord(54, 0) = 1.621296376281917e-03; pointCoord(54, 1) =4.023497179739645e-01; pointCoord(54, 2) =7.289580822874776e-06;
pointCoord(55, 0) = 2.016212963762819e-01; pointCoord(55, 1) =4.023497179739645e-01; pointCoord(55, 2) =7.289580822874776e-06;
pointCoord(56, 0) = 4.016212963762820e-01; pointCoord(56, 1) =4.023497179739645e-01; pointCoord(56, 2) =7.289580822874776e-06;
pointCoord(57, 0) = 1.621296376281917e-03; pointCoord(57, 1) =6.023497179739644e-01; pointCoord(57, 2) =7.289580822874776e-06;
pointCoord(58, 0) = 2.016212963762819e-01; pointCoord(58, 1) =6.023497179739644e-01; pointCoord(58, 2) =7.289580822874776e-06;
pointCoord(59, 0) = 1.621296376281917e-03; pointCoord(59, 1) =8.023497179739645e-01; pointCoord(59, 2) =7.289580822874776e-06;
pointCoord(60, 0) = 2.349717973964455e-03; pointCoord(60, 1) =1.621296376281917e-03; pointCoord(60, 2) =7.289580822874776e-06;
pointCoord(61, 0) = 2.023497179739645e-01; pointCoord(61, 1) =1.621296376281917e-03; pointCoord(61, 2) =7.289580822874776e-06;
pointCoord(62, 0) = 4.023497179739645e-01; pointCoord(62, 1) =1.621296376281917e-03; pointCoord(62, 2) =7.289580822874776e-06;
pointCoord(63, 0) = 6.023497179739644e-01; pointCoord(63, 1) =1.621296376281917e-03; pointCoord(63, 2) =7.289580822874776e-06;
pointCoord(64, 0) = 8.023497179739645e-01; pointCoord(64, 1) =1.621296376281917e-03; pointCoord(64, 2) =7.289580822874776e-06;
pointCoord(65, 0) = 2.349717973964455e-03; pointCoord(65, 1) =2.016212963762819e-01; pointCoord(65, 2) =7.289580822874776e-06;
pointCoord(66, 0) = 2.023497179739645e-01; pointCoord(66, 1) =2.016212963762819e-01; pointCoord(66, 2) =7.289580822874776e-06;
pointCoord(67, 0) = 4.023497179739645e-01; pointCoord(67, 1) =2.016212963762819e-01; pointCoord(67, 2) =7.289580822874776e-06;
pointCoord(68, 0) = 6.023497179739644e-01; pointCoord(68, 1) =2.016212963762819e-01; pointCoord(68, 2) =7.289580822874776e-06;
pointCoord(69, 0) = 2.349717973964455e-03; pointCoord(69, 1) =4.016212963762820e-01; pointCoord(69, 2) =7.289580822874776e-06;
pointCoord(70, 0) = 2.023497179739645e-01; pointCoord(70, 1) =4.016212963762820e-01; pointCoord(70, 2) =7.289580822874776e-06;
pointCoord(71, 0) = 4.023497179739645e-01; pointCoord(71, 1) =4.016212963762820e-01; pointCoord(71, 2) =7.289580822874776e-06;
pointCoord(72, 0) = 2.349717973964455e-03; pointCoord(72, 1) =6.016212963762819e-01; pointCoord(72, 2) =7.289580822874776e-06;
pointCoord(73, 0) = 2.023497179739645e-01; pointCoord(73, 1) =6.016212963762819e-01; pointCoord(73, 2) =7.289580822874776e-06;
pointCoord(74, 0) = 2.349717973964455e-03; pointCoord(74, 1) =8.016212963762820e-01; pointCoord(74, 2) =7.289580822874776e-06;
pointCoord(75, 0) = 3.028955545771836e-03; pointCoord(75, 1) =9.420588044745352e-04; pointCoord(75, 2) =6.305189409073112e-06;
pointCoord(76, 0) = 2.030289555457719e-01; pointCoord(76, 1) =9.420588044745352e-04; pointCoord(76, 2) =6.305189409073112e-06;
pointCoord(77, 0) = 4.030289555457718e-01; pointCoord(77, 1) =9.420588044745352e-04; pointCoord(77, 2) =6.305189409073112e-06;
pointCoord(78, 0) = 6.030289555457718e-01; pointCoord(78, 1) =9.420588044745352e-04; pointCoord(78, 2) =6.305189409073112e-06;
pointCoord(79, 0) = 8.030289555457719e-01; pointCoord(79, 1) =9.420588044745352e-04; pointCoord(79, 2) =6.305189409073112e-06;
pointCoord(80, 0) = 3.028955545771836e-03; pointCoord(80, 1) =2.009420588044745e-01; pointCoord(80, 2) =6.305189409073112e-06;
pointCoord(81, 0) = 2.030289555457719e-01; pointCoord(81, 1) =2.009420588044745e-01; pointCoord(81, 2) =6.305189409073112e-06;
pointCoord(82, 0) = 4.030289555457718e-01; pointCoord(82, 1) =2.009420588044745e-01; pointCoord(82, 2) =6.305189409073112e-06;
pointCoord(83, 0) = 6.030289555457718e-01; pointCoord(83, 1) =2.009420588044745e-01; pointCoord(83, 2) =6.305189409073112e-06;
pointCoord(84, 0) = 3.028955545771836e-03; pointCoord(84, 1) =4.009420588044745e-01; pointCoord(84, 2) =6.305189409073112e-06;
pointCoord(85, 0) = 2.030289555457719e-01; pointCoord(85, 1) =4.009420588044745e-01; pointCoord(85, 2) =6.305189409073112e-06;
pointCoord(86, 0) = 4.030289555457718e-01; pointCoord(86, 1) =4.009420588044745e-01; pointCoord(86, 2) =6.305189409073112e-06;
pointCoord(87, 0) = 3.028955545771836e-03; pointCoord(87, 1) =6.009420588044745e-01; pointCoord(87, 2) =6.305189409073112e-06;
pointCoord(88, 0) = 2.030289555457719e-01; pointCoord(88, 1) =6.009420588044745e-01; pointCoord(88, 2) =6.305189409073112e-06;
pointCoord(89, 0) = 3.028955545771836e-03; pointCoord(89, 1) =8.009420588044746e-01; pointCoord(89, 2) =6.305189409073112e-06;
pointCoord(90, 0) = 3.567294182208055e-03; pointCoord(90, 1) =4.037201680383164e-04; pointCoord(90, 2) =4.469636080836934e-06;
pointCoord(91, 0) = 2.035672941822081e-01; pointCoord(91, 1) =4.037201680383164e-04; pointCoord(91, 2) =4.469636080836934e-06;
pointCoord(92, 0) = 4.035672941822081e-01; pointCoord(92, 1) =4.037201680383164e-04; pointCoord(92, 2) =4.469636080836934e-06;
pointCoord(93, 0) = 6.035672941822080e-01; pointCoord(93, 1) =4.037201680383164e-04; pointCoord(93, 2) =4.469636080836934e-06;
pointCoord(94, 0) = 8.035672941822081e-01; pointCoord(94, 1) =4.037201680383164e-04; pointCoord(94, 2) =4.469636080836934e-06;
pointCoord(95, 0) = 3.567294182208055e-03; pointCoord(95, 1) =2.004037201680383e-01; pointCoord(95, 2) =4.469636080836934e-06;
pointCoord(96, 0) = 2.035672941822081e-01; pointCoord(96, 1) =2.004037201680383e-01; pointCoord(96, 2) =4.469636080836934e-06;
pointCoord(97, 0) = 4.035672941822081e-01; pointCoord(97, 1) =2.004037201680383e-01; pointCoord(97, 2) =4.469636080836934e-06;
pointCoord(98, 0) = 6.035672941822080e-01; pointCoord(98, 1) =2.004037201680383e-01; pointCoord(98, 2) =4.469636080836934e-06;
pointCoord(99, 0) = 3.567294182208055e-03; pointCoord(99, 1) =4.004037201680383e-01; pointCoord(99, 2) =4.469636080836934e-06;
pointCoord(100, 0) = 2.035672941822081e-01; pointCoord(100, 1) =4.004037201680383e-01; pointCoord(100, 2) =4.469636080836934e-06;
pointCoord(101, 0) = 4.035672941822081e-01; pointCoord(101, 1) =4.004037201680383e-01; pointCoord(101, 2) =4.469636080836934e-06;
pointCoord(102, 0) = 3.567294182208055e-03; pointCoord(102, 1) =6.004037201680383e-01; pointCoord(102, 2) =4.469636080836934e-06;
pointCoord(103, 0) = 2.035672941822081e-01; pointCoord(103, 1) =6.004037201680383e-01; pointCoord(103, 2) =4.469636080836934e-06;
pointCoord(104, 0) = 3.567294182208055e-03; pointCoord(104, 1) =8.004037201680384e-01; pointCoord(104, 2) =4.469636080836934e-06;
pointCoord(105, 0) = 3.892169575397059e-03; pointCoord(105, 1) =7.884477484931305e-05; pointCoord(105, 2) =2.034592200391275e-06;
pointCoord(106, 0) = 2.038921695753971e-01; pointCoord(106, 1) =7.884477484931305e-05; pointCoord(106, 2) =2.034592200391275e-06;
pointCoord(107, 0) = 4.038921695753971e-01; pointCoord(107, 1) =7.884477484931305e-05; pointCoord(107, 2) =2.034592200391275e-06;
pointCoord(108, 0) = 6.038921695753970e-01; pointCoord(108, 1) =7.884477484931305e-05; pointCoord(108, 2) =2.034592200391275e-06;
pointCoord(109, 0) = 8.038921695753971e-01; pointCoord(109, 1) =7.884477484931305e-05; pointCoord(109, 2) =2.034592200391275e-06;
pointCoord(110, 0) = 3.892169575397059e-03; pointCoord(110, 1) =2.000788447748493e-01; pointCoord(110, 2) =2.034592200391275e-06;
pointCoord(111, 0) = 2.038921695753971e-01; pointCoord(111, 1) =2.000788447748493e-01; pointCoord(111, 2) =2.034592200391275e-06;
pointCoord(112, 0) = 4.038921695753971e-01; pointCoord(112, 1) =2.000788447748493e-01; pointCoord(112, 2) =2.034592200391275e-06;
pointCoord(113, 0) = 6.038921695753970e-01; pointCoord(113, 1) =2.000788447748493e-01; pointCoord(113, 2) =2.034592200391275e-06;
pointCoord(114, 0) = 3.892169575397059e-03; pointCoord(114, 1) =4.000788447748493e-01; pointCoord(114, 2) =2.034592200391275e-06;
pointCoord(115, 0) = 2.038921695753971e-01; pointCoord(115, 1) =4.000788447748493e-01; pointCoord(115, 2) =2.034592200391275e-06;
pointCoord(116, 0) = 4.038921695753971e-01; pointCoord(116, 1) =4.000788447748493e-01; pointCoord(116, 2) =2.034592200391275e-06;
pointCoord(117, 0) = 3.892169575397059e-03; pointCoord(117, 1) =6.000788447748493e-01; pointCoord(117, 2) =2.034592200391275e-06;
pointCoord(118, 0) = 2.038921695753971e-01; pointCoord(118, 1) =6.000788447748493e-01; pointCoord(118, 2) =2.034592200391275e-06;
pointCoord(119, 0) = 3.892169575397059e-03; pointCoord(119, 1) =8.000788447748494e-01; pointCoord(119, 2) =2.034592200391275e-06;
pointCoord(120, 0) = 4.037201680383164e-04; pointCoord(120, 1) =1.992963209059901e-02; pointCoord(120, 2) =2.288651636172858e-05;
pointCoord(121, 0) = 2.004037201680383e-01; pointCoord(121, 1) =1.992963209059901e-02; pointCoord(121, 2) =2.288651636172858e-05;
pointCoord(122, 0) = 4.004037201680383e-01; pointCoord(122, 1) =1.992963209059901e-02; pointCoord(122, 2) =2.288651636172858e-05;
pointCoord(123, 0) = 6.004037201680383e-01; pointCoord(123, 1) =1.992963209059901e-02; pointCoord(123, 2) =2.288651636172858e-05;
pointCoord(124, 0) = 8.004037201680384e-01; pointCoord(124, 1) =1.992963209059901e-02; pointCoord(124, 2) =2.288651636172858e-05;
pointCoord(125, 0) = 4.037201680383164e-04; pointCoord(125, 1) =2.199296320905990e-01; pointCoord(125, 2) =2.288651636172858e-05;
pointCoord(126, 0) = 2.004037201680383e-01; pointCoord(126, 1) =2.199296320905990e-01; pointCoord(126, 2) =2.288651636172858e-05;
pointCoord(127, 0) = 4.004037201680383e-01; pointCoord(127, 1) =2.199296320905990e-01; pointCoord(127, 2) =2.288651636172858e-05;
pointCoord(128, 0) = 6.004037201680383e-01; pointCoord(128, 1) =2.199296320905990e-01; pointCoord(128, 2) =2.288651636172858e-05;
pointCoord(129, 0) = 4.037201680383164e-04; pointCoord(129, 1) =4.199296320905990e-01; pointCoord(129, 2) =2.288651636172858e-05;
pointCoord(130, 0) = 2.004037201680383e-01; pointCoord(130, 1) =4.199296320905990e-01; pointCoord(130, 2) =2.288651636172858e-05;
pointCoord(131, 0) = 4.004037201680383e-01; pointCoord(131, 1) =4.199296320905990e-01; pointCoord(131, 2) =2.288651636172858e-05;
pointCoord(132, 0) = 4.037201680383164e-04; pointCoord(132, 1) =6.199296320905990e-01; pointCoord(132, 2) =2.288651636172858e-05;
pointCoord(133, 0) = 2.004037201680383e-01; pointCoord(133, 1) =6.199296320905990e-01; pointCoord(133, 2) =2.288651636172858e-05;
pointCoord(134, 0) = 4.037201680383164e-04; pointCoord(134, 1) =8.199296320905991e-01; pointCoord(134, 2) =2.288651636172858e-05;
pointCoord(135, 0) = 2.067226070369159e-03; pointCoord(135, 1) =1.826612618826817e-02; pointCoord(135, 2) =5.027759335525545e-05;
pointCoord(136, 0) = 2.020672260703692e-01; pointCoord(136, 1) =1.826612618826817e-02; pointCoord(136, 2) =5.027759335525545e-05;
pointCoord(137, 0) = 4.020672260703692e-01; pointCoord(137, 1) =1.826612618826817e-02; pointCoord(137, 2) =5.027759335525545e-05;
pointCoord(138, 0) = 6.020672260703691e-01; pointCoord(138, 1) =1.826612618826817e-02; pointCoord(138, 2) =5.027759335525545e-05;
pointCoord(139, 0) = 8.020672260703692e-01; pointCoord(139, 1) =1.826612618826817e-02; pointCoord(139, 2) =5.027759335525545e-05;
pointCoord(140, 0) = 2.067226070369159e-03; pointCoord(140, 1) =2.182661261882682e-01; pointCoord(140, 2) =5.027759335525545e-05;
pointCoord(141, 0) = 2.020672260703692e-01; pointCoord(141, 1) =2.182661261882682e-01; pointCoord(141, 2) =5.027759335525545e-05;
pointCoord(142, 0) = 4.020672260703692e-01; pointCoord(142, 1) =2.182661261882682e-01; pointCoord(142, 2) =5.027759335525545e-05;
pointCoord(143, 0) = 6.020672260703691e-01; pointCoord(143, 1) =2.182661261882682e-01; pointCoord(143, 2) =5.027759335525545e-05;
pointCoord(144, 0) = 2.067226070369159e-03; pointCoord(144, 1) =4.182661261882682e-01; pointCoord(144, 2) =5.027759335525545e-05;
pointCoord(145, 0) = 2.020672260703692e-01; pointCoord(145, 1) =4.182661261882682e-01; pointCoord(145, 2) =5.027759335525545e-05;
pointCoord(146, 0) = 4.020672260703692e-01; pointCoord(146, 1) =4.182661261882682e-01; pointCoord(146, 2) =5.027759335525545e-05;
pointCoord(147, 0) = 2.067226070369159e-03; pointCoord(147, 1) =6.182661261882682e-01; pointCoord(147, 2) =5.027759335525545e-05;
pointCoord(148, 0) = 2.020672260703692e-01; pointCoord(148, 1) =6.182661261882682e-01; pointCoord(148, 2) =5.027759335525545e-05;
pointCoord(149, 0) = 2.067226070369159e-03; pointCoord(149, 1) =8.182661261882682e-01; pointCoord(149, 2) =5.027759335525545e-05;
pointCoord(150, 0) = 4.823758322239011e-03; pointCoord(150, 1) =1.550959393639832e-02; pointCoord(150, 2) =7.092518124604937e-05;
pointCoord(151, 0) = 2.048237583222390e-01; pointCoord(151, 1) =1.550959393639832e-02; pointCoord(151, 2) =7.092518124604937e-05;
pointCoord(152, 0) = 4.048237583222390e-01; pointCoord(152, 1) =1.550959393639832e-02; pointCoord(152, 2) =7.092518124604937e-05;
pointCoord(153, 0) = 6.048237583222390e-01; pointCoord(153, 1) =1.550959393639832e-02; pointCoord(153, 2) =7.092518124604937e-05;
pointCoord(154, 0) = 8.048237583222391e-01; pointCoord(154, 1) =1.550959393639832e-02; pointCoord(154, 2) =7.092518124604937e-05;
pointCoord(155, 0) = 4.823758322239011e-03; pointCoord(155, 1) =2.155095939363983e-01; pointCoord(155, 2) =7.092518124604937e-05;
pointCoord(156, 0) = 2.048237583222390e-01; pointCoord(156, 1) =2.155095939363983e-01; pointCoord(156, 2) =7.092518124604937e-05;
pointCoord(157, 0) = 4.048237583222390e-01; pointCoord(157, 1) =2.155095939363983e-01; pointCoord(157, 2) =7.092518124604937e-05;
pointCoord(158, 0) = 6.048237583222390e-01; pointCoord(158, 1) =2.155095939363983e-01; pointCoord(158, 2) =7.092518124604937e-05;
pointCoord(159, 0) = 4.823758322239011e-03; pointCoord(159, 1) =4.155095939363984e-01; pointCoord(159, 2) =7.092518124604937e-05;
pointCoord(160, 0) = 2.048237583222390e-01; pointCoord(160, 1) =4.155095939363984e-01; pointCoord(160, 2) =7.092518124604937e-05;
pointCoord(161, 0) = 4.048237583222390e-01; pointCoord(161, 1) =4.155095939363984e-01; pointCoord(161, 2) =7.092518124604937e-05;
pointCoord(162, 0) = 4.823758322239011e-03; pointCoord(162, 1) =6.155095939363983e-01; pointCoord(162, 2) =7.092518124604937e-05;
pointCoord(163, 0) = 2.048237583222390e-01; pointCoord(163, 1) =6.155095939363983e-01; pointCoord(163, 2) =7.092518124604937e-05;
pointCoord(164, 0) = 4.823758322239011e-03; pointCoord(164, 1) =8.155095939363983e-01; pointCoord(164, 2) =7.092518124604937e-05;
pointCoord(165, 0) = 8.301755528168038e-03; pointCoord(165, 1) =1.203159673046929e-02; pointCoord(165, 2) =8.199830449599807e-05;
pointCoord(166, 0) = 2.083017555281680e-01; pointCoord(166, 1) =1.203159673046929e-02; pointCoord(166, 2) =8.199830449599807e-05;
pointCoord(167, 0) = 4.083017555281681e-01; pointCoord(167, 1) =1.203159673046929e-02; pointCoord(167, 2) =8.199830449599807e-05;
pointCoord(168, 0) = 6.083017555281680e-01; pointCoord(168, 1) =1.203159673046929e-02; pointCoord(168, 2) =8.199830449599807e-05;
pointCoord(169, 0) = 8.083017555281681e-01; pointCoord(169, 1) =1.203159673046929e-02; pointCoord(169, 2) =8.199830449599807e-05;
pointCoord(170, 0) = 8.301755528168038e-03; pointCoord(170, 1) =2.120315967304693e-01; pointCoord(170, 2) =8.199830449599807e-05;
pointCoord(171, 0) = 2.083017555281680e-01; pointCoord(171, 1) =2.120315967304693e-01; pointCoord(171, 2) =8.199830449599807e-05;
pointCoord(172, 0) = 4.083017555281681e-01; pointCoord(172, 1) =2.120315967304693e-01; pointCoord(172, 2) =8.199830449599807e-05;
pointCoord(173, 0) = 6.083017555281680e-01; pointCoord(173, 1) =2.120315967304693e-01; pointCoord(173, 2) =8.199830449599807e-05;
pointCoord(174, 0) = 8.301755528168038e-03; pointCoord(174, 1) =4.120315967304693e-01; pointCoord(174, 2) =8.199830449599807e-05;
pointCoord(175, 0) = 2.083017555281680e-01; pointCoord(175, 1) =4.120315967304693e-01; pointCoord(175, 2) =8.199830449599807e-05;
pointCoord(176, 0) = 4.083017555281681e-01; pointCoord(176, 1) =4.120315967304693e-01; pointCoord(176, 2) =8.199830449599807e-05;
pointCoord(177, 0) = 8.301755528168038e-03; pointCoord(177, 1) =6.120315967304693e-01; pointCoord(177, 2) =8.199830449599807e-05;
pointCoord(178, 0) = 2.083017555281680e-01; pointCoord(178, 1) =6.120315967304693e-01; pointCoord(178, 2) =8.199830449599807e-05;
pointCoord(179, 0) = 8.301755528168038e-03; pointCoord(179, 1) =8.120315967304693e-01; pointCoord(179, 2) =8.199830449599807e-05;
pointCoord(180, 0) = 1.203159673046929e-02; pointCoord(180, 1) =8.301755528168038e-03; pointCoord(180, 2) =8.199830449599807e-05;
pointCoord(181, 0) = 2.120315967304693e-01; pointCoord(181, 1) =8.301755528168038e-03; pointCoord(181, 2) =8.199830449599807e-05;
pointCoord(182, 0) = 4.120315967304693e-01; pointCoord(182, 1) =8.301755528168038e-03; pointCoord(182, 2) =8.199830449599807e-05;
pointCoord(183, 0) = 6.120315967304693e-01; pointCoord(183, 1) =8.301755528168038e-03; pointCoord(183, 2) =8.199830449599807e-05;
pointCoord(184, 0) = 8.120315967304693e-01; pointCoord(184, 1) =8.301755528168038e-03; pointCoord(184, 2) =8.199830449599807e-05;
pointCoord(185, 0) = 1.203159673046929e-02; pointCoord(185, 1) =2.083017555281680e-01; pointCoord(185, 2) =8.199830449599807e-05;
pointCoord(186, 0) = 2.120315967304693e-01; pointCoord(186, 1) =2.083017555281680e-01; pointCoord(186, 2) =8.199830449599807e-05;
pointCoord(187, 0) = 4.120315967304693e-01; pointCoord(187, 1) =2.083017555281680e-01; pointCoord(187, 2) =8.199830449599807e-05;
pointCoord(188, 0) = 6.120315967304693e-01; pointCoord(188, 1) =2.083017555281680e-01; pointCoord(188, 2) =8.199830449599807e-05;
pointCoord(189, 0) = 1.203159673046929e-02; pointCoord(189, 1) =4.083017555281681e-01; pointCoord(189, 2) =8.199830449599807e-05;
pointCoord(190, 0) = 2.120315967304693e-01; pointCoord(190, 1) =4.083017555281681e-01; pointCoord(190, 2) =8.199830449599807e-05;
pointCoord(191, 0) = 4.120315967304693e-01; pointCoord(191, 1) =4.083017555281681e-01; pointCoord(191, 2) =8.199830449599807e-05;
pointCoord(192, 0) = 1.203159673046929e-02; pointCoord(192, 1) =6.083017555281680e-01; pointCoord(192, 2) =8.199830449599807e-05;
pointCoord(193, 0) = 2.120315967304693e-01; pointCoord(193, 1) =6.083017555281680e-01; pointCoord(193, 2) =8.199830449599807e-05;
pointCoord(194, 0) = 1.203159673046929e-02; pointCoord(194, 1) =8.083017555281681e-01; pointCoord(194, 2) =8.199830449599807e-05;
pointCoord(195, 0) = 1.550959393639832e-02; pointCoord(195, 1) =4.823758322239011e-03; pointCoord(195, 2) =7.092518124604937e-05;
pointCoord(196, 0) = 2.155095939363983e-01; pointCoord(196, 1) =4.823758322239011e-03; pointCoord(196, 2) =7.092518124604937e-05;
pointCoord(197, 0) = 4.155095939363984e-01; pointCoord(197, 1) =4.823758322239011e-03; pointCoord(197, 2) =7.092518124604937e-05;
pointCoord(198, 0) = 6.155095939363983e-01; pointCoord(198, 1) =4.823758322239011e-03; pointCoord(198, 2) =7.092518124604937e-05;
pointCoord(199, 0) = 8.155095939363983e-01; pointCoord(199, 1) =4.823758322239011e-03; pointCoord(199, 2) =7.092518124604937e-05;
pointCoord(200, 0) = 1.550959393639832e-02; pointCoord(200, 1) =2.048237583222390e-01; pointCoord(200, 2) =7.092518124604937e-05;
pointCoord(201, 0) = 2.155095939363983e-01; pointCoord(201, 1) =2.048237583222390e-01; pointCoord(201, 2) =7.092518124604937e-05;
pointCoord(202, 0) = 4.155095939363984e-01; pointCoord(202, 1) =2.048237583222390e-01; pointCoord(202, 2) =7.092518124604937e-05;
pointCoord(203, 0) = 6.155095939363983e-01; pointCoord(203, 1) =2.048237583222390e-01; pointCoord(203, 2) =7.092518124604937e-05;
pointCoord(204, 0) = 1.550959393639832e-02; pointCoord(204, 1) =4.048237583222390e-01; pointCoord(204, 2) =7.092518124604937e-05;
pointCoord(205, 0) = 2.155095939363983e-01; pointCoord(205, 1) =4.048237583222390e-01; pointCoord(205, 2) =7.092518124604937e-05;
pointCoord(206, 0) = 4.155095939363984e-01; pointCoord(206, 1) =4.048237583222390e-01; pointCoord(206, 2) =7.092518124604937e-05;
pointCoord(207, 0) = 1.550959393639832e-02; pointCoord(207, 1) =6.048237583222390e-01; pointCoord(207, 2) =7.092518124604937e-05;
pointCoord(208, 0) = 2.155095939363983e-01; pointCoord(208, 1) =6.048237583222390e-01; pointCoord(208, 2) =7.092518124604937e-05;
pointCoord(209, 0) = 1.550959393639832e-02; pointCoord(209, 1) =8.048237583222391e-01; pointCoord(209, 2) =7.092518124604937e-05;
pointCoord(210, 0) = 1.826612618826817e-02; pointCoord(210, 1) =2.067226070369159e-03; pointCoord(210, 2) =5.027759335525545e-05;
pointCoord(211, 0) = 2.182661261882682e-01; pointCoord(211, 1) =2.067226070369159e-03; pointCoord(211, 2) =5.027759335525545e-05;
pointCoord(212, 0) = 4.182661261882682e-01; pointCoord(212, 1) =2.067226070369159e-03; pointCoord(212, 2) =5.027759335525545e-05;
pointCoord(213, 0) = 6.182661261882682e-01; pointCoord(213, 1) =2.067226070369159e-03; pointCoord(213, 2) =5.027759335525545e-05;
pointCoord(214, 0) = 8.182661261882682e-01; pointCoord(214, 1) =2.067226070369159e-03; pointCoord(214, 2) =5.027759335525545e-05;
pointCoord(215, 0) = 1.826612618826817e-02; pointCoord(215, 1) =2.020672260703692e-01; pointCoord(215, 2) =5.027759335525545e-05;
pointCoord(216, 0) = 2.182661261882682e-01; pointCoord(216, 1) =2.020672260703692e-01; pointCoord(216, 2) =5.027759335525545e-05;
pointCoord(217, 0) = 4.182661261882682e-01; pointCoord(217, 1) =2.020672260703692e-01; pointCoord(217, 2) =5.027759335525545e-05;
pointCoord(218, 0) = 6.182661261882682e-01; pointCoord(218, 1) =2.020672260703692e-01; pointCoord(218, 2) =5.027759335525545e-05;
pointCoord(219, 0) = 1.826612618826817e-02; pointCoord(219, 1) =4.020672260703692e-01; pointCoord(219, 2) =5.027759335525545e-05;
pointCoord(220, 0) = 2.182661261882682e-01; pointCoord(220, 1) =4.020672260703692e-01; pointCoord(220, 2) =5.027759335525545e-05;
pointCoord(221, 0) = 4.182661261882682e-01; pointCoord(221, 1) =4.020672260703692e-01; pointCoord(221, 2) =5.027759335525545e-05;
pointCoord(222, 0) = 1.826612618826817e-02; pointCoord(222, 1) =6.020672260703691e-01; pointCoord(222, 2) =5.027759335525545e-05;
pointCoord(223, 0) = 2.182661261882682e-01; pointCoord(223, 1) =6.020672260703691e-01; pointCoord(223, 2) =5.027759335525545e-05;
pointCoord(224, 0) = 1.826612618826817e-02; pointCoord(224, 1) =8.020672260703692e-01; pointCoord(224, 2) =5.027759335525545e-05;
pointCoord(225, 0) = 1.992963209059901e-02; pointCoord(225, 1) =4.037201680383164e-04; pointCoord(225, 2) =2.288651636172858e-05;
pointCoord(226, 0) = 2.199296320905990e-01; pointCoord(226, 1) =4.037201680383164e-04; pointCoord(226, 2) =2.288651636172858e-05;
pointCoord(227, 0) = 4.199296320905990e-01; pointCoord(227, 1) =4.037201680383164e-04; pointCoord(227, 2) =2.288651636172858e-05;
pointCoord(228, 0) = 6.199296320905990e-01; pointCoord(228, 1) =4.037201680383164e-04; pointCoord(228, 2) =2.288651636172858e-05;
pointCoord(229, 0) = 8.199296320905991e-01; pointCoord(229, 1) =4.037201680383164e-04; pointCoord(229, 2) =2.288651636172858e-05;
pointCoord(230, 0) = 1.992963209059901e-02; pointCoord(230, 1) =2.004037201680383e-01; pointCoord(230, 2) =2.288651636172858e-05;
pointCoord(231, 0) = 2.199296320905990e-01; pointCoord(231, 1) =2.004037201680383e-01; pointCoord(231, 2) =2.288651636172858e-05;
pointCoord(232, 0) = 4.199296320905990e-01; pointCoord(232, 1) =2.004037201680383e-01; pointCoord(232, 2) =2.288651636172858e-05;
pointCoord(233, 0) = 6.199296320905990e-01; pointCoord(233, 1) =2.004037201680383e-01; pointCoord(233, 2) =2.288651636172858e-05;
pointCoord(234, 0) = 1.992963209059901e-02; pointCoord(234, 1) =4.004037201680383e-01; pointCoord(234, 2) =2.288651636172858e-05;
pointCoord(235, 0) = 2.199296320905990e-01; pointCoord(235, 1) =4.004037201680383e-01; pointCoord(235, 2) =2.288651636172858e-05;
pointCoord(236, 0) = 4.199296320905990e-01; pointCoord(236, 1) =4.004037201680383e-01; pointCoord(236, 2) =2.288651636172858e-05;
pointCoord(237, 0) = 1.992963209059901e-02; pointCoord(237, 1) =6.004037201680383e-01; pointCoord(237, 2) =2.288651636172858e-05;
pointCoord(238, 0) = 2.199296320905990e-01; pointCoord(238, 1) =6.004037201680383e-01; pointCoord(238, 2) =2.288651636172858e-05;
pointCoord(239, 0) = 1.992963209059901e-02; pointCoord(239, 1) =8.004037201680384e-01; pointCoord(239, 2) =2.288651636172858e-05;
pointCoord(240, 0) = 9.420588044745352e-04; pointCoord(240, 1) =4.650470020389257e-02; pointCoord(240, 2) =7.533611717515951e-05;
pointCoord(241, 0) = 2.009420588044745e-01; pointCoord(241, 1) =4.650470020389257e-02; pointCoord(241, 2) =7.533611717515951e-05;
pointCoord(242, 0) = 4.009420588044745e-01; pointCoord(242, 1) =4.650470020389257e-02; pointCoord(242, 2) =7.533611717515951e-05;
pointCoord(243, 0) = 6.009420588044745e-01; pointCoord(243, 1) =4.650470020389257e-02; pointCoord(243, 2) =7.533611717515951e-05;
pointCoord(244, 0) = 8.009420588044746e-01; pointCoord(244, 1) =4.650470020389257e-02; pointCoord(244, 2) =7.533611717515951e-05;
pointCoord(245, 0) = 9.420588044745352e-04; pointCoord(245, 1) =2.465047002038926e-01; pointCoord(245, 2) =7.533611717515951e-05;
pointCoord(246, 0) = 2.009420588044745e-01; pointCoord(246, 1) =2.465047002038926e-01; pointCoord(246, 2) =7.533611717515951e-05;
pointCoord(247, 0) = 4.009420588044745e-01; pointCoord(247, 1) =2.465047002038926e-01; pointCoord(247, 2) =7.533611717515951e-05;
pointCoord(248, 0) = 6.009420588044745e-01; pointCoord(248, 1) =2.465047002038926e-01; pointCoord(248, 2) =7.533611717515951e-05;
pointCoord(249, 0) = 9.420588044745352e-04; pointCoord(249, 1) =4.465047002038926e-01; pointCoord(249, 2) =7.533611717515951e-05;
pointCoord(250, 0) = 2.009420588044745e-01; pointCoord(250, 1) =4.465047002038926e-01; pointCoord(250, 2) =7.533611717515951e-05;
pointCoord(251, 0) = 4.009420588044745e-01; pointCoord(251, 1) =4.465047002038926e-01; pointCoord(251, 2) =7.533611717515951e-05;
pointCoord(252, 0) = 9.420588044745352e-04; pointCoord(252, 1) =6.465047002038925e-01; pointCoord(252, 2) =7.533611717515951e-05;
pointCoord(253, 0) = 2.009420588044745e-01; pointCoord(253, 1) =6.465047002038925e-01; pointCoord(253, 2) =7.533611717515951e-05;
pointCoord(254, 0) = 9.420588044745352e-04; pointCoord(254, 1) =8.465047002038926e-01; pointCoord(254, 2) =7.533611717515951e-05;
pointCoord(255, 0) = 4.823758322239011e-03; pointCoord(255, 1) =4.262300068612809e-02; pointCoord(255, 2) =1.655000090197417e-04;
pointCoord(256, 0) = 2.048237583222390e-01; pointCoord(256, 1) =4.262300068612809e-02; pointCoord(256, 2) =1.655000090197417e-04;
pointCoord(257, 0) = 4.048237583222390e-01; pointCoord(257, 1) =4.262300068612809e-02; pointCoord(257, 2) =1.655000090197417e-04;
pointCoord(258, 0) = 6.048237583222390e-01; pointCoord(258, 1) =4.262300068612809e-02; pointCoord(258, 2) =1.655000090197417e-04;
pointCoord(259, 0) = 8.048237583222391e-01; pointCoord(259, 1) =4.262300068612809e-02; pointCoord(259, 2) =1.655000090197417e-04;
pointCoord(260, 0) = 4.823758322239011e-03; pointCoord(260, 1) =2.426230006861281e-01; pointCoord(260, 2) =1.655000090197417e-04;
pointCoord(261, 0) = 2.048237583222390e-01; pointCoord(261, 1) =2.426230006861281e-01; pointCoord(261, 2) =1.655000090197417e-04;
pointCoord(262, 0) = 4.048237583222390e-01; pointCoord(262, 1) =2.426230006861281e-01; pointCoord(262, 2) =1.655000090197417e-04;
pointCoord(263, 0) = 6.048237583222390e-01; pointCoord(263, 1) =2.426230006861281e-01; pointCoord(263, 2) =1.655000090197417e-04;
pointCoord(264, 0) = 4.823758322239011e-03; pointCoord(264, 1) =4.426230006861281e-01; pointCoord(264, 2) =1.655000090197417e-04;
pointCoord(265, 0) = 2.048237583222390e-01; pointCoord(265, 1) =4.426230006861281e-01; pointCoord(265, 2) =1.655000090197417e-04;
pointCoord(266, 0) = 4.048237583222390e-01; pointCoord(266, 1) =4.426230006861281e-01; pointCoord(266, 2) =1.655000090197417e-04;
pointCoord(267, 0) = 4.823758322239011e-03; pointCoord(267, 1) =6.426230006861281e-01; pointCoord(267, 2) =1.655000090197417e-04;
pointCoord(268, 0) = 2.048237583222390e-01; pointCoord(268, 1) =6.426230006861281e-01; pointCoord(268, 2) =1.655000090197417e-04;
pointCoord(269, 0) = 4.823758322239011e-03; pointCoord(269, 1) =8.426230006861282e-01; pointCoord(269, 2) =1.655000090197417e-04;
pointCoord(270, 0) = 1.125597470199032e-02; pointCoord(270, 1) =3.619078430637677e-02; pointCoord(270, 2) =2.334661894615331e-04;
pointCoord(271, 0) = 2.112559747019903e-01; pointCoord(271, 1) =3.619078430637677e-02; pointCoord(271, 2) =2.334661894615331e-04;
pointCoord(272, 0) = 4.112559747019903e-01; pointCoord(272, 1) =3.619078430637677e-02; pointCoord(272, 2) =2.334661894615331e-04;
pointCoord(273, 0) = 6.112559747019903e-01; pointCoord(273, 1) =3.619078430637677e-02; pointCoord(273, 2) =2.334661894615331e-04;
pointCoord(274, 0) = 8.112559747019904e-01; pointCoord(274, 1) =3.619078430637677e-02; pointCoord(274, 2) =2.334661894615331e-04;
pointCoord(275, 0) = 1.125597470199032e-02; pointCoord(275, 1) =2.361907843063768e-01; pointCoord(275, 2) =2.334661894615331e-04;
pointCoord(276, 0) = 2.112559747019903e-01; pointCoord(276, 1) =2.361907843063768e-01; pointCoord(276, 2) =2.334661894615331e-04;
pointCoord(277, 0) = 4.112559747019903e-01; pointCoord(277, 1) =2.361907843063768e-01; pointCoord(277, 2) =2.334661894615331e-04;
pointCoord(278, 0) = 6.112559747019903e-01; pointCoord(278, 1) =2.361907843063768e-01; pointCoord(278, 2) =2.334661894615331e-04;
pointCoord(279, 0) = 1.125597470199032e-02; pointCoord(279, 1) =4.361907843063768e-01; pointCoord(279, 2) =2.334661894615331e-04;
pointCoord(280, 0) = 2.112559747019903e-01; pointCoord(280, 1) =4.361907843063768e-01; pointCoord(280, 2) =2.334661894615331e-04;
pointCoord(281, 0) = 4.112559747019903e-01; pointCoord(281, 1) =4.361907843063768e-01; pointCoord(281, 2) =2.334661894615331e-04;
pointCoord(282, 0) = 1.125597470199032e-02; pointCoord(282, 1) =6.361907843063768e-01; pointCoord(282, 2) =2.334661894615331e-04;
pointCoord(283, 0) = 2.112559747019903e-01; pointCoord(283, 1) =6.361907843063768e-01; pointCoord(283, 2) =2.334661894615331e-04;
pointCoord(284, 0) = 1.125597470199032e-02; pointCoord(284, 1) =8.361907843063768e-01; pointCoord(284, 2) =2.334661894615331e-04;
pointCoord(285, 0) = 1.937168986604502e-02; pointCoord(285, 1) =2.807506914232209e-02; pointCoord(285, 2) =2.699158656581297e-04;
pointCoord(286, 0) = 2.193716898660450e-01; pointCoord(286, 1) =2.807506914232209e-02; pointCoord(286, 2) =2.699158656581297e-04;
pointCoord(287, 0) = 4.193716898660451e-01; pointCoord(287, 1) =2.807506914232209e-02; pointCoord(287, 2) =2.699158656581297e-04;
pointCoord(288, 0) = 6.193716898660450e-01; pointCoord(288, 1) =2.807506914232209e-02; pointCoord(288, 2) =2.699158656581297e-04;
pointCoord(289, 0) = 8.193716898660450e-01; pointCoord(289, 1) =2.807506914232209e-02; pointCoord(289, 2) =2.699158656581297e-04;
pointCoord(290, 0) = 1.937168986604502e-02; pointCoord(290, 1) =2.280750691423221e-01; pointCoord(290, 2) =2.699158656581297e-04;
pointCoord(291, 0) = 2.193716898660450e-01; pointCoord(291, 1) =2.280750691423221e-01; pointCoord(291, 2) =2.699158656581297e-04;
pointCoord(292, 0) = 4.193716898660451e-01; pointCoord(292, 1) =2.280750691423221e-01; pointCoord(292, 2) =2.699158656581297e-04;
pointCoord(293, 0) = 6.193716898660450e-01; pointCoord(293, 1) =2.280750691423221e-01; pointCoord(293, 2) =2.699158656581297e-04;
pointCoord(294, 0) = 1.937168986604502e-02; pointCoord(294, 1) =4.280750691423221e-01; pointCoord(294, 2) =2.699158656581297e-04;
pointCoord(295, 0) = 2.193716898660450e-01; pointCoord(295, 1) =4.280750691423221e-01; pointCoord(295, 2) =2.699158656581297e-04;
pointCoord(296, 0) = 4.193716898660451e-01; pointCoord(296, 1) =4.280750691423221e-01; pointCoord(296, 2) =2.699158656581297e-04;
pointCoord(297, 0) = 1.937168986604502e-02; pointCoord(297, 1) =6.280750691423220e-01; pointCoord(297, 2) =2.699158656581297e-04;
pointCoord(298, 0) = 2.193716898660450e-01; pointCoord(298, 1) =6.280750691423220e-01; pointCoord(298, 2) =2.699158656581297e-04;
pointCoord(299, 0) = 1.937168986604502e-02; pointCoord(299, 1) =8.280750691423221e-01; pointCoord(299, 2) =2.699158656581297e-04;
pointCoord(300, 0) = 2.807506914232209e-02; pointCoord(300, 1) =1.937168986604502e-02; pointCoord(300, 2) =2.699158656581297e-04;
pointCoord(301, 0) = 2.280750691423221e-01; pointCoord(301, 1) =1.937168986604502e-02; pointCoord(301, 2) =2.699158656581297e-04;
pointCoord(302, 0) = 4.280750691423221e-01; pointCoord(302, 1) =1.937168986604502e-02; pointCoord(302, 2) =2.699158656581297e-04;
pointCoord(303, 0) = 6.280750691423220e-01; pointCoord(303, 1) =1.937168986604502e-02; pointCoord(303, 2) =2.699158656581297e-04;
pointCoord(304, 0) = 8.280750691423221e-01; pointCoord(304, 1) =1.937168986604502e-02; pointCoord(304, 2) =2.699158656581297e-04;
pointCoord(305, 0) = 2.807506914232209e-02; pointCoord(305, 1) =2.193716898660450e-01; pointCoord(305, 2) =2.699158656581297e-04;
pointCoord(306, 0) = 2.280750691423221e-01; pointCoord(306, 1) =2.193716898660450e-01; pointCoord(306, 2) =2.699158656581297e-04;
pointCoord(307, 0) = 4.280750691423221e-01; pointCoord(307, 1) =2.193716898660450e-01; pointCoord(307, 2) =2.699158656581297e-04;
pointCoord(308, 0) = 6.280750691423220e-01; pointCoord(308, 1) =2.193716898660450e-01; pointCoord(308, 2) =2.699158656581297e-04;
pointCoord(309, 0) = 2.807506914232209e-02; pointCoord(309, 1) =4.193716898660451e-01; pointCoord(309, 2) =2.699158656581297e-04;
pointCoord(310, 0) = 2.280750691423221e-01; pointCoord(310, 1) =4.193716898660451e-01; pointCoord(310, 2) =2.699158656581297e-04;
pointCoord(311, 0) = 4.280750691423221e-01; pointCoord(311, 1) =4.193716898660451e-01; pointCoord(311, 2) =2.699158656581297e-04;
pointCoord(312, 0) = 2.807506914232209e-02; pointCoord(312, 1) =6.193716898660450e-01; pointCoord(312, 2) =2.699158656581297e-04;
pointCoord(313, 0) = 2.280750691423221e-01; pointCoord(313, 1) =6.193716898660450e-01; pointCoord(313, 2) =2.699158656581297e-04;
pointCoord(314, 0) = 2.807506914232209e-02; pointCoord(314, 1) =8.193716898660450e-01; pointCoord(314, 2) =2.699158656581297e-04;
pointCoord(315, 0) = 3.619078430637677e-02; pointCoord(315, 1) =1.125597470199032e-02; pointCoord(315, 2) =2.334661894615331e-04;
pointCoord(316, 0) = 2.361907843063768e-01; pointCoord(316, 1) =1.125597470199032e-02; pointCoord(316, 2) =2.334661894615331e-04;
pointCoord(317, 0) = 4.361907843063768e-01; pointCoord(317, 1) =1.125597470199032e-02; pointCoord(317, 2) =2.334661894615331e-04;
pointCoord(318, 0) = 6.361907843063768e-01; pointCoord(318, 1) =1.125597470199032e-02; pointCoord(318, 2) =2.334661894615331e-04;
pointCoord(319, 0) = 8.361907843063768e-01; pointCoord(319, 1) =1.125597470199032e-02; pointCoord(319, 2) =2.334661894615331e-04;
pointCoord(320, 0) = 3.619078430637677e-02; pointCoord(320, 1) =2.112559747019903e-01; pointCoord(320, 2) =2.334661894615331e-04;
pointCoord(321, 0) = 2.361907843063768e-01; pointCoord(321, 1) =2.112559747019903e-01; pointCoord(321, 2) =2.334661894615331e-04;
pointCoord(322, 0) = 4.361907843063768e-01; pointCoord(322, 1) =2.112559747019903e-01; pointCoord(322, 2) =2.334661894615331e-04;
pointCoord(323, 0) = 6.361907843063768e-01; pointCoord(323, 1) =2.112559747019903e-01; pointCoord(323, 2) =2.334661894615331e-04;
pointCoord(324, 0) = 3.619078430637677e-02; pointCoord(324, 1) =4.112559747019903e-01; pointCoord(324, 2) =2.334661894615331e-04;
pointCoord(325, 0) = 2.361907843063768e-01; pointCoord(325, 1) =4.112559747019903e-01; pointCoord(325, 2) =2.334661894615331e-04;
pointCoord(326, 0) = 4.361907843063768e-01; pointCoord(326, 1) =4.112559747019903e-01; pointCoord(326, 2) =2.334661894615331e-04;
pointCoord(327, 0) = 3.619078430637677e-02; pointCoord(327, 1) =6.112559747019903e-01; pointCoord(327, 2) =2.334661894615331e-04;
pointCoord(328, 0) = 2.361907843063768e-01; pointCoord(328, 1) =6.112559747019903e-01; pointCoord(328, 2) =2.334661894615331e-04;
pointCoord(329, 0) = 3.619078430637677e-02; pointCoord(329, 1) =8.112559747019904e-01; pointCoord(329, 2) =2.334661894615331e-04;
pointCoord(330, 0) = 4.262300068612809e-02; pointCoord(330, 1) =4.823758322239011e-03; pointCoord(330, 2) =1.655000090197417e-04;
pointCoord(331, 0) = 2.426230006861281e-01; pointCoord(331, 1) =4.823758322239011e-03; pointCoord(331, 2) =1.655000090197417e-04;
pointCoord(332, 0) = 4.426230006861281e-01; pointCoord(332, 1) =4.823758322239011e-03; pointCoord(332, 2) =1.655000090197417e-04;
pointCoord(333, 0) = 6.426230006861281e-01; pointCoord(333, 1) =4.823758322239011e-03; pointCoord(333, 2) =1.655000090197417e-04;
pointCoord(334, 0) = 8.426230006861282e-01; pointCoord(334, 1) =4.823758322239011e-03; pointCoord(334, 2) =1.655000090197417e-04;
pointCoord(335, 0) = 4.262300068612809e-02; pointCoord(335, 1) =2.048237583222390e-01; pointCoord(335, 2) =1.655000090197417e-04;
pointCoord(336, 0) = 2.426230006861281e-01; pointCoord(336, 1) =2.048237583222390e-01; pointCoord(336, 2) =1.655000090197417e-04;
pointCoord(337, 0) = 4.426230006861281e-01; pointCoord(337, 1) =2.048237583222390e-01; pointCoord(337, 2) =1.655000090197417e-04;
pointCoord(338, 0) = 6.426230006861281e-01; pointCoord(338, 1) =2.048237583222390e-01; pointCoord(338, 2) =1.655000090197417e-04;
pointCoord(339, 0) = 4.262300068612809e-02; pointCoord(339, 1) =4.048237583222390e-01; pointCoord(339, 2) =1.655000090197417e-04;
pointCoord(340, 0) = 2.426230006861281e-01; pointCoord(340, 1) =4.048237583222390e-01; pointCoord(340, 2) =1.655000090197417e-04;
pointCoord(341, 0) = 4.426230006861281e-01; pointCoord(341, 1) =4.048237583222390e-01; pointCoord(341, 2) =1.655000090197417e-04;
pointCoord(342, 0) = 4.262300068612809e-02; pointCoord(342, 1) =6.048237583222390e-01; pointCoord(342, 2) =1.655000090197417e-04;
pointCoord(343, 0) = 2.426230006861281e-01; pointCoord(343, 1) =6.048237583222390e-01; pointCoord(343, 2) =1.655000090197417e-04;
pointCoord(344, 0) = 4.262300068612809e-02; pointCoord(344, 1) =8.048237583222391e-01; pointCoord(344, 2) =1.655000090197417e-04;
pointCoord(345, 0) = 4.650470020389257e-02; pointCoord(345, 1) =9.420588044745352e-04; pointCoord(345, 2) =7.533611717515951e-05;
pointCoord(346, 0) = 2.465047002038926e-01; pointCoord(346, 1) =9.420588044745352e-04; pointCoord(346, 2) =7.533611717515951e-05;
pointCoord(347, 0) = 4.465047002038926e-01; pointCoord(347, 1) =9.420588044745352e-04; pointCoord(347, 2) =7.533611717515951e-05;
pointCoord(348, 0) = 6.465047002038925e-01; pointCoord(348, 1) =9.420588044745352e-04; pointCoord(348, 2) =7.533611717515951e-05;
pointCoord(349, 0) = 8.465047002038926e-01; pointCoord(349, 1) =9.420588044745352e-04; pointCoord(349, 2) =7.533611717515951e-05;
pointCoord(350, 0) = 4.650470020389257e-02; pointCoord(350, 1) =2.009420588044745e-01; pointCoord(350, 2) =7.533611717515951e-05;
pointCoord(351, 0) = 2.465047002038926e-01; pointCoord(351, 1) =2.009420588044745e-01; pointCoord(351, 2) =7.533611717515951e-05;
pointCoord(352, 0) = 4.465047002038926e-01; pointCoord(352, 1) =2.009420588044745e-01; pointCoord(352, 2) =7.533611717515951e-05;
pointCoord(353, 0) = 6.465047002038925e-01; pointCoord(353, 1) =2.009420588044745e-01; pointCoord(353, 2) =7.533611717515951e-05;
pointCoord(354, 0) = 4.650470020389257e-02; pointCoord(354, 1) =4.009420588044745e-01; pointCoord(354, 2) =7.533611717515951e-05;
pointCoord(355, 0) = 2.465047002038926e-01; pointCoord(355, 1) =4.009420588044745e-01; pointCoord(355, 2) =7.533611717515951e-05;
pointCoord(356, 0) = 4.465047002038926e-01; pointCoord(356, 1) =4.009420588044745e-01; pointCoord(356, 2) =7.533611717515951e-05;
pointCoord(357, 0) = 4.650470020389257e-02; pointCoord(357, 1) =6.009420588044745e-01; pointCoord(357, 2) =7.533611717515951e-05;
pointCoord(358, 0) = 2.465047002038926e-01; pointCoord(358, 1) =6.009420588044745e-01; pointCoord(358, 2) =7.533611717515951e-05;
pointCoord(359, 0) = 4.650470020389257e-02; pointCoord(359, 1) =8.009420588044746e-01; pointCoord(359, 2) =7.533611717515951e-05;
pointCoord(360, 0) = 1.621296376281917e-03; pointCoord(360, 1) =8.003523937415311e-02; pointCoord(360, 2) =1.498966925243746e-04;
pointCoord(361, 0) = 2.016212963762819e-01; pointCoord(361, 1) =8.003523937415311e-02; pointCoord(361, 2) =1.498966925243746e-04;
pointCoord(362, 0) = 4.016212963762820e-01; pointCoord(362, 1) =8.003523937415311e-02; pointCoord(362, 2) =1.498966925243746e-04;
pointCoord(363, 0) = 6.016212963762819e-01; pointCoord(363, 1) =8.003523937415311e-02; pointCoord(363, 2) =1.498966925243746e-04;
pointCoord(364, 0) = 8.016212963762820e-01; pointCoord(364, 1) =8.003523937415311e-02; pointCoord(364, 2) =1.498966925243746e-04;
pointCoord(365, 0) = 1.621296376281917e-03; pointCoord(365, 1) =2.800352393741531e-01; pointCoord(365, 2) =1.498966925243746e-04;
pointCoord(366, 0) = 2.016212963762819e-01; pointCoord(366, 1) =2.800352393741531e-01; pointCoord(366, 2) =1.498966925243746e-04;
pointCoord(367, 0) = 4.016212963762820e-01; pointCoord(367, 1) =2.800352393741531e-01; pointCoord(367, 2) =1.498966925243746e-04;
pointCoord(368, 0) = 6.016212963762819e-01; pointCoord(368, 1) =2.800352393741531e-01; pointCoord(368, 2) =1.498966925243746e-04;
pointCoord(369, 0) = 1.621296376281917e-03; pointCoord(369, 1) =4.800352393741532e-01; pointCoord(369, 2) =1.498966925243746e-04;
pointCoord(370, 0) = 2.016212963762819e-01; pointCoord(370, 1) =4.800352393741532e-01; pointCoord(370, 2) =1.498966925243746e-04;
pointCoord(371, 0) = 4.016212963762820e-01; pointCoord(371, 1) =4.800352393741532e-01; pointCoord(371, 2) =1.498966925243746e-04;
pointCoord(372, 0) = 1.621296376281917e-03; pointCoord(372, 1) =6.800352393741531e-01; pointCoord(372, 2) =1.498966925243746e-04;
pointCoord(373, 0) = 2.016212963762819e-01; pointCoord(373, 1) =6.800352393741531e-01; pointCoord(373, 2) =1.498966925243746e-04;
pointCoord(374, 0) = 1.621296376281917e-03; pointCoord(374, 1) =8.800352393741532e-01; pointCoord(374, 2) =1.498966925243746e-04;
pointCoord(375, 0) = 8.301755528168038e-03; pointCoord(375, 1) =7.335478022226698e-02; pointCoord(375, 2) =3.292962910091858e-04;
pointCoord(376, 0) = 2.083017555281680e-01; pointCoord(376, 1) =7.335478022226698e-02; pointCoord(376, 2) =3.292962910091858e-04;
pointCoord(377, 0) = 4.083017555281681e-01; pointCoord(377, 1) =7.335478022226698e-02; pointCoord(377, 2) =3.292962910091858e-04;
pointCoord(378, 0) = 6.083017555281680e-01; pointCoord(378, 1) =7.335478022226698e-02; pointCoord(378, 2) =3.292962910091858e-04;
pointCoord(379, 0) = 8.083017555281681e-01; pointCoord(379, 1) =7.335478022226698e-02; pointCoord(379, 2) =3.292962910091858e-04;
pointCoord(380, 0) = 8.301755528168038e-03; pointCoord(380, 1) =2.733547802222670e-01; pointCoord(380, 2) =3.292962910091858e-04;
pointCoord(381, 0) = 2.083017555281680e-01; pointCoord(381, 1) =2.733547802222670e-01; pointCoord(381, 2) =3.292962910091858e-04;
pointCoord(382, 0) = 4.083017555281681e-01; pointCoord(382, 1) =2.733547802222670e-01; pointCoord(382, 2) =3.292962910091858e-04;
pointCoord(383, 0) = 6.083017555281680e-01; pointCoord(383, 1) =2.733547802222670e-01; pointCoord(383, 2) =3.292962910091858e-04;
pointCoord(384, 0) = 8.301755528168038e-03; pointCoord(384, 1) =4.733547802222670e-01; pointCoord(384, 2) =3.292962910091858e-04;
pointCoord(385, 0) = 2.083017555281680e-01; pointCoord(385, 1) =4.733547802222670e-01; pointCoord(385, 2) =3.292962910091858e-04;
pointCoord(386, 0) = 4.083017555281681e-01; pointCoord(386, 1) =4.733547802222670e-01; pointCoord(386, 2) =3.292962910091858e-04;
pointCoord(387, 0) = 8.301755528168038e-03; pointCoord(387, 1) =6.733547802222669e-01; pointCoord(387, 2) =3.292962910091858e-04;
pointCoord(388, 0) = 2.083017555281680e-01; pointCoord(388, 1) =6.733547802222669e-01; pointCoord(388, 2) =3.292962910091858e-04;
pointCoord(389, 0) = 8.301755528168038e-03; pointCoord(389, 1) =8.733547802222670e-01; pointCoord(389, 2) =3.292962910091858e-04;
pointCoord(390, 0) = 1.937168986604502e-02; pointCoord(390, 1) =6.228484588439000e-02; pointCoord(390, 2) =4.645289793099657e-04;
pointCoord(391, 0) = 2.193716898660450e-01; pointCoord(391, 1) =6.228484588439000e-02; pointCoord(391, 2) =4.645289793099657e-04;
pointCoord(392, 0) = 4.193716898660451e-01; pointCoord(392, 1) =6.228484588439000e-02; pointCoord(392, 2) =4.645289793099657e-04;
pointCoord(393, 0) = 6.193716898660450e-01; pointCoord(393, 1) =6.228484588439000e-02; pointCoord(393, 2) =4.645289793099657e-04;
pointCoord(394, 0) = 8.193716898660450e-01; pointCoord(394, 1) =6.228484588439000e-02; pointCoord(394, 2) =4.645289793099657e-04;
pointCoord(395, 0) = 1.937168986604502e-02; pointCoord(395, 1) =2.622848458843900e-01; pointCoord(395, 2) =4.645289793099657e-04;
pointCoord(396, 0) = 2.193716898660450e-01; pointCoord(396, 1) =2.622848458843900e-01; pointCoord(396, 2) =4.645289793099657e-04;
pointCoord(397, 0) = 4.193716898660451e-01; pointCoord(397, 1) =2.622848458843900e-01; pointCoord(397, 2) =4.645289793099657e-04;
pointCoord(398, 0) = 6.193716898660450e-01; pointCoord(398, 1) =2.622848458843900e-01; pointCoord(398, 2) =4.645289793099657e-04;
pointCoord(399, 0) = 1.937168986604502e-02; pointCoord(399, 1) =4.622848458843900e-01; pointCoord(399, 2) =4.645289793099657e-04;
pointCoord(400, 0) = 2.193716898660450e-01; pointCoord(400, 1) =4.622848458843900e-01; pointCoord(400, 2) =4.645289793099657e-04;
pointCoord(401, 0) = 4.193716898660451e-01; pointCoord(401, 1) =4.622848458843900e-01; pointCoord(401, 2) =4.645289793099657e-04;
pointCoord(402, 0) = 1.937168986604502e-02; pointCoord(402, 1) =6.622848458843900e-01; pointCoord(402, 2) =4.645289793099657e-04;
pointCoord(403, 0) = 2.193716898660450e-01; pointCoord(403, 1) =6.622848458843900e-01; pointCoord(403, 2) =4.645289793099657e-04;
pointCoord(404, 0) = 1.937168986604502e-02; pointCoord(404, 1) =8.622848458843900e-01; pointCoord(404, 2) =4.645289793099657e-04;
pointCoord(405, 0) = 3.333894915381037e-02; pointCoord(405, 1) =4.831758659662465e-02; pointCoord(405, 2) =5.370531033333869e-04;
pointCoord(406, 0) = 2.333389491538104e-01; pointCoord(406, 1) =4.831758659662465e-02; pointCoord(406, 2) =5.370531033333869e-04;
pointCoord(407, 0) = 4.333389491538104e-01; pointCoord(407, 1) =4.831758659662465e-02; pointCoord(407, 2) =5.370531033333869e-04;
pointCoord(408, 0) = 6.333389491538104e-01; pointCoord(408, 1) =4.831758659662465e-02; pointCoord(408, 2) =5.370531033333869e-04;
pointCoord(409, 0) = 8.333389491538105e-01; pointCoord(409, 1) =4.831758659662465e-02; pointCoord(409, 2) =5.370531033333869e-04;
pointCoord(410, 0) = 3.333894915381037e-02; pointCoord(410, 1) =2.483175865966247e-01; pointCoord(410, 2) =5.370531033333869e-04;
pointCoord(411, 0) = 2.333389491538104e-01; pointCoord(411, 1) =2.483175865966247e-01; pointCoord(411, 2) =5.370531033333869e-04;
pointCoord(412, 0) = 4.333389491538104e-01; pointCoord(412, 1) =2.483175865966247e-01; pointCoord(412, 2) =5.370531033333869e-04;
pointCoord(413, 0) = 6.333389491538104e-01; pointCoord(413, 1) =2.483175865966247e-01; pointCoord(413, 2) =5.370531033333869e-04;
pointCoord(414, 0) = 3.333894915381037e-02; pointCoord(414, 1) =4.483175865966247e-01; pointCoord(414, 2) =5.370531033333869e-04;
pointCoord(415, 0) = 2.333389491538104e-01; pointCoord(415, 1) =4.483175865966247e-01; pointCoord(415, 2) =5.370531033333869e-04;
pointCoord(416, 0) = 4.333389491538104e-01; pointCoord(416, 1) =4.483175865966247e-01; pointCoord(416, 2) =5.370531033333869e-04;
pointCoord(417, 0) = 3.333894915381037e-02; pointCoord(417, 1) =6.483175865966246e-01; pointCoord(417, 2) =5.370531033333869e-04;
pointCoord(418, 0) = 2.333389491538104e-01; pointCoord(418, 1) =6.483175865966246e-01; pointCoord(418, 2) =5.370531033333869e-04;
pointCoord(419, 0) = 3.333894915381037e-02; pointCoord(419, 1) =8.483175865966247e-01; pointCoord(419, 2) =5.370531033333869e-04;
pointCoord(420, 0) = 4.831758659662465e-02; pointCoord(420, 1) =3.333894915381037e-02; pointCoord(420, 2) =5.370531033333869e-04;
pointCoord(421, 0) = 2.483175865966247e-01; pointCoord(421, 1) =3.333894915381037e-02; pointCoord(421, 2) =5.370531033333869e-04;
pointCoord(422, 0) = 4.483175865966247e-01; pointCoord(422, 1) =3.333894915381037e-02; pointCoord(422, 2) =5.370531033333869e-04;
pointCoord(423, 0) = 6.483175865966246e-01; pointCoord(423, 1) =3.333894915381037e-02; pointCoord(423, 2) =5.370531033333869e-04;
pointCoord(424, 0) = 8.483175865966247e-01; pointCoord(424, 1) =3.333894915381037e-02; pointCoord(424, 2) =5.370531033333869e-04;
pointCoord(425, 0) = 4.831758659662465e-02; pointCoord(425, 1) =2.333389491538104e-01; pointCoord(425, 2) =5.370531033333869e-04;
pointCoord(426, 0) = 2.483175865966247e-01; pointCoord(426, 1) =2.333389491538104e-01; pointCoord(426, 2) =5.370531033333869e-04;
pointCoord(427, 0) = 4.483175865966247e-01; pointCoord(427, 1) =2.333389491538104e-01; pointCoord(427, 2) =5.370531033333869e-04;
pointCoord(428, 0) = 6.483175865966246e-01; pointCoord(428, 1) =2.333389491538104e-01; pointCoord(428, 2) =5.370531033333869e-04;
pointCoord(429, 0) = 4.831758659662465e-02; pointCoord(429, 1) =4.333389491538104e-01; pointCoord(429, 2) =5.370531033333869e-04;
pointCoord(430, 0) = 2.483175865966247e-01; pointCoord(430, 1) =4.333389491538104e-01; pointCoord(430, 2) =5.370531033333869e-04;
pointCoord(431, 0) = 4.483175865966247e-01; pointCoord(431, 1) =4.333389491538104e-01; pointCoord(431, 2) =5.370531033333869e-04;
pointCoord(432, 0) = 4.831758659662465e-02; pointCoord(432, 1) =6.333389491538104e-01; pointCoord(432, 2) =5.370531033333869e-04;
pointCoord(433, 0) = 2.483175865966247e-01; pointCoord(433, 1) =6.333389491538104e-01; pointCoord(433, 2) =5.370531033333869e-04;
pointCoord(434, 0) = 4.831758659662465e-02; pointCoord(434, 1) =8.333389491538105e-01; pointCoord(434, 2) =5.370531033333869e-04;
pointCoord(435, 0) = 6.228484588439000e-02; pointCoord(435, 1) =1.937168986604502e-02; pointCoord(435, 2) =4.645289793099657e-04;
pointCoord(436, 0) = 2.622848458843900e-01; pointCoord(436, 1) =1.937168986604502e-02; pointCoord(436, 2) =4.645289793099657e-04;
pointCoord(437, 0) = 4.622848458843900e-01; pointCoord(437, 1) =1.937168986604502e-02; pointCoord(437, 2) =4.645289793099657e-04;
pointCoord(438, 0) = 6.622848458843900e-01; pointCoord(438, 1) =1.937168986604502e-02; pointCoord(438, 2) =4.645289793099657e-04;
pointCoord(439, 0) = 8.622848458843900e-01; pointCoord(439, 1) =1.937168986604502e-02; pointCoord(439, 2) =4.645289793099657e-04;
pointCoord(440, 0) = 6.228484588439000e-02; pointCoord(440, 1) =2.193716898660450e-01; pointCoord(440, 2) =4.645289793099657e-04;
pointCoord(441, 0) = 2.622848458843900e-01; pointCoord(441, 1) =2.193716898660450e-01; pointCoord(441, 2) =4.645289793099657e-04;
pointCoord(442, 0) = 4.622848458843900e-01; pointCoord(442, 1) =2.193716898660450e-01; pointCoord(442, 2) =4.645289793099657e-04;
pointCoord(443, 0) = 6.622848458843900e-01; pointCoord(443, 1) =2.193716898660450e-01; pointCoord(443, 2) =4.645289793099657e-04;
pointCoord(444, 0) = 6.228484588439000e-02; pointCoord(444, 1) =4.193716898660451e-01; pointCoord(444, 2) =4.645289793099657e-04;
pointCoord(445, 0) = 2.622848458843900e-01; pointCoord(445, 1) =4.193716898660451e-01; pointCoord(445, 2) =4.645289793099657e-04;
pointCoord(446, 0) = 4.622848458843900e-01; pointCoord(446, 1) =4.193716898660451e-01; pointCoord(446, 2) =4.645289793099657e-04;
pointCoord(447, 0) = 6.228484588439000e-02; pointCoord(447, 1) =6.193716898660450e-01; pointCoord(447, 2) =4.645289793099657e-04;
pointCoord(448, 0) = 2.622848458843900e-01; pointCoord(448, 1) =6.193716898660450e-01; pointCoord(448, 2) =4.645289793099657e-04;
pointCoord(449, 0) = 6.228484588439000e-02; pointCoord(449, 1) =8.193716898660450e-01; pointCoord(449, 2) =4.645289793099657e-04;
pointCoord(450, 0) = 7.335478022226698e-02; pointCoord(450, 1) =8.301755528168038e-03; pointCoord(450, 2) =3.292962910091858e-04;
pointCoord(451, 0) = 2.733547802222670e-01; pointCoord(451, 1) =8.301755528168038e-03; pointCoord(451, 2) =3.292962910091858e-04;
pointCoord(452, 0) = 4.733547802222670e-01; pointCoord(452, 1) =8.301755528168038e-03; pointCoord(452, 2) =3.292962910091858e-04;
pointCoord(453, 0) = 6.733547802222669e-01; pointCoord(453, 1) =8.301755528168038e-03; pointCoord(453, 2) =3.292962910091858e-04;
pointCoord(454, 0) = 8.733547802222670e-01; pointCoord(454, 1) =8.301755528168038e-03; pointCoord(454, 2) =3.292962910091858e-04;
pointCoord(455, 0) = 7.335478022226698e-02; pointCoord(455, 1) =2.083017555281680e-01; pointCoord(455, 2) =3.292962910091858e-04;
pointCoord(456, 0) = 2.733547802222670e-01; pointCoord(456, 1) =2.083017555281680e-01; pointCoord(456, 2) =3.292962910091858e-04;
pointCoord(457, 0) = 4.733547802222670e-01; pointCoord(457, 1) =2.083017555281680e-01; pointCoord(457, 2) =3.292962910091858e-04;
pointCoord(458, 0) = 6.733547802222669e-01; pointCoord(458, 1) =2.083017555281680e-01; pointCoord(458, 2) =3.292962910091858e-04;
pointCoord(459, 0) = 7.335478022226698e-02; pointCoord(459, 1) =4.083017555281681e-01; pointCoord(459, 2) =3.292962910091858e-04;
pointCoord(460, 0) = 2.733547802222670e-01; pointCoord(460, 1) =4.083017555281681e-01; pointCoord(460, 2) =3.292962910091858e-04;
pointCoord(461, 0) = 4.733547802222670e-01; pointCoord(461, 1) =4.083017555281681e-01; pointCoord(461, 2) =3.292962910091858e-04;
pointCoord(462, 0) = 7.335478022226698e-02; pointCoord(462, 1) =6.083017555281680e-01; pointCoord(462, 2) =3.292962910091858e-04;
pointCoord(463, 0) = 2.733547802222670e-01; pointCoord(463, 1) =6.083017555281680e-01; pointCoord(463, 2) =3.292962910091858e-04;
pointCoord(464, 0) = 7.335478022226698e-02; pointCoord(464, 1) =8.083017555281681e-01; pointCoord(464, 2) =3.292962910091858e-04;
pointCoord(465, 0) = 8.003523937415311e-02; pointCoord(465, 1) =1.621296376281917e-03; pointCoord(465, 2) =1.498966925243746e-04;
pointCoord(466, 0) = 2.800352393741531e-01; pointCoord(466, 1) =1.621296376281917e-03; pointCoord(466, 2) =1.498966925243746e-04;
pointCoord(467, 0) = 4.800352393741532e-01; pointCoord(467, 1) =1.621296376281917e-03; pointCoord(467, 2) =1.498966925243746e-04;
pointCoord(468, 0) = 6.800352393741531e-01; pointCoord(468, 1) =1.621296376281917e-03; pointCoord(468, 2) =1.498966925243746e-04;
pointCoord(469, 0) = 8.800352393741532e-01; pointCoord(469, 1) =1.621296376281917e-03; pointCoord(469, 2) =1.498966925243746e-04;
pointCoord(470, 0) = 8.003523937415311e-02; pointCoord(470, 1) =2.016212963762819e-01; pointCoord(470, 2) =1.498966925243746e-04;
pointCoord(471, 0) = 2.800352393741531e-01; pointCoord(471, 1) =2.016212963762819e-01; pointCoord(471, 2) =1.498966925243746e-04;
pointCoord(472, 0) = 4.800352393741532e-01; pointCoord(472, 1) =2.016212963762819e-01; pointCoord(472, 2) =1.498966925243746e-04;
pointCoord(473, 0) = 6.800352393741531e-01; pointCoord(473, 1) =2.016212963762819e-01; pointCoord(473, 2) =1.498966925243746e-04;
pointCoord(474, 0) = 8.003523937415311e-02; pointCoord(474, 1) =4.016212963762820e-01; pointCoord(474, 2) =1.498966925243746e-04;
pointCoord(475, 0) = 2.800352393741531e-01; pointCoord(475, 1) =4.016212963762820e-01; pointCoord(475, 2) =1.498966925243746e-04;
pointCoord(476, 0) = 4.800352393741532e-01; pointCoord(476, 1) =4.016212963762820e-01; pointCoord(476, 2) =1.498966925243746e-04;
pointCoord(477, 0) = 8.003523937415311e-02; pointCoord(477, 1) =6.016212963762819e-01; pointCoord(477, 2) =1.498966925243746e-04;
pointCoord(478, 0) = 2.800352393741531e-01; pointCoord(478, 1) =6.016212963762819e-01; pointCoord(478, 2) =1.498966925243746e-04;
pointCoord(479, 0) = 8.003523937415311e-02; pointCoord(479, 1) =8.016212963762820e-01; pointCoord(479, 2) =1.498966925243746e-04;
pointCoord(480, 0) = 2.349717973964455e-03; pointCoord(480, 1) =1.159937462756005e-01; pointCoord(480, 2) =2.172427927521020e-04;
pointCoord(481, 0) = 2.023497179739645e-01; pointCoord(481, 1) =1.159937462756005e-01; pointCoord(481, 2) =2.172427927521020e-04;
pointCoord(482, 0) = 4.023497179739645e-01; pointCoord(482, 1) =1.159937462756005e-01; pointCoord(482, 2) =2.172427927521020e-04;
pointCoord(483, 0) = 6.023497179739644e-01; pointCoord(483, 1) =1.159937462756005e-01; pointCoord(483, 2) =2.172427927521020e-04;
pointCoord(484, 0) = 8.023497179739645e-01; pointCoord(484, 1) =1.159937462756005e-01; pointCoord(484, 2) =2.172427927521020e-04;
pointCoord(485, 0) = 2.349717973964455e-03; pointCoord(485, 1) =3.159937462756005e-01; pointCoord(485, 2) =2.172427927521020e-04;
pointCoord(486, 0) = 2.023497179739645e-01; pointCoord(486, 1) =3.159937462756005e-01; pointCoord(486, 2) =2.172427927521020e-04;
pointCoord(487, 0) = 4.023497179739645e-01; pointCoord(487, 1) =3.159937462756005e-01; pointCoord(487, 2) =2.172427927521020e-04;
pointCoord(488, 0) = 6.023497179739644e-01; pointCoord(488, 1) =3.159937462756005e-01; pointCoord(488, 2) =2.172427927521020e-04;
pointCoord(489, 0) = 2.349717973964455e-03; pointCoord(489, 1) =5.159937462756006e-01; pointCoord(489, 2) =2.172427927521020e-04;
pointCoord(490, 0) = 2.023497179739645e-01; pointCoord(490, 1) =5.159937462756006e-01; pointCoord(490, 2) =2.172427927521020e-04;
pointCoord(491, 0) = 4.023497179739645e-01; pointCoord(491, 1) =5.159937462756006e-01; pointCoord(491, 2) =2.172427927521020e-04;
pointCoord(492, 0) = 2.349717973964455e-03; pointCoord(492, 1) =7.159937462756005e-01; pointCoord(492, 2) =2.172427927521020e-04;
pointCoord(493, 0) = 2.023497179739645e-01; pointCoord(493, 1) =7.159937462756005e-01; pointCoord(493, 2) =2.172427927521020e-04;
pointCoord(494, 0) = 2.349717973964455e-03; pointCoord(494, 1) =9.159937462756006e-01; pointCoord(494, 2) =2.172427927521020e-04;
pointCoord(495, 0) = 1.203159673046929e-02; pointCoord(495, 1) =1.063118675190957e-01; pointCoord(495, 2) =4.772436582622514e-04;
pointCoord(496, 0) = 2.120315967304693e-01; pointCoord(496, 1) =1.063118675190957e-01; pointCoord(496, 2) =4.772436582622514e-04;
pointCoord(497, 0) = 4.120315967304693e-01; pointCoord(497, 1) =1.063118675190957e-01; pointCoord(497, 2) =4.772436582622514e-04;
pointCoord(498, 0) = 6.120315967304693e-01; pointCoord(498, 1) =1.063118675190957e-01; pointCoord(498, 2) =4.772436582622514e-04;
pointCoord(499, 0) = 8.120315967304693e-01; pointCoord(499, 1) =1.063118675190957e-01; pointCoord(499, 2) =4.772436582622514e-04;
pointCoord(500, 0) = 1.203159673046929e-02; pointCoord(500, 1) =3.063118675190957e-01; pointCoord(500, 2) =4.772436582622514e-04;
pointCoord(501, 0) = 2.120315967304693e-01; pointCoord(501, 1) =3.063118675190957e-01; pointCoord(501, 2) =4.772436582622514e-04;
pointCoord(502, 0) = 4.120315967304693e-01; pointCoord(502, 1) =3.063118675190957e-01; pointCoord(502, 2) =4.772436582622514e-04;
pointCoord(503, 0) = 6.120315967304693e-01; pointCoord(503, 1) =3.063118675190957e-01; pointCoord(503, 2) =4.772436582622514e-04;
pointCoord(504, 0) = 1.203159673046929e-02; pointCoord(504, 1) =5.063118675190957e-01; pointCoord(504, 2) =4.772436582622514e-04;
pointCoord(505, 0) = 2.120315967304693e-01; pointCoord(505, 1) =5.063118675190957e-01; pointCoord(505, 2) =4.772436582622514e-04;
pointCoord(506, 0) = 4.120315967304693e-01; pointCoord(506, 1) =5.063118675190957e-01; pointCoord(506, 2) =4.772436582622514e-04;
pointCoord(507, 0) = 1.203159673046929e-02; pointCoord(507, 1) =7.063118675190957e-01; pointCoord(507, 2) =4.772436582622514e-04;
pointCoord(508, 0) = 2.120315967304693e-01; pointCoord(508, 1) =7.063118675190957e-01; pointCoord(508, 2) =4.772436582622514e-04;
pointCoord(509, 0) = 1.203159673046929e-02; pointCoord(509, 1) =9.063118675190958e-01; pointCoord(509, 2) =4.772436582622514e-04;
pointCoord(510, 0) = 2.807506914232209e-02; pointCoord(510, 1) =9.026839510724288e-02; pointCoord(510, 2) =6.732341526693159e-04;
pointCoord(511, 0) = 2.280750691423221e-01; pointCoord(511, 1) =9.026839510724288e-02; pointCoord(511, 2) =6.732341526693159e-04;
pointCoord(512, 0) = 4.280750691423221e-01; pointCoord(512, 1) =9.026839510724288e-02; pointCoord(512, 2) =6.732341526693159e-04;
pointCoord(513, 0) = 6.280750691423220e-01; pointCoord(513, 1) =9.026839510724288e-02; pointCoord(513, 2) =6.732341526693159e-04;
pointCoord(514, 0) = 8.280750691423221e-01; pointCoord(514, 1) =9.026839510724288e-02; pointCoord(514, 2) =6.732341526693159e-04;
pointCoord(515, 0) = 2.807506914232209e-02; pointCoord(515, 1) =2.902683951072429e-01; pointCoord(515, 2) =6.732341526693159e-04;
pointCoord(516, 0) = 2.280750691423221e-01; pointCoord(516, 1) =2.902683951072429e-01; pointCoord(516, 2) =6.732341526693159e-04;
pointCoord(517, 0) = 4.280750691423221e-01; pointCoord(517, 1) =2.902683951072429e-01; pointCoord(517, 2) =6.732341526693159e-04;
pointCoord(518, 0) = 6.280750691423220e-01; pointCoord(518, 1) =2.902683951072429e-01; pointCoord(518, 2) =6.732341526693159e-04;
pointCoord(519, 0) = 2.807506914232209e-02; pointCoord(519, 1) =4.902683951072429e-01; pointCoord(519, 2) =6.732341526693159e-04;
pointCoord(520, 0) = 2.280750691423221e-01; pointCoord(520, 1) =4.902683951072429e-01; pointCoord(520, 2) =6.732341526693159e-04;
pointCoord(521, 0) = 4.280750691423221e-01; pointCoord(521, 1) =4.902683951072429e-01; pointCoord(521, 2) =6.732341526693159e-04;
pointCoord(522, 0) = 2.807506914232209e-02; pointCoord(522, 1) =6.902683951072428e-01; pointCoord(522, 2) =6.732341526693159e-04;
pointCoord(523, 0) = 2.280750691423221e-01; pointCoord(523, 1) =6.902683951072428e-01; pointCoord(523, 2) =6.732341526693159e-04;
pointCoord(524, 0) = 2.807506914232209e-02; pointCoord(524, 1) =8.902683951072430e-01; pointCoord(524, 2) =6.732341526693159e-04;
pointCoord(525, 0) = 4.831758659662465e-02; pointCoord(525, 1) =7.002587765294030e-02; pointCoord(525, 2) =7.783421639230390e-04;
pointCoord(526, 0) = 2.483175865966247e-01; pointCoord(526, 1) =7.002587765294030e-02; pointCoord(526, 2) =7.783421639230390e-04;
pointCoord(527, 0) = 4.483175865966247e-01; pointCoord(527, 1) =7.002587765294030e-02; pointCoord(527, 2) =7.783421639230390e-04;
pointCoord(528, 0) = 6.483175865966246e-01; pointCoord(528, 1) =7.002587765294030e-02; pointCoord(528, 2) =7.783421639230390e-04;
pointCoord(529, 0) = 8.483175865966247e-01; pointCoord(529, 1) =7.002587765294030e-02; pointCoord(529, 2) =7.783421639230390e-04;
pointCoord(530, 0) = 4.831758659662465e-02; pointCoord(530, 1) =2.700258776529403e-01; pointCoord(530, 2) =7.783421639230390e-04;
pointCoord(531, 0) = 2.483175865966247e-01; pointCoord(531, 1) =2.700258776529403e-01; pointCoord(531, 2) =7.783421639230390e-04;
pointCoord(532, 0) = 4.483175865966247e-01; pointCoord(532, 1) =2.700258776529403e-01; pointCoord(532, 2) =7.783421639230390e-04;
pointCoord(533, 0) = 6.483175865966246e-01; pointCoord(533, 1) =2.700258776529403e-01; pointCoord(533, 2) =7.783421639230390e-04;
pointCoord(534, 0) = 4.831758659662465e-02; pointCoord(534, 1) =4.700258776529403e-01; pointCoord(534, 2) =7.783421639230390e-04;
pointCoord(535, 0) = 2.483175865966247e-01; pointCoord(535, 1) =4.700258776529403e-01; pointCoord(535, 2) =7.783421639230390e-04;
pointCoord(536, 0) = 4.483175865966247e-01; pointCoord(536, 1) =4.700258776529403e-01; pointCoord(536, 2) =7.783421639230390e-04;
pointCoord(537, 0) = 4.831758659662465e-02; pointCoord(537, 1) =6.700258776529403e-01; pointCoord(537, 2) =7.783421639230390e-04;
pointCoord(538, 0) = 2.483175865966247e-01; pointCoord(538, 1) =6.700258776529403e-01; pointCoord(538, 2) =7.783421639230390e-04;
pointCoord(539, 0) = 4.831758659662465e-02; pointCoord(539, 1) =8.700258776529404e-01; pointCoord(539, 2) =7.783421639230390e-04;
pointCoord(540, 0) = 7.002587765294030e-02; pointCoord(540, 1) =4.831758659662465e-02; pointCoord(540, 2) =7.783421639230390e-04;
pointCoord(541, 0) = 2.700258776529403e-01; pointCoord(541, 1) =4.831758659662465e-02; pointCoord(541, 2) =7.783421639230390e-04;
pointCoord(542, 0) = 4.700258776529403e-01; pointCoord(542, 1) =4.831758659662465e-02; pointCoord(542, 2) =7.783421639230390e-04;
pointCoord(543, 0) = 6.700258776529403e-01; pointCoord(543, 1) =4.831758659662465e-02; pointCoord(543, 2) =7.783421639230390e-04;
pointCoord(544, 0) = 8.700258776529404e-01; pointCoord(544, 1) =4.831758659662465e-02; pointCoord(544, 2) =7.783421639230390e-04;
pointCoord(545, 0) = 7.002587765294030e-02; pointCoord(545, 1) =2.483175865966247e-01; pointCoord(545, 2) =7.783421639230390e-04;
pointCoord(546, 0) = 2.700258776529403e-01; pointCoord(546, 1) =2.483175865966247e-01; pointCoord(546, 2) =7.783421639230390e-04;
pointCoord(547, 0) = 4.700258776529403e-01; pointCoord(547, 1) =2.483175865966247e-01; pointCoord(547, 2) =7.783421639230390e-04;
pointCoord(548, 0) = 6.700258776529403e-01; pointCoord(548, 1) =2.483175865966247e-01; pointCoord(548, 2) =7.783421639230390e-04;
pointCoord(549, 0) = 7.002587765294030e-02; pointCoord(549, 1) =4.483175865966247e-01; pointCoord(549, 2) =7.783421639230390e-04;
pointCoord(550, 0) = 2.700258776529403e-01; pointCoord(550, 1) =4.483175865966247e-01; pointCoord(550, 2) =7.783421639230390e-04;
pointCoord(551, 0) = 4.700258776529403e-01; pointCoord(551, 1) =4.483175865966247e-01; pointCoord(551, 2) =7.783421639230390e-04;
pointCoord(552, 0) = 7.002587765294030e-02; pointCoord(552, 1) =6.483175865966246e-01; pointCoord(552, 2) =7.783421639230390e-04;
pointCoord(553, 0) = 2.700258776529403e-01; pointCoord(553, 1) =6.483175865966246e-01; pointCoord(553, 2) =7.783421639230390e-04;
pointCoord(554, 0) = 7.002587765294030e-02; pointCoord(554, 1) =8.483175865966247e-01; pointCoord(554, 2) =7.783421639230390e-04;
pointCoord(555, 0) = 9.026839510724288e-02; pointCoord(555, 1) =2.807506914232209e-02; pointCoord(555, 2) =6.732341526693159e-04;
pointCoord(556, 0) = 2.902683951072429e-01; pointCoord(556, 1) =2.807506914232209e-02; pointCoord(556, 2) =6.732341526693159e-04;
pointCoord(557, 0) = 4.902683951072429e-01; pointCoord(557, 1) =2.807506914232209e-02; pointCoord(557, 2) =6.732341526693159e-04;
pointCoord(558, 0) = 6.902683951072428e-01; pointCoord(558, 1) =2.807506914232209e-02; pointCoord(558, 2) =6.732341526693159e-04;
pointCoord(559, 0) = 8.902683951072430e-01; pointCoord(559, 1) =2.807506914232209e-02; pointCoord(559, 2) =6.732341526693159e-04;
pointCoord(560, 0) = 9.026839510724288e-02; pointCoord(560, 1) =2.280750691423221e-01; pointCoord(560, 2) =6.732341526693159e-04;
pointCoord(561, 0) = 2.902683951072429e-01; pointCoord(561, 1) =2.280750691423221e-01; pointCoord(561, 2) =6.732341526693159e-04;
pointCoord(562, 0) = 4.902683951072429e-01; pointCoord(562, 1) =2.280750691423221e-01; pointCoord(562, 2) =6.732341526693159e-04;
pointCoord(563, 0) = 6.902683951072428e-01; pointCoord(563, 1) =2.280750691423221e-01; pointCoord(563, 2) =6.732341526693159e-04;
pointCoord(564, 0) = 9.026839510724288e-02; pointCoord(564, 1) =4.280750691423221e-01; pointCoord(564, 2) =6.732341526693159e-04;
pointCoord(565, 0) = 2.902683951072429e-01; pointCoord(565, 1) =4.280750691423221e-01; pointCoord(565, 2) =6.732341526693159e-04;
pointCoord(566, 0) = 4.902683951072429e-01; pointCoord(566, 1) =4.280750691423221e-01; pointCoord(566, 2) =6.732341526693159e-04;
pointCoord(567, 0) = 9.026839510724288e-02; pointCoord(567, 1) =6.280750691423220e-01; pointCoord(567, 2) =6.732341526693159e-04;
pointCoord(568, 0) = 2.902683951072429e-01; pointCoord(568, 1) =6.280750691423220e-01; pointCoord(568, 2) =6.732341526693159e-04;
pointCoord(569, 0) = 9.026839510724288e-02; pointCoord(569, 1) =8.280750691423221e-01; pointCoord(569, 2) =6.732341526693159e-04;
pointCoord(570, 0) = 1.063118675190957e-01; pointCoord(570, 1) =1.203159673046929e-02; pointCoord(570, 2) =4.772436582622514e-04;
pointCoord(571, 0) = 3.063118675190957e-01; pointCoord(571, 1) =1.203159673046929e-02; pointCoord(571, 2) =4.772436582622514e-04;
pointCoord(572, 0) = 5.063118675190957e-01; pointCoord(572, 1) =1.203159673046929e-02; pointCoord(572, 2) =4.772436582622514e-04;
pointCoord(573, 0) = 7.063118675190957e-01; pointCoord(573, 1) =1.203159673046929e-02; pointCoord(573, 2) =4.772436582622514e-04;
pointCoord(574, 0) = 9.063118675190958e-01; pointCoord(574, 1) =1.203159673046929e-02; pointCoord(574, 2) =4.772436582622514e-04;
pointCoord(575, 0) = 1.063118675190957e-01; pointCoord(575, 1) =2.120315967304693e-01; pointCoord(575, 2) =4.772436582622514e-04;
pointCoord(576, 0) = 3.063118675190957e-01; pointCoord(576, 1) =2.120315967304693e-01; pointCoord(576, 2) =4.772436582622514e-04;
pointCoord(577, 0) = 5.063118675190957e-01; pointCoord(577, 1) =2.120315967304693e-01; pointCoord(577, 2) =4.772436582622514e-04;
pointCoord(578, 0) = 7.063118675190957e-01; pointCoord(578, 1) =2.120315967304693e-01; pointCoord(578, 2) =4.772436582622514e-04;
pointCoord(579, 0) = 1.063118675190957e-01; pointCoord(579, 1) =4.120315967304693e-01; pointCoord(579, 2) =4.772436582622514e-04;
pointCoord(580, 0) = 3.063118675190957e-01; pointCoord(580, 1) =4.120315967304693e-01; pointCoord(580, 2) =4.772436582622514e-04;
pointCoord(581, 0) = 5.063118675190957e-01; pointCoord(581, 1) =4.120315967304693e-01; pointCoord(581, 2) =4.772436582622514e-04;
pointCoord(582, 0) = 1.063118675190957e-01; pointCoord(582, 1) =6.120315967304693e-01; pointCoord(582, 2) =4.772436582622514e-04;
pointCoord(583, 0) = 3.063118675190957e-01; pointCoord(583, 1) =6.120315967304693e-01; pointCoord(583, 2) =4.772436582622514e-04;
pointCoord(584, 0) = 1.063118675190957e-01; pointCoord(584, 1) =8.120315967304693e-01; pointCoord(584, 2) =4.772436582622514e-04;
pointCoord(585, 0) = 1.159937462756005e-01; pointCoord(585, 1) =2.349717973964455e-03; pointCoord(585, 2) =2.172427927521020e-04;
pointCoord(586, 0) = 3.159937462756005e-01; pointCoord(586, 1) =2.349717973964455e-03; pointCoord(586, 2) =2.172427927521020e-04;
pointCoord(587, 0) = 5.159937462756006e-01; pointCoord(587, 1) =2.349717973964455e-03; pointCoord(587, 2) =2.172427927521020e-04;
pointCoord(588, 0) = 7.159937462756005e-01; pointCoord(588, 1) =2.349717973964455e-03; pointCoord(588, 2) =2.172427927521020e-04;
pointCoord(589, 0) = 9.159937462756006e-01; pointCoord(589, 1) =2.349717973964455e-03; pointCoord(589, 2) =2.172427927521020e-04;
pointCoord(590, 0) = 1.159937462756005e-01; pointCoord(590, 1) =2.023497179739645e-01; pointCoord(590, 2) =2.172427927521020e-04;
pointCoord(591, 0) = 3.159937462756005e-01; pointCoord(591, 1) =2.023497179739645e-01; pointCoord(591, 2) =2.172427927521020e-04;
pointCoord(592, 0) = 5.159937462756006e-01; pointCoord(592, 1) =2.023497179739645e-01; pointCoord(592, 2) =2.172427927521020e-04;
pointCoord(593, 0) = 7.159937462756005e-01; pointCoord(593, 1) =2.023497179739645e-01; pointCoord(593, 2) =2.172427927521020e-04;
pointCoord(594, 0) = 1.159937462756005e-01; pointCoord(594, 1) =4.023497179739645e-01; pointCoord(594, 2) =2.172427927521020e-04;
pointCoord(595, 0) = 3.159937462756005e-01; pointCoord(595, 1) =4.023497179739645e-01; pointCoord(595, 2) =2.172427927521020e-04;
pointCoord(596, 0) = 5.159937462756006e-01; pointCoord(596, 1) =4.023497179739645e-01; pointCoord(596, 2) =2.172427927521020e-04;
pointCoord(597, 0) = 1.159937462756005e-01; pointCoord(597, 1) =6.023497179739644e-01; pointCoord(597, 2) =2.172427927521020e-04;
pointCoord(598, 0) = 3.159937462756005e-01; pointCoord(598, 1) =6.023497179739644e-01; pointCoord(598, 2) =2.172427927521020e-04;
pointCoord(599, 0) = 1.159937462756005e-01; pointCoord(599, 1) =8.023497179739645e-01; pointCoord(599, 2) =2.172427927521020e-04;
pointCoord(600, 0) = 3.028955545771836e-03; pointCoord(600, 1) =1.495242854458611e-01; pointCoord(600, 2) =2.422245286926613e-04;
pointCoord(601, 0) = 2.030289555457719e-01; pointCoord(601, 1) =1.495242854458611e-01; pointCoord(601, 2) =2.422245286926613e-04;
pointCoord(602, 0) = 4.030289555457718e-01; pointCoord(602, 1) =1.495242854458611e-01; pointCoord(602, 2) =2.422245286926613e-04;
pointCoord(603, 0) = 6.030289555457718e-01; pointCoord(603, 1) =1.495242854458611e-01; pointCoord(603, 2) =2.422245286926613e-04;
pointCoord(604, 0) = 8.030289555457719e-01; pointCoord(604, 1) =1.495242854458611e-01; pointCoord(604, 2) =2.422245286926613e-04;
pointCoord(605, 0) = 3.028955545771836e-03; pointCoord(605, 1) =3.495242854458611e-01; pointCoord(605, 2) =2.422245286926613e-04;
pointCoord(606, 0) = 2.030289555457719e-01; pointCoord(606, 1) =3.495242854458611e-01; pointCoord(606, 2) =2.422245286926613e-04;
pointCoord(607, 0) = 4.030289555457718e-01; pointCoord(607, 1) =3.495242854458611e-01; pointCoord(607, 2) =2.422245286926613e-04;
pointCoord(608, 0) = 6.030289555457718e-01; pointCoord(608, 1) =3.495242854458611e-01; pointCoord(608, 2) =2.422245286926613e-04;
pointCoord(609, 0) = 3.028955545771836e-03; pointCoord(609, 1) =5.495242854458611e-01; pointCoord(609, 2) =2.422245286926613e-04;
pointCoord(610, 0) = 2.030289555457719e-01; pointCoord(610, 1) =5.495242854458611e-01; pointCoord(610, 2) =2.422245286926613e-04;
pointCoord(611, 0) = 4.030289555457718e-01; pointCoord(611, 1) =5.495242854458611e-01; pointCoord(611, 2) =2.422245286926613e-04;
pointCoord(612, 0) = 3.028955545771836e-03; pointCoord(612, 1) =7.495242854458610e-01; pointCoord(612, 2) =2.422245286926613e-04;
pointCoord(613, 0) = 2.030289555457719e-01; pointCoord(613, 1) =7.495242854458610e-01; pointCoord(613, 2) =2.422245286926613e-04;
pointCoord(614, 0) = 3.028955545771836e-03; pointCoord(614, 1) =9.495242854458611e-01; pointCoord(614, 2) =2.422245286926613e-04;
pointCoord(615, 0) = 1.550959393639832e-02; pointCoord(615, 1) =1.370436470552346e-01; pointCoord(615, 2) =5.321240752324881e-04;
pointCoord(616, 0) = 2.155095939363983e-01; pointCoord(616, 1) =1.370436470552346e-01; pointCoord(616, 2) =5.321240752324881e-04;
pointCoord(617, 0) = 4.155095939363984e-01; pointCoord(617, 1) =1.370436470552346e-01; pointCoord(617, 2) =5.321240752324881e-04;
pointCoord(618, 0) = 6.155095939363983e-01; pointCoord(618, 1) =1.370436470552346e-01; pointCoord(618, 2) =5.321240752324881e-04;
pointCoord(619, 0) = 8.155095939363983e-01; pointCoord(619, 1) =1.370436470552346e-01; pointCoord(619, 2) =5.321240752324881e-04;
pointCoord(620, 0) = 1.550959393639832e-02; pointCoord(620, 1) =3.370436470552346e-01; pointCoord(620, 2) =5.321240752324881e-04;
pointCoord(621, 0) = 2.155095939363983e-01; pointCoord(621, 1) =3.370436470552346e-01; pointCoord(621, 2) =5.321240752324881e-04;
pointCoord(622, 0) = 4.155095939363984e-01; pointCoord(622, 1) =3.370436470552346e-01; pointCoord(622, 2) =5.321240752324881e-04;
pointCoord(623, 0) = 6.155095939363983e-01; pointCoord(623, 1) =3.370436470552346e-01; pointCoord(623, 2) =5.321240752324881e-04;
pointCoord(624, 0) = 1.550959393639832e-02; pointCoord(624, 1) =5.370436470552347e-01; pointCoord(624, 2) =5.321240752324881e-04;
pointCoord(625, 0) = 2.155095939363983e-01; pointCoord(625, 1) =5.370436470552347e-01; pointCoord(625, 2) =5.321240752324881e-04;
pointCoord(626, 0) = 4.155095939363984e-01; pointCoord(626, 1) =5.370436470552347e-01; pointCoord(626, 2) =5.321240752324881e-04;
pointCoord(627, 0) = 1.550959393639832e-02; pointCoord(627, 1) =7.370436470552346e-01; pointCoord(627, 2) =5.321240752324881e-04;
pointCoord(628, 0) = 2.155095939363983e-01; pointCoord(628, 1) =7.370436470552346e-01; pointCoord(628, 2) =5.321240752324881e-04;
pointCoord(629, 0) = 1.550959393639832e-02; pointCoord(629, 1) =9.370436470552346e-01; pointCoord(629, 2) =5.321240752324881e-04;
pointCoord(630, 0) = 3.619078430637677e-02; pointCoord(630, 1) =1.163624566852561e-01; pointCoord(630, 2) =7.506524072180082e-04;
pointCoord(631, 0) = 2.361907843063768e-01; pointCoord(631, 1) =1.163624566852561e-01; pointCoord(631, 2) =7.506524072180082e-04;
pointCoord(632, 0) = 4.361907843063768e-01; pointCoord(632, 1) =1.163624566852561e-01; pointCoord(632, 2) =7.506524072180082e-04;
pointCoord(633, 0) = 6.361907843063768e-01; pointCoord(633, 1) =1.163624566852561e-01; pointCoord(633, 2) =7.506524072180082e-04;
pointCoord(634, 0) = 8.361907843063768e-01; pointCoord(634, 1) =1.163624566852561e-01; pointCoord(634, 2) =7.506524072180082e-04;
pointCoord(635, 0) = 3.619078430637677e-02; pointCoord(635, 1) =3.163624566852561e-01; pointCoord(635, 2) =7.506524072180082e-04;
pointCoord(636, 0) = 2.361907843063768e-01; pointCoord(636, 1) =3.163624566852561e-01; pointCoord(636, 2) =7.506524072180082e-04;
pointCoord(637, 0) = 4.361907843063768e-01; pointCoord(637, 1) =3.163624566852561e-01; pointCoord(637, 2) =7.506524072180082e-04;
pointCoord(638, 0) = 6.361907843063768e-01; pointCoord(638, 1) =3.163624566852561e-01; pointCoord(638, 2) =7.506524072180082e-04;
pointCoord(639, 0) = 3.619078430637677e-02; pointCoord(639, 1) =5.163624566852562e-01; pointCoord(639, 2) =7.506524072180082e-04;
pointCoord(640, 0) = 2.361907843063768e-01; pointCoord(640, 1) =5.163624566852562e-01; pointCoord(640, 2) =7.506524072180082e-04;
pointCoord(641, 0) = 4.361907843063768e-01; pointCoord(641, 1) =5.163624566852562e-01; pointCoord(641, 2) =7.506524072180082e-04;
pointCoord(642, 0) = 3.619078430637677e-02; pointCoord(642, 1) =7.163624566852561e-01; pointCoord(642, 2) =7.506524072180082e-04;
pointCoord(643, 0) = 2.361907843063768e-01; pointCoord(643, 1) =7.163624566852561e-01; pointCoord(643, 2) =7.506524072180082e-04;
pointCoord(644, 0) = 3.619078430637677e-02; pointCoord(644, 1) =9.163624566852562e-01; pointCoord(644, 2) =7.506524072180082e-04;
pointCoord(645, 0) = 6.228484588439000e-02; pointCoord(645, 1) =9.026839510724288e-02; pointCoord(645, 2) =8.678472663211519e-04;
pointCoord(646, 0) = 2.622848458843900e-01; pointCoord(646, 1) =9.026839510724288e-02; pointCoord(646, 2) =8.678472663211519e-04;
pointCoord(647, 0) = 4.622848458843900e-01; pointCoord(647, 1) =9.026839510724288e-02; pointCoord(647, 2) =8.678472663211519e-04;
pointCoord(648, 0) = 6.622848458843900e-01; pointCoord(648, 1) =9.026839510724288e-02; pointCoord(648, 2) =8.678472663211519e-04;
pointCoord(649, 0) = 8.622848458843900e-01; pointCoord(649, 1) =9.026839510724288e-02; pointCoord(649, 2) =8.678472663211519e-04;
pointCoord(650, 0) = 6.228484588439000e-02; pointCoord(650, 1) =2.902683951072429e-01; pointCoord(650, 2) =8.678472663211519e-04;
pointCoord(651, 0) = 2.622848458843900e-01; pointCoord(651, 1) =2.902683951072429e-01; pointCoord(651, 2) =8.678472663211519e-04;
pointCoord(652, 0) = 4.622848458843900e-01; pointCoord(652, 1) =2.902683951072429e-01; pointCoord(652, 2) =8.678472663211519e-04;
pointCoord(653, 0) = 6.622848458843900e-01; pointCoord(653, 1) =2.902683951072429e-01; pointCoord(653, 2) =8.678472663211519e-04;
pointCoord(654, 0) = 6.228484588439000e-02; pointCoord(654, 1) =4.902683951072429e-01; pointCoord(654, 2) =8.678472663211519e-04;
pointCoord(655, 0) = 2.622848458843900e-01; pointCoord(655, 1) =4.902683951072429e-01; pointCoord(655, 2) =8.678472663211519e-04;
pointCoord(656, 0) = 4.622848458843900e-01; pointCoord(656, 1) =4.902683951072429e-01; pointCoord(656, 2) =8.678472663211519e-04;
pointCoord(657, 0) = 6.228484588439000e-02; pointCoord(657, 1) =6.902683951072428e-01; pointCoord(657, 2) =8.678472663211519e-04;
pointCoord(658, 0) = 2.622848458843900e-01; pointCoord(658, 1) =6.902683951072428e-01; pointCoord(658, 2) =8.678472663211519e-04;
pointCoord(659, 0) = 6.228484588439000e-02; pointCoord(659, 1) =8.902683951072430e-01; pointCoord(659, 2) =8.678472663211519e-04;
pointCoord(660, 0) = 9.026839510724288e-02; pointCoord(660, 1) =6.228484588439000e-02; pointCoord(660, 2) =8.678472663211519e-04;
pointCoord(661, 0) = 2.902683951072429e-01; pointCoord(661, 1) =6.228484588439000e-02; pointCoord(661, 2) =8.678472663211519e-04;
pointCoord(662, 0) = 4.902683951072429e-01; pointCoord(662, 1) =6.228484588439000e-02; pointCoord(662, 2) =8.678472663211519e-04;
pointCoord(663, 0) = 6.902683951072428e-01; pointCoord(663, 1) =6.228484588439000e-02; pointCoord(663, 2) =8.678472663211519e-04;
pointCoord(664, 0) = 8.902683951072430e-01; pointCoord(664, 1) =6.228484588439000e-02; pointCoord(664, 2) =8.678472663211519e-04;
pointCoord(665, 0) = 9.026839510724288e-02; pointCoord(665, 1) =2.622848458843900e-01; pointCoord(665, 2) =8.678472663211519e-04;
pointCoord(666, 0) = 2.902683951072429e-01; pointCoord(666, 1) =2.622848458843900e-01; pointCoord(666, 2) =8.678472663211519e-04;
pointCoord(667, 0) = 4.902683951072429e-01; pointCoord(667, 1) =2.622848458843900e-01; pointCoord(667, 2) =8.678472663211519e-04;
pointCoord(668, 0) = 6.902683951072428e-01; pointCoord(668, 1) =2.622848458843900e-01; pointCoord(668, 2) =8.678472663211519e-04;
pointCoord(669, 0) = 9.026839510724288e-02; pointCoord(669, 1) =4.622848458843900e-01; pointCoord(669, 2) =8.678472663211519e-04;
pointCoord(670, 0) = 2.902683951072429e-01; pointCoord(670, 1) =4.622848458843900e-01; pointCoord(670, 2) =8.678472663211519e-04;
pointCoord(671, 0) = 4.902683951072429e-01; pointCoord(671, 1) =4.622848458843900e-01; pointCoord(671, 2) =8.678472663211519e-04;
pointCoord(672, 0) = 9.026839510724288e-02; pointCoord(672, 1) =6.622848458843900e-01; pointCoord(672, 2) =8.678472663211519e-04;
pointCoord(673, 0) = 2.902683951072429e-01; pointCoord(673, 1) =6.622848458843900e-01; pointCoord(673, 2) =8.678472663211519e-04;
pointCoord(674, 0) = 9.026839510724288e-02; pointCoord(674, 1) =8.622848458843900e-01; pointCoord(674, 2) =8.678472663211519e-04;
pointCoord(675, 0) = 1.163624566852561e-01; pointCoord(675, 1) =3.619078430637677e-02; pointCoord(675, 2) =7.506524072180082e-04;
pointCoord(676, 0) = 3.163624566852561e-01; pointCoord(676, 1) =3.619078430637677e-02; pointCoord(676, 2) =7.506524072180082e-04;
pointCoord(677, 0) = 5.163624566852562e-01; pointCoord(677, 1) =3.619078430637677e-02; pointCoord(677, 2) =7.506524072180082e-04;
pointCoord(678, 0) = 7.163624566852561e-01; pointCoord(678, 1) =3.619078430637677e-02; pointCoord(678, 2) =7.506524072180082e-04;
pointCoord(679, 0) = 9.163624566852562e-01; pointCoord(679, 1) =3.619078430637677e-02; pointCoord(679, 2) =7.506524072180082e-04;
pointCoord(680, 0) = 1.163624566852561e-01; pointCoord(680, 1) =2.361907843063768e-01; pointCoord(680, 2) =7.506524072180082e-04;
pointCoord(681, 0) = 3.163624566852561e-01; pointCoord(681, 1) =2.361907843063768e-01; pointCoord(681, 2) =7.506524072180082e-04;
pointCoord(682, 0) = 5.163624566852562e-01; pointCoord(682, 1) =2.361907843063768e-01; pointCoord(682, 2) =7.506524072180082e-04;
pointCoord(683, 0) = 7.163624566852561e-01; pointCoord(683, 1) =2.361907843063768e-01; pointCoord(683, 2) =7.506524072180082e-04;
pointCoord(684, 0) = 1.163624566852561e-01; pointCoord(684, 1) =4.361907843063768e-01; pointCoord(684, 2) =7.506524072180082e-04;
pointCoord(685, 0) = 3.163624566852561e-01; pointCoord(685, 1) =4.361907843063768e-01; pointCoord(685, 2) =7.506524072180082e-04;
pointCoord(686, 0) = 5.163624566852562e-01; pointCoord(686, 1) =4.361907843063768e-01; pointCoord(686, 2) =7.506524072180082e-04;
pointCoord(687, 0) = 1.163624566852561e-01; pointCoord(687, 1) =6.361907843063768e-01; pointCoord(687, 2) =7.506524072180082e-04;
pointCoord(688, 0) = 3.163624566852561e-01; pointCoord(688, 1) =6.361907843063768e-01; pointCoord(688, 2) =7.506524072180082e-04;
pointCoord(689, 0) = 1.163624566852561e-01; pointCoord(689, 1) =8.361907843063768e-01; pointCoord(689, 2) =7.506524072180082e-04;
pointCoord(690, 0) = 1.370436470552346e-01; pointCoord(690, 1) =1.550959393639832e-02; pointCoord(690, 2) =5.321240752324881e-04;
pointCoord(691, 0) = 3.370436470552346e-01; pointCoord(691, 1) =1.550959393639832e-02; pointCoord(691, 2) =5.321240752324881e-04;
pointCoord(692, 0) = 5.370436470552347e-01; pointCoord(692, 1) =1.550959393639832e-02; pointCoord(692, 2) =5.321240752324881e-04;
pointCoord(693, 0) = 7.370436470552346e-01; pointCoord(693, 1) =1.550959393639832e-02; pointCoord(693, 2) =5.321240752324881e-04;
pointCoord(694, 0) = 9.370436470552346e-01; pointCoord(694, 1) =1.550959393639832e-02; pointCoord(694, 2) =5.321240752324881e-04;
pointCoord(695, 0) = 1.370436470552346e-01; pointCoord(695, 1) =2.155095939363983e-01; pointCoord(695, 2) =5.321240752324881e-04;
pointCoord(696, 0) = 3.370436470552346e-01; pointCoord(696, 1) =2.155095939363983e-01; pointCoord(696, 2) =5.321240752324881e-04;
pointCoord(697, 0) = 5.370436470552347e-01; pointCoord(697, 1) =2.155095939363983e-01; pointCoord(697, 2) =5.321240752324881e-04;
pointCoord(698, 0) = 7.370436470552346e-01; pointCoord(698, 1) =2.155095939363983e-01; pointCoord(698, 2) =5.321240752324881e-04;
pointCoord(699, 0) = 1.370436470552346e-01; pointCoord(699, 1) =4.155095939363984e-01; pointCoord(699, 2) =5.321240752324881e-04;
pointCoord(700, 0) = 3.370436470552346e-01; pointCoord(700, 1) =4.155095939363984e-01; pointCoord(700, 2) =5.321240752324881e-04;
pointCoord(701, 0) = 5.370436470552347e-01; pointCoord(701, 1) =4.155095939363984e-01; pointCoord(701, 2) =5.321240752324881e-04;
pointCoord(702, 0) = 1.370436470552346e-01; pointCoord(702, 1) =6.155095939363983e-01; pointCoord(702, 2) =5.321240752324881e-04;
pointCoord(703, 0) = 3.370436470552346e-01; pointCoord(703, 1) =6.155095939363983e-01; pointCoord(703, 2) =5.321240752324881e-04;
pointCoord(704, 0) = 1.370436470552346e-01; pointCoord(704, 1) =8.155095939363983e-01; pointCoord(704, 2) =5.321240752324881e-04;
pointCoord(705, 0) = 1.495242854458611e-01; pointCoord(705, 1) =3.028955545771836e-03; pointCoord(705, 2) =2.422245286926613e-04;
pointCoord(706, 0) = 3.495242854458611e-01; pointCoord(706, 1) =3.028955545771836e-03; pointCoord(706, 2) =2.422245286926613e-04;
pointCoord(707, 0) = 5.495242854458611e-01; pointCoord(707, 1) =3.028955545771836e-03; pointCoord(707, 2) =2.422245286926613e-04;
pointCoord(708, 0) = 7.495242854458610e-01; pointCoord(708, 1) =3.028955545771836e-03; pointCoord(708, 2) =2.422245286926613e-04;
pointCoord(709, 0) = 9.495242854458611e-01; pointCoord(709, 1) =3.028955545771836e-03; pointCoord(709, 2) =2.422245286926613e-04;
pointCoord(710, 0) = 1.495242854458611e-01; pointCoord(710, 1) =2.030289555457719e-01; pointCoord(710, 2) =2.422245286926613e-04;
pointCoord(711, 0) = 3.495242854458611e-01; pointCoord(711, 1) =2.030289555457719e-01; pointCoord(711, 2) =2.422245286926613e-04;
pointCoord(712, 0) = 5.495242854458611e-01; pointCoord(712, 1) =2.030289555457719e-01; pointCoord(712, 2) =2.422245286926613e-04;
pointCoord(713, 0) = 7.495242854458610e-01; pointCoord(713, 1) =2.030289555457719e-01; pointCoord(713, 2) =2.422245286926613e-04;
pointCoord(714, 0) = 1.495242854458611e-01; pointCoord(714, 1) =4.030289555457718e-01; pointCoord(714, 2) =2.422245286926613e-04;
pointCoord(715, 0) = 3.495242854458611e-01; pointCoord(715, 1) =4.030289555457718e-01; pointCoord(715, 2) =2.422245286926613e-04;
pointCoord(716, 0) = 5.495242854458611e-01; pointCoord(716, 1) =4.030289555457718e-01; pointCoord(716, 2) =2.422245286926613e-04;
pointCoord(717, 0) = 1.495242854458611e-01; pointCoord(717, 1) =6.030289555457718e-01; pointCoord(717, 2) =2.422245286926613e-04;
pointCoord(718, 0) = 3.495242854458611e-01; pointCoord(718, 1) =6.030289555457718e-01; pointCoord(718, 2) =2.422245286926613e-04;
pointCoord(719, 0) = 1.495242854458611e-01; pointCoord(719, 1) =8.030289555457719e-01; pointCoord(719, 2) =2.422245286926613e-04;
pointCoord(720, 0) = 3.567294182208055e-03; pointCoord(720, 1) =1.760993535591546e-01; pointCoord(720, 2) =2.022265498028209e-04;
pointCoord(721, 0) = 2.035672941822081e-01; pointCoord(721, 1) =1.760993535591546e-01; pointCoord(721, 2) =2.022265498028209e-04;
pointCoord(722, 0) = 4.035672941822081e-01; pointCoord(722, 1) =1.760993535591546e-01; pointCoord(722, 2) =2.022265498028209e-04;
pointCoord(723, 0) = 6.035672941822080e-01; pointCoord(723, 1) =1.760993535591546e-01; pointCoord(723, 2) =2.022265498028209e-04;
pointCoord(724, 0) = 8.035672941822081e-01; pointCoord(724, 1) =1.760993535591546e-01; pointCoord(724, 2) =2.022265498028209e-04;
pointCoord(725, 0) = 3.567294182208055e-03; pointCoord(725, 1) =3.760993535591546e-01; pointCoord(725, 2) =2.022265498028209e-04;
pointCoord(726, 0) = 2.035672941822081e-01; pointCoord(726, 1) =3.760993535591546e-01; pointCoord(726, 2) =2.022265498028209e-04;
pointCoord(727, 0) = 4.035672941822081e-01; pointCoord(727, 1) =3.760993535591546e-01; pointCoord(727, 2) =2.022265498028209e-04;
pointCoord(728, 0) = 6.035672941822080e-01; pointCoord(728, 1) =3.760993535591546e-01; pointCoord(728, 2) =2.022265498028209e-04;
pointCoord(729, 0) = 3.567294182208055e-03; pointCoord(729, 1) =5.760993535591546e-01; pointCoord(729, 2) =2.022265498028209e-04;
pointCoord(730, 0) = 2.035672941822081e-01; pointCoord(730, 1) =5.760993535591546e-01; pointCoord(730, 2) =2.022265498028209e-04;
pointCoord(731, 0) = 4.035672941822081e-01; pointCoord(731, 1) =5.760993535591546e-01; pointCoord(731, 2) =2.022265498028209e-04;
pointCoord(732, 0) = 3.567294182208055e-03; pointCoord(732, 1) =7.760993535591546e-01; pointCoord(732, 2) =2.022265498028209e-04;
pointCoord(733, 0) = 2.035672941822081e-01; pointCoord(733, 1) =7.760993535591546e-01; pointCoord(733, 2) =2.022265498028209e-04;
pointCoord(734, 0) = 3.567294182208055e-03; pointCoord(734, 1) =9.760993535591547e-01; pointCoord(734, 2) =2.022265498028209e-04;
pointCoord(735, 0) = 1.826612618826817e-02; pointCoord(735, 1) =1.614005215530945e-01; pointCoord(735, 2) =4.442556514902738e-04;
pointCoord(736, 0) = 2.182661261882682e-01; pointCoord(736, 1) =1.614005215530945e-01; pointCoord(736, 2) =4.442556514902738e-04;
pointCoord(737, 0) = 4.182661261882682e-01; pointCoord(737, 1) =1.614005215530945e-01; pointCoord(737, 2) =4.442556514902738e-04;
pointCoord(738, 0) = 6.182661261882682e-01; pointCoord(738, 1) =1.614005215530945e-01; pointCoord(738, 2) =4.442556514902738e-04;
pointCoord(739, 0) = 8.182661261882682e-01; pointCoord(739, 1) =1.614005215530945e-01; pointCoord(739, 2) =4.442556514902738e-04;
pointCoord(740, 0) = 1.826612618826817e-02; pointCoord(740, 1) =3.614005215530945e-01; pointCoord(740, 2) =4.442556514902738e-04;
pointCoord(741, 0) = 2.182661261882682e-01; pointCoord(741, 1) =3.614005215530945e-01; pointCoord(741, 2) =4.442556514902738e-04;
pointCoord(742, 0) = 4.182661261882682e-01; pointCoord(742, 1) =3.614005215530945e-01; pointCoord(742, 2) =4.442556514902738e-04;
pointCoord(743, 0) = 6.182661261882682e-01; pointCoord(743, 1) =3.614005215530945e-01; pointCoord(743, 2) =4.442556514902738e-04;
pointCoord(744, 0) = 1.826612618826817e-02; pointCoord(744, 1) =5.614005215530945e-01; pointCoord(744, 2) =4.442556514902738e-04;
pointCoord(745, 0) = 2.182661261882682e-01; pointCoord(745, 1) =5.614005215530945e-01; pointCoord(745, 2) =4.442556514902738e-04;
pointCoord(746, 0) = 4.182661261882682e-01; pointCoord(746, 1) =5.614005215530945e-01; pointCoord(746, 2) =4.442556514902738e-04;
pointCoord(747, 0) = 1.826612618826817e-02; pointCoord(747, 1) =7.614005215530945e-01; pointCoord(747, 2) =4.442556514902738e-04;
pointCoord(748, 0) = 2.182661261882682e-01; pointCoord(748, 1) =7.614005215530945e-01; pointCoord(748, 2) =4.442556514902738e-04;
pointCoord(749, 0) = 1.826612618826817e-02; pointCoord(749, 1) =9.614005215530945e-01; pointCoord(749, 2) =4.442556514902738e-04;
pointCoord(750, 0) = 4.262300068612809e-02; pointCoord(750, 1) =1.370436470552346e-01; pointCoord(750, 2) =6.266989030061803e-04;
pointCoord(751, 0) = 2.426230006861281e-01; pointCoord(751, 1) =1.370436470552346e-01; pointCoord(751, 2) =6.266989030061803e-04;
pointCoord(752, 0) = 4.426230006861281e-01; pointCoord(752, 1) =1.370436470552346e-01; pointCoord(752, 2) =6.266989030061803e-04;
pointCoord(753, 0) = 6.426230006861281e-01; pointCoord(753, 1) =1.370436470552346e-01; pointCoord(753, 2) =6.266989030061803e-04;
pointCoord(754, 0) = 8.426230006861282e-01; pointCoord(754, 1) =1.370436470552346e-01; pointCoord(754, 2) =6.266989030061803e-04;
pointCoord(755, 0) = 4.262300068612809e-02; pointCoord(755, 1) =3.370436470552346e-01; pointCoord(755, 2) =6.266989030061803e-04;
pointCoord(756, 0) = 2.426230006861281e-01; pointCoord(756, 1) =3.370436470552346e-01; pointCoord(756, 2) =6.266989030061803e-04;
pointCoord(757, 0) = 4.426230006861281e-01; pointCoord(757, 1) =3.370436470552346e-01; pointCoord(757, 2) =6.266989030061803e-04;
pointCoord(758, 0) = 6.426230006861281e-01; pointCoord(758, 1) =3.370436470552346e-01; pointCoord(758, 2) =6.266989030061803e-04;
pointCoord(759, 0) = 4.262300068612809e-02; pointCoord(759, 1) =5.370436470552347e-01; pointCoord(759, 2) =6.266989030061803e-04;
pointCoord(760, 0) = 2.426230006861281e-01; pointCoord(760, 1) =5.370436470552347e-01; pointCoord(760, 2) =6.266989030061803e-04;
pointCoord(761, 0) = 4.426230006861281e-01; pointCoord(761, 1) =5.370436470552347e-01; pointCoord(761, 2) =6.266989030061803e-04;
pointCoord(762, 0) = 4.262300068612809e-02; pointCoord(762, 1) =7.370436470552346e-01; pointCoord(762, 2) =6.266989030061803e-04;
pointCoord(763, 0) = 2.426230006861281e-01; pointCoord(763, 1) =7.370436470552346e-01; pointCoord(763, 2) =6.266989030061803e-04;
pointCoord(764, 0) = 4.262300068612809e-02; pointCoord(764, 1) =9.370436470552346e-01; pointCoord(764, 2) =6.266989030061803e-04;
pointCoord(765, 0) = 7.335478022226698e-02; pointCoord(765, 1) =1.063118675190957e-01; pointCoord(765, 2) =7.245416447754392e-04;
pointCoord(766, 0) = 2.733547802222670e-01; pointCoord(766, 1) =1.063118675190957e-01; pointCoord(766, 2) =7.245416447754392e-04;
pointCoord(767, 0) = 4.733547802222670e-01; pointCoord(767, 1) =1.063118675190957e-01; pointCoord(767, 2) =7.245416447754392e-04;
pointCoord(768, 0) = 6.733547802222669e-01; pointCoord(768, 1) =1.063118675190957e-01; pointCoord(768, 2) =7.245416447754392e-04;
pointCoord(769, 0) = 8.733547802222670e-01; pointCoord(769, 1) =1.063118675190957e-01; pointCoord(769, 2) =7.245416447754392e-04;
pointCoord(770, 0) = 7.335478022226698e-02; pointCoord(770, 1) =3.063118675190957e-01; pointCoord(770, 2) =7.245416447754392e-04;
pointCoord(771, 0) = 2.733547802222670e-01; pointCoord(771, 1) =3.063118675190957e-01; pointCoord(771, 2) =7.245416447754392e-04;
pointCoord(772, 0) = 4.733547802222670e-01; pointCoord(772, 1) =3.063118675190957e-01; pointCoord(772, 2) =7.245416447754392e-04;
pointCoord(773, 0) = 6.733547802222669e-01; pointCoord(773, 1) =3.063118675190957e-01; pointCoord(773, 2) =7.245416447754392e-04;
pointCoord(774, 0) = 7.335478022226698e-02; pointCoord(774, 1) =5.063118675190957e-01; pointCoord(774, 2) =7.245416447754392e-04;
pointCoord(775, 0) = 2.733547802222670e-01; pointCoord(775, 1) =5.063118675190957e-01; pointCoord(775, 2) =7.245416447754392e-04;
pointCoord(776, 0) = 4.733547802222670e-01; pointCoord(776, 1) =5.063118675190957e-01; pointCoord(776, 2) =7.245416447754392e-04;
pointCoord(777, 0) = 7.335478022226698e-02; pointCoord(777, 1) =7.063118675190957e-01; pointCoord(777, 2) =7.245416447754392e-04;
pointCoord(778, 0) = 2.733547802222670e-01; pointCoord(778, 1) =7.063118675190957e-01; pointCoord(778, 2) =7.245416447754392e-04;
pointCoord(779, 0) = 7.335478022226698e-02; pointCoord(779, 1) =9.063118675190958e-01; pointCoord(779, 2) =7.245416447754392e-04;
pointCoord(780, 0) = 1.063118675190957e-01; pointCoord(780, 1) =7.335478022226698e-02; pointCoord(780, 2) =7.245416447754392e-04;
pointCoord(781, 0) = 3.063118675190957e-01; pointCoord(781, 1) =7.335478022226698e-02; pointCoord(781, 2) =7.245416447754392e-04;
pointCoord(782, 0) = 5.063118675190957e-01; pointCoord(782, 1) =7.335478022226698e-02; pointCoord(782, 2) =7.245416447754392e-04;
pointCoord(783, 0) = 7.063118675190957e-01; pointCoord(783, 1) =7.335478022226698e-02; pointCoord(783, 2) =7.245416447754392e-04;
pointCoord(784, 0) = 9.063118675190958e-01; pointCoord(784, 1) =7.335478022226698e-02; pointCoord(784, 2) =7.245416447754392e-04;
pointCoord(785, 0) = 1.063118675190957e-01; pointCoord(785, 1) =2.733547802222670e-01; pointCoord(785, 2) =7.245416447754392e-04;
pointCoord(786, 0) = 3.063118675190957e-01; pointCoord(786, 1) =2.733547802222670e-01; pointCoord(786, 2) =7.245416447754392e-04;
pointCoord(787, 0) = 5.063118675190957e-01; pointCoord(787, 1) =2.733547802222670e-01; pointCoord(787, 2) =7.245416447754392e-04;
pointCoord(788, 0) = 7.063118675190957e-01; pointCoord(788, 1) =2.733547802222670e-01; pointCoord(788, 2) =7.245416447754392e-04;
pointCoord(789, 0) = 1.063118675190957e-01; pointCoord(789, 1) =4.733547802222670e-01; pointCoord(789, 2) =7.245416447754392e-04;
pointCoord(790, 0) = 3.063118675190957e-01; pointCoord(790, 1) =4.733547802222670e-01; pointCoord(790, 2) =7.245416447754392e-04;
pointCoord(791, 0) = 5.063118675190957e-01; pointCoord(791, 1) =4.733547802222670e-01; pointCoord(791, 2) =7.245416447754392e-04;
pointCoord(792, 0) = 1.063118675190957e-01; pointCoord(792, 1) =6.733547802222669e-01; pointCoord(792, 2) =7.245416447754392e-04;
pointCoord(793, 0) = 3.063118675190957e-01; pointCoord(793, 1) =6.733547802222669e-01; pointCoord(793, 2) =7.245416447754392e-04;
pointCoord(794, 0) = 1.063118675190957e-01; pointCoord(794, 1) =8.733547802222670e-01; pointCoord(794, 2) =7.245416447754392e-04;
pointCoord(795, 0) = 1.370436470552346e-01; pointCoord(795, 1) =4.262300068612809e-02; pointCoord(795, 2) =6.266989030061803e-04;
pointCoord(796, 0) = 3.370436470552346e-01; pointCoord(796, 1) =4.262300068612809e-02; pointCoord(796, 2) =6.266989030061803e-04;
pointCoord(797, 0) = 5.370436470552347e-01; pointCoord(797, 1) =4.262300068612809e-02; pointCoord(797, 2) =6.266989030061803e-04;
pointCoord(798, 0) = 7.370436470552346e-01; pointCoord(798, 1) =4.262300068612809e-02; pointCoord(798, 2) =6.266989030061803e-04;
pointCoord(799, 0) = 9.370436470552346e-01; pointCoord(799, 1) =4.262300068612809e-02; pointCoord(799, 2) =6.266989030061803e-04;
pointCoord(800, 0) = 1.370436470552346e-01; pointCoord(800, 1) =2.426230006861281e-01; pointCoord(800, 2) =6.266989030061803e-04;
pointCoord(801, 0) = 3.370436470552346e-01; pointCoord(801, 1) =2.426230006861281e-01; pointCoord(801, 2) =6.266989030061803e-04;
pointCoord(802, 0) = 5.370436470552347e-01; pointCoord(802, 1) =2.426230006861281e-01; pointCoord(802, 2) =6.266989030061803e-04;
pointCoord(803, 0) = 7.370436470552346e-01; pointCoord(803, 1) =2.426230006861281e-01; pointCoord(803, 2) =6.266989030061803e-04;
pointCoord(804, 0) = 1.370436470552346e-01; pointCoord(804, 1) =4.426230006861281e-01; pointCoord(804, 2) =6.266989030061803e-04;
pointCoord(805, 0) = 3.370436470552346e-01; pointCoord(805, 1) =4.426230006861281e-01; pointCoord(805, 2) =6.266989030061803e-04;
pointCoord(806, 0) = 5.370436470552347e-01; pointCoord(806, 1) =4.426230006861281e-01; pointCoord(806, 2) =6.266989030061803e-04;
pointCoord(807, 0) = 1.370436470552346e-01; pointCoord(807, 1) =6.426230006861281e-01; pointCoord(807, 2) =6.266989030061803e-04;
pointCoord(808, 0) = 3.370436470552346e-01; pointCoord(808, 1) =6.426230006861281e-01; pointCoord(808, 2) =6.266989030061803e-04;
pointCoord(809, 0) = 1.370436470552346e-01; pointCoord(809, 1) =8.426230006861282e-01; pointCoord(809, 2) =6.266989030061803e-04;
pointCoord(810, 0) = 1.614005215530945e-01; pointCoord(810, 1) =1.826612618826817e-02; pointCoord(810, 2) =4.442556514902738e-04;
pointCoord(811, 0) = 3.614005215530945e-01; pointCoord(811, 1) =1.826612618826817e-02; pointCoord(811, 2) =4.442556514902738e-04;
pointCoord(812, 0) = 5.614005215530945e-01; pointCoord(812, 1) =1.826612618826817e-02; pointCoord(812, 2) =4.442556514902738e-04;
pointCoord(813, 0) = 7.614005215530945e-01; pointCoord(813, 1) =1.826612618826817e-02; pointCoord(813, 2) =4.442556514902738e-04;
pointCoord(814, 0) = 9.614005215530945e-01; pointCoord(814, 1) =1.826612618826817e-02; pointCoord(814, 2) =4.442556514902738e-04;
pointCoord(815, 0) = 1.614005215530945e-01; pointCoord(815, 1) =2.182661261882682e-01; pointCoord(815, 2) =4.442556514902738e-04;
pointCoord(816, 0) = 3.614005215530945e-01; pointCoord(816, 1) =2.182661261882682e-01; pointCoord(816, 2) =4.442556514902738e-04;
pointCoord(817, 0) = 5.614005215530945e-01; pointCoord(817, 1) =2.182661261882682e-01; pointCoord(817, 2) =4.442556514902738e-04;
pointCoord(818, 0) = 7.614005215530945e-01; pointCoord(818, 1) =2.182661261882682e-01; pointCoord(818, 2) =4.442556514902738e-04;
pointCoord(819, 0) = 1.614005215530945e-01; pointCoord(819, 1) =4.182661261882682e-01; pointCoord(819, 2) =4.442556514902738e-04;
pointCoord(820, 0) = 3.614005215530945e-01; pointCoord(820, 1) =4.182661261882682e-01; pointCoord(820, 2) =4.442556514902738e-04;
pointCoord(821, 0) = 5.614005215530945e-01; pointCoord(821, 1) =4.182661261882682e-01; pointCoord(821, 2) =4.442556514902738e-04;
pointCoord(822, 0) = 1.614005215530945e-01; pointCoord(822, 1) =6.182661261882682e-01; pointCoord(822, 2) =4.442556514902738e-04;
pointCoord(823, 0) = 3.614005215530945e-01; pointCoord(823, 1) =6.182661261882682e-01; pointCoord(823, 2) =4.442556514902738e-04;
pointCoord(824, 0) = 1.614005215530945e-01; pointCoord(824, 1) =8.182661261882682e-01; pointCoord(824, 2) =4.442556514902738e-04;
pointCoord(825, 0) = 1.760993535591546e-01; pointCoord(825, 1) =3.567294182208055e-03; pointCoord(825, 2) =2.022265498028209e-04;
pointCoord(826, 0) = 3.760993535591546e-01; pointCoord(826, 1) =3.567294182208055e-03; pointCoord(826, 2) =2.022265498028209e-04;
pointCoord(827, 0) = 5.760993535591546e-01; pointCoord(827, 1) =3.567294182208055e-03; pointCoord(827, 2) =2.022265498028209e-04;
pointCoord(828, 0) = 7.760993535591546e-01; pointCoord(828, 1) =3.567294182208055e-03; pointCoord(828, 2) =2.022265498028209e-04;
pointCoord(829, 0) = 9.760993535591547e-01; pointCoord(829, 1) =3.567294182208055e-03; pointCoord(829, 2) =2.022265498028209e-04;
pointCoord(830, 0) = 1.760993535591546e-01; pointCoord(830, 1) =2.035672941822081e-01; pointCoord(830, 2) =2.022265498028209e-04;
pointCoord(831, 0) = 3.760993535591546e-01; pointCoord(831, 1) =2.035672941822081e-01; pointCoord(831, 2) =2.022265498028209e-04;
pointCoord(832, 0) = 5.760993535591546e-01; pointCoord(832, 1) =2.035672941822081e-01; pointCoord(832, 2) =2.022265498028209e-04;
pointCoord(833, 0) = 7.760993535591546e-01; pointCoord(833, 1) =2.035672941822081e-01; pointCoord(833, 2) =2.022265498028209e-04;
pointCoord(834, 0) = 1.760993535591546e-01; pointCoord(834, 1) =4.035672941822081e-01; pointCoord(834, 2) =2.022265498028209e-04;
pointCoord(835, 0) = 3.760993535591546e-01; pointCoord(835, 1) =4.035672941822081e-01; pointCoord(835, 2) =2.022265498028209e-04;
pointCoord(836, 0) = 5.760993535591546e-01; pointCoord(836, 1) =4.035672941822081e-01; pointCoord(836, 2) =2.022265498028209e-04;
pointCoord(837, 0) = 1.760993535591546e-01; pointCoord(837, 1) =6.035672941822080e-01; pointCoord(837, 2) =2.022265498028209e-04;
pointCoord(838, 0) = 3.760993535591546e-01; pointCoord(838, 1) =6.035672941822080e-01; pointCoord(838, 2) =2.022265498028209e-04;
pointCoord(839, 0) = 1.760993535591546e-01; pointCoord(839, 1) =8.035672941822081e-01; pointCoord(839, 2) =2.022265498028209e-04;
pointCoord(840, 0) = 3.892169575397059e-03; pointCoord(840, 1) =1.921368160743566e-01; pointCoord(840, 2) =1.004375733945300e-04;
pointCoord(841, 0) = 2.038921695753971e-01; pointCoord(841, 1) =1.921368160743566e-01; pointCoord(841, 2) =1.004375733945300e-04;
pointCoord(842, 0) = 4.038921695753971e-01; pointCoord(842, 1) =1.921368160743566e-01; pointCoord(842, 2) =1.004375733945300e-04;
pointCoord(843, 0) = 6.038921695753970e-01; pointCoord(843, 1) =1.921368160743566e-01; pointCoord(843, 2) =1.004375733945300e-04;
pointCoord(844, 0) = 8.038921695753971e-01; pointCoord(844, 1) =1.921368160743566e-01; pointCoord(844, 2) =1.004375733945300e-04;
pointCoord(845, 0) = 3.892169575397059e-03; pointCoord(845, 1) =3.921368160743566e-01; pointCoord(845, 2) =1.004375733945300e-04;
pointCoord(846, 0) = 2.038921695753971e-01; pointCoord(846, 1) =3.921368160743566e-01; pointCoord(846, 2) =1.004375733945300e-04;
pointCoord(847, 0) = 4.038921695753971e-01; pointCoord(847, 1) =3.921368160743566e-01; pointCoord(847, 2) =1.004375733945300e-04;
pointCoord(848, 0) = 6.038921695753970e-01; pointCoord(848, 1) =3.921368160743566e-01; pointCoord(848, 2) =1.004375733945300e-04;
pointCoord(849, 0) = 3.892169575397059e-03; pointCoord(849, 1) =5.921368160743566e-01; pointCoord(849, 2) =1.004375733945300e-04;
pointCoord(850, 0) = 2.038921695753971e-01; pointCoord(850, 1) =5.921368160743566e-01; pointCoord(850, 2) =1.004375733945300e-04;
pointCoord(851, 0) = 4.038921695753971e-01; pointCoord(851, 1) =5.921368160743566e-01; pointCoord(851, 2) =1.004375733945300e-04;
pointCoord(852, 0) = 3.892169575397059e-03; pointCoord(852, 1) =7.921368160743566e-01; pointCoord(852, 2) =1.004375733945300e-04;
pointCoord(853, 0) = 2.038921695753971e-01; pointCoord(853, 1) =7.921368160743566e-01; pointCoord(853, 2) =1.004375733945300e-04;
pointCoord(854, 0) = 3.892169575397059e-03; pointCoord(854, 1) =9.921368160743567e-01; pointCoord(854, 2) =1.004375733945300e-04;
pointCoord(855, 0) = 1.992963209059901e-02; pointCoord(855, 1) =1.760993535591546e-01; pointCoord(855, 2) =2.206434300837125e-04;
pointCoord(856, 0) = 2.199296320905990e-01; pointCoord(856, 1) =1.760993535591546e-01; pointCoord(856, 2) =2.206434300837125e-04;
pointCoord(857, 0) = 4.199296320905990e-01; pointCoord(857, 1) =1.760993535591546e-01; pointCoord(857, 2) =2.206434300837125e-04;
pointCoord(858, 0) = 6.199296320905990e-01; pointCoord(858, 1) =1.760993535591546e-01; pointCoord(858, 2) =2.206434300837125e-04;
pointCoord(859, 0) = 8.199296320905991e-01; pointCoord(859, 1) =1.760993535591546e-01; pointCoord(859, 2) =2.206434300837125e-04;
pointCoord(860, 0) = 1.992963209059901e-02; pointCoord(860, 1) =3.760993535591546e-01; pointCoord(860, 2) =2.206434300837125e-04;
pointCoord(861, 0) = 2.199296320905990e-01; pointCoord(861, 1) =3.760993535591546e-01; pointCoord(861, 2) =2.206434300837125e-04;
pointCoord(862, 0) = 4.199296320905990e-01; pointCoord(862, 1) =3.760993535591546e-01; pointCoord(862, 2) =2.206434300837125e-04;
pointCoord(863, 0) = 6.199296320905990e-01; pointCoord(863, 1) =3.760993535591546e-01; pointCoord(863, 2) =2.206434300837125e-04;
pointCoord(864, 0) = 1.992963209059901e-02; pointCoord(864, 1) =5.760993535591546e-01; pointCoord(864, 2) =2.206434300837125e-04;
pointCoord(865, 0) = 2.199296320905990e-01; pointCoord(865, 1) =5.760993535591546e-01; pointCoord(865, 2) =2.206434300837125e-04;
pointCoord(866, 0) = 4.199296320905990e-01; pointCoord(866, 1) =5.760993535591546e-01; pointCoord(866, 2) =2.206434300837125e-04;
pointCoord(867, 0) = 1.992963209059901e-02; pointCoord(867, 1) =7.760993535591546e-01; pointCoord(867, 2) =2.206434300837125e-04;
pointCoord(868, 0) = 2.199296320905990e-01; pointCoord(868, 1) =7.760993535591546e-01; pointCoord(868, 2) =2.206434300837125e-04;
pointCoord(869, 0) = 1.992963209059901e-02; pointCoord(869, 1) =9.760993535591547e-01; pointCoord(869, 2) =2.206434300837125e-04;
pointCoord(870, 0) = 4.650470020389257e-02; pointCoord(870, 1) =1.495242854458611e-01; pointCoord(870, 2) =3.112554564587477e-04;
pointCoord(871, 0) = 2.465047002038926e-01; pointCoord(871, 1) =1.495242854458611e-01; pointCoord(871, 2) =3.112554564587477e-04;
pointCoord(872, 0) = 4.465047002038926e-01; pointCoord(872, 1) =1.495242854458611e-01; pointCoord(872, 2) =3.112554564587477e-04;
pointCoord(873, 0) = 6.465047002038925e-01; pointCoord(873, 1) =1.495242854458611e-01; pointCoord(873, 2) =3.112554564587477e-04;
pointCoord(874, 0) = 8.465047002038926e-01; pointCoord(874, 1) =1.495242854458611e-01; pointCoord(874, 2) =3.112554564587477e-04;
pointCoord(875, 0) = 4.650470020389257e-02; pointCoord(875, 1) =3.495242854458611e-01; pointCoord(875, 2) =3.112554564587477e-04;
pointCoord(876, 0) = 2.465047002038926e-01; pointCoord(876, 1) =3.495242854458611e-01; pointCoord(876, 2) =3.112554564587477e-04;
pointCoord(877, 0) = 4.465047002038926e-01; pointCoord(877, 1) =3.495242854458611e-01; pointCoord(877, 2) =3.112554564587477e-04;
pointCoord(878, 0) = 6.465047002038925e-01; pointCoord(878, 1) =3.495242854458611e-01; pointCoord(878, 2) =3.112554564587477e-04;
pointCoord(879, 0) = 4.650470020389257e-02; pointCoord(879, 1) =5.495242854458611e-01; pointCoord(879, 2) =3.112554564587477e-04;
pointCoord(880, 0) = 2.465047002038926e-01; pointCoord(880, 1) =5.495242854458611e-01; pointCoord(880, 2) =3.112554564587477e-04;
pointCoord(881, 0) = 4.465047002038926e-01; pointCoord(881, 1) =5.495242854458611e-01; pointCoord(881, 2) =3.112554564587477e-04;
pointCoord(882, 0) = 4.650470020389257e-02; pointCoord(882, 1) =7.495242854458610e-01; pointCoord(882, 2) =3.112554564587477e-04;
pointCoord(883, 0) = 2.465047002038926e-01; pointCoord(883, 1) =7.495242854458610e-01; pointCoord(883, 2) =3.112554564587477e-04;
pointCoord(884, 0) = 4.650470020389257e-02; pointCoord(884, 1) =9.495242854458611e-01; pointCoord(884, 2) =3.112554564587477e-04;
pointCoord(885, 0) = 8.003523937415311e-02; pointCoord(885, 1) =1.159937462756005e-01; pointCoord(885, 2) =3.598499044536019e-04;
pointCoord(886, 0) = 2.800352393741531e-01; pointCoord(886, 1) =1.159937462756005e-01; pointCoord(886, 2) =3.598499044536019e-04;
pointCoord(887, 0) = 4.800352393741532e-01; pointCoord(887, 1) =1.159937462756005e-01; pointCoord(887, 2) =3.598499044536019e-04;
pointCoord(888, 0) = 6.800352393741531e-01; pointCoord(888, 1) =1.159937462756005e-01; pointCoord(888, 2) =3.598499044536019e-04;
pointCoord(889, 0) = 8.800352393741532e-01; pointCoord(889, 1) =1.159937462756005e-01; pointCoord(889, 2) =3.598499044536019e-04;
pointCoord(890, 0) = 8.003523937415311e-02; pointCoord(890, 1) =3.159937462756005e-01; pointCoord(890, 2) =3.598499044536019e-04;
pointCoord(891, 0) = 2.800352393741531e-01; pointCoord(891, 1) =3.159937462756005e-01; pointCoord(891, 2) =3.598499044536019e-04;
pointCoord(892, 0) = 4.800352393741532e-01; pointCoord(892, 1) =3.159937462756005e-01; pointCoord(892, 2) =3.598499044536019e-04;
pointCoord(893, 0) = 6.800352393741531e-01; pointCoord(893, 1) =3.159937462756005e-01; pointCoord(893, 2) =3.598499044536019e-04;
pointCoord(894, 0) = 8.003523937415311e-02; pointCoord(894, 1) =5.159937462756006e-01; pointCoord(894, 2) =3.598499044536019e-04;
pointCoord(895, 0) = 2.800352393741531e-01; pointCoord(895, 1) =5.159937462756006e-01; pointCoord(895, 2) =3.598499044536019e-04;
pointCoord(896, 0) = 4.800352393741532e-01; pointCoord(896, 1) =5.159937462756006e-01; pointCoord(896, 2) =3.598499044536019e-04;
pointCoord(897, 0) = 8.003523937415311e-02; pointCoord(897, 1) =7.159937462756005e-01; pointCoord(897, 2) =3.598499044536019e-04;
pointCoord(898, 0) = 2.800352393741531e-01; pointCoord(898, 1) =7.159937462756005e-01; pointCoord(898, 2) =3.598499044536019e-04;
pointCoord(899, 0) = 8.003523937415311e-02; pointCoord(899, 1) =9.159937462756006e-01; pointCoord(899, 2) =3.598499044536019e-04;
pointCoord(900, 0) = 1.159937462756005e-01; pointCoord(900, 1) =8.003523937415311e-02; pointCoord(900, 2) =3.598499044536019e-04;
pointCoord(901, 0) = 3.159937462756005e-01; pointCoord(901, 1) =8.003523937415311e-02; pointCoord(901, 2) =3.598499044536019e-04;
pointCoord(902, 0) = 5.159937462756006e-01; pointCoord(902, 1) =8.003523937415311e-02; pointCoord(902, 2) =3.598499044536019e-04;
pointCoord(903, 0) = 7.159937462756005e-01; pointCoord(903, 1) =8.003523937415311e-02; pointCoord(903, 2) =3.598499044536019e-04;
pointCoord(904, 0) = 9.159937462756006e-01; pointCoord(904, 1) =8.003523937415311e-02; pointCoord(904, 2) =3.598499044536019e-04;
pointCoord(905, 0) = 1.159937462756005e-01; pointCoord(905, 1) =2.800352393741531e-01; pointCoord(905, 2) =3.598499044536019e-04;
pointCoord(906, 0) = 3.159937462756005e-01; pointCoord(906, 1) =2.800352393741531e-01; pointCoord(906, 2) =3.598499044536019e-04;
pointCoord(907, 0) = 5.159937462756006e-01; pointCoord(907, 1) =2.800352393741531e-01; pointCoord(907, 2) =3.598499044536019e-04;
pointCoord(908, 0) = 7.159937462756005e-01; pointCoord(908, 1) =2.800352393741531e-01; pointCoord(908, 2) =3.598499044536019e-04;
pointCoord(909, 0) = 1.159937462756005e-01; pointCoord(909, 1) =4.800352393741532e-01; pointCoord(909, 2) =3.598499044536019e-04;
pointCoord(910, 0) = 3.159937462756005e-01; pointCoord(910, 1) =4.800352393741532e-01; pointCoord(910, 2) =3.598499044536019e-04;
pointCoord(911, 0) = 5.159937462756006e-01; pointCoord(911, 1) =4.800352393741532e-01; pointCoord(911, 2) =3.598499044536019e-04;
pointCoord(912, 0) = 1.159937462756005e-01; pointCoord(912, 1) =6.800352393741531e-01; pointCoord(912, 2) =3.598499044536019e-04;
pointCoord(913, 0) = 3.159937462756005e-01; pointCoord(913, 1) =6.800352393741531e-01; pointCoord(913, 2) =3.598499044536019e-04;
pointCoord(914, 0) = 1.159937462756005e-01; pointCoord(914, 1) =8.800352393741532e-01; pointCoord(914, 2) =3.598499044536019e-04;
pointCoord(915, 0) = 1.495242854458611e-01; pointCoord(915, 1) =4.650470020389257e-02; pointCoord(915, 2) =3.112554564587477e-04;
pointCoord(916, 0) = 3.495242854458611e-01; pointCoord(916, 1) =4.650470020389257e-02; pointCoord(916, 2) =3.112554564587477e-04;
pointCoord(917, 0) = 5.495242854458611e-01; pointCoord(917, 1) =4.650470020389257e-02; pointCoord(917, 2) =3.112554564587477e-04;
pointCoord(918, 0) = 7.495242854458610e-01; pointCoord(918, 1) =4.650470020389257e-02; pointCoord(918, 2) =3.112554564587477e-04;
pointCoord(919, 0) = 9.495242854458611e-01; pointCoord(919, 1) =4.650470020389257e-02; pointCoord(919, 2) =3.112554564587477e-04;
pointCoord(920, 0) = 1.495242854458611e-01; pointCoord(920, 1) =2.465047002038926e-01; pointCoord(920, 2) =3.112554564587477e-04;
pointCoord(921, 0) = 3.495242854458611e-01; pointCoord(921, 1) =2.465047002038926e-01; pointCoord(921, 2) =3.112554564587477e-04;
pointCoord(922, 0) = 5.495242854458611e-01; pointCoord(922, 1) =2.465047002038926e-01; pointCoord(922, 2) =3.112554564587477e-04;
pointCoord(923, 0) = 7.495242854458610e-01; pointCoord(923, 1) =2.465047002038926e-01; pointCoord(923, 2) =3.112554564587477e-04;
pointCoord(924, 0) = 1.495242854458611e-01; pointCoord(924, 1) =4.465047002038926e-01; pointCoord(924, 2) =3.112554564587477e-04;
pointCoord(925, 0) = 3.495242854458611e-01; pointCoord(925, 1) =4.465047002038926e-01; pointCoord(925, 2) =3.112554564587477e-04;
pointCoord(926, 0) = 5.495242854458611e-01; pointCoord(926, 1) =4.465047002038926e-01; pointCoord(926, 2) =3.112554564587477e-04;
pointCoord(927, 0) = 1.495242854458611e-01; pointCoord(927, 1) =6.465047002038925e-01; pointCoord(927, 2) =3.112554564587477e-04;
pointCoord(928, 0) = 3.495242854458611e-01; pointCoord(928, 1) =6.465047002038925e-01; pointCoord(928, 2) =3.112554564587477e-04;
pointCoord(929, 0) = 1.495242854458611e-01; pointCoord(929, 1) =8.465047002038926e-01; pointCoord(929, 2) =3.112554564587477e-04;
pointCoord(930, 0) = 1.760993535591546e-01; pointCoord(930, 1) =1.992963209059901e-02; pointCoord(930, 2) =2.206434300837125e-04;
pointCoord(931, 0) = 3.760993535591546e-01; pointCoord(931, 1) =1.992963209059901e-02; pointCoord(931, 2) =2.206434300837125e-04;
pointCoord(932, 0) = 5.760993535591546e-01; pointCoord(932, 1) =1.992963209059901e-02; pointCoord(932, 2) =2.206434300837125e-04;
pointCoord(933, 0) = 7.760993535591546e-01; pointCoord(933, 1) =1.992963209059901e-02; pointCoord(933, 2) =2.206434300837125e-04;
pointCoord(934, 0) = 9.760993535591547e-01; pointCoord(934, 1) =1.992963209059901e-02; pointCoord(934, 2) =2.206434300837125e-04;
pointCoord(935, 0) = 1.760993535591546e-01; pointCoord(935, 1) =2.199296320905990e-01; pointCoord(935, 2) =2.206434300837125e-04;
pointCoord(936, 0) = 3.760993535591546e-01; pointCoord(936, 1) =2.199296320905990e-01; pointCoord(936, 2) =2.206434300837125e-04;
pointCoord(937, 0) = 5.760993535591546e-01; pointCoord(937, 1) =2.199296320905990e-01; pointCoord(937, 2) =2.206434300837125e-04;
pointCoord(938, 0) = 7.760993535591546e-01; pointCoord(938, 1) =2.199296320905990e-01; pointCoord(938, 2) =2.206434300837125e-04;
pointCoord(939, 0) = 1.760993535591546e-01; pointCoord(939, 1) =4.199296320905990e-01; pointCoord(939, 2) =2.206434300837125e-04;
pointCoord(940, 0) = 3.760993535591546e-01; pointCoord(940, 1) =4.199296320905990e-01; pointCoord(940, 2) =2.206434300837125e-04;
pointCoord(941, 0) = 5.760993535591546e-01; pointCoord(941, 1) =4.199296320905990e-01; pointCoord(941, 2) =2.206434300837125e-04;
pointCoord(942, 0) = 1.760993535591546e-01; pointCoord(942, 1) =6.199296320905990e-01; pointCoord(942, 2) =2.206434300837125e-04;
pointCoord(943, 0) = 3.760993535591546e-01; pointCoord(943, 1) =6.199296320905990e-01; pointCoord(943, 2) =2.206434300837125e-04;
pointCoord(944, 0) = 1.760993535591546e-01; pointCoord(944, 1) =8.199296320905991e-01; pointCoord(944, 2) =2.206434300837125e-04;
pointCoord(945, 0) = 1.921368160743566e-01; pointCoord(945, 1) =3.892169575397059e-03; pointCoord(945, 2) =1.004375733945300e-04;
pointCoord(946, 0) = 3.921368160743566e-01; pointCoord(946, 1) =3.892169575397059e-03; pointCoord(946, 2) =1.004375733945300e-04;
pointCoord(947, 0) = 5.921368160743566e-01; pointCoord(947, 1) =3.892169575397059e-03; pointCoord(947, 2) =1.004375733945300e-04;
pointCoord(948, 0) = 7.921368160743566e-01; pointCoord(948, 1) =3.892169575397059e-03; pointCoord(948, 2) =1.004375733945300e-04;
pointCoord(949, 0) = 9.921368160743567e-01; pointCoord(949, 1) =3.892169575397059e-03; pointCoord(949, 2) =1.004375733945300e-04;
pointCoord(950, 0) = 1.921368160743566e-01; pointCoord(950, 1) =2.038921695753971e-01; pointCoord(950, 2) =1.004375733945300e-04;
pointCoord(951, 0) = 3.921368160743566e-01; pointCoord(951, 1) =2.038921695753971e-01; pointCoord(951, 2) =1.004375733945300e-04;
pointCoord(952, 0) = 5.921368160743566e-01; pointCoord(952, 1) =2.038921695753971e-01; pointCoord(952, 2) =1.004375733945300e-04;
pointCoord(953, 0) = 7.921368160743566e-01; pointCoord(953, 1) =2.038921695753971e-01; pointCoord(953, 2) =1.004375733945300e-04;
pointCoord(954, 0) = 1.921368160743566e-01; pointCoord(954, 1) =4.038921695753971e-01; pointCoord(954, 2) =1.004375733945300e-04;
pointCoord(955, 0) = 3.921368160743566e-01; pointCoord(955, 1) =4.038921695753971e-01; pointCoord(955, 2) =1.004375733945300e-04;
pointCoord(956, 0) = 5.921368160743566e-01; pointCoord(956, 1) =4.038921695753971e-01; pointCoord(956, 2) =1.004375733945300e-04;
pointCoord(957, 0) = 1.921368160743566e-01; pointCoord(957, 1) =6.038921695753970e-01; pointCoord(957, 2) =1.004375733945300e-04;
pointCoord(958, 0) = 3.921368160743566e-01; pointCoord(958, 1) =6.038921695753970e-01; pointCoord(958, 2) =1.004375733945300e-04;
pointCoord(959, 0) = 1.921368160743566e-01; pointCoord(959, 1) =8.038921695753971e-01; pointCoord(959, 2) =1.004375733945300e-04;
pointCoord(960, 0) = 1.999211552251507e-01; pointCoord(960, 1) =1.961078304246029e-01; pointCoord(960, 2) =2.034592200391275e-06;
pointCoord(961, 0) = 3.999211552251507e-01; pointCoord(961, 1) =1.961078304246029e-01; pointCoord(961, 2) =2.034592200391275e-06;
pointCoord(962, 0) = 5.999211552251507e-01; pointCoord(962, 1) =1.961078304246029e-01; pointCoord(962, 2) =2.034592200391275e-06;
pointCoord(963, 0) = 7.999211552251507e-01; pointCoord(963, 1) =1.961078304246029e-01; pointCoord(963, 2) =2.034592200391275e-06;
pointCoord(964, 0) = 1.999211552251507e-01; pointCoord(964, 1) =3.961078304246030e-01; pointCoord(964, 2) =2.034592200391275e-06;
pointCoord(965, 0) = 3.999211552251507e-01; pointCoord(965, 1) =3.961078304246030e-01; pointCoord(965, 2) =2.034592200391275e-06;
pointCoord(966, 0) = 5.999211552251507e-01; pointCoord(966, 1) =3.961078304246030e-01; pointCoord(966, 2) =2.034592200391275e-06;
pointCoord(967, 0) = 1.999211552251507e-01; pointCoord(967, 1) =5.961078304246029e-01; pointCoord(967, 2) =2.034592200391275e-06;
pointCoord(968, 0) = 3.999211552251507e-01; pointCoord(968, 1) =5.961078304246029e-01; pointCoord(968, 2) =2.034592200391275e-06;
pointCoord(969, 0) = 1.999211552251507e-01; pointCoord(969, 1) =7.961078304246030e-01; pointCoord(969, 2) =2.034592200391275e-06;
pointCoord(970, 0) = 1.995962798319617e-01; pointCoord(970, 1) =1.964327058177920e-01; pointCoord(970, 2) =4.469636080836934e-06;
pointCoord(971, 0) = 3.995962798319617e-01; pointCoord(971, 1) =1.964327058177920e-01; pointCoord(971, 2) =4.469636080836934e-06;
pointCoord(972, 0) = 5.995962798319616e-01; pointCoord(972, 1) =1.964327058177920e-01; pointCoord(972, 2) =4.469636080836934e-06;
pointCoord(973, 0) = 7.995962798319617e-01; pointCoord(973, 1) =1.964327058177920e-01; pointCoord(973, 2) =4.469636080836934e-06;
pointCoord(974, 0) = 1.995962798319617e-01; pointCoord(974, 1) =3.964327058177919e-01; pointCoord(974, 2) =4.469636080836934e-06;
pointCoord(975, 0) = 3.995962798319617e-01; pointCoord(975, 1) =3.964327058177919e-01; pointCoord(975, 2) =4.469636080836934e-06;
pointCoord(976, 0) = 5.995962798319616e-01; pointCoord(976, 1) =3.964327058177919e-01; pointCoord(976, 2) =4.469636080836934e-06;
pointCoord(977, 0) = 1.995962798319617e-01; pointCoord(977, 1) =5.964327058177920e-01; pointCoord(977, 2) =4.469636080836934e-06;
pointCoord(978, 0) = 3.995962798319617e-01; pointCoord(978, 1) =5.964327058177920e-01; pointCoord(978, 2) =4.469636080836934e-06;
pointCoord(979, 0) = 1.995962798319617e-01; pointCoord(979, 1) =7.964327058177920e-01; pointCoord(979, 2) =4.469636080836934e-06;
pointCoord(980, 0) = 1.990579411955255e-01; pointCoord(980, 1) =1.969710444542282e-01; pointCoord(980, 2) =6.305189409073112e-06;
pointCoord(981, 0) = 3.990579411955255e-01; pointCoord(981, 1) =1.969710444542282e-01; pointCoord(981, 2) =6.305189409073112e-06;
pointCoord(982, 0) = 5.990579411955255e-01; pointCoord(982, 1) =1.969710444542282e-01; pointCoord(982, 2) =6.305189409073112e-06;
pointCoord(983, 0) = 7.990579411955255e-01; pointCoord(983, 1) =1.969710444542282e-01; pointCoord(983, 2) =6.305189409073112e-06;
pointCoord(984, 0) = 1.990579411955255e-01; pointCoord(984, 1) =3.969710444542282e-01; pointCoord(984, 2) =6.305189409073112e-06;
pointCoord(985, 0) = 3.990579411955255e-01; pointCoord(985, 1) =3.969710444542282e-01; pointCoord(985, 2) =6.305189409073112e-06;
pointCoord(986, 0) = 5.990579411955255e-01; pointCoord(986, 1) =3.969710444542282e-01; pointCoord(986, 2) =6.305189409073112e-06;
pointCoord(987, 0) = 1.990579411955255e-01; pointCoord(987, 1) =5.969710444542281e-01; pointCoord(987, 2) =6.305189409073112e-06;
pointCoord(988, 0) = 3.990579411955255e-01; pointCoord(988, 1) =5.969710444542281e-01; pointCoord(988, 2) =6.305189409073112e-06;
pointCoord(989, 0) = 1.990579411955255e-01; pointCoord(989, 1) =7.969710444542282e-01; pointCoord(989, 2) =6.305189409073112e-06;
pointCoord(990, 0) = 1.983787036237181e-01; pointCoord(990, 1) =1.976502820260355e-01; pointCoord(990, 2) =7.289580822874776e-06;
pointCoord(991, 0) = 3.983787036237181e-01; pointCoord(991, 1) =1.976502820260355e-01; pointCoord(991, 2) =7.289580822874776e-06;
pointCoord(992, 0) = 5.983787036237180e-01; pointCoord(992, 1) =1.976502820260355e-01; pointCoord(992, 2) =7.289580822874776e-06;
pointCoord(993, 0) = 7.983787036237181e-01; pointCoord(993, 1) =1.976502820260355e-01; pointCoord(993, 2) =7.289580822874776e-06;
pointCoord(994, 0) = 1.983787036237181e-01; pointCoord(994, 1) =3.976502820260356e-01; pointCoord(994, 2) =7.289580822874776e-06;
pointCoord(995, 0) = 3.983787036237181e-01; pointCoord(995, 1) =3.976502820260356e-01; pointCoord(995, 2) =7.289580822874776e-06;
pointCoord(996, 0) = 5.983787036237180e-01; pointCoord(996, 1) =3.976502820260356e-01; pointCoord(996, 2) =7.289580822874776e-06;
pointCoord(997, 0) = 1.983787036237181e-01; pointCoord(997, 1) =5.976502820260355e-01; pointCoord(997, 2) =7.289580822874776e-06;
pointCoord(998, 0) = 3.983787036237181e-01; pointCoord(998, 1) =5.976502820260355e-01; pointCoord(998, 2) =7.289580822874776e-06;
pointCoord(999, 0) = 1.983787036237181e-01; pointCoord(999, 1) =7.976502820260356e-01; pointCoord(999, 2) =7.289580822874776e-06;
pointCoord(1000, 0) = 1.976502820260355e-01; pointCoord(1000, 1) =1.983787036237181e-01; pointCoord(1000, 2) =7.289580822874776e-06;
pointCoord(1001, 0) = 3.976502820260356e-01; pointCoord(1001, 1) =1.983787036237181e-01; pointCoord(1001, 2) =7.289580822874776e-06;
pointCoord(1002, 0) = 5.976502820260355e-01; pointCoord(1002, 1) =1.983787036237181e-01; pointCoord(1002, 2) =7.289580822874776e-06;
pointCoord(1003, 0) = 7.976502820260356e-01; pointCoord(1003, 1) =1.983787036237181e-01; pointCoord(1003, 2) =7.289580822874776e-06;
pointCoord(1004, 0) = 1.976502820260355e-01; pointCoord(1004, 1) =3.983787036237181e-01; pointCoord(1004, 2) =7.289580822874776e-06;
pointCoord(1005, 0) = 3.976502820260356e-01; pointCoord(1005, 1) =3.983787036237181e-01; pointCoord(1005, 2) =7.289580822874776e-06;
pointCoord(1006, 0) = 5.976502820260355e-01; pointCoord(1006, 1) =3.983787036237181e-01; pointCoord(1006, 2) =7.289580822874776e-06;
pointCoord(1007, 0) = 1.976502820260355e-01; pointCoord(1007, 1) =5.983787036237180e-01; pointCoord(1007, 2) =7.289580822874776e-06;
pointCoord(1008, 0) = 3.976502820260356e-01; pointCoord(1008, 1) =5.983787036237180e-01; pointCoord(1008, 2) =7.289580822874776e-06;
pointCoord(1009, 0) = 1.976502820260355e-01; pointCoord(1009, 1) =7.983787036237181e-01; pointCoord(1009, 2) =7.289580822874776e-06;
pointCoord(1010, 0) = 1.969710444542282e-01; pointCoord(1010, 1) =1.990579411955255e-01; pointCoord(1010, 2) =6.305189409073112e-06;
pointCoord(1011, 0) = 3.969710444542282e-01; pointCoord(1011, 1) =1.990579411955255e-01; pointCoord(1011, 2) =6.305189409073112e-06;
pointCoord(1012, 0) = 5.969710444542281e-01; pointCoord(1012, 1) =1.990579411955255e-01; pointCoord(1012, 2) =6.305189409073112e-06;
pointCoord(1013, 0) = 7.969710444542282e-01; pointCoord(1013, 1) =1.990579411955255e-01; pointCoord(1013, 2) =6.305189409073112e-06;
pointCoord(1014, 0) = 1.969710444542282e-01; pointCoord(1014, 1) =3.990579411955255e-01; pointCoord(1014, 2) =6.305189409073112e-06;
pointCoord(1015, 0) = 3.969710444542282e-01; pointCoord(1015, 1) =3.990579411955255e-01; pointCoord(1015, 2) =6.305189409073112e-06;
pointCoord(1016, 0) = 5.969710444542281e-01; pointCoord(1016, 1) =3.990579411955255e-01; pointCoord(1016, 2) =6.305189409073112e-06;
pointCoord(1017, 0) = 1.969710444542282e-01; pointCoord(1017, 1) =5.990579411955255e-01; pointCoord(1017, 2) =6.305189409073112e-06;
pointCoord(1018, 0) = 3.969710444542282e-01; pointCoord(1018, 1) =5.990579411955255e-01; pointCoord(1018, 2) =6.305189409073112e-06;
pointCoord(1019, 0) = 1.969710444542282e-01; pointCoord(1019, 1) =7.990579411955255e-01; pointCoord(1019, 2) =6.305189409073112e-06;
pointCoord(1020, 0) = 1.964327058177920e-01; pointCoord(1020, 1) =1.995962798319617e-01; pointCoord(1020, 2) =4.469636080836934e-06;
pointCoord(1021, 0) = 3.964327058177919e-01; pointCoord(1021, 1) =1.995962798319617e-01; pointCoord(1021, 2) =4.469636080836934e-06;
pointCoord(1022, 0) = 5.964327058177920e-01; pointCoord(1022, 1) =1.995962798319617e-01; pointCoord(1022, 2) =4.469636080836934e-06;
pointCoord(1023, 0) = 7.964327058177920e-01; pointCoord(1023, 1) =1.995962798319617e-01; pointCoord(1023, 2) =4.469636080836934e-06;
pointCoord(1024, 0) = 1.964327058177920e-01; pointCoord(1024, 1) =3.995962798319617e-01; pointCoord(1024, 2) =4.469636080836934e-06;
pointCoord(1025, 0) = 3.964327058177919e-01; pointCoord(1025, 1) =3.995962798319617e-01; pointCoord(1025, 2) =4.469636080836934e-06;
pointCoord(1026, 0) = 5.964327058177920e-01; pointCoord(1026, 1) =3.995962798319617e-01; pointCoord(1026, 2) =4.469636080836934e-06;
pointCoord(1027, 0) = 1.964327058177920e-01; pointCoord(1027, 1) =5.995962798319616e-01; pointCoord(1027, 2) =4.469636080836934e-06;
pointCoord(1028, 0) = 3.964327058177919e-01; pointCoord(1028, 1) =5.995962798319616e-01; pointCoord(1028, 2) =4.469636080836934e-06;
pointCoord(1029, 0) = 1.964327058177920e-01; pointCoord(1029, 1) =7.995962798319617e-01; pointCoord(1029, 2) =4.469636080836934e-06;
pointCoord(1030, 0) = 1.961078304246029e-01; pointCoord(1030, 1) =1.999211552251507e-01; pointCoord(1030, 2) =2.034592200391275e-06;
pointCoord(1031, 0) = 3.961078304246030e-01; pointCoord(1031, 1) =1.999211552251507e-01; pointCoord(1031, 2) =2.034592200391275e-06;
pointCoord(1032, 0) = 5.961078304246029e-01; pointCoord(1032, 1) =1.999211552251507e-01; pointCoord(1032, 2) =2.034592200391275e-06;
pointCoord(1033, 0) = 7.961078304246030e-01; pointCoord(1033, 1) =1.999211552251507e-01; pointCoord(1033, 2) =2.034592200391275e-06;
pointCoord(1034, 0) = 1.961078304246029e-01; pointCoord(1034, 1) =3.999211552251507e-01; pointCoord(1034, 2) =2.034592200391275e-06;
pointCoord(1035, 0) = 3.961078304246030e-01; pointCoord(1035, 1) =3.999211552251507e-01; pointCoord(1035, 2) =2.034592200391275e-06;
pointCoord(1036, 0) = 5.961078304246029e-01; pointCoord(1036, 1) =3.999211552251507e-01; pointCoord(1036, 2) =2.034592200391275e-06;
pointCoord(1037, 0) = 1.961078304246029e-01; pointCoord(1037, 1) =5.999211552251507e-01; pointCoord(1037, 2) =2.034592200391275e-06;
pointCoord(1038, 0) = 3.961078304246030e-01; pointCoord(1038, 1) =5.999211552251507e-01; pointCoord(1038, 2) =2.034592200391275e-06;
pointCoord(1039, 0) = 1.961078304246029e-01; pointCoord(1039, 1) =7.999211552251507e-01; pointCoord(1039, 2) =2.034592200391275e-06;
pointCoord(1040, 0) = 1.995962798319617e-01; pointCoord(1040, 1) =1.800703679094010e-01; pointCoord(1040, 2) =2.288651636172858e-05;
pointCoord(1041, 0) = 3.995962798319617e-01; pointCoord(1041, 1) =1.800703679094010e-01; pointCoord(1041, 2) =2.288651636172858e-05;
pointCoord(1042, 0) = 5.995962798319616e-01; pointCoord(1042, 1) =1.800703679094010e-01; pointCoord(1042, 2) =2.288651636172858e-05;
pointCoord(1043, 0) = 7.995962798319617e-01; pointCoord(1043, 1) =1.800703679094010e-01; pointCoord(1043, 2) =2.288651636172858e-05;
pointCoord(1044, 0) = 1.995962798319617e-01; pointCoord(1044, 1) =3.800703679094010e-01; pointCoord(1044, 2) =2.288651636172858e-05;
pointCoord(1045, 0) = 3.995962798319617e-01; pointCoord(1045, 1) =3.800703679094010e-01; pointCoord(1045, 2) =2.288651636172858e-05;
pointCoord(1046, 0) = 5.995962798319616e-01; pointCoord(1046, 1) =3.800703679094010e-01; pointCoord(1046, 2) =2.288651636172858e-05;
pointCoord(1047, 0) = 1.995962798319617e-01; pointCoord(1047, 1) =5.800703679094009e-01; pointCoord(1047, 2) =2.288651636172858e-05;
pointCoord(1048, 0) = 3.995962798319617e-01; pointCoord(1048, 1) =5.800703679094009e-01; pointCoord(1048, 2) =2.288651636172858e-05;
pointCoord(1049, 0) = 1.995962798319617e-01; pointCoord(1049, 1) =7.800703679094010e-01; pointCoord(1049, 2) =2.288651636172858e-05;
pointCoord(1050, 0) = 1.979327739296309e-01; pointCoord(1050, 1) =1.817338738117318e-01; pointCoord(1050, 2) =5.027759335525545e-05;
pointCoord(1051, 0) = 3.979327739296309e-01; pointCoord(1051, 1) =1.817338738117318e-01; pointCoord(1051, 2) =5.027759335525545e-05;
pointCoord(1052, 0) = 5.979327739296308e-01; pointCoord(1052, 1) =1.817338738117318e-01; pointCoord(1052, 2) =5.027759335525545e-05;
pointCoord(1053, 0) = 7.979327739296309e-01; pointCoord(1053, 1) =1.817338738117318e-01; pointCoord(1053, 2) =5.027759335525545e-05;
pointCoord(1054, 0) = 1.979327739296309e-01; pointCoord(1054, 1) =3.817338738117318e-01; pointCoord(1054, 2) =5.027759335525545e-05;
pointCoord(1055, 0) = 3.979327739296309e-01; pointCoord(1055, 1) =3.817338738117318e-01; pointCoord(1055, 2) =5.027759335525545e-05;
pointCoord(1056, 0) = 5.979327739296308e-01; pointCoord(1056, 1) =3.817338738117318e-01; pointCoord(1056, 2) =5.027759335525545e-05;
pointCoord(1057, 0) = 1.979327739296309e-01; pointCoord(1057, 1) =5.817338738117318e-01; pointCoord(1057, 2) =5.027759335525545e-05;
pointCoord(1058, 0) = 3.979327739296309e-01; pointCoord(1058, 1) =5.817338738117318e-01; pointCoord(1058, 2) =5.027759335525545e-05;
pointCoord(1059, 0) = 1.979327739296309e-01; pointCoord(1059, 1) =7.817338738117319e-01; pointCoord(1059, 2) =5.027759335525545e-05;
pointCoord(1060, 0) = 1.951762416777610e-01; pointCoord(1060, 1) =1.844904060636017e-01; pointCoord(1060, 2) =7.092518124604937e-05;
pointCoord(1061, 0) = 3.951762416777610e-01; pointCoord(1061, 1) =1.844904060636017e-01; pointCoord(1061, 2) =7.092518124604937e-05;
pointCoord(1062, 0) = 5.951762416777610e-01; pointCoord(1062, 1) =1.844904060636017e-01; pointCoord(1062, 2) =7.092518124604937e-05;
pointCoord(1063, 0) = 7.951762416777610e-01; pointCoord(1063, 1) =1.844904060636017e-01; pointCoord(1063, 2) =7.092518124604937e-05;
pointCoord(1064, 0) = 1.951762416777610e-01; pointCoord(1064, 1) =3.844904060636017e-01; pointCoord(1064, 2) =7.092518124604937e-05;
pointCoord(1065, 0) = 3.951762416777610e-01; pointCoord(1065, 1) =3.844904060636017e-01; pointCoord(1065, 2) =7.092518124604937e-05;
pointCoord(1066, 0) = 5.951762416777610e-01; pointCoord(1066, 1) =3.844904060636017e-01; pointCoord(1066, 2) =7.092518124604937e-05;
pointCoord(1067, 0) = 1.951762416777610e-01; pointCoord(1067, 1) =5.844904060636017e-01; pointCoord(1067, 2) =7.092518124604937e-05;
pointCoord(1068, 0) = 3.951762416777610e-01; pointCoord(1068, 1) =5.844904060636017e-01; pointCoord(1068, 2) =7.092518124604937e-05;
pointCoord(1069, 0) = 1.951762416777610e-01; pointCoord(1069, 1) =7.844904060636018e-01; pointCoord(1069, 2) =7.092518124604937e-05;
pointCoord(1070, 0) = 1.916982444718320e-01; pointCoord(1070, 1) =1.879684032695307e-01; pointCoord(1070, 2) =8.199830449599807e-05;
pointCoord(1071, 0) = 3.916982444718320e-01; pointCoord(1071, 1) =1.879684032695307e-01; pointCoord(1071, 2) =8.199830449599807e-05;
pointCoord(1072, 0) = 5.916982444718319e-01; pointCoord(1072, 1) =1.879684032695307e-01; pointCoord(1072, 2) =8.199830449599807e-05;
pointCoord(1073, 0) = 7.916982444718320e-01; pointCoord(1073, 1) =1.879684032695307e-01; pointCoord(1073, 2) =8.199830449599807e-05;
pointCoord(1074, 0) = 1.916982444718320e-01; pointCoord(1074, 1) =3.879684032695307e-01; pointCoord(1074, 2) =8.199830449599807e-05;
pointCoord(1075, 0) = 3.916982444718320e-01; pointCoord(1075, 1) =3.879684032695307e-01; pointCoord(1075, 2) =8.199830449599807e-05;
pointCoord(1076, 0) = 5.916982444718319e-01; pointCoord(1076, 1) =3.879684032695307e-01; pointCoord(1076, 2) =8.199830449599807e-05;
pointCoord(1077, 0) = 1.916982444718320e-01; pointCoord(1077, 1) =5.879684032695307e-01; pointCoord(1077, 2) =8.199830449599807e-05;
pointCoord(1078, 0) = 3.916982444718320e-01; pointCoord(1078, 1) =5.879684032695307e-01; pointCoord(1078, 2) =8.199830449599807e-05;
pointCoord(1079, 0) = 1.916982444718320e-01; pointCoord(1079, 1) =7.879684032695308e-01; pointCoord(1079, 2) =8.199830449599807e-05;
pointCoord(1080, 0) = 1.879684032695307e-01; pointCoord(1080, 1) =1.916982444718320e-01; pointCoord(1080, 2) =8.199830449599807e-05;
pointCoord(1081, 0) = 3.879684032695307e-01; pointCoord(1081, 1) =1.916982444718320e-01; pointCoord(1081, 2) =8.199830449599807e-05;
pointCoord(1082, 0) = 5.879684032695307e-01; pointCoord(1082, 1) =1.916982444718320e-01; pointCoord(1082, 2) =8.199830449599807e-05;
pointCoord(1083, 0) = 7.879684032695308e-01; pointCoord(1083, 1) =1.916982444718320e-01; pointCoord(1083, 2) =8.199830449599807e-05;
pointCoord(1084, 0) = 1.879684032695307e-01; pointCoord(1084, 1) =3.916982444718320e-01; pointCoord(1084, 2) =8.199830449599807e-05;
pointCoord(1085, 0) = 3.879684032695307e-01; pointCoord(1085, 1) =3.916982444718320e-01; pointCoord(1085, 2) =8.199830449599807e-05;
pointCoord(1086, 0) = 5.879684032695307e-01; pointCoord(1086, 1) =3.916982444718320e-01; pointCoord(1086, 2) =8.199830449599807e-05;
pointCoord(1087, 0) = 1.879684032695307e-01; pointCoord(1087, 1) =5.916982444718319e-01; pointCoord(1087, 2) =8.199830449599807e-05;
pointCoord(1088, 0) = 3.879684032695307e-01; pointCoord(1088, 1) =5.916982444718319e-01; pointCoord(1088, 2) =8.199830449599807e-05;
pointCoord(1089, 0) = 1.879684032695307e-01; pointCoord(1089, 1) =7.916982444718320e-01; pointCoord(1089, 2) =8.199830449599807e-05;
pointCoord(1090, 0) = 1.844904060636017e-01; pointCoord(1090, 1) =1.951762416777610e-01; pointCoord(1090, 2) =7.092518124604937e-05;
pointCoord(1091, 0) = 3.844904060636017e-01; pointCoord(1091, 1) =1.951762416777610e-01; pointCoord(1091, 2) =7.092518124604937e-05;
pointCoord(1092, 0) = 5.844904060636017e-01; pointCoord(1092, 1) =1.951762416777610e-01; pointCoord(1092, 2) =7.092518124604937e-05;
pointCoord(1093, 0) = 7.844904060636018e-01; pointCoord(1093, 1) =1.951762416777610e-01; pointCoord(1093, 2) =7.092518124604937e-05;
pointCoord(1094, 0) = 1.844904060636017e-01; pointCoord(1094, 1) =3.951762416777610e-01; pointCoord(1094, 2) =7.092518124604937e-05;
pointCoord(1095, 0) = 3.844904060636017e-01; pointCoord(1095, 1) =3.951762416777610e-01; pointCoord(1095, 2) =7.092518124604937e-05;
pointCoord(1096, 0) = 5.844904060636017e-01; pointCoord(1096, 1) =3.951762416777610e-01; pointCoord(1096, 2) =7.092518124604937e-05;
pointCoord(1097, 0) = 1.844904060636017e-01; pointCoord(1097, 1) =5.951762416777610e-01; pointCoord(1097, 2) =7.092518124604937e-05;
pointCoord(1098, 0) = 3.844904060636017e-01; pointCoord(1098, 1) =5.951762416777610e-01; pointCoord(1098, 2) =7.092518124604937e-05;
pointCoord(1099, 0) = 1.844904060636017e-01; pointCoord(1099, 1) =7.951762416777610e-01; pointCoord(1099, 2) =7.092518124604937e-05;
pointCoord(1100, 0) = 1.817338738117318e-01; pointCoord(1100, 1) =1.979327739296309e-01; pointCoord(1100, 2) =5.027759335525545e-05;
pointCoord(1101, 0) = 3.817338738117318e-01; pointCoord(1101, 1) =1.979327739296309e-01; pointCoord(1101, 2) =5.027759335525545e-05;
pointCoord(1102, 0) = 5.817338738117318e-01; pointCoord(1102, 1) =1.979327739296309e-01; pointCoord(1102, 2) =5.027759335525545e-05;
pointCoord(1103, 0) = 7.817338738117319e-01; pointCoord(1103, 1) =1.979327739296309e-01; pointCoord(1103, 2) =5.027759335525545e-05;
pointCoord(1104, 0) = 1.817338738117318e-01; pointCoord(1104, 1) =3.979327739296309e-01; pointCoord(1104, 2) =5.027759335525545e-05;
pointCoord(1105, 0) = 3.817338738117318e-01; pointCoord(1105, 1) =3.979327739296309e-01; pointCoord(1105, 2) =5.027759335525545e-05;
pointCoord(1106, 0) = 5.817338738117318e-01; pointCoord(1106, 1) =3.979327739296309e-01; pointCoord(1106, 2) =5.027759335525545e-05;
pointCoord(1107, 0) = 1.817338738117318e-01; pointCoord(1107, 1) =5.979327739296308e-01; pointCoord(1107, 2) =5.027759335525545e-05;
pointCoord(1108, 0) = 3.817338738117318e-01; pointCoord(1108, 1) =5.979327739296308e-01; pointCoord(1108, 2) =5.027759335525545e-05;
pointCoord(1109, 0) = 1.817338738117318e-01; pointCoord(1109, 1) =7.979327739296309e-01; pointCoord(1109, 2) =5.027759335525545e-05;
pointCoord(1110, 0) = 1.800703679094010e-01; pointCoord(1110, 1) =1.995962798319617e-01; pointCoord(1110, 2) =2.288651636172858e-05;
pointCoord(1111, 0) = 3.800703679094010e-01; pointCoord(1111, 1) =1.995962798319617e-01; pointCoord(1111, 2) =2.288651636172858e-05;
pointCoord(1112, 0) = 5.800703679094009e-01; pointCoord(1112, 1) =1.995962798319617e-01; pointCoord(1112, 2) =2.288651636172858e-05;
pointCoord(1113, 0) = 7.800703679094010e-01; pointCoord(1113, 1) =1.995962798319617e-01; pointCoord(1113, 2) =2.288651636172858e-05;
pointCoord(1114, 0) = 1.800703679094010e-01; pointCoord(1114, 1) =3.995962798319617e-01; pointCoord(1114, 2) =2.288651636172858e-05;
pointCoord(1115, 0) = 3.800703679094010e-01; pointCoord(1115, 1) =3.995962798319617e-01; pointCoord(1115, 2) =2.288651636172858e-05;
pointCoord(1116, 0) = 5.800703679094009e-01; pointCoord(1116, 1) =3.995962798319617e-01; pointCoord(1116, 2) =2.288651636172858e-05;
pointCoord(1117, 0) = 1.800703679094010e-01; pointCoord(1117, 1) =5.995962798319616e-01; pointCoord(1117, 2) =2.288651636172858e-05;
pointCoord(1118, 0) = 3.800703679094010e-01; pointCoord(1118, 1) =5.995962798319616e-01; pointCoord(1118, 2) =2.288651636172858e-05;
pointCoord(1119, 0) = 1.800703679094010e-01; pointCoord(1119, 1) =7.995962798319617e-01; pointCoord(1119, 2) =2.288651636172858e-05;
pointCoord(1120, 0) = 1.990579411955255e-01; pointCoord(1120, 1) =1.534952997961074e-01; pointCoord(1120, 2) =7.533611717515951e-05;
pointCoord(1121, 0) = 3.990579411955255e-01; pointCoord(1121, 1) =1.534952997961074e-01; pointCoord(1121, 2) =7.533611717515951e-05;
pointCoord(1122, 0) = 5.990579411955255e-01; pointCoord(1122, 1) =1.534952997961074e-01; pointCoord(1122, 2) =7.533611717515951e-05;
pointCoord(1123, 0) = 7.990579411955255e-01; pointCoord(1123, 1) =1.534952997961074e-01; pointCoord(1123, 2) =7.533611717515951e-05;
pointCoord(1124, 0) = 1.990579411955255e-01; pointCoord(1124, 1) =3.534952997961074e-01; pointCoord(1124, 2) =7.533611717515951e-05;
pointCoord(1125, 0) = 3.990579411955255e-01; pointCoord(1125, 1) =3.534952997961074e-01; pointCoord(1125, 2) =7.533611717515951e-05;
pointCoord(1126, 0) = 5.990579411955255e-01; pointCoord(1126, 1) =3.534952997961074e-01; pointCoord(1126, 2) =7.533611717515951e-05;
pointCoord(1127, 0) = 1.990579411955255e-01; pointCoord(1127, 1) =5.534952997961075e-01; pointCoord(1127, 2) =7.533611717515951e-05;
pointCoord(1128, 0) = 3.990579411955255e-01; pointCoord(1128, 1) =5.534952997961075e-01; pointCoord(1128, 2) =7.533611717515951e-05;
pointCoord(1129, 0) = 1.990579411955255e-01; pointCoord(1129, 1) =7.534952997961075e-01; pointCoord(1129, 2) =7.533611717515951e-05;
pointCoord(1130, 0) = 1.951762416777610e-01; pointCoord(1130, 1) =1.573769993138719e-01; pointCoord(1130, 2) =1.655000090197417e-04;
pointCoord(1131, 0) = 3.951762416777610e-01; pointCoord(1131, 1) =1.573769993138719e-01; pointCoord(1131, 2) =1.655000090197417e-04;
pointCoord(1132, 0) = 5.951762416777610e-01; pointCoord(1132, 1) =1.573769993138719e-01; pointCoord(1132, 2) =1.655000090197417e-04;
pointCoord(1133, 0) = 7.951762416777610e-01; pointCoord(1133, 1) =1.573769993138719e-01; pointCoord(1133, 2) =1.655000090197417e-04;
pointCoord(1134, 0) = 1.951762416777610e-01; pointCoord(1134, 1) =3.573769993138720e-01; pointCoord(1134, 2) =1.655000090197417e-04;
pointCoord(1135, 0) = 3.951762416777610e-01; pointCoord(1135, 1) =3.573769993138720e-01; pointCoord(1135, 2) =1.655000090197417e-04;
pointCoord(1136, 0) = 5.951762416777610e-01; pointCoord(1136, 1) =3.573769993138720e-01; pointCoord(1136, 2) =1.655000090197417e-04;
pointCoord(1137, 0) = 1.951762416777610e-01; pointCoord(1137, 1) =5.573769993138719e-01; pointCoord(1137, 2) =1.655000090197417e-04;
pointCoord(1138, 0) = 3.951762416777610e-01; pointCoord(1138, 1) =5.573769993138719e-01; pointCoord(1138, 2) =1.655000090197417e-04;
pointCoord(1139, 0) = 1.951762416777610e-01; pointCoord(1139, 1) =7.573769993138719e-01; pointCoord(1139, 2) =1.655000090197417e-04;
pointCoord(1140, 0) = 1.887440252980097e-01; pointCoord(1140, 1) =1.638092156936232e-01; pointCoord(1140, 2) =2.334661894615331e-04;
pointCoord(1141, 0) = 3.887440252980097e-01; pointCoord(1141, 1) =1.638092156936232e-01; pointCoord(1141, 2) =2.334661894615331e-04;
pointCoord(1142, 0) = 5.887440252980096e-01; pointCoord(1142, 1) =1.638092156936232e-01; pointCoord(1142, 2) =2.334661894615331e-04;
pointCoord(1143, 0) = 7.887440252980097e-01; pointCoord(1143, 1) =1.638092156936232e-01; pointCoord(1143, 2) =2.334661894615331e-04;
pointCoord(1144, 0) = 1.887440252980097e-01; pointCoord(1144, 1) =3.638092156936232e-01; pointCoord(1144, 2) =2.334661894615331e-04;
pointCoord(1145, 0) = 3.887440252980097e-01; pointCoord(1145, 1) =3.638092156936232e-01; pointCoord(1145, 2) =2.334661894615331e-04;
pointCoord(1146, 0) = 5.887440252980096e-01; pointCoord(1146, 1) =3.638092156936232e-01; pointCoord(1146, 2) =2.334661894615331e-04;
pointCoord(1147, 0) = 1.887440252980097e-01; pointCoord(1147, 1) =5.638092156936232e-01; pointCoord(1147, 2) =2.334661894615331e-04;
pointCoord(1148, 0) = 3.887440252980097e-01; pointCoord(1148, 1) =5.638092156936232e-01; pointCoord(1148, 2) =2.334661894615331e-04;
pointCoord(1149, 0) = 1.887440252980097e-01; pointCoord(1149, 1) =7.638092156936233e-01; pointCoord(1149, 2) =2.334661894615331e-04;
pointCoord(1150, 0) = 1.806283101339550e-01; pointCoord(1150, 1) =1.719249308576779e-01; pointCoord(1150, 2) =2.699158656581297e-04;
pointCoord(1151, 0) = 3.806283101339550e-01; pointCoord(1151, 1) =1.719249308576779e-01; pointCoord(1151, 2) =2.699158656581297e-04;
pointCoord(1152, 0) = 5.806283101339550e-01; pointCoord(1152, 1) =1.719249308576779e-01; pointCoord(1152, 2) =2.699158656581297e-04;
pointCoord(1153, 0) = 7.806283101339551e-01; pointCoord(1153, 1) =1.719249308576779e-01; pointCoord(1153, 2) =2.699158656581297e-04;
pointCoord(1154, 0) = 1.806283101339550e-01; pointCoord(1154, 1) =3.719249308576779e-01; pointCoord(1154, 2) =2.699158656581297e-04;
pointCoord(1155, 0) = 3.806283101339550e-01; pointCoord(1155, 1) =3.719249308576779e-01; pointCoord(1155, 2) =2.699158656581297e-04;
pointCoord(1156, 0) = 5.806283101339550e-01; pointCoord(1156, 1) =3.719249308576779e-01; pointCoord(1156, 2) =2.699158656581297e-04;
pointCoord(1157, 0) = 1.806283101339550e-01; pointCoord(1157, 1) =5.719249308576779e-01; pointCoord(1157, 2) =2.699158656581297e-04;
pointCoord(1158, 0) = 3.806283101339550e-01; pointCoord(1158, 1) =5.719249308576779e-01; pointCoord(1158, 2) =2.699158656581297e-04;
pointCoord(1159, 0) = 1.806283101339550e-01; pointCoord(1159, 1) =7.719249308576780e-01; pointCoord(1159, 2) =2.699158656581297e-04;
pointCoord(1160, 0) = 1.719249308576779e-01; pointCoord(1160, 1) =1.806283101339550e-01; pointCoord(1160, 2) =2.699158656581297e-04;
pointCoord(1161, 0) = 3.719249308576779e-01; pointCoord(1161, 1) =1.806283101339550e-01; pointCoord(1161, 2) =2.699158656581297e-04;
pointCoord(1162, 0) = 5.719249308576779e-01; pointCoord(1162, 1) =1.806283101339550e-01; pointCoord(1162, 2) =2.699158656581297e-04;
pointCoord(1163, 0) = 7.719249308576780e-01; pointCoord(1163, 1) =1.806283101339550e-01; pointCoord(1163, 2) =2.699158656581297e-04;
pointCoord(1164, 0) = 1.719249308576779e-01; pointCoord(1164, 1) =3.806283101339550e-01; pointCoord(1164, 2) =2.699158656581297e-04;
pointCoord(1165, 0) = 3.719249308576779e-01; pointCoord(1165, 1) =3.806283101339550e-01; pointCoord(1165, 2) =2.699158656581297e-04;
pointCoord(1166, 0) = 5.719249308576779e-01; pointCoord(1166, 1) =3.806283101339550e-01; pointCoord(1166, 2) =2.699158656581297e-04;
pointCoord(1167, 0) = 1.719249308576779e-01; pointCoord(1167, 1) =5.806283101339550e-01; pointCoord(1167, 2) =2.699158656581297e-04;
pointCoord(1168, 0) = 3.719249308576779e-01; pointCoord(1168, 1) =5.806283101339550e-01; pointCoord(1168, 2) =2.699158656581297e-04;
pointCoord(1169, 0) = 1.719249308576779e-01; pointCoord(1169, 1) =7.806283101339551e-01; pointCoord(1169, 2) =2.699158656581297e-04;
pointCoord(1170, 0) = 1.638092156936232e-01; pointCoord(1170, 1) =1.887440252980097e-01; pointCoord(1170, 2) =2.334661894615331e-04;
pointCoord(1171, 0) = 3.638092156936232e-01; pointCoord(1171, 1) =1.887440252980097e-01; pointCoord(1171, 2) =2.334661894615331e-04;
pointCoord(1172, 0) = 5.638092156936232e-01; pointCoord(1172, 1) =1.887440252980097e-01; pointCoord(1172, 2) =2.334661894615331e-04;
pointCoord(1173, 0) = 7.638092156936233e-01; pointCoord(1173, 1) =1.887440252980097e-01; pointCoord(1173, 2) =2.334661894615331e-04;
pointCoord(1174, 0) = 1.638092156936232e-01; pointCoord(1174, 1) =3.887440252980097e-01; pointCoord(1174, 2) =2.334661894615331e-04;
pointCoord(1175, 0) = 3.638092156936232e-01; pointCoord(1175, 1) =3.887440252980097e-01; pointCoord(1175, 2) =2.334661894615331e-04;
pointCoord(1176, 0) = 5.638092156936232e-01; pointCoord(1176, 1) =3.887440252980097e-01; pointCoord(1176, 2) =2.334661894615331e-04;
pointCoord(1177, 0) = 1.638092156936232e-01; pointCoord(1177, 1) =5.887440252980096e-01; pointCoord(1177, 2) =2.334661894615331e-04;
pointCoord(1178, 0) = 3.638092156936232e-01; pointCoord(1178, 1) =5.887440252980096e-01; pointCoord(1178, 2) =2.334661894615331e-04;
pointCoord(1179, 0) = 1.638092156936232e-01; pointCoord(1179, 1) =7.887440252980097e-01; pointCoord(1179, 2) =2.334661894615331e-04;
pointCoord(1180, 0) = 1.573769993138719e-01; pointCoord(1180, 1) =1.951762416777610e-01; pointCoord(1180, 2) =1.655000090197417e-04;
pointCoord(1181, 0) = 3.573769993138720e-01; pointCoord(1181, 1) =1.951762416777610e-01; pointCoord(1181, 2) =1.655000090197417e-04;
pointCoord(1182, 0) = 5.573769993138719e-01; pointCoord(1182, 1) =1.951762416777610e-01; pointCoord(1182, 2) =1.655000090197417e-04;
pointCoord(1183, 0) = 7.573769993138719e-01; pointCoord(1183, 1) =1.951762416777610e-01; pointCoord(1183, 2) =1.655000090197417e-04;
pointCoord(1184, 0) = 1.573769993138719e-01; pointCoord(1184, 1) =3.951762416777610e-01; pointCoord(1184, 2) =1.655000090197417e-04;
pointCoord(1185, 0) = 3.573769993138720e-01; pointCoord(1185, 1) =3.951762416777610e-01; pointCoord(1185, 2) =1.655000090197417e-04;
pointCoord(1186, 0) = 5.573769993138719e-01; pointCoord(1186, 1) =3.951762416777610e-01; pointCoord(1186, 2) =1.655000090197417e-04;
pointCoord(1187, 0) = 1.573769993138719e-01; pointCoord(1187, 1) =5.951762416777610e-01; pointCoord(1187, 2) =1.655000090197417e-04;
pointCoord(1188, 0) = 3.573769993138720e-01; pointCoord(1188, 1) =5.951762416777610e-01; pointCoord(1188, 2) =1.655000090197417e-04;
pointCoord(1189, 0) = 1.573769993138719e-01; pointCoord(1189, 1) =7.951762416777610e-01; pointCoord(1189, 2) =1.655000090197417e-04;
pointCoord(1190, 0) = 1.534952997961074e-01; pointCoord(1190, 1) =1.990579411955255e-01; pointCoord(1190, 2) =7.533611717515951e-05;
pointCoord(1191, 0) = 3.534952997961074e-01; pointCoord(1191, 1) =1.990579411955255e-01; pointCoord(1191, 2) =7.533611717515951e-05;
pointCoord(1192, 0) = 5.534952997961075e-01; pointCoord(1192, 1) =1.990579411955255e-01; pointCoord(1192, 2) =7.533611717515951e-05;
pointCoord(1193, 0) = 7.534952997961075e-01; pointCoord(1193, 1) =1.990579411955255e-01; pointCoord(1193, 2) =7.533611717515951e-05;
pointCoord(1194, 0) = 1.534952997961074e-01; pointCoord(1194, 1) =3.990579411955255e-01; pointCoord(1194, 2) =7.533611717515951e-05;
pointCoord(1195, 0) = 3.534952997961074e-01; pointCoord(1195, 1) =3.990579411955255e-01; pointCoord(1195, 2) =7.533611717515951e-05;
pointCoord(1196, 0) = 5.534952997961075e-01; pointCoord(1196, 1) =3.990579411955255e-01; pointCoord(1196, 2) =7.533611717515951e-05;
pointCoord(1197, 0) = 1.534952997961074e-01; pointCoord(1197, 1) =5.990579411955255e-01; pointCoord(1197, 2) =7.533611717515951e-05;
pointCoord(1198, 0) = 3.534952997961074e-01; pointCoord(1198, 1) =5.990579411955255e-01; pointCoord(1198, 2) =7.533611717515951e-05;
pointCoord(1199, 0) = 1.534952997961074e-01; pointCoord(1199, 1) =7.990579411955255e-01; pointCoord(1199, 2) =7.533611717515951e-05;
pointCoord(1200, 0) = 1.983787036237181e-01; pointCoord(1200, 1) =1.199647606258469e-01; pointCoord(1200, 2) =1.498966925243746e-04;
pointCoord(1201, 0) = 3.983787036237181e-01; pointCoord(1201, 1) =1.199647606258469e-01; pointCoord(1201, 2) =1.498966925243746e-04;
pointCoord(1202, 0) = 5.983787036237180e-01; pointCoord(1202, 1) =1.199647606258469e-01; pointCoord(1202, 2) =1.498966925243746e-04;
pointCoord(1203, 0) = 7.983787036237181e-01; pointCoord(1203, 1) =1.199647606258469e-01; pointCoord(1203, 2) =1.498966925243746e-04;
pointCoord(1204, 0) = 1.983787036237181e-01; pointCoord(1204, 1) =3.199647606258469e-01; pointCoord(1204, 2) =1.498966925243746e-04;
pointCoord(1205, 0) = 3.983787036237181e-01; pointCoord(1205, 1) =3.199647606258469e-01; pointCoord(1205, 2) =1.498966925243746e-04;
pointCoord(1206, 0) = 5.983787036237180e-01; pointCoord(1206, 1) =3.199647606258469e-01; pointCoord(1206, 2) =1.498966925243746e-04;
pointCoord(1207, 0) = 1.983787036237181e-01; pointCoord(1207, 1) =5.199647606258468e-01; pointCoord(1207, 2) =1.498966925243746e-04;
pointCoord(1208, 0) = 3.983787036237181e-01; pointCoord(1208, 1) =5.199647606258468e-01; pointCoord(1208, 2) =1.498966925243746e-04;
pointCoord(1209, 0) = 1.983787036237181e-01; pointCoord(1209, 1) =7.199647606258469e-01; pointCoord(1209, 2) =1.498966925243746e-04;
pointCoord(1210, 0) = 1.916982444718320e-01; pointCoord(1210, 1) =1.266452197777330e-01; pointCoord(1210, 2) =3.292962910091858e-04;
pointCoord(1211, 0) = 3.916982444718320e-01; pointCoord(1211, 1) =1.266452197777330e-01; pointCoord(1211, 2) =3.292962910091858e-04;
pointCoord(1212, 0) = 5.916982444718319e-01; pointCoord(1212, 1) =1.266452197777330e-01; pointCoord(1212, 2) =3.292962910091858e-04;
pointCoord(1213, 0) = 7.916982444718320e-01; pointCoord(1213, 1) =1.266452197777330e-01; pointCoord(1213, 2) =3.292962910091858e-04;
pointCoord(1214, 0) = 1.916982444718320e-01; pointCoord(1214, 1) =3.266452197777330e-01; pointCoord(1214, 2) =3.292962910091858e-04;
pointCoord(1215, 0) = 3.916982444718320e-01; pointCoord(1215, 1) =3.266452197777330e-01; pointCoord(1215, 2) =3.292962910091858e-04;
pointCoord(1216, 0) = 5.916982444718319e-01; pointCoord(1216, 1) =3.266452197777330e-01; pointCoord(1216, 2) =3.292962910091858e-04;
pointCoord(1217, 0) = 1.916982444718320e-01; pointCoord(1217, 1) =5.266452197777330e-01; pointCoord(1217, 2) =3.292962910091858e-04;
pointCoord(1218, 0) = 3.916982444718320e-01; pointCoord(1218, 1) =5.266452197777330e-01; pointCoord(1218, 2) =3.292962910091858e-04;
pointCoord(1219, 0) = 1.916982444718320e-01; pointCoord(1219, 1) =7.266452197777331e-01; pointCoord(1219, 2) =3.292962910091858e-04;
pointCoord(1220, 0) = 1.806283101339550e-01; pointCoord(1220, 1) =1.377151541156100e-01; pointCoord(1220, 2) =4.645289793099657e-04;
pointCoord(1221, 0) = 3.806283101339550e-01; pointCoord(1221, 1) =1.377151541156100e-01; pointCoord(1221, 2) =4.645289793099657e-04;
pointCoord(1222, 0) = 5.806283101339550e-01; pointCoord(1222, 1) =1.377151541156100e-01; pointCoord(1222, 2) =4.645289793099657e-04;
pointCoord(1223, 0) = 7.806283101339551e-01; pointCoord(1223, 1) =1.377151541156100e-01; pointCoord(1223, 2) =4.645289793099657e-04;
pointCoord(1224, 0) = 1.806283101339550e-01; pointCoord(1224, 1) =3.377151541156100e-01; pointCoord(1224, 2) =4.645289793099657e-04;
pointCoord(1225, 0) = 3.806283101339550e-01; pointCoord(1225, 1) =3.377151541156100e-01; pointCoord(1225, 2) =4.645289793099657e-04;
pointCoord(1226, 0) = 5.806283101339550e-01; pointCoord(1226, 1) =3.377151541156100e-01; pointCoord(1226, 2) =4.645289793099657e-04;
pointCoord(1227, 0) = 1.806283101339550e-01; pointCoord(1227, 1) =5.377151541156100e-01; pointCoord(1227, 2) =4.645289793099657e-04;
pointCoord(1228, 0) = 3.806283101339550e-01; pointCoord(1228, 1) =5.377151541156100e-01; pointCoord(1228, 2) =4.645289793099657e-04;
pointCoord(1229, 0) = 1.806283101339550e-01; pointCoord(1229, 1) =7.377151541156101e-01; pointCoord(1229, 2) =4.645289793099657e-04;
pointCoord(1230, 0) = 1.666610508461897e-01; pointCoord(1230, 1) =1.516824134033753e-01; pointCoord(1230, 2) =5.370531033333869e-04;
pointCoord(1231, 0) = 3.666610508461897e-01; pointCoord(1231, 1) =1.516824134033753e-01; pointCoord(1231, 2) =5.370531033333869e-04;
pointCoord(1232, 0) = 5.666610508461896e-01; pointCoord(1232, 1) =1.516824134033753e-01; pointCoord(1232, 2) =5.370531033333869e-04;
pointCoord(1233, 0) = 7.666610508461896e-01; pointCoord(1233, 1) =1.516824134033753e-01; pointCoord(1233, 2) =5.370531033333869e-04;
pointCoord(1234, 0) = 1.666610508461897e-01; pointCoord(1234, 1) =3.516824134033754e-01; pointCoord(1234, 2) =5.370531033333869e-04;
pointCoord(1235, 0) = 3.666610508461897e-01; pointCoord(1235, 1) =3.516824134033754e-01; pointCoord(1235, 2) =5.370531033333869e-04;
pointCoord(1236, 0) = 5.666610508461896e-01; pointCoord(1236, 1) =3.516824134033754e-01; pointCoord(1236, 2) =5.370531033333869e-04;
pointCoord(1237, 0) = 1.666610508461897e-01; pointCoord(1237, 1) =5.516824134033753e-01; pointCoord(1237, 2) =5.370531033333869e-04;
pointCoord(1238, 0) = 3.666610508461897e-01; pointCoord(1238, 1) =5.516824134033753e-01; pointCoord(1238, 2) =5.370531033333869e-04;
pointCoord(1239, 0) = 1.666610508461897e-01; pointCoord(1239, 1) =7.516824134033754e-01; pointCoord(1239, 2) =5.370531033333869e-04;
pointCoord(1240, 0) = 1.516824134033753e-01; pointCoord(1240, 1) =1.666610508461897e-01; pointCoord(1240, 2) =5.370531033333869e-04;
pointCoord(1241, 0) = 3.516824134033754e-01; pointCoord(1241, 1) =1.666610508461897e-01; pointCoord(1241, 2) =5.370531033333869e-04;
pointCoord(1242, 0) = 5.516824134033753e-01; pointCoord(1242, 1) =1.666610508461897e-01; pointCoord(1242, 2) =5.370531033333869e-04;
pointCoord(1243, 0) = 7.516824134033754e-01; pointCoord(1243, 1) =1.666610508461897e-01; pointCoord(1243, 2) =5.370531033333869e-04;
pointCoord(1244, 0) = 1.516824134033753e-01; pointCoord(1244, 1) =3.666610508461897e-01; pointCoord(1244, 2) =5.370531033333869e-04;
pointCoord(1245, 0) = 3.516824134033754e-01; pointCoord(1245, 1) =3.666610508461897e-01; pointCoord(1245, 2) =5.370531033333869e-04;
pointCoord(1246, 0) = 5.516824134033753e-01; pointCoord(1246, 1) =3.666610508461897e-01; pointCoord(1246, 2) =5.370531033333869e-04;
pointCoord(1247, 0) = 1.516824134033753e-01; pointCoord(1247, 1) =5.666610508461896e-01; pointCoord(1247, 2) =5.370531033333869e-04;
pointCoord(1248, 0) = 3.516824134033754e-01; pointCoord(1248, 1) =5.666610508461896e-01; pointCoord(1248, 2) =5.370531033333869e-04;
pointCoord(1249, 0) = 1.516824134033753e-01; pointCoord(1249, 1) =7.666610508461896e-01; pointCoord(1249, 2) =5.370531033333869e-04;
pointCoord(1250, 0) = 1.377151541156100e-01; pointCoord(1250, 1) =1.806283101339550e-01; pointCoord(1250, 2) =4.645289793099657e-04;
pointCoord(1251, 0) = 3.377151541156100e-01; pointCoord(1251, 1) =1.806283101339550e-01; pointCoord(1251, 2) =4.645289793099657e-04;
pointCoord(1252, 0) = 5.377151541156100e-01; pointCoord(1252, 1) =1.806283101339550e-01; pointCoord(1252, 2) =4.645289793099657e-04;
pointCoord(1253, 0) = 7.377151541156101e-01; pointCoord(1253, 1) =1.806283101339550e-01; pointCoord(1253, 2) =4.645289793099657e-04;
pointCoord(1254, 0) = 1.377151541156100e-01; pointCoord(1254, 1) =3.806283101339550e-01; pointCoord(1254, 2) =4.645289793099657e-04;
pointCoord(1255, 0) = 3.377151541156100e-01; pointCoord(1255, 1) =3.806283101339550e-01; pointCoord(1255, 2) =4.645289793099657e-04;
pointCoord(1256, 0) = 5.377151541156100e-01; pointCoord(1256, 1) =3.806283101339550e-01; pointCoord(1256, 2) =4.645289793099657e-04;
pointCoord(1257, 0) = 1.377151541156100e-01; pointCoord(1257, 1) =5.806283101339550e-01; pointCoord(1257, 2) =4.645289793099657e-04;
pointCoord(1258, 0) = 3.377151541156100e-01; pointCoord(1258, 1) =5.806283101339550e-01; pointCoord(1258, 2) =4.645289793099657e-04;
pointCoord(1259, 0) = 1.377151541156100e-01; pointCoord(1259, 1) =7.806283101339551e-01; pointCoord(1259, 2) =4.645289793099657e-04;
pointCoord(1260, 0) = 1.266452197777330e-01; pointCoord(1260, 1) =1.916982444718320e-01; pointCoord(1260, 2) =3.292962910091858e-04;
pointCoord(1261, 0) = 3.266452197777330e-01; pointCoord(1261, 1) =1.916982444718320e-01; pointCoord(1261, 2) =3.292962910091858e-04;
pointCoord(1262, 0) = 5.266452197777330e-01; pointCoord(1262, 1) =1.916982444718320e-01; pointCoord(1262, 2) =3.292962910091858e-04;
pointCoord(1263, 0) = 7.266452197777331e-01; pointCoord(1263, 1) =1.916982444718320e-01; pointCoord(1263, 2) =3.292962910091858e-04;
pointCoord(1264, 0) = 1.266452197777330e-01; pointCoord(1264, 1) =3.916982444718320e-01; pointCoord(1264, 2) =3.292962910091858e-04;
pointCoord(1265, 0) = 3.266452197777330e-01; pointCoord(1265, 1) =3.916982444718320e-01; pointCoord(1265, 2) =3.292962910091858e-04;
pointCoord(1266, 0) = 5.266452197777330e-01; pointCoord(1266, 1) =3.916982444718320e-01; pointCoord(1266, 2) =3.292962910091858e-04;
pointCoord(1267, 0) = 1.266452197777330e-01; pointCoord(1267, 1) =5.916982444718319e-01; pointCoord(1267, 2) =3.292962910091858e-04;
pointCoord(1268, 0) = 3.266452197777330e-01; pointCoord(1268, 1) =5.916982444718319e-01; pointCoord(1268, 2) =3.292962910091858e-04;
pointCoord(1269, 0) = 1.266452197777330e-01; pointCoord(1269, 1) =7.916982444718320e-01; pointCoord(1269, 2) =3.292962910091858e-04;
pointCoord(1270, 0) = 1.199647606258469e-01; pointCoord(1270, 1) =1.983787036237181e-01; pointCoord(1270, 2) =1.498966925243746e-04;
pointCoord(1271, 0) = 3.199647606258469e-01; pointCoord(1271, 1) =1.983787036237181e-01; pointCoord(1271, 2) =1.498966925243746e-04;
pointCoord(1272, 0) = 5.199647606258468e-01; pointCoord(1272, 1) =1.983787036237181e-01; pointCoord(1272, 2) =1.498966925243746e-04;
pointCoord(1273, 0) = 7.199647606258469e-01; pointCoord(1273, 1) =1.983787036237181e-01; pointCoord(1273, 2) =1.498966925243746e-04;
pointCoord(1274, 0) = 1.199647606258469e-01; pointCoord(1274, 1) =3.983787036237181e-01; pointCoord(1274, 2) =1.498966925243746e-04;
pointCoord(1275, 0) = 3.199647606258469e-01; pointCoord(1275, 1) =3.983787036237181e-01; pointCoord(1275, 2) =1.498966925243746e-04;
pointCoord(1276, 0) = 5.199647606258468e-01; pointCoord(1276, 1) =3.983787036237181e-01; pointCoord(1276, 2) =1.498966925243746e-04;
pointCoord(1277, 0) = 1.199647606258469e-01; pointCoord(1277, 1) =5.983787036237180e-01; pointCoord(1277, 2) =1.498966925243746e-04;
pointCoord(1278, 0) = 3.199647606258469e-01; pointCoord(1278, 1) =5.983787036237180e-01; pointCoord(1278, 2) =1.498966925243746e-04;
pointCoord(1279, 0) = 1.199647606258469e-01; pointCoord(1279, 1) =7.983787036237181e-01; pointCoord(1279, 2) =1.498966925243746e-04;
pointCoord(1280, 0) = 1.976502820260355e-01; pointCoord(1280, 1) =8.400625372439949e-02; pointCoord(1280, 2) =2.172427927521020e-04;
pointCoord(1281, 0) = 3.976502820260356e-01; pointCoord(1281, 1) =8.400625372439949e-02; pointCoord(1281, 2) =2.172427927521020e-04;
pointCoord(1282, 0) = 5.976502820260355e-01; pointCoord(1282, 1) =8.400625372439949e-02; pointCoord(1282, 2) =2.172427927521020e-04;
pointCoord(1283, 0) = 7.976502820260356e-01; pointCoord(1283, 1) =8.400625372439949e-02; pointCoord(1283, 2) =2.172427927521020e-04;
pointCoord(1284, 0) = 1.976502820260355e-01; pointCoord(1284, 1) =2.840062537243995e-01; pointCoord(1284, 2) =2.172427927521020e-04;
pointCoord(1285, 0) = 3.976502820260356e-01; pointCoord(1285, 1) =2.840062537243995e-01; pointCoord(1285, 2) =2.172427927521020e-04;
pointCoord(1286, 0) = 5.976502820260355e-01; pointCoord(1286, 1) =2.840062537243995e-01; pointCoord(1286, 2) =2.172427927521020e-04;
pointCoord(1287, 0) = 1.976502820260355e-01; pointCoord(1287, 1) =4.840062537243994e-01; pointCoord(1287, 2) =2.172427927521020e-04;
pointCoord(1288, 0) = 3.976502820260356e-01; pointCoord(1288, 1) =4.840062537243994e-01; pointCoord(1288, 2) =2.172427927521020e-04;
pointCoord(1289, 0) = 1.976502820260355e-01; pointCoord(1289, 1) =6.840062537243995e-01; pointCoord(1289, 2) =2.172427927521020e-04;
pointCoord(1290, 0) = 1.879684032695307e-01; pointCoord(1290, 1) =9.368813248090434e-02; pointCoord(1290, 2) =4.772436582622514e-04;
pointCoord(1291, 0) = 3.879684032695307e-01; pointCoord(1291, 1) =9.368813248090434e-02; pointCoord(1291, 2) =4.772436582622514e-04;
pointCoord(1292, 0) = 5.879684032695307e-01; pointCoord(1292, 1) =9.368813248090434e-02; pointCoord(1292, 2) =4.772436582622514e-04;
pointCoord(1293, 0) = 7.879684032695308e-01; pointCoord(1293, 1) =9.368813248090434e-02; pointCoord(1293, 2) =4.772436582622514e-04;
pointCoord(1294, 0) = 1.879684032695307e-01; pointCoord(1294, 1) =2.936881324809044e-01; pointCoord(1294, 2) =4.772436582622514e-04;
pointCoord(1295, 0) = 3.879684032695307e-01; pointCoord(1295, 1) =2.936881324809044e-01; pointCoord(1295, 2) =4.772436582622514e-04;
pointCoord(1296, 0) = 5.879684032695307e-01; pointCoord(1296, 1) =2.936881324809044e-01; pointCoord(1296, 2) =4.772436582622514e-04;
pointCoord(1297, 0) = 1.879684032695307e-01; pointCoord(1297, 1) =4.936881324809043e-01; pointCoord(1297, 2) =4.772436582622514e-04;
pointCoord(1298, 0) = 3.879684032695307e-01; pointCoord(1298, 1) =4.936881324809043e-01; pointCoord(1298, 2) =4.772436582622514e-04;
pointCoord(1299, 0) = 1.879684032695307e-01; pointCoord(1299, 1) =6.936881324809043e-01; pointCoord(1299, 2) =4.772436582622514e-04;
pointCoord(1300, 0) = 1.719249308576779e-01; pointCoord(1300, 1) =1.097316048927571e-01; pointCoord(1300, 2) =6.732341526693159e-04;
pointCoord(1301, 0) = 3.719249308576779e-01; pointCoord(1301, 1) =1.097316048927571e-01; pointCoord(1301, 2) =6.732341526693159e-04;
pointCoord(1302, 0) = 5.719249308576779e-01; pointCoord(1302, 1) =1.097316048927571e-01; pointCoord(1302, 2) =6.732341526693159e-04;
pointCoord(1303, 0) = 7.719249308576780e-01; pointCoord(1303, 1) =1.097316048927571e-01; pointCoord(1303, 2) =6.732341526693159e-04;
pointCoord(1304, 0) = 1.719249308576779e-01; pointCoord(1304, 1) =3.097316048927571e-01; pointCoord(1304, 2) =6.732341526693159e-04;
pointCoord(1305, 0) = 3.719249308576779e-01; pointCoord(1305, 1) =3.097316048927571e-01; pointCoord(1305, 2) =6.732341526693159e-04;
pointCoord(1306, 0) = 5.719249308576779e-01; pointCoord(1306, 1) =3.097316048927571e-01; pointCoord(1306, 2) =6.732341526693159e-04;
pointCoord(1307, 0) = 1.719249308576779e-01; pointCoord(1307, 1) =5.097316048927572e-01; pointCoord(1307, 2) =6.732341526693159e-04;
pointCoord(1308, 0) = 3.719249308576779e-01; pointCoord(1308, 1) =5.097316048927572e-01; pointCoord(1308, 2) =6.732341526693159e-04;
pointCoord(1309, 0) = 1.719249308576779e-01; pointCoord(1309, 1) =7.097316048927571e-01; pointCoord(1309, 2) =6.732341526693159e-04;
pointCoord(1310, 0) = 1.516824134033753e-01; pointCoord(1310, 1) =1.299741223470597e-01; pointCoord(1310, 2) =7.783421639230390e-04;
pointCoord(1311, 0) = 3.516824134033754e-01; pointCoord(1311, 1) =1.299741223470597e-01; pointCoord(1311, 2) =7.783421639230390e-04;
pointCoord(1312, 0) = 5.516824134033753e-01; pointCoord(1312, 1) =1.299741223470597e-01; pointCoord(1312, 2) =7.783421639230390e-04;
pointCoord(1313, 0) = 7.516824134033754e-01; pointCoord(1313, 1) =1.299741223470597e-01; pointCoord(1313, 2) =7.783421639230390e-04;
pointCoord(1314, 0) = 1.516824134033753e-01; pointCoord(1314, 1) =3.299741223470597e-01; pointCoord(1314, 2) =7.783421639230390e-04;
pointCoord(1315, 0) = 3.516824134033754e-01; pointCoord(1315, 1) =3.299741223470597e-01; pointCoord(1315, 2) =7.783421639230390e-04;
pointCoord(1316, 0) = 5.516824134033753e-01; pointCoord(1316, 1) =3.299741223470597e-01; pointCoord(1316, 2) =7.783421639230390e-04;
pointCoord(1317, 0) = 1.516824134033753e-01; pointCoord(1317, 1) =5.299741223470597e-01; pointCoord(1317, 2) =7.783421639230390e-04;
pointCoord(1318, 0) = 3.516824134033754e-01; pointCoord(1318, 1) =5.299741223470597e-01; pointCoord(1318, 2) =7.783421639230390e-04;
pointCoord(1319, 0) = 1.516824134033753e-01; pointCoord(1319, 1) =7.299741223470597e-01; pointCoord(1319, 2) =7.783421639230390e-04;
pointCoord(1320, 0) = 1.299741223470597e-01; pointCoord(1320, 1) =1.516824134033753e-01; pointCoord(1320, 2) =7.783421639230390e-04;
pointCoord(1321, 0) = 3.299741223470597e-01; pointCoord(1321, 1) =1.516824134033753e-01; pointCoord(1321, 2) =7.783421639230390e-04;
pointCoord(1322, 0) = 5.299741223470597e-01; pointCoord(1322, 1) =1.516824134033753e-01; pointCoord(1322, 2) =7.783421639230390e-04;
pointCoord(1323, 0) = 7.299741223470597e-01; pointCoord(1323, 1) =1.516824134033753e-01; pointCoord(1323, 2) =7.783421639230390e-04;
pointCoord(1324, 0) = 1.299741223470597e-01; pointCoord(1324, 1) =3.516824134033754e-01; pointCoord(1324, 2) =7.783421639230390e-04;
pointCoord(1325, 0) = 3.299741223470597e-01; pointCoord(1325, 1) =3.516824134033754e-01; pointCoord(1325, 2) =7.783421639230390e-04;
pointCoord(1326, 0) = 5.299741223470597e-01; pointCoord(1326, 1) =3.516824134033754e-01; pointCoord(1326, 2) =7.783421639230390e-04;
pointCoord(1327, 0) = 1.299741223470597e-01; pointCoord(1327, 1) =5.516824134033753e-01; pointCoord(1327, 2) =7.783421639230390e-04;
pointCoord(1328, 0) = 3.299741223470597e-01; pointCoord(1328, 1) =5.516824134033753e-01; pointCoord(1328, 2) =7.783421639230390e-04;
pointCoord(1329, 0) = 1.299741223470597e-01; pointCoord(1329, 1) =7.516824134033754e-01; pointCoord(1329, 2) =7.783421639230390e-04;
pointCoord(1330, 0) = 1.097316048927571e-01; pointCoord(1330, 1) =1.719249308576779e-01; pointCoord(1330, 2) =6.732341526693159e-04;
pointCoord(1331, 0) = 3.097316048927571e-01; pointCoord(1331, 1) =1.719249308576779e-01; pointCoord(1331, 2) =6.732341526693159e-04;
pointCoord(1332, 0) = 5.097316048927572e-01; pointCoord(1332, 1) =1.719249308576779e-01; pointCoord(1332, 2) =6.732341526693159e-04;
pointCoord(1333, 0) = 7.097316048927571e-01; pointCoord(1333, 1) =1.719249308576779e-01; pointCoord(1333, 2) =6.732341526693159e-04;
pointCoord(1334, 0) = 1.097316048927571e-01; pointCoord(1334, 1) =3.719249308576779e-01; pointCoord(1334, 2) =6.732341526693159e-04;
pointCoord(1335, 0) = 3.097316048927571e-01; pointCoord(1335, 1) =3.719249308576779e-01; pointCoord(1335, 2) =6.732341526693159e-04;
pointCoord(1336, 0) = 5.097316048927572e-01; pointCoord(1336, 1) =3.719249308576779e-01; pointCoord(1336, 2) =6.732341526693159e-04;
pointCoord(1337, 0) = 1.097316048927571e-01; pointCoord(1337, 1) =5.719249308576779e-01; pointCoord(1337, 2) =6.732341526693159e-04;
pointCoord(1338, 0) = 3.097316048927571e-01; pointCoord(1338, 1) =5.719249308576779e-01; pointCoord(1338, 2) =6.732341526693159e-04;
pointCoord(1339, 0) = 1.097316048927571e-01; pointCoord(1339, 1) =7.719249308576780e-01; pointCoord(1339, 2) =6.732341526693159e-04;
pointCoord(1340, 0) = 9.368813248090434e-02; pointCoord(1340, 1) =1.879684032695307e-01; pointCoord(1340, 2) =4.772436582622514e-04;
pointCoord(1341, 0) = 2.936881324809044e-01; pointCoord(1341, 1) =1.879684032695307e-01; pointCoord(1341, 2) =4.772436582622514e-04;
pointCoord(1342, 0) = 4.936881324809043e-01; pointCoord(1342, 1) =1.879684032695307e-01; pointCoord(1342, 2) =4.772436582622514e-04;
pointCoord(1343, 0) = 6.936881324809043e-01; pointCoord(1343, 1) =1.879684032695307e-01; pointCoord(1343, 2) =4.772436582622514e-04;
pointCoord(1344, 0) = 9.368813248090434e-02; pointCoord(1344, 1) =3.879684032695307e-01; pointCoord(1344, 2) =4.772436582622514e-04;
pointCoord(1345, 0) = 2.936881324809044e-01; pointCoord(1345, 1) =3.879684032695307e-01; pointCoord(1345, 2) =4.772436582622514e-04;
pointCoord(1346, 0) = 4.936881324809043e-01; pointCoord(1346, 1) =3.879684032695307e-01; pointCoord(1346, 2) =4.772436582622514e-04;
pointCoord(1347, 0) = 9.368813248090434e-02; pointCoord(1347, 1) =5.879684032695307e-01; pointCoord(1347, 2) =4.772436582622514e-04;
pointCoord(1348, 0) = 2.936881324809044e-01; pointCoord(1348, 1) =5.879684032695307e-01; pointCoord(1348, 2) =4.772436582622514e-04;
pointCoord(1349, 0) = 9.368813248090434e-02; pointCoord(1349, 1) =7.879684032695308e-01; pointCoord(1349, 2) =4.772436582622514e-04;
pointCoord(1350, 0) = 8.400625372439949e-02; pointCoord(1350, 1) =1.976502820260355e-01; pointCoord(1350, 2) =2.172427927521020e-04;
pointCoord(1351, 0) = 2.840062537243995e-01; pointCoord(1351, 1) =1.976502820260355e-01; pointCoord(1351, 2) =2.172427927521020e-04;
pointCoord(1352, 0) = 4.840062537243994e-01; pointCoord(1352, 1) =1.976502820260355e-01; pointCoord(1352, 2) =2.172427927521020e-04;
pointCoord(1353, 0) = 6.840062537243995e-01; pointCoord(1353, 1) =1.976502820260355e-01; pointCoord(1353, 2) =2.172427927521020e-04;
pointCoord(1354, 0) = 8.400625372439949e-02; pointCoord(1354, 1) =3.976502820260356e-01; pointCoord(1354, 2) =2.172427927521020e-04;
pointCoord(1355, 0) = 2.840062537243995e-01; pointCoord(1355, 1) =3.976502820260356e-01; pointCoord(1355, 2) =2.172427927521020e-04;
pointCoord(1356, 0) = 4.840062537243994e-01; pointCoord(1356, 1) =3.976502820260356e-01; pointCoord(1356, 2) =2.172427927521020e-04;
pointCoord(1357, 0) = 8.400625372439949e-02; pointCoord(1357, 1) =5.976502820260355e-01; pointCoord(1357, 2) =2.172427927521020e-04;
pointCoord(1358, 0) = 2.840062537243995e-01; pointCoord(1358, 1) =5.976502820260355e-01; pointCoord(1358, 2) =2.172427927521020e-04;
pointCoord(1359, 0) = 8.400625372439949e-02; pointCoord(1359, 1) =7.976502820260356e-01; pointCoord(1359, 2) =2.172427927521020e-04;
pointCoord(1360, 0) = 1.969710444542282e-01; pointCoord(1360, 1) =5.047571455413893e-02; pointCoord(1360, 2) =2.422245286926613e-04;
pointCoord(1361, 0) = 3.969710444542282e-01; pointCoord(1361, 1) =5.047571455413893e-02; pointCoord(1361, 2) =2.422245286926613e-04;
pointCoord(1362, 0) = 5.969710444542281e-01; pointCoord(1362, 1) =5.047571455413893e-02; pointCoord(1362, 2) =2.422245286926613e-04;
pointCoord(1363, 0) = 7.969710444542282e-01; pointCoord(1363, 1) =5.047571455413893e-02; pointCoord(1363, 2) =2.422245286926613e-04;
pointCoord(1364, 0) = 1.969710444542282e-01; pointCoord(1364, 1) =2.504757145541390e-01; pointCoord(1364, 2) =2.422245286926613e-04;
pointCoord(1365, 0) = 3.969710444542282e-01; pointCoord(1365, 1) =2.504757145541390e-01; pointCoord(1365, 2) =2.422245286926613e-04;
pointCoord(1366, 0) = 5.969710444542281e-01; pointCoord(1366, 1) =2.504757145541390e-01; pointCoord(1366, 2) =2.422245286926613e-04;
pointCoord(1367, 0) = 1.969710444542282e-01; pointCoord(1367, 1) =4.504757145541389e-01; pointCoord(1367, 2) =2.422245286926613e-04;
pointCoord(1368, 0) = 3.969710444542282e-01; pointCoord(1368, 1) =4.504757145541389e-01; pointCoord(1368, 2) =2.422245286926613e-04;
pointCoord(1369, 0) = 1.969710444542282e-01; pointCoord(1369, 1) =6.504757145541390e-01; pointCoord(1369, 2) =2.422245286926613e-04;
pointCoord(1370, 0) = 1.844904060636017e-01; pointCoord(1370, 1) =6.295635294476543e-02; pointCoord(1370, 2) =5.321240752324881e-04;
pointCoord(1371, 0) = 3.844904060636017e-01; pointCoord(1371, 1) =6.295635294476543e-02; pointCoord(1371, 2) =5.321240752324881e-04;
pointCoord(1372, 0) = 5.844904060636017e-01; pointCoord(1372, 1) =6.295635294476543e-02; pointCoord(1372, 2) =5.321240752324881e-04;
pointCoord(1373, 0) = 7.844904060636018e-01; pointCoord(1373, 1) =6.295635294476543e-02; pointCoord(1373, 2) =5.321240752324881e-04;
pointCoord(1374, 0) = 1.844904060636017e-01; pointCoord(1374, 1) =2.629563529447654e-01; pointCoord(1374, 2) =5.321240752324881e-04;
pointCoord(1375, 0) = 3.844904060636017e-01; pointCoord(1375, 1) =2.629563529447654e-01; pointCoord(1375, 2) =5.321240752324881e-04;
pointCoord(1376, 0) = 5.844904060636017e-01; pointCoord(1376, 1) =2.629563529447654e-01; pointCoord(1376, 2) =5.321240752324881e-04;
pointCoord(1377, 0) = 1.844904060636017e-01; pointCoord(1377, 1) =4.629563529447654e-01; pointCoord(1377, 2) =5.321240752324881e-04;
pointCoord(1378, 0) = 3.844904060636017e-01; pointCoord(1378, 1) =4.629563529447654e-01; pointCoord(1378, 2) =5.321240752324881e-04;
pointCoord(1379, 0) = 1.844904060636017e-01; pointCoord(1379, 1) =6.629563529447655e-01; pointCoord(1379, 2) =5.321240752324881e-04;
pointCoord(1380, 0) = 1.638092156936232e-01; pointCoord(1380, 1) =8.363754331474391e-02; pointCoord(1380, 2) =7.506524072180082e-04;
pointCoord(1381, 0) = 3.638092156936232e-01; pointCoord(1381, 1) =8.363754331474391e-02; pointCoord(1381, 2) =7.506524072180082e-04;
pointCoord(1382, 0) = 5.638092156936232e-01; pointCoord(1382, 1) =8.363754331474391e-02; pointCoord(1382, 2) =7.506524072180082e-04;
pointCoord(1383, 0) = 7.638092156936233e-01; pointCoord(1383, 1) =8.363754331474391e-02; pointCoord(1383, 2) =7.506524072180082e-04;
pointCoord(1384, 0) = 1.638092156936232e-01; pointCoord(1384, 1) =2.836375433147439e-01; pointCoord(1384, 2) =7.506524072180082e-04;
pointCoord(1385, 0) = 3.638092156936232e-01; pointCoord(1385, 1) =2.836375433147439e-01; pointCoord(1385, 2) =7.506524072180082e-04;
pointCoord(1386, 0) = 5.638092156936232e-01; pointCoord(1386, 1) =2.836375433147439e-01; pointCoord(1386, 2) =7.506524072180082e-04;
pointCoord(1387, 0) = 1.638092156936232e-01; pointCoord(1387, 1) =4.836375433147438e-01; pointCoord(1387, 2) =7.506524072180082e-04;
pointCoord(1388, 0) = 3.638092156936232e-01; pointCoord(1388, 1) =4.836375433147438e-01; pointCoord(1388, 2) =7.506524072180082e-04;
pointCoord(1389, 0) = 1.638092156936232e-01; pointCoord(1389, 1) =6.836375433147439e-01; pointCoord(1389, 2) =7.506524072180082e-04;
pointCoord(1390, 0) = 1.377151541156100e-01; pointCoord(1390, 1) =1.097316048927571e-01; pointCoord(1390, 2) =8.678472663211519e-04;
pointCoord(1391, 0) = 3.377151541156100e-01; pointCoord(1391, 1) =1.097316048927571e-01; pointCoord(1391, 2) =8.678472663211519e-04;
pointCoord(1392, 0) = 5.377151541156100e-01; pointCoord(1392, 1) =1.097316048927571e-01; pointCoord(1392, 2) =8.678472663211519e-04;
pointCoord(1393, 0) = 7.377151541156101e-01; pointCoord(1393, 1) =1.097316048927571e-01; pointCoord(1393, 2) =8.678472663211519e-04;
pointCoord(1394, 0) = 1.377151541156100e-01; pointCoord(1394, 1) =3.097316048927571e-01; pointCoord(1394, 2) =8.678472663211519e-04;
pointCoord(1395, 0) = 3.377151541156100e-01; pointCoord(1395, 1) =3.097316048927571e-01; pointCoord(1395, 2) =8.678472663211519e-04;
pointCoord(1396, 0) = 5.377151541156100e-01; pointCoord(1396, 1) =3.097316048927571e-01; pointCoord(1396, 2) =8.678472663211519e-04;
pointCoord(1397, 0) = 1.377151541156100e-01; pointCoord(1397, 1) =5.097316048927572e-01; pointCoord(1397, 2) =8.678472663211519e-04;
pointCoord(1398, 0) = 3.377151541156100e-01; pointCoord(1398, 1) =5.097316048927572e-01; pointCoord(1398, 2) =8.678472663211519e-04;
pointCoord(1399, 0) = 1.377151541156100e-01; pointCoord(1399, 1) =7.097316048927571e-01; pointCoord(1399, 2) =8.678472663211519e-04;
pointCoord(1400, 0) = 1.097316048927571e-01; pointCoord(1400, 1) =1.377151541156100e-01; pointCoord(1400, 2) =8.678472663211519e-04;
pointCoord(1401, 0) = 3.097316048927571e-01; pointCoord(1401, 1) =1.377151541156100e-01; pointCoord(1401, 2) =8.678472663211519e-04;
pointCoord(1402, 0) = 5.097316048927572e-01; pointCoord(1402, 1) =1.377151541156100e-01; pointCoord(1402, 2) =8.678472663211519e-04;
pointCoord(1403, 0) = 7.097316048927571e-01; pointCoord(1403, 1) =1.377151541156100e-01; pointCoord(1403, 2) =8.678472663211519e-04;
pointCoord(1404, 0) = 1.097316048927571e-01; pointCoord(1404, 1) =3.377151541156100e-01; pointCoord(1404, 2) =8.678472663211519e-04;
pointCoord(1405, 0) = 3.097316048927571e-01; pointCoord(1405, 1) =3.377151541156100e-01; pointCoord(1405, 2) =8.678472663211519e-04;
pointCoord(1406, 0) = 5.097316048927572e-01; pointCoord(1406, 1) =3.377151541156100e-01; pointCoord(1406, 2) =8.678472663211519e-04;
pointCoord(1407, 0) = 1.097316048927571e-01; pointCoord(1407, 1) =5.377151541156100e-01; pointCoord(1407, 2) =8.678472663211519e-04;
pointCoord(1408, 0) = 3.097316048927571e-01; pointCoord(1408, 1) =5.377151541156100e-01; pointCoord(1408, 2) =8.678472663211519e-04;
pointCoord(1409, 0) = 1.097316048927571e-01; pointCoord(1409, 1) =7.377151541156101e-01; pointCoord(1409, 2) =8.678472663211519e-04;
pointCoord(1410, 0) = 8.363754331474391e-02; pointCoord(1410, 1) =1.638092156936232e-01; pointCoord(1410, 2) =7.506524072180082e-04;
pointCoord(1411, 0) = 2.836375433147439e-01; pointCoord(1411, 1) =1.638092156936232e-01; pointCoord(1411, 2) =7.506524072180082e-04;
pointCoord(1412, 0) = 4.836375433147438e-01; pointCoord(1412, 1) =1.638092156936232e-01; pointCoord(1412, 2) =7.506524072180082e-04;
pointCoord(1413, 0) = 6.836375433147439e-01; pointCoord(1413, 1) =1.638092156936232e-01; pointCoord(1413, 2) =7.506524072180082e-04;
pointCoord(1414, 0) = 8.363754331474391e-02; pointCoord(1414, 1) =3.638092156936232e-01; pointCoord(1414, 2) =7.506524072180082e-04;
pointCoord(1415, 0) = 2.836375433147439e-01; pointCoord(1415, 1) =3.638092156936232e-01; pointCoord(1415, 2) =7.506524072180082e-04;
pointCoord(1416, 0) = 4.836375433147438e-01; pointCoord(1416, 1) =3.638092156936232e-01; pointCoord(1416, 2) =7.506524072180082e-04;
pointCoord(1417, 0) = 8.363754331474391e-02; pointCoord(1417, 1) =5.638092156936232e-01; pointCoord(1417, 2) =7.506524072180082e-04;
pointCoord(1418, 0) = 2.836375433147439e-01; pointCoord(1418, 1) =5.638092156936232e-01; pointCoord(1418, 2) =7.506524072180082e-04;
pointCoord(1419, 0) = 8.363754331474391e-02; pointCoord(1419, 1) =7.638092156936233e-01; pointCoord(1419, 2) =7.506524072180082e-04;
pointCoord(1420, 0) = 6.295635294476543e-02; pointCoord(1420, 1) =1.844904060636017e-01; pointCoord(1420, 2) =5.321240752324881e-04;
pointCoord(1421, 0) = 2.629563529447654e-01; pointCoord(1421, 1) =1.844904060636017e-01; pointCoord(1421, 2) =5.321240752324881e-04;
pointCoord(1422, 0) = 4.629563529447654e-01; pointCoord(1422, 1) =1.844904060636017e-01; pointCoord(1422, 2) =5.321240752324881e-04;
pointCoord(1423, 0) = 6.629563529447655e-01; pointCoord(1423, 1) =1.844904060636017e-01; pointCoord(1423, 2) =5.321240752324881e-04;
pointCoord(1424, 0) = 6.295635294476543e-02; pointCoord(1424, 1) =3.844904060636017e-01; pointCoord(1424, 2) =5.321240752324881e-04;
pointCoord(1425, 0) = 2.629563529447654e-01; pointCoord(1425, 1) =3.844904060636017e-01; pointCoord(1425, 2) =5.321240752324881e-04;
pointCoord(1426, 0) = 4.629563529447654e-01; pointCoord(1426, 1) =3.844904060636017e-01; pointCoord(1426, 2) =5.321240752324881e-04;
pointCoord(1427, 0) = 6.295635294476543e-02; pointCoord(1427, 1) =5.844904060636017e-01; pointCoord(1427, 2) =5.321240752324881e-04;
pointCoord(1428, 0) = 2.629563529447654e-01; pointCoord(1428, 1) =5.844904060636017e-01; pointCoord(1428, 2) =5.321240752324881e-04;
pointCoord(1429, 0) = 6.295635294476543e-02; pointCoord(1429, 1) =7.844904060636018e-01; pointCoord(1429, 2) =5.321240752324881e-04;
pointCoord(1430, 0) = 5.047571455413893e-02; pointCoord(1430, 1) =1.969710444542282e-01; pointCoord(1430, 2) =2.422245286926613e-04;
pointCoord(1431, 0) = 2.504757145541390e-01; pointCoord(1431, 1) =1.969710444542282e-01; pointCoord(1431, 2) =2.422245286926613e-04;
pointCoord(1432, 0) = 4.504757145541389e-01; pointCoord(1432, 1) =1.969710444542282e-01; pointCoord(1432, 2) =2.422245286926613e-04;
pointCoord(1433, 0) = 6.504757145541390e-01; pointCoord(1433, 1) =1.969710444542282e-01; pointCoord(1433, 2) =2.422245286926613e-04;
pointCoord(1434, 0) = 5.047571455413893e-02; pointCoord(1434, 1) =3.969710444542282e-01; pointCoord(1434, 2) =2.422245286926613e-04;
pointCoord(1435, 0) = 2.504757145541390e-01; pointCoord(1435, 1) =3.969710444542282e-01; pointCoord(1435, 2) =2.422245286926613e-04;
pointCoord(1436, 0) = 4.504757145541389e-01; pointCoord(1436, 1) =3.969710444542282e-01; pointCoord(1436, 2) =2.422245286926613e-04;
pointCoord(1437, 0) = 5.047571455413893e-02; pointCoord(1437, 1) =5.969710444542281e-01; pointCoord(1437, 2) =2.422245286926613e-04;
pointCoord(1438, 0) = 2.504757145541390e-01; pointCoord(1438, 1) =5.969710444542281e-01; pointCoord(1438, 2) =2.422245286926613e-04;
pointCoord(1439, 0) = 5.047571455413893e-02; pointCoord(1439, 1) =7.969710444542282e-01; pointCoord(1439, 2) =2.422245286926613e-04;
pointCoord(1440, 0) = 1.964327058177920e-01; pointCoord(1440, 1) =2.390064644084539e-02; pointCoord(1440, 2) =2.022265498028209e-04;
pointCoord(1441, 0) = 3.964327058177919e-01; pointCoord(1441, 1) =2.390064644084539e-02; pointCoord(1441, 2) =2.022265498028209e-04;
pointCoord(1442, 0) = 5.964327058177920e-01; pointCoord(1442, 1) =2.390064644084539e-02; pointCoord(1442, 2) =2.022265498028209e-04;
pointCoord(1443, 0) = 7.964327058177920e-01; pointCoord(1443, 1) =2.390064644084539e-02; pointCoord(1443, 2) =2.022265498028209e-04;
pointCoord(1444, 0) = 1.964327058177920e-01; pointCoord(1444, 1) =2.239006464408454e-01; pointCoord(1444, 2) =2.022265498028209e-04;
pointCoord(1445, 0) = 3.964327058177919e-01; pointCoord(1445, 1) =2.239006464408454e-01; pointCoord(1445, 2) =2.022265498028209e-04;
pointCoord(1446, 0) = 5.964327058177920e-01; pointCoord(1446, 1) =2.239006464408454e-01; pointCoord(1446, 2) =2.022265498028209e-04;
pointCoord(1447, 0) = 1.964327058177920e-01; pointCoord(1447, 1) =4.239006464408454e-01; pointCoord(1447, 2) =2.022265498028209e-04;
pointCoord(1448, 0) = 3.964327058177919e-01; pointCoord(1448, 1) =4.239006464408454e-01; pointCoord(1448, 2) =2.022265498028209e-04;
pointCoord(1449, 0) = 1.964327058177920e-01; pointCoord(1449, 1) =6.239006464408454e-01; pointCoord(1449, 2) =2.022265498028209e-04;
pointCoord(1450, 0) = 1.817338738117318e-01; pointCoord(1450, 1) =3.859947844690551e-02; pointCoord(1450, 2) =4.442556514902738e-04;
pointCoord(1451, 0) = 3.817338738117318e-01; pointCoord(1451, 1) =3.859947844690551e-02; pointCoord(1451, 2) =4.442556514902738e-04;
pointCoord(1452, 0) = 5.817338738117318e-01; pointCoord(1452, 1) =3.859947844690551e-02; pointCoord(1452, 2) =4.442556514902738e-04;
pointCoord(1453, 0) = 7.817338738117319e-01; pointCoord(1453, 1) =3.859947844690551e-02; pointCoord(1453, 2) =4.442556514902738e-04;
pointCoord(1454, 0) = 1.817338738117318e-01; pointCoord(1454, 1) =2.385994784469055e-01; pointCoord(1454, 2) =4.442556514902738e-04;
pointCoord(1455, 0) = 3.817338738117318e-01; pointCoord(1455, 1) =2.385994784469055e-01; pointCoord(1455, 2) =4.442556514902738e-04;
pointCoord(1456, 0) = 5.817338738117318e-01; pointCoord(1456, 1) =2.385994784469055e-01; pointCoord(1456, 2) =4.442556514902738e-04;
pointCoord(1457, 0) = 1.817338738117318e-01; pointCoord(1457, 1) =4.385994784469055e-01; pointCoord(1457, 2) =4.442556514902738e-04;
pointCoord(1458, 0) = 3.817338738117318e-01; pointCoord(1458, 1) =4.385994784469055e-01; pointCoord(1458, 2) =4.442556514902738e-04;
pointCoord(1459, 0) = 1.817338738117318e-01; pointCoord(1459, 1) =6.385994784469056e-01; pointCoord(1459, 2) =4.442556514902738e-04;
pointCoord(1460, 0) = 1.573769993138719e-01; pointCoord(1460, 1) =6.295635294476543e-02; pointCoord(1460, 2) =6.266989030061803e-04;
pointCoord(1461, 0) = 3.573769993138720e-01; pointCoord(1461, 1) =6.295635294476543e-02; pointCoord(1461, 2) =6.266989030061803e-04;
pointCoord(1462, 0) = 5.573769993138719e-01; pointCoord(1462, 1) =6.295635294476543e-02; pointCoord(1462, 2) =6.266989030061803e-04;
pointCoord(1463, 0) = 7.573769993138719e-01; pointCoord(1463, 1) =6.295635294476543e-02; pointCoord(1463, 2) =6.266989030061803e-04;
pointCoord(1464, 0) = 1.573769993138719e-01; pointCoord(1464, 1) =2.629563529447654e-01; pointCoord(1464, 2) =6.266989030061803e-04;
pointCoord(1465, 0) = 3.573769993138720e-01; pointCoord(1465, 1) =2.629563529447654e-01; pointCoord(1465, 2) =6.266989030061803e-04;
pointCoord(1466, 0) = 5.573769993138719e-01; pointCoord(1466, 1) =2.629563529447654e-01; pointCoord(1466, 2) =6.266989030061803e-04;
pointCoord(1467, 0) = 1.573769993138719e-01; pointCoord(1467, 1) =4.629563529447654e-01; pointCoord(1467, 2) =6.266989030061803e-04;
pointCoord(1468, 0) = 3.573769993138720e-01; pointCoord(1468, 1) =4.629563529447654e-01; pointCoord(1468, 2) =6.266989030061803e-04;
pointCoord(1469, 0) = 1.573769993138719e-01; pointCoord(1469, 1) =6.629563529447655e-01; pointCoord(1469, 2) =6.266989030061803e-04;
pointCoord(1470, 0) = 1.266452197777330e-01; pointCoord(1470, 1) =9.368813248090434e-02; pointCoord(1470, 2) =7.245416447754392e-04;
pointCoord(1471, 0) = 3.266452197777330e-01; pointCoord(1471, 1) =9.368813248090434e-02; pointCoord(1471, 2) =7.245416447754392e-04;
pointCoord(1472, 0) = 5.266452197777330e-01; pointCoord(1472, 1) =9.368813248090434e-02; pointCoord(1472, 2) =7.245416447754392e-04;
pointCoord(1473, 0) = 7.266452197777331e-01; pointCoord(1473, 1) =9.368813248090434e-02; pointCoord(1473, 2) =7.245416447754392e-04;
pointCoord(1474, 0) = 1.266452197777330e-01; pointCoord(1474, 1) =2.936881324809044e-01; pointCoord(1474, 2) =7.245416447754392e-04;
pointCoord(1475, 0) = 3.266452197777330e-01; pointCoord(1475, 1) =2.936881324809044e-01; pointCoord(1475, 2) =7.245416447754392e-04;
pointCoord(1476, 0) = 5.266452197777330e-01; pointCoord(1476, 1) =2.936881324809044e-01; pointCoord(1476, 2) =7.245416447754392e-04;
pointCoord(1477, 0) = 1.266452197777330e-01; pointCoord(1477, 1) =4.936881324809043e-01; pointCoord(1477, 2) =7.245416447754392e-04;
pointCoord(1478, 0) = 3.266452197777330e-01; pointCoord(1478, 1) =4.936881324809043e-01; pointCoord(1478, 2) =7.245416447754392e-04;
pointCoord(1479, 0) = 1.266452197777330e-01; pointCoord(1479, 1) =6.936881324809043e-01; pointCoord(1479, 2) =7.245416447754392e-04;
pointCoord(1480, 0) = 9.368813248090434e-02; pointCoord(1480, 1) =1.266452197777330e-01; pointCoord(1480, 2) =7.245416447754392e-04;
pointCoord(1481, 0) = 2.936881324809044e-01; pointCoord(1481, 1) =1.266452197777330e-01; pointCoord(1481, 2) =7.245416447754392e-04;
pointCoord(1482, 0) = 4.936881324809043e-01; pointCoord(1482, 1) =1.266452197777330e-01; pointCoord(1482, 2) =7.245416447754392e-04;
pointCoord(1483, 0) = 6.936881324809043e-01; pointCoord(1483, 1) =1.266452197777330e-01; pointCoord(1483, 2) =7.245416447754392e-04;
pointCoord(1484, 0) = 9.368813248090434e-02; pointCoord(1484, 1) =3.266452197777330e-01; pointCoord(1484, 2) =7.245416447754392e-04;
pointCoord(1485, 0) = 2.936881324809044e-01; pointCoord(1485, 1) =3.266452197777330e-01; pointCoord(1485, 2) =7.245416447754392e-04;
pointCoord(1486, 0) = 4.936881324809043e-01; pointCoord(1486, 1) =3.266452197777330e-01; pointCoord(1486, 2) =7.245416447754392e-04;
pointCoord(1487, 0) = 9.368813248090434e-02; pointCoord(1487, 1) =5.266452197777330e-01; pointCoord(1487, 2) =7.245416447754392e-04;
pointCoord(1488, 0) = 2.936881324809044e-01; pointCoord(1488, 1) =5.266452197777330e-01; pointCoord(1488, 2) =7.245416447754392e-04;
pointCoord(1489, 0) = 9.368813248090434e-02; pointCoord(1489, 1) =7.266452197777331e-01; pointCoord(1489, 2) =7.245416447754392e-04;
pointCoord(1490, 0) = 6.295635294476543e-02; pointCoord(1490, 1) =1.573769993138719e-01; pointCoord(1490, 2) =6.266989030061803e-04;
pointCoord(1491, 0) = 2.629563529447654e-01; pointCoord(1491, 1) =1.573769993138719e-01; pointCoord(1491, 2) =6.266989030061803e-04;
pointCoord(1492, 0) = 4.629563529447654e-01; pointCoord(1492, 1) =1.573769993138719e-01; pointCoord(1492, 2) =6.266989030061803e-04;
pointCoord(1493, 0) = 6.629563529447655e-01; pointCoord(1493, 1) =1.573769993138719e-01; pointCoord(1493, 2) =6.266989030061803e-04;
pointCoord(1494, 0) = 6.295635294476543e-02; pointCoord(1494, 1) =3.573769993138720e-01; pointCoord(1494, 2) =6.266989030061803e-04;
pointCoord(1495, 0) = 2.629563529447654e-01; pointCoord(1495, 1) =3.573769993138720e-01; pointCoord(1495, 2) =6.266989030061803e-04;
pointCoord(1496, 0) = 4.629563529447654e-01; pointCoord(1496, 1) =3.573769993138720e-01; pointCoord(1496, 2) =6.266989030061803e-04;
pointCoord(1497, 0) = 6.295635294476543e-02; pointCoord(1497, 1) =5.573769993138719e-01; pointCoord(1497, 2) =6.266989030061803e-04;
pointCoord(1498, 0) = 2.629563529447654e-01; pointCoord(1498, 1) =5.573769993138719e-01; pointCoord(1498, 2) =6.266989030061803e-04;
pointCoord(1499, 0) = 6.295635294476543e-02; pointCoord(1499, 1) =7.573769993138719e-01; pointCoord(1499, 2) =6.266989030061803e-04;
pointCoord(1500, 0) = 3.859947844690551e-02; pointCoord(1500, 1) =1.817338738117318e-01; pointCoord(1500, 2) =4.442556514902738e-04;
pointCoord(1501, 0) = 2.385994784469055e-01; pointCoord(1501, 1) =1.817338738117318e-01; pointCoord(1501, 2) =4.442556514902738e-04;
pointCoord(1502, 0) = 4.385994784469055e-01; pointCoord(1502, 1) =1.817338738117318e-01; pointCoord(1502, 2) =4.442556514902738e-04;
pointCoord(1503, 0) = 6.385994784469056e-01; pointCoord(1503, 1) =1.817338738117318e-01; pointCoord(1503, 2) =4.442556514902738e-04;
pointCoord(1504, 0) = 3.859947844690551e-02; pointCoord(1504, 1) =3.817338738117318e-01; pointCoord(1504, 2) =4.442556514902738e-04;
pointCoord(1505, 0) = 2.385994784469055e-01; pointCoord(1505, 1) =3.817338738117318e-01; pointCoord(1505, 2) =4.442556514902738e-04;
pointCoord(1506, 0) = 4.385994784469055e-01; pointCoord(1506, 1) =3.817338738117318e-01; pointCoord(1506, 2) =4.442556514902738e-04;
pointCoord(1507, 0) = 3.859947844690551e-02; pointCoord(1507, 1) =5.817338738117318e-01; pointCoord(1507, 2) =4.442556514902738e-04;
pointCoord(1508, 0) = 2.385994784469055e-01; pointCoord(1508, 1) =5.817338738117318e-01; pointCoord(1508, 2) =4.442556514902738e-04;
pointCoord(1509, 0) = 3.859947844690551e-02; pointCoord(1509, 1) =7.817338738117319e-01; pointCoord(1509, 2) =4.442556514902738e-04;
pointCoord(1510, 0) = 2.390064644084539e-02; pointCoord(1510, 1) =1.964327058177920e-01; pointCoord(1510, 2) =2.022265498028209e-04;
pointCoord(1511, 0) = 2.239006464408454e-01; pointCoord(1511, 1) =1.964327058177920e-01; pointCoord(1511, 2) =2.022265498028209e-04;
pointCoord(1512, 0) = 4.239006464408454e-01; pointCoord(1512, 1) =1.964327058177920e-01; pointCoord(1512, 2) =2.022265498028209e-04;
pointCoord(1513, 0) = 6.239006464408454e-01; pointCoord(1513, 1) =1.964327058177920e-01; pointCoord(1513, 2) =2.022265498028209e-04;
pointCoord(1514, 0) = 2.390064644084539e-02; pointCoord(1514, 1) =3.964327058177919e-01; pointCoord(1514, 2) =2.022265498028209e-04;
pointCoord(1515, 0) = 2.239006464408454e-01; pointCoord(1515, 1) =3.964327058177919e-01; pointCoord(1515, 2) =2.022265498028209e-04;
pointCoord(1516, 0) = 4.239006464408454e-01; pointCoord(1516, 1) =3.964327058177919e-01; pointCoord(1516, 2) =2.022265498028209e-04;
pointCoord(1517, 0) = 2.390064644084539e-02; pointCoord(1517, 1) =5.964327058177920e-01; pointCoord(1517, 2) =2.022265498028209e-04;
pointCoord(1518, 0) = 2.239006464408454e-01; pointCoord(1518, 1) =5.964327058177920e-01; pointCoord(1518, 2) =2.022265498028209e-04;
pointCoord(1519, 0) = 2.390064644084539e-02; pointCoord(1519, 1) =7.964327058177920e-01; pointCoord(1519, 2) =2.022265498028209e-04;
pointCoord(1520, 0) = 1.961078304246029e-01; pointCoord(1520, 1) =7.863183925643419e-03; pointCoord(1520, 2) =1.004375733945300e-04;
pointCoord(1521, 0) = 3.961078304246030e-01; pointCoord(1521, 1) =7.863183925643419e-03; pointCoord(1521, 2) =1.004375733945300e-04;
pointCoord(1522, 0) = 5.961078304246029e-01; pointCoord(1522, 1) =7.863183925643419e-03; pointCoord(1522, 2) =1.004375733945300e-04;
pointCoord(1523, 0) = 7.961078304246030e-01; pointCoord(1523, 1) =7.863183925643419e-03; pointCoord(1523, 2) =1.004375733945300e-04;
pointCoord(1524, 0) = 1.961078304246029e-01; pointCoord(1524, 1) =2.078631839256434e-01; pointCoord(1524, 2) =1.004375733945300e-04;
pointCoord(1525, 0) = 3.961078304246030e-01; pointCoord(1525, 1) =2.078631839256434e-01; pointCoord(1525, 2) =1.004375733945300e-04;
pointCoord(1526, 0) = 5.961078304246029e-01; pointCoord(1526, 1) =2.078631839256434e-01; pointCoord(1526, 2) =1.004375733945300e-04;
pointCoord(1527, 0) = 1.961078304246029e-01; pointCoord(1527, 1) =4.078631839256434e-01; pointCoord(1527, 2) =1.004375733945300e-04;
pointCoord(1528, 0) = 3.961078304246030e-01; pointCoord(1528, 1) =4.078631839256434e-01; pointCoord(1528, 2) =1.004375733945300e-04;
pointCoord(1529, 0) = 1.961078304246029e-01; pointCoord(1529, 1) =6.078631839256434e-01; pointCoord(1529, 2) =1.004375733945300e-04;
pointCoord(1530, 0) = 1.800703679094010e-01; pointCoord(1530, 1) =2.390064644084539e-02; pointCoord(1530, 2) =2.206434300837125e-04;
pointCoord(1531, 0) = 3.800703679094010e-01; pointCoord(1531, 1) =2.390064644084539e-02; pointCoord(1531, 2) =2.206434300837125e-04;
pointCoord(1532, 0) = 5.800703679094009e-01; pointCoord(1532, 1) =2.390064644084539e-02; pointCoord(1532, 2) =2.206434300837125e-04;
pointCoord(1533, 0) = 7.800703679094010e-01; pointCoord(1533, 1) =2.390064644084539e-02; pointCoord(1533, 2) =2.206434300837125e-04;
pointCoord(1534, 0) = 1.800703679094010e-01; pointCoord(1534, 1) =2.239006464408454e-01; pointCoord(1534, 2) =2.206434300837125e-04;
pointCoord(1535, 0) = 3.800703679094010e-01; pointCoord(1535, 1) =2.239006464408454e-01; pointCoord(1535, 2) =2.206434300837125e-04;
pointCoord(1536, 0) = 5.800703679094009e-01; pointCoord(1536, 1) =2.239006464408454e-01; pointCoord(1536, 2) =2.206434300837125e-04;
pointCoord(1537, 0) = 1.800703679094010e-01; pointCoord(1537, 1) =4.239006464408454e-01; pointCoord(1537, 2) =2.206434300837125e-04;
pointCoord(1538, 0) = 3.800703679094010e-01; pointCoord(1538, 1) =4.239006464408454e-01; pointCoord(1538, 2) =2.206434300837125e-04;
pointCoord(1539, 0) = 1.800703679094010e-01; pointCoord(1539, 1) =6.239006464408454e-01; pointCoord(1539, 2) =2.206434300837125e-04;
pointCoord(1540, 0) = 1.534952997961074e-01; pointCoord(1540, 1) =5.047571455413893e-02; pointCoord(1540, 2) =3.112554564587477e-04;
pointCoord(1541, 0) = 3.534952997961074e-01; pointCoord(1541, 1) =5.047571455413893e-02; pointCoord(1541, 2) =3.112554564587477e-04;
pointCoord(1542, 0) = 5.534952997961075e-01; pointCoord(1542, 1) =5.047571455413893e-02; pointCoord(1542, 2) =3.112554564587477e-04;
pointCoord(1543, 0) = 7.534952997961075e-01; pointCoord(1543, 1) =5.047571455413893e-02; pointCoord(1543, 2) =3.112554564587477e-04;
pointCoord(1544, 0) = 1.534952997961074e-01; pointCoord(1544, 1) =2.504757145541390e-01; pointCoord(1544, 2) =3.112554564587477e-04;
pointCoord(1545, 0) = 3.534952997961074e-01; pointCoord(1545, 1) =2.504757145541390e-01; pointCoord(1545, 2) =3.112554564587477e-04;
pointCoord(1546, 0) = 5.534952997961075e-01; pointCoord(1546, 1) =2.504757145541390e-01; pointCoord(1546, 2) =3.112554564587477e-04;
pointCoord(1547, 0) = 1.534952997961074e-01; pointCoord(1547, 1) =4.504757145541389e-01; pointCoord(1547, 2) =3.112554564587477e-04;
pointCoord(1548, 0) = 3.534952997961074e-01; pointCoord(1548, 1) =4.504757145541389e-01; pointCoord(1548, 2) =3.112554564587477e-04;
pointCoord(1549, 0) = 1.534952997961074e-01; pointCoord(1549, 1) =6.504757145541390e-01; pointCoord(1549, 2) =3.112554564587477e-04;
pointCoord(1550, 0) = 1.199647606258469e-01; pointCoord(1550, 1) =8.400625372439949e-02; pointCoord(1550, 2) =3.598499044536019e-04;
pointCoord(1551, 0) = 3.199647606258469e-01; pointCoord(1551, 1) =8.400625372439949e-02; pointCoord(1551, 2) =3.598499044536019e-04;
pointCoord(1552, 0) = 5.199647606258468e-01; pointCoord(1552, 1) =8.400625372439949e-02; pointCoord(1552, 2) =3.598499044536019e-04;
pointCoord(1553, 0) = 7.199647606258469e-01; pointCoord(1553, 1) =8.400625372439949e-02; pointCoord(1553, 2) =3.598499044536019e-04;
pointCoord(1554, 0) = 1.199647606258469e-01; pointCoord(1554, 1) =2.840062537243995e-01; pointCoord(1554, 2) =3.598499044536019e-04;
pointCoord(1555, 0) = 3.199647606258469e-01; pointCoord(1555, 1) =2.840062537243995e-01; pointCoord(1555, 2) =3.598499044536019e-04;
pointCoord(1556, 0) = 5.199647606258468e-01; pointCoord(1556, 1) =2.840062537243995e-01; pointCoord(1556, 2) =3.598499044536019e-04;
pointCoord(1557, 0) = 1.199647606258469e-01; pointCoord(1557, 1) =4.840062537243994e-01; pointCoord(1557, 2) =3.598499044536019e-04;
pointCoord(1558, 0) = 3.199647606258469e-01; pointCoord(1558, 1) =4.840062537243994e-01; pointCoord(1558, 2) =3.598499044536019e-04;
pointCoord(1559, 0) = 1.199647606258469e-01; pointCoord(1559, 1) =6.840062537243995e-01; pointCoord(1559, 2) =3.598499044536019e-04;
pointCoord(1560, 0) = 8.400625372439949e-02; pointCoord(1560, 1) =1.199647606258469e-01; pointCoord(1560, 2) =3.598499044536019e-04;
pointCoord(1561, 0) = 2.840062537243995e-01; pointCoord(1561, 1) =1.199647606258469e-01; pointCoord(1561, 2) =3.598499044536019e-04;
pointCoord(1562, 0) = 4.840062537243994e-01; pointCoord(1562, 1) =1.199647606258469e-01; pointCoord(1562, 2) =3.598499044536019e-04;
pointCoord(1563, 0) = 6.840062537243995e-01; pointCoord(1563, 1) =1.199647606258469e-01; pointCoord(1563, 2) =3.598499044536019e-04;
pointCoord(1564, 0) = 8.400625372439949e-02; pointCoord(1564, 1) =3.199647606258469e-01; pointCoord(1564, 2) =3.598499044536019e-04;
pointCoord(1565, 0) = 2.840062537243995e-01; pointCoord(1565, 1) =3.199647606258469e-01; pointCoord(1565, 2) =3.598499044536019e-04;
pointCoord(1566, 0) = 4.840062537243994e-01; pointCoord(1566, 1) =3.199647606258469e-01; pointCoord(1566, 2) =3.598499044536019e-04;
pointCoord(1567, 0) = 8.400625372439949e-02; pointCoord(1567, 1) =5.199647606258468e-01; pointCoord(1567, 2) =3.598499044536019e-04;
pointCoord(1568, 0) = 2.840062537243995e-01; pointCoord(1568, 1) =5.199647606258468e-01; pointCoord(1568, 2) =3.598499044536019e-04;
pointCoord(1569, 0) = 8.400625372439949e-02; pointCoord(1569, 1) =7.199647606258469e-01; pointCoord(1569, 2) =3.598499044536019e-04;
pointCoord(1570, 0) = 5.047571455413893e-02; pointCoord(1570, 1) =1.534952997961074e-01; pointCoord(1570, 2) =3.112554564587477e-04;
pointCoord(1571, 0) = 2.504757145541390e-01; pointCoord(1571, 1) =1.534952997961074e-01; pointCoord(1571, 2) =3.112554564587477e-04;
pointCoord(1572, 0) = 4.504757145541389e-01; pointCoord(1572, 1) =1.534952997961074e-01; pointCoord(1572, 2) =3.112554564587477e-04;
pointCoord(1573, 0) = 6.504757145541390e-01; pointCoord(1573, 1) =1.534952997961074e-01; pointCoord(1573, 2) =3.112554564587477e-04;
pointCoord(1574, 0) = 5.047571455413893e-02; pointCoord(1574, 1) =3.534952997961074e-01; pointCoord(1574, 2) =3.112554564587477e-04;
pointCoord(1575, 0) = 2.504757145541390e-01; pointCoord(1575, 1) =3.534952997961074e-01; pointCoord(1575, 2) =3.112554564587477e-04;
pointCoord(1576, 0) = 4.504757145541389e-01; pointCoord(1576, 1) =3.534952997961074e-01; pointCoord(1576, 2) =3.112554564587477e-04;
pointCoord(1577, 0) = 5.047571455413893e-02; pointCoord(1577, 1) =5.534952997961075e-01; pointCoord(1577, 2) =3.112554564587477e-04;
pointCoord(1578, 0) = 2.504757145541390e-01; pointCoord(1578, 1) =5.534952997961075e-01; pointCoord(1578, 2) =3.112554564587477e-04;
pointCoord(1579, 0) = 5.047571455413893e-02; pointCoord(1579, 1) =7.534952997961075e-01; pointCoord(1579, 2) =3.112554564587477e-04;
pointCoord(1580, 0) = 2.390064644084539e-02; pointCoord(1580, 1) =1.800703679094010e-01; pointCoord(1580, 2) =2.206434300837125e-04;
pointCoord(1581, 0) = 2.239006464408454e-01; pointCoord(1581, 1) =1.800703679094010e-01; pointCoord(1581, 2) =2.206434300837125e-04;
pointCoord(1582, 0) = 4.239006464408454e-01; pointCoord(1582, 1) =1.800703679094010e-01; pointCoord(1582, 2) =2.206434300837125e-04;
pointCoord(1583, 0) = 6.239006464408454e-01; pointCoord(1583, 1) =1.800703679094010e-01; pointCoord(1583, 2) =2.206434300837125e-04;
pointCoord(1584, 0) = 2.390064644084539e-02; pointCoord(1584, 1) =3.800703679094010e-01; pointCoord(1584, 2) =2.206434300837125e-04;
pointCoord(1585, 0) = 2.239006464408454e-01; pointCoord(1585, 1) =3.800703679094010e-01; pointCoord(1585, 2) =2.206434300837125e-04;
pointCoord(1586, 0) = 4.239006464408454e-01; pointCoord(1586, 1) =3.800703679094010e-01; pointCoord(1586, 2) =2.206434300837125e-04;
pointCoord(1587, 0) = 2.390064644084539e-02; pointCoord(1587, 1) =5.800703679094009e-01; pointCoord(1587, 2) =2.206434300837125e-04;
pointCoord(1588, 0) = 2.239006464408454e-01; pointCoord(1588, 1) =5.800703679094009e-01; pointCoord(1588, 2) =2.206434300837125e-04;
pointCoord(1589, 0) = 2.390064644084539e-02; pointCoord(1589, 1) =7.800703679094010e-01; pointCoord(1589, 2) =2.206434300837125e-04;
pointCoord(1590, 0) = 7.863183925643419e-03; pointCoord(1590, 1) =1.961078304246029e-01; pointCoord(1590, 2) =1.004375733945300e-04;
pointCoord(1591, 0) = 2.078631839256434e-01; pointCoord(1591, 1) =1.961078304246029e-01; pointCoord(1591, 2) =1.004375733945300e-04;
pointCoord(1592, 0) = 4.078631839256434e-01; pointCoord(1592, 1) =1.961078304246029e-01; pointCoord(1592, 2) =1.004375733945300e-04;
pointCoord(1593, 0) = 6.078631839256434e-01; pointCoord(1593, 1) =1.961078304246029e-01; pointCoord(1593, 2) =1.004375733945300e-04;
pointCoord(1594, 0) = 7.863183925643419e-03; pointCoord(1594, 1) =3.961078304246030e-01; pointCoord(1594, 2) =1.004375733945300e-04;
pointCoord(1595, 0) = 2.078631839256434e-01; pointCoord(1595, 1) =3.961078304246030e-01; pointCoord(1595, 2) =1.004375733945300e-04;
pointCoord(1596, 0) = 4.078631839256434e-01; pointCoord(1596, 1) =3.961078304246030e-01; pointCoord(1596, 2) =1.004375733945300e-04;
pointCoord(1597, 0) = 7.863183925643419e-03; pointCoord(1597, 1) =5.961078304246029e-01; pointCoord(1597, 2) =1.004375733945300e-04;
pointCoord(1598, 0) = 2.078631839256434e-01; pointCoord(1598, 1) =5.961078304246029e-01; pointCoord(1598, 2) =1.004375733945300e-04;
pointCoord(1599, 0) = 7.863183925643419e-03; pointCoord(1599, 1) =7.961078304246030e-01; pointCoord(1599, 2) =1.004375733945300e-04;
    }
    else if(nH == 625)
    {
        pointCoord(0, 0) = 4.401110654046414e-04; pointCoord(0, 1) =8.941904340728963e-03; pointCoord(0, 2) =2.633266629202920e-05;
pointCoord(1, 0) = 2.004401110654047e-01; pointCoord(1, 1) =8.941904340728963e-03; pointCoord(1, 2) =2.633266629202920e-05;
pointCoord(2, 0) = 4.004401110654047e-01; pointCoord(2, 1) =8.941904340728963e-03; pointCoord(2, 2) =2.633266629202920e-05;
pointCoord(3, 0) = 6.004401110654046e-01; pointCoord(3, 1) =8.941904340728963e-03; pointCoord(3, 2) =2.633266629202920e-05;
pointCoord(4, 0) = 8.004401110654047e-01; pointCoord(4, 1) =8.941904340728963e-03; pointCoord(4, 2) =2.633266629202920e-05;
pointCoord(5, 0) = 4.401110654046414e-04; pointCoord(5, 1) =2.089419043407290e-01; pointCoord(5, 2) =2.633266629202920e-05;
pointCoord(6, 0) = 2.004401110654047e-01; pointCoord(6, 1) =2.089419043407290e-01; pointCoord(6, 2) =2.633266629202920e-05;
pointCoord(7, 0) = 4.004401110654047e-01; pointCoord(7, 1) =2.089419043407290e-01; pointCoord(7, 2) =2.633266629202920e-05;
pointCoord(8, 0) = 6.004401110654046e-01; pointCoord(8, 1) =2.089419043407290e-01; pointCoord(8, 2) =2.633266629202920e-05;
pointCoord(9, 0) = 4.401110654046414e-04; pointCoord(9, 1) =4.089419043407290e-01; pointCoord(9, 2) =2.633266629202920e-05;
pointCoord(10, 0) = 2.004401110654047e-01; pointCoord(10, 1) =4.089419043407290e-01; pointCoord(10, 2) =2.633266629202920e-05;
pointCoord(11, 0) = 4.004401110654047e-01; pointCoord(11, 1) =4.089419043407290e-01; pointCoord(11, 2) =2.633266629202920e-05;
pointCoord(12, 0) = 4.401110654046414e-04; pointCoord(12, 1) =6.089419043407289e-01; pointCoord(12, 2) =2.633266629202920e-05;
pointCoord(13, 0) = 2.004401110654047e-01; pointCoord(13, 1) =6.089419043407289e-01; pointCoord(13, 2) =2.633266629202920e-05;
pointCoord(14, 0) = 4.401110654046414e-04; pointCoord(14, 1) =8.089419043407290e-01; pointCoord(14, 2) =2.633266629202920e-05;
pointCoord(15, 0) = 2.165044021495976e-03; pointCoord(15, 1) =7.216971384637627e-03; pointCoord(15, 2) =5.319602735277754e-05;
pointCoord(16, 0) = 2.021650440214960e-01; pointCoord(16, 1) =7.216971384637627e-03; pointCoord(16, 2) =5.319602735277754e-05;
pointCoord(17, 0) = 4.021650440214960e-01; pointCoord(17, 1) =7.216971384637627e-03; pointCoord(17, 2) =5.319602735277754e-05;
pointCoord(18, 0) = 6.021650440214960e-01; pointCoord(18, 1) =7.216971384637627e-03; pointCoord(18, 2) =5.319602735277754e-05;
pointCoord(19, 0) = 8.021650440214960e-01; pointCoord(19, 1) =7.216971384637627e-03; pointCoord(19, 2) =5.319602735277754e-05;
pointCoord(20, 0) = 2.165044021495976e-03; pointCoord(20, 1) =2.072169713846376e-01; pointCoord(20, 2) =5.319602735277754e-05;
pointCoord(21, 0) = 2.021650440214960e-01; pointCoord(21, 1) =2.072169713846376e-01; pointCoord(21, 2) =5.319602735277754e-05;
pointCoord(22, 0) = 4.021650440214960e-01; pointCoord(22, 1) =2.072169713846376e-01; pointCoord(22, 2) =5.319602735277754e-05;
pointCoord(23, 0) = 6.021650440214960e-01; pointCoord(23, 1) =2.072169713846376e-01; pointCoord(23, 2) =5.319602735277754e-05;
pointCoord(24, 0) = 2.165044021495976e-03; pointCoord(24, 1) =4.072169713846376e-01; pointCoord(24, 2) =5.319602735277754e-05;
pointCoord(25, 0) = 2.021650440214960e-01; pointCoord(25, 1) =4.072169713846376e-01; pointCoord(25, 2) =5.319602735277754e-05;
pointCoord(26, 0) = 4.021650440214960e-01; pointCoord(26, 1) =4.072169713846376e-01; pointCoord(26, 2) =5.319602735277754e-05;
pointCoord(27, 0) = 2.165044021495976e-03; pointCoord(27, 1) =6.072169713846376e-01; pointCoord(27, 2) =5.319602735277754e-05;
pointCoord(28, 0) = 2.021650440214960e-01; pointCoord(28, 1) =6.072169713846376e-01; pointCoord(28, 2) =5.319602735277754e-05;
pointCoord(29, 0) = 2.165044021495976e-03; pointCoord(29, 1) =8.072169713846377e-01; pointCoord(29, 2) =5.319602735277754e-05;
pointCoord(30, 0) = 4.691007703066801e-03; pointCoord(30, 1) =4.691007703066801e-03; pointCoord(30, 2) =6.322778128282771e-05;
pointCoord(31, 0) = 2.046910077030668e-01; pointCoord(31, 1) =4.691007703066801e-03; pointCoord(31, 2) =6.322778128282771e-05;
pointCoord(32, 0) = 4.046910077030668e-01; pointCoord(32, 1) =4.691007703066801e-03; pointCoord(32, 2) =6.322778128282771e-05;
pointCoord(33, 0) = 6.046910077030668e-01; pointCoord(33, 1) =4.691007703066801e-03; pointCoord(33, 2) =6.322778128282771e-05;
pointCoord(34, 0) = 8.046910077030669e-01; pointCoord(34, 1) =4.691007703066801e-03; pointCoord(34, 2) =6.322778128282771e-05;
pointCoord(35, 0) = 4.691007703066801e-03; pointCoord(35, 1) =2.046910077030668e-01; pointCoord(35, 2) =6.322778128282771e-05;
pointCoord(36, 0) = 2.046910077030668e-01; pointCoord(36, 1) =2.046910077030668e-01; pointCoord(36, 2) =6.322778128282771e-05;
pointCoord(37, 0) = 4.046910077030668e-01; pointCoord(37, 1) =2.046910077030668e-01; pointCoord(37, 2) =6.322778128282771e-05;
pointCoord(38, 0) = 6.046910077030668e-01; pointCoord(38, 1) =2.046910077030668e-01; pointCoord(38, 2) =6.322778128282771e-05;
pointCoord(39, 0) = 4.691007703066801e-03; pointCoord(39, 1) =4.046910077030668e-01; pointCoord(39, 2) =6.322778128282771e-05;
pointCoord(40, 0) = 2.046910077030668e-01; pointCoord(40, 1) =4.046910077030668e-01; pointCoord(40, 2) =6.322778128282771e-05;
pointCoord(41, 0) = 4.046910077030668e-01; pointCoord(41, 1) =4.046910077030668e-01; pointCoord(41, 2) =6.322778128282771e-05;
pointCoord(42, 0) = 4.691007703066801e-03; pointCoord(42, 1) =6.046910077030668e-01; pointCoord(42, 2) =6.322778128282771e-05;
pointCoord(43, 0) = 2.046910077030668e-01; pointCoord(43, 1) =6.046910077030668e-01; pointCoord(43, 2) =6.322778128282771e-05;
pointCoord(44, 0) = 4.691007703066801e-03; pointCoord(44, 1) =8.046910077030669e-01; pointCoord(44, 2) =6.322778128282771e-05;
pointCoord(45, 0) = 7.216971384637627e-03; pointCoord(45, 1) =2.165044021495976e-03; pointCoord(45, 2) =5.319602735277754e-05;
pointCoord(46, 0) = 2.072169713846376e-01; pointCoord(46, 1) =2.165044021495976e-03; pointCoord(46, 2) =5.319602735277754e-05;
pointCoord(47, 0) = 4.072169713846376e-01; pointCoord(47, 1) =2.165044021495976e-03; pointCoord(47, 2) =5.319602735277754e-05;
pointCoord(48, 0) = 6.072169713846376e-01; pointCoord(48, 1) =2.165044021495976e-03; pointCoord(48, 2) =5.319602735277754e-05;
pointCoord(49, 0) = 8.072169713846377e-01; pointCoord(49, 1) =2.165044021495976e-03; pointCoord(49, 2) =5.319602735277754e-05;
pointCoord(50, 0) = 7.216971384637627e-03; pointCoord(50, 1) =2.021650440214960e-01; pointCoord(50, 2) =5.319602735277754e-05;
pointCoord(51, 0) = 2.072169713846376e-01; pointCoord(51, 1) =2.021650440214960e-01; pointCoord(51, 2) =5.319602735277754e-05;
pointCoord(52, 0) = 4.072169713846376e-01; pointCoord(52, 1) =2.021650440214960e-01; pointCoord(52, 2) =5.319602735277754e-05;
pointCoord(53, 0) = 6.072169713846376e-01; pointCoord(53, 1) =2.021650440214960e-01; pointCoord(53, 2) =5.319602735277754e-05;
pointCoord(54, 0) = 7.216971384637627e-03; pointCoord(54, 1) =4.021650440214960e-01; pointCoord(54, 2) =5.319602735277754e-05;
pointCoord(55, 0) = 2.072169713846376e-01; pointCoord(55, 1) =4.021650440214960e-01; pointCoord(55, 2) =5.319602735277754e-05;
pointCoord(56, 0) = 4.072169713846376e-01; pointCoord(56, 1) =4.021650440214960e-01; pointCoord(56, 2) =5.319602735277754e-05;
pointCoord(57, 0) = 7.216971384637627e-03; pointCoord(57, 1) =6.021650440214960e-01; pointCoord(57, 2) =5.319602735277754e-05;
pointCoord(58, 0) = 2.072169713846376e-01; pointCoord(58, 1) =6.021650440214960e-01; pointCoord(58, 2) =5.319602735277754e-05;
pointCoord(59, 0) = 7.216971384637627e-03; pointCoord(59, 1) =8.021650440214960e-01; pointCoord(59, 2) =5.319602735277754e-05;
pointCoord(60, 0) = 8.941904340728963e-03; pointCoord(60, 1) =4.401110654046414e-04; pointCoord(60, 2) =2.633266629202920e-05;
pointCoord(61, 0) = 2.089419043407290e-01; pointCoord(61, 1) =4.401110654046414e-04; pointCoord(61, 2) =2.633266629202920e-05;
pointCoord(62, 0) = 4.089419043407290e-01; pointCoord(62, 1) =4.401110654046414e-04; pointCoord(62, 2) =2.633266629202920e-05;
pointCoord(63, 0) = 6.089419043407289e-01; pointCoord(63, 1) =4.401110654046414e-04; pointCoord(63, 2) =2.633266629202920e-05;
pointCoord(64, 0) = 8.089419043407290e-01; pointCoord(64, 1) =4.401110654046414e-04; pointCoord(64, 2) =2.633266629202920e-05;
pointCoord(65, 0) = 8.941904340728963e-03; pointCoord(65, 1) =2.004401110654047e-01; pointCoord(65, 2) =2.633266629202920e-05;
pointCoord(66, 0) = 2.089419043407290e-01; pointCoord(66, 1) =2.004401110654047e-01; pointCoord(66, 2) =2.633266629202920e-05;
pointCoord(67, 0) = 4.089419043407290e-01; pointCoord(67, 1) =2.004401110654047e-01; pointCoord(67, 2) =2.633266629202920e-05;
pointCoord(68, 0) = 6.089419043407289e-01; pointCoord(68, 1) =2.004401110654047e-01; pointCoord(68, 2) =2.633266629202920e-05;
pointCoord(69, 0) = 8.941904340728963e-03; pointCoord(69, 1) =4.004401110654047e-01; pointCoord(69, 2) =2.633266629202920e-05;
pointCoord(70, 0) = 2.089419043407290e-01; pointCoord(70, 1) =4.004401110654047e-01; pointCoord(70, 2) =2.633266629202920e-05;
pointCoord(71, 0) = 4.089419043407290e-01; pointCoord(71, 1) =4.004401110654047e-01; pointCoord(71, 2) =2.633266629202920e-05;
pointCoord(72, 0) = 8.941904340728963e-03; pointCoord(72, 1) =6.004401110654046e-01; pointCoord(72, 2) =2.633266629202920e-05;
pointCoord(73, 0) = 2.089419043407290e-01; pointCoord(73, 1) =6.004401110654046e-01; pointCoord(73, 2) =2.633266629202920e-05;
pointCoord(74, 0) = 8.941904340728963e-03; pointCoord(74, 1) =8.004401110654047e-01; pointCoord(74, 2) =2.633266629202920e-05;
pointCoord(75, 0) = 2.165044021495976e-03; pointCoord(75, 1) =4.398802496793571e-02; pointCoord(75, 2) =2.616879011700777e-04;
pointCoord(76, 0) = 2.021650440214960e-01; pointCoord(76, 1) =4.398802496793571e-02; pointCoord(76, 2) =2.616879011700777e-04;
pointCoord(77, 0) = 4.021650440214960e-01; pointCoord(77, 1) =4.398802496793571e-02; pointCoord(77, 2) =2.616879011700777e-04;
pointCoord(78, 0) = 6.021650440214960e-01; pointCoord(78, 1) =4.398802496793571e-02; pointCoord(78, 2) =2.616879011700777e-04;
pointCoord(79, 0) = 8.021650440214960e-01; pointCoord(79, 1) =4.398802496793571e-02; pointCoord(79, 2) =2.616879011700777e-04;
pointCoord(80, 0) = 2.165044021495976e-03; pointCoord(80, 1) =2.439880249679357e-01; pointCoord(80, 2) =2.616879011700777e-04;
pointCoord(81, 0) = 2.021650440214960e-01; pointCoord(81, 1) =2.439880249679357e-01; pointCoord(81, 2) =2.616879011700777e-04;
pointCoord(82, 0) = 4.021650440214960e-01; pointCoord(82, 1) =2.439880249679357e-01; pointCoord(82, 2) =2.616879011700777e-04;
pointCoord(83, 0) = 6.021650440214960e-01; pointCoord(83, 1) =2.439880249679357e-01; pointCoord(83, 2) =2.616879011700777e-04;
pointCoord(84, 0) = 2.165044021495976e-03; pointCoord(84, 1) =4.439880249679357e-01; pointCoord(84, 2) =2.616879011700777e-04;
pointCoord(85, 0) = 2.021650440214960e-01; pointCoord(85, 1) =4.439880249679357e-01; pointCoord(85, 2) =2.616879011700777e-04;
pointCoord(86, 0) = 4.021650440214960e-01; pointCoord(86, 1) =4.439880249679357e-01; pointCoord(86, 2) =2.616879011700777e-04;
pointCoord(87, 0) = 2.165044021495976e-03; pointCoord(87, 1) =6.439880249679357e-01; pointCoord(87, 2) =2.616879011700777e-04;
pointCoord(88, 0) = 2.021650440214960e-01; pointCoord(88, 1) =6.439880249679357e-01; pointCoord(88, 2) =2.616879011700777e-04;
pointCoord(89, 0) = 2.165044021495976e-03; pointCoord(89, 1) =8.439880249679358e-01; pointCoord(89, 2) =2.616879011700777e-04;
pointCoord(90, 0) = 1.065052888571620e-02; pointCoord(90, 1) =3.550254010371548e-02; pointCoord(90, 2) =5.286497232810855e-04;
pointCoord(91, 0) = 2.106505288857162e-01; pointCoord(91, 1) =3.550254010371548e-02; pointCoord(91, 2) =5.286497232810855e-04;
pointCoord(92, 0) = 4.106505288857162e-01; pointCoord(92, 1) =3.550254010371548e-02; pointCoord(92, 2) =5.286497232810855e-04;
pointCoord(93, 0) = 6.106505288857161e-01; pointCoord(93, 1) =3.550254010371548e-02; pointCoord(93, 2) =5.286497232810855e-04;
pointCoord(94, 0) = 8.106505288857162e-01; pointCoord(94, 1) =3.550254010371548e-02; pointCoord(94, 2) =5.286497232810855e-04;
pointCoord(95, 0) = 1.065052888571620e-02; pointCoord(95, 1) =2.355025401037155e-01; pointCoord(95, 2) =5.286497232810855e-04;
pointCoord(96, 0) = 2.106505288857162e-01; pointCoord(96, 1) =2.355025401037155e-01; pointCoord(96, 2) =5.286497232810855e-04;
pointCoord(97, 0) = 4.106505288857162e-01; pointCoord(97, 1) =2.355025401037155e-01; pointCoord(97, 2) =5.286497232810855e-04;
pointCoord(98, 0) = 6.106505288857161e-01; pointCoord(98, 1) =2.355025401037155e-01; pointCoord(98, 2) =5.286497232810855e-04;
pointCoord(99, 0) = 1.065052888571620e-02; pointCoord(99, 1) =4.355025401037155e-01; pointCoord(99, 2) =5.286497232810855e-04;
pointCoord(100, 0) = 2.106505288857162e-01; pointCoord(100, 1) =4.355025401037155e-01; pointCoord(100, 2) =5.286497232810855e-04;
pointCoord(101, 0) = 4.106505288857162e-01; pointCoord(101, 1) =4.355025401037155e-01; pointCoord(101, 2) =5.286497232810855e-04;
pointCoord(102, 0) = 1.065052888571620e-02; pointCoord(102, 1) =6.355025401037154e-01; pointCoord(102, 2) =5.286497232810855e-04;
pointCoord(103, 0) = 2.106505288857162e-01; pointCoord(103, 1) =6.355025401037154e-01; pointCoord(103, 2) =5.286497232810855e-04;
pointCoord(104, 0) = 1.065052888571620e-02; pointCoord(104, 1) =8.355025401037155e-01; pointCoord(104, 2) =5.286497232810855e-04;
pointCoord(105, 0) = 2.307653449471584e-02; pointCoord(105, 1) =2.307653449471584e-02; pointCoord(105, 2) =6.283429560853968e-04;
pointCoord(106, 0) = 2.230765344947158e-01; pointCoord(106, 1) =2.307653449471584e-02; pointCoord(106, 2) =6.283429560853968e-04;
pointCoord(107, 0) = 4.230765344947159e-01; pointCoord(107, 1) =2.307653449471584e-02; pointCoord(107, 2) =6.283429560853968e-04;
pointCoord(108, 0) = 6.230765344947158e-01; pointCoord(108, 1) =2.307653449471584e-02; pointCoord(108, 2) =6.283429560853968e-04;
pointCoord(109, 0) = 8.230765344947159e-01; pointCoord(109, 1) =2.307653449471584e-02; pointCoord(109, 2) =6.283429560853968e-04;
pointCoord(110, 0) = 2.307653449471584e-02; pointCoord(110, 1) =2.230765344947158e-01; pointCoord(110, 2) =6.283429560853968e-04;
pointCoord(111, 0) = 2.230765344947158e-01; pointCoord(111, 1) =2.230765344947158e-01; pointCoord(111, 2) =6.283429560853968e-04;
pointCoord(112, 0) = 4.230765344947159e-01; pointCoord(112, 1) =2.230765344947158e-01; pointCoord(112, 2) =6.283429560853968e-04;
pointCoord(113, 0) = 6.230765344947158e-01; pointCoord(113, 1) =2.230765344947158e-01; pointCoord(113, 2) =6.283429560853968e-04;
pointCoord(114, 0) = 2.307653449471584e-02; pointCoord(114, 1) =4.230765344947159e-01; pointCoord(114, 2) =6.283429560853968e-04;
pointCoord(115, 0) = 2.230765344947158e-01; pointCoord(115, 1) =4.230765344947159e-01; pointCoord(115, 2) =6.283429560853968e-04;
pointCoord(116, 0) = 4.230765344947159e-01; pointCoord(116, 1) =4.230765344947159e-01; pointCoord(116, 2) =6.283429560853968e-04;
pointCoord(117, 0) = 2.307653449471584e-02; pointCoord(117, 1) =6.230765344947158e-01; pointCoord(117, 2) =6.283429560853968e-04;
pointCoord(118, 0) = 2.230765344947158e-01; pointCoord(118, 1) =6.230765344947158e-01; pointCoord(118, 2) =6.283429560853968e-04;
pointCoord(119, 0) = 2.307653449471584e-02; pointCoord(119, 1) =8.230765344947159e-01; pointCoord(119, 2) =6.283429560853968e-04;
pointCoord(120, 0) = 3.550254010371548e-02; pointCoord(120, 1) =1.065052888571620e-02; pointCoord(120, 2) =5.286497232810855e-04;
pointCoord(121, 0) = 2.355025401037155e-01; pointCoord(121, 1) =1.065052888571620e-02; pointCoord(121, 2) =5.286497232810855e-04;
pointCoord(122, 0) = 4.355025401037155e-01; pointCoord(122, 1) =1.065052888571620e-02; pointCoord(122, 2) =5.286497232810855e-04;
pointCoord(123, 0) = 6.355025401037154e-01; pointCoord(123, 1) =1.065052888571620e-02; pointCoord(123, 2) =5.286497232810855e-04;
pointCoord(124, 0) = 8.355025401037155e-01; pointCoord(124, 1) =1.065052888571620e-02; pointCoord(124, 2) =5.286497232810855e-04;
pointCoord(125, 0) = 3.550254010371548e-02; pointCoord(125, 1) =2.106505288857162e-01; pointCoord(125, 2) =5.286497232810855e-04;
pointCoord(126, 0) = 2.355025401037155e-01; pointCoord(126, 1) =2.106505288857162e-01; pointCoord(126, 2) =5.286497232810855e-04;
pointCoord(127, 0) = 4.355025401037155e-01; pointCoord(127, 1) =2.106505288857162e-01; pointCoord(127, 2) =5.286497232810855e-04;
pointCoord(128, 0) = 6.355025401037154e-01; pointCoord(128, 1) =2.106505288857162e-01; pointCoord(128, 2) =5.286497232810855e-04;
pointCoord(129, 0) = 3.550254010371548e-02; pointCoord(129, 1) =4.106505288857162e-01; pointCoord(129, 2) =5.286497232810855e-04;
pointCoord(130, 0) = 2.355025401037155e-01; pointCoord(130, 1) =4.106505288857162e-01; pointCoord(130, 2) =5.286497232810855e-04;
pointCoord(131, 0) = 4.355025401037155e-01; pointCoord(131, 1) =4.106505288857162e-01; pointCoord(131, 2) =5.286497232810855e-04;
pointCoord(132, 0) = 3.550254010371548e-02; pointCoord(132, 1) =6.106505288857161e-01; pointCoord(132, 2) =5.286497232810855e-04;
pointCoord(133, 0) = 2.355025401037155e-01; pointCoord(133, 1) =6.106505288857161e-01; pointCoord(133, 2) =5.286497232810855e-04;
pointCoord(134, 0) = 3.550254010371548e-02; pointCoord(134, 1) =8.106505288857162e-01; pointCoord(134, 2) =5.286497232810855e-04;
pointCoord(135, 0) = 4.398802496793571e-02; pointCoord(135, 1) =2.165044021495976e-03; pointCoord(135, 2) =2.616879011700777e-04;
pointCoord(136, 0) = 2.439880249679357e-01; pointCoord(136, 1) =2.165044021495976e-03; pointCoord(136, 2) =2.616879011700777e-04;
pointCoord(137, 0) = 4.439880249679357e-01; pointCoord(137, 1) =2.165044021495976e-03; pointCoord(137, 2) =2.616879011700777e-04;
pointCoord(138, 0) = 6.439880249679357e-01; pointCoord(138, 1) =2.165044021495976e-03; pointCoord(138, 2) =2.616879011700777e-04;
pointCoord(139, 0) = 8.439880249679358e-01; pointCoord(139, 1) =2.165044021495976e-03; pointCoord(139, 2) =2.616879011700777e-04;
pointCoord(140, 0) = 4.398802496793571e-02; pointCoord(140, 1) =2.021650440214960e-01; pointCoord(140, 2) =2.616879011700777e-04;
pointCoord(141, 0) = 2.439880249679357e-01; pointCoord(141, 1) =2.021650440214960e-01; pointCoord(141, 2) =2.616879011700777e-04;
pointCoord(142, 0) = 4.439880249679357e-01; pointCoord(142, 1) =2.021650440214960e-01; pointCoord(142, 2) =2.616879011700777e-04;
pointCoord(143, 0) = 6.439880249679357e-01; pointCoord(143, 1) =2.021650440214960e-01; pointCoord(143, 2) =2.616879011700777e-04;
pointCoord(144, 0) = 4.398802496793571e-02; pointCoord(144, 1) =4.021650440214960e-01; pointCoord(144, 2) =2.616879011700777e-04;
pointCoord(145, 0) = 2.439880249679357e-01; pointCoord(145, 1) =4.021650440214960e-01; pointCoord(145, 2) =2.616879011700777e-04;
pointCoord(146, 0) = 4.439880249679357e-01; pointCoord(146, 1) =4.021650440214960e-01; pointCoord(146, 2) =2.616879011700777e-04;
pointCoord(147, 0) = 4.398802496793571e-02; pointCoord(147, 1) =6.021650440214960e-01; pointCoord(147, 2) =2.616879011700777e-04;
pointCoord(148, 0) = 2.439880249679357e-01; pointCoord(148, 1) =6.021650440214960e-01; pointCoord(148, 2) =2.616879011700777e-04;
pointCoord(149, 0) = 4.398802496793571e-02; pointCoord(149, 1) =8.021650440214960e-01; pointCoord(149, 2) =2.616879011700777e-04;
pointCoord(150, 0) = 4.691007703066801e-03; pointCoord(150, 1) =9.530899229693320e-02; pointCoord(150, 2) =6.739253619376046e-04;
pointCoord(151, 0) = 2.046910077030668e-01; pointCoord(151, 1) =9.530899229693320e-02; pointCoord(151, 2) =6.739253619376046e-04;
pointCoord(152, 0) = 4.046910077030668e-01; pointCoord(152, 1) =9.530899229693320e-02; pointCoord(152, 2) =6.739253619376046e-04;
pointCoord(153, 0) = 6.046910077030668e-01; pointCoord(153, 1) =9.530899229693320e-02; pointCoord(153, 2) =6.739253619376046e-04;
pointCoord(154, 0) = 8.046910077030669e-01; pointCoord(154, 1) =9.530899229693320e-02; pointCoord(154, 2) =6.739253619376046e-04;
pointCoord(155, 0) = 4.691007703066801e-03; pointCoord(155, 1) =2.953089922969332e-01; pointCoord(155, 2) =6.739253619376046e-04;
pointCoord(156, 0) = 2.046910077030668e-01; pointCoord(156, 1) =2.953089922969332e-01; pointCoord(156, 2) =6.739253619376046e-04;
pointCoord(157, 0) = 4.046910077030668e-01; pointCoord(157, 1) =2.953089922969332e-01; pointCoord(157, 2) =6.739253619376046e-04;
pointCoord(158, 0) = 6.046910077030668e-01; pointCoord(158, 1) =2.953089922969332e-01; pointCoord(158, 2) =6.739253619376046e-04;
pointCoord(159, 0) = 4.691007703066801e-03; pointCoord(159, 1) =4.953089922969332e-01; pointCoord(159, 2) =6.739253619376046e-04;
pointCoord(160, 0) = 2.046910077030668e-01; pointCoord(160, 1) =4.953089922969332e-01; pointCoord(160, 2) =6.739253619376046e-04;
pointCoord(161, 0) = 4.046910077030668e-01; pointCoord(161, 1) =4.953089922969332e-01; pointCoord(161, 2) =6.739253619376046e-04;
pointCoord(162, 0) = 4.691007703066801e-03; pointCoord(162, 1) =6.953089922969332e-01; pointCoord(162, 2) =6.739253619376046e-04;
pointCoord(163, 0) = 2.046910077030668e-01; pointCoord(163, 1) =6.953089922969332e-01; pointCoord(163, 2) =6.739253619376046e-04;
pointCoord(164, 0) = 4.691007703066801e-03; pointCoord(164, 1) =8.953089922969333e-01; pointCoord(164, 2) =6.739253619376046e-04;
pointCoord(165, 0) = 2.307653449471584e-02; pointCoord(165, 1) =7.692346550528414e-02; pointCoord(165, 2) =1.361432662753753e-03;
pointCoord(166, 0) = 2.230765344947158e-01; pointCoord(166, 1) =7.692346550528414e-02; pointCoord(166, 2) =1.361432662753753e-03;
pointCoord(167, 0) = 4.230765344947159e-01; pointCoord(167, 1) =7.692346550528414e-02; pointCoord(167, 2) =1.361432662753753e-03;
pointCoord(168, 0) = 6.230765344947158e-01; pointCoord(168, 1) =7.692346550528414e-02; pointCoord(168, 2) =1.361432662753753e-03;
pointCoord(169, 0) = 8.230765344947159e-01; pointCoord(169, 1) =7.692346550528414e-02; pointCoord(169, 2) =1.361432662753753e-03;
pointCoord(170, 0) = 2.307653449471584e-02; pointCoord(170, 1) =2.769234655052841e-01; pointCoord(170, 2) =1.361432662753753e-03;
pointCoord(171, 0) = 2.230765344947158e-01; pointCoord(171, 1) =2.769234655052841e-01; pointCoord(171, 2) =1.361432662753753e-03;
pointCoord(172, 0) = 4.230765344947159e-01; pointCoord(172, 1) =2.769234655052841e-01; pointCoord(172, 2) =1.361432662753753e-03;
pointCoord(173, 0) = 6.230765344947158e-01; pointCoord(173, 1) =2.769234655052841e-01; pointCoord(173, 2) =1.361432662753753e-03;
pointCoord(174, 0) = 2.307653449471584e-02; pointCoord(174, 1) =4.769234655052842e-01; pointCoord(174, 2) =1.361432662753753e-03;
pointCoord(175, 0) = 2.230765344947158e-01; pointCoord(175, 1) =4.769234655052842e-01; pointCoord(175, 2) =1.361432662753753e-03;
pointCoord(176, 0) = 4.230765344947159e-01; pointCoord(176, 1) =4.769234655052842e-01; pointCoord(176, 2) =1.361432662753753e-03;
pointCoord(177, 0) = 2.307653449471584e-02; pointCoord(177, 1) =6.769234655052841e-01; pointCoord(177, 2) =1.361432662753753e-03;
pointCoord(178, 0) = 2.230765344947158e-01; pointCoord(178, 1) =6.769234655052841e-01; pointCoord(178, 2) =1.361432662753753e-03;
pointCoord(179, 0) = 2.307653449471584e-02; pointCoord(179, 1) =8.769234655052842e-01; pointCoord(179, 2) =1.361432662753753e-03;
pointCoord(180, 0) = 5.000000000000000e-02; pointCoord(180, 1) =5.000000000000000e-02; pointCoord(180, 2) =1.618172839506173e-03;
pointCoord(181, 0) = 2.500000000000000e-01; pointCoord(181, 1) =5.000000000000000e-02; pointCoord(181, 2) =1.618172839506173e-03;
pointCoord(182, 0) = 4.500000000000000e-01; pointCoord(182, 1) =5.000000000000000e-02; pointCoord(182, 2) =1.618172839506173e-03;
pointCoord(183, 0) = 6.500000000000000e-01; pointCoord(183, 1) =5.000000000000000e-02; pointCoord(183, 2) =1.618172839506173e-03;
pointCoord(184, 0) = 8.500000000000001e-01; pointCoord(184, 1) =5.000000000000000e-02; pointCoord(184, 2) =1.618172839506173e-03;
pointCoord(185, 0) = 5.000000000000000e-02; pointCoord(185, 1) =2.500000000000000e-01; pointCoord(185, 2) =1.618172839506173e-03;
pointCoord(186, 0) = 2.500000000000000e-01; pointCoord(186, 1) =2.500000000000000e-01; pointCoord(186, 2) =1.618172839506173e-03;
pointCoord(187, 0) = 4.500000000000000e-01; pointCoord(187, 1) =2.500000000000000e-01; pointCoord(187, 2) =1.618172839506173e-03;
pointCoord(188, 0) = 6.500000000000000e-01; pointCoord(188, 1) =2.500000000000000e-01; pointCoord(188, 2) =1.618172839506173e-03;
pointCoord(189, 0) = 5.000000000000000e-02; pointCoord(189, 1) =4.500000000000000e-01; pointCoord(189, 2) =1.618172839506173e-03;
pointCoord(190, 0) = 2.500000000000000e-01; pointCoord(190, 1) =4.500000000000000e-01; pointCoord(190, 2) =1.618172839506173e-03;
pointCoord(191, 0) = 4.500000000000000e-01; pointCoord(191, 1) =4.500000000000000e-01; pointCoord(191, 2) =1.618172839506173e-03;
pointCoord(192, 0) = 5.000000000000000e-02; pointCoord(192, 1) =6.500000000000000e-01; pointCoord(192, 2) =1.618172839506173e-03;
pointCoord(193, 0) = 2.500000000000000e-01; pointCoord(193, 1) =6.500000000000000e-01; pointCoord(193, 2) =1.618172839506173e-03;
pointCoord(194, 0) = 5.000000000000000e-02; pointCoord(194, 1) =8.500000000000001e-01; pointCoord(194, 2) =1.618172839506173e-03;
pointCoord(195, 0) = 7.692346550528414e-02; pointCoord(195, 1) =2.307653449471584e-02; pointCoord(195, 2) =1.361432662753753e-03;
pointCoord(196, 0) = 2.769234655052841e-01; pointCoord(196, 1) =2.307653449471584e-02; pointCoord(196, 2) =1.361432662753753e-03;
pointCoord(197, 0) = 4.769234655052842e-01; pointCoord(197, 1) =2.307653449471584e-02; pointCoord(197, 2) =1.361432662753753e-03;
pointCoord(198, 0) = 6.769234655052841e-01; pointCoord(198, 1) =2.307653449471584e-02; pointCoord(198, 2) =1.361432662753753e-03;
pointCoord(199, 0) = 8.769234655052842e-01; pointCoord(199, 1) =2.307653449471584e-02; pointCoord(199, 2) =1.361432662753753e-03;
pointCoord(200, 0) = 7.692346550528414e-02; pointCoord(200, 1) =2.230765344947158e-01; pointCoord(200, 2) =1.361432662753753e-03;
pointCoord(201, 0) = 2.769234655052841e-01; pointCoord(201, 1) =2.230765344947158e-01; pointCoord(201, 2) =1.361432662753753e-03;
pointCoord(202, 0) = 4.769234655052842e-01; pointCoord(202, 1) =2.230765344947158e-01; pointCoord(202, 2) =1.361432662753753e-03;
pointCoord(203, 0) = 6.769234655052841e-01; pointCoord(203, 1) =2.230765344947158e-01; pointCoord(203, 2) =1.361432662753753e-03;
pointCoord(204, 0) = 7.692346550528414e-02; pointCoord(204, 1) =4.230765344947159e-01; pointCoord(204, 2) =1.361432662753753e-03;
pointCoord(205, 0) = 2.769234655052841e-01; pointCoord(205, 1) =4.230765344947159e-01; pointCoord(205, 2) =1.361432662753753e-03;
pointCoord(206, 0) = 4.769234655052842e-01; pointCoord(206, 1) =4.230765344947159e-01; pointCoord(206, 2) =1.361432662753753e-03;
pointCoord(207, 0) = 7.692346550528414e-02; pointCoord(207, 1) =6.230765344947158e-01; pointCoord(207, 2) =1.361432662753753e-03;
pointCoord(208, 0) = 2.769234655052841e-01; pointCoord(208, 1) =6.230765344947158e-01; pointCoord(208, 2) =1.361432662753753e-03;
pointCoord(209, 0) = 7.692346550528414e-02; pointCoord(209, 1) =8.230765344947159e-01; pointCoord(209, 2) =1.361432662753753e-03;
pointCoord(210, 0) = 9.530899229693320e-02; pointCoord(210, 1) =4.691007703066801e-03; pointCoord(210, 2) =6.739253619376046e-04;
pointCoord(211, 0) = 2.953089922969332e-01; pointCoord(211, 1) =4.691007703066801e-03; pointCoord(211, 2) =6.739253619376046e-04;
pointCoord(212, 0) = 4.953089922969332e-01; pointCoord(212, 1) =4.691007703066801e-03; pointCoord(212, 2) =6.739253619376046e-04;
pointCoord(213, 0) = 6.953089922969332e-01; pointCoord(213, 1) =4.691007703066801e-03; pointCoord(213, 2) =6.739253619376046e-04;
pointCoord(214, 0) = 8.953089922969333e-01; pointCoord(214, 1) =4.691007703066801e-03; pointCoord(214, 2) =6.739253619376046e-04;
pointCoord(215, 0) = 9.530899229693320e-02; pointCoord(215, 1) =2.046910077030668e-01; pointCoord(215, 2) =6.739253619376046e-04;
pointCoord(216, 0) = 2.953089922969332e-01; pointCoord(216, 1) =2.046910077030668e-01; pointCoord(216, 2) =6.739253619376046e-04;
pointCoord(217, 0) = 4.953089922969332e-01; pointCoord(217, 1) =2.046910077030668e-01; pointCoord(217, 2) =6.739253619376046e-04;
pointCoord(218, 0) = 6.953089922969332e-01; pointCoord(218, 1) =2.046910077030668e-01; pointCoord(218, 2) =6.739253619376046e-04;
pointCoord(219, 0) = 9.530899229693320e-02; pointCoord(219, 1) =4.046910077030668e-01; pointCoord(219, 2) =6.739253619376046e-04;
pointCoord(220, 0) = 2.953089922969332e-01; pointCoord(220, 1) =4.046910077030668e-01; pointCoord(220, 2) =6.739253619376046e-04;
pointCoord(221, 0) = 4.953089922969332e-01; pointCoord(221, 1) =4.046910077030668e-01; pointCoord(221, 2) =6.739253619376046e-04;
pointCoord(222, 0) = 9.530899229693320e-02; pointCoord(222, 1) =6.046910077030668e-01; pointCoord(222, 2) =6.739253619376046e-04;
pointCoord(223, 0) = 2.953089922969332e-01; pointCoord(223, 1) =6.046910077030668e-01; pointCoord(223, 2) =6.739253619376046e-04;
pointCoord(224, 0) = 9.530899229693320e-02; pointCoord(224, 1) =8.046910077030669e-01; pointCoord(224, 2) =6.739253619376046e-04;
pointCoord(225, 0) = 7.216971384637627e-03; pointCoord(225, 1) =1.466299596259307e-01; pointCoord(225, 2) =8.723120988299224e-04;
pointCoord(226, 0) = 2.072169713846376e-01; pointCoord(226, 1) =1.466299596259307e-01; pointCoord(226, 2) =8.723120988299224e-04;
pointCoord(227, 0) = 4.072169713846376e-01; pointCoord(227, 1) =1.466299596259307e-01; pointCoord(227, 2) =8.723120988299224e-04;
pointCoord(228, 0) = 6.072169713846376e-01; pointCoord(228, 1) =1.466299596259307e-01; pointCoord(228, 2) =8.723120988299224e-04;
pointCoord(229, 0) = 8.072169713846377e-01; pointCoord(229, 1) =1.466299596259307e-01; pointCoord(229, 2) =8.723120988299224e-04;
pointCoord(230, 0) = 7.216971384637627e-03; pointCoord(230, 1) =3.466299596259307e-01; pointCoord(230, 2) =8.723120988299224e-04;
pointCoord(231, 0) = 2.072169713846376e-01; pointCoord(231, 1) =3.466299596259307e-01; pointCoord(231, 2) =8.723120988299224e-04;
pointCoord(232, 0) = 4.072169713846376e-01; pointCoord(232, 1) =3.466299596259307e-01; pointCoord(232, 2) =8.723120988299224e-04;
pointCoord(233, 0) = 6.072169713846376e-01; pointCoord(233, 1) =3.466299596259307e-01; pointCoord(233, 2) =8.723120988299224e-04;
pointCoord(234, 0) = 7.216971384637627e-03; pointCoord(234, 1) =5.466299596259308e-01; pointCoord(234, 2) =8.723120988299224e-04;
pointCoord(235, 0) = 2.072169713846376e-01; pointCoord(235, 1) =5.466299596259308e-01; pointCoord(235, 2) =8.723120988299224e-04;
pointCoord(236, 0) = 4.072169713846376e-01; pointCoord(236, 1) =5.466299596259308e-01; pointCoord(236, 2) =8.723120988299224e-04;
pointCoord(237, 0) = 7.216971384637627e-03; pointCoord(237, 1) =7.466299596259307e-01; pointCoord(237, 2) =8.723120988299224e-04;
pointCoord(238, 0) = 2.072169713846376e-01; pointCoord(238, 1) =7.466299596259307e-01; pointCoord(238, 2) =8.723120988299224e-04;
pointCoord(239, 0) = 7.216971384637627e-03; pointCoord(239, 1) =9.466299596259307e-01; pointCoord(239, 2) =8.723120988299224e-04;
pointCoord(240, 0) = 3.550254010371548e-02; pointCoord(240, 1) =1.183443909068528e-01; pointCoord(240, 2) =1.762204318958826e-03;
pointCoord(241, 0) = 2.355025401037155e-01; pointCoord(241, 1) =1.183443909068528e-01; pointCoord(241, 2) =1.762204318958826e-03;
pointCoord(242, 0) = 4.355025401037155e-01; pointCoord(242, 1) =1.183443909068528e-01; pointCoord(242, 2) =1.762204318958826e-03;
pointCoord(243, 0) = 6.355025401037154e-01; pointCoord(243, 1) =1.183443909068528e-01; pointCoord(243, 2) =1.762204318958826e-03;
pointCoord(244, 0) = 8.355025401037155e-01; pointCoord(244, 1) =1.183443909068528e-01; pointCoord(244, 2) =1.762204318958826e-03;
pointCoord(245, 0) = 3.550254010371548e-02; pointCoord(245, 1) =3.183443909068528e-01; pointCoord(245, 2) =1.762204318958826e-03;
pointCoord(246, 0) = 2.355025401037155e-01; pointCoord(246, 1) =3.183443909068528e-01; pointCoord(246, 2) =1.762204318958826e-03;
pointCoord(247, 0) = 4.355025401037155e-01; pointCoord(247, 1) =3.183443909068528e-01; pointCoord(247, 2) =1.762204318958826e-03;
pointCoord(248, 0) = 6.355025401037154e-01; pointCoord(248, 1) =3.183443909068528e-01; pointCoord(248, 2) =1.762204318958826e-03;
pointCoord(249, 0) = 3.550254010371548e-02; pointCoord(249, 1) =5.183443909068528e-01; pointCoord(249, 2) =1.762204318958826e-03;
pointCoord(250, 0) = 2.355025401037155e-01; pointCoord(250, 1) =5.183443909068528e-01; pointCoord(250, 2) =1.762204318958826e-03;
pointCoord(251, 0) = 4.355025401037155e-01; pointCoord(251, 1) =5.183443909068528e-01; pointCoord(251, 2) =1.762204318958826e-03;
pointCoord(252, 0) = 3.550254010371548e-02; pointCoord(252, 1) =7.183443909068528e-01; pointCoord(252, 2) =1.762204318958826e-03;
pointCoord(253, 0) = 2.355025401037155e-01; pointCoord(253, 1) =7.183443909068528e-01; pointCoord(253, 2) =1.762204318958826e-03;
pointCoord(254, 0) = 3.550254010371548e-02; pointCoord(254, 1) =9.183443909068528e-01; pointCoord(254, 2) =1.762204318958826e-03;
pointCoord(255, 0) = 7.692346550528414e-02; pointCoord(255, 1) =7.692346550528414e-02; pointCoord(255, 2) =2.094522369422110e-03;
pointCoord(256, 0) = 2.769234655052841e-01; pointCoord(256, 1) =7.692346550528414e-02; pointCoord(256, 2) =2.094522369422110e-03;
pointCoord(257, 0) = 4.769234655052842e-01; pointCoord(257, 1) =7.692346550528414e-02; pointCoord(257, 2) =2.094522369422110e-03;
pointCoord(258, 0) = 6.769234655052841e-01; pointCoord(258, 1) =7.692346550528414e-02; pointCoord(258, 2) =2.094522369422110e-03;
pointCoord(259, 0) = 8.769234655052842e-01; pointCoord(259, 1) =7.692346550528414e-02; pointCoord(259, 2) =2.094522369422110e-03;
pointCoord(260, 0) = 7.692346550528414e-02; pointCoord(260, 1) =2.769234655052841e-01; pointCoord(260, 2) =2.094522369422110e-03;
pointCoord(261, 0) = 2.769234655052841e-01; pointCoord(261, 1) =2.769234655052841e-01; pointCoord(261, 2) =2.094522369422110e-03;
pointCoord(262, 0) = 4.769234655052842e-01; pointCoord(262, 1) =2.769234655052841e-01; pointCoord(262, 2) =2.094522369422110e-03;
pointCoord(263, 0) = 6.769234655052841e-01; pointCoord(263, 1) =2.769234655052841e-01; pointCoord(263, 2) =2.094522369422110e-03;
pointCoord(264, 0) = 7.692346550528414e-02; pointCoord(264, 1) =4.769234655052842e-01; pointCoord(264, 2) =2.094522369422110e-03;
pointCoord(265, 0) = 2.769234655052841e-01; pointCoord(265, 1) =4.769234655052842e-01; pointCoord(265, 2) =2.094522369422110e-03;
pointCoord(266, 0) = 4.769234655052842e-01; pointCoord(266, 1) =4.769234655052842e-01; pointCoord(266, 2) =2.094522369422110e-03;
pointCoord(267, 0) = 7.692346550528414e-02; pointCoord(267, 1) =6.769234655052841e-01; pointCoord(267, 2) =2.094522369422110e-03;
pointCoord(268, 0) = 2.769234655052841e-01; pointCoord(268, 1) =6.769234655052841e-01; pointCoord(268, 2) =2.094522369422110e-03;
pointCoord(269, 0) = 7.692346550528414e-02; pointCoord(269, 1) =8.769234655052842e-01; pointCoord(269, 2) =2.094522369422110e-03;
pointCoord(270, 0) = 1.183443909068528e-01; pointCoord(270, 1) =3.550254010371548e-02; pointCoord(270, 2) =1.762204318958826e-03;
pointCoord(271, 0) = 3.183443909068528e-01; pointCoord(271, 1) =3.550254010371548e-02; pointCoord(271, 2) =1.762204318958826e-03;
pointCoord(272, 0) = 5.183443909068528e-01; pointCoord(272, 1) =3.550254010371548e-02; pointCoord(272, 2) =1.762204318958826e-03;
pointCoord(273, 0) = 7.183443909068528e-01; pointCoord(273, 1) =3.550254010371548e-02; pointCoord(273, 2) =1.762204318958826e-03;
pointCoord(274, 0) = 9.183443909068528e-01; pointCoord(274, 1) =3.550254010371548e-02; pointCoord(274, 2) =1.762204318958826e-03;
pointCoord(275, 0) = 1.183443909068528e-01; pointCoord(275, 1) =2.355025401037155e-01; pointCoord(275, 2) =1.762204318958826e-03;
pointCoord(276, 0) = 3.183443909068528e-01; pointCoord(276, 1) =2.355025401037155e-01; pointCoord(276, 2) =1.762204318958826e-03;
pointCoord(277, 0) = 5.183443909068528e-01; pointCoord(277, 1) =2.355025401037155e-01; pointCoord(277, 2) =1.762204318958826e-03;
pointCoord(278, 0) = 7.183443909068528e-01; pointCoord(278, 1) =2.355025401037155e-01; pointCoord(278, 2) =1.762204318958826e-03;
pointCoord(279, 0) = 1.183443909068528e-01; pointCoord(279, 1) =4.355025401037155e-01; pointCoord(279, 2) =1.762204318958826e-03;
pointCoord(280, 0) = 3.183443909068528e-01; pointCoord(280, 1) =4.355025401037155e-01; pointCoord(280, 2) =1.762204318958826e-03;
pointCoord(281, 0) = 5.183443909068528e-01; pointCoord(281, 1) =4.355025401037155e-01; pointCoord(281, 2) =1.762204318958826e-03;
pointCoord(282, 0) = 1.183443909068528e-01; pointCoord(282, 1) =6.355025401037154e-01; pointCoord(282, 2) =1.762204318958826e-03;
pointCoord(283, 0) = 3.183443909068528e-01; pointCoord(283, 1) =6.355025401037154e-01; pointCoord(283, 2) =1.762204318958826e-03;
pointCoord(284, 0) = 1.183443909068528e-01; pointCoord(284, 1) =8.355025401037155e-01; pointCoord(284, 2) =1.762204318958826e-03;
pointCoord(285, 0) = 1.466299596259307e-01; pointCoord(285, 1) =7.216971384637627e-03; pointCoord(285, 2) =8.723120988299224e-04;
pointCoord(286, 0) = 3.466299596259307e-01; pointCoord(286, 1) =7.216971384637627e-03; pointCoord(286, 2) =8.723120988299224e-04;
pointCoord(287, 0) = 5.466299596259308e-01; pointCoord(287, 1) =7.216971384637627e-03; pointCoord(287, 2) =8.723120988299224e-04;
pointCoord(288, 0) = 7.466299596259307e-01; pointCoord(288, 1) =7.216971384637627e-03; pointCoord(288, 2) =8.723120988299224e-04;
pointCoord(289, 0) = 9.466299596259307e-01; pointCoord(289, 1) =7.216971384637627e-03; pointCoord(289, 2) =8.723120988299224e-04;
pointCoord(290, 0) = 1.466299596259307e-01; pointCoord(290, 1) =2.072169713846376e-01; pointCoord(290, 2) =8.723120988299224e-04;
pointCoord(291, 0) = 3.466299596259307e-01; pointCoord(291, 1) =2.072169713846376e-01; pointCoord(291, 2) =8.723120988299224e-04;
pointCoord(292, 0) = 5.466299596259308e-01; pointCoord(292, 1) =2.072169713846376e-01; pointCoord(292, 2) =8.723120988299224e-04;
pointCoord(293, 0) = 7.466299596259307e-01; pointCoord(293, 1) =2.072169713846376e-01; pointCoord(293, 2) =8.723120988299224e-04;
pointCoord(294, 0) = 1.466299596259307e-01; pointCoord(294, 1) =4.072169713846376e-01; pointCoord(294, 2) =8.723120988299224e-04;
pointCoord(295, 0) = 3.466299596259307e-01; pointCoord(295, 1) =4.072169713846376e-01; pointCoord(295, 2) =8.723120988299224e-04;
pointCoord(296, 0) = 5.466299596259308e-01; pointCoord(296, 1) =4.072169713846376e-01; pointCoord(296, 2) =8.723120988299224e-04;
pointCoord(297, 0) = 1.466299596259307e-01; pointCoord(297, 1) =6.072169713846376e-01; pointCoord(297, 2) =8.723120988299224e-04;
pointCoord(298, 0) = 3.466299596259307e-01; pointCoord(298, 1) =6.072169713846376e-01; pointCoord(298, 2) =8.723120988299224e-04;
pointCoord(299, 0) = 1.466299596259307e-01; pointCoord(299, 1) =8.072169713846377e-01; pointCoord(299, 2) =8.723120988299224e-04;
pointCoord(300, 0) = 8.941904340728963e-03; pointCoord(300, 1) =1.816760802531374e-01; pointCoord(300, 2) =5.350108223322572e-04;
pointCoord(301, 0) = 2.089419043407290e-01; pointCoord(301, 1) =1.816760802531374e-01; pointCoord(301, 2) =5.350108223322572e-04;
pointCoord(302, 0) = 4.089419043407290e-01; pointCoord(302, 1) =1.816760802531374e-01; pointCoord(302, 2) =5.350108223322572e-04;
pointCoord(303, 0) = 6.089419043407289e-01; pointCoord(303, 1) =1.816760802531374e-01; pointCoord(303, 2) =5.350108223322572e-04;
pointCoord(304, 0) = 8.089419043407290e-01; pointCoord(304, 1) =1.816760802531374e-01; pointCoord(304, 2) =5.350108223322572e-04;
pointCoord(305, 0) = 8.941904340728963e-03; pointCoord(305, 1) =3.816760802531374e-01; pointCoord(305, 2) =5.350108223322572e-04;
pointCoord(306, 0) = 2.089419043407290e-01; pointCoord(306, 1) =3.816760802531374e-01; pointCoord(306, 2) =5.350108223322572e-04;
pointCoord(307, 0) = 4.089419043407290e-01; pointCoord(307, 1) =3.816760802531374e-01; pointCoord(307, 2) =5.350108223322572e-04;
pointCoord(308, 0) = 6.089419043407289e-01; pointCoord(308, 1) =3.816760802531374e-01; pointCoord(308, 2) =5.350108223322572e-04;
pointCoord(309, 0) = 8.941904340728963e-03; pointCoord(309, 1) =5.816760802531374e-01; pointCoord(309, 2) =5.350108223322572e-04;
pointCoord(310, 0) = 2.089419043407290e-01; pointCoord(310, 1) =5.816760802531374e-01; pointCoord(310, 2) =5.350108223322572e-04;
pointCoord(311, 0) = 4.089419043407290e-01; pointCoord(311, 1) =5.816760802531374e-01; pointCoord(311, 2) =5.350108223322572e-04;
pointCoord(312, 0) = 8.941904340728963e-03; pointCoord(312, 1) =7.816760802531374e-01; pointCoord(312, 2) =5.350108223322572e-04;
pointCoord(313, 0) = 2.089419043407290e-01; pointCoord(313, 1) =7.816760802531374e-01; pointCoord(313, 2) =5.350108223322572e-04;
pointCoord(314, 0) = 8.941904340728963e-03; pointCoord(314, 1) =9.816760802531375e-01; pointCoord(314, 2) =5.350108223322572e-04;
pointCoord(315, 0) = 4.398802496793571e-02; pointCoord(315, 1) =1.466299596259307e-01; pointCoord(315, 2) =1.080803972647223e-03;
pointCoord(316, 0) = 2.439880249679357e-01; pointCoord(316, 1) =1.466299596259307e-01; pointCoord(316, 2) =1.080803972647223e-03;
pointCoord(317, 0) = 4.439880249679357e-01; pointCoord(317, 1) =1.466299596259307e-01; pointCoord(317, 2) =1.080803972647223e-03;
pointCoord(318, 0) = 6.439880249679357e-01; pointCoord(318, 1) =1.466299596259307e-01; pointCoord(318, 2) =1.080803972647223e-03;
pointCoord(319, 0) = 8.439880249679358e-01; pointCoord(319, 1) =1.466299596259307e-01; pointCoord(319, 2) =1.080803972647223e-03;
pointCoord(320, 0) = 4.398802496793571e-02; pointCoord(320, 1) =3.466299596259307e-01; pointCoord(320, 2) =1.080803972647223e-03;
pointCoord(321, 0) = 2.439880249679357e-01; pointCoord(321, 1) =3.466299596259307e-01; pointCoord(321, 2) =1.080803972647223e-03;
pointCoord(322, 0) = 4.439880249679357e-01; pointCoord(322, 1) =3.466299596259307e-01; pointCoord(322, 2) =1.080803972647223e-03;
pointCoord(323, 0) = 6.439880249679357e-01; pointCoord(323, 1) =3.466299596259307e-01; pointCoord(323, 2) =1.080803972647223e-03;
pointCoord(324, 0) = 4.398802496793571e-02; pointCoord(324, 1) =5.466299596259308e-01; pointCoord(324, 2) =1.080803972647223e-03;
pointCoord(325, 0) = 2.439880249679357e-01; pointCoord(325, 1) =5.466299596259308e-01; pointCoord(325, 2) =1.080803972647223e-03;
pointCoord(326, 0) = 4.439880249679357e-01; pointCoord(326, 1) =5.466299596259308e-01; pointCoord(326, 2) =1.080803972647223e-03;
pointCoord(327, 0) = 4.398802496793571e-02; pointCoord(327, 1) =7.466299596259307e-01; pointCoord(327, 2) =1.080803972647223e-03;
pointCoord(328, 0) = 2.439880249679357e-01; pointCoord(328, 1) =7.466299596259307e-01; pointCoord(328, 2) =1.080803972647223e-03;
pointCoord(329, 0) = 4.398802496793571e-02; pointCoord(329, 1) =9.466299596259307e-01; pointCoord(329, 2) =1.080803972647223e-03;
pointCoord(330, 0) = 9.530899229693320e-02; pointCoord(330, 1) =9.530899229693320e-02; pointCoord(330, 2) =1.284622942592381e-03;
pointCoord(331, 0) = 2.953089922969332e-01; pointCoord(331, 1) =9.530899229693320e-02; pointCoord(331, 2) =1.284622942592381e-03;
pointCoord(332, 0) = 4.953089922969332e-01; pointCoord(332, 1) =9.530899229693320e-02; pointCoord(332, 2) =1.284622942592381e-03;
pointCoord(333, 0) = 6.953089922969332e-01; pointCoord(333, 1) =9.530899229693320e-02; pointCoord(333, 2) =1.284622942592381e-03;
pointCoord(334, 0) = 8.953089922969333e-01; pointCoord(334, 1) =9.530899229693320e-02; pointCoord(334, 2) =1.284622942592381e-03;
pointCoord(335, 0) = 9.530899229693320e-02; pointCoord(335, 1) =2.953089922969332e-01; pointCoord(335, 2) =1.284622942592381e-03;
pointCoord(336, 0) = 2.953089922969332e-01; pointCoord(336, 1) =2.953089922969332e-01; pointCoord(336, 2) =1.284622942592381e-03;
pointCoord(337, 0) = 4.953089922969332e-01; pointCoord(337, 1) =2.953089922969332e-01; pointCoord(337, 2) =1.284622942592381e-03;
pointCoord(338, 0) = 6.953089922969332e-01; pointCoord(338, 1) =2.953089922969332e-01; pointCoord(338, 2) =1.284622942592381e-03;
pointCoord(339, 0) = 9.530899229693320e-02; pointCoord(339, 1) =4.953089922969332e-01; pointCoord(339, 2) =1.284622942592381e-03;
pointCoord(340, 0) = 2.953089922969332e-01; pointCoord(340, 1) =4.953089922969332e-01; pointCoord(340, 2) =1.284622942592381e-03;
pointCoord(341, 0) = 4.953089922969332e-01; pointCoord(341, 1) =4.953089922969332e-01; pointCoord(341, 2) =1.284622942592381e-03;
pointCoord(342, 0) = 9.530899229693320e-02; pointCoord(342, 1) =6.953089922969332e-01; pointCoord(342, 2) =1.284622942592381e-03;
pointCoord(343, 0) = 2.953089922969332e-01; pointCoord(343, 1) =6.953089922969332e-01; pointCoord(343, 2) =1.284622942592381e-03;
pointCoord(344, 0) = 9.530899229693320e-02; pointCoord(344, 1) =8.953089922969333e-01; pointCoord(344, 2) =1.284622942592381e-03;
pointCoord(345, 0) = 1.466299596259307e-01; pointCoord(345, 1) =4.398802496793571e-02; pointCoord(345, 2) =1.080803972647223e-03;
pointCoord(346, 0) = 3.466299596259307e-01; pointCoord(346, 1) =4.398802496793571e-02; pointCoord(346, 2) =1.080803972647223e-03;
pointCoord(347, 0) = 5.466299596259308e-01; pointCoord(347, 1) =4.398802496793571e-02; pointCoord(347, 2) =1.080803972647223e-03;
pointCoord(348, 0) = 7.466299596259307e-01; pointCoord(348, 1) =4.398802496793571e-02; pointCoord(348, 2) =1.080803972647223e-03;
pointCoord(349, 0) = 9.466299596259307e-01; pointCoord(349, 1) =4.398802496793571e-02; pointCoord(349, 2) =1.080803972647223e-03;
pointCoord(350, 0) = 1.466299596259307e-01; pointCoord(350, 1) =2.439880249679357e-01; pointCoord(350, 2) =1.080803972647223e-03;
pointCoord(351, 0) = 3.466299596259307e-01; pointCoord(351, 1) =2.439880249679357e-01; pointCoord(351, 2) =1.080803972647223e-03;
pointCoord(352, 0) = 5.466299596259308e-01; pointCoord(352, 1) =2.439880249679357e-01; pointCoord(352, 2) =1.080803972647223e-03;
pointCoord(353, 0) = 7.466299596259307e-01; pointCoord(353, 1) =2.439880249679357e-01; pointCoord(353, 2) =1.080803972647223e-03;
pointCoord(354, 0) = 1.466299596259307e-01; pointCoord(354, 1) =4.439880249679357e-01; pointCoord(354, 2) =1.080803972647223e-03;
pointCoord(355, 0) = 3.466299596259307e-01; pointCoord(355, 1) =4.439880249679357e-01; pointCoord(355, 2) =1.080803972647223e-03;
pointCoord(356, 0) = 5.466299596259308e-01; pointCoord(356, 1) =4.439880249679357e-01; pointCoord(356, 2) =1.080803972647223e-03;
pointCoord(357, 0) = 1.466299596259307e-01; pointCoord(357, 1) =6.439880249679357e-01; pointCoord(357, 2) =1.080803972647223e-03;
pointCoord(358, 0) = 3.466299596259307e-01; pointCoord(358, 1) =6.439880249679357e-01; pointCoord(358, 2) =1.080803972647223e-03;
pointCoord(359, 0) = 1.466299596259307e-01; pointCoord(359, 1) =8.439880249679358e-01; pointCoord(359, 2) =1.080803972647223e-03;
pointCoord(360, 0) = 1.816760802531374e-01; pointCoord(360, 1) =8.941904340728963e-03; pointCoord(360, 2) =5.350108223322572e-04;
pointCoord(361, 0) = 3.816760802531374e-01; pointCoord(361, 1) =8.941904340728963e-03; pointCoord(361, 2) =5.350108223322572e-04;
pointCoord(362, 0) = 5.816760802531374e-01; pointCoord(362, 1) =8.941904340728963e-03; pointCoord(362, 2) =5.350108223322572e-04;
pointCoord(363, 0) = 7.816760802531374e-01; pointCoord(363, 1) =8.941904340728963e-03; pointCoord(363, 2) =5.350108223322572e-04;
pointCoord(364, 0) = 9.816760802531375e-01; pointCoord(364, 1) =8.941904340728963e-03; pointCoord(364, 2) =5.350108223322572e-04;
pointCoord(365, 0) = 1.816760802531374e-01; pointCoord(365, 1) =2.089419043407290e-01; pointCoord(365, 2) =5.350108223322572e-04;
pointCoord(366, 0) = 3.816760802531374e-01; pointCoord(366, 1) =2.089419043407290e-01; pointCoord(366, 2) =5.350108223322572e-04;
pointCoord(367, 0) = 5.816760802531374e-01; pointCoord(367, 1) =2.089419043407290e-01; pointCoord(367, 2) =5.350108223322572e-04;
pointCoord(368, 0) = 7.816760802531374e-01; pointCoord(368, 1) =2.089419043407290e-01; pointCoord(368, 2) =5.350108223322572e-04;
pointCoord(369, 0) = 1.816760802531374e-01; pointCoord(369, 1) =4.089419043407290e-01; pointCoord(369, 2) =5.350108223322572e-04;
pointCoord(370, 0) = 3.816760802531374e-01; pointCoord(370, 1) =4.089419043407290e-01; pointCoord(370, 2) =5.350108223322572e-04;
pointCoord(371, 0) = 5.816760802531374e-01; pointCoord(371, 1) =4.089419043407290e-01; pointCoord(371, 2) =5.350108223322572e-04;
pointCoord(372, 0) = 1.816760802531374e-01; pointCoord(372, 1) =6.089419043407289e-01; pointCoord(372, 2) =5.350108223322572e-04;
pointCoord(373, 0) = 3.816760802531374e-01; pointCoord(373, 1) =6.089419043407289e-01; pointCoord(373, 2) =5.350108223322572e-04;
pointCoord(374, 0) = 1.816760802531374e-01; pointCoord(374, 1) =8.089419043407290e-01; pointCoord(374, 2) =5.350108223322572e-04;
pointCoord(375, 0) = 1.995598889345954e-01; pointCoord(375, 1) =1.910580956592710e-01; pointCoord(375, 2) =2.633266629202920e-05;
pointCoord(376, 0) = 3.995598889345954e-01; pointCoord(376, 1) =1.910580956592710e-01; pointCoord(376, 2) =2.633266629202920e-05;
pointCoord(377, 0) = 5.995598889345953e-01; pointCoord(377, 1) =1.910580956592710e-01; pointCoord(377, 2) =2.633266629202920e-05;
pointCoord(378, 0) = 7.995598889345954e-01; pointCoord(378, 1) =1.910580956592710e-01; pointCoord(378, 2) =2.633266629202920e-05;
pointCoord(379, 0) = 1.995598889345954e-01; pointCoord(379, 1) =3.910580956592711e-01; pointCoord(379, 2) =2.633266629202920e-05;
pointCoord(380, 0) = 3.995598889345954e-01; pointCoord(380, 1) =3.910580956592711e-01; pointCoord(380, 2) =2.633266629202920e-05;
pointCoord(381, 0) = 5.995598889345953e-01; pointCoord(381, 1) =3.910580956592711e-01; pointCoord(381, 2) =2.633266629202920e-05;
pointCoord(382, 0) = 1.995598889345954e-01; pointCoord(382, 1) =5.910580956592710e-01; pointCoord(382, 2) =2.633266629202920e-05;
pointCoord(383, 0) = 3.995598889345954e-01; pointCoord(383, 1) =5.910580956592710e-01; pointCoord(383, 2) =2.633266629202920e-05;
pointCoord(384, 0) = 1.995598889345954e-01; pointCoord(384, 1) =7.910580956592711e-01; pointCoord(384, 2) =2.633266629202920e-05;
pointCoord(385, 0) = 1.978349559785040e-01; pointCoord(385, 1) =1.927830286153624e-01; pointCoord(385, 2) =5.319602735277754e-05;
pointCoord(386, 0) = 3.978349559785040e-01; pointCoord(386, 1) =1.927830286153624e-01; pointCoord(386, 2) =5.319602735277754e-05;
pointCoord(387, 0) = 5.978349559785040e-01; pointCoord(387, 1) =1.927830286153624e-01; pointCoord(387, 2) =5.319602735277754e-05;
pointCoord(388, 0) = 7.978349559785041e-01; pointCoord(388, 1) =1.927830286153624e-01; pointCoord(388, 2) =5.319602735277754e-05;
pointCoord(389, 0) = 1.978349559785040e-01; pointCoord(389, 1) =3.927830286153624e-01; pointCoord(389, 2) =5.319602735277754e-05;
pointCoord(390, 0) = 3.978349559785040e-01; pointCoord(390, 1) =3.927830286153624e-01; pointCoord(390, 2) =5.319602735277754e-05;
pointCoord(391, 0) = 5.978349559785040e-01; pointCoord(391, 1) =3.927830286153624e-01; pointCoord(391, 2) =5.319602735277754e-05;
pointCoord(392, 0) = 1.978349559785040e-01; pointCoord(392, 1) =5.927830286153624e-01; pointCoord(392, 2) =5.319602735277754e-05;
pointCoord(393, 0) = 3.978349559785040e-01; pointCoord(393, 1) =5.927830286153624e-01; pointCoord(393, 2) =5.319602735277754e-05;
pointCoord(394, 0) = 1.978349559785040e-01; pointCoord(394, 1) =7.927830286153624e-01; pointCoord(394, 2) =5.319602735277754e-05;
pointCoord(395, 0) = 1.953089922969332e-01; pointCoord(395, 1) =1.953089922969332e-01; pointCoord(395, 2) =6.322778128282771e-05;
pointCoord(396, 0) = 3.953089922969332e-01; pointCoord(396, 1) =1.953089922969332e-01; pointCoord(396, 2) =6.322778128282771e-05;
pointCoord(397, 0) = 5.953089922969331e-01; pointCoord(397, 1) =1.953089922969332e-01; pointCoord(397, 2) =6.322778128282771e-05;
pointCoord(398, 0) = 7.953089922969332e-01; pointCoord(398, 1) =1.953089922969332e-01; pointCoord(398, 2) =6.322778128282771e-05;
pointCoord(399, 0) = 1.953089922969332e-01; pointCoord(399, 1) =3.953089922969332e-01; pointCoord(399, 2) =6.322778128282771e-05;
pointCoord(400, 0) = 3.953089922969332e-01; pointCoord(400, 1) =3.953089922969332e-01; pointCoord(400, 2) =6.322778128282771e-05;
pointCoord(401, 0) = 5.953089922969331e-01; pointCoord(401, 1) =3.953089922969332e-01; pointCoord(401, 2) =6.322778128282771e-05;
pointCoord(402, 0) = 1.953089922969332e-01; pointCoord(402, 1) =5.953089922969331e-01; pointCoord(402, 2) =6.322778128282771e-05;
pointCoord(403, 0) = 3.953089922969332e-01; pointCoord(403, 1) =5.953089922969331e-01; pointCoord(403, 2) =6.322778128282771e-05;
pointCoord(404, 0) = 1.953089922969332e-01; pointCoord(404, 1) =7.953089922969332e-01; pointCoord(404, 2) =6.322778128282771e-05;
pointCoord(405, 0) = 1.927830286153624e-01; pointCoord(405, 1) =1.978349559785040e-01; pointCoord(405, 2) =5.319602735277754e-05;
pointCoord(406, 0) = 3.927830286153624e-01; pointCoord(406, 1) =1.978349559785040e-01; pointCoord(406, 2) =5.319602735277754e-05;
pointCoord(407, 0) = 5.927830286153624e-01; pointCoord(407, 1) =1.978349559785040e-01; pointCoord(407, 2) =5.319602735277754e-05;
pointCoord(408, 0) = 7.927830286153624e-01; pointCoord(408, 1) =1.978349559785040e-01; pointCoord(408, 2) =5.319602735277754e-05;
pointCoord(409, 0) = 1.927830286153624e-01; pointCoord(409, 1) =3.978349559785040e-01; pointCoord(409, 2) =5.319602735277754e-05;
pointCoord(410, 0) = 3.927830286153624e-01; pointCoord(410, 1) =3.978349559785040e-01; pointCoord(410, 2) =5.319602735277754e-05;
pointCoord(411, 0) = 5.927830286153624e-01; pointCoord(411, 1) =3.978349559785040e-01; pointCoord(411, 2) =5.319602735277754e-05;
pointCoord(412, 0) = 1.927830286153624e-01; pointCoord(412, 1) =5.978349559785040e-01; pointCoord(412, 2) =5.319602735277754e-05;
pointCoord(413, 0) = 3.927830286153624e-01; pointCoord(413, 1) =5.978349559785040e-01; pointCoord(413, 2) =5.319602735277754e-05;
pointCoord(414, 0) = 1.927830286153624e-01; pointCoord(414, 1) =7.978349559785041e-01; pointCoord(414, 2) =5.319602735277754e-05;
pointCoord(415, 0) = 1.910580956592710e-01; pointCoord(415, 1) =1.995598889345954e-01; pointCoord(415, 2) =2.633266629202920e-05;
pointCoord(416, 0) = 3.910580956592711e-01; pointCoord(416, 1) =1.995598889345954e-01; pointCoord(416, 2) =2.633266629202920e-05;
pointCoord(417, 0) = 5.910580956592710e-01; pointCoord(417, 1) =1.995598889345954e-01; pointCoord(417, 2) =2.633266629202920e-05;
pointCoord(418, 0) = 7.910580956592711e-01; pointCoord(418, 1) =1.995598889345954e-01; pointCoord(418, 2) =2.633266629202920e-05;
pointCoord(419, 0) = 1.910580956592710e-01; pointCoord(419, 1) =3.995598889345954e-01; pointCoord(419, 2) =2.633266629202920e-05;
pointCoord(420, 0) = 3.910580956592711e-01; pointCoord(420, 1) =3.995598889345954e-01; pointCoord(420, 2) =2.633266629202920e-05;
pointCoord(421, 0) = 5.910580956592710e-01; pointCoord(421, 1) =3.995598889345954e-01; pointCoord(421, 2) =2.633266629202920e-05;
pointCoord(422, 0) = 1.910580956592710e-01; pointCoord(422, 1) =5.995598889345953e-01; pointCoord(422, 2) =2.633266629202920e-05;
pointCoord(423, 0) = 3.910580956592711e-01; pointCoord(423, 1) =5.995598889345953e-01; pointCoord(423, 2) =2.633266629202920e-05;
pointCoord(424, 0) = 1.910580956592710e-01; pointCoord(424, 1) =7.995598889345954e-01; pointCoord(424, 2) =2.633266629202920e-05;
pointCoord(425, 0) = 1.978349559785040e-01; pointCoord(425, 1) =1.560119750320643e-01; pointCoord(425, 2) =2.616879011700777e-04;
pointCoord(426, 0) = 3.978349559785040e-01; pointCoord(426, 1) =1.560119750320643e-01; pointCoord(426, 2) =2.616879011700777e-04;
pointCoord(427, 0) = 5.978349559785040e-01; pointCoord(427, 1) =1.560119750320643e-01; pointCoord(427, 2) =2.616879011700777e-04;
pointCoord(428, 0) = 7.978349559785041e-01; pointCoord(428, 1) =1.560119750320643e-01; pointCoord(428, 2) =2.616879011700777e-04;
pointCoord(429, 0) = 1.978349559785040e-01; pointCoord(429, 1) =3.560119750320643e-01; pointCoord(429, 2) =2.616879011700777e-04;
pointCoord(430, 0) = 3.978349559785040e-01; pointCoord(430, 1) =3.560119750320643e-01; pointCoord(430, 2) =2.616879011700777e-04;
pointCoord(431, 0) = 5.978349559785040e-01; pointCoord(431, 1) =3.560119750320643e-01; pointCoord(431, 2) =2.616879011700777e-04;
pointCoord(432, 0) = 1.978349559785040e-01; pointCoord(432, 1) =5.560119750320642e-01; pointCoord(432, 2) =2.616879011700777e-04;
pointCoord(433, 0) = 3.978349559785040e-01; pointCoord(433, 1) =5.560119750320642e-01; pointCoord(433, 2) =2.616879011700777e-04;
pointCoord(434, 0) = 1.978349559785040e-01; pointCoord(434, 1) =7.560119750320643e-01; pointCoord(434, 2) =2.616879011700777e-04;
pointCoord(435, 0) = 1.893494711142838e-01; pointCoord(435, 1) =1.644974598962845e-01; pointCoord(435, 2) =5.286497232810855e-04;
pointCoord(436, 0) = 3.893494711142838e-01; pointCoord(436, 1) =1.644974598962845e-01; pointCoord(436, 2) =5.286497232810855e-04;
pointCoord(437, 0) = 5.893494711142838e-01; pointCoord(437, 1) =1.644974598962845e-01; pointCoord(437, 2) =5.286497232810855e-04;
pointCoord(438, 0) = 7.893494711142839e-01; pointCoord(438, 1) =1.644974598962845e-01; pointCoord(438, 2) =5.286497232810855e-04;
pointCoord(439, 0) = 1.893494711142838e-01; pointCoord(439, 1) =3.644974598962846e-01; pointCoord(439, 2) =5.286497232810855e-04;
pointCoord(440, 0) = 3.893494711142838e-01; pointCoord(440, 1) =3.644974598962846e-01; pointCoord(440, 2) =5.286497232810855e-04;
pointCoord(441, 0) = 5.893494711142838e-01; pointCoord(441, 1) =3.644974598962846e-01; pointCoord(441, 2) =5.286497232810855e-04;
pointCoord(442, 0) = 1.893494711142838e-01; pointCoord(442, 1) =5.644974598962845e-01; pointCoord(442, 2) =5.286497232810855e-04;
pointCoord(443, 0) = 3.893494711142838e-01; pointCoord(443, 1) =5.644974598962845e-01; pointCoord(443, 2) =5.286497232810855e-04;
pointCoord(444, 0) = 1.893494711142838e-01; pointCoord(444, 1) =7.644974598962846e-01; pointCoord(444, 2) =5.286497232810855e-04;
pointCoord(445, 0) = 1.769234655052842e-01; pointCoord(445, 1) =1.769234655052842e-01; pointCoord(445, 2) =6.283429560853968e-04;
pointCoord(446, 0) = 3.769234655052842e-01; pointCoord(446, 1) =1.769234655052842e-01; pointCoord(446, 2) =6.283429560853968e-04;
pointCoord(447, 0) = 5.769234655052842e-01; pointCoord(447, 1) =1.769234655052842e-01; pointCoord(447, 2) =6.283429560853968e-04;
pointCoord(448, 0) = 7.769234655052842e-01; pointCoord(448, 1) =1.769234655052842e-01; pointCoord(448, 2) =6.283429560853968e-04;
pointCoord(449, 0) = 1.769234655052842e-01; pointCoord(449, 1) =3.769234655052842e-01; pointCoord(449, 2) =6.283429560853968e-04;
pointCoord(450, 0) = 3.769234655052842e-01; pointCoord(450, 1) =3.769234655052842e-01; pointCoord(450, 2) =6.283429560853968e-04;
pointCoord(451, 0) = 5.769234655052842e-01; pointCoord(451, 1) =3.769234655052842e-01; pointCoord(451, 2) =6.283429560853968e-04;
pointCoord(452, 0) = 1.769234655052842e-01; pointCoord(452, 1) =5.769234655052842e-01; pointCoord(452, 2) =6.283429560853968e-04;
pointCoord(453, 0) = 3.769234655052842e-01; pointCoord(453, 1) =5.769234655052842e-01; pointCoord(453, 2) =6.283429560853968e-04;
pointCoord(454, 0) = 1.769234655052842e-01; pointCoord(454, 1) =7.769234655052842e-01; pointCoord(454, 2) =6.283429560853968e-04;
pointCoord(455, 0) = 1.644974598962845e-01; pointCoord(455, 1) =1.893494711142838e-01; pointCoord(455, 2) =5.286497232810855e-04;
pointCoord(456, 0) = 3.644974598962846e-01; pointCoord(456, 1) =1.893494711142838e-01; pointCoord(456, 2) =5.286497232810855e-04;
pointCoord(457, 0) = 5.644974598962845e-01; pointCoord(457, 1) =1.893494711142838e-01; pointCoord(457, 2) =5.286497232810855e-04;
pointCoord(458, 0) = 7.644974598962846e-01; pointCoord(458, 1) =1.893494711142838e-01; pointCoord(458, 2) =5.286497232810855e-04;
pointCoord(459, 0) = 1.644974598962845e-01; pointCoord(459, 1) =3.893494711142838e-01; pointCoord(459, 2) =5.286497232810855e-04;
pointCoord(460, 0) = 3.644974598962846e-01; pointCoord(460, 1) =3.893494711142838e-01; pointCoord(460, 2) =5.286497232810855e-04;
pointCoord(461, 0) = 5.644974598962845e-01; pointCoord(461, 1) =3.893494711142838e-01; pointCoord(461, 2) =5.286497232810855e-04;
pointCoord(462, 0) = 1.644974598962845e-01; pointCoord(462, 1) =5.893494711142838e-01; pointCoord(462, 2) =5.286497232810855e-04;
pointCoord(463, 0) = 3.644974598962846e-01; pointCoord(463, 1) =5.893494711142838e-01; pointCoord(463, 2) =5.286497232810855e-04;
pointCoord(464, 0) = 1.644974598962845e-01; pointCoord(464, 1) =7.893494711142839e-01; pointCoord(464, 2) =5.286497232810855e-04;
pointCoord(465, 0) = 1.560119750320643e-01; pointCoord(465, 1) =1.978349559785040e-01; pointCoord(465, 2) =2.616879011700777e-04;
pointCoord(466, 0) = 3.560119750320643e-01; pointCoord(466, 1) =1.978349559785040e-01; pointCoord(466, 2) =2.616879011700777e-04;
pointCoord(467, 0) = 5.560119750320642e-01; pointCoord(467, 1) =1.978349559785040e-01; pointCoord(467, 2) =2.616879011700777e-04;
pointCoord(468, 0) = 7.560119750320643e-01; pointCoord(468, 1) =1.978349559785040e-01; pointCoord(468, 2) =2.616879011700777e-04;
pointCoord(469, 0) = 1.560119750320643e-01; pointCoord(469, 1) =3.978349559785040e-01; pointCoord(469, 2) =2.616879011700777e-04;
pointCoord(470, 0) = 3.560119750320643e-01; pointCoord(470, 1) =3.978349559785040e-01; pointCoord(470, 2) =2.616879011700777e-04;
pointCoord(471, 0) = 5.560119750320642e-01; pointCoord(471, 1) =3.978349559785040e-01; pointCoord(471, 2) =2.616879011700777e-04;
pointCoord(472, 0) = 1.560119750320643e-01; pointCoord(472, 1) =5.978349559785040e-01; pointCoord(472, 2) =2.616879011700777e-04;
pointCoord(473, 0) = 3.560119750320643e-01; pointCoord(473, 1) =5.978349559785040e-01; pointCoord(473, 2) =2.616879011700777e-04;
pointCoord(474, 0) = 1.560119750320643e-01; pointCoord(474, 1) =7.978349559785041e-01; pointCoord(474, 2) =2.616879011700777e-04;
pointCoord(475, 0) = 1.953089922969332e-01; pointCoord(475, 1) =1.046910077030668e-01; pointCoord(475, 2) =6.739253619376046e-04;
pointCoord(476, 0) = 3.953089922969332e-01; pointCoord(476, 1) =1.046910077030668e-01; pointCoord(476, 2) =6.739253619376046e-04;
pointCoord(477, 0) = 5.953089922969331e-01; pointCoord(477, 1) =1.046910077030668e-01; pointCoord(477, 2) =6.739253619376046e-04;
pointCoord(478, 0) = 7.953089922969332e-01; pointCoord(478, 1) =1.046910077030668e-01; pointCoord(478, 2) =6.739253619376046e-04;
pointCoord(479, 0) = 1.953089922969332e-01; pointCoord(479, 1) =3.046910077030668e-01; pointCoord(479, 2) =6.739253619376046e-04;
pointCoord(480, 0) = 3.953089922969332e-01; pointCoord(480, 1) =3.046910077030668e-01; pointCoord(480, 2) =6.739253619376046e-04;
pointCoord(481, 0) = 5.953089922969331e-01; pointCoord(481, 1) =3.046910077030668e-01; pointCoord(481, 2) =6.739253619376046e-04;
pointCoord(482, 0) = 1.953089922969332e-01; pointCoord(482, 1) =5.046910077030667e-01; pointCoord(482, 2) =6.739253619376046e-04;
pointCoord(483, 0) = 3.953089922969332e-01; pointCoord(483, 1) =5.046910077030667e-01; pointCoord(483, 2) =6.739253619376046e-04;
pointCoord(484, 0) = 1.953089922969332e-01; pointCoord(484, 1) =7.046910077030668e-01; pointCoord(484, 2) =6.739253619376046e-04;
pointCoord(485, 0) = 1.769234655052842e-01; pointCoord(485, 1) =1.230765344947159e-01; pointCoord(485, 2) =1.361432662753753e-03;
pointCoord(486, 0) = 3.769234655052842e-01; pointCoord(486, 1) =1.230765344947159e-01; pointCoord(486, 2) =1.361432662753753e-03;
pointCoord(487, 0) = 5.769234655052842e-01; pointCoord(487, 1) =1.230765344947159e-01; pointCoord(487, 2) =1.361432662753753e-03;
pointCoord(488, 0) = 7.769234655052842e-01; pointCoord(488, 1) =1.230765344947159e-01; pointCoord(488, 2) =1.361432662753753e-03;
pointCoord(489, 0) = 1.769234655052842e-01; pointCoord(489, 1) =3.230765344947159e-01; pointCoord(489, 2) =1.361432662753753e-03;
pointCoord(490, 0) = 3.769234655052842e-01; pointCoord(490, 1) =3.230765344947159e-01; pointCoord(490, 2) =1.361432662753753e-03;
pointCoord(491, 0) = 5.769234655052842e-01; pointCoord(491, 1) =3.230765344947159e-01; pointCoord(491, 2) =1.361432662753753e-03;
pointCoord(492, 0) = 1.769234655052842e-01; pointCoord(492, 1) =5.230765344947158e-01; pointCoord(492, 2) =1.361432662753753e-03;
pointCoord(493, 0) = 3.769234655052842e-01; pointCoord(493, 1) =5.230765344947158e-01; pointCoord(493, 2) =1.361432662753753e-03;
pointCoord(494, 0) = 1.769234655052842e-01; pointCoord(494, 1) =7.230765344947159e-01; pointCoord(494, 2) =1.361432662753753e-03;
pointCoord(495, 0) = 1.500000000000000e-01; pointCoord(495, 1) =1.500000000000000e-01; pointCoord(495, 2) =1.618172839506173e-03;
pointCoord(496, 0) = 3.500000000000000e-01; pointCoord(496, 1) =1.500000000000000e-01; pointCoord(496, 2) =1.618172839506173e-03;
pointCoord(497, 0) = 5.499999999999999e-01; pointCoord(497, 1) =1.500000000000000e-01; pointCoord(497, 2) =1.618172839506173e-03;
pointCoord(498, 0) = 7.500000000000000e-01; pointCoord(498, 1) =1.500000000000000e-01; pointCoord(498, 2) =1.618172839506173e-03;
pointCoord(499, 0) = 1.500000000000000e-01; pointCoord(499, 1) =3.500000000000000e-01; pointCoord(499, 2) =1.618172839506173e-03;
pointCoord(500, 0) = 3.500000000000000e-01; pointCoord(500, 1) =3.500000000000000e-01; pointCoord(500, 2) =1.618172839506173e-03;
pointCoord(501, 0) = 5.499999999999999e-01; pointCoord(501, 1) =3.500000000000000e-01; pointCoord(501, 2) =1.618172839506173e-03;
pointCoord(502, 0) = 1.500000000000000e-01; pointCoord(502, 1) =5.499999999999999e-01; pointCoord(502, 2) =1.618172839506173e-03;
pointCoord(503, 0) = 3.500000000000000e-01; pointCoord(503, 1) =5.499999999999999e-01; pointCoord(503, 2) =1.618172839506173e-03;
pointCoord(504, 0) = 1.500000000000000e-01; pointCoord(504, 1) =7.500000000000000e-01; pointCoord(504, 2) =1.618172839506173e-03;
pointCoord(505, 0) = 1.230765344947159e-01; pointCoord(505, 1) =1.769234655052842e-01; pointCoord(505, 2) =1.361432662753753e-03;
pointCoord(506, 0) = 3.230765344947159e-01; pointCoord(506, 1) =1.769234655052842e-01; pointCoord(506, 2) =1.361432662753753e-03;
pointCoord(507, 0) = 5.230765344947158e-01; pointCoord(507, 1) =1.769234655052842e-01; pointCoord(507, 2) =1.361432662753753e-03;
pointCoord(508, 0) = 7.230765344947159e-01; pointCoord(508, 1) =1.769234655052842e-01; pointCoord(508, 2) =1.361432662753753e-03;
pointCoord(509, 0) = 1.230765344947159e-01; pointCoord(509, 1) =3.769234655052842e-01; pointCoord(509, 2) =1.361432662753753e-03;
pointCoord(510, 0) = 3.230765344947159e-01; pointCoord(510, 1) =3.769234655052842e-01; pointCoord(510, 2) =1.361432662753753e-03;
pointCoord(511, 0) = 5.230765344947158e-01; pointCoord(511, 1) =3.769234655052842e-01; pointCoord(511, 2) =1.361432662753753e-03;
pointCoord(512, 0) = 1.230765344947159e-01; pointCoord(512, 1) =5.769234655052842e-01; pointCoord(512, 2) =1.361432662753753e-03;
pointCoord(513, 0) = 3.230765344947159e-01; pointCoord(513, 1) =5.769234655052842e-01; pointCoord(513, 2) =1.361432662753753e-03;
pointCoord(514, 0) = 1.230765344947159e-01; pointCoord(514, 1) =7.769234655052842e-01; pointCoord(514, 2) =1.361432662753753e-03;
pointCoord(515, 0) = 1.046910077030668e-01; pointCoord(515, 1) =1.953089922969332e-01; pointCoord(515, 2) =6.739253619376046e-04;
pointCoord(516, 0) = 3.046910077030668e-01; pointCoord(516, 1) =1.953089922969332e-01; pointCoord(516, 2) =6.739253619376046e-04;
pointCoord(517, 0) = 5.046910077030667e-01; pointCoord(517, 1) =1.953089922969332e-01; pointCoord(517, 2) =6.739253619376046e-04;
pointCoord(518, 0) = 7.046910077030668e-01; pointCoord(518, 1) =1.953089922969332e-01; pointCoord(518, 2) =6.739253619376046e-04;
pointCoord(519, 0) = 1.046910077030668e-01; pointCoord(519, 1) =3.953089922969332e-01; pointCoord(519, 2) =6.739253619376046e-04;
pointCoord(520, 0) = 3.046910077030668e-01; pointCoord(520, 1) =3.953089922969332e-01; pointCoord(520, 2) =6.739253619376046e-04;
pointCoord(521, 0) = 5.046910077030667e-01; pointCoord(521, 1) =3.953089922969332e-01; pointCoord(521, 2) =6.739253619376046e-04;
pointCoord(522, 0) = 1.046910077030668e-01; pointCoord(522, 1) =5.953089922969331e-01; pointCoord(522, 2) =6.739253619376046e-04;
pointCoord(523, 0) = 3.046910077030668e-01; pointCoord(523, 1) =5.953089922969331e-01; pointCoord(523, 2) =6.739253619376046e-04;
pointCoord(524, 0) = 1.046910077030668e-01; pointCoord(524, 1) =7.953089922969332e-01; pointCoord(524, 2) =6.739253619376046e-04;
pointCoord(525, 0) = 1.927830286153624e-01; pointCoord(525, 1) =5.337004037406934e-02; pointCoord(525, 2) =8.723120988299224e-04;
pointCoord(526, 0) = 3.927830286153624e-01; pointCoord(526, 1) =5.337004037406934e-02; pointCoord(526, 2) =8.723120988299224e-04;
pointCoord(527, 0) = 5.927830286153624e-01; pointCoord(527, 1) =5.337004037406934e-02; pointCoord(527, 2) =8.723120988299224e-04;
pointCoord(528, 0) = 7.927830286153624e-01; pointCoord(528, 1) =5.337004037406934e-02; pointCoord(528, 2) =8.723120988299224e-04;
pointCoord(529, 0) = 1.927830286153624e-01; pointCoord(529, 1) =2.533700403740693e-01; pointCoord(529, 2) =8.723120988299224e-04;
pointCoord(530, 0) = 3.927830286153624e-01; pointCoord(530, 1) =2.533700403740693e-01; pointCoord(530, 2) =8.723120988299224e-04;
pointCoord(531, 0) = 5.927830286153624e-01; pointCoord(531, 1) =2.533700403740693e-01; pointCoord(531, 2) =8.723120988299224e-04;
pointCoord(532, 0) = 1.927830286153624e-01; pointCoord(532, 1) =4.533700403740693e-01; pointCoord(532, 2) =8.723120988299224e-04;
pointCoord(533, 0) = 3.927830286153624e-01; pointCoord(533, 1) =4.533700403740693e-01; pointCoord(533, 2) =8.723120988299224e-04;
pointCoord(534, 0) = 1.927830286153624e-01; pointCoord(534, 1) =6.533700403740694e-01; pointCoord(534, 2) =8.723120988299224e-04;
pointCoord(535, 0) = 1.644974598962845e-01; pointCoord(535, 1) =8.165560909314720e-02; pointCoord(535, 2) =1.762204318958826e-03;
pointCoord(536, 0) = 3.644974598962846e-01; pointCoord(536, 1) =8.165560909314720e-02; pointCoord(536, 2) =1.762204318958826e-03;
pointCoord(537, 0) = 5.644974598962845e-01; pointCoord(537, 1) =8.165560909314720e-02; pointCoord(537, 2) =1.762204318958826e-03;
pointCoord(538, 0) = 7.644974598962846e-01; pointCoord(538, 1) =8.165560909314720e-02; pointCoord(538, 2) =1.762204318958826e-03;
pointCoord(539, 0) = 1.644974598962845e-01; pointCoord(539, 1) =2.816556090931472e-01; pointCoord(539, 2) =1.762204318958826e-03;
pointCoord(540, 0) = 3.644974598962846e-01; pointCoord(540, 1) =2.816556090931472e-01; pointCoord(540, 2) =1.762204318958826e-03;
pointCoord(541, 0) = 5.644974598962845e-01; pointCoord(541, 1) =2.816556090931472e-01; pointCoord(541, 2) =1.762204318958826e-03;
pointCoord(542, 0) = 1.644974598962845e-01; pointCoord(542, 1) =4.816556090931472e-01; pointCoord(542, 2) =1.762204318958826e-03;
pointCoord(543, 0) = 3.644974598962846e-01; pointCoord(543, 1) =4.816556090931472e-01; pointCoord(543, 2) =1.762204318958826e-03;
pointCoord(544, 0) = 1.644974598962845e-01; pointCoord(544, 1) =6.816556090931473e-01; pointCoord(544, 2) =1.762204318958826e-03;
pointCoord(545, 0) = 1.230765344947159e-01; pointCoord(545, 1) =1.230765344947159e-01; pointCoord(545, 2) =2.094522369422110e-03;
pointCoord(546, 0) = 3.230765344947159e-01; pointCoord(546, 1) =1.230765344947159e-01; pointCoord(546, 2) =2.094522369422110e-03;
pointCoord(547, 0) = 5.230765344947158e-01; pointCoord(547, 1) =1.230765344947159e-01; pointCoord(547, 2) =2.094522369422110e-03;
pointCoord(548, 0) = 7.230765344947159e-01; pointCoord(548, 1) =1.230765344947159e-01; pointCoord(548, 2) =2.094522369422110e-03;
pointCoord(549, 0) = 1.230765344947159e-01; pointCoord(549, 1) =3.230765344947159e-01; pointCoord(549, 2) =2.094522369422110e-03;
pointCoord(550, 0) = 3.230765344947159e-01; pointCoord(550, 1) =3.230765344947159e-01; pointCoord(550, 2) =2.094522369422110e-03;
pointCoord(551, 0) = 5.230765344947158e-01; pointCoord(551, 1) =3.230765344947159e-01; pointCoord(551, 2) =2.094522369422110e-03;
pointCoord(552, 0) = 1.230765344947159e-01; pointCoord(552, 1) =5.230765344947158e-01; pointCoord(552, 2) =2.094522369422110e-03;
pointCoord(553, 0) = 3.230765344947159e-01; pointCoord(553, 1) =5.230765344947158e-01; pointCoord(553, 2) =2.094522369422110e-03;
pointCoord(554, 0) = 1.230765344947159e-01; pointCoord(554, 1) =7.230765344947159e-01; pointCoord(554, 2) =2.094522369422110e-03;
pointCoord(555, 0) = 8.165560909314720e-02; pointCoord(555, 1) =1.644974598962845e-01; pointCoord(555, 2) =1.762204318958826e-03;
pointCoord(556, 0) = 2.816556090931472e-01; pointCoord(556, 1) =1.644974598962845e-01; pointCoord(556, 2) =1.762204318958826e-03;
pointCoord(557, 0) = 4.816556090931472e-01; pointCoord(557, 1) =1.644974598962845e-01; pointCoord(557, 2) =1.762204318958826e-03;
pointCoord(558, 0) = 6.816556090931473e-01; pointCoord(558, 1) =1.644974598962845e-01; pointCoord(558, 2) =1.762204318958826e-03;
pointCoord(559, 0) = 8.165560909314720e-02; pointCoord(559, 1) =3.644974598962846e-01; pointCoord(559, 2) =1.762204318958826e-03;
pointCoord(560, 0) = 2.816556090931472e-01; pointCoord(560, 1) =3.644974598962846e-01; pointCoord(560, 2) =1.762204318958826e-03;
pointCoord(561, 0) = 4.816556090931472e-01; pointCoord(561, 1) =3.644974598962846e-01; pointCoord(561, 2) =1.762204318958826e-03;
pointCoord(562, 0) = 8.165560909314720e-02; pointCoord(562, 1) =5.644974598962845e-01; pointCoord(562, 2) =1.762204318958826e-03;
pointCoord(563, 0) = 2.816556090931472e-01; pointCoord(563, 1) =5.644974598962845e-01; pointCoord(563, 2) =1.762204318958826e-03;
pointCoord(564, 0) = 8.165560909314720e-02; pointCoord(564, 1) =7.644974598962846e-01; pointCoord(564, 2) =1.762204318958826e-03;
pointCoord(565, 0) = 5.337004037406934e-02; pointCoord(565, 1) =1.927830286153624e-01; pointCoord(565, 2) =8.723120988299224e-04;
pointCoord(566, 0) = 2.533700403740693e-01; pointCoord(566, 1) =1.927830286153624e-01; pointCoord(566, 2) =8.723120988299224e-04;
pointCoord(567, 0) = 4.533700403740693e-01; pointCoord(567, 1) =1.927830286153624e-01; pointCoord(567, 2) =8.723120988299224e-04;
pointCoord(568, 0) = 6.533700403740694e-01; pointCoord(568, 1) =1.927830286153624e-01; pointCoord(568, 2) =8.723120988299224e-04;
pointCoord(569, 0) = 5.337004037406934e-02; pointCoord(569, 1) =3.927830286153624e-01; pointCoord(569, 2) =8.723120988299224e-04;
pointCoord(570, 0) = 2.533700403740693e-01; pointCoord(570, 1) =3.927830286153624e-01; pointCoord(570, 2) =8.723120988299224e-04;
pointCoord(571, 0) = 4.533700403740693e-01; pointCoord(571, 1) =3.927830286153624e-01; pointCoord(571, 2) =8.723120988299224e-04;
pointCoord(572, 0) = 5.337004037406934e-02; pointCoord(572, 1) =5.927830286153624e-01; pointCoord(572, 2) =8.723120988299224e-04;
pointCoord(573, 0) = 2.533700403740693e-01; pointCoord(573, 1) =5.927830286153624e-01; pointCoord(573, 2) =8.723120988299224e-04;
pointCoord(574, 0) = 5.337004037406934e-02; pointCoord(574, 1) =7.927830286153624e-01; pointCoord(574, 2) =8.723120988299224e-04;
pointCoord(575, 0) = 1.910580956592710e-01; pointCoord(575, 1) =1.832391974686259e-02; pointCoord(575, 2) =5.350108223322572e-04;
pointCoord(576, 0) = 3.910580956592711e-01; pointCoord(576, 1) =1.832391974686259e-02; pointCoord(576, 2) =5.350108223322572e-04;
pointCoord(577, 0) = 5.910580956592710e-01; pointCoord(577, 1) =1.832391974686259e-02; pointCoord(577, 2) =5.350108223322572e-04;
pointCoord(578, 0) = 7.910580956592711e-01; pointCoord(578, 1) =1.832391974686259e-02; pointCoord(578, 2) =5.350108223322572e-04;
pointCoord(579, 0) = 1.910580956592710e-01; pointCoord(579, 1) =2.183239197468626e-01; pointCoord(579, 2) =5.350108223322572e-04;
pointCoord(580, 0) = 3.910580956592711e-01; pointCoord(580, 1) =2.183239197468626e-01; pointCoord(580, 2) =5.350108223322572e-04;
pointCoord(581, 0) = 5.910580956592710e-01; pointCoord(581, 1) =2.183239197468626e-01; pointCoord(581, 2) =5.350108223322572e-04;
pointCoord(582, 0) = 1.910580956592710e-01; pointCoord(582, 1) =4.183239197468626e-01; pointCoord(582, 2) =5.350108223322572e-04;
pointCoord(583, 0) = 3.910580956592711e-01; pointCoord(583, 1) =4.183239197468626e-01; pointCoord(583, 2) =5.350108223322572e-04;
pointCoord(584, 0) = 1.910580956592710e-01; pointCoord(584, 1) =6.183239197468626e-01; pointCoord(584, 2) =5.350108223322572e-04;
pointCoord(585, 0) = 1.560119750320643e-01; pointCoord(585, 1) =5.337004037406934e-02; pointCoord(585, 2) =1.080803972647223e-03;
pointCoord(586, 0) = 3.560119750320643e-01; pointCoord(586, 1) =5.337004037406934e-02; pointCoord(586, 2) =1.080803972647223e-03;
pointCoord(587, 0) = 5.560119750320642e-01; pointCoord(587, 1) =5.337004037406934e-02; pointCoord(587, 2) =1.080803972647223e-03;
pointCoord(588, 0) = 7.560119750320643e-01; pointCoord(588, 1) =5.337004037406934e-02; pointCoord(588, 2) =1.080803972647223e-03;
pointCoord(589, 0) = 1.560119750320643e-01; pointCoord(589, 1) =2.533700403740693e-01; pointCoord(589, 2) =1.080803972647223e-03;
pointCoord(590, 0) = 3.560119750320643e-01; pointCoord(590, 1) =2.533700403740693e-01; pointCoord(590, 2) =1.080803972647223e-03;
pointCoord(591, 0) = 5.560119750320642e-01; pointCoord(591, 1) =2.533700403740693e-01; pointCoord(591, 2) =1.080803972647223e-03;
pointCoord(592, 0) = 1.560119750320643e-01; pointCoord(592, 1) =4.533700403740693e-01; pointCoord(592, 2) =1.080803972647223e-03;
pointCoord(593, 0) = 3.560119750320643e-01; pointCoord(593, 1) =4.533700403740693e-01; pointCoord(593, 2) =1.080803972647223e-03;
pointCoord(594, 0) = 1.560119750320643e-01; pointCoord(594, 1) =6.533700403740694e-01; pointCoord(594, 2) =1.080803972647223e-03;
pointCoord(595, 0) = 1.046910077030668e-01; pointCoord(595, 1) =1.046910077030668e-01; pointCoord(595, 2) =1.284622942592381e-03;
pointCoord(596, 0) = 3.046910077030668e-01; pointCoord(596, 1) =1.046910077030668e-01; pointCoord(596, 2) =1.284622942592381e-03;
pointCoord(597, 0) = 5.046910077030667e-01; pointCoord(597, 1) =1.046910077030668e-01; pointCoord(597, 2) =1.284622942592381e-03;
pointCoord(598, 0) = 7.046910077030668e-01; pointCoord(598, 1) =1.046910077030668e-01; pointCoord(598, 2) =1.284622942592381e-03;
pointCoord(599, 0) = 1.046910077030668e-01; pointCoord(599, 1) =3.046910077030668e-01; pointCoord(599, 2) =1.284622942592381e-03;
pointCoord(600, 0) = 3.046910077030668e-01; pointCoord(600, 1) =3.046910077030668e-01; pointCoord(600, 2) =1.284622942592381e-03;
pointCoord(601, 0) = 5.046910077030667e-01; pointCoord(601, 1) =3.046910077030668e-01; pointCoord(601, 2) =1.284622942592381e-03;
pointCoord(602, 0) = 1.046910077030668e-01; pointCoord(602, 1) =5.046910077030667e-01; pointCoord(602, 2) =1.284622942592381e-03;
pointCoord(603, 0) = 3.046910077030668e-01; pointCoord(603, 1) =5.046910077030667e-01; pointCoord(603, 2) =1.284622942592381e-03;
pointCoord(604, 0) = 1.046910077030668e-01; pointCoord(604, 1) =7.046910077030668e-01; pointCoord(604, 2) =1.284622942592381e-03;
pointCoord(605, 0) = 5.337004037406934e-02; pointCoord(605, 1) =1.560119750320643e-01; pointCoord(605, 2) =1.080803972647223e-03;
pointCoord(606, 0) = 2.533700403740693e-01; pointCoord(606, 1) =1.560119750320643e-01; pointCoord(606, 2) =1.080803972647223e-03;
pointCoord(607, 0) = 4.533700403740693e-01; pointCoord(607, 1) =1.560119750320643e-01; pointCoord(607, 2) =1.080803972647223e-03;
pointCoord(608, 0) = 6.533700403740694e-01; pointCoord(608, 1) =1.560119750320643e-01; pointCoord(608, 2) =1.080803972647223e-03;
pointCoord(609, 0) = 5.337004037406934e-02; pointCoord(609, 1) =3.560119750320643e-01; pointCoord(609, 2) =1.080803972647223e-03;
pointCoord(610, 0) = 2.533700403740693e-01; pointCoord(610, 1) =3.560119750320643e-01; pointCoord(610, 2) =1.080803972647223e-03;
pointCoord(611, 0) = 4.533700403740693e-01; pointCoord(611, 1) =3.560119750320643e-01; pointCoord(611, 2) =1.080803972647223e-03;
pointCoord(612, 0) = 5.337004037406934e-02; pointCoord(612, 1) =5.560119750320642e-01; pointCoord(612, 2) =1.080803972647223e-03;
pointCoord(613, 0) = 2.533700403740693e-01; pointCoord(613, 1) =5.560119750320642e-01; pointCoord(613, 2) =1.080803972647223e-03;
pointCoord(614, 0) = 5.337004037406934e-02; pointCoord(614, 1) =7.560119750320643e-01; pointCoord(614, 2) =1.080803972647223e-03;
pointCoord(615, 0) = 1.832391974686259e-02; pointCoord(615, 1) =1.910580956592710e-01; pointCoord(615, 2) =5.350108223322572e-04;
pointCoord(616, 0) = 2.183239197468626e-01; pointCoord(616, 1) =1.910580956592710e-01; pointCoord(616, 2) =5.350108223322572e-04;
pointCoord(617, 0) = 4.183239197468626e-01; pointCoord(617, 1) =1.910580956592710e-01; pointCoord(617, 2) =5.350108223322572e-04;
pointCoord(618, 0) = 6.183239197468626e-01; pointCoord(618, 1) =1.910580956592710e-01; pointCoord(618, 2) =5.350108223322572e-04;
pointCoord(619, 0) = 1.832391974686259e-02; pointCoord(619, 1) =3.910580956592711e-01; pointCoord(619, 2) =5.350108223322572e-04;
pointCoord(620, 0) = 2.183239197468626e-01; pointCoord(620, 1) =3.910580956592711e-01; pointCoord(620, 2) =5.350108223322572e-04;
pointCoord(621, 0) = 4.183239197468626e-01; pointCoord(621, 1) =3.910580956592711e-01; pointCoord(621, 2) =5.350108223322572e-04;
pointCoord(622, 0) = 1.832391974686259e-02; pointCoord(622, 1) =5.910580956592710e-01; pointCoord(622, 2) =5.350108223322572e-04;
pointCoord(623, 0) = 2.183239197468626e-01; pointCoord(623, 1) =5.910580956592710e-01; pointCoord(623, 2) =5.350108223322572e-04;
pointCoord(624, 0) = 1.832391974686259e-02; pointCoord(624, 1) =7.910580956592711e-01; pointCoord(624, 2) =5.350108223322572e-04;
    }
    else
    {
        std::cout << "SELECT A VALID NUMBER FOR THE AMOUNT OF HAMMER POINTS. THE SELECTED NUMBER IS" << nH << "." << std::endl;
    }

    return pointCoord;
    //

    // return pointCoord;
}

bounded_matrix<double, 2, 2> Element::referenceJacobianMatrix(const double &xsi1, const double &xsi2)
{
    matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);
    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    for (int i = 0; i < connection_.size(); i++)
    {
        bounded_vector<double, 2> initialCoord = connection_[i]->getInitialCoordinate();
        dx1_dxsi1 += initialCoord(0) * dphi_dxsi(0, i);
        dx1_dxsi2 += initialCoord(0) * dphi_dxsi(1, i);
        dx2_dxsi1 += initialCoord(1) * dphi_dxsi(0, i);
        dx2_dxsi2 += initialCoord(1) * dphi_dxsi(1, i);
    }

    bounded_matrix<double, 2, 2> referenceJacobianMatrix;
    referenceJacobianMatrix(0, 0) = dx1_dxsi1;
    referenceJacobianMatrix(1, 0) = dx2_dxsi1;
    referenceJacobianMatrix(0, 1) = dx1_dxsi2;
    referenceJacobianMatrix(1, 1) = dx2_dxsi2;

    return referenceJacobianMatrix;
}

bounded_matrix<double, 2, 2> Element::currentJacobianMatrix(const double &xsi1, const double &xsi2)
{
    matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);
    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    for (int i = 0; i < connection_.size(); i++)
    {
        bounded_vector<double, 2> currentCoord = connection_[i]->getCurrentCoordinate();
        dx1_dxsi1 += currentCoord(0) * dphi_dxsi(0, i);
        dx1_dxsi2 += currentCoord(0) * dphi_dxsi(1, i);
        dx2_dxsi1 += currentCoord(1) * dphi_dxsi(0, i);
        dx2_dxsi2 += currentCoord(1) * dphi_dxsi(1, i);
    }

    bounded_matrix<double, 2, 2> currentJacobianMatrix;
    currentJacobianMatrix(0, 0) = dx1_dxsi1;
    currentJacobianMatrix(1, 0) = dx2_dxsi1;
    currentJacobianMatrix(0, 1) = dx1_dxsi2;
    currentJacobianMatrix(1, 1) = dx2_dxsi2;

    return currentJacobianMatrix;
}

double Element::jacobianDeterminant(const bounded_matrix<double, 2, 2> &jacobianMatrix)
{
    return (jacobianMatrix(0, 0) * jacobianMatrix(1, 1) - jacobianMatrix(0, 1) * jacobianMatrix(1, 0));
}

std::pair<vector<double>, matrix<double>> Element::elementContributions(const std::string &ep, const std::string &typeAnalyze,
                                                                        const int &step, const int &numberOfStep, const double &deltat,
                                                                        const double &beta, const double &gama,
                                                                        const int &hammerPoints, const int &hammerPointsBlendZone)
{
    std::vector<Cell *> cellAux;
    std::vector<int> freedom;
    int freedomDegree = connection_.size();
    for (Cell *cell : incidenceCell_)
    {
        bool teste = false;
        for (int i = 0; i < cellAux.size(); i++)
        {
            if (cell == cellAux[i])
            {
                teste = true;
            }
        }
        if (teste == false)
        {
            cellAux.push_back(cell);
            freedom.push_back(freedomDegree);
            freedomDegree += cell->getControlPoints().size();
        }
    }

    vector<double> rhs(2 * freedomDegree, 0.0);
    matrix<double> tangent(2 * freedomDegree, 2 * freedomDegree, 0.0);

    matrix<double> domainIntegrationPoints;
    if (insideBlendZone_.size() > 0) //the element has at least one point in blend zone;
    {
        domainIntegrationPoints = hammerQuadrature(hammerPointsBlendZone);
    }
    else
    {
        domainIntegrationPoints = hammerQuadrature(hammerPoints);
    }

    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);
    if (typeAnalyze == "STATIC")
    {
        density = 0.0;
    }
    double thickness = mesh_->getThickness();

    int auxiliar = 0;

    if (insideBlendZone_.size() > 0) //the element has at least one hammer point in blending zone
    {
        for (int ih = 0; ih < hammerPointsBlendZone; ih++)
        {
            if (insideBlendZone_[ih] == false)
            {
                double xsi1 = domainIntegrationPoints(ih, 0);
                double xsi2 = domainIntegrationPoints(ih, 1);
                double weight = domainIntegrationPoints(ih, 2);

                vector<double> phi = domainShapeFunction(xsi1, xsi2);
                matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);  //row = direction, column = node
                bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
                double j0 = jacobianDeterminant(A0);
                bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
                A0I(0, 0) = A0(1, 1) / j0;
                A0I(1, 1) = A0(0, 0) / j0;
                A0I(0, 1) = -A0(0, 1) / j0;
                A0I(1, 0) = -A0(1, 0) / j0;
                bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2); //current configuration map
                bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                     //current deformation gradient
                identity_matrix<double> I(2);                                        //identity matrix
                bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);   //current green strain tensor

                bounded_matrix<double, 2, 2> S; //second piola kirchhoff stress tensor

                // bounded_vector<double, 2> coordPo;
                // coordPo(0) = 0.0;
                // coordPo(1) = 0.0;
                // int ik = 0;
                // for (Node *n : connection_)
                // {
                //     coordPo += n->getInitialCoordinate() * phi(ik++);
                // }

                if (ep == "EPD")
                {
                    S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
                    S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
                    S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                    S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
                }
                else
                {
                    S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
                    S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
                    S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                    S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
                }

                bounded_matrix<double, 2, 2> dA_dy;
                bounded_matrix<double, 2, 2> dA_dy2;

                for (int i = 0; i < connection_.size(); i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (j == 0)
                        {
                            dA_dy(0, 0) = dphi_dxsi(0, i);
                            dA_dy(0, 1) = dphi_dxsi(1, i);
                            dA_dy(1, 0) = 0.0;
                            dA_dy(1, 1) = 0.0;
                        }
                        else
                        {
                            dA_dy(1, 0) = dphi_dxsi(0, i);
                            dA_dy(1, 1) = dphi_dxsi(1, i);
                            dA_dy(0, 0) = 0.0;
                            dA_dy(0, 1) = 0.0;
                        }

                        bounded_matrix<double, 2, 2> mat1 = prod(trans(A0I), trans(dA_dy));
                        bounded_matrix<double, 2, 2> mat2 = prod(dA_dy, A0I);

                        bounded_matrix<double, 2, 2> dE_dy = 0.5 * (prod(mat1, Ac) + prod(trans(Ac), mat2)); //first derivative of E regarding i,j

                        double r = dE_dy(0, 0) * S(0, 0) + dE_dy(1, 1) * S(1, 1) + dE_dy(0, 1) * S(0, 1) + dE_dy(1, 0) * S(1, 0); //internal force

                        double accel = 0.0;

                        if (typeAnalyze == "DYNAMIC")
                        {
                            for (int m = 0; m < connection_.size(); m++)
                            {
                                accel += phi(m) * connection_[m]->getCurrentAcceleration()(j);
                            }
                        }

                        double m = density * phi(i) * accel; //inertial force

                        double shape = 0.0;
                        // if (j == 0)
                        // {
                        //     // shape = cos(3.1415926535897932384626433 * coordPo(0) / 10.0) * functions.first(i);
                        //     if (coordPo(0) <= 10.0 or coordPo(0) >= 20.0)
                        //     {
                        //         shape = 0.0;
                        //     }
                        //     else
                        //     {
                        //         shape = sin(3.1415926535897932384626433 / 10.0 * (coordPo(0) - 10.0)) * phi(i);
                        //     }
                        // }
                        rhs(2 * i + j) -= (r + m - shape) * weight * j0 * thickness;

                        for (int k = 0; k < connection_.size(); k++)
                        {
                            for (int l = 0; l < 2; l++)
                            {
                                if (l == 0)
                                {
                                    dA_dy2(0, 0) = dphi_dxsi(0, k);
                                    dA_dy2(0, 1) = dphi_dxsi(1, k);
                                    dA_dy2(1, 0) = 0.0;
                                    dA_dy2(1, 1) = 0.0;
                                }
                                else
                                {
                                    dA_dy2(1, 0) = dphi_dxsi(0, k);
                                    dA_dy2(1, 1) = dphi_dxsi(1, k);
                                    dA_dy2(0, 0) = 0.0;
                                    dA_dy2(0, 1) = 0.0;
                                }

                                bounded_matrix<double, 2, 2> mat3 = prod(trans(A0I), trans(dA_dy2));
                                bounded_matrix<double, 2, 2> mat4 = prod(dA_dy2, A0I);

                                bounded_matrix<double, 2, 2> dE_dy2 = 0.5 * (prod(mat3, Ac) + prod(trans(Ac), mat4)); //first derivative of E regarding k,l
                                bounded_matrix<double, 2, 2> dE_dy3 = 0.5 * (prod(mat1, mat4) + prod(mat3, mat2));    //second derivative of E regarding i,j,k,l

                                bounded_matrix<double, 2, 2> dS_dy;
                                if (ep == "EPD")
                                {
                                    dS_dy(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(0, 0) + poisson * dE_dy(1, 1)));
                                    dS_dy(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(1, 1) + poisson * dE_dy(0, 0)));
                                    dS_dy(1, 0) = (young / (1.0 + poisson)) * dE_dy(1, 0);
                                    dS_dy(0, 1) = (young / (1.0 + poisson)) * dE_dy(0, 1);
                                }
                                else
                                {
                                    dS_dy(0, 0) = (young / (1.0 - poisson * poisson)) * (dE_dy(0, 0) + poisson * dE_dy(1, 1));
                                    dS_dy(1, 1) = (young / (1.0 - poisson * poisson)) * (dE_dy(1, 1) + poisson * dE_dy(0, 0));
                                    dS_dy(1, 0) = (young / (1.0 + poisson)) * dE_dy(1, 0);
                                    dS_dy(0, 1) = (young / (1.0 + poisson)) * dE_dy(0, 1);
                                }

                                double v = dE_dy3(0, 0) * S(0, 0) + dE_dy3(1, 1) * S(1, 1) + dE_dy3(0, 1) * S(0, 1) + dE_dy3(1, 0) * S(1, 0) +                //second part of equation 5.88
                                           dE_dy2(0, 0) * dS_dy(0, 0) + dE_dy2(1, 1) * dS_dy(1, 1) + dE_dy2(0, 1) * dS_dy(0, 1) + dE_dy2(1, 0) * dS_dy(1, 0); //viscosity and pressure contribution

                                tangent(2 * i + j, 2 * k + l) += v * weight * j0 * thickness;
                                if (j == l)
                                {
                                    double mm = density * phi(i) * phi(k); //mass contribution
                                    tangent(2 * i + j, 2 * k + l) += (mm / (beta * deltat * deltat)) * weight * j0 * thickness;
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                //os próximos 4 intens foram calculados em GlobalSolid::incidenceLocalxGlobal() e enviados para cada elemento;
                //os vetores utilizados abaixo possuem tamanho igual a quantidade de pontos de hammer que estão dentro da blend zone, por isso o uso da variável "auxiliar"...
                Cell *cell = incidenceCell_[auxiliar];                                   //célula em que o ponto de hammer está inserido
                const bounded_vector<double, 2> xsiGlobal = xsiIncidenceCell_[auxiliar]; //coordenadas adimensionais da célula que correspondem aos pontos de hammer nas coordenadas globais;
                const double blend = bValue_[auxiliar];                                  //valor da blend function no ponto de hammer, em que a distância assinalada foi interpolada;
                const bounded_vector<double, 2> dblend_dxsi = db_dxsiValue_[auxiliar];   //derivada de b em relação a xsiLocal = db_dDAhammer * dDAhammer_dxsiLocal, sendo DAhammer a distância assinalada interpolada;

                //Pontos de Hammer para integração do elemento finito triangular
                const double xsi1 = domainIntegrationPoints(ih, 0);
                const double xsi2 = domainIntegrationPoints(ih, 1);
                const double weight = domainIntegrationPoints(ih, 2);

                vector<double> phiLocal = domainShapeFunction(xsi1, xsi2);                 //funções de forma calculada no ponto de hammer
                matrix<double> dphi_dxsiLocal = domainDerivativeShapeFunction(xsi1, xsi2); //derivadas das funções de forma em relação a xsi local (direção, nó);

                //GLOBAL
                std::vector<ControlPoint *> cps = cell->getControlPoints();
                int nglobal = cps.size(); //quantidade de pontos de controle que definem a célula em que o ponto de hammer está inserido
                vector<double> wpc(nglobal);
                bounded_vector<int, 2> inc;
                inc = cps[nglobal - 1]->getINC();
                for (int i = 0; i < nglobal; i++)
                {
                    wpc(i) = cps[i]->getWeight();
                }                                                                                                                       //conectividade da célula
                const std::pair<vector<double>, matrix<double>> functionsGlobal = cell->shapeFunctionAndDerivates(xsiGlobal, wpc, inc); //functionsGlobal.first(i) são as funções de forma globais calculadas em xsiGlobal que houve incidência;
                                                                                                                                        //functionsGlobal.second(j, i) são as derivadas das funções de forma globais em relação a xsi global;

                bounded_matrix<double, 2, 2> Jlocal = referenceJacobianMatrix(xsi1, xsi2);
                bounded_matrix<double, 2, 2> JGlobalI = inverseMatrix(cell->referenceJacobianMatrix(functionsGlobal.second));
                bounded_matrix<double, 2, 2> M = trans(prod(JGlobalI, Jlocal)); //para transforma dPhiGlobal_dXsiGlobal em dPhiGlobal_dXsiLocal

                int nlocal = connection_.size();               //poderia estar de fora...
                int nnum = nlocal + nglobal;                   //novo número funções de forma que definem o elemento na blend zone
                vector<double> phiBlended(nnum, 0.0);          //novas funções de forma que definem o elemento
                matrix<double> dphi_dxsiBlended(2, nnum, 0.0); //derivadas das novas funções de forma em relação a xsi local

                //novas funções de forma relacionadas as funções do elemento local
                for (int i = 0; i < nlocal; i++)
                {
                    phiBlended(i) = (1.0 - blend) * phiLocal(i);
                    dphi_dxsiBlended(0, i) = -dblend_dxsi(0) * phiLocal(i) + (1.0 - blend) * dphi_dxsiLocal(0, i);
                    dphi_dxsiBlended(1, i) = -dblend_dxsi(1) * phiLocal(i) + (1.0 - blend) * dphi_dxsiLocal(1, i);
                }

                //novas funções de forma relacionadas relacionadas as funções do global
                for (int i = 0; i < nglobal; i++)
                {
                    phiBlended(nlocal + i) = blend * functionsGlobal.first(i);

                    bounded_vector<double, 2> dPhiG_dXsiG, dPhiG_dxsiL;
                    dPhiG_dXsiG(0) = functionsGlobal.second(0, i); //dPhiGlobal_dxsiGlobal1
                    dPhiG_dXsiG(1) = functionsGlobal.second(1, i); //dPhiGlobal_dxsiGlobal2
                    dPhiG_dxsiL = prod(M, dPhiG_dXsiG);

                    dphi_dxsiBlended(0, nlocal + i) = dblend_dxsi(0) * functionsGlobal.first(i) + blend * dPhiG_dxsiL(0);
                    dphi_dxsiBlended(1, nlocal + i) = dblend_dxsi(1) * functionsGlobal.first(i) + blend * dPhiG_dxsiL(1);
                }

                int aux;
                for (int iaux = 0; iaux < cellAux.size(); iaux++)
                {
                    if (cell == cellAux[iaux])
                    {
                        aux = freedom[iaux]; //para auxiliar na contribuição da matriz e vetor quando o elemento triangular incide sobre mais de uma célula...
                        break;
                    }
                }

                bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrixBlendZone(dphi_dxsiBlended, cps); //mapeamento da configuração inicial
                bounded_matrix<double, 2, 2> A1 = currentJacobianMatrixBlendZone(dphi_dxsiBlended, cps);   //mapeamento da configuração atual
                bounded_matrix<double, 2, 2> A0I = inverseMatrix(A0);                                      //inverse initial configuration map
                bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                                           //gradiente da função mudança de configuração
                identity_matrix<double> I(2);                                                              //identity matrix
                bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);                         //current green strain tensor
                bounded_matrix<double, 2, 2> S;                                                            //second piola kirchhoff stress tensor
                double j0 = jacobianDeterminant(Jlocal);                                                   //para multiplicar com o peso

                // bounded_vector<double, 2> coordPo;
                // coordPo(0) = 0.0;
                // coordPo(1) = 0.0;
                // int ik = 0;
                // for (Node *n : connection_)
                // {
                //     coordPo += n->getInitialCoordinate() * phiBlended(ik++);
                // }
                // for (ControlPoint *n : cps)
                // {
                //     coordPo += n->getInitialCoordinate() * phiBlended(ik++);
                // }

                if (ep == "EPD")
                {
                    S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
                    S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
                    S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                    S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
                }
                else
                {
                    S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
                    S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
                    S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                    S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
                }

                bounded_matrix<double, 2, 2> dA_dy;  //deriavada de A1 em relação a y
                bounded_matrix<double, 2, 2> dA_dy2; //deriavada de A1 em relação a y

                for (int i = 0; i < nnum; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (j == 0)
                        {
                            dA_dy(0, 0) = dphi_dxsiBlended(0, i);
                            dA_dy(0, 1) = dphi_dxsiBlended(1, i);
                            dA_dy(1, 0) = 0.0;
                            dA_dy(1, 1) = 0.0;
                        }
                        else
                        {
                            dA_dy(1, 0) = dphi_dxsiBlended(0, i);
                            dA_dy(1, 1) = dphi_dxsiBlended(1, i);
                            dA_dy(0, 0) = 0.0;
                            dA_dy(0, 1) = 0.0;
                        }

                        bounded_matrix<double, 2, 2> mat1 = prod(trans(A0I), trans(dA_dy));
                        bounded_matrix<double, 2, 2> mat2 = prod(dA_dy, A0I);

                        bounded_matrix<double, 2, 2> dE_dy = 0.5 * (prod(mat1, Ac) + prod(trans(Ac), mat2)); //first derivative of E regarding i,j

                        double r = dE_dy(0, 0) * S(0, 0) + dE_dy(1, 1) * S(1, 1) + dE_dy(0, 1) * S(0, 1) + dE_dy(1, 0) * S(1, 0); //internal force

                        double accel = 0.0;

                        if (typeAnalyze == "DYNAMIC")
                        {
                            for (int m = 0; m < nlocal; m++)
                            {
                                accel += phiBlended(m) * connection_[m]->getCurrentAcceleration()(j);
                            }
                            for (int m = 0; m < nglobal; m++)
                            {
                                accel += phiBlended(m + nlocal) * cps[m]->getCurrentAcceleration()(j);
                            }
                        }

                        double m = density * phiBlended(i) * accel; //inertial force

                        double shape = 0.0;
                        // if (j == 0)
                        // {
                        //     if (coordPo(0) <= 30.0 / 3.0 or coordPo(0) >= 2.0 * 30.0 / 3.0)
                        //     {
                        //         shape = 0.0;
                        //     }
                        //     else
                        //     {
                        //         shape = sin(3.1415926535897932384626433 / 10.0 * (coordPo(0) - 10.0)) * phiBlended(i);
                        //     }
                        // }

                        if (i < nlocal)
                        {
                            rhs(2 * i + j) -= (r + m - shape) * weight * j0 * thickness;
                        }
                        else if (i >= nlocal)
                        {
                            rhs(2 * (i + aux - nlocal) + j) -= (r + m - shape) * weight * j0 * thickness;
                        }

                        for (int k = 0; k < nnum; k++)
                        {
                            for (int l = 0; l < 2; l++)
                            {
                                if (l == 0)
                                {
                                    dA_dy2(0, 0) = dphi_dxsiBlended(0, k);
                                    dA_dy2(0, 1) = dphi_dxsiBlended(1, k);
                                    dA_dy2(1, 0) = 0.0;
                                    dA_dy2(1, 1) = 0.0;
                                }
                                else
                                {
                                    dA_dy2(1, 0) = dphi_dxsiBlended(0, k);
                                    dA_dy2(1, 1) = dphi_dxsiBlended(1, k);
                                    dA_dy2(0, 0) = 0.0;
                                    dA_dy2(0, 1) = 0.0;
                                }

                                bounded_matrix<double, 2, 2> mat3 = prod(trans(A0I), trans(dA_dy2));
                                bounded_matrix<double, 2, 2> mat4 = prod(dA_dy2, A0I);

                                bounded_matrix<double, 2, 2> dE_dy2 = 0.5 * (prod(mat3, Ac) + prod(trans(Ac), mat4)); //first derivative of E regarding k,l
                                bounded_matrix<double, 2, 2> dE_dy3 = 0.5 * (prod(mat1, mat4) + prod(mat3, mat2));    //second derivative of E regarding i,j,k,l

                                bounded_matrix<double, 2, 2> dS_dy;
                                if (ep == "EPD")
                                {
                                    dS_dy(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(0, 0) + poisson * dE_dy(1, 1)));
                                    dS_dy(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(1, 1) + poisson * dE_dy(0, 0)));
                                    dS_dy(1, 0) = (young / (1.0 + poisson)) * dE_dy(1, 0);
                                    dS_dy(0, 1) = (young / (1.0 + poisson)) * dE_dy(0, 1);
                                }
                                else
                                {
                                    dS_dy(0, 0) = (young / (1.0 - poisson * poisson)) * (dE_dy(0, 0) + poisson * dE_dy(1, 1));
                                    dS_dy(1, 1) = (young / (1.0 - poisson * poisson)) * (dE_dy(1, 1) + poisson * dE_dy(0, 0));
                                    dS_dy(1, 0) = (young / (1.0 + poisson)) * dE_dy(1, 0);
                                    dS_dy(0, 1) = (young / (1.0 + poisson)) * dE_dy(0, 1);
                                }

                                double v = dE_dy3(0, 0) * S(0, 0) + dE_dy3(1, 1) * S(1, 1) + dE_dy3(0, 1) * S(0, 1) + dE_dy3(1, 0) * S(1, 0) +                //second part of equation 5.88
                                           dE_dy2(0, 0) * dS_dy(0, 0) + dE_dy2(1, 1) * dS_dy(1, 1) + dE_dy2(0, 1) * dS_dy(0, 1) + dE_dy2(1, 0) * dS_dy(1, 0); //viscosity and pressure contribution

                                if (i < nlocal and k < nlocal)
                                {
                                    tangent(2 * i + j, 2 * k + l) += v * weight * j0 * thickness;
                                }
                                else if (i >= nlocal and k < nlocal)
                                {
                                    tangent(2 * (aux + i - nlocal) + j, 2 * k + l) += v * weight * j0 * thickness;
                                }
                                else if (i < nlocal and k >= nlocal)
                                {
                                    tangent(2 * i + j, 2 * (aux + k - nlocal) + l) += v * weight * j0 * thickness;
                                }
                                else if (i >= nlocal and k >= nlocal)
                                {
                                    tangent(2 * (aux + i - nlocal) + j, 2 * (aux + k - nlocal) + l) += v * weight * j0 * thickness;
                                }

                                if (j == l and typeAnalyze == "DYNAMIC")
                                {
                                    double mm = density * phiBlended(i) * phiBlended(k); //mass contribution

                                    if (i < nlocal and k < nlocal)
                                    {
                                        tangent(2 * i + j, 2 * k + l) += (mm / (beta * deltat * deltat)) * weight * j0 * thickness;
                                    }
                                    else if (i >= nlocal and k < nlocal)
                                    {
                                        tangent(2 * (aux + i - nlocal) + j, 2 * k + l) += (mm / (beta * deltat * deltat)) * weight * j0 * thickness;
                                    }
                                    else if (i < nlocal and k >= nlocal)
                                    {
                                        tangent(2 * i + j, 2 * (aux + k - nlocal) + l) += (mm / (beta * deltat * deltat)) * weight * j0 * thickness;
                                    }
                                    else if (i >= nlocal and k >= nlocal)
                                    {
                                        tangent(2 * (aux + i - nlocal) + j, 2 * (aux + k - nlocal) + l) += (mm / (beta * deltat * deltat)) * weight * j0 * thickness;
                                    }
                                }
                            }
                        }
                    }
                }
                auxiliar += 1;
            }
        }
    }
    else //the element doesn't have no hammer point in blending zone
    {
        for (int ih = 0; ih < hammerPoints; ih++)
        {
            double xsi1 = domainIntegrationPoints(ih, 0);
            double xsi2 = domainIntegrationPoints(ih, 1);
            double weight = domainIntegrationPoints(ih, 2);

            vector<double> phi = domainShapeFunction(xsi1, xsi2);
            matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);  //row = direction, column = node
            bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
            double j0 = jacobianDeterminant(A0);
            bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
            A0I(0, 0) = A0(1, 1) / j0;
            A0I(1, 1) = A0(0, 0) / j0;
            A0I(0, 1) = -A0(0, 1) / j0;
            A0I(1, 0) = -A0(1, 0) / j0;
            bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2); //current configuration map
            bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                     //current deformation gradient
            identity_matrix<double> I(2);                                        //identity matrix
            bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);   //current green strain tensor

            bounded_matrix<double, 2, 2> S; //second piola kirchhoff stress tensor

            // bounded_vector<double, 2> coordPo;
            // coordPo(0) = 0.0;
            // coordPo(1) = 0.0;
            // int ik = 0;
            // for (Node *n : connection_)
            // {
            //     coordPo += n->getInitialCoordinate() * phi(ik++);
            // }

            if (ep == "EPD")
            {
                S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
                S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
                S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
            }
            else
            {
                S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
                S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
                S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
                S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
            }

            bounded_matrix<double, 2, 2> dA_dy;
            bounded_matrix<double, 2, 2> dA_dy2;

            for (int i = 0; i < connection_.size(); i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    if (j == 0)
                    {
                        dA_dy(0, 0) = dphi_dxsi(0, i);
                        dA_dy(0, 1) = dphi_dxsi(1, i);
                        dA_dy(1, 0) = 0.0;
                        dA_dy(1, 1) = 0.0;
                    }
                    else
                    {
                        dA_dy(1, 0) = dphi_dxsi(0, i);
                        dA_dy(1, 1) = dphi_dxsi(1, i);
                        dA_dy(0, 0) = 0.0;
                        dA_dy(0, 1) = 0.0;
                    }

                    bounded_matrix<double, 2, 2> mat1 = prod(trans(A0I), trans(dA_dy));
                    bounded_matrix<double, 2, 2> mat2 = prod(dA_dy, A0I);

                    bounded_matrix<double, 2, 2> dE_dy = 0.5 * (prod(mat1, Ac) + prod(trans(Ac), mat2)); //first derivative of E regarding i,j

                    double r = dE_dy(0, 0) * S(0, 0) + dE_dy(1, 1) * S(1, 1) + dE_dy(0, 1) * S(0, 1) + dE_dy(1, 0) * S(1, 0); //internal force

                    double accel = 0.0;

                    if (typeAnalyze == "DYNAMIC")
                    {
                        for (int m = 0; m < connection_.size(); m++)
                        {
                            accel += phi(m) * connection_[m]->getCurrentAcceleration()(j);
                        }
                    }

                    double m = density * phi(i) * accel; //inertial force

                    double shape = 0.0;
                    // if (j == 0)
                    // {
                    //     if (coordPo(0) <= 30.0 / 3.0 or coordPo(0) >= 2.0 * 30.0 / 3.0)
                    //     {
                    //         shape = 0.0;
                    //     }
                    //     else
                    //     {
                    //         shape = sin(3.1415926535897932384626433 / 10.0 * (coordPo(0) - 10.0)) * phi(i);

                    //     }
                    // }

                    rhs(2 * i + j) -= (r + m - shape) * weight * j0 * thickness;

                    for (int k = 0; k < connection_.size(); k++)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            if (l == 0)
                            {
                                dA_dy2(0, 0) = dphi_dxsi(0, k);
                                dA_dy2(0, 1) = dphi_dxsi(1, k);
                                dA_dy2(1, 0) = 0.0;
                                dA_dy2(1, 1) = 0.0;
                            }
                            else
                            {
                                dA_dy2(1, 0) = dphi_dxsi(0, k);
                                dA_dy2(1, 1) = dphi_dxsi(1, k);
                                dA_dy2(0, 0) = 0.0;
                                dA_dy2(0, 1) = 0.0;
                            }

                            bounded_matrix<double, 2, 2> mat3 = prod(trans(A0I), trans(dA_dy2));
                            bounded_matrix<double, 2, 2> mat4 = prod(dA_dy2, A0I);

                            bounded_matrix<double, 2, 2> dE_dy2 = 0.5 * (prod(mat3, Ac) + prod(trans(Ac), mat4)); //first derivative of E regarding k,l
                            bounded_matrix<double, 2, 2> dE_dy3 = 0.5 * (prod(mat1, mat4) + prod(mat3, mat2));    //second derivative of E regarding i,j,k,l

                            bounded_matrix<double, 2, 2> dS_dy;
                            if (ep == "EPD")
                            {
                                dS_dy(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(0, 0) + poisson * dE_dy(1, 1)));
                                dS_dy(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(1, 1) + poisson * dE_dy(0, 0)));
                                dS_dy(1, 0) = (young / (1.0 + poisson)) * dE_dy(1, 0);
                                dS_dy(0, 1) = (young / (1.0 + poisson)) * dE_dy(0, 1);
                            }
                            else
                            {
                                dS_dy(0, 0) = (young / (1.0 - poisson * poisson)) * (dE_dy(0, 0) + poisson * dE_dy(1, 1));
                                dS_dy(1, 1) = (young / (1.0 - poisson * poisson)) * (dE_dy(1, 1) + poisson * dE_dy(0, 0));
                                dS_dy(1, 0) = (young / (1.0 + poisson)) * dE_dy(1, 0);
                                dS_dy(0, 1) = (young / (1.0 + poisson)) * dE_dy(0, 1);
                            }

                            double v = dE_dy3(0, 0) * S(0, 0) + dE_dy3(1, 1) * S(1, 1) + dE_dy3(0, 1) * S(0, 1) + dE_dy3(1, 0) * S(1, 0) +                //second part of equation 5.88
                                       dE_dy2(0, 0) * dS_dy(0, 0) + dE_dy2(1, 1) * dS_dy(1, 1) + dE_dy2(0, 1) * dS_dy(0, 1) + dE_dy2(1, 0) * dS_dy(1, 0); //viscosity and pressure contribution

                            tangent(2 * i + j, 2 * k + l) += v * weight * j0 * thickness;
                            if (j == l)
                            {
                                double mm = density * phi(i) * phi(k); //mass contribution
                                tangent(2 * i + j, 2 * k + l) += (mm / (beta * deltat * deltat)) * weight * j0 * thickness;
                            }
                        }
                    }
                }
            }
        }
    }

    return std::make_pair(rhs, tangent);
}

matrix<double> Element::massMatrix(const int &hammerPoints, const int &hammerPointsBlendZone)
{
    std::vector<Cell *> cellAux;
    std::vector<int> freedom;
    int freedomDegree = connection_.size();
    for (Cell *cell : incidenceCell_)
    {
        bool teste = false;
        for (int i = 0; i < cellAux.size(); i++)
        {
            if (cell == cellAux[i])
            {
                teste = true;
            }
        }
        if (teste == false)
        {
            cellAux.push_back(cell);
            freedom.push_back(freedomDegree);
            freedomDegree += cell->getControlPoints().size();
        }
    }
    int auxiliar = 0; //auxiliar pegar distancia e cell

    matrix<double> mass(2 * freedomDegree, 2 * freedomDegree, 0.0);

    matrix<double> domainIntegrationPoints;
    if (insideBlendZone_.size() > 0) //the element has at least one point in blend zone;
    {
        domainIntegrationPoints = hammerQuadrature(hammerPointsBlendZone);
    }
    else
    {
        domainIntegrationPoints = hammerQuadrature(hammerPoints);
    }

    double auxcalc = mesh_->getMaterial()->getDensity() * mesh_->getThickness();
    double resul;
    int nlocal = connection_.size();

    if (insideBlendZone_.size() > 0)
    {
        for (int ih = 0; ih < hammerPointsBlendZone; ih++)
        {
            if (insideBlendZone_[ih] == false)
            {
                double xsi1 = domainIntegrationPoints(ih, 0);
                double xsi2 = domainIntegrationPoints(ih, 1);
                double weight = domainIntegrationPoints(ih, 2);

                vector<double> phi = domainShapeFunction(xsi1, xsi2);

                for (int i = 0; i < nlocal; i++)
                {
                    for (int k = 0; k < nlocal; k++)
                    {
                        resul = auxcalc * phi(i) * phi(k) * weight;
                        mass(2 * i, 2 * k) += resul;
                        mass(2 * i + 1, 2 * k + 1) += resul;
                    }
                }
            }
            else
            {
                Cell *cell = incidenceCell_[auxiliar];
                std::vector<ControlPoint *> cps = cell->getControlPoints();
                bounded_vector<double, 2> xsiGlobal = xsiIncidenceCell_[auxiliar];
                double blend = bValue_[auxiliar];

                const double xsi1 = domainIntegrationPoints(ih, 0);
                const double xsi2 = domainIntegrationPoints(ih, 1);
                const double weight = domainIntegrationPoints(ih, 2);

                vector<double> phiLocal = domainShapeFunction(xsi1, xsi2);

                vector<double> wpc2(cps.size());
                bounded_vector<int, 2> INC_;
                INC_ = cps[cps.size() - 1]->getINC();
                for (int i = 0; i < cps.size(); i++)
                {
                    wpc2(i) = cps[i]->getWeight();
                }
                vector<double> phiGlobal = cell->shapeFunction(xsiGlobal, wpc2, INC_);

                int nglobal = cps.size();
                int nnum = nlocal + nglobal;

                vector<double> phiBlended(nnum, 0.0);
                for (int i = 0; i < nlocal; i++)
                {
                    phiBlended(i) = (1.0 - blend) * phiLocal(i);
                }
                for (int i = 0; i < nglobal; i++)
                {
                    phiBlended(nlocal + i) = blend * phiGlobal(i);
                }

                int aux;
                for (int iaux = 0; iaux < cellAux.size(); iaux++)
                {
                    if (cell == cellAux[iaux])
                    {
                        aux = freedom[iaux]; //para auxiliar na contribuição da matriz e vetor
                    }
                }

                for (int i = 0; i < nnum; i++)
                {
                    for (int j = 0; j < nnum; j++)
                    {
                        resul = auxcalc * phiBlended(i) * phiBlended(j) * weight;

                        if (i < nlocal and j < nlocal)
                        {
                            mass(2 * i, 2 * j) += resul;
                            mass(2 * i + 1, 2 * j + 1) += resul;
                        }
                        else if (i >= nlocal and j < nlocal)
                        {
                            mass(2 * (aux + i - nlocal), 2 * j) += resul;
                            mass(2 * (aux + i - nlocal) + 1, 2 * j + 1) += resul;
                        }
                        else if (i < nlocal and j >= nlocal)
                        {
                            mass(2 * i, 2 * (aux + j - nlocal)) += resul;
                            mass(2 * i + 1, 2 * (aux + j - nlocal) + 1) += resul;
                        }
                        else if (i >= nlocal and j >= nlocal)
                        {
                            mass(2 * (aux + i - nlocal), 2 * (aux + j - nlocal)) += resul;
                            mass(2 * (aux + i - nlocal) + 1, 2 * (aux + j - nlocal) + 1) += resul;
                        }
                    }
                }

                auxiliar += 1;
            }
        }
    }
    else
    {
        for (int ih = 0; ih < hammerPoints; ih++)
        {
            double xsi1 = domainIntegrationPoints(ih, 0);
            double xsi2 = domainIntegrationPoints(ih, 1);
            double weight = domainIntegrationPoints(ih, 2);

            vector<double> phi = domainShapeFunction(xsi1, xsi2);

            for (int i = 0; i < nlocal; i++)
            {
                for (int k = 0; k < nlocal; k++)
                {
                    resul = auxcalc * phi(i) * phi(k) * weight;
                    mass(2 * i, 2 * k) += resul;
                    mass(2 * i + 1, 2 * k + 1) += resul;
                }
            }
        }
    }

    return mass;
}

bounded_vector<double, 2> Element::calculateGlobalCoordinate(const bounded_vector<double, 2> &qxsi)
{
    vector<double> phi;
    phi = domainShapeFunction(qxsi(0), qxsi(1));

    bounded_vector<double, 2> coordIP;
    coordIP(0) = 0.0;
    coordIP(1) = 0.0;

    for (int cp = 0; cp < connection_.size(); cp++) //calculando coordenadas globais dos pontos de integração
    {
        bounded_vector<double, 2> coordinateCP = connection_[cp]->getCurrentCoordinate();
        coordIP(0) += phi(cp) * coordinateCP(0);
        coordIP(1) += phi(cp) * coordinateCP(1);
    }

    return coordIP;
}

void Element::setIncidenceOnGlobal(const std::vector<bool> &inside, std::vector<Cell *> cells, const std::vector<bounded_vector<double, 2>> xsis,
                                   const std::vector<double> &bValue, const std::vector<bounded_vector<double, 2>> &db_dxsiValue)
{
    insideBlendZone_ = inside;
    incidenceCell_ = cells;
    xsiIncidenceCell_ = xsis;
    bValue_ = bValue;
    db_dxsiValue_ = db_dxsiValue;
}

// void Element::setDistanceFromFEBoundary(const std::vector<double> &distance)
// {
//     distanceFE_ = distance;
// }

// std::vector<double> Element::getDistanceFromFEBoundary()
// {
//     return distanceFE_;
// }

// double Element::blendFunction(const double &distance, const double &blendZoneThickness)
// {
//     double blend = 2.0 * pow((distance / blendZoneThickness), 3.0) - 3.0 * pow((distance / blendZoneThickness), 2.0) + 1.0;
//     return blend;
// }

std::vector<int> Element::getFreedomDegree()
{
    std::vector<int> indexNodeAndCP;

    for (Node *no : connection_)
    {
        indexNodeAndCP.push_back(no->getIndex());
    }

    std::vector<Cell *> cellAux;
    for (Cell *cell : incidenceCell_)
    {
        bool teste = false;
        for (int i = 0; i < cellAux.size(); i++)
        {
            if (cell == cellAux[i])
            {
                teste = true;
            }
        }
        if (teste == false)
        {
            cellAux.push_back(cell);
        }
    }
    for (Cell *cell : cellAux)
    {
        std::vector<ControlPoint *> cps = cell->getControlPoints();
        for (ControlPoint *cp : cps)
        {
            indexNodeAndCP.push_back(cp->getIndex());
        }
    }

    return indexNodeAndCP;
}

bounded_matrix<double, 2, 2> Element::referenceJacobianMatrixBlendZone(const matrix<double> &dphi_dxsiBlended, std::vector<ControlPoint *> cpsGLobal)
{
    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    int nLocal = connection_.size();
    for (int i = 0; i < nLocal; i++)
    {
        bounded_vector<double, 2> initialCoord = connection_[i]->getInitialCoordinate();
        dx1_dxsi1 += initialCoord(0) * dphi_dxsiBlended(0, i);
        dx1_dxsi2 += initialCoord(0) * dphi_dxsiBlended(1, i);
        dx2_dxsi1 += initialCoord(1) * dphi_dxsiBlended(0, i);
        dx2_dxsi2 += initialCoord(1) * dphi_dxsiBlended(1, i);
    }

    int nGlobal = cpsGLobal.size();
    // vector<double> wpc(nGlobal);
    // for (int i = 0; i < nGlobal; i++)
    // {
    //     wpc(i) = cpsGLobal[i]->getWeight();
    // }

    for (int i = 0; i < nGlobal; i++)
    {
        bounded_vector<double, 2> initialCoord = cpsGLobal[i]->getInitialCoordinate();
        dx1_dxsi1 += initialCoord(0) * dphi_dxsiBlended(0, i + nLocal);
        dx1_dxsi2 += initialCoord(0) * dphi_dxsiBlended(1, i + nLocal);
        dx2_dxsi1 += initialCoord(1) * dphi_dxsiBlended(0, i + nLocal);
        dx2_dxsi2 += initialCoord(1) * dphi_dxsiBlended(1, i + nLocal);
    }

    bounded_matrix<double, 2, 2> referenceJacobianMatrix;
    referenceJacobianMatrix(0, 0) = dx1_dxsi1;
    referenceJacobianMatrix(1, 0) = dx2_dxsi1;
    referenceJacobianMatrix(0, 1) = dx1_dxsi2;
    referenceJacobianMatrix(1, 1) = dx2_dxsi2;

    return referenceJacobianMatrix;
}

bounded_matrix<double, 2, 2> Element::currentJacobianMatrixBlendZone(const matrix<double> &dphi_dxsiBlended, std::vector<ControlPoint *> cpsGLobal)
{
    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    int nLocal = connection_.size();
    for (int i = 0; i < nLocal; i++)
    {
        bounded_vector<double, 2> currentCoord = connection_[i]->getCurrentCoordinate();
        dx1_dxsi1 += currentCoord(0) * dphi_dxsiBlended(0, i);
        dx1_dxsi2 += currentCoord(0) * dphi_dxsiBlended(1, i);
        dx2_dxsi1 += currentCoord(1) * dphi_dxsiBlended(0, i);
        dx2_dxsi2 += currentCoord(1) * dphi_dxsiBlended(1, i);
    }

    int nGlobal = cpsGLobal.size();
    // vector<double> wpc(nGlobal);
    // for (int i = 0; i < nGlobal; i++)
    // {
    //     wpc(i) = cpsGLobal[i]->getWeight();
    // }
    for (int i = 0; i < nGlobal; i++)
    {
        bounded_vector<double, 2> currentCoord = cpsGLobal[i]->getCurrentCoordinate();
        dx1_dxsi1 += currentCoord(0) * dphi_dxsiBlended(0, i + nLocal);
        dx1_dxsi2 += currentCoord(0) * dphi_dxsiBlended(1, i + nLocal);
        dx2_dxsi1 += currentCoord(1) * dphi_dxsiBlended(0, i + nLocal);
        dx2_dxsi2 += currentCoord(1) * dphi_dxsiBlended(1, i + nLocal);
    }

    bounded_matrix<double, 2, 2> referenceJacobianMatrix;
    referenceJacobianMatrix(0, 0) = dx1_dxsi1;
    referenceJacobianMatrix(1, 0) = dx2_dxsi1;
    referenceJacobianMatrix(0, 1) = dx1_dxsi2;
    referenceJacobianMatrix(1, 1) = dx2_dxsi2;

    return referenceJacobianMatrix;
}

std::vector<bool> Element::ipInsideBlendZone()
{
    return insideBlendZone_;
}

std::vector<bounded_vector<double, 2>> Element::xsiIncidenceCell()
{
    return xsiIncidenceCell_;
}

std::vector<Cell *> Element::incidenceCell()
{
    return incidenceCell_;
}

std::vector<double> Element::getBlendValue()
{
    return bValue_;
}

bounded_matrix<double, 2, 2> Element::inverseMatrix(const bounded_matrix<double, 2, 2> &matrix)
{
    bounded_matrix<double, 2, 2> inverse;
    double detinv = 1.0 / (matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0));

    inverse(0, 0) = detinv * matrix(1, 1);
    inverse(1, 0) = -detinv * matrix(1, 0);
    inverse(0, 1) = -detinv * matrix(0, 1);
    inverse(1, 1) = detinv * matrix(0, 0);

    return inverse;
}

vector<double> Element::diagonalMass(const int &hammerPoints, const int &hammerPointsBlendZone)
{
    std::vector<Cell *> cellAux;
    std::vector<int> freedom;
    int freedomDegree = connection_.size();
    for (Cell *cell : incidenceCell_)
    {
        bool teste = false;
        for (int i = 0; i < cellAux.size(); i++)
        {
            if (cell == cellAux[i])
            {
                teste = true;
            }
        }
        if (teste == false)
        {
            cellAux.push_back(cell);
            freedom.push_back(freedomDegree);
            freedomDegree += cell->getControlPoints().size();
        }
    }

    matrix<double> domainIntegrationPoints;
    if (insideBlendZone_.size() > 0) //the element has at least one point in blend zone;
    {
        domainIntegrationPoints = hammerQuadrature(hammerPointsBlendZone);
    }
    else
    {
        domainIntegrationPoints = hammerQuadrature(hammerPoints);
    }

    int auxiliar = 0; //auxiliar pegar distancia e cell

    vector<double> mass(freedomDegree, 0.0);
    if (insideBlendZone_.size() > 0)
    {
        for (int ih = 0; ih < hammerPointsBlendZone; ih++)
        {
            if (insideBlendZone_[ih] == false)
            {
                double xsi1 = domainIntegrationPoints(ih, 0);
                double xsi2 = domainIntegrationPoints(ih, 1);
                double weight = domainIntegrationPoints(ih, 2);

                vector<double> phi = domainShapeFunction(xsi1, xsi2);

                bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
                double j0 = jacobianDeterminant(A0);

                for (int i = 0; i < connection_.size(); i++)
                {
                    mass(i) += phi(i) * phi(i) * weight * j0;
                }
            }
            else
            {
                Cell *cell = incidenceCell_[auxiliar];
                std::vector<ControlPoint *> cps = cell->getControlPoints();
                bounded_vector<double, 2> xsiGlobal = xsiIncidenceCell_[auxiliar];
                double blend = bValue_[auxiliar];

                const double xsi1 = domainIntegrationPoints(ih, 0);
                const double xsi2 = domainIntegrationPoints(ih, 1);
                const double weight = domainIntegrationPoints(ih, 2);

                bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
                double j0 = jacobianDeterminant(A0);

                vector<double> phiLocal = domainShapeFunction(xsi1, xsi2);
                for (int i = 0; i < phiLocal.size(); i++)
                {
                    double phiblend = (1.0 - blend) * phiLocal(i);
                    mass(i) += phiblend * phiblend * weight;
                }

                vector<double> wpc2(cps.size());
                bounded_vector<int, 2> INC_;
                INC_ = cps[cps.size() - 1]->getINC();
                for (int i = 0; i < cps.size(); i++)
                {
                    wpc2(i) = cps[i]->getWeight();
                }
                vector<double> phiGlobal = cell->shapeFunction(xsiGlobal, wpc2, INC_);

                int aux;
                for (int iaux = 0; iaux < cellAux.size(); iaux++)
                {
                    if (cell == cellAux[iaux])
                    {
                        aux = freedom[iaux]; //para auxiliar na contribuição da matriz e vetor
                    }
                }

                for (int i = 0; i < phiGlobal.size(); i++)
                {
                    double phiblend = blend * phiGlobal(i);
                    mass(aux + i) += phiblend * phiblend * weight * j0;
                }

                auxiliar += 1;
            }
        }
    }
    else
    {
        for (int ih = 0; ih < hammerPoints; ih++)
        {
            double xsi1 = domainIntegrationPoints(ih, 0);
            double xsi2 = domainIntegrationPoints(ih, 1);
            double weight = domainIntegrationPoints(ih, 2);

            bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
            double j0 = jacobianDeterminant(A0);

            vector<double> phi = domainShapeFunction(xsi1, xsi2);

            for (int i = 0; i < connection_.size(); i++)
            {
                mass(i) += phi(i) * phi(i) * weight * j0;
            }
        }
    }

    return mass;
}

bounded_vector<double, 4> Element::getCauchyStress(const bounded_vector<double, 2> &qxsi, const std::string &ep)
{
    bounded_vector<double, 4> cauchStress;

    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);
    double xsi1 = qxsi(0);
    double xsi2 = qxsi(1);

    bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
    double j0 = jacobianDeterminant(A0);
    bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
    A0I(0, 0) = A0(1, 1) / j0;
    A0I(1, 1) = A0(0, 0) / j0;
    A0I(0, 1) = -A0(0, 1) / j0;
    A0I(1, 0) = -A0(1, 0) / j0;
    bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2); //current configuration map
    bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                     //current deformation gradient
    identity_matrix<double> I(2);                                        //identity matrix
    bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);   //current green strain tensor

    // Ec(0, 0) = Ac(0, 0) - 1.0;
    // Ec(1, 1) = Ac(1, 1) - 1.0;
    // Ec(0, 1) = 0.5 * (Ac(0, 1) + Ac(1, 0));
    // Ec(1, 0) = 0.5 * (Ac(0, 1) + Ac(1, 0));

    bounded_matrix<double, 2, 2> S; //second piola kirchhoff stress tensor

    if (ep == "EPD")
    {
        S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
        S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
        S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
        S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
    }
    else
    {
        S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
        S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
        S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
        S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
    }

    bounded_matrix<double, 2, 2> sigma;
    double jac = jacobianDeterminant(Ac);
    bounded_matrix<double, 2, 2> mat1;
    mat1 = prod(Ac, S);
    sigma = (1.0 / jac) * prod(mat1, trans(Ac));

    cauchStress(0) = sigma(0, 0);
    cauchStress(1) = sigma(1, 1);
    cauchStress(3) = sigma(0, 1);

    if (ep == "EPD")
    {
        cauchStress(2) = poisson * (cauchStress(0) + cauchStress(1));
    }
    else //EPT
    {
        cauchStress(2) = 0.0;
    }
    return cauchStress;
}

bounded_vector<double, 4> Element::getCauchyStressBlendZone(const bounded_vector<double, 2> &qxsiLocal, const std::string &ep, Cell *cell, const bounded_vector<double, 2> &xsiGlobal, const bounded_vector<double, 2> &blendFunctions)
{
    bounded_vector<double, 4> cauchStress;

    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);

    double blend = blendFunctions(0);
    const double xsi1 = qxsiLocal(0);
    const double xsi2 = qxsiLocal(1);

    vector<double> phiLocal = domainShapeFunction(xsi1, xsi2);                 //funções de forma calculada no ponto de hammer
    matrix<double> dphi_dxsiLocal = domainDerivativeShapeFunction(xsi1, xsi2); //derivadas das funções de forma em relação a xsi local (direção, nó);

    bounded_vector<double, 2> dDA_dxsi;
    dDA_dxsi(0) = 0.0;
    dDA_dxsi(1) = 0.0;
    for (int i = 0; i < connection_.size(); i++)
    {
        double DAnode = connection_[i]->getDistanceToBoundary();
        dDA_dxsi(0) += dphi_dxsiLocal(0, i) * DAnode; //dDAhammer_dxsiLocal1
        dDA_dxsi(1) += dphi_dxsiLocal(1, i) * DAnode; //dDAhammer_dxsiLocal1
    }

    bounded_vector<double, 2> dblend_dxsi;
    dblend_dxsi(0) = blendFunctions(1) * dDA_dxsi(0); //(db_dxsiLocal1)
    dblend_dxsi(1) = blendFunctions(1) * dDA_dxsi(1); //(db_dxsiLocal2)

    //GLOBAL
    std::vector<ControlPoint *> cps = cell->getControlPoints();
    int nglobal = cps.size(); //quantidade de pontos de controle que definem a célula em que o ponto de hammer está inserido
    vector<double> wpc(nglobal);
    bounded_vector<int, 2> inc;
    inc = cps[nglobal - 1]->getINC();
    for (int i = 0; i < nglobal; i++)
    {
        wpc(i) = cps[i]->getWeight();
    }                                                                                                                       //conectividade da célula
    const std::pair<vector<double>, matrix<double>> functionsGlobal = cell->shapeFunctionAndDerivates(xsiGlobal, wpc, inc); //functionsGlobal.first(i) são as funções de forma globais calculadas em xsiGlobal que houve incidência;
                                                                                                                            //functionsGlobal.second(j, i) são as derivadas das funções de forma globais em relação a xsi global;

    bounded_matrix<double, 2, 2> Jlocal = referenceJacobianMatrix(xsi1, xsi2);
    bounded_matrix<double, 2, 2> JGlobalI = inverseMatrix(cell->referenceJacobianMatrix(functionsGlobal.second));
    bounded_matrix<double, 2, 2> M = trans(prod(JGlobalI, Jlocal)); //para transforma dPhiGlobal_dXsiGlobal em dPhiGlobal_dXsiLocal

    int nlocal = connection_.size();
    int nnum = nlocal + nglobal;
    matrix<double> dphi_dxsiBlended(2, nnum, 0.0);

    for (int i = 0; i < nlocal; i++)
    {
        dphi_dxsiBlended(0, i) = -dblend_dxsi(0) * phiLocal(i) + (1.0 - blend) * dphi_dxsiLocal(0, i);
        dphi_dxsiBlended(1, i) = -dblend_dxsi(1) * phiLocal(i) + (1.0 - blend) * dphi_dxsiLocal(1, i);
    }

    for (int i = 0; i < nglobal; i++)
    {
        bounded_vector<double, 2> dPhiG_dXsiG, dPhiG_dxsiL;
        dPhiG_dXsiG(0) = functionsGlobal.second(0, i); //dPhiGlobal_dxsiGlobal1
        dPhiG_dXsiG(1) = functionsGlobal.second(1, i); //dPhiGlobal_dxsiGlobal2
        dPhiG_dxsiL = prod(M, dPhiG_dXsiG);

        dphi_dxsiBlended(0, nlocal + i) = dblend_dxsi(0) * functionsGlobal.first(i) + blend * dPhiG_dxsiL(0);
        dphi_dxsiBlended(1, nlocal + i) = dblend_dxsi(1) * functionsGlobal.first(i) + blend * dPhiG_dxsiL(1);
    }

    bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrixBlendZone(dphi_dxsiBlended, cps); //mapeamento da configuração inicial
    bounded_matrix<double, 2, 2> A1 = currentJacobianMatrixBlendZone(dphi_dxsiBlended, cps);   //mapeamento da configuração atual
    bounded_matrix<double, 2, 2> A0I = inverseMatrix(A0);                                      //inverse initial configuration map
    bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                                           //gradiente da função mudança de configuração
    identity_matrix<double> I(2);                                                              //identity matrix
    bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);                         //current green strain tensor
    bounded_matrix<double, 2, 2> S;                                                            //second piola kirchhoff stress tensor
    double j0 = jacobianDeterminant(Jlocal);                                                   //para multiplicar com o peso

    if (ep == "EPD")
    {
        S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
        S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
        S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
        S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
    }
    else
    {
        S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
        S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
        S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
        S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
    }

    //double jac = jacobianDeterminant(prod(currentJacobianMatrix(xsi1, xsi2), inverseMatrix(Jlocal)));

    // A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
    // //double j0 = jacobianDeterminant(A0);
    // A0I = inverseMatrix(A0); //inverse initial configuration map
    // A1 = currentJacobianMatrix(xsi1, xsi2); //current configuration map
    // Ac = prod(A1, A0I);                     //current deformation gradient

    bounded_matrix<double, 2, 2> sigma;
    double jac = jacobianDeterminant(Ac);
    bounded_matrix<double, 2, 2> mat1;
    mat1 = prod(Ac, S);
    sigma = (1.0 / jac) * prod(mat1, trans(Ac));

    cauchStress(0) = sigma(0, 0);
    cauchStress(1) = sigma(1, 1);
    cauchStress(3) = sigma(0, 1);

    if (ep == "EPD")
    {
        cauchStress(2) = poisson * (cauchStress(0) + cauchStress(1));
    }
    else //EPT
    {
        cauchStress(2) = 0.0;
    }

    return cauchStress;
}

std::vector<Node *> Element::getSideNodes(const int &side)
{
    std::vector<Node *> nodes;
    if (side == 0)
    {
        nodes.push_back(connection_[0]);
        nodes.push_back(connection_[1]);
        if (order_ == 2)
        {
            nodes.push_back(connection_[3]);
        }
        else if (order_ == 3)
        {
            nodes.push_back(connection_[3]);
            nodes.push_back(connection_[4]);
        }
    }
    else if (side == 1)
    {
        nodes.push_back(connection_[1]);
        nodes.push_back(connection_[2]);
        if (order_ == 2)
        {
            nodes.push_back(connection_[4]);
        }
        else if (order_ == 3)
        {
            nodes.push_back(connection_[5]);
            nodes.push_back(connection_[6]);
        }
    }
    else if (side == 2)
    {
        nodes.push_back(connection_[2]);
        nodes.push_back(connection_[0]);
        if (order_ == 2)
        {
            nodes.push_back(connection_[5]);
        }
        else if (order_ == 3)
        {
            nodes.push_back(connection_[7]);
            nodes.push_back(connection_[8]);
        }
    }

    return nodes;
}

bounded_vector<double, 2> Element::contributionJ_Integral4(const int &gaussPoints, const int &side, const std::string &ep, const double &rotation, Node *tipNode, const bool &crack)
{
    bounded_vector<double, 2> jIntegral;
    jIntegral(0) = 0.0;
    jIntegral(1) = 0.0;

    int nnos = connection_.size();
    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);

    double mi = young / (2.0 * (1.0 + poisson));
    double youngLine;

    double k;

    if (ep == "EPD")
    {
        k = 3.0 - 4.0 * poisson;
        youngLine = young / (1.0 - poisson * poisson);
    }
    else
    {
        k = (3.0 - poisson) / (1.0 + poisson);
        youngLine = young;
    }

    std::vector<Node *> conecOfBoundary = getSideNodes(side);
    BoundaryElement *bound = new BoundaryElement(0, conecOfBoundary);

    matrix<double> gaussIntegrationPoints(gaussPoints, 2);
    gaussIntegrationPoints = bound->boundaryIsoQuadrature(gaussPoints);

    bounded_matrix<double, 2, 2> matrixRotation;
    matrixRotation(0, 0) = cos(rotation);
    matrixRotation(0, 1) = sin(rotation);
    matrixRotation(1, 0) = -sin(rotation);
    matrixRotation(1, 1) = cos(rotation);

    bounded_vector<double, 2> tipCurrentCoordinate = tipNode->getCurrentCoordinate();
    //bounded_vector<double, 2> tipInitialCoordinate = tipNode->getInitialCoordinate();

    for (int iq = 0; iq < gaussPoints; iq++)
    {
        double xsi = gaussIntegrationPoints(iq, 0);
        double weight = gaussIntegrationPoints(iq, 1);

        matrix<double> functions = bound->shapeFunctionsAndDerivates(xsi); //phi, phi', phi''
        bounded_vector<double, 2> tangent, coordP;
        tangent(0) = 0.0;
        tangent(1) = 0.0;
        coordP(0) = 0.0;
        coordP(1) = 0.0;

        for (int ih = 0; ih < conecOfBoundary.size(); ih++)
        {
            tangent += functions(ih, 1) * conecOfBoundary[ih]->getCurrentCoordinate();
            coordP += functions(ih, 0) * conecOfBoundary[ih]->getCurrentCoordinate();
        }

        double jacobian = norm_2(tangent);

        bounded_vector<double, 2> normalaux;
        normalaux(0) = tangent(1) / jacobian;
        normalaux(1) = -tangent(0) / jacobian;

        bounded_vector<double, 2> normal = prod(matrixRotation, normalaux);

        bounded_vector<double, 2> tipToPointAux = coordP - tipCurrentCoordinate;

        double radius = norm_2(tipToPointAux);

        bounded_vector<double, 2> tipToPoint = prod(matrixRotation, tipToPointAux);

        double teta = atan2(tipToPoint(1), tipToPoint(0));

        double xsi1;
        double xsi2;

        if (side == 0)
        {
            xsi1 = (-1.0 * xsi + 1.0) / (2.0);
            xsi2 = (xsi + 1.0) / (2.0);
        }
        else if (side == 1)
        {
            xsi1 = 0.0;
            xsi2 = (-1.0 * xsi + 1.0) / (2.0);
        }
        else if (side == 2)
        {
            xsi1 = (xsi + 1.0) / (2.0);
            xsi2 = 0.0;
        }

        vector<double> phi = domainShapeFunction(xsi1, xsi2);
        matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

        bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
        bounded_matrix<double, 2, 2> A0I = inverseMatrix(A0);                  //inverse initial configuration maP
        bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2);   //current configuration map
        bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                       //current deformation gradient
        identity_matrix<double> I(2);                                          //identity matrix
        bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);     //current green strain tensor

        bounded_matrix<double, 2, 2> S; //second piola kirchhoff stress tensor

        if (ep == "EPD")
        {
            S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }
        else
        {
            S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }

        bounded_matrix<double, 2, 2> sigmaaux;
        double jac = jacobianDeterminant(Ac);
        bounded_matrix<double, 2, 2> mat1;
        mat1 = prod(Ac, S);
        sigmaaux = (1.0 / jac) * prod(mat1, trans(Ac));

        bounded_matrix<double, 2, 2> sigma;
        mat1 = prod(matrixRotation, sigmaaux);
        sigma = prod(mat1, trans(matrixRotation));

        double sigma33;

        if (ep == "EPD")
        {
            sigma33 = poisson * (sigmaaux(0, 0) + sigmaaux(1, 1));
        }
        else
        {
            sigma33 = 0.0;
        }
        bounded_matrix<double, 2, 2> epsolon, epsolonaux;

        epsolonaux(0, 0) = (1.0 / young) * (sigmaaux(0, 0) - poisson * (sigmaaux(1, 1) + sigma33));
        epsolonaux(1, 1) = (1.0 / young) * (sigmaaux(1, 1) - poisson * (sigmaaux(0, 0) + sigma33));
        epsolonaux(0, 1) = ((1.0 + poisson) / young) * sigmaaux(0, 1);
        epsolonaux(1, 0) = ((1.0 + poisson) / young) * sigmaaux(1, 0);

        mat1 = prod(matrixRotation, epsolonaux);
        epsolon = prod(mat1, trans(matrixRotation));

        bounded_vector<double, 2> force = prod(sigma, normal);

        double du1_dxsi1, du1_dxsi2, du2_dxsi1, du2_dxsi2;
        du1_dxsi1 = 0.0;
        du1_dxsi2 = 0.0;

        du2_dxsi1 = 0.0;
        du2_dxsi2 = 0.0;

        for (int in = 0; in < nnos; in++)
        {
            bounded_vector<double, 2> deslocaux = connection_[in]->getCurrentCoordinate() - connection_[in]->getInitialCoordinate();
            bounded_vector<double, 2> desloc = prod(matrixRotation, deslocaux);

            du1_dxsi1 += dphi_dxsi(0, in) * (desloc(0));
            du1_dxsi2 += dphi_dxsi(1, in) * (desloc(0));

            du2_dxsi1 += dphi_dxsi(0, in) * (desloc(1));
            du2_dxsi2 += dphi_dxsi(1, in) * (desloc(1));
        }

        double du1_dy1, du2_dy1, du1_dy2, du2_dy2;

        bounded_matrix<double, 2, 2> A1I = inverseMatrix(A1);

        du1_dy1 = du1_dxsi1 * A1I(0, 0) + du1_dxsi2 * A1I(1, 0);
        du2_dy1 = du2_dxsi1 * A1I(0, 0) + du2_dxsi2 * A1I(1, 0);

        du1_dy2 = du1_dxsi1 * A1I(0, 1) + du1_dxsi2 * A1I(1, 1);
        du2_dy2 = du2_dxsi1 * A1I(0, 1) + du2_dxsi2 * A1I(1, 1);

        double du1_dy1Chapeu, du2_dy1Chapeu; //, du1_dy2Chapeu, du2_dy2Chapeu;

        du1_dy1Chapeu = du1_dy1 * cos(rotation) + du1_dy2 * sin(rotation);

        du2_dy1Chapeu = du2_dy1 * cos(rotation) + du2_dy2 * sin(rotation);

        //AUXILIARY STATES K1 = 1 E K2 = 0;
        const double pi = 3.14159265359;
        sigmaaux(0, 0) = (1.0 / (sqrt(2.0 * pi * radius))) * (cos(0.5 * teta) * (1.0 - sin(0.5 * teta) * sin(1.5 * teta)));
        sigmaaux(1, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * (cos(0.5 * teta) * (1.0 + sin(0.5 * teta) * sin(1.5 * teta)));
        sigmaaux(0, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * sin(0.5 * teta) * cos(0.5 * teta) * cos(1.5 * teta);
        sigmaaux(1, 0) = sigmaaux(0, 1);
        if (ep == "EPD")
        {
            sigma33 = poisson * (sigmaaux(0, 0) + sigmaaux(1, 1));
        }
        else
        {
            sigma33 = 0.0;
        }

        epsolonaux(0, 0) = (1.0 / young) * (sigmaaux(0, 0) - poisson * (sigmaaux(1, 1) + sigma33));
        epsolonaux(1, 1) = (1.0 / young) * (sigmaaux(1, 1) - poisson * (sigmaaux(0, 0) + sigma33));
        epsolonaux(0, 1) = ((1.0 + poisson) / young) * sigmaaux(0, 1);
        epsolonaux(1, 0) = ((1.0 + poisson) / young) * sigmaaux(1, 0);

        double du1aux_dr, du2aux_dr, du1aux_dteta, du2aux_dteta;

        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (cos(0.5 * teta) * (k - cos(teta)));

        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k - cos(teta)));

        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (-0.5 * sin(0.5 * teta) * (k - cos(teta)) + sin(teta) * cos(0.5 * teta));

        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k - cos(teta)) + sin(teta) * sin(0.5 * teta));

        double du1aux_dy1chapeu, du2aux_dy1chapeu;

        du1aux_dy1chapeu = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;

        du2aux_dy1chapeu = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        double we = (sigma(0, 0) * epsolonaux(0, 0) + sigma(0, 1) * epsolonaux(0, 1) + sigma(1, 0) * epsolonaux(1, 0) + sigma(1, 1) * epsolonaux(1, 1));

        double we1 = (sigmaaux(0, 0) * epsolon(0, 0) + sigmaaux(0, 1) * epsolon(0, 1) + sigmaaux(1, 0) * epsolon(1, 0) + sigmaaux(1, 1) * epsolon(1, 1));

        bounded_vector<double, 2> forceAux = prod(sigmaaux, normal);

        if (crack == true)
        {
            forceAux(0) = 0.0;
            forceAux(1) = 0.0;
            force(0) = 0.0;
            force(1) = 0.0;
        }

        //double w = (sigma(0, 0) * epsolon(0, 0) + sigma(0, 1) * epsolon(0, 1) + sigma(1, 0) * epsolon(1, 0) + sigma(1, 1) * epsolon(1, 1));

        //jaIntegral(0) += (0.5 * w * normal(0) - force(0) * du1_dy1Chapeu - force(1) * du2_dy1Chapeu) * weight * jacobian;

        jIntegral(0) += (0.5 * youngLine * (0.5 * (we + we1) * normal(0) - force(0) * du1aux_dy1chapeu - force(1) * du2aux_dy1chapeu - forceAux(0) * du1_dy1Chapeu - forceAux(1) * du2_dy1Chapeu)) * weight * jacobian;

        //AUXILIARY STATE K1 = 0 E K2 = 1
        sigmaaux(0, 0) = (1.0 / (sqrt(2.0 * pi * radius))) * (-sin(0.5 * teta) * (2.0 + cos(0.5 * teta) * cos(1.5 * teta)));
        sigmaaux(1, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * (sin(0.5 * teta) * cos(0.5 * teta) * cos(1.5 * teta));
        sigmaaux(0, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * (cos(0.5 * teta) * (1.0 - sin(0.5 * teta) * sin(1.5 * teta)));
        sigmaaux(1, 0) = sigmaaux(0, 1);
        if (ep == "EPD")
        {
            sigma33 = poisson * (sigmaaux(0, 0) + sigmaaux(1, 1));
        }
        else
        {
            sigma33 = 0.0;
        }

        epsolonaux(0, 0) = (1.0 / young) * (sigmaaux(0, 0) - poisson * (sigmaaux(1, 1) + sigma33));
        epsolonaux(1, 1) = (1.0 / young) * (sigmaaux(1, 1) - poisson * (sigmaaux(0, 0) + sigma33));
        epsolonaux(0, 1) = ((1.0 + poisson) / young) * sigmaaux(0, 1);
        epsolonaux(1, 0) = ((1.0 + poisson) / young) * sigmaaux(1, 0);

        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k + 2.0 + cos(teta)));

        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (-cos(0.5 * teta) * (k - 2.0 + cos(teta)));
        ///falta daqui para baixo
        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k + 2.0 + cos(teta)) - sin(teta) * sin(0.5 * teta));

        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * sin(0.5 * teta) * (k - 2.0 + cos(teta)) + sin(teta) * cos(0.5 * teta));

        du1aux_dy1chapeu = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;

        du2aux_dy1chapeu = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        we = (sigma(0, 0) * epsolonaux(0, 0) + sigma(0, 1) * epsolonaux(0, 1) + sigma(1, 0) * epsolonaux(1, 0) + sigma(1, 1) * epsolonaux(1, 1));

        we1 = (sigmaaux(0, 0) * epsolon(0, 0) + sigmaaux(0, 1) * epsolon(0, 1) + sigmaaux(1, 0) * epsolon(1, 0) + sigmaaux(1, 1) * epsolon(1, 1));

        forceAux = prod(sigmaaux, normal);

        if (crack == true)
        {
            forceAux(0) = 0.0;
            forceAux(1) = 0.0;
            force(0) = 0.0;
            force(1) = 0.0;
        }

        jIntegral(1) += (0.5 * youngLine * (0.5 * (we + we1) * normal(0) - force(0) * du1aux_dy1chapeu - force(1) * du2aux_dy1chapeu - forceAux(0) * du1_dy1Chapeu - forceAux(1) * du2_dy1Chapeu)) * weight * jacobian;
    }
    delete bound;
    return jIntegral;
}

bounded_vector<double, 2> Element::contributionJ_IntegralInitial(const int &gaussPoints, const int &sideaux, const std::string &ep, const double &rotation, const bounded_vector<double, 2> &tipCoordinates)
{
    int side;
    if (sideaux >= 10)
    {
        side = sideaux - 11;
    }
    else
    {
        side = sideaux;
    }

    bounded_vector<double, 2> jIntegral;
    jIntegral(0) = 0.0;
    jIntegral(1) = 0.0;

    int nnos = connection_.size();
    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);
    double thickness = mesh_->getThickness();

    double mi = young / (2.0 * (1.0 + poisson));
    double youngLine;

    double k;

    if (ep == "EPD")
    {
        k = 3.0 - 4.0 * poisson;
        youngLine = thickness * young / (1.0 - poisson * poisson);
    }
    else
    {
        k = (3.0 - poisson) / (1.0 + poisson);
        youngLine = thickness * young;
    }

    std::vector<Node *> conecOfBoundary = getSideNodes(side);
    BoundaryElement *bound = new BoundaryElement(0, conecOfBoundary);

    matrix<double> gaussIntegrationPoints(gaussPoints, 2);
    gaussIntegrationPoints = bound->boundaryIsoQuadrature(gaussPoints);

    bounded_matrix<double, 2, 2> matrixRotation;
    matrixRotation(0, 0) = cos(rotation);
    matrixRotation(0, 1) = sin(rotation);
    matrixRotation(1, 0) = -sin(rotation);
    matrixRotation(1, 1) = cos(rotation);

    for (int iq = 0; iq < gaussPoints; iq++)
    {
        const double xsi = gaussIntegrationPoints(iq, 0);
        const double weight = gaussIntegrationPoints(iq, 1);

        matrix<double> functions = bound->shapeFunctionsAndDerivates(xsi); //phi, phi', phi''
        bounded_vector<double, 2> tangent, coordP;
        tangent(0) = 0.0;
        tangent(1) = 0.0;
        coordP(0) = 0.0;
        coordP(1) = 0.0;

        for (int ih = 0; ih < conecOfBoundary.size(); ih++)
        {
            tangent += functions(ih, 1) * conecOfBoundary[ih]->getInitialCoordinate();
            coordP += functions(ih, 0) * conecOfBoundary[ih]->getInitialCoordinate();
        }

        const double jacobian = norm_2(tangent);

        bounded_vector<double, 2> normalaux;
        normalaux(0) = tangent(1) / jacobian;
        normalaux(1) = -tangent(0) / jacobian;

        const bounded_vector<double, 2> normal = prod(matrixRotation, normalaux);

        bounded_vector<double, 2> tipToPointAux = coordP - tipCoordinates;
        const double radius = norm_2(tipToPointAux);

        bounded_vector<double, 2> tipToPoint = prod(matrixRotation, tipToPointAux);

        const double teta = atan2(tipToPoint(1), tipToPoint(0));

        double xsi1;
        double xsi2;

        if (side == 0)
        {
            xsi1 = (-1.0 * xsi + 1.0) / (2.0);
            xsi2 = (xsi + 1.0) / (2.0);
        }
        else if (side == 1)
        {
            xsi1 = 0.0;
            xsi2 = (-1.0 * xsi + 1.0) / (2.0);
        }
        else if (side == 2)
        {
            xsi1 = (xsi + 1.0) / (2.0);
            xsi2 = 0.0;
        }
        else
        {
            std::cout << "ERROR!" << std::endl;
        }

        vector<double> phi = domainShapeFunction(xsi1, xsi2);
        matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);

        bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2);
        bounded_matrix<double, 2, 2> A0I = inverseMatrix(A0);
        bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2);
        bounded_matrix<double, 2, 2> Ac_aux = prod(A1, A0I); //dyi / dxj

        bounded_matrix<double, 2, 2> Ac, mataux;

        mataux = prod(matrixRotation, Ac_aux);
        Ac = prod(mataux, trans(matrixRotation)); //gradiente nos eixos da ponta da fissura

        const identity_matrix<double> I(2);
        bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);
        bounded_matrix<double, 2, 2> S;

        if (ep == "EPD")
        {
            S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }
        else
        {
            S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }

        bounded_matrix<double, 2, 2> Pt;
        Pt = prod(Ac, S);
        bounded_vector<double, 2> force = prod(Pt, normal);

        const double pi = 3.14159265358979323846;

        //AUXILIARY STATES K1 = 1 E K2 = 0;
        double du1aux_dr, du2aux_dr, du1aux_dteta, du2aux_dteta;

        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (cos(0.5 * teta) * (k - cos(teta)));
        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k - cos(teta)));

        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (-0.5 * sin(0.5 * teta) * (k - cos(teta)) + sin(teta) * cos(0.5 * teta));
        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k - cos(teta)) + sin(teta) * sin(0.5 * teta));

        double du1aux_dx1c, du2aux_dx1c, du1aux_dx2c, du2aux_dx2c;

        du1aux_dx1c = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;
        du2aux_dx1c = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        du1aux_dx2c = sin(teta) * du1aux_dr + (cos(teta) / radius) * du1aux_dteta;
        du2aux_dx2c = sin(teta) * du2aux_dr + (cos(teta) / radius) * du2aux_dteta;

        bounded_matrix<double, 2, 2> Acaux, Ecaux, Saux, Ptaux;

        Acaux(0, 0) = du1aux_dx1c + 1.0;
        Acaux(0, 1) = du1aux_dx2c;
        Acaux(1, 0) = du2aux_dx1c;
        Acaux(1, 1) = du2aux_dx2c + 1.0;

        Ecaux = 0.5 * (prod(trans(Acaux), Acaux) - I);

        if (ep == "EPD")
        {
            Saux(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }
        else
        {
            Saux(0, 0) = (young / (1.0 - poisson * poisson)) * (Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / (1.0 - poisson * poisson)) * (Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }

        Ptaux = prod(Acaux, Saux);

        bounded_vector<double, 2> forceAux = prod(Ptaux, normal);

        double we = (S(0, 0) * Ecaux(0, 0) + S(0, 1) * Ecaux(0, 1) + S(1, 0) * Ecaux(1, 0) + S(1, 1) * Ecaux(1, 1));
        double we1 = (Saux(0, 0) * Ec(0, 0) + Saux(0, 1) * Ec(0, 1) + Saux(1, 0) * Ec(1, 0) + Saux(1, 1) * Ec(1, 1));

        jIntegral(0) += (0.5 * youngLine * (0.5 * (we + we1) * normal(0) - force(0) * (Acaux(0, 0) - 1.0) - force(1) * Acaux(1, 0) - forceAux(0) * (Ac(0, 0) - 1.0) - forceAux(1) * Ac(1, 0))) * weight * jacobian;

        //AUXILIARY STATE K1 = 0 E K2 = 1
        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k + 2.0 + cos(teta)));
        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (-cos(0.5 * teta) * (k - 2.0 + cos(teta)));

        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k + 2.0 + cos(teta)) - sin(teta) * sin(0.5 * teta));
        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * sin(0.5 * teta) * (k - 2.0 + cos(teta)) + sin(teta) * cos(0.5 * teta));

        du1aux_dx1c = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;
        du2aux_dx1c = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        du1aux_dx2c = sin(teta) * du1aux_dr + (cos(teta) / radius) * du1aux_dteta;
        du2aux_dx2c = sin(teta) * du2aux_dr + (cos(teta) / radius) * du2aux_dteta;

        Acaux(0, 0) = du1aux_dx1c + 1.0;
        Acaux(0, 1) = du1aux_dx2c;
        Acaux(1, 0) = du2aux_dx1c;
        Acaux(1, 1) = du2aux_dx2c + 1.0;

        Ecaux = 0.5 * (prod(trans(Acaux), Acaux) - I);

        if (ep == "EPD")
        {
            Saux(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }
        else
        {
            Saux(0, 0) = (young / (1.0 - poisson * poisson)) * (Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / (1.0 - poisson * poisson)) * (Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }

        Ptaux = prod(Acaux, Saux);

        forceAux = prod(Ptaux, normal);

        we = (S(0, 0) * Ecaux(0, 0) + S(0, 1) * Ecaux(0, 1) + S(1, 0) * Ecaux(1, 0) + S(1, 1) * Ecaux(1, 1));
        we1 = (Saux(0, 0) * Ec(0, 0) + Saux(0, 1) * Ec(0, 1) + Saux(1, 0) * Ec(1, 0) + Saux(1, 1) * Ec(1, 1));

        jIntegral(1) += (0.5 * youngLine * (0.5 * (we + we1) * normal(0) - force(0) * (Acaux(0, 0) - 1.0) - force(1) * Acaux(1, 0) - forceAux(0) * (Ac(0, 0) - 1.0) - forceAux(1) * Ac(1, 0))) * weight * jacobian;
    }
    delete bound;
    return jIntegral;
}

bounded_vector<double, 2> Element::contributionJ_Integral5(const int &gaussPoints, const int &side, const std::string &ep, const double &rotation, Node *tipNode, const bool &crack)
{
    bounded_vector<double, 2> jIntegral;
    jIntegral(0) = 0.0;
    jIntegral(1) = 0.0;

    int nnos = connection_.size();
    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);

    double mi = young / (2.0 * (1.0 + poisson));
    double youngLine;

    double k;

    if (ep == "EPD")
    {
        k = 3.0 - 4.0 * poisson;
        youngLine = young / (1.0 - poisson * poisson);
    }
    else
    {
        k = (3.0 - poisson) / (1.0 + poisson);
        youngLine = young;
    }

    std::vector<Node *> conecOfBoundary = getSideNodes(side);
    BoundaryElement *bound = new BoundaryElement(0, conecOfBoundary);

    matrix<double> gaussIntegrationPoints(gaussPoints, 2);
    gaussIntegrationPoints = bound->boundaryIsoQuadrature(gaussPoints);

    bounded_vector<double, 2> tipCurrentCoordinate = tipNode->getCurrentCoordinate();

    for (int iq = 0; iq < gaussPoints; iq++)
    {
        double xsi = gaussIntegrationPoints(iq, 0);
        double weight = gaussIntegrationPoints(iq, 1);

        matrix<double> functions = bound->shapeFunctionsAndDerivates(xsi); //phi, phi', phi''

        double xsi1;
        double xsi2;

        if (side == 0)
        {
            xsi1 = (-1.0 * xsi + 1.0) / (2.0);
            xsi2 = (xsi + 1.0) / (2.0);
        }
        else if (side == 1)
        {
            xsi1 = 0.0;
            xsi2 = (-1.0 * xsi + 1.0) / (2.0);
        }
        else if (side == 2)
        {
            xsi1 = (xsi + 1.0) / (2.0);
            xsi2 = 0.0;
        }

        vector<double> phi = domainShapeFunction(xsi1, xsi2);
        matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

        bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
        bounded_matrix<double, 2, 2> A0I = inverseMatrix(A0);                  //inverse initial configuration maP
        bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2);   //current configuration map
        bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                       //current deformation gradient

        bounded_vector<double, 2> crackDirectionInitial, crackDirectionCurrent;
        crackDirectionInitial(0) = cos(rotation);
        crackDirectionInitial(1) = sin(rotation);

        crackDirectionCurrent = prod(Ac, crackDirectionInitial);

        double alfa = atan2(crackDirectionCurrent(1), crackDirectionCurrent(0));

        bounded_matrix<double, 2, 2> matrixRotation;
        matrixRotation(0, 0) = cos(alfa);
        matrixRotation(0, 1) = sin(alfa);
        matrixRotation(1, 0) = -sin(alfa);
        matrixRotation(1, 1) = cos(alfa);

        bounded_vector<double, 2> tangent, coordP;
        tangent(0) = 0.0;
        tangent(1) = 0.0;
        coordP(0) = 0.0;
        coordP(1) = 0.0;

        for (int ih = 0; ih < conecOfBoundary.size(); ih++)
        {
            tangent += functions(ih, 1) * conecOfBoundary[ih]->getCurrentCoordinate();
            coordP += functions(ih, 0) * conecOfBoundary[ih]->getCurrentCoordinate();
        }

        double jacobian = norm_2(tangent);

        bounded_vector<double, 2> normalaux;
        normalaux(0) = tangent(1) / jacobian;
        normalaux(1) = -tangent(0) / jacobian; //normal em x1 e x2

        bounded_vector<double, 2> normal = prod(matrixRotation, normalaux); //normal em y1c e y2c

        bounded_vector<double, 2> tipToPointAux = coordP - tipCurrentCoordinate;

        double radius = norm_2(tipToPointAux);

        bounded_vector<double, 2> tipToPoint = prod(matrixRotation, tipToPointAux);

        double teta = atan2(tipToPoint(1), tipToPoint(0));

        identity_matrix<double> I(2);                                      //identity matrix
        bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor

        bounded_matrix<double, 2, 2> S; //second piola kirchhoff stress tensor

        if (ep == "EPD")
        {
            S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }
        else
        {
            S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }

        bounded_matrix<double, 2, 2> sigmaaux;
        double jac = jacobianDeterminant(Ac);
        bounded_matrix<double, 2, 2> mat1;
        mat1 = prod(Ac, S);
        sigmaaux = (1.0 / jac) * prod(mat1, trans(Ac));

        bounded_matrix<double, 2, 2> sigma;
        mat1 = prod(matrixRotation, sigmaaux);
        sigma = prod(mat1, trans(matrixRotation));

        double sigma33;

        if (ep == "EPD")
        {
            sigma33 = poisson * (sigmaaux(0, 0) + sigmaaux(1, 1));
        }
        else
        {
            sigma33 = 0.0;
        }
        bounded_matrix<double, 2, 2> epsolon, epsolonaux;

        epsolonaux(0, 0) = (1.0 / young) * (sigmaaux(0, 0) - poisson * (sigmaaux(1, 1) + sigma33));
        epsolonaux(1, 1) = (1.0 / young) * (sigmaaux(1, 1) - poisson * (sigmaaux(0, 0) + sigma33));
        epsolonaux(0, 1) = ((1.0 + poisson) / young) * sigmaaux(0, 1);
        epsolonaux(1, 0) = ((1.0 + poisson) / young) * sigmaaux(1, 0);

        mat1 = prod(matrixRotation, epsolonaux);
        epsolon = prod(mat1, trans(matrixRotation));

        bounded_vector<double, 2> force = prod(sigma, normal);

        double du1_dxsi1, du1_dxsi2, du2_dxsi1, du2_dxsi2;
        du1_dxsi1 = 0.0;
        du1_dxsi2 = 0.0;

        du2_dxsi1 = 0.0;
        du2_dxsi2 = 0.0;

        for (int in = 0; in < nnos; in++)
        {
            bounded_vector<double, 2> deslocaux = connection_[in]->getCurrentCoordinate() - connection_[in]->getInitialCoordinate();
            bounded_vector<double, 2> desloc = prod(matrixRotation, deslocaux);

            du1_dxsi1 += dphi_dxsi(0, in) * (desloc(0));
            du1_dxsi2 += dphi_dxsi(1, in) * (desloc(0));

            du2_dxsi1 += dphi_dxsi(0, in) * (desloc(1));
            du2_dxsi2 += dphi_dxsi(1, in) * (desloc(1));
        }

        double du1_dy1, du2_dy1, du1_dy2, du2_dy2;

        bounded_matrix<double, 2, 2> A1I = inverseMatrix(A1);

        du1_dy1 = du1_dxsi1 * A1I(0, 0) + du1_dxsi2 * A1I(1, 0);
        du2_dy1 = du2_dxsi1 * A1I(0, 0) + du2_dxsi2 * A1I(1, 0);

        du1_dy2 = du1_dxsi1 * A1I(0, 1) + du1_dxsi2 * A1I(1, 1);
        du2_dy2 = du2_dxsi1 * A1I(0, 1) + du2_dxsi2 * A1I(1, 1);

        double du1_dy1Chapeu, du2_dy1Chapeu; //, du1_dy2Chapeu, du2_dy2Chapeu;

        du1_dy1Chapeu = du1_dy1 * cos(rotation) + du1_dy2 * sin(rotation);

        du2_dy1Chapeu = du2_dy1 * cos(rotation) + du2_dy2 * sin(rotation);

        //AUXILIARY STATES K1 = 1 E K2 = 0;
        const double pi = 3.14159265359;
        sigmaaux(0, 0) = (1.0 / (sqrt(2.0 * pi * radius))) * (cos(0.5 * teta) * (1.0 - sin(0.5 * teta) * sin(1.5 * teta)));
        sigmaaux(1, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * (cos(0.5 * teta) * (1.0 + sin(0.5 * teta) * sin(1.5 * teta)));
        sigmaaux(0, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * sin(0.5 * teta) * cos(0.5 * teta) * cos(1.5 * teta);
        sigmaaux(1, 0) = sigmaaux(0, 1);
        if (ep == "EPD")
        {
            sigma33 = poisson * (sigmaaux(0, 0) + sigmaaux(1, 1));
        }
        else
        {
            sigma33 = 0.0;
        }

        epsolonaux(0, 0) = (1.0 / young) * (sigmaaux(0, 0) - poisson * (sigmaaux(1, 1) + sigma33));
        epsolonaux(1, 1) = (1.0 / young) * (sigmaaux(1, 1) - poisson * (sigmaaux(0, 0) + sigma33));
        epsolonaux(0, 1) = ((1.0 + poisson) / young) * sigmaaux(0, 1);
        epsolonaux(1, 0) = ((1.0 + poisson) / young) * sigmaaux(1, 0);

        double du1aux_dr, du2aux_dr, du1aux_dteta, du2aux_dteta;

        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (cos(0.5 * teta) * (k - cos(teta)));

        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k - cos(teta)));

        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (-0.5 * sin(0.5 * teta) * (k - cos(teta)) + sin(teta) * cos(0.5 * teta));

        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k - cos(teta)) + sin(teta) * sin(0.5 * teta));

        double du1aux_dy1chapeu, du2aux_dy1chapeu;

        du1aux_dy1chapeu = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;

        du2aux_dy1chapeu = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        double we = (sigma(0, 0) * epsolonaux(0, 0) + sigma(0, 1) * epsolonaux(0, 1) + sigma(1, 0) * epsolonaux(1, 0) + sigma(1, 1) * epsolonaux(1, 1));

        double we1 = (sigmaaux(0, 0) * epsolon(0, 0) + sigmaaux(0, 1) * epsolon(0, 1) + sigmaaux(1, 0) * epsolon(1, 0) + sigmaaux(1, 1) * epsolon(1, 1));

        bounded_vector<double, 2> forceAux = prod(sigmaaux, normal);

        if (crack == true)
        {
            forceAux(0) = 0.0;
            forceAux(1) = 0.0;
            force(0) = 0.0;
            force(1) = 0.0;
        }

        //double w = (sigma(0, 0) * epsolon(0, 0) + sigma(0, 1) * epsolon(0, 1) + sigma(1, 0) * epsolon(1, 0) + sigma(1, 1) * epsolon(1, 1));

        //jaIntegral(0) += (0.5 * w * normal(0) - force(0) * du1_dy1Chapeu - force(1) * du2_dy1Chapeu) * weight * jacobian;

        jIntegral(0) += (0.5 * youngLine * (0.5 * (we + we1) * normal(0) - force(0) * du1aux_dy1chapeu - force(1) * du2aux_dy1chapeu - forceAux(0) * du1_dy1Chapeu - forceAux(1) * du2_dy1Chapeu)) * weight * jacobian;

        //AUXILIARY STATE K1 = 0 E K2 = 1
        sigmaaux(0, 0) = (1.0 / (sqrt(2.0 * pi * radius))) * (-sin(0.5 * teta) * (2.0 + cos(0.5 * teta) * cos(1.5 * teta)));
        sigmaaux(1, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * (sin(0.5 * teta) * cos(0.5 * teta) * cos(1.5 * teta));
        sigmaaux(0, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * (cos(0.5 * teta) * (1.0 - sin(0.5 * teta) * sin(1.5 * teta)));
        sigmaaux(1, 0) = sigmaaux(0, 1);
        if (ep == "EPD")
        {
            sigma33 = poisson * (sigmaaux(0, 0) + sigmaaux(1, 1));
        }
        else
        {
            sigma33 = 0.0;
        }

        epsolonaux(0, 0) = (1.0 / young) * (sigmaaux(0, 0) - poisson * (sigmaaux(1, 1) + sigma33));
        epsolonaux(1, 1) = (1.0 / young) * (sigmaaux(1, 1) - poisson * (sigmaaux(0, 0) + sigma33));
        epsolonaux(0, 1) = ((1.0 + poisson) / young) * sigmaaux(0, 1);
        epsolonaux(1, 0) = ((1.0 + poisson) / young) * sigmaaux(1, 0);

        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k + 2.0 + cos(teta)));

        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (-cos(0.5 * teta) * (k - 2.0 + cos(teta)));
        ///falta daqui para baixo
        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k + 2.0 + cos(teta)) - sin(teta) * sin(0.5 * teta));

        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * sin(0.5 * teta) * (k - 2.0 + cos(teta)) + sin(teta) * cos(0.5 * teta));

        du1aux_dy1chapeu = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;

        du2aux_dy1chapeu = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        we = (sigma(0, 0) * epsolonaux(0, 0) + sigma(0, 1) * epsolonaux(0, 1) + sigma(1, 0) * epsolonaux(1, 0) + sigma(1, 1) * epsolonaux(1, 1));

        we1 = (sigmaaux(0, 0) * epsolon(0, 0) + sigmaaux(0, 1) * epsolon(0, 1) + sigmaaux(1, 0) * epsolon(1, 0) + sigmaaux(1, 1) * epsolon(1, 1));

        forceAux = prod(sigmaaux, normal);

        if (crack == true)
        {
            forceAux(0) = 0.0;
            forceAux(1) = 0.0;
            force(0) = 0.0;
            force(1) = 0.0;
        }

        jIntegral(1) += (0.5 * youngLine * (0.5 * (we + we1) * normal(0) - force(0) * du1aux_dy1chapeu - force(1) * du2aux_dy1chapeu - forceAux(0) * du1_dy1Chapeu - forceAux(1) * du2_dy1Chapeu)) * weight * jacobian;
    }
    delete bound;
    return jIntegral;
}

bounded_vector<double, 2> Element::contributionJ_IntegralFromRice(const int &gaussPoints, const int &sideaux, const std::string &ep, const double &rotation, Node *tipNode)
{
    int side;
    if (sideaux >= 10)
    {
        side = sideaux - 11;
    }
    else
    {
        side = sideaux;
    }

    bounded_vector<double, 2> jIntegral;
    jIntegral(0) = 0.0;
    jIntegral(1) = 0.0;

    int nnos = connection_.size();
    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);

    double mi = young / (2.0 * (1.0 + poisson));
    double youngLine;
    double thickness = mesh_->getThickness();

    double k;

    if (ep == "EPD")
    {
        k = 3.0 - 4.0 * poisson;
        youngLine = thickness * young / (1.0 - poisson * poisson);
    }
    else
    {
        k = (3.0 - poisson) / (1.0 + poisson);
        youngLine = thickness * young;
    }

    std::vector<Node *> conecOfBoundary = getSideNodes(side);
    BoundaryElement *bound = new BoundaryElement(0, conecOfBoundary);

    matrix<double> gaussIntegrationPoints(gaussPoints, 2);
    gaussIntegrationPoints = bound->boundaryIsoQuadrature(gaussPoints);

    bounded_matrix<double, 2, 2> matrixRotation;
    matrixRotation(0, 0) = cos(rotation);
    matrixRotation(0, 1) = sin(rotation);
    matrixRotation(1, 0) = -sin(rotation);
    matrixRotation(1, 1) = cos(rotation);

    // bounded_vector<double, 2> tipCurrentCoordinate = tipNode->getCurrentCoordinate();
    // bounded_vector<double, 2> tipInitialCoordinate = tipNode->getInitialCoordinate();

    for (int iq = 0; iq < gaussPoints; iq++)
    {
        double xsi = gaussIntegrationPoints(iq, 0);
        double weight = gaussIntegrationPoints(iq, 1);

        matrix<double> functions = bound->shapeFunctionsAndDerivates(xsi); //phi, phi', phi''
        bounded_vector<double, 2> tangent, coordP;
        tangent(0) = 0.0;
        tangent(1) = 0.0;
        coordP(0) = 0.0;
        coordP(1) = 0.0;

        for (int ih = 0; ih < conecOfBoundary.size(); ih++)
        {
            tangent += functions(ih, 1) * conecOfBoundary[ih]->getInitialCoordinate();
            coordP += functions(ih, 0) * conecOfBoundary[ih]->getInitialCoordinate();
        }

        double jacobian = norm_2(tangent);

        bounded_vector<double, 2> normalaux;
        normalaux(0) = tangent(1) / jacobian;
        normalaux(1) = -tangent(0) / jacobian;

        bounded_vector<double, 2> normal = prod(matrixRotation, normalaux);

        bounded_vector<double, 2> tipToPointAux = coordP - tipNode->getInitialCoordinate();

        double radius = norm_2(tipToPointAux);

        bounded_vector<double, 2> tipToPoint = prod(matrixRotation, tipToPointAux);

        double teta = atan2(tipToPoint(1), tipToPoint(0));

        double xsi1;
        double xsi2;

        if (side == 0)
        {
            xsi1 = (-1.0 * xsi + 1.0) / (2.0);
            xsi2 = (xsi + 1.0) / (2.0);
        }
        else if (side == 1)
        {
            xsi1 = 0.0;
            xsi2 = (-1.0 * xsi + 1.0) / (2.0);
        }
        else if (side == 2)
        {
            xsi1 = (xsi + 1.0) / (2.0);
            xsi2 = 0.0;
        }
        else
        {
            std::cout << "ERROR!" << std::endl;
        }

        vector<double> phi = domainShapeFunction(xsi1, xsi2);
        matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);

        bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2);
        bounded_matrix<double, 2, 2> A0I = inverseMatrix(A0);
        bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2);
        bounded_matrix<double, 2, 2> Ac_aux = prod(A1, A0I); //dyi / dxj

        bounded_matrix<double, 2, 2> Ac, mataux;

        mataux = prod(matrixRotation, Ac_aux);
        Ac = prod(mataux, trans(matrixRotation)); //gradiente nos eixos da ponta da fissura

        const identity_matrix<double> I(2); //identity matrix
        bounded_matrix<double, 2, 2> Ec;    // = 0.5 * (prod(trans(Ac), Ac) - I);     //current green strain tensor

        Ec(0, 0) = Ac(0, 0) - 1.0;
        Ec(1, 1) = Ac(1, 1) - 1.0;
        Ec(0, 1) = 0.5 * (Ac(0, 1) + Ac(1, 0));
        Ec(1, 0) = 0.5 * (Ac(0, 1) + Ac(1, 0));

        bounded_matrix<double, 2, 2> S; //second piola kirchhoff stress tensor

        if (ep == "EPD")
        {
            S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }
        else
        {
            S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }

        bounded_vector<double, 2> force = prod(S, normal);

        //AUXILIARY STATES K1 = 1 E K2 = 0;
        const double pi = 3.14159265359;
        bounded_matrix<double, 2, 2> sigmaaux, epsolonaux;
        sigmaaux(0, 0) = (1.0 / (sqrt(2.0 * pi * radius))) * (cos(0.5 * teta) * (1.0 - sin(0.5 * teta) * sin(1.5 * teta)));
        sigmaaux(1, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * (cos(0.5 * teta) * (1.0 + sin(0.5 * teta) * sin(1.5 * teta)));
        sigmaaux(0, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * sin(0.5 * teta) * cos(0.5 * teta) * cos(1.5 * teta);
        sigmaaux(1, 0) = sigmaaux(0, 1);
        double sigma33;
        if (ep == "EPD")
        {
            sigma33 = poisson * (sigmaaux(0, 0) + sigmaaux(1, 1));
        }
        else
        {
            sigma33 = 0.0;
        }

        epsolonaux(0, 0) = (1.0 / young) * (sigmaaux(0, 0) - poisson * (sigmaaux(1, 1) + sigma33));
        epsolonaux(1, 1) = (1.0 / young) * (sigmaaux(1, 1) - poisson * (sigmaaux(0, 0) + sigma33));
        epsolonaux(0, 1) = ((1.0 + poisson) / young) * sigmaaux(0, 1);
        epsolonaux(1, 0) = ((1.0 + poisson) / young) * sigmaaux(1, 0);

        double du1aux_dr, du2aux_dr, du1aux_dteta, du2aux_dteta;

        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (cos(0.5 * teta) * (k - cos(teta)));

        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k - cos(teta)));

        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (-0.5 * sin(0.5 * teta) * (k - cos(teta)) + sin(teta) * cos(0.5 * teta));

        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k - cos(teta)) + sin(teta) * sin(0.5 * teta));

        double du1aux_dy1chapeu, du2aux_dy1chapeu;

        du1aux_dy1chapeu = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;

        du2aux_dy1chapeu = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        double we = (S(0, 0) * epsolonaux(0, 0) + S(0, 1) * epsolonaux(0, 1) + S(1, 0) * epsolonaux(1, 0) + S(1, 1) * epsolonaux(1, 1));

        double we1 = (sigmaaux(0, 0) * Ec(0, 0) + sigmaaux(0, 1) * Ec(0, 1) + sigmaaux(1, 0) * Ec(1, 0) + sigmaaux(1, 1) * Ec(1, 1));

        bounded_vector<double, 2> forceAux = prod(sigmaaux, normal);

        jIntegral(0) += (0.5 * youngLine * (0.5 * (we + we1) * normal(0) - force(0) * du1aux_dy1chapeu - force(1) * du2aux_dy1chapeu - forceAux(0) * (Ac(0, 0) - 1.0) - forceAux(1) * Ac(1, 0))) * weight * jacobian;

        //AUXILIARY STATE K1 = 0 E K2 = 1
        sigmaaux(0, 0) = (1.0 / (sqrt(2.0 * pi * radius))) * (-sin(0.5 * teta) * (2.0 + cos(0.5 * teta) * cos(1.5 * teta)));
        sigmaaux(1, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * (sin(0.5 * teta) * cos(0.5 * teta) * cos(1.5 * teta));
        sigmaaux(0, 1) = (1.0 / (sqrt(2.0 * pi * radius))) * (cos(0.5 * teta) * (1.0 - sin(0.5 * teta) * sin(1.5 * teta)));
        sigmaaux(1, 0) = sigmaaux(0, 1);
        if (ep == "EPD")
        {
            sigma33 = poisson * (sigmaaux(0, 0) + sigmaaux(1, 1));
        }
        else
        {
            sigma33 = 0.0;
        }

        epsolonaux(0, 0) = (1.0 / young) * (sigmaaux(0, 0) - poisson * (sigmaaux(1, 1) + sigma33));
        epsolonaux(1, 1) = (1.0 / young) * (sigmaaux(1, 1) - poisson * (sigmaaux(0, 0) + sigma33));
        epsolonaux(0, 1) = ((1.0 + poisson) / young) * sigmaaux(0, 1);
        epsolonaux(1, 0) = ((1.0 + poisson) / young) * sigmaaux(1, 0);

        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k + 2.0 + cos(teta)));

        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (-cos(0.5 * teta) * (k - 2.0 + cos(teta)));
        ///falta daqui para baixo
        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k + 2.0 + cos(teta)) - sin(teta) * sin(0.5 * teta));

        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * sin(0.5 * teta) * (k - 2.0 + cos(teta)) + sin(teta) * cos(0.5 * teta));

        du1aux_dy1chapeu = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;

        du2aux_dy1chapeu = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        we = (S(0, 0) * epsolonaux(0, 0) + S(0, 1) * epsolonaux(0, 1) + S(1, 0) * epsolonaux(1, 0) + S(1, 1) * epsolonaux(1, 1));

        we1 = (sigmaaux(0, 0) * Ec(0, 0) + sigmaaux(0, 1) * Ec(0, 1) + sigmaaux(1, 0) * Ec(1, 0) + sigmaaux(1, 1) * Ec(1, 1));

        forceAux = prod(sigmaaux, normal);

        jIntegral(1) += (0.5 * youngLine * (0.5 * (we + we1) * normal(0) - force(0) * du1aux_dy1chapeu - force(1) * du2aux_dy1chapeu - forceAux(0) * (Ac(0, 0) - 1.0) - forceAux(1) * Ac(1, 0))) * weight * jacobian;

        // if (crack == true)
        // {
        //     forceAux(0) = 0.0;
        //     forceAux(1) = 0.0;
        //     force(0) = 0.0;
        //     force(1) = 0.0;
        // }

        // jIntegral(1) += (0.5 * youngLine * (0.5 * (we + we1) * normal(0) - force(0) * du1aux_dy1chapeu - force(1) * du2aux_dy1chapeu - forceAux(0) * du1_dy1Chapeu - forceAux(1) * du2_dy1Chapeu)) * weight * jacobian;
    }
    delete bound;
    return jIntegral;
}

bounded_vector<double, 2> Element::domainJ_IntegralInitial(const int &nHammer, const std::string &ep, const std::string &type, const double &Jintegral_radius, const double &rotation, const bounded_vector<double, 2> &tipCoordinates)
{
    bounded_vector<double, 2> jIntegral;
    jIntegral(0) = 0.0;
    jIntegral(1) = 0.0;

    int nnos = connection_.size();

    //Material properties
    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);
    double thickness = mesh_->getThickness();
    double mi = young / (2.0 * (1.0 + poisson));
    double youngLine;
    double k;
    if (ep == "EPD")
    {
        k = 3.0 - 4.0 * poisson;
        youngLine = thickness * young / (1.0 - poisson * poisson);
    }
    else
    {
        k = (3.0 - poisson) / (1.0 + poisson);
        youngLine = thickness * young;
    }

    matrix<double> domainIntegrationPoints = hammerQuadrature(nHammer);

    bounded_matrix<double, 2, 2> matrixRotation;
    matrixRotation(0, 0) = cos(rotation);
    matrixRotation(0, 1) = sin(rotation);
    matrixRotation(1, 0) = -sin(rotation);
    matrixRotation(1, 1) = cos(rotation);

    for (int ihh = 0; ihh < nHammer; ihh++)
    {

        double xsi1 = domainIntegrationPoints(ihh, 0);
        double xsi2 = domainIntegrationPoints(ihh, 1);
        double weight = domainIntegrationPoints(ihh, 2);

        vector<double> phi = domainShapeFunction(xsi1, xsi2);
        matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node
        matrix<double> d2phi_d2xsi = domainSecondDerivativeShapeFunction(xsi1, xsi2);

        bounded_matrix<double, 2, 2> dA1_dxsi1, dA1_dxsi2, A0loc, A1loc, mataux;
        dA1_dxsi1(0, 0) = 0.0;
        dA1_dxsi1(0, 1) = 0.0;
        dA1_dxsi1(1, 0) = 0.0;
        dA1_dxsi1(1, 1) = 0.0;

        dA1_dxsi2(0, 0) = 0.0;
        dA1_dxsi2(0, 1) = 0.0;
        dA1_dxsi2(1, 0) = 0.0;
        dA1_dxsi2(1, 1) = 0.0;

        A0loc(0, 0) = 0.0;
        A0loc(0, 1) = 0.0;
        A0loc(1, 0) = 0.0;
        A0loc(1, 1) = 0.0;

        A1loc(0, 0) = 0.0;
        A1loc(0, 1) = 0.0;
        A1loc(1, 0) = 0.0;
        A1loc(1, 1) = 0.0;

        for (int in = 0; in < nnos; in++)
        {
            bounded_vector<double, 2> coordNode = prod(matrixRotation, connection_[in]->getCurrentCoordinate());

            bounded_vector<double, 2> initialNode = prod(matrixRotation, connection_[in]->getInitialCoordinate());

            dA1_dxsi1(0, 0) += d2phi_d2xsi(0, in) * coordNode(0);
            dA1_dxsi1(0, 1) += d2phi_d2xsi(1, in) * coordNode(0);
            dA1_dxsi1(1, 0) += d2phi_d2xsi(0, in) * coordNode(1);
            dA1_dxsi1(1, 1) += d2phi_d2xsi(1, in) * coordNode(1);

            dA1_dxsi2(0, 0) += d2phi_d2xsi(1, in) * coordNode(0);
            dA1_dxsi2(0, 1) += d2phi_d2xsi(2, in) * coordNode(0);
            dA1_dxsi2(1, 0) += d2phi_d2xsi(1, in) * coordNode(1);
            dA1_dxsi2(1, 1) += d2phi_d2xsi(2, in) * coordNode(1);

            A0loc(0, 0) += initialNode(0) * dphi_dxsi(0, in);
            A0loc(0, 1) += initialNode(0) * dphi_dxsi(1, in);
            A0loc(1, 0) += initialNode(1) * dphi_dxsi(0, in);
            A0loc(1, 1) += initialNode(1) * dphi_dxsi(1, in);

            // A1loc(0, 0) += coordNode(0) * dphi_dxsi(0, in);
            // A1loc(0, 1) += coordNode(0) * dphi_dxsi(1, in);
            // A1loc(1, 0) += coordNode(1) * dphi_dxsi(0, in);
            // A1loc(1, 1) += coordNode(1) * dphi_dxsi(1, in);
        }

        bounded_matrix<double, 2, 2> A0glob = referenceJacobianMatrix(xsi1, xsi2);
        double j0 = jacobianDeterminant(A0glob);
        bounded_matrix<double, 2, 2> A0Iglob = inverseMatrix(A0glob);
        bounded_matrix<double, 2, 2> A1glob = currentJacobianMatrix(xsi1, xsi2);
        bounded_matrix<double, 2, 2> Acglob = prod(A1glob, A0Iglob); //dyi / dxj

        bounded_matrix<double, 2, 2> Ac; // Aclocteste;

        //Aclocteste = prod(A1loc, inverseMatrix(A0loc));

        mataux = prod(matrixRotation, Acglob);
        Ac = prod(mataux, trans(matrixRotation)); //gradiente nos eixos da ponta da fissura

        // std::cout << Aclocteste(0, 0) / Ac(0, 0) << " " << Aclocteste(0, 1) / Ac(0, 1) << " " << Aclocteste(1, 0) / Ac(1, 0) << " " << Aclocteste(1, 1) / Ac(1, 1) << std::endl;

        const identity_matrix<double> I(2);
        bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);
        bounded_matrix<double, 2, 2> S;

        if (ep == "EPD")
        {
            S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }
        else
        {
            S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }

        bounded_matrix<double, 2, 2> Pt;
        Pt = prod(Ac, S);

        bounded_matrix<double, 2, 2> dAc_dxsi1, dAc_dxsi2, A0Iloc;

        A0Iloc = inverseMatrix(A0loc);
        dAc_dxsi1 = prod(dA1_dxsi1, A0Iloc);
        dAc_dxsi2 = prod(dA1_dxsi2, A0Iloc);

        bounded_vector<double, 2> a;
        a(0) = Ac(0, 0) - 1.0;
        a(1) = Ac(1, 0);

        bounded_matrix<double, 2, 2> grad_a, dEc_dx1, dAc_dx1;

        grad_a(0, 0) = dAc_dxsi1(0, 0) * A0Iloc(0, 0) + dAc_dxsi2(0, 0) * A0Iloc(1, 0); //ok
        grad_a(0, 1) = dAc_dxsi1(0, 0) * A0Iloc(0, 1) + dAc_dxsi2(0, 0) * A0Iloc(1, 1);
        grad_a(1, 0) = dAc_dxsi1(1, 0) * A0Iloc(0, 0) + dAc_dxsi2(1, 0) * A0Iloc(1, 0); //ok
        grad_a(1, 1) = dAc_dxsi1(1, 0) * A0Iloc(0, 1) + dAc_dxsi2(1, 0) * A0Iloc(1, 1);

        dAc_dx1 = dAc_dxsi1 * A0Iloc(0, 0) + dAc_dxsi2 * A0Iloc(1, 0);

        dEc_dx1 = 0.5 * (prod(trans(dAc_dx1), Ac) + prod(trans(Ac), dAc_dx1));

        //////////////////////////////////
        //Q FUNCTION
        //////////////////////////////////
        double radius = 0.0;
        bounded_vector<double, 2> dradius_dxsi, coordP;
        coordP(0) = 0.0;
        coordP(1) = 0.0;
        dradius_dxsi(0) = 0.0;
        dradius_dxsi(1) = 0.0;

        for (int i = 0; i < connection_.size(); i++)
        {
            coordP += phi(i) * connection_[i]->getInitialCoordinate();

            const double radiusNODE = norm_2(connection_[i]->getInitialCoordinate() - tipCoordinates);

            dradius_dxsi(0) += dphi_dxsi(0, i) * radiusNODE; //dDAhammer_dxsiLocal1
            dradius_dxsi(1) += dphi_dxsi(1, i) * radiusNODE; //dDAhammer_dxsiLocal1

            radius += phi(i) * radiusNODE;
        }

        bounded_vector<double, 2> tipToPointAux = coordP - tipCoordinates;
        bounded_vector<double, 2> tipToPoint = prod(matrixRotation, tipToPointAux);
        const double teta = atan2(tipToPoint(1), tipToPoint(0));

        const double aux = radius / Jintegral_radius;
        const double q = 2.0 * aux * aux * aux - 3.0 * aux * aux + 1.0;
        const double dq_dradius = (6.0 / Jintegral_radius) * aux * aux - (6.0 / Jintegral_radius) * aux;

        bounded_vector<double, 2> dq_dxsi;
        dq_dxsi(0) = dq_dradius * dradius_dxsi(0);
        dq_dxsi(1) = dq_dradius * dradius_dxsi(1);

        bounded_vector<double, 2> grad_q;

        grad_q(0) = dq_dxsi(0) * A0Iloc(0, 0) + dq_dxsi(1) * A0Iloc(1, 0);
        grad_q(1) = dq_dxsi(0) * A0Iloc(0, 1) + dq_dxsi(1) * A0Iloc(1, 1);

        //////////////////////////////////
        //AUXILIARY STATES K1 = 1 E K2 = 0;
        //////////////////////////////////

        const double pi = 3.14159265358979323846;
        const double meio_teta = 0.5 * teta;

        double dl_dr, d2l_dr2, dm_dtheta, d2m_dtetha2, dn_dtheta, d2n_dtheta2, K1, K2, l, m, n;

        K1 = 1.0;
        K2 = 0.0;

        l = (sqrt(radius / (2.0 * pi)) / (2.0 * mi));
        dl_dr = sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius);
        d2l_dr2 = -sqrt(2.0 * pi * radius) / (16.0 * mi * pi * radius * radius);

        m = K1 * cos(meio_teta) * (k - cos(teta)) + K2 * sin(meio_teta) * (k + 2.0 + cos(teta));
        dm_dtheta = K1 * (-0.5 * sin(meio_teta) * (k - cos(teta)) + sin(teta) * cos(meio_teta)) + K2 * (0.5 * cos(meio_teta) * (k + 2.0 + cos(teta)) - sin(teta) * sin(meio_teta));
        d2m_dtetha2 = K1 * (-cos(meio_teta) * 0.25 * (k - cos(teta)) - sin(meio_teta) * sin(teta) + cos(meio_teta) * cos(teta)) + K2 * (-sin(meio_teta) * 0.25 * (k + 2.0 + cos(teta)) - sin(teta) * cos(meio_teta) - cos(teta) * sin(meio_teta));

        n = K1 * sin(meio_teta) * (k - cos(teta)) - K2 * cos(meio_teta) * (k - 2.0 + cos(teta));
        dn_dtheta = K1 * (0.5 * cos(meio_teta) * (k - cos(teta)) + sin(teta) * sin(meio_teta)) + K2 * (0.5 * sin(meio_teta) * (k - 2.0 + cos(teta)) + sin(teta) * cos(meio_teta));
        d2n_dtheta2 = K1 * (-sin(meio_teta) * 0.25 * (k - cos(teta)) + cos(meio_teta) * sin(teta) + sin(meio_teta) * cos(teta)) - K2 * (-cos(meio_teta) * 0.25 * (k - 2.0 + cos(teta)) + sin(meio_teta) * sin(teta) - cos(teta) * cos(meio_teta));

        double du1aux_dr, du2aux_dr, du1aux_dteta, du2aux_dteta, d2u1aux_dr2, d2u1aux_drdteta, d2u1aux_dteta2, d2u2aux_dr2, d2u2aux_drdteta, d2u2aux_dteta2;

        //FIRST DERIVATES OF U1
        du1aux_dr = dl_dr * m;
        du1aux_dteta = l * dm_dtheta;

        //SECOND DERIVATES OF U1
        d2u1aux_dr2 = d2l_dr2 * m;
        d2u1aux_drdteta = dl_dr * dm_dtheta;
        d2u1aux_dteta2 = l * d2m_dtetha2;

        //FIRST DERIVATES OF U2
        du2aux_dr = dl_dr * n;
        du2aux_dteta = l * dn_dtheta;

        //SECOND DERIVATES OF U2
        d2u2aux_dr2 = d2l_dr2 * n;
        d2u2aux_drdteta = dl_dr * dn_dtheta;
        d2u2aux_dteta2 = l * d2n_dtheta2;

        //FIRST DERIVATES: dui_dxj
        double du1aux_dx1, du2aux_dx1, du1aux_dx2, du2aux_dx2;
        du1aux_dx1 = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;
        du1aux_dx2 = sin(teta) * du1aux_dr + (cos(teta) / radius) * du1aux_dteta;

        du2aux_dx1 = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;
        du2aux_dx2 = sin(teta) * du2aux_dr + (cos(teta) / radius) * du2aux_dteta;

        //SECOND DERIVATES: d2ui_dxjdxk
        double d2u1aux_dx1dx1, d2u1aux_dx1dx2, d2u1aux_dx2dx2, d2u2aux_dx1dx1, d2u2aux_dx1dx2, d2u2aux_dx2dx2;

        d2u1aux_dx1dx1 = cos(teta) * cos(teta) * d2u1aux_dr2 + 2.0 * sin(teta) * cos(teta) / (radius * radius) * du1aux_dteta - 2.0 * sin(teta) * cos(teta) / radius * d2u1aux_drdteta + sin(teta) * sin(teta) * (du1aux_dr / radius + d2u1aux_dteta2 / (radius * radius));
        d2u1aux_dx1dx2 = sin(teta) * (cos(teta) * d2u1aux_dr2 + sin(teta) / (radius * radius) * du1aux_dteta - sin(teta) / radius * d2u1aux_drdteta) + cos(teta) / radius * (-sin(teta) * du1aux_dr + cos(teta) * d2u1aux_drdteta - cos(teta) / radius * du1aux_dteta - sin(teta) / radius * d2u1aux_dteta2);
        d2u1aux_dx2dx2 = sin(teta) * sin(teta) * d2u1aux_dr2 - 2.0 * sin(teta) * cos(teta) / (radius * radius) * du1aux_dteta + 2.0 * sin(teta) * cos(teta) / radius * d2u1aux_drdteta + cos(teta) * cos(teta) * (du1aux_dr / radius + d2u1aux_dteta2 / (radius * radius));

        d2u2aux_dx1dx1 = cos(teta) * cos(teta) * d2u2aux_dr2 + 2.0 * sin(teta) * cos(teta) / (radius * radius) * du2aux_dteta - 2.0 * sin(teta) * cos(teta) / radius * d2u2aux_drdteta + sin(teta) * sin(teta) * (du2aux_dr / radius + d2u2aux_dteta2 / (radius * radius));
        d2u2aux_dx1dx2 = sin(teta) * (cos(teta) * d2u2aux_dr2 + sin(teta) / (radius * radius) * du2aux_dteta - sin(teta) / radius * d2u2aux_drdteta) + cos(teta) / radius * (-sin(teta) * du2aux_dr + cos(teta) * d2u2aux_drdteta - cos(teta) / radius * du2aux_dteta - sin(teta) / radius * d2u2aux_dteta2);
        d2u2aux_dx2dx2 = sin(teta) * sin(teta) * d2u2aux_dr2 - 2.0 * sin(teta) * cos(teta) / (radius * radius) * du2aux_dteta + 2.0 * sin(teta) * cos(teta) / radius * d2u2aux_drdteta + cos(teta) * cos(teta) * (du2aux_dr / radius + d2u2aux_dteta2 / (radius * radius));

        bounded_matrix<double, 2, 2> Acaux, Ecaux, Saux, Ptaux;

        Acaux(0, 0) = du1aux_dx1 + 1.0;
        Acaux(0, 1) = du1aux_dx2;
        Acaux(1, 0) = du2aux_dx1;
        Acaux(1, 1) = du2aux_dx2 + 1.0;

        Ecaux = 0.5 * (prod(trans(Acaux), Acaux) - I);

        if (ep == "EPD")
        {
            Saux(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }
        else
        {
            Saux(0, 0) = (young / (1.0 - poisson * poisson)) * (Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / (1.0 - poisson * poisson)) * (Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }

        Ptaux = prod(Acaux, Saux);

        // bounded_matrix<double, 2, 2> Accoupled, Scoupled, Ptcouped, Ptsoma;

        // Accoupled = Ac + Acaux - I;
        // Scoupled = S + Saux;

        // Ptcouped = prod(Accoupled, Scoupled);
        // Ptsoma = Pt + Ptaux;

        // std::cout <<std::scientific << Ptsoma(0, 0) / Ptcouped(0, 0) - 1.0 << " " << Ptsoma(0, 1) / Ptcouped(0, 1) - 1.0 << " " << Ptsoma(1, 0) / Ptcouped(1, 0) - 1.0 << " " << Ptsoma(1, 1) / Ptcouped(1, 1) - 1.0 << " " << std::endl;

        bounded_vector<double, 2> divPtaux, aaux;
        bounded_matrix<double, 2, 2> dAcaux_dx1, dAcaux_dx2, dEcaux_dx1, dEcaux_dx2, grad_aaux, dSaux_dx1, dSaux_dx2, dPtaux_dx1, dPtaux_dx2;

        aaux(0) = du1aux_dx1;
        aaux(1) = du2aux_dx1;

        dAcaux_dx1(0, 0) = d2u1aux_dx1dx1;
        dAcaux_dx1(0, 1) = d2u1aux_dx1dx2;
        dAcaux_dx1(1, 0) = d2u2aux_dx1dx1;
        dAcaux_dx1(1, 1) = d2u2aux_dx1dx2;

        dAcaux_dx2(0, 0) = d2u1aux_dx1dx2;
        dAcaux_dx2(0, 1) = d2u1aux_dx2dx2;
        dAcaux_dx2(1, 0) = d2u2aux_dx1dx2;
        dAcaux_dx2(1, 1) = d2u2aux_dx2dx2;

        dEcaux_dx1 = 0.5 * (prod(trans(dAcaux_dx1), Acaux) + prod(trans(Acaux), dAcaux_dx1));
        dEcaux_dx2 = 0.5 * (prod(trans(dAcaux_dx2), Acaux) + prod(trans(Acaux), dAcaux_dx2));

        if (ep == "EPD")
        {
            dSaux_dx1(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx1(0, 0) + poisson * dEcaux_dx1(1, 1)));
            dSaux_dx1(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx1(1, 1) + poisson * dEcaux_dx1(0, 0)));
            dSaux_dx1(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx1(1, 0);
            dSaux_dx1(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx1(0, 1);

            dSaux_dx2(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx2(0, 0) + poisson * dEcaux_dx2(1, 1)));
            dSaux_dx2(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx2(1, 1) + poisson * dEcaux_dx2(0, 0)));
            dSaux_dx2(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx2(1, 0);
            dSaux_dx2(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx2(0, 1);
        }
        else
        {
            dSaux_dx1(0, 0) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx1(0, 0) + poisson * dEcaux_dx1(1, 1));
            dSaux_dx1(1, 1) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx1(1, 1) + poisson * dEcaux_dx1(0, 0));
            dSaux_dx1(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx1(1, 0);
            dSaux_dx1(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx1(0, 1);

            dSaux_dx2(0, 0) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx2(0, 0) + poisson * dEcaux_dx2(1, 1));
            dSaux_dx2(1, 1) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx2(1, 1) + poisson * dEcaux_dx2(0, 0));
            dSaux_dx2(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx2(1, 0);
            dSaux_dx2(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx2(0, 1);
        }

        dPtaux_dx1 = prod(dAcaux_dx1, Saux) + prod(Acaux, dSaux_dx1);
        dPtaux_dx2 = prod(dAcaux_dx2, Saux) + prod(Acaux, dSaux_dx2);

        divPtaux(0) = dPtaux_dx1(0, 0) + dPtaux_dx2(0, 1);
        divPtaux(1) = dPtaux_dx1(1, 0) + dPtaux_dx2(1, 1);

        grad_aaux(0, 0) = dAcaux_dx1(0, 0);
        grad_aaux(1, 0) = dAcaux_dx1(1, 0);
        grad_aaux(0, 1) = dAcaux_dx2(0, 0);
        grad_aaux(1, 1) = dAcaux_dx2(1, 0);

        bounded_vector<double, 2> term1, term2, term3, unitX1, accel;

        unitX1(0) = 1.0;
        unitX1(1) = 0.0;

        double term4, term5, term6, term7, term8, term9;

        term1 = prod(trans(Pt), aaux);

        term2 = prod(trans(Ptaux), a);

        term3 = -0.5 * (tensorContraction(S, Ecaux) + tensorContraction(Saux, Ec)) * unitX1;

        if (type == "STATIC")
        {
            term4 = 0.0;
        }
        else
        {
            accel = 0.0;
            for (int i = 0; i < connection_.size(); i++)
            {
                accel += phi(i) * connection_[i]->getCurrentAcceleration();
            }
            term4 = density * (inner_prod(accel, aaux));
        }

        term5 = inner_prod(divPtaux, a);

        term6 = tensorContraction(Pt, grad_aaux);

        term7 = tensorContraction(Ptaux, grad_a);

        term8 = -1.0 * tensorContraction(S, dEcaux_dx1);

        term9 = -1.0 * tensorContraction(Saux, dEc_dx1);

        // std::cout << 0.5 * youngLine * term4 << " " << 0.5 * youngLine * term5 << " " << 0.5 * youngLine * term6 << " " << 0.5 * youngLine * term7 << " " << 0.5 * youngLine * term8 << " " << 0.5 * youngLine * term9 << std::endl;

        // term5 = 0.0;

        // term4 = 0.0;

        // term6 = 0.0;

        // term7 = 0.0;

        // term8 = 0.0;

        // term9 = 0.0;

        //q = 1.0;

        jIntegral(0) += 0.5 * youngLine * (inner_prod((term1 + term2 + term3), grad_q) + (term4 + term5 + term6 + term7 + term8 + term9) * q) * weight * j0 * thickness;

        //////////////////////////////////
        //AUXILIARY STATE K1 = 0 E K2 = 1
        //////////////////////////////////

        K1 = 0.0;
        K2 = 1.0;

        l = (sqrt(radius / (2.0 * pi)) / (2.0 * mi));
        dl_dr = sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius);
        d2l_dr2 = -sqrt(2.0 * pi * radius) / (16.0 * mi * pi * radius * radius);

        m = K1 * cos(meio_teta) * (k - cos(teta)) + K2 * sin(meio_teta) * (k + 2.0 + cos(teta));
        dm_dtheta = K1 * (-0.5 * sin(meio_teta) * (k - cos(teta)) + sin(teta) * cos(meio_teta)) + K2 * (0.5 * cos(meio_teta) * (k + 2.0 + cos(teta)) - sin(teta) * sin(meio_teta));
        d2m_dtetha2 = K1 * (-cos(meio_teta) * 0.25 * (k - cos(teta)) - sin(meio_teta) * sin(teta) + cos(meio_teta) * cos(teta)) + K2 * (-sin(meio_teta) * 0.25 * (k + 2.0 + cos(teta)) - sin(teta) * cos(meio_teta) - cos(teta) * sin(meio_teta));

        n = K1 * sin(meio_teta) * (k - cos(teta)) - K2 * cos(meio_teta) * (k - 2.0 + cos(teta));
        dn_dtheta = K1 * (0.5 * cos(meio_teta) * (k - cos(teta)) + sin(teta) * sin(meio_teta)) + K2 * (0.5 * sin(meio_teta) * (k - 2.0 + cos(teta)) + sin(teta) * cos(meio_teta));
        d2n_dtheta2 = K1 * (-sin(meio_teta) * 0.25 * (k - cos(teta)) + cos(meio_teta) * sin(teta) + sin(meio_teta) * cos(teta)) - K2 * (-cos(meio_teta) * 0.25 * (k - 2.0 + cos(teta)) + sin(meio_teta) * sin(teta) - cos(teta) * cos(meio_teta));

        //FIRST DERIVATES OF U1
        du1aux_dr = dl_dr * m;
        du1aux_dteta = l * dm_dtheta;

        //SECOND DERIVATES OF U1
        d2u1aux_dr2 = d2l_dr2 * m;
        d2u1aux_drdteta = dl_dr * dm_dtheta;
        d2u1aux_dteta2 = l * d2m_dtetha2;

        //FIRST DERIVATES OF U2
        du2aux_dr = dl_dr * n;
        du2aux_dteta = l * dn_dtheta;

        //SECOND DERIVATES OF U2
        d2u2aux_dr2 = d2l_dr2 * n;
        d2u2aux_drdteta = dl_dr * dn_dtheta;
        d2u2aux_dteta2 = l * d2n_dtheta2;

        //FIRST DERIVATES: dui_dxj
        du1aux_dx1 = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;
        du1aux_dx2 = sin(teta) * du1aux_dr + (cos(teta) / radius) * du1aux_dteta;

        du2aux_dx1 = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;
        du2aux_dx2 = sin(teta) * du2aux_dr + (cos(teta) / radius) * du2aux_dteta;

        //SECOND DERIVATES: d2ui_dxjdxk
        d2u1aux_dx1dx1 = cos(teta) * cos(teta) * d2u1aux_dr2 + 2.0 * sin(teta) * cos(teta) / (radius * radius) * du1aux_dteta - 2.0 * sin(teta) * cos(teta) / radius * d2u1aux_drdteta + sin(teta) * sin(teta) * (du1aux_dr / radius + d2u1aux_dteta2 / (radius * radius));
        d2u1aux_dx1dx2 = sin(teta) * (cos(teta) * d2u1aux_dr2 + sin(teta) / (radius * radius) * du1aux_dteta - sin(teta) / radius * d2u1aux_drdteta) + cos(teta) / radius * (-sin(teta) * du1aux_dr + cos(teta) * d2u1aux_drdteta - cos(teta) / radius * du1aux_dteta - sin(teta) / radius * d2u1aux_dteta2);
        d2u1aux_dx2dx2 = sin(teta) * sin(teta) * d2u1aux_dr2 - 2.0 * sin(teta) * cos(teta) / (radius * radius) * du1aux_dteta + 2.0 * sin(teta) * cos(teta) / radius * d2u1aux_drdteta + cos(teta) * cos(teta) * (du1aux_dr / radius + d2u1aux_dteta2 / (radius * radius));

        d2u2aux_dx1dx1 = cos(teta) * cos(teta) * d2u2aux_dr2 + 2.0 * sin(teta) * cos(teta) / (radius * radius) * du2aux_dteta - 2.0 * sin(teta) * cos(teta) / radius * d2u2aux_drdteta + sin(teta) * sin(teta) * (du2aux_dr / radius + d2u2aux_dteta2 / (radius * radius));
        d2u2aux_dx1dx2 = sin(teta) * (cos(teta) * d2u2aux_dr2 + sin(teta) / (radius * radius) * du2aux_dteta - sin(teta) / radius * d2u2aux_drdteta) + cos(teta) / radius * (-sin(teta) * du2aux_dr + cos(teta) * d2u2aux_drdteta - cos(teta) / radius * du2aux_dteta - sin(teta) / radius * d2u2aux_dteta2);
        d2u2aux_dx2dx2 = sin(teta) * sin(teta) * d2u2aux_dr2 - 2.0 * sin(teta) * cos(teta) / (radius * radius) * du2aux_dteta + 2.0 * sin(teta) * cos(teta) / radius * d2u2aux_drdteta + cos(teta) * cos(teta) * (du2aux_dr / radius + d2u2aux_dteta2 / (radius * radius));

        Acaux(0, 0) = du1aux_dx1 + 1.0;
        Acaux(0, 1) = du1aux_dx2;
        Acaux(1, 0) = du2aux_dx1;
        Acaux(1, 1) = du2aux_dx2 + 1.0;

        Ecaux = 0.5 * (prod(trans(Acaux), Acaux) - I);

        if (ep == "EPD")
        {
            Saux(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }
        else
        {
            Saux(0, 0) = (young / (1.0 - poisson * poisson)) * (Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / (1.0 - poisson * poisson)) * (Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }

        Ptaux = prod(Acaux, Saux);

        aaux(0) = du1aux_dx1;
        aaux(1) = du2aux_dx1;

        dAcaux_dx1(0, 0) = d2u1aux_dx1dx1;
        dAcaux_dx1(0, 1) = d2u1aux_dx1dx2;
        dAcaux_dx1(1, 0) = d2u2aux_dx1dx1;
        dAcaux_dx1(1, 1) = d2u2aux_dx1dx2;

        dAcaux_dx2(0, 0) = d2u1aux_dx1dx2;
        dAcaux_dx2(0, 1) = d2u1aux_dx2dx2;
        dAcaux_dx2(1, 0) = d2u2aux_dx1dx2;
        dAcaux_dx2(1, 1) = d2u2aux_dx2dx2;

        dEcaux_dx1 = 0.5 * (prod(trans(dAcaux_dx1), Acaux) + prod(trans(Acaux), dAcaux_dx1));
        dEcaux_dx2 = 0.5 * (prod(trans(dAcaux_dx2), Acaux) + prod(trans(Acaux), dAcaux_dx2));

        if (ep == "EPD")
        {
            dSaux_dx1(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx1(0, 0) + poisson * dEcaux_dx1(1, 1)));
            dSaux_dx1(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx1(1, 1) + poisson * dEcaux_dx1(0, 0)));
            dSaux_dx1(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx1(1, 0);
            dSaux_dx1(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx1(0, 1);

            dSaux_dx2(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx2(0, 0) + poisson * dEcaux_dx2(1, 1)));
            dSaux_dx2(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx2(1, 1) + poisson * dEcaux_dx2(0, 0)));
            dSaux_dx2(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx2(1, 0);
            dSaux_dx2(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx2(0, 1);
        }
        else
        {
            dSaux_dx1(0, 0) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx1(0, 0) + poisson * dEcaux_dx1(1, 1));
            dSaux_dx1(1, 1) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx1(1, 1) + poisson * dEcaux_dx1(0, 0));
            dSaux_dx1(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx1(1, 0);
            dSaux_dx1(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx1(0, 1);

            dSaux_dx2(0, 0) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx2(0, 0) + poisson * dEcaux_dx2(1, 1));
            dSaux_dx2(1, 1) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx2(1, 1) + poisson * dEcaux_dx2(0, 0));
            dSaux_dx2(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx2(1, 0);
            dSaux_dx2(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx2(0, 1);
        }

        dPtaux_dx1 = prod(dAcaux_dx1, Saux) + prod(Acaux, dSaux_dx1);
        dPtaux_dx2 = prod(dAcaux_dx2, Saux) + prod(Acaux, dSaux_dx2);

        divPtaux(0) = dPtaux_dx1(0, 0) + dPtaux_dx2(0, 1);
        divPtaux(1) = dPtaux_dx1(1, 0) + dPtaux_dx2(1, 1);

        grad_aaux(0, 0) = dAcaux_dx1(0, 0);
        grad_aaux(1, 0) = dAcaux_dx1(1, 0);
        grad_aaux(0, 1) = dAcaux_dx2(0, 0);
        grad_aaux(1, 1) = dAcaux_dx2(1, 0);

        unitX1(0) = 1.0;
        unitX1(1) = 0.0;

        term1 = prod(Pt, aaux);

        term2 = prod(Ptaux, a);

        term3 = -1.0 * tensorContraction(S, Ecaux) * unitX1;

        if (type == "STATIC")
        {
            term4 = 0.0;
        }
        else
        {
            accel = 0.0;
            for (int i = 0; i < connection_.size(); i++)
            {
                accel += phi(i) * connection_[i]->getCurrentAcceleration();
            }
            term4 = density * (inner_prod(accel, aaux));
        }

        term5 = inner_prod(divPtaux, a);

        term6 = tensorContraction(Pt, grad_aaux);

        term7 = tensorContraction(Ptaux, grad_a);

        term8 = -1.0 * tensorContraction(S, dEcaux_dx1);

        term9 = -1.0 * tensorContraction(Saux, dEc_dx1);

        term5 = 0.0;
        term4 = 0.0;

        term6 = 0.0;

        term7 = 0.0;

        term8 = 0.0;

        term9 = 0.0;

        jIntegral(1) += 0.5 * youngLine * (inner_prod((term1 + term2 + term3), grad_q) + (term4 + term5 + term6 + term7 + term8 + term9) * q) * weight * j0 * thickness;
    }
    return jIntegral;
}

bounded_vector<double, 2> Element::domainJ_IntegralInitial2(const int &nHammer, const std::string &ep, const std::string &type, const double &Jintegral_radius, const double &rotation, const bounded_vector<double, 2> &tipCoordinates)
{
    bounded_vector<double, 2> jIntegral;
    jIntegral(0) = 0.0;
    jIntegral(1) = 0.0;

    int nnos = connection_.size();

    //Material properties
    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);
    double thickness = mesh_->getThickness();
    double mi = young / (2.0 * (1.0 + poisson));
    double youngLine;
    double k;
    if (ep == "EPD")
    {
        k = 3.0 - 4.0 * poisson;
        youngLine = young / (1.0 - poisson * poisson);
    }
    else
    {
        k = (3.0 - poisson) / (1.0 + poisson);
        youngLine = young;
    }

    matrix<double> domainIntegrationPoints = hammerQuadrature(nHammer);

    bounded_matrix<double, 2, 2> matrixRotation;
    matrixRotation(0, 0) = cos(rotation);
    matrixRotation(0, 1) = sin(rotation);
    matrixRotation(1, 0) = -sin(rotation);
    matrixRotation(1, 1) = cos(rotation);

    for (int ihh = 0; ihh < nHammer; ihh++)
    {

        double xsi1 = domainIntegrationPoints(ihh, 0);
        double xsi2 = domainIntegrationPoints(ihh, 1);
        double weight = domainIntegrationPoints(ihh, 2);

        vector<double> phi = domainShapeFunction(xsi1, xsi2);
        matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node
        matrix<double> d2phi_d2xsi = domainSecondDerivativeShapeFunction(xsi1, xsi2);

        bounded_matrix<double, 2, 2> dA1_dxsi1, dA1_dxsi2, A0loc, mataux; //A1loc,
        dA1_dxsi1(0, 0) = 0.0;
        dA1_dxsi1(0, 1) = 0.0;
        dA1_dxsi1(1, 0) = 0.0;
        dA1_dxsi1(1, 1) = 0.0;

        dA1_dxsi2(0, 0) = 0.0;
        dA1_dxsi2(0, 1) = 0.0;
        dA1_dxsi2(1, 0) = 0.0;
        dA1_dxsi2(1, 1) = 0.0;

        A0loc(0, 0) = 0.0;
        A0loc(0, 1) = 0.0;
        A0loc(1, 0) = 0.0;
        A0loc(1, 1) = 0.0;

        // A1loc(0, 0) = 0.0;
        // A1loc(0, 1) = 0.0;
        // A1loc(1, 0) = 0.0;
        // A1loc(1, 1) = 0.0;

        for (int in = 0; in < nnos; in++)
        {
            bounded_vector<double, 2> coordNode = prod(matrixRotation, connection_[in]->getCurrentCoordinate());

            bounded_vector<double, 2> initialNode = prod(matrixRotation, connection_[in]->getInitialCoordinate());

            dA1_dxsi1(0, 0) += d2phi_d2xsi(0, in) * coordNode(0);
            dA1_dxsi1(0, 1) += d2phi_d2xsi(1, in) * coordNode(0);
            dA1_dxsi1(1, 0) += d2phi_d2xsi(0, in) * coordNode(1);
            dA1_dxsi1(1, 1) += d2phi_d2xsi(1, in) * coordNode(1);

            dA1_dxsi2(0, 0) += d2phi_d2xsi(1, in) * coordNode(0);
            dA1_dxsi2(0, 1) += d2phi_d2xsi(2, in) * coordNode(0);
            dA1_dxsi2(1, 0) += d2phi_d2xsi(1, in) * coordNode(1);
            dA1_dxsi2(1, 1) += d2phi_d2xsi(2, in) * coordNode(1);

            A0loc(0, 0) += initialNode(0) * dphi_dxsi(0, in);
            A0loc(0, 1) += initialNode(0) * dphi_dxsi(1, in);
            A0loc(1, 0) += initialNode(1) * dphi_dxsi(0, in);
            A0loc(1, 1) += initialNode(1) * dphi_dxsi(1, in);

            // A1loc(0, 0) += coordNode(0) * dphi_dxsi(0, in);
            // A1loc(0, 1) += coordNode(0) * dphi_dxsi(1, in);
            // A1loc(1, 0) += coordNode(1) * dphi_dxsi(0, in);
            // A1loc(1, 1) += coordNode(1) * dphi_dxsi(1, in);
        }

        bounded_matrix<double, 2, 2> A0glob = referenceJacobianMatrix(xsi1, xsi2);
        double j0 = jacobianDeterminant(A0glob);
        bounded_matrix<double, 2, 2> A0Iglob = inverseMatrix(A0glob);
        bounded_matrix<double, 2, 2> A1glob = currentJacobianMatrix(xsi1, xsi2);
        bounded_matrix<double, 2, 2> Acglob = prod(A1glob, A0Iglob); //dyi / dxj

        bounded_matrix<double, 2, 2> Ac; // Aclocteste;

        //Aclocteste = prod(A1loc, inverseMatrix(A0loc));

        mataux = prod(matrixRotation, Acglob);
        Ac = prod(mataux, trans(matrixRotation)); //gradiente nos eixos da ponta da fissura

        // std::cout << Aclocteste(0, 0) / Ac(0, 0) << " " << Aclocteste(0, 1) / Ac(0, 1) << " " << Aclocteste(1, 0) / Ac(1, 0) << " " << Aclocteste(1, 1) / Ac(1, 1) << std::endl;

        const identity_matrix<double> I(2);
        bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);
        bounded_matrix<double, 2, 2> S;

        if (ep == "EPD")
        {
            S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }
        else
        {
            S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
            S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
            S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
            S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
        }

        bounded_matrix<double, 2, 2> Pt;
        Pt = prod(Ac, S);

        bounded_matrix<double, 2, 2> dAc_dxsi1, dAc_dxsi2, A0Iloc;

        A0Iloc = inverseMatrix(A0loc);
        dAc_dxsi1 = prod(dA1_dxsi1, A0Iloc);
        dAc_dxsi2 = prod(dA1_dxsi2, A0Iloc);

        bounded_vector<double, 2> a;
        a(0) = Ac(0, 0) - 1.0;
        a(1) = Ac(1, 0);

        bounded_matrix<double, 2, 2> grad_a, dEc_dx1, dAc_dx1;

        grad_a(0, 0) = dAc_dxsi1(0, 0) * A0Iloc(0, 0) + dAc_dxsi2(0, 0) * A0Iloc(1, 0); //ok
        grad_a(0, 1) = dAc_dxsi1(0, 0) * A0Iloc(0, 1) + dAc_dxsi2(0, 0) * A0Iloc(1, 1);
        grad_a(1, 0) = dAc_dxsi1(1, 0) * A0Iloc(0, 0) + dAc_dxsi2(1, 0) * A0Iloc(1, 0); //ok
        grad_a(1, 1) = dAc_dxsi1(1, 0) * A0Iloc(0, 1) + dAc_dxsi2(1, 0) * A0Iloc(1, 1);

        dAc_dx1 = dAc_dxsi1 * A0Iloc(0, 0) + dAc_dxsi2 * A0Iloc(1, 0); //é igual a grad_a

        dEc_dx1 = 0.5 * (prod(trans(dAc_dx1), Ac) + prod(trans(Ac), dAc_dx1));

        //////////////////////////////////
        //Q FUNCTION
        //////////////////////////////////
        double radius = 0.0;
        bounded_vector<double, 2> dradius_dxsi, coordP;
        coordP(0) = 0.0;
        coordP(1) = 0.0;
        dradius_dxsi(0) = 0.0;
        dradius_dxsi(1) = 0.0;

        for (int in = 0; in < connection_.size(); in++)
        {
            coordP += phi(in) * connection_[in]->getInitialCoordinate();

            const double radiusNODE = norm_2(connection_[in]->getInitialCoordinate() - tipCoordinates);

            dradius_dxsi(0) += dphi_dxsi(0, in) * radiusNODE; //dDAhammer_dxsiLocal1
            dradius_dxsi(1) += dphi_dxsi(1, in) * radiusNODE; //dDAhammer_dxsiLocal1

            radius += phi(in) * radiusNODE;
        }

        bounded_vector<double, 2> tipToPointAux = coordP - tipCoordinates;
        bounded_vector<double, 2> tipToPoint = prod(matrixRotation, tipToPointAux);
        const double teta = atan2(tipToPoint(1), tipToPoint(0));

        const double aux = radius / Jintegral_radius;
        const double q = 2.0 * aux * aux * aux - 3.0 * aux * aux + 1.0;
        const double dq_dradius = (6.0 / Jintegral_radius) * aux * aux - (6.0 / Jintegral_radius) * aux;

        // const double q = 1.0 - radius / Jintegral_radius;
        // const double dq_dradius = -1.0 / Jintegral_radius;

        bounded_vector<double, 2> dq_dxsi;
        dq_dxsi(0) = dq_dradius * dradius_dxsi(0);
        dq_dxsi(1) = dq_dradius * dradius_dxsi(1);

        bounded_vector<double, 2> grad_q;

        grad_q(0) = dq_dxsi(0) * A0Iloc(0, 0) + dq_dxsi(1) * A0Iloc(1, 0);
        grad_q(1) = dq_dxsi(0) * A0Iloc(0, 1) + dq_dxsi(1) * A0Iloc(1, 1);

        //////////////////////////////////
        //AUXILIARY STATES K1 = 1 E K2 = 0;
        //////////////////////////////////

        const double pi = 3.14159265358979323846;
        const double meio_teta = 0.5 * teta;

        double dl_dr, d2l_dr2, dm_dtheta, d2m_dtetha2, dn_dtheta, d2n_dtheta2, K1, K2, l, m, n;

        K1 = 1.0;
        K2 = 0.0;

        l = (sqrt(radius / (2.0 * pi)) / (2.0 * mi));
        dl_dr = sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius);
        d2l_dr2 = -sqrt(2.0 * pi * radius) / (16.0 * mi * pi * radius * radius);

        m = K1 * cos(meio_teta) * (k - cos(teta)) + K2 * sin(meio_teta) * (k + 2.0 + cos(teta));
        dm_dtheta = K1 * (-0.5 * sin(meio_teta) * (k - cos(teta)) + sin(teta) * cos(meio_teta)) + K2 * (0.5 * cos(meio_teta) * (k + 2.0 + cos(teta)) - sin(teta) * sin(meio_teta));
        d2m_dtetha2 = K1 * (-cos(meio_teta) * 0.25 * (k - cos(teta)) - sin(meio_teta) * sin(teta) + cos(meio_teta) * cos(teta)) + K2 * (-sin(meio_teta) * 0.25 * (k + 2.0 + cos(teta)) - sin(teta) * cos(meio_teta) - cos(teta) * sin(meio_teta));

        n = K1 * sin(meio_teta) * (k - cos(teta)) - K2 * cos(meio_teta) * (k - 2.0 + cos(teta));
        dn_dtheta = K1 * (0.5 * cos(meio_teta) * (k - cos(teta)) + sin(teta) * sin(meio_teta)) + K2 * (0.5 * sin(meio_teta) * (k - 2.0 + cos(teta)) + sin(teta) * cos(meio_teta));
        d2n_dtheta2 = K1 * (-sin(meio_teta) * 0.25 * (k - cos(teta)) + cos(meio_teta) * sin(teta) + sin(meio_teta) * cos(teta)) - K2 * (-cos(meio_teta) * 0.25 * (k - 2.0 + cos(teta)) + sin(meio_teta) * sin(teta) - cos(teta) * cos(meio_teta));

        double du1aux_dr, du2aux_dr, du1aux_dteta, du2aux_dteta, d2u1aux_dr2, d2u1aux_drdteta, d2u1aux_dteta2, d2u2aux_dr2, d2u2aux_drdteta, d2u2aux_dteta2;

        //FIRST DERIVATES OF U1
        du1aux_dr = dl_dr * m;
        du1aux_dteta = l * dm_dtheta;

        //SECOND DERIVATES OF U1
        d2u1aux_dr2 = d2l_dr2 * m;
        d2u1aux_drdteta = dl_dr * dm_dtheta;
        d2u1aux_dteta2 = l * d2m_dtetha2;

        //FIRST DERIVATES OF U2
        du2aux_dr = dl_dr * n;
        du2aux_dteta = l * dn_dtheta;

        //SECOND DERIVATES OF U2
        d2u2aux_dr2 = d2l_dr2 * n;
        d2u2aux_drdteta = dl_dr * dn_dtheta;
        d2u2aux_dteta2 = l * d2n_dtheta2;

        //FIRST DERIVATES: dui_dxj
        double du1aux_dx1, du2aux_dx1, du1aux_dx2, du2aux_dx2;
        du1aux_dx1 = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;
        du1aux_dx2 = sin(teta) * du1aux_dr + (cos(teta) / radius) * du1aux_dteta;

        du2aux_dx1 = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;
        du2aux_dx2 = sin(teta) * du2aux_dr + (cos(teta) / radius) * du2aux_dteta;

        //SECOND DERIVATES: d2ui_dxjdxk
        double d2u1aux_dx1dx1, d2u1aux_dx1dx2, d2u1aux_dx2dx2, d2u2aux_dx1dx1, d2u2aux_dx1dx2, d2u2aux_dx2dx2;

        d2u1aux_dx1dx1 = cos(teta) * cos(teta) * d2u1aux_dr2 + 2.0 * sin(teta) * cos(teta) / (radius * radius) * du1aux_dteta - 2.0 * sin(teta) * cos(teta) / radius * d2u1aux_drdteta + sin(teta) * sin(teta) * (du1aux_dr / radius + d2u1aux_dteta2 / (radius * radius));
        d2u1aux_dx1dx2 = sin(teta) * (cos(teta) * d2u1aux_dr2 + sin(teta) / (radius * radius) * du1aux_dteta - sin(teta) / radius * d2u1aux_drdteta) + cos(teta) / radius * (-sin(teta) * du1aux_dr + cos(teta) * d2u1aux_drdteta - cos(teta) / radius * du1aux_dteta - sin(teta) / radius * d2u1aux_dteta2);
        d2u1aux_dx2dx2 = sin(teta) * sin(teta) * d2u1aux_dr2 - 2.0 * sin(teta) * cos(teta) / (radius * radius) * du1aux_dteta + 2.0 * sin(teta) * cos(teta) / radius * d2u1aux_drdteta + cos(teta) * cos(teta) * (du1aux_dr / radius + d2u1aux_dteta2 / (radius * radius));

        d2u2aux_dx1dx1 = cos(teta) * cos(teta) * d2u2aux_dr2 + 2.0 * sin(teta) * cos(teta) / (radius * radius) * du2aux_dteta - 2.0 * sin(teta) * cos(teta) / radius * d2u2aux_drdteta + sin(teta) * sin(teta) * (du2aux_dr / radius + d2u2aux_dteta2 / (radius * radius));
        d2u2aux_dx1dx2 = sin(teta) * (cos(teta) * d2u2aux_dr2 + sin(teta) / (radius * radius) * du2aux_dteta - sin(teta) / radius * d2u2aux_drdteta) + cos(teta) / radius * (-sin(teta) * du2aux_dr + cos(teta) * d2u2aux_drdteta - cos(teta) / radius * du2aux_dteta - sin(teta) / radius * d2u2aux_dteta2);
        d2u2aux_dx2dx2 = sin(teta) * sin(teta) * d2u2aux_dr2 - 2.0 * sin(teta) * cos(teta) / (radius * radius) * du2aux_dteta + 2.0 * sin(teta) * cos(teta) / radius * d2u2aux_drdteta + cos(teta) * cos(teta) * (du2aux_dr / radius + d2u2aux_dteta2 / (radius * radius));

        bounded_matrix<double, 2, 2> Acaux, Ecaux, Saux, Ptaux;

        Acaux(0, 0) = du1aux_dx1 + 1.0;
        Acaux(0, 1) = du1aux_dx2;
        Acaux(1, 0) = du2aux_dx1;
        Acaux(1, 1) = du2aux_dx2 + 1.0;

        Ecaux = 0.5 * (prod(trans(Acaux), Acaux) - I);

        if (ep == "EPD")
        {
            Saux(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }
        else
        {
            Saux(0, 0) = (young / (1.0 - poisson * poisson)) * (Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / (1.0 - poisson * poisson)) * (Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }

        Ptaux = prod(Acaux, Saux);

        // bounded_matrix<double, 2, 2> Accoupled, Scoupled, Ptcouped, Ptsoma;

        // Accoupled = Ac + Acaux - I;
        // Scoupled = S + Saux;

        // Ptcouped = prod(Accoupled, Scoupled);
        // Ptsoma = Pt + Ptaux;

        // std::cout <<std::scientific << Ptsoma(0, 0) / Ptcouped(0, 0) - 1.0 << " " << Ptsoma(0, 1) / Ptcouped(0, 1) - 1.0 << " " << Ptsoma(1, 0) / Ptcouped(1, 0) - 1.0 << " " << Ptsoma(1, 1) / Ptcouped(1, 1) - 1.0 << " " << std::endl;

        bounded_vector<double, 2> aaux;                                                         //divPtaux
        bounded_matrix<double, 2, 2> dAcaux_dx1, dAcaux_dx2, dEcaux_dx1, dEcaux_dx2, grad_aaux; // dSaux_dx1, dSaux_dx2, dPtaux_dx1, dPtaux_dx2;

        aaux(0) = du1aux_dx1;
        aaux(1) = du2aux_dx1;

        dAcaux_dx1(0, 0) = d2u1aux_dx1dx1;
        dAcaux_dx1(0, 1) = d2u1aux_dx1dx2;
        dAcaux_dx1(1, 0) = d2u2aux_dx1dx1;
        dAcaux_dx1(1, 1) = d2u2aux_dx1dx2;

        dAcaux_dx2(0, 0) = d2u1aux_dx1dx2;
        dAcaux_dx2(0, 1) = d2u1aux_dx2dx2;
        dAcaux_dx2(1, 0) = d2u2aux_dx1dx2;
        dAcaux_dx2(1, 1) = d2u2aux_dx2dx2;

        dEcaux_dx1 = 0.5 * (prod(trans(dAcaux_dx1), Acaux) + prod(trans(Acaux), dAcaux_dx1));
        dEcaux_dx2 = 0.5 * (prod(trans(dAcaux_dx2), Acaux) + prod(trans(Acaux), dAcaux_dx2));

        // if (ep == "EPD")
        // {
        //     dSaux_dx1(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx1(0, 0) + poisson * dEcaux_dx1(1, 1)));
        //     dSaux_dx1(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx1(1, 1) + poisson * dEcaux_dx1(0, 0)));
        //     dSaux_dx1(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx1(1, 0);
        //     dSaux_dx1(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx1(0, 1);

        //     dSaux_dx2(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx2(0, 0) + poisson * dEcaux_dx2(1, 1)));
        //     dSaux_dx2(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx2(1, 1) + poisson * dEcaux_dx2(0, 0)));
        //     dSaux_dx2(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx2(1, 0);
        //     dSaux_dx2(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx2(0, 1);
        // }
        // else
        // {
        //     dSaux_dx1(0, 0) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx1(0, 0) + poisson * dEcaux_dx1(1, 1));
        //     dSaux_dx1(1, 1) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx1(1, 1) + poisson * dEcaux_dx1(0, 0));
        //     dSaux_dx1(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx1(1, 0);
        //     dSaux_dx1(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx1(0, 1);

        //     dSaux_dx2(0, 0) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx2(0, 0) + poisson * dEcaux_dx2(1, 1));
        //     dSaux_dx2(1, 1) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx2(1, 1) + poisson * dEcaux_dx2(0, 0));
        //     dSaux_dx2(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx2(1, 0);
        //     dSaux_dx2(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx2(0, 1);
        // }

        // dPtaux_dx1 = prod(dAcaux_dx1, Saux) + prod(Acaux, dSaux_dx1);
        // dPtaux_dx2 = prod(dAcaux_dx2, Saux) + prod(Acaux, dSaux_dx2);

        // divPtaux(0) = dPtaux_dx1(0, 0) + dPtaux_dx2(0, 1);
        // divPtaux(1) = dPtaux_dx1(1, 0) + dPtaux_dx2(1, 1);

        grad_aaux(0, 0) = dAcaux_dx1(0, 0);
        grad_aaux(1, 0) = dAcaux_dx1(1, 0);
        grad_aaux(0, 1) = dAcaux_dx2(0, 0);
        grad_aaux(1, 1) = dAcaux_dx2(1, 0);

        bounded_vector<double, 2> term1, term2, term3, unitX1, accel;

        unitX1(0) = 1.0;
        unitX1(1) = 0.0;

        double term4, term5, term6, term7, term8;

        term1 = prod(trans(Pt), aaux);

        term2 = prod(trans(Ptaux), a);

        term3 = -0.5 * (tensorContraction(S, Ecaux) + tensorContraction(Saux, Ec)) * unitX1;

        if (type == "STATIC")
        {
            term4 = 0.0;
        }
        else
        {
            accel(0) = 0.0;
            accel(1) = 0.0;
            for (int ij = 0; ij < connection_.size(); ij++)
            {
                accel += phi(ij) * prod(matrixRotation, connection_[ij]->getCurrentAcceleration());
            }
            term4 = density * (inner_prod(accel, aaux));
        }

        term5 = tensorContraction(Pt, grad_aaux);

        term6 = tensorContraction(Ptaux, grad_a);

        term7 = -1.0 * tensorContraction(S, dEcaux_dx1);

        term8 = -1.0 * tensorContraction(Saux, dEc_dx1);

        jIntegral(0) += 0.5 * youngLine * (inner_prod((term1 + term2 + term3), grad_q) + (term4 + term5 + term6 + term7 + term8) * q) * weight * j0 * thickness;

        //////////////////////////////////
        //AUXILIARY STATE K1 = 0 E K2 = 1
        //////////////////////////////////

        K1 = 0.0;
        K2 = 1.0;

        l = (sqrt(radius / (2.0 * pi)) / (2.0 * mi));
        dl_dr = sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius);
        d2l_dr2 = -sqrt(2.0 * pi * radius) / (16.0 * mi * pi * radius * radius);

        m = K1 * cos(meio_teta) * (k - cos(teta)) + K2 * sin(meio_teta) * (k + 2.0 + cos(teta));
        dm_dtheta = K1 * (-0.5 * sin(meio_teta) * (k - cos(teta)) + sin(teta) * cos(meio_teta)) + K2 * (0.5 * cos(meio_teta) * (k + 2.0 + cos(teta)) - sin(teta) * sin(meio_teta));
        d2m_dtetha2 = K1 * (-cos(meio_teta) * 0.25 * (k - cos(teta)) - sin(meio_teta) * sin(teta) + cos(meio_teta) * cos(teta)) + K2 * (-sin(meio_teta) * 0.25 * (k + 2.0 + cos(teta)) - sin(teta) * cos(meio_teta) - cos(teta) * sin(meio_teta));

        n = K1 * sin(meio_teta) * (k - cos(teta)) - K2 * cos(meio_teta) * (k - 2.0 + cos(teta));
        dn_dtheta = K1 * (0.5 * cos(meio_teta) * (k - cos(teta)) + sin(teta) * sin(meio_teta)) + K2 * (0.5 * sin(meio_teta) * (k - 2.0 + cos(teta)) + sin(teta) * cos(meio_teta));
        d2n_dtheta2 = K1 * (-sin(meio_teta) * 0.25 * (k - cos(teta)) + cos(meio_teta) * sin(teta) + sin(meio_teta) * cos(teta)) - K2 * (-cos(meio_teta) * 0.25 * (k - 2.0 + cos(teta)) + sin(meio_teta) * sin(teta) - cos(teta) * cos(meio_teta));

        //FIRST DERIVATES OF U1
        du1aux_dr = dl_dr * m;
        du1aux_dteta = l * dm_dtheta;

        //SECOND DERIVATES OF U1
        d2u1aux_dr2 = d2l_dr2 * m;
        d2u1aux_drdteta = dl_dr * dm_dtheta;
        d2u1aux_dteta2 = l * d2m_dtetha2;

        //FIRST DERIVATES OF U2
        du2aux_dr = dl_dr * n;
        du2aux_dteta = l * dn_dtheta;

        //SECOND DERIVATES OF U2
        d2u2aux_dr2 = d2l_dr2 * n;
        d2u2aux_drdteta = dl_dr * dn_dtheta;
        d2u2aux_dteta2 = l * d2n_dtheta2;

        //FIRST DERIVATES: dui_dxj
        du1aux_dx1 = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;
        du1aux_dx2 = sin(teta) * du1aux_dr + (cos(teta) / radius) * du1aux_dteta;

        du2aux_dx1 = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;
        du2aux_dx2 = sin(teta) * du2aux_dr + (cos(teta) / radius) * du2aux_dteta;

        //SECOND DERIVATES: d2ui_dxjdxk
        d2u1aux_dx1dx1 = cos(teta) * cos(teta) * d2u1aux_dr2 + 2.0 * sin(teta) * cos(teta) / (radius * radius) * du1aux_dteta - 2.0 * sin(teta) * cos(teta) / radius * d2u1aux_drdteta + sin(teta) * sin(teta) * (du1aux_dr / radius + d2u1aux_dteta2 / (radius * radius));
        d2u1aux_dx1dx2 = sin(teta) * (cos(teta) * d2u1aux_dr2 + sin(teta) / (radius * radius) * du1aux_dteta - sin(teta) / radius * d2u1aux_drdteta) + cos(teta) / radius * (-sin(teta) * du1aux_dr + cos(teta) * d2u1aux_drdteta - cos(teta) / radius * du1aux_dteta - sin(teta) / radius * d2u1aux_dteta2);
        d2u1aux_dx2dx2 = sin(teta) * sin(teta) * d2u1aux_dr2 - 2.0 * sin(teta) * cos(teta) / (radius * radius) * du1aux_dteta + 2.0 * sin(teta) * cos(teta) / radius * d2u1aux_drdteta + cos(teta) * cos(teta) * (du1aux_dr / radius + d2u1aux_dteta2 / (radius * radius));

        d2u2aux_dx1dx1 = cos(teta) * cos(teta) * d2u2aux_dr2 + 2.0 * sin(teta) * cos(teta) / (radius * radius) * du2aux_dteta - 2.0 * sin(teta) * cos(teta) / radius * d2u2aux_drdteta + sin(teta) * sin(teta) * (du2aux_dr / radius + d2u2aux_dteta2 / (radius * radius));
        d2u2aux_dx1dx2 = sin(teta) * (cos(teta) * d2u2aux_dr2 + sin(teta) / (radius * radius) * du2aux_dteta - sin(teta) / radius * d2u2aux_drdteta) + cos(teta) / radius * (-sin(teta) * du2aux_dr + cos(teta) * d2u2aux_drdteta - cos(teta) / radius * du2aux_dteta - sin(teta) / radius * d2u2aux_dteta2);
        d2u2aux_dx2dx2 = sin(teta) * sin(teta) * d2u2aux_dr2 - 2.0 * sin(teta) * cos(teta) / (radius * radius) * du2aux_dteta + 2.0 * sin(teta) * cos(teta) / radius * d2u2aux_drdteta + cos(teta) * cos(teta) * (du2aux_dr / radius + d2u2aux_dteta2 / (radius * radius));

        Acaux(0, 0) = du1aux_dx1 + 1.0;
        Acaux(0, 1) = du1aux_dx2;
        Acaux(1, 0) = du2aux_dx1;
        Acaux(1, 1) = du2aux_dx2 + 1.0;

        Ecaux = 0.5 * (prod(trans(Acaux), Acaux) - I);

        if (ep == "EPD")
        {
            Saux(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }
        else
        {
            Saux(0, 0) = (young / (1.0 - poisson * poisson)) * (Ecaux(0, 0) + poisson * Ecaux(1, 1));
            Saux(1, 1) = (young / (1.0 - poisson * poisson)) * (Ecaux(1, 1) + poisson * Ecaux(0, 0));
            Saux(1, 0) = (young / (1.0 + poisson)) * Ecaux(1, 0);
            Saux(0, 1) = (young / (1.0 + poisson)) * Ecaux(0, 1);
        }

        Ptaux = prod(Acaux, Saux);

        // bounded_matrix<double, 2, 2> Accoupled, Scoupled, Ptcouped, Ptsoma;

        // Accoupled = Ac + Acaux - I;
        // Scoupled = S + Saux;

        // Ptcouped = prod(Accoupled, Scoupled);
        // Ptsoma = Pt + Ptaux;

        // std::cout <<std::scientific << Ptsoma(0, 0) / Ptcouped(0, 0) - 1.0 << " " << Ptsoma(0, 1) / Ptcouped(0, 1) - 1.0 << " " << Ptsoma(1, 0) / Ptcouped(1, 0) - 1.0 << " " << Ptsoma(1, 1) / Ptcouped(1, 1) - 1.0 << " " << std::endl;

        aaux(0) = du1aux_dx1;
        aaux(1) = du2aux_dx1;

        dAcaux_dx1(0, 0) = d2u1aux_dx1dx1;
        dAcaux_dx1(0, 1) = d2u1aux_dx1dx2;
        dAcaux_dx1(1, 0) = d2u2aux_dx1dx1;
        dAcaux_dx1(1, 1) = d2u2aux_dx1dx2;

        dAcaux_dx2(0, 0) = d2u1aux_dx1dx2;
        dAcaux_dx2(0, 1) = d2u1aux_dx2dx2;
        dAcaux_dx2(1, 0) = d2u2aux_dx1dx2;
        dAcaux_dx2(1, 1) = d2u2aux_dx2dx2;

        dEcaux_dx1 = 0.5 * (prod(trans(dAcaux_dx1), Acaux) + prod(trans(Acaux), dAcaux_dx1));
        dEcaux_dx2 = 0.5 * (prod(trans(dAcaux_dx2), Acaux) + prod(trans(Acaux), dAcaux_dx2));

        // if (ep == "EPD")
        // {
        //     dSaux_dx1(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx1(0, 0) + poisson * dEcaux_dx1(1, 1)));
        //     dSaux_dx1(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx1(1, 1) + poisson * dEcaux_dx1(0, 0)));
        //     dSaux_dx1(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx1(1, 0);
        //     dSaux_dx1(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx1(0, 1);

        //     dSaux_dx2(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx2(0, 0) + poisson * dEcaux_dx2(1, 1)));
        //     dSaux_dx2(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dEcaux_dx2(1, 1) + poisson * dEcaux_dx2(0, 0)));
        //     dSaux_dx2(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx2(1, 0);
        //     dSaux_dx2(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx2(0, 1);
        // }
        // else
        // {
        //     dSaux_dx1(0, 0) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx1(0, 0) + poisson * dEcaux_dx1(1, 1));
        //     dSaux_dx1(1, 1) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx1(1, 1) + poisson * dEcaux_dx1(0, 0));
        //     dSaux_dx1(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx1(1, 0);
        //     dSaux_dx1(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx1(0, 1);

        //     dSaux_dx2(0, 0) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx2(0, 0) + poisson * dEcaux_dx2(1, 1));
        //     dSaux_dx2(1, 1) = (young / (1.0 - poisson * poisson)) * (dEcaux_dx2(1, 1) + poisson * dEcaux_dx2(0, 0));
        //     dSaux_dx2(1, 0) = (young / (1.0 + poisson)) * dEcaux_dx2(1, 0);
        //     dSaux_dx2(0, 1) = (young / (1.0 + poisson)) * dEcaux_dx2(0, 1);
        // }

        // dPtaux_dx1 = prod(dAcaux_dx1, Saux) + prod(Acaux, dSaux_dx1);
        // dPtaux_dx2 = prod(dAcaux_dx2, Saux) + prod(Acaux, dSaux_dx2);

        // divPtaux(0) = dPtaux_dx1(0, 0) + dPtaux_dx2(0, 1);
        // divPtaux(1) = dPtaux_dx1(1, 0) + dPtaux_dx2(1, 1);

        grad_aaux(0, 0) = dAcaux_dx1(0, 0);
        grad_aaux(1, 0) = dAcaux_dx1(1, 0);
        grad_aaux(0, 1) = dAcaux_dx2(0, 0);
        grad_aaux(1, 1) = dAcaux_dx2(1, 0);

        unitX1(0) = 1.0;
        unitX1(1) = 0.0;

        term1 = prod(trans(Pt), aaux);

        term2 = prod(trans(Ptaux), a);

        term3 = -0.5 * (tensorContraction(S, Ecaux) + tensorContraction(Saux, Ec)) * unitX1;

        if (type == "STATIC")
        {
            term4 = 0.0;
        }
        else
        {
            accel(0) = 0.0;
            accel(1) = 0.0;
            for (int ij = 0; ij < connection_.size(); ij++)
            {
                accel += phi(ij) * connection_[ij]->getCurrentAcceleration();
            }
            term4 = density * (inner_prod(accel, aaux));
        }

        term5 = tensorContraction(Pt, grad_aaux);

        term6 = tensorContraction(Ptaux, grad_a);

        term7 = -1.0 * tensorContraction(S, dEcaux_dx1);

        term8 = -1.0 * tensorContraction(Saux, dEc_dx1);

        jIntegral(1) += 0.5 * youngLine * (inner_prod((term1 + term2 + term3), grad_q)) * weight * j0 * thickness;
    }
    return jIntegral;
}

double Element::tensorContraction(const bounded_matrix<double, 2, 2> &tensorA, const bounded_matrix<double, 2, 2> &tensorB)
{
    double scalar = 0.0;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            scalar += tensorA(i, j) * tensorB(i, j);
        }
    }

    return scalar;
}

bounded_vector<double, 2> Element::contributionDynamicJ_IntegralInitial(const int &hammerPoints, const std::string &ep, const double &rotation, const bounded_vector<double, 2> &tipCoordinates)
{
    bounded_vector<double, 2> jIntegral;
    jIntegral(0) = 0.0;
    jIntegral(1) = 0.0;

    int nnos = connection_.size();
    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);
    double thickness = mesh_->getThickness();

    double mi = young / (2.0 * (1.0 + poisson));
    double youngLine;

    double k;

    if (ep == "EPD")
    {
        k = 3.0 - 4.0 * poisson;
        youngLine = density * thickness * young / (1.0 - poisson * poisson);
    }
    else
    {
        k = (3.0 - poisson) / (1.0 + poisson);
        youngLine = density * thickness * young;
    }

    matrix<double> domainIntegrationPoints = hammerQuadrature(hammerPoints);

    bounded_matrix<double, 2, 2> matrixRotation;
    matrixRotation(0, 0) = cos(rotation);
    matrixRotation(0, 1) = sin(rotation);
    matrixRotation(1, 0) = -sin(rotation);
    matrixRotation(1, 1) = cos(rotation);

    const double pi = 3.14159265358979323846;

    for (int ihh = 0; ihh < hammerPoints; ihh++)
    {
        double xsi1 = domainIntegrationPoints(ihh, 0);
        double xsi2 = domainIntegrationPoints(ihh, 1);
        double weight = domainIntegrationPoints(ihh, 2);

        vector<double> phi = domainShapeFunction(xsi1, xsi2);
        matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2); //row = direction, column = node

        bounded_vector<double, 2> accelaux, coordP;
        accelaux(0) = 0.0;
        accelaux(1) = 0.0;
        coordP(0) = 0.0;
        coordP(1) = 0.0;

        for (int in = 0; in < nnos; in++)
        {
            accelaux += phi(in) * connection_[in]->getCurrentAcceleration();
            coordP += phi(in) * connection_[in]->getInitialCoordinate();
        }

        bounded_vector<double, 2> tipToPointAux = coordP - tipCoordinates;
        bounded_vector<double, 2> tipToPoint = prod(matrixRotation, tipToPointAux);
        const double radius = norm_2(tipToPointAux);
        const double teta = atan2(tipToPoint(1), tipToPoint(0));

        // bounded_vector<double, 2> accel = prod(matrixRotation, accelaux);

        bounded_vector<double, 2> accel = accelaux;

        bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2);
        double j0 = jacobianDeterminant(A0);
        // bounded_matrix<double, 2, 2> A0I = inverseMatrix(A0);
        bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2);
        // bounded_matrix<double, 2, 2> Ac_aux = prod(A1, A0I); //dyi / dxj

        // bounded_matrix<double, 2, 2> Ac, mataux;

        // mataux = prod(matrixRotation, Ac_aux);
        // Ac = prod(mataux, trans(matrixRotation)); //gradiente nos eixos da ponta da fissura

        //AUXILIARY STATE K1 = 1 E K2 = 0
        double du1aux_dr, du2aux_dr, du1aux_dteta, du2aux_dteta;

        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (cos(0.5 * teta) * (k - cos(teta)));
        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k - cos(teta)));

        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (-0.5 * sin(0.5 * teta) * (k - cos(teta)) + sin(teta) * cos(0.5 * teta));
        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k - cos(teta)) + sin(teta) * sin(0.5 * teta));

        double du1aux_dx1c, du2aux_dx1c;

        du1aux_dx1c = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;
        du2aux_dx1c = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        jIntegral(0) += density * (accel(0) * (A1(0, 0) - 1.0) + accel(1) * A1(1, 0)) * j0 * weight;

        //std::cout<<density * (accel(0) * A1(0, 0) + accel(1) * A1(1, 0)) * j0 * weight<<std::endl;

        // jIntegral(0) += 0.5 * youngLine * (accel(0) * du1aux_dx1c + accel(1) * du2aux_dx1c) * j0 * weight;

        // std::cout << "AUXILIAR I: " << " " << accel(0) << " " << accel(1) << " " << du1aux_dx1c << " " << du2aux_dx1c << std::endl;

        //AUXILIARY STATE K1 = 0 E K2 = 1
        du1aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (sin(0.5 * teta) * (k + 2.0 + cos(teta)));
        du2aux_dr = (sqrt(2.0 * pi * radius) / (8.0 * mi * pi * radius)) * (-cos(0.5 * teta) * (k - 2.0 + cos(teta)));

        du1aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * cos(0.5 * teta) * (k + 2.0 + cos(teta)) - sin(teta) * sin(0.5 * teta));
        du2aux_dteta = (sqrt(radius / (2.0 * pi)) / (2.0 * mi)) * (0.5 * sin(0.5 * teta) * (k - 2.0 + cos(teta)) + sin(teta) * cos(0.5 * teta));

        du1aux_dx1c = cos(teta) * du1aux_dr - (sin(teta) / radius) * du1aux_dteta;
        du2aux_dx1c = cos(teta) * du2aux_dr - (sin(teta) / radius) * du2aux_dteta;

        // std::cout << "AUXILIAR II: " <<  " " << accel(0) << " " << accel(1) << " " << du1aux_dx1c << " " << du2aux_dx1c << std::endl;

        //jIntegral(1) += 0.5 * youngLine * (accel(0) * du1aux_dx1c + accel(1) * du2aux_dx1c) * j0 * weight;
    }
    //  std::cout << jIntegral(0) << std::endl;
    return jIntegral;
}

bounded_vector<double, 2> Element::errorL2(const int &hammerPoints, const int &hammerPointsBlendZone)
{
    double sigmaA = -8.0, sigmaB = 6.0;
    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);

    double errorL2 = 0.0, errorH1 = 0.0, K1, K2;

    K1 = (1.0 - poisson * poisson) / young; //PARA EPD
    K2 = poisson / (1.0 - poisson);

    double a = 2.0; //raio interno
    double b = 6.0; //raio externos
    double C = (a * a * b * b * (sigmaA - sigmaB)) / (b * b - a * a);
    double b2 = (sigmaB * b * b - sigmaA * a * a) / (b * b - a * a);

    double radius, theta;

    int auxiliar = 0; //auxiliar pegar distancia e cell

    matrix<double> domainIntegrationPoints;
    if (insideBlendZone_.size() > 0) //the element has at least one point in blend zone;
    {
        domainIntegrationPoints = hammerQuadrature(hammerPointsBlendZone);
    }
    else
    {
        domainIntegrationPoints = hammerQuadrature(hammerPoints);
    }

    int nlocal = connection_.size();

    if (insideBlendZone_.size() > 0) //the element has at least one point in blend zone;
    {
        for (int ih = 0; ih < hammerPointsBlendZone; ih++)
        {
            if (insideBlendZone_[ih] == false)
            {
                double xsi1 = domainIntegrationPoints(ih, 0);
                double xsi2 = domainIntegrationPoints(ih, 1);
                double weight = domainIntegrationPoints(ih, 2);

                vector<double> phi = domainShapeFunction(xsi1, xsi2);

                bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
                double j0 = jacobianDeterminant(A0);
                bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
                A0I(0, 0) = A0(1, 1) / j0;
                A0I(1, 1) = A0(0, 0) / j0;
                A0I(0, 1) = -A0(0, 1) / j0;
                A0I(1, 0) = -A0(1, 0) / j0;
                bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2); //current configuration map
                bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                     //current deformation gradient
                identity_matrix<double> I(2);                                        //identity matrix

                bounded_vector<double, 2> initCoord, numDisp;
                initCoord(0) = 0.0;
                initCoord(1) = 0.0;
                numDisp(0) = 0.0;
                numDisp(1) = 0.0;

                for (int i = 0; i < nlocal; i++)
                {
                    initCoord += phi(i) * connection_[i]->getInitialCoordinate();
                    numDisp += phi(i) * connection_[i]->getCurrentDisplacement();
                }

                radius = norm_2(initCoord);
                theta = atan2(initCoord(1), initCoord(0));
                double urr = K1 * (b2 * (1.0 - K2) * radius - C / radius * (1.0 + K2));
                bounded_vector<double, 2> dif;
                dif(0) = cos(theta) * urr - numDisp(0);
                dif(1) = sin(theta) * urr - numDisp(1);

                // double l = 30.0;
                // double S = 1.0;
                // double pi = 3.1415926535897932384626433;
                // double ux; // = -1.0 / (young * S) * (-l * l / (pi * pi) * cos(pi * initCoord(0) / l) - 2.0 * l * initCoord(0) / (pi * pi) + l * l / (pi * pi));
                // if (initCoord(0) <= l / 3.0)
                // {
                //     ux = 1.0 / (young * S) * l * initCoord(0) / (3.0 * pi);
                // }
                // else if (initCoord(0) >= 2.0 * l / 3.0)
                // {
                //     ux = -1.0 / (young * S) * l * initCoord(0) / (3.0 * pi) + 1.0 / (young * S) * l * l / (3.0 * pi);
                // }
                // else
                // {
                //     ux = 1.0 / (young * S) * (l * l / (3.0 * 3.0 * pi * pi) * sin(pi / (l / 3.0) * (initCoord(0) - l / 3.0)) + l * l / (9.0 * pi));
                //     // std::cout << "NAO DEVERIA ENTRAR AQUI PARA O PROBLEMA EM QUESTÃO\n"
                //     //           << std::endl;
                // }
                // dif(0) = ux - numDisp(0);
                // dif(1) = 0.0 - numDisp(1);

                errorL2 += (dif(0) * dif(0) + dif(1) * dif(1)) * j0 * weight;

                bounded_matrix<double, 2, 2> gradU_num = Ac - I;

                double dur_dr = K1 * (b2 * (1.0 - K2) + C / (radius * radius) * (1.0 + K2));
                double du1_dr = cos(theta) * dur_dr;
                double du2_dr = sin(theta) * dur_dr;
                double du1_dtheta = -sin(theta) * urr;
                double du2_dtheta = cos(theta) * urr;

                bounded_matrix<double, 2, 2> gradU_teo;
                gradU_teo(0, 0) = cos(theta) * du1_dr - sin(theta) / radius * du1_dtheta;
                gradU_teo(0, 1) = sin(theta) * du1_dr + cos(theta) / radius * du1_dtheta;
                gradU_teo(1, 0) = cos(theta) * du2_dr - sin(theta) / radius * du2_dtheta;
                gradU_teo(1, 1) = sin(theta) * du2_dr + cos(theta) / radius * du2_dtheta;

                // gradU_teo(0, 0) = -1.0 / (young * S) * (l / pi * sin(pi * initCoord(0) / l) - 2.0 * l / (pi * pi));
                // gradU_teo(0, 1) = 0.0;
                // gradU_teo(1, 0) = 0.0;
                // gradU_teo(1, 1) = 0.0;

                // if (initCoord(0) <= l / 3.0)
                // {
                //     gradU_teo(0, 0) = 1.0 / (young * S) * l / (3.0 * pi);
                // }
                // else if (initCoord(0) >= 2.0 * l / 3.0)
                // {
                //     gradU_teo(0, 0) = -1.0 / (young * S) * l / (3.0 * pi);
                // }
                // else
                // {
                //     gradU_teo(0, 0) = 1.0 / (young * S) * cos(pi / (l / 3.0) * (initCoord(0) - l / 3.0)) * l / (3.0 * pi);
                // }

                bounded_matrix<double, 2, 2> gradDif = gradU_teo - gradU_num;

                errorH1 += (gradDif(0, 0) * gradDif(0, 0) + gradDif(0, 1) * gradDif(0, 1) + gradDif(1, 0) * gradDif(1, 0) + gradDif(1, 1) * gradDif(1, 1)) * j0 * weight;
            }
            else
            {
                Cell *cell = incidenceCell_[auxiliar];                                   //célula em que o ponto de hammer está inserido
                const bounded_vector<double, 2> xsiGlobal = xsiIncidenceCell_[auxiliar]; //coordenadas adimensionais da célula que correspondem aos pontos de hammer nas coordenadas globais;
                const double blend = bValue_[auxiliar];                                  //valor da blend function no ponto de hammer, em que a distância assinalada foi interpolada;
                const bounded_vector<double, 2> dblend_dxsi = db_dxsiValue_[auxiliar];   //derivada de b em relação a xsiLocal = db_dDAhammer * dDAhammer_dxsiLocal, sendo DAhammer a distância assinalada interpolada;

                //Pontos de Hammer para integração do elemento finito triangular
                const double xsi1 = domainIntegrationPoints(ih, 0);
                const double xsi2 = domainIntegrationPoints(ih, 1);
                const double weight = domainIntegrationPoints(ih, 2);

                vector<double> phiLocal = domainShapeFunction(xsi1, xsi2);                 //funções de forma calculada no ponto de hammer
                matrix<double> dphi_dxsiLocal = domainDerivativeShapeFunction(xsi1, xsi2); //derivadas das funções de forma em relação a xsi local (direção, nó);

                //GLOBAL
                std::vector<ControlPoint *> cps = cell->getControlPoints();
                int nglobal = cps.size(); //quantidade de pontos de controle que definem a célula em que o ponto de hammer está inserido
                vector<double> wpc(nglobal);
                bounded_vector<int, 2> inc;
                inc = cps[nglobal - 1]->getINC();
                for (int i = 0; i < nglobal; i++)
                {
                    wpc(i) = cps[i]->getWeight();
                }                                                                                                                       //conectividade da célula
                const std::pair<vector<double>, matrix<double>> functionsGlobal = cell->shapeFunctionAndDerivates(xsiGlobal, wpc, inc); //functionsGlobal.first(i) são as funções de forma globais calculadas em xsiGlobal que houve incidência;
                                                                                                                                        //functionsGlobal.second(j, i) são as derivadas das funções de forma globais em relação a xsi global;

                bounded_matrix<double, 2, 2> Jlocal = referenceJacobianMatrix(xsi1, xsi2);
                bounded_matrix<double, 2, 2> JGlobalI = inverseMatrix(cell->referenceJacobianMatrix(functionsGlobal.second));
                bounded_matrix<double, 2, 2> M = trans(prod(JGlobalI, Jlocal)); //para transforma dPhiGlobal_dXsiGlobal em dPhiGlobal_dXsiLocal

                int nlocal = connection_.size();               //poderia estar de fora...
                int nnum = nlocal + nglobal;                   //novo número funções de forma que definem o elemento na blend zone
                vector<double> phiBlended(nnum, 0.0);          //novas funções de forma que definem o elemento
                matrix<double> dphi_dxsiBlended(2, nnum, 0.0); //derivadas das novas funções de forma em relação a xsi local

                //novas funções de forma relacionadas as funções do elemento local
                for (int i = 0; i < nlocal; i++)
                {
                    phiBlended(i) = (1.0 - blend) * phiLocal(i);
                    dphi_dxsiBlended(0, i) = -dblend_dxsi(0) * phiLocal(i) + (1.0 - blend) * dphi_dxsiLocal(0, i);
                    dphi_dxsiBlended(1, i) = -dblend_dxsi(1) * phiLocal(i) + (1.0 - blend) * dphi_dxsiLocal(1, i);
                }

                //novas funções de forma relacionadas relacionadas as funções do global
                for (int i = 0; i < nglobal; i++)
                {
                    phiBlended(nlocal + i) = blend * functionsGlobal.first(i);

                    bounded_vector<double, 2> dPhiG_dXsiG, dPhiG_dxsiL;
                    dPhiG_dXsiG(0) = functionsGlobal.second(0, i); //dPhiGlobal_dxsiGlobal1
                    dPhiG_dXsiG(1) = functionsGlobal.second(1, i); //dPhiGlobal_dxsiGlobal2
                    dPhiG_dxsiL = prod(M, dPhiG_dXsiG);

                    dphi_dxsiBlended(0, nlocal + i) = dblend_dxsi(0) * functionsGlobal.first(i) + blend * dPhiG_dxsiL(0);
                    dphi_dxsiBlended(1, nlocal + i) = dblend_dxsi(1) * functionsGlobal.first(i) + blend * dPhiG_dxsiL(1);
                }

                bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrixBlendZone(dphi_dxsiBlended, cps); //mapeamento da configuração inicial
                bounded_matrix<double, 2, 2> A1 = currentJacobianMatrixBlendZone(dphi_dxsiBlended, cps);   //mapeamento da configuração atual
                bounded_matrix<double, 2, 2> A0I = inverseMatrix(A0);                                      //inverse initial configuration map
                bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                                           //gradiente da função mudança de configuração
                identity_matrix<double> I(2);                                                              //identity matrix
                double j0 = jacobianDeterminant(Jlocal);

                bounded_vector<double, 2> initCoord, numDisp;
                initCoord(0) = 0.0;
                initCoord(1) = 0.0;
                numDisp(0) = 0.0;
                numDisp(1) = 0.0;
                for (int i = 0; i < nlocal; i++)
                {
                    initCoord += phiBlended(i) * connection_[i]->getInitialCoordinate();
                    numDisp += phiBlended(i) * connection_[i]->getCurrentDisplacement();
                }
                for (int i = 0; i < nglobal; i++)
                {
                    initCoord += phiBlended(i + nlocal) * cps[i]->getInitialCoordinate();
                    numDisp += phiBlended(i + nlocal) * cps[i]->getCurrentDisplacement();
                }

                radius = sqrt(initCoord(0) * initCoord(0) + initCoord(1) * initCoord(1));
                theta = atan2(initCoord(1), initCoord(0));

                double urr = K1 * (b2 * (1.0 - K2) * radius - C / radius * (1.0 + K2));

                bounded_vector<double, 2> dif;
                dif(0) = cos(theta) * urr - numDisp(0);
                dif(1) = sin(theta) * urr - numDisp(1);

                // double l = 30.0;
                // double S = 1.0;
                // double pi = 3.1415926535897932384626433;
                // double ux; // = -1.0 / (young * S) * (-l * l / (pi * pi) * cos(pi * initCoord(0) / l) - 2.0 * l * initCoord(0) / (pi * pi) + l * l / (pi * pi));

                // if (initCoord(0) <= l / 3.0)
                // {
                //     ux = 1.0 / (young * S) * l * initCoord(0) / (3.0 * pi);
                // }
                // else if (initCoord(0) >= 2.0 * l / 3.0)
                // {
                //     ux = -1.0 / (young * S) * l * initCoord(0) / (3.0 * pi) + 1.0 / (young * S) * l * l / (3.0 * pi);
                // }
                // else
                // {
                //     ux = 1.0 / (young * S) * (l * l / (3.0 * 3.0 * pi * pi) * sin(pi / (l / 3.0) * (initCoord(0) - l / 3.0)) + l * l / (9.0 * pi));
                //     // std::cout << "NAO DEVERIA ENTRAR AQUI PARA O PROBLEMA EM QUESTÃO\n"
                //     //           << std::endl;
                // }

                // dif(0) = ux - numDisp(0);
                // dif(1) = 0.0 - numDisp(1);

                errorL2 += (dif(0) * dif(0) + dif(1) * dif(1)) * j0 * weight;

                bounded_matrix<double, 2, 2> gradU_num = Ac - I;

                double dur_dr = K1 * (b2 * (1.0 - K2) + C / (radius * radius) * (1.0 + K2));
                double du1_dr = cos(theta) * dur_dr;
                double du2_dr = sin(theta) * dur_dr;
                double du1_dtheta = -sin(theta) * urr;
                double du2_dtheta = cos(theta) * urr;

                bounded_matrix<double, 2, 2> gradU_teo;
                gradU_teo(0, 0) = cos(theta) * du1_dr - sin(theta) / radius * du1_dtheta;
                gradU_teo(0, 1) = sin(theta) * du1_dr + cos(theta) / radius * du1_dtheta;
                gradU_teo(1, 0) = cos(theta) * du2_dr - sin(theta) / radius * du2_dtheta;
                gradU_teo(1, 1) = sin(theta) * du2_dr + cos(theta) / radius * du2_dtheta;

                // gradU_teo(0, 0) = -1.0 / (young * S) * (l / pi * sin(pi * initCoord(0) / l) - 2.0 * l / (pi * pi));
                // gradU_teo(0, 1) = 0.0;
                // gradU_teo(1, 0) = 0.0;
                // gradU_teo(1, 1) = 0.0;

                // if (initCoord(0) <= l / 3.0)
                // {
                //     gradU_teo(0, 0) = 1.0 / (young * S) * l / (3.0 * pi);
                // }
                // else if (initCoord(0) >= 2.0 * l / 3.0)
                // {
                //     gradU_teo(0, 0) = -1.0 / (young * S) * l / (3.0 * pi);
                // }
                // else
                // {
                //     gradU_teo(0, 0) = 1.0 / (young * S) * cos(pi / (l / 3.0) * (initCoord(0) - l / 3.0)) * l / (3.0 * pi);
                // }

                bounded_matrix<double, 2, 2> gradDif = gradU_teo - gradU_num;

                errorH1 += (gradDif(0, 0) * gradDif(0, 0) + gradDif(0, 1) * gradDif(0, 1) + gradDif(1, 0) * gradDif(1, 0) + gradDif(1, 1) * gradDif(1, 1)) * j0 * weight;

                auxiliar++;
            }
        }
    }
    else
    {
        for (int ih = 0; ih < hammerPoints; ih++)
        {
            double xsi1 = domainIntegrationPoints(ih, 0);
            double xsi2 = domainIntegrationPoints(ih, 1);
            double weight = domainIntegrationPoints(ih, 2);

            vector<double> phi = domainShapeFunction(xsi1, xsi2);

            bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi1, xsi2); //initial configuration map
            double j0 = jacobianDeterminant(A0);
            bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
            A0I(0, 0) = A0(1, 1) / j0;
            A0I(1, 1) = A0(0, 0) / j0;
            A0I(0, 1) = -A0(0, 1) / j0;
            A0I(1, 0) = -A0(1, 0) / j0;
            bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi1, xsi2); //current configuration map
            bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                     //current deformation gradient
            identity_matrix<double> I(2);                                        //identity matrix

            bounded_vector<double, 2> initCoord, numDisp;
            initCoord(0) = 0.0;
            initCoord(1) = 0.0;
            numDisp(0) = 0.0;
            numDisp(1) = 0.0;

            for (int i = 0; i < nlocal; i++)
            {
                initCoord += phi(i) * connection_[i]->getInitialCoordinate();
                numDisp += phi(i) * connection_[i]->getCurrentDisplacement();
            }

            radius = norm_2(initCoord);
            theta = atan2(initCoord(1), initCoord(0));
            double urr = K1 * (b2 * (1.0 - K2) * radius - C / radius * (1.0 + K2));
            bounded_vector<double, 2> dif;
            dif(0) = cos(theta) * urr - numDisp(0);
            dif(1) = sin(theta) * urr - numDisp(1);

            // double l = 30.0;
            // double S = 1.0;
            // double pi = 3.1415926535897932384626433;
            // double ux; // = -1.0 / (young * S) * (-l * l / (pi * pi) * cos(pi * initCoord(0) / l) - 2.0 * l * initCoord(0) / (pi * pi) + l * l / (pi * pi));

            // if (initCoord(0) <= l / 3.0)
            // {
            //     ux = 1.0 / (young * S) * l * initCoord(0) / (3.0 * pi);
            // }
            // else if (initCoord(0) >= 2.0 * l / 3.0)
            // {
            //     ux = -1.0 / (young * S) * l * initCoord(0) / (3.0 * pi) + 1.0 / (young * S) * l * l / (3.0 * pi);
            // }
            // else
            // {
            //     ux = 1.0 / (young * S) * (l * l / (3.0 * 3.0 * pi * pi) * sin(pi / (l / 3.0) * (initCoord(0) - l / 3.0)) + l * l / (9.0 * pi));
            //     // std::cout << "NAO DEVERIA ENTRAR AQUI PARA O PROBLEMA EM QUESTÃO\n"
            //     //           << std::endl;
            // }

            // dif(0) = ux - numDisp(0);
            // dif(1) = 0.0 - numDisp(1);

            errorL2 += (dif(0) * dif(0) + dif(1) * dif(1)) * j0 * weight;

            bounded_matrix<double, 2, 2> gradU_num = Ac - I;

            double dur_dr = K1 * (b2 * (1.0 - K2) + C / (radius * radius) * (1.0 + K2));
            double du1_dr = cos(theta) * dur_dr;
            double du2_dr = sin(theta) * dur_dr;
            double du1_dtheta = -sin(theta) * urr;
            double du2_dtheta = cos(theta) * urr;

            bounded_matrix<double, 2, 2> gradU_teo;
            gradU_teo(0, 0) = cos(theta) * du1_dr - sin(theta) / radius * du1_dtheta;
            gradU_teo(0, 1) = sin(theta) * du1_dr + cos(theta) / radius * du1_dtheta;
            gradU_teo(1, 0) = cos(theta) * du2_dr - sin(theta) / radius * du2_dtheta;
            gradU_teo(1, 1) = sin(theta) * du2_dr + cos(theta) / radius * du2_dtheta;

            // gradU_teo(0, 0) = -1.0 / (young * S) * (l / pi * sin(pi * initCoord(0) / l) - 2.0 * l / (pi * pi));
            // gradU_teo(0, 1) = 0.0;
            // gradU_teo(1, 0) = 0.0;
            // gradU_teo(1, 1) = 0.0;

            // if (initCoord(0) <= l / 3.0)
            // {
            //     gradU_teo(0, 0) = 1.0 / (young * S) * l / (3.0 * pi);
            // }
            // else if (initCoord(0) >= 2.0 * l / 3.0)
            // {
            //     gradU_teo(0, 0) = -1.0 / (young * S) * l / (3.0 * pi);
            // }
            // else
            // {
            //     gradU_teo(0, 0) = 1.0 / (young * S) * cos(pi / (l / 3.0) * (initCoord(0) - l / 3.0)) * l / (3.0 * pi);
            // }

            bounded_matrix<double, 2, 2> gradDif = gradU_teo - gradU_num;

            errorH1 += (gradDif(0, 0) * gradDif(0, 0) + gradDif(0, 1) * gradDif(0, 1) + gradDif(1, 0) * gradDif(1, 0) + gradDif(1, 1) * gradDif(1, 1)) * j0 * weight;
        }
    }
    bounded_vector<double, 2> norms;
    norms(0) = errorL2;
    norms(1) = errorH1;
    return norms;
}
