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
    else
    {
        std::cout << "SELECT A VALID NUMBER FOR THE AMOUNT OF HAMMER POINTS." << std::endl;
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

    if (insideBlendZone_.size() > 0) //the element has at least one hammer point in blend zone
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

                        double shape = (shapeForce_(j) * step / numberOfStep / thickness) * phi(i);

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

                        double shape = (shapeForce_(j) * step) / (thickness * numberOfStep) * phiBlended(i);

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
    else //the element doesn't have no hammer point in blend zone
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

                    double shape = (shapeForce_(j) * step / numberOfStep / thickness) * phi(i);

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

                for (int i = 0; i < connection_.size(); i++)
                {
                    mass(i) += phi(i) * phi(i) * weight;
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
                    mass(aux + i) += phiblend * phiblend * weight;
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

            for (int i = 0; i < connection_.size(); i++)
            {
                mass(i) += phi(i) * phi(i) * weight;
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

bounded_vector<double, 2> Element::contributionJ_IntegralInitial(const int &gaussPoints, const int &side, const std::string &ep, const double &rotation, const bounded_vector<double, 2> &tipCoordinates)
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

        const double pi = 3.14159265359;

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

bounded_vector<double, 2> Element::contributionJ_IntegralFromRice(const int &gaussPoints, const int &side, const std::string &ep, const double &rotation, Node *tipNode)
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