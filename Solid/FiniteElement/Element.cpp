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
    std::string elementType = mesh->getElementType();
    if (elementType == "T3")
    {
        order_ = 1;
    }
    else if (elementType == "T6")
    {
        order_ = 2;
    }
    else if (elementType == "T10")
    {
        order_ = 3;
    }

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

matrix<double> Element::hammerQuadrature()
{
    matrix<double> pointCoord(112, 3);

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

    return pointCoord;
    // pointCoord(0, 0) = 0.166666666666667;
    // pointCoord(1, 0) = 0.050643253661728;
    // pointCoord(2, 0) = 0.398713492676543;
    // pointCoord(3, 0) = 0.050643253661728;
    // pointCoord(4, 0) = 0.235071032052557;
    // pointCoord(5, 0) = 0.235071032052557;
    // pointCoord(6, 0) = 0.029857935894885;
    // pointCoord(7, 0) = 0.333333333333333;
    // pointCoord(8, 0) = 0.101286507323457;
    // pointCoord(9, 0) = 0.449356746338272;
    // pointCoord(10, 0) = 0.449356746338272;
    // pointCoord(11, 0) = 0.264928967947442;
    // pointCoord(12, 0) = 0.470142064105115;
    // pointCoord(13, 0) = 0.264928967947442;
    // pointCoord(14, 0) = 0.666666666666667;
    // pointCoord(15, 0) = 0.550643253661728;
    // pointCoord(16, 0) = 0.898713492676543;
    // pointCoord(17, 0) = 0.550643253661728;
    // pointCoord(18, 0) = 0.735071032052557;
    // pointCoord(19, 0) = 0.735071032052558;
    // pointCoord(20, 0) = 0.529857935894885;
    // pointCoord(21, 0) = 0.166666666666667;
    // pointCoord(22, 0) = 0.050643253661728;
    // pointCoord(23, 0) = 0.398713492676543;
    // pointCoord(24, 0) = 0.050643253661728;
    // pointCoord(25, 0) = 0.235071032052557;
    // pointCoord(26, 0) = 0.235071032052557;
    // pointCoord(27, 0) = 0.029857935894885;

    // pointCoord(0, 1) = 0.166666666666667;
    // pointCoord(1, 1) = 0.050643253661729;
    // pointCoord(2, 1) = 0.050643253661729;
    // pointCoord(3, 1) = 0.398713492676544;
    // pointCoord(4, 1) = 0.029857935894885;
    // pointCoord(5, 1) = 0.235071032052557;
    // pointCoord(6, 1) = 0.235071032052557;
    // pointCoord(7, 1) = 0.333333333333333;
    // pointCoord(8, 1) = 0.449356746338272;
    // pointCoord(9, 1) = 0.101286507323457;
    // pointCoord(10, 1) = 0.449356746338272;
    // pointCoord(11, 1) = 0.264928967947442;
    // pointCoord(12, 1) = 0.264928967947442;
    // pointCoord(13, 1) = 0.470142064105115;
    // pointCoord(14, 1) = 0.166666666666667;
    // pointCoord(15, 1) = 0.050643253661729;
    // pointCoord(16, 1) = 0.050643253661729;
    // pointCoord(17, 1) = 0.398713492676544;
    // pointCoord(18, 1) = 0.029857935894885;
    // pointCoord(19, 1) = 0.235071032052557;
    // pointCoord(20, 1) = 0.235071032052557;
    // pointCoord(21, 1) = 0.666666666666667;
    // pointCoord(22, 1) = 0.550643253661729;
    // pointCoord(23, 1) = 0.550643253661729;
    // pointCoord(24, 1) = 0.898713492676544;
    // pointCoord(25, 1) = 0.529857935894885;
    // pointCoord(26, 1) = 0.735071032052558;
    // pointCoord(27, 1) = 0.735071032052557;

    // pointCoord(0, 2) = 0.0281250000000000;
    // pointCoord(1, 2) = 0.0157423975681034;
    // pointCoord(2, 2) = 0.0157423975681034;
    // pointCoord(3, 2) = 0.0157423975681034;
    // pointCoord(4, 2) = 0.0165492690985632;
    // pointCoord(5, 2) = 0.0165492690985632;
    // pointCoord(6, 2) = 0.0165492690985632;
    // pointCoord(7, 2) = 0.0281250000000000;
    // pointCoord(8, 2) = 0.0157423975681034;
    // pointCoord(9, 2) = 0.0157423975681034;
    // pointCoord(10, 2) = 0.0157423975681034;
    // pointCoord(11, 2) = 0.0165492690985632;
    // pointCoord(12, 2) = 0.0165492690985632;
    // pointCoord(13, 2) = 0.0165492690985632;
    // pointCoord(14, 2) = 0.0281250000000000;
    // pointCoord(15, 2) = 0.0157423975681034;
    // pointCoord(16, 2) = 0.0157423975681034;
    // pointCoord(17, 2) = 0.0157423975681034;
    // pointCoord(18, 2) = 0.0165492690985632;
    // pointCoord(19, 2) = 0.0165492690985632;
    // pointCoord(20, 2) = 0.0165492690985632;
    // pointCoord(21, 2) = 0.0281250000000000;
    // pointCoord(22, 2) = 0.0157423975681034;
    // pointCoord(23, 2) = 0.0157423975681034;
    // pointCoord(24, 2) = 0.0157423975681034;
    // pointCoord(25, 2) = 0.0165492690985632;
    // pointCoord(26, 2) = 0.0165492690985632;
    // pointCoord(27, 2) = 0.0165492690985632;

    // return pointCoord;

    // if (order_ == 1)
    // {
    //     matrix<double> aux(1, 3, 0.0);
    //     aux(0, 0) = 0.333333333333333; //xsi1
    //     aux(0, 1) = 0.333333333333333; //xsi2
    //     aux(0, 2) = 0.5;               //weight
    //     hammer = aux;
    // }
    // else if (order_ == 2)
    // {
    //     matrix<double> aux(4, 3, 0.0);
    //     aux(0, 0) = 0.333333333333333;
    //     aux(0, 1) = 0.333333333333333;
    //     aux(0, 2) = -0.281250000000000;

    //     aux(1, 0) = 0.6;
    //     aux(1, 1) = 0.2;
    //     aux(1, 2) = -0.260416666666667;

    //     aux(2, 0) = 0.2;
    //     aux(2, 1) = 0.6;
    //     aux(2, 2) = -0.260416666666667;

    //     aux(3, 0) = 0.2;
    //     aux(3, 1) = 0.2;
    //     aux(3, 2) = -0.260416666666667;
    //     hammer = aux;
    // }
    // else if (order_ == 3 or order_ == 2)
    // {
    //     matrix<double> aux(7, 3, 0.0);
    //     aux(0, 0) = 0.333333333333333;
    //     aux(0, 1) = 0.333333333333333;
    //     aux(0, 2) = 0.112500000000000;

    //     aux(1, 0) = 0.797426985353087;
    //     aux(1, 1) = 0.101286507323456;
    //     aux(1, 2) = 0.062969590272414;

    //     aux(2, 0) = 0.101286507323456;
    //     aux(2, 1) = 0.797426985353087;
    //     aux(2, 2) = 0.062969590272414;

    //     aux(3, 0) = 0.101286507323456;
    //     aux(3, 1) = 0.101286507323456;
    //     aux(3, 2) = 0.062969590272414;

    //     aux(4, 0) = 0.470142064105115;
    //     aux(4, 1) = 0.470142064105115;
    //     aux(4, 2) = 0.066197076394253;

    //     aux(5, 0) = 0.059715871789770;
    //     aux(5, 1) = 0.470142064105115;
    //     aux(5, 2) = 0.066197076394253;

    //     aux(6, 0) = 0.470142064105115;
    //     aux(6, 1) = 0.059715871789770;
    //     aux(6, 2) = 0.066197076394253;
    //     hammer = aux;
    // }
    // matrix<double> hammer(12, 3, 0.0);

    // hammer(0, 0) = 0.501426509658179;
    // hammer(0, 1) = 0.249286745170910;
    // hammer(0, 2) = 0.05839313786319;

    // hammer(1, 0) = 0.249286745170910;
    // hammer(1, 1) = 0.249286745170910;
    // hammer(1, 2) = 0.05839313786319;

    // hammer(2, 0) = 0.249286745170910;
    // hammer(2, 1) = 0.501426509658179;
    // hammer(2, 2) = 0.05839313786319;

    // hammer(3, 0) = 0.873821971016996;
    // hammer(3, 1) = 0.063089014491502;
    // hammer(3, 2) = 0.025422453185104;

    // hammer(4, 0) = 0.063089014491502;
    // hammer(4, 1) = 0.063089014491502;
    // hammer(4, 2) = 0.025422453185104;

    // hammer(5, 0) = 0.063089014491502;
    // hammer(5, 1) = 0.873821971016996;
    // hammer(5, 2) = 0.025422453185104;

    // hammer(6, 0) = 0.053145049844816;
    // hammer(6, 1) = 0.310352451033785;
    // hammer(6, 2) = 0.041425537809187;

    // hammer(7, 0) = 0.310352451033785;
    // hammer(7, 1) = 0.636502499121399;
    // hammer(7, 2) = 0.041425537809187;

    // hammer(8, 0) = 0.636502499121399;
    // hammer(8, 1) = 0.053145049844816;
    // hammer(8, 2) = 0.041425537809187;

    // hammer(9, 0) = 0.310352451033785;
    // hammer(9, 1) = 0.053145049844816;
    // hammer(9, 2) = 0.041425537809187;

    // hammer(10, 0) = 0.636502499121399;
    // hammer(10, 1) = 0.310352451033785;
    // hammer(10, 2) = 0.041425537809187;

    // hammer(11, 0) = 0.053145049844816;
    // hammer(11, 1) = 0.636502499121399;
    // hammer(11, 2) = 0.041425537809187;

    // return hammer;
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
                                                                        const double &beta, const double &gama)
{
    std::vector<Cell *> cellAux; //vetor para identificar a quantidade de células diferentes em que o elemento está incidindo
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

    //a variável freedomDegree é = nº de nós do elemento finito + nº pontos de controle * (quantidade de células diferentes que estão sob o elemento finito triangular);

    int auxiliar = 0; //auxiliar pegar distancia e cell

    vector<double> rhs(2 * freedomDegree, 0.0);
    matrix<double> tangent(2 * freedomDegree, 2 * freedomDegree, 0.0);
    matrix<double> domainIntegrationPoints_ = hammerQuadrature(); 
    double young, poisson, density;
    mesh_->getMaterial()->setProperties(young, poisson, density);
    if (typeAnalyze == "STATIC")
    {
        density = 0.0;
    }
    double thickness = mesh_->getThickness();

    for (int ih = 0; ih < domainIntegrationPoints_.size1(); ih++)
    {
        if (insideBlendZone_[ih] == false) //o ponto de hammer não está na blend zone; o vetor insideBlendZone_ possui tamanho igual ao número de pontos de hammer que estão sendo utilizados;
        {
            // std::cout<<index_<<std::endl;
            double xsi1 = domainIntegrationPoints_(ih, 0);
            double xsi2 = domainIntegrationPoints_(ih, 1);
            double weight = domainIntegrationPoints_(ih, 2);

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
            // double j1 = jacobianDeterminant(A1);
            // bounded_matrix<double, 2, 2> A1I; //inverse current configuration map
            // A1I(0, 0) = A1(1, 1) / j1;
            // A1I(1, 1) = A1(0, 0) / j1;
            // A1I(1, 0) = -A1(1, 0) / j1;
            // A1I(0, 1) = -A1(0, 1) / j1;
            bounded_matrix<double, 2, 2> Ac = prod(A1, A0I); //current deformation gradient
            // double jac = jacobianDeterminant(Ac);
            // bounded_matrix<double, 2, 2> AcI; //inverse current deformation gradient - ACHO QUE NÃO É NECESSÁRIO(REVER)
            // AcI(0, 0) = Ac(1, 1) / jac;
            // AcI(1, 1) = Ac(0, 0) / jac;
            // AcI(1, 0) = -Ac(1, 0) / jac;
            // AcI(0, 1) = -Ac(0, 1) / jac;
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
                    double shape = 0.0;

                    if (typeAnalyze == "STATIC")
                    {
                        for (int m = 0; m < connection_.size(); m++)
                        {
                            shape += phi(m); //* (1.0 * shapeForce_(j) * step / (1.0 * numberOfStep));
                        }
                    }
                    else
                    {
                        for (int m = 0; m < connection_.size(); m++)
                        {
                            accel += phi(m) * connection_[m]->getCurrentAcceleration()(j);
                            shape += phi(m);
                        }
                    }

                    double m = density * phi(i) * accel; //inertial force

                    shape = (shape * shapeForce_(j) * step / numberOfStep / thickness) * phi(i);

                    // // bounded_vector<double, 2> b;
                    // // b(0) = 0.0;
                    // // b(1) = phi(i) * density * gravity_; //domain force

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
                            dS_dy(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(0, 0) + poisson * dE_dy(1, 1)));
                            dS_dy(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(1, 1) + poisson * dE_dy(0, 0)));
                            dS_dy(1, 0) = (young / (1.0 + poisson)) * dE_dy(1, 0);
                            dS_dy(0, 1) = (young / (1.0 + poisson)) * dE_dy(0, 1);

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
        else //insideBlendZone_[ih] == true, logo, o ponto de hammer está na blend zone
        {
            //os próximos 4 intens foram calculados em GlobalSolid::incidenceLocalxGlobal() e enviados para cada elemento;
            //os vetores utilizados abaixo possuem tamanho igual a quantidade de pontos de hammer que estão dentro da blend zone, por isso o uso da variável "auxiliar"...
            Cell *cell = incidenceCell_[auxiliar];                                   //célula em que o ponto de hammer está inserido
            const bounded_vector<double, 2> xsiGlobal = xsiIncidenceCell_[auxiliar]; //coordenadas adimensionais da célula que correspondem aos pontos de hammer nas coordenadas globais;
            const double blend = bValue_[auxiliar];                                  //valor da blend function no ponto de hammer, em que a distância assinalada foi interpolada;
            const bounded_vector<double, 2> dblend_dxsi = db_dxsiValue_[auxiliar];   //derivada de b em relação a xsiLocal = db_dDAhammer * dDAhammer_dxsiLocal, sendo DAhammer a distância assinalada interpolada;

            //Pontos de Hammer para integração do elemento finito triangular
            const double xsi1 = domainIntegrationPoints_(ih, 0);
            const double xsi2 = domainIntegrationPoints_(ih, 1);
            const double weight = domainIntegrationPoints_(ih, 2);

            vector<double> phiLocal = domainShapeFunction(xsi1, xsi2);                 //funções de forma calculada no ponto de hammer
            matrix<double> dphi_dxsiLocal = domainDerivativeShapeFunction(xsi1, xsi2); //derivadas das funções de forma em relação a xsi local (direção, nó);

            //GLOBAL
            std::vector<ControlPoint *> cps = cell->getControlPoints();                                                   //conectividade da célula
            const std::pair<vector<double>, matrix<double>> functionsGlobal = cell->shapeFunctionAndDerivates(xsiGlobal); //functionsGlobal.first(i) são as funções de forma globais calculadas em xsiGlobal que houve incidência;
                                                                                                                          //functionsGlobal.second(j, i) são as derivadas das funções de forma globais em relação a xsi global;

            bounded_matrix<double, 2, 2> Jlocal = referenceJacobianMatrix(xsi1, xsi2);
            bounded_matrix<double, 2, 2> JGlobalI = inverseMatrix(cell->referenceJacobianMatrix(functionsGlobal.second));
            bounded_matrix<double, 2, 2> M = trans(prod(JGlobalI, Jlocal)); //para transforma dPhiGlobal_dXsiGlobal em dPhiGlobal_dXsiLocal

            //Teste para ver se dPhiGlobal_dXsiLocal estava certo...
            // bounded_matrix<double, 2, 2> teste;
            // teste(0, 0) = 0.0;
            // teste(0, 1) = 0.0;
            // teste(1, 0) = 0.0;
            // teste(1, 1) = 0.0;
            // for (int i = 0; i < nglobal; i++)
            // {
            //     bounded_vector<double, 2> dPhiG_dXsiG, dPhiG_dxsiL, coord;
            //     dPhiG_dXsiG(0) = functionsGlobal.second(0, i); //dPhiGlobal_dxsiGlobal1
            //     dPhiG_dXsiG(1) = functionsGlobal.second(1, i); //dPhiGlobal_dxsiGlobal2
            //     dPhiG_dxsiL = prod(M, dPhiG_dXsiG);
            //     coord = cps[i]->getInitialCoordinate();
            //     teste(0, 0) += coord(0) * dPhiG_dxsiL(0);
            //     teste(0, 1) += coord(0) * dPhiG_dxsiL(1);
            //     teste(1, 0) += coord(1) * dPhiG_dxsiL(0);
            //     teste(1, 1) += coord(1) * dPhiG_dxsiL(1);
            // }
            // std::cout.precision(5);
            // std::cout << std::fixed;
            // std::cout << Jlocal(0, 0) - teste(0, 0) << " " << Jlocal(0, 1) - teste(0, 1) << " " << Jlocal(1, 0) - teste(1, 0) << " " << Jlocal(1, 1) - teste(1, 1) << std::endl;

            int nlocal = connection_.size(); //poderia estar de fora...
            int nglobal = cps.size(); //quantidade de pontos de controle que definem a célula em que o ponto de hammer está inserido
            int nnum = nlocal + nglobal; //novo número funções de forma que definem o elemento na blend zone
            vector<double> phiBlended(nnum, 0.0); //novas funções de forma que definem o elemento
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
            bounded_matrix<double, 2, 2> A0I = inverseMatrix(A0); //inverse initial configuration map
            bounded_matrix<double, 2, 2> Ac = prod(A1, A0I); //gradiente da função mudança de configuração
            identity_matrix<double> I(2);                                      //identity matrix
            bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor
            bounded_matrix<double, 2, 2> S;                                    //second piola kirchhoff stress tensor
            double j0 = jacobianDeterminant(Jlocal);                           //para multiplicar com o peso

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

            bounded_matrix<double, 2, 2> dA_dy; //deriavada de A1 em relação a y
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
                            dS_dy(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(0, 0) + poisson * dE_dy(1, 1)));
                            dS_dy(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) * ((1.0 - poisson) * dE_dy(1, 1) + poisson * dE_dy(0, 0)));
                            dS_dy(1, 0) = (young / (1.0 + poisson)) * dE_dy(1, 0);
                            dS_dy(0, 1) = (young / (1.0 + poisson)) * dE_dy(0, 1);

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

    return std::make_pair(rhs, tangent);
}

void Element::StressCalculate(const std::string &ep)
{
    if (order_ == 2)
    {
        double young, poisson, density;
        mesh_->getMaterial()->setProperties(young, poisson, density);
        for (int i = 0; i < 6; i++)
        {
            bounded_matrix<double, 6, 2> xsi;
            xsi(0, 0) = 1.0;
            xsi(0, 1) = 0.0;

            xsi(1, 0) = 0.0;
            xsi(1, 1) = 1.0;

            xsi(2, 0) = 0.0;
            xsi(2, 1) = 0.0;

            xsi(3, 0) = 0.5;
            xsi(3, 1) = 0.5;

            xsi(4, 0) = 0.0;
            xsi(4, 1) = 0.5;

            xsi(5, 0) = 0.5;
            xsi(5, 1) = 0.0;

            bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi(i, 0), xsi(i, 1)); //initial configuration map
            double j0 = jacobianDeterminant(A0);
            bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
            A0I(0, 0) = A0(1, 1) / j0;
            A0I(1, 1) = A0(0, 0) / j0;
            A0I(0, 1) = -A0(0, 1) / j0;
            A0I(1, 0) = -A0(1, 0) / j0;
            bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi(i, 0), xsi(i, 1));
            bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                   //current deformation gradient
            identity_matrix<double> I(2);                                      //identity matrix
            bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor
            bounded_matrix<double, 2, 2> S;                                    //second piola kirchhoff stress tensor

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

            bounded_vector<double, 4> cauchStress;

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

            connection_[i]->setStressState(cauchStress);
        }
    }
    else if (order_ == 1)
    {
        double young, poisson, density;
        mesh_->getMaterial()->setProperties(young, poisson, density);
        for (int i = 0; i < 3; i++)
        {
            bounded_matrix<double, 3, 2> xsi;
            xsi(0, 0) = 1.0;
            xsi(0, 1) = 0.0;

            xsi(1, 0) = 0.0;
            xsi(1, 1) = 1.0;

            xsi(2, 0) = 0.0;
            xsi(2, 1) = 0.0;

            bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi(i, 0), xsi(i, 1)); //initial configuration map
            double j0 = jacobianDeterminant(A0);
            bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
            A0I(0, 0) = A0(1, 1) / j0;
            A0I(1, 1) = A0(0, 0) / j0;
            A0I(0, 1) = -A0(0, 1) / j0;
            A0I(1, 0) = -A0(1, 0) / j0;
            bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi(i, 0), xsi(i, 1));
            bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                   //current deformation gradient
            identity_matrix<double> I(2);                                      //identity matrix
            bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor
            bounded_matrix<double, 2, 2> S;                                    //second piola kirchhoff stress tensor

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

            bounded_vector<double, 4> cauchStress;

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

            connection_[i]->setStressState(cauchStress);
        }
    }
    else if (order_ == 3)
    {
        double young, poisson, density;
        mesh_->getMaterial()->setProperties(young, poisson, density);
        for (int i = 0; i < 10; i++)
        {
            bounded_matrix<double, 10, 2> xsi;
            xsi(0, 0) = 1.0;
            xsi(0, 1) = 0.0;

            xsi(1, 0) = 0.0;
            xsi(1, 1) = 1.0;

            xsi(2, 0) = 0.0;
            xsi(2, 1) = 0.0;

            xsi(3, 0) = 0.666666666666667;
            xsi(3, 1) = 0.333333333333333;

            xsi(4, 0) = 0.333333333333333;
            xsi(4, 1) = 0.666666666666667;

            xsi(5, 0) = 0.0;
            xsi(5, 1) = 0.666666666666667;

            xsi(6, 0) = 0.0;
            xsi(6, 1) = 0.333333333333333;

            xsi(7, 0) = 0.333333333333333;
            xsi(7, 1) = 0.0;

            xsi(8, 0) = 0.666666666666667;
            xsi(8, 1) = 0.0;

            xsi(9, 0) = 0.333333333333333;
            xsi(9, 1) = 0.333333333333333;

            bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi(i, 0), xsi(i, 1)); //initial configuration map
            double j0 = jacobianDeterminant(A0);
            bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
            A0I(0, 0) = A0(1, 1) / j0;
            A0I(1, 1) = A0(0, 0) / j0;
            A0I(0, 1) = -A0(0, 1) / j0;
            A0I(1, 0) = -A0(1, 0) / j0;
            bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi(i, 0), xsi(i, 1));
            bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                   //current deformation gradient
            identity_matrix<double> I(2);                                      //identity matrix
            bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor
            bounded_matrix<double, 2, 2> S;                                    //second piola kirchhoff stress tensor

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

            bounded_vector<double, 4> cauchStress;

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

            connection_[i]->setStressState(cauchStress);
        }
    }
}

matrix<double> Element::massMatrix()
{
    matrix<double> mass(2 * connection_.size(), 2 * connection_.size(), 0.0);
    matrix<double> domainIntegrationPoints_ = hammerQuadrature();
    double auxiliar = mesh_->getMaterial()->getDensity() * mesh_->getThickness();
    double resul;

    for (int ih = 0; ih < domainIntegrationPoints_.size1(); ih++)
    {
        double xsi1 = domainIntegrationPoints_(ih, 0);
        double xsi2 = domainIntegrationPoints_(ih, 1);
        double weight = domainIntegrationPoints_(ih, 2);

        vector<double> phi = domainShapeFunction(xsi1, xsi2);

        for (int i = 0; i < connection_.size(); i++)
        {
            for (int k = 0; k < connection_.size(); k++)
            {
                resul = auxiliar * phi(i) * phi(k);
                mass(2 * i, 2 * k) += resul;
                mass(2 * i + 1, 2 * k + 1) += resul;
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
    vector<double> wpc(nGlobal);
    for (int i = 0; i < nGlobal; i++)
    {
        wpc(i) = cpsGLobal[i]->getWeight();
    }

    for (int i = 0; i < nGlobal; i++)
    {
        bounded_vector<double, 2> initialCoord = cpsGLobal[i]->getInitialCoordinate();
        dx1_dxsi1 += initialCoord(0) * dphi_dxsiBlended(0, i + nLocal) / wpc(i);
        dx1_dxsi2 += initialCoord(0) * dphi_dxsiBlended(1, i + nLocal) / wpc(i);
        dx2_dxsi1 += initialCoord(1) * dphi_dxsiBlended(0, i + nLocal) / wpc(i);
        dx2_dxsi2 += initialCoord(1) * dphi_dxsiBlended(1, i + nLocal) / wpc(i);
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
    vector<double> wpc(nGlobal);
    for (int i = 0; i < nGlobal; i++)
    {
        wpc(i) = cpsGLobal[i]->getWeight();
    }
    for (int i = 0; i < nGlobal; i++)
    {
        bounded_vector<double, 2> currentCoord = cpsGLobal[i]->getCurrentCoordinate();
        dx1_dxsi1 += currentCoord(0) * dphi_dxsiBlended(0, i + nLocal) / wpc(i);
        dx1_dxsi2 += currentCoord(0) * dphi_dxsiBlended(1, i + nLocal) / wpc(i);
        dx2_dxsi1 += currentCoord(1) * dphi_dxsiBlended(0, i + nLocal) / wpc(i);
        dx2_dxsi2 += currentCoord(1) * dphi_dxsiBlended(1, i + nLocal) / wpc(i);
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

vector<double> Element::diagonalMass()
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

    vector<double> mass(freedomDegree, 0.0);
    matrix<double> domainIntegrationPoints_ = hammerQuadrature();
    double auxcalc = mesh_->getMaterial()->getDensity() * mesh_->getThickness();

    for (int ih = 0; ih < domainIntegrationPoints_.size1(); ih++)
    {
        if (insideBlendZone_[ih] == false)
        {
            double xsi1 = domainIntegrationPoints_(ih, 0);
            double xsi2 = domainIntegrationPoints_(ih, 1);
            double weight = domainIntegrationPoints_(ih, 2);

            vector<double> phi = domainShapeFunction(xsi1, xsi2);

            for (int i = 0; i < connection_.size(); i++)
            {
                mass(i) += auxcalc * phi(i) * phi(i) * weight;
            }
        }
        else
        {
            Cell *cell = incidenceCell_[auxiliar];
            std::vector<ControlPoint *> cps = cell->getControlPoints();
            bounded_vector<double, 2> xsiGlobal = xsiIncidenceCell_[auxiliar];
            double blend = bValue_[auxiliar];

            const double xsi1 = domainIntegrationPoints_(ih, 0);
            const double xsi2 = domainIntegrationPoints_(ih, 1);
            const double weight = domainIntegrationPoints_(ih, 2);

            vector<double> phiLocal = domainShapeFunction(xsi1, xsi2);
            for (int i = 0; i < phiLocal.size(); i++)
            {
                double phiblend = (1.0 - blend) * phiLocal(i);
                mass(i) += auxcalc * phiblend * phiblend * weight;
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
                mass(aux + i) += auxcalc * phiblend * phiblend * weight;
            }

            auxiliar += 1;
        }
    }
    return mass;
}