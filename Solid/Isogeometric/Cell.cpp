#include "Cell.h"

Cell::Cell() {}

Cell::Cell(const int &index,
           Patch *patch,
           const std::vector<ControlPoint *> &controlPoints)
//    const bounded_vector<int, 2> &order,
//    Material *material,
//    const double &thickness)
{
    index_ = index;
    controlPoints_ = controlPoints;
    // material_ = material;
    // thickness_ = thickness;
    patch_ = patch;
    // order_ = order;
    shapeForce_(0) = 0.0;
    shapeForce_(1) = 0.0;
}

Cell::~Cell() {}

int Cell::getIndex()
{
    return index_;
}

Patch *Cell::getPatch()
{
    return patch_;
}

std::vector<ControlPoint *> Cell::getControlPoints()
{
    return controlPoints_;
}

ControlPoint *Cell::getControlPoint(const int &index)
{
    return controlPoints_[index];
}

// Material *Cell::getMaterial()
// {
//     return material_;
// }

void Cell::setShapeForce(const bounded_vector<double, 2> &shapeForce)
{
    shapeForce_ = shapeForce;
}

vector<double> Cell::shapeFunction(const bounded_vector<double, 2> &qxsi,
                                   const vector<double> &wpc,
                                   const bounded_vector<double, 2> inc)
{
    int degm = patch_->getDegree(0);
    int degn = patch_->getDegree(1);
    int npcm = patch_->getNpc_Dir(0);
    int npcn = patch_->getNpc_Dir(1);

    // int degm = 2;
    // int degn = 2;
    // int npcm = 121;
    // int npcn = 91;

    int dimu = degm + npcm + 1;
    int dimv = degn + npcn + 1;
    int numLocalBF = (degm + 1) * (degn + 1);

    vector<double> uknot(dimu);
    vector<double> vknot(dimv);
    vector<double> uLeft(degm + 1);
    vector<double> uRight(degm + 1);
    vector<double> vLeft(degn + 1);
    vector<double> vRight(degn + 1);
    vector<double> phiL1(degm + 1);
    vector<double> phiL2(degn + 1);
    vector<double> phit(numLocalBF);
    matrix<double> uBF(degm + 1, degm + 1); // stores de base functions in each direction
    matrix<double> vBF(degn + 1, degn + 1); // stores de base functions in each direction

    //vector<double> wpc(numLocalBF);
    vector<double> phiIso(numLocalBF);

    double saved, temp;

    uknot = patch_->getKnotVectorU();
    vknot = patch_->getKnotVectorV();

    // indexes from index space where a function related with first local point of element start
    int uind = inc(0);
    int vind = inc(1);

    // parametric coordinates that defines the element
    double u1 = uknot(uind);
    double u2 = uknot(uind + 1);
    double v1 = vknot(vind);
    double v2 = vknot(vind + 1);

    // relates integration space with the parametric one
    const double xsi1 = ((u2 - u1) * qxsi(0) + (u2 + u1)) * 0.5;
    const double xsi2 = ((v2 - v1) * qxsi(1) + (v2 + v1)) * 0.5;

    // for (int i = 0; i < numLocalBF; i++)
    // {

    //     wpc(i) = wpcs(i);
    // }

    uBF(0, 0) = 1.0;

    for (int j = 1; j < (degm + 1); j++)
    {
        uLeft(j) = xsi1 - uknot(uind + 1 - j);

        uRight(j) = uknot(uind + j) - xsi1;
        saved = 0.;

        for (int r = 0; r < j; r++)
        {
            //lower triangle (Piegl and Tiller pg. 68-71)
            uBF(j, r) = uRight(r + 1) + uLeft(j - r);
            temp = uBF(r, j - 1) / uBF(j, r);

            //upper triangle
            uBF(r, j) = saved + uRight(r + 1) * temp;
            saved = uLeft(j - r) * temp;
        }

        uBF(j, j) = saved;
    }

    for (int i = 0; i < degm + 1; i++)
    {
        phiL1(i) = uBF(i, degm);
    }

    vBF(0, 0) = 1.0;

    for (int j = 1; j < (degn + 1); j++)
    {
        vLeft(j) = xsi2 - vknot(vind + 1 - j);
        vRight(j) = vknot(vind + j) - xsi2;
        saved = 0.;

        for (int r = 0; r < j; r++)
        {
            //lower triangle (Piegl and Tiller pg. 68-71)
            vBF(j, r) = vRight(r + 1) + vLeft(j - r);
            temp = vBF(r, j - 1) / vBF(j, r);

            //upper triangle
            vBF(r, j) = saved + vRight(r + 1) * temp;
            saved = vLeft(j - r) * temp;
        }

        vBF(j, j) = saved;
    }

    for (int i = 0; i < degn + 1; i++)
    {
        phiL2(i) = vBF(i, degn);
    }

    int index = 0;
    for (int j = 0; j < degn + 1; j++)
    {
        for (int k = 0; k < degm + 1; k++)
        {
            phit(index) = phiL1(k) * phiL2(j);
            index++;
        }
    }

    double sum = 0.0;

    for (int i = 0; i < numLocalBF; i++)
    {
        sum += phit(i) * wpc(i);
    }

    for (int i = 0; i < numLocalBF; i++)
    {
        phiIso(i) = phit(i) * wpc(i) / sum;
    }
    return phiIso;
}
std::pair<vector<double>, matrix<double>> Cell::shapeFunctionAndDerivates(const bounded_vector<double, 2> &qxsi,
                                                                          const vector<double> &wpc,
                                                                          const bounded_vector<double, 2> inc)
{
    int degm = patch_->getDegree(0);
    int degn = patch_->getDegree(1);
    int npcm = patch_->getNpc_Dir(0);
    int npcn = patch_->getNpc_Dir(1);

    // int degm = 2;
    // int degn = 2;
    // int npcm = 121;
    // int npcn = 91;

    int dimu = degm + npcm + 1;
    int dimv = degn + npcn + 1;
    int numLocalBF = (degm + 1) * (degn + 1);

    vector<double> uknot(dimu);
    vector<double> vknot(dimv);
    vector<double> uLeft(degm + 1);
    vector<double> uRight(degm + 1);
    vector<double> vLeft(degn + 1);
    vector<double> vRight(degn + 1);
    vector<double> phiL1(degm + 1);
    vector<double> phiL2(degn + 1);
    vector<double> phit(numLocalBF);
    vector<double> dphiL1(degm + 1);
    vector<double> dphiL2(degn + 1);
    matrix<double> dphit(2, numLocalBF);
    matrix<double> aux1(2, degm);           // stores the two lines more recently computed
    matrix<double> aux2(2, degn);           // stores the two lines more recently computed
    matrix<double> uBF(degm + 1, degm + 1); // stores de base functions in each direction
    matrix<double> vBF(degn + 1, degn + 1); // stores de base functions in each direction

    //vector<double> wpc(numLocalBF);
    vector<double> phiIso(numLocalBF);
    matrix<double> dphiIso(2, numLocalBF);

    double saved, temp;

    uknot = patch_->getKnotVectorU();
    vknot = patch_->getKnotVectorV();

    // indexes from index space where a function related with first local point of element start
    int uind = inc(0);
    int vind = inc(1);

    // parametric coordinates that defines the element
    double u1 = uknot(uind);
    double u2 = uknot(uind + 1);
    double v1 = vknot(vind);
    double v2 = vknot(vind + 1);

    // relates integration space with the parametric one
    const double xsi1 = ((u2 - u1) * qxsi(0) + (u2 + u1)) * 0.5;
    const double xsi2 = ((v2 - v1) * qxsi(1) + (v2 + v1)) * 0.5;

    // for (int i = 0; i < numLocalBF; i++)
    // {

    //     wpc(i) = wpcs(i);
    // }

    uBF(0, 0) = 1.0;

    for (int j = 1; j < (degm + 1); j++)
    {
        uLeft(j) = xsi1 - uknot(uind + 1 - j);
        uRight(j) = uknot(uind + j) - xsi1;
        saved = 0.;

        for (int r = 0; r < j; r++)
        {
            //lower triangle (Piegl and Tiller pg. 68-71)
            uBF(j, r) = uRight(r + 1) + uLeft(j - r);
            temp = uBF(r, j - 1) / uBF(j, r);

            //upper triangle
            uBF(r, j) = saved + uRight(r + 1) * temp;
            saved = uLeft(j - r) * temp;
        }

        uBF(j, j) = saved;
    }

    for (int i = 0; i < degm + 1; i++)
    {
        phiL1(i) = uBF(i, degm);
    }

    vBF(0, 0) = 1.0;

    for (int j = 1; j < (degn + 1); j++)
    {
        vLeft(j) = xsi2 - vknot(vind + 1 - j);
        vRight(j) = vknot(vind + j) - xsi2;
        saved = 0.;

        for (int r = 0; r < j; r++)
        {
            //lower triangle (Piegl and Tiller pg. 68-71)
            vBF(j, r) = vRight(r + 1) + vLeft(j - r);
            temp = vBF(r, j - 1) / vBF(j, r);

            //upper triangle
            vBF(r, j) = saved + vRight(r + 1) * temp;
            saved = vLeft(j - r) * temp;
        }

        vBF(j, j) = saved;
    }

    for (int i = 0; i < degn + 1; i++)
    {
        phiL2(i) = vBF(i, degn);
    }

    int index = 0;
    for (int j = 0; j < degn + 1; j++)
    {
        for (int k = 0; k < degm + 1; k++)
        {
            phit(index) = phiL1(k) * phiL2(j);
            index++;
        }
    }

    double sum = 0.0;

    for (int i = 0; i < numLocalBF; i++)
    {
        sum += phit(i) * wpc(i);
    }

    for (int i = 0; i < numLocalBF; i++)
    {
        phiIso(i) = phit(i) * wpc(i) / sum;
    }

    //DERIVATES

    int s1, s2, rk, pk, j1, j2, cor, k;
    double d;

    for (int r = 0; r < (degm + 1); r++)
    {
        s1 = 0;
        s2 = 1;
        aux1(0, 0) = 1.0;
        k = 1;
        d = 0.0;
        rk = r - k;
        pk = degm - k;

        if (r >= k)
        {
            aux1(s2, 0) = aux1(s1, 0) / uBF(pk + 1, rk);
            d = aux1(s2, 0) * uBF(rk, pk);
        }

        if (rk >= -1)
        {
            j1 = 1;
        }
        else
        {
            j1 = -rk;
        }

        if ((r - 1) <= pk)
        {
            j2 = k - 1;
        }
        else
        {
            j2 = degm - r;
        }

        for (int j = j1; j <= j2; j++)
        {
            aux1(s2, j) = (aux1(s1, j) - aux1(s1, j - 1)) / uBF(pk + 1, rk + j);

            d = d + aux1(s2, j) * uBF(rk + j, pk);
        }

        if (r <= pk)
        {
            aux1(s2, k) = -aux1(s1, k - 1) / uBF(pk + 1, r);
            d = d + aux1(s2, k) * uBF(r, pk);
        }

        dphiL1(r) = d;

        int j = s1;
        s1 = s2;
        s2 = j;
    }

    for (int i = 0; i < (degm + 1); i++)
    {
        dphiL1(i) = dphiL1(i) * degm;
    }

    // derivatives v direction;
    for (int r = 0; r < (degn + 1); r++)
    {
        s1 = 0;
        s2 = 1;
        aux2(0, 0) = 1.0;
        k = 1;
        d = 0.0;
        rk = r - k;
        pk = degn - k;

        if (r >= k)
        {
            aux2(s2, 0) = aux2(s1, 0) / vBF(pk + 1, rk);
            d = aux2(s2, 0) * vBF(rk, pk);
        }

        if (rk >= -1)
        {
            j1 = 1;
        }
        else
        {
            j1 = -rk;
        }

        if ((r - 1) <= pk)
        {
            j2 = k - 1;
        }
        else
        {
            j2 = degn - r;
        }

        for (int j = j1; j <= j2; j++)
        {
            aux2(s2, j) = (aux2(s1, j) - aux2(s1, j - 1)) / vBF(pk + 1, rk + j);

            d = d + aux2(s2, j) * vBF(rk + j, pk);
        }

        if (r <= pk)
        {
            aux2(s2, k) = -aux2(s1, k - 1) / vBF(pk + 1, r);
            d = d + aux2(s2, k) * vBF(r, pk);
        }

        dphiL2(r) = d;

        int j = s1;
        s1 = s2;
        s2 = j;
    }

    for (int i = 0; i < (degn + 1); i++)
    {
        dphiL2(i) = dphiL2(i) * degn;
    }

    index = 0;
    for (int j = 0; j < (degn + 1); j++)
    {
        for (int k = 0; k < (degm + 1); k++)
        {
            dphit(0, index) = dphiL1(k) * phiL2(j);
            dphit(1, index) = phiL1(k) * dphiL2(j);

            index++;
        }
    }

    sum = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int i = 0; i < numLocalBF; i++)
    {
        sum += phit(i) * wpc(i);
        sum1 += dphit(0, i) * wpc(i);
        sum2 += dphit(1, i) * wpc(i);
    }

    bounded_vector<int, 2> auxiliar = patch_->getSpanNumber();

    for (int i = 0; i < numLocalBF; i++)
    {
        dphiIso(0, i) = (1.0 / (sum * sum * 2 * auxiliar(0))) * (dphit(0, i) * wpc(i) * sum - phit(i) * wpc(i) * sum1);
        dphiIso(1, i) = (1.0 / (sum * sum * 2 * auxiliar(1))) * (dphit(1, i) * wpc(i) * sum - phit(i) * wpc(i) * sum2);
    }

    // for (int i = 0; i < numLocalBF; i++)
    // {
    //     dphiIso(0, i) = (1/394)*dphiIso(0, i);
    //     dphiIso(1, i) = (1/14)*dphiIso(0, i);
    // }

    return std::make_pair(phiIso, dphiIso);
};

bounded_matrix<double, 2, 2> Cell::referenceJacobianMatrix(const matrix<double> &dphi_dxsi, const vector<double> &wpc)
{
    //matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);
    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    for (int i = 0; i < controlPoints_.size(); i++)
    {
        bounded_vector<double, 2> initialCoord = controlPoints_[i]->getInitialCoordinate();
        dx1_dxsi1 += initialCoord(0) * dphi_dxsi(0, i) / wpc(i);
        dx1_dxsi2 += initialCoord(0) * dphi_dxsi(1, i) / wpc(i);
        dx2_dxsi1 += initialCoord(1) * dphi_dxsi(0, i) / wpc(i);
        dx2_dxsi2 += initialCoord(1) * dphi_dxsi(1, i) / wpc(i);
    }

    bounded_matrix<double, 2, 2> referenceJacobianMatrix;
    referenceJacobianMatrix(0, 0) = dx1_dxsi1;
    referenceJacobianMatrix(1, 0) = dx2_dxsi1;
    referenceJacobianMatrix(0, 1) = dx1_dxsi2;
    referenceJacobianMatrix(1, 1) = dx2_dxsi2;

    return referenceJacobianMatrix;
}

bounded_matrix<double, 2, 2> Cell::currentJacobianMatrix(const matrix<double> &dphi_dxsi, const vector<double> &wpc)
{
    //matrix<double> dphi_dxsi = domainDerivativeShapeFunction(xsi1, xsi2);
    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    for (int i = 0; i < controlPoints_.size(); i++)
    {
        bounded_vector<double, 2> currentCoord = controlPoints_[i]->getCurrentCoordinate();
        dx1_dxsi1 += currentCoord(0) * dphi_dxsi(0, i) / wpc(i);
        dx1_dxsi2 += currentCoord(0) * dphi_dxsi(1, i) / wpc(i);
        dx2_dxsi1 += currentCoord(1) * dphi_dxsi(0, i) / wpc(i);
        dx2_dxsi2 += currentCoord(1) * dphi_dxsi(1, i) / wpc(i);
    }

    bounded_matrix<double, 2, 2> currentJacobianMatrix;
    currentJacobianMatrix(0, 0) = dx1_dxsi1;
    currentJacobianMatrix(1, 0) = dx2_dxsi1;
    currentJacobianMatrix(0, 1) = dx1_dxsi2;
    currentJacobianMatrix(1, 1) = dx2_dxsi2;

    return currentJacobianMatrix;
}

double Cell::jacobianDeterminant(const bounded_matrix<double, 2, 2> &jacobianMatrix)
{
    return (jacobianMatrix(0, 0) * jacobianMatrix(1, 1) - jacobianMatrix(0, 1) * jacobianMatrix(1, 0));
}

matrix<double> Cell::isoQuadrature()
{
    matrix<double> pointCoordIso(16, 3); //xsi1, xsi2, weight

    // pointCoordIso(0, 0) = 0.;
    // pointCoordIso(0, 1) = 0.;
    // pointCoordIso(0, 2) = 0.790123456790124;

    // pointCoordIso(1, 0) = 0.774596669241483;
    // pointCoordIso(1, 1) = 0.000000000000000;
    // pointCoordIso(1, 2) = 0.493827160493828;

    // pointCoordIso(2, 0) = -0.774596669241483;
    // pointCoordIso(2, 1) = 0.;
    // pointCoordIso(2, 2) = 0.493827160493828;

    // pointCoordIso(3, 0) = 0.;
    // pointCoordIso(3, 1) = 0.774596669241483;
    // pointCoordIso(3, 2) = 0.493827160493828;

    // pointCoordIso(4, 0) = 0.774596669241483;
    // pointCoordIso(4, 1) = 0.774596669241483;
    // pointCoordIso(4, 2) = 0.308641975308642;

    // pointCoordIso(5, 0) = -0.774596669241483;
    // pointCoordIso(5, 1) = 0.774596669241483;
    // pointCoordIso(5, 2) = 0.308641975308642;

    // pointCoordIso(6, 0) = 0.;
    // pointCoordIso(6, 1) = -0.774596669241483;
    // pointCoordIso(6, 2) = 0.493827160493828;

    // pointCoordIso(7, 0) = 0.774596669241483;
    // pointCoordIso(7, 1) = -0.774596669241483;
    // pointCoordIso(7, 2) = 0.308641975308642;

    // pointCoordIso(8, 0) = -0.774596669241483;
    // pointCoordIso(8, 1) = -0.774596669241483;
    // pointCoordIso(8, 2) = 0.308641975308642;

    pointCoordIso(0, 0) = -0.861136311594053;
    pointCoordIso(1, 0) = -0.861136311594053;
    pointCoordIso(2, 0) = -0.861136311594053;
    pointCoordIso(3, 0) = -0.861136311594053;
    pointCoordIso(4, 0) = 0.861136311594053;
    pointCoordIso(5, 0) = 0.861136311594053;
    pointCoordIso(6, 0) = 0.861136311594053;
    pointCoordIso(7, 0) = 0.861136311594053;
    pointCoordIso(8, 0) = -0.339981043584856;
    pointCoordIso(9, 0) = -0.339981043584856;
    pointCoordIso(10, 0) = -0.339981043584856;
    pointCoordIso(11, 0) = -0.339981043584856;
    pointCoordIso(12, 0) = 0.339981043584856;
    pointCoordIso(13, 0) = 0.339981043584856;
    pointCoordIso(14, 0) = 0.339981043584856;
    pointCoordIso(15, 0) = 0.339981043584856;

    pointCoordIso(0, 1) = -0.861136311594053;
    pointCoordIso(1, 1) = 0.861136311594053;
    pointCoordIso(2, 1) = -0.339981043584856;
    pointCoordIso(3, 1) = 0.339981043584856;
    pointCoordIso(4, 1) = -0.861136311594053;
    pointCoordIso(5, 1) = 0.861136311594053;
    pointCoordIso(6, 1) = -0.339981043584856;
    pointCoordIso(7, 1) = 0.339981043584856;
    pointCoordIso(8, 1) = -0.861136311594053;
    pointCoordIso(9, 1) = 0.861136311594053;
    pointCoordIso(10, 1) = -0.339981043584856;
    pointCoordIso(11, 1) = 0.339981043584856;
    pointCoordIso(12, 1) = -0.861136311594053;
    pointCoordIso(13, 1) = 0.861136311594053;
    pointCoordIso(14, 1) = -0.339981043584856;
    pointCoordIso(15, 1) = 0.339981043584856;

    pointCoordIso(0, 2) = 0.121002993285602;
    pointCoordIso(1, 2) = 0.121002993285602;
    pointCoordIso(2, 2) = 0.226851851851852;
    pointCoordIso(3, 2) = 0.226851851851852;
    pointCoordIso(4, 2) = 0.121002993285602;
    pointCoordIso(5, 2) = 0.121002993285602;
    pointCoordIso(6, 2) = 0.226851851851852;
    pointCoordIso(7, 2) = 0.226851851851852;
    pointCoordIso(8, 2) = 0.226851851851852;
    pointCoordIso(9, 2) = 0.226851851851852;
    pointCoordIso(10, 2) = 0.425293303010694;
    pointCoordIso(11, 2) = 0.425293303010694;
    pointCoordIso(12, 2) = 0.226851851851852;
    pointCoordIso(13, 2) = 0.226851851851852;
    pointCoordIso(14, 2) = 0.425293303010694;
    pointCoordIso(15, 2) = 0.425293303010694;

    // pointCoordIso(0, 0) = 0.000000000000000;
    // pointCoordIso(1, 0) = 0.000000000000000;
    // pointCoordIso(2, 0) = 0.000000000000000;
    // pointCoordIso(3, 0) = 0.000000000000000;
    // pointCoordIso(4, 0) = 0.000000000000000;
    // pointCoordIso(5, 0) = 0.538469310195683;
    // pointCoordIso(6, 0) = 0.538469310195683;
    // pointCoordIso(7, 0) = 0.538469310195683;
    // pointCoordIso(8, 0) = 0.538469310195683;
    // pointCoordIso(9, 0) = 0.538469310195683;
    // pointCoordIso(10, 0) = -0.538469310195683;
    // pointCoordIso(11, 0) = -0.538469310195683;
    // pointCoordIso(12, 0) = -0.538469310195683;
    // pointCoordIso(13, 0) = -0.538469310195683;
    // pointCoordIso(14, 0) = -0.538469310195683;
    // pointCoordIso(15, 0) = 0.906179845938664;
    // pointCoordIso(16, 0) = 0.906179845938664;
    // pointCoordIso(17, 0) = 0.906179845938664;
    // pointCoordIso(18, 0) = 0.906179845938664;
    // pointCoordIso(19, 0) = 0.906179845938664;
    // pointCoordIso(20, 0) = -0.906179845938664;
    // pointCoordIso(21, 0) = -0.906179845938664;
    // pointCoordIso(22, 0) = -0.906179845938664;
    // pointCoordIso(23, 0) = -0.906179845938664;
    // pointCoordIso(24, 0) = -0.906179845938664;

    // pointCoordIso(0, 1) = 0.000000000000000;
    // pointCoordIso(1, 1) = 0.538469310195683;
    // pointCoordIso(2, 1) = -0.538469310195683;
    // pointCoordIso(3, 1) = 0.906179845938664;
    // pointCoordIso(4, 1) = -0.906179845938664;
    // pointCoordIso(5, 1) = 0.000000000000000;
    // pointCoordIso(6, 1) = 0.538469310195683;
    // pointCoordIso(7, 1) = -0.538469310195683;
    // pointCoordIso(8, 1) = 0.906179845938664;
    // pointCoordIso(9, 1) = -0.906179845938664;
    // pointCoordIso(10, 1) = 0.000000000000000;
    // pointCoordIso(11, 1) = 0.538469310195683;
    // pointCoordIso(12, 1) = -0.538469310195683;
    // pointCoordIso(13, 1) = 0.906179845938664;
    // pointCoordIso(14, 1) = -0.906179845938664;
    // pointCoordIso(15, 1) = 0.000000000000000;
    // pointCoordIso(16, 1) = 0.538469310195683;
    // pointCoordIso(17, 1) = -0.538469310195683;
    // pointCoordIso(18, 1) = 0.906179845938664;
    // pointCoordIso(19, 1) = -0.906179845938664;
    // pointCoordIso(20, 1) = 0.000000000000000;
    // pointCoordIso(21, 1) = 0.538469310195683;
    // pointCoordIso(22, 1) = -0.538469310195683;
    // pointCoordIso(23, 1) = 0.906179845938664;
    // pointCoordIso(24, 1) = -0.906179845938664;

    // pointCoordIso(0, 2) = 0.323634567901235;
    // pointCoordIso(1, 2) = 0.27228653255075;
    // pointCoordIso(2, 2) = 0.27228653255075;
    // pointCoordIso(3, 2) = 0.134785072387521;
    // pointCoordIso(4, 2) = 0.134785072387521;
    // pointCoordIso(5, 2) = 0.272286532550750;
    // pointCoordIso(6, 2) = 0.229085404223991;
    // pointCoordIso(7, 2) = 0.229085404223991;
    // pointCoordIso(8, 2) = 0.1134;
    // pointCoordIso(9, 2) = 0.1134;
    // pointCoordIso(10, 2) = 0.272286532550750;
    // pointCoordIso(11, 2) = 0.229085404223991;
    // pointCoordIso(12, 2) = 0.229085404223991;
    // pointCoordIso(13, 2) = 0.1134;
    // pointCoordIso(14, 2) = 0.1134;
    // pointCoordIso(15, 2) = 0.134785072387521;
    // pointCoordIso(16, 2) = 0.1134;
    // pointCoordIso(17, 2) = 0.1134;
    // pointCoordIso(18, 2) = 0.056134348862429;
    // pointCoordIso(19, 2) = 0.056134348862429;
    // pointCoordIso(20, 2) = 0.134785072387521;
    // pointCoordIso(21, 2) = 0.1134;
    // pointCoordIso(22, 2) = 0.1134;
    // pointCoordIso(23, 2) = 0.056134348862429;
    // pointCoordIso(24, 2) = 0.056134348862429;

    return pointCoordIso;
}

std::pair<vector<double>, matrix<double>> Cell::cellContributions(const std::string &ep, const std::string &typeAnalyze,
                                                                  const int &step, const int &numberOfStep, const double &deltat,
                                                                  const double &beta, const double &gama)
{
    int npc = controlPoints_.size();

    vector<double> rhs(2 * npc, 0.0);
    matrix<double> tangent(2 * npc, 2 * npc, 0.0);
    vector<double> wpc(npc);
    bounded_vector<int, 2> inc_;
    matrix<double> domainIntegrationPoints_ = isoQuadrature();

    inc_ = controlPoints_[npc - 1]->getINC();

    for (int i = 0; i < npc; i++)
    {
        wpc(i) = controlPoints_[i]->getWeight();
    }

    double young, poisson, density;
    patch_->getMaterial()->setProperties(young, poisson, density);
    if (typeAnalyze == "STATIC")
    {
        density = 0.0;
    }
    double thickness = patch_->getThickness();

    for (int ih = 0; ih < domainIntegrationPoints_.size1(); ih++)
    {

        double xsi1 = domainIntegrationPoints_(ih, 0);
        double xsi2 = domainIntegrationPoints_(ih, 1);
        double weight = domainIntegrationPoints_(ih, 2);
        bounded_vector<double, 2> qxsi;
        qxsi(0) = xsi1;
        qxsi(1) = xsi2;

        std::pair<vector<double>, matrix<double>> functions;

        functions = shapeFunctionAndDerivates(qxsi, wpc, inc_);

        bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(functions.second, wpc); //initial configuration map
        double j0 = jacobianDeterminant(A0);
        bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
        A0I(0, 0) = A0(1, 1) / j0;
        A0I(1, 1) = A0(0, 0) / j0;
        A0I(0, 1) = -A0(0, 1) / j0;
        A0I(1, 0) = -A0(1, 0) / j0;

        bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(functions.second, wpc); //current configuration map
        double j1 = jacobianDeterminant(A1);
        bounded_matrix<double, 2, 2> A1I; //inverse current configuration map
        A1I(0, 0) = A1(1, 1) / j1;
        A1I(1, 1) = A1(0, 0) / j1;
        A1I(1, 0) = -A1(1, 0) / j1;
        A1I(0, 1) = -A1(0, 1) / j1;
        bounded_matrix<double, 2, 2> Ac = prod(A1, A0I); //current deformation gradient
        double jac = jacobianDeterminant(Ac);
        bounded_matrix<double, 2, 2> AcI; //inverse current deformation gradient
        AcI(0, 0) = Ac(1, 1) / jac;
        AcI(1, 1) = Ac(0, 0) / jac;
        AcI(1, 0) = -Ac(1, 0) / jac;
        AcI(0, 1) = -Ac(0, 1) / jac;
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

        for (int i = 0; i < controlPoints_.size(); i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (j == 0)
                {
                    dA_dy(0, 0) = functions.second(0, i);
                    dA_dy(0, 1) = functions.second(1, i);
                    dA_dy(1, 0) = 0.0;
                    dA_dy(1, 1) = 0.0;
                }
                else
                {
                    dA_dy(1, 0) = functions.second(0, i);
                    dA_dy(1, 1) = functions.second(1, i);
                    dA_dy(0, 0) = 0.0;
                    dA_dy(0, 1) = 0.0;
                }

                bounded_matrix<double, 2, 2> mat1 = prod(trans(A0I), trans(dA_dy));
                bounded_matrix<double, 2, 2> mat2 = prod(dA_dy, A0I);

                bounded_matrix<double, 2, 2> dE_dy = 0.5 * (prod(mat1, Ac) + prod(trans(Ac), mat2)); //first derivative of E regarding i,j

                double r = dE_dy(0, 0) * S(0, 0) + dE_dy(1, 1) * S(1, 1) + dE_dy(0, 1) * S(0, 1) + dE_dy(1, 0) * S(1, 0); //internal force

                double accel = 0.0;
                double shape = 0.0;

                //shape = phi(i) * shape;

                if (typeAnalyze == "STATIC")
                {
                    for (int m = 0; m < controlPoints_.size(); m++)
                    {
                        shape += functions.first(m); //* (1.0 * shapeForce_(j) * step / (1.0 * numberOfStep));
                    }
                }
                else
                {
                    for (int m = 0; m < controlPoints_.size(); m++)
                    {
                        accel += functions.first(m) * controlPoints_[m]->getCurrentAcceleration()(j);
                        shape += functions.first(m);
                    }
                }

                double m = density * functions.first(i) * accel; //inertial force

                shape = (shape * shapeForce_(j) * step / numberOfStep / thickness) * functions.first(i);

                rhs(2 * i + j) -= (r + m - shape) * weight * j0 * thickness;

                for (int k = 0; k < controlPoints_.size(); k++)
                {
                    for (int l = 0; l < 2; l++)
                    {
                        if (l == 0)
                        {
                            dA_dy2(0, 0) = functions.second(0, k);
                            dA_dy2(0, 1) = functions.second(1, k);
                            dA_dy2(1, 0) = 0.0;
                            dA_dy2(1, 1) = 0.0;
                        }
                        else
                        {
                            dA_dy2(1, 0) = functions.second(0, k);
                            dA_dy2(1, 1) = functions.second(1, k);
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

                        double v = dE_dy3(0, 0) * S(0, 0) + dE_dy3(1, 1) * S(1, 1) + dE_dy3(0, 1) * S(0, 1) + dE_dy3(1, 0) * S(1, 0) + //second part of equation 5.88
                                   dE_dy2(0, 0) * dS_dy(0, 0) + dE_dy2(1, 1) * dS_dy(1, 1) + dE_dy2(0, 1) * dS_dy(0, 1) + dE_dy2(1, 0) * dS_dy(1, 0);

                        tangent(2 * i + j, 2 * k + l) += v * weight * j0 * thickness;
                        if (j == l)
                        {
                            double mm = density * functions.first(i) * functions.first(k); //mass contribution
                            tangent(2 * i + j, 2 * k + l) += (mm / (beta * deltat * deltat)) * weight * j0 * thickness;
                        }
                    }
                }
            }
        }
    }
    return std::make_pair(rhs, tangent);
}

bounded_vector<double, 4> Cell::getCauchStress(const bounded_vector<double, 2> &qxsi, const std::string &ep) //(SIGMA11, SIGMA22, SIGMA33, SIGMA12)
{
    bounded_vector<double, 4> cauchStress;

    int npc = controlPoints_.size();
    vector<double> wpc(npc);
    bounded_vector<int, 2> inc_;

    inc_ = controlPoints_[npc - 1]->getINC();

    for (int i = 0; i < npc; i++)
    {
        wpc(i) = controlPoints_[i]->getWeight();
    }

    double young, poisson, density;
    patch_->getMaterial()->setProperties(young, poisson, density);

    std::pair<vector<double>, matrix<double>> functions;
    functions = shapeFunctionAndDerivates(qxsi, wpc, inc_); //retorna as funções de formas e derivadas calculadas nos pontos qxsi

    bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(functions.second, wpc); //initial configuration map
    double j0 = jacobianDeterminant(A0);
    bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
    A0I(0, 0) = A0(1, 1) / j0;
    A0I(1, 1) = A0(0, 0) / j0;
    A0I(0, 1) = -A0(0, 1) / j0;
    A0I(1, 0) = -A0(1, 0) / j0;

    bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(functions.second, wpc); //current configuration map
    bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                                //current deformation gradient
    identity_matrix<double> I(2);                                                   //identity matrix
    bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);              //current green strain tensor
    bounded_matrix<double, 2, 2> S;                                                 //second piola kirchhoff stress tensor

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

bounded_vector<double, 4> Cell::getGreen(const bounded_vector<double, 2> &qxsi, const std::string &ep) //(SIGMA11, SIGMA22, SIGMA33, SIGMA12)
{
    bounded_vector<double, 4> green;

    int npc = controlPoints_.size();
    vector<double> wpc(npc);
    bounded_vector<int, 2> inc_;

    inc_ = controlPoints_[npc - 1]->getINC();

    for (int i = 0; i < npc; i++)
    {
        wpc(i) = controlPoints_[i]->getWeight();
    }

    double young, poisson, density;
    patch_->getMaterial()->setProperties(young, poisson, density);

    std::pair<vector<double>, matrix<double>> functions;
    functions = shapeFunctionAndDerivates(qxsi, wpc, inc_); //retorna as funções de formas e derivadas calculadas nos pontos qxsi

    bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(functions.second, wpc); //initial configuration map
    double j0 = jacobianDeterminant(A0);
    bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
    A0I(0, 0) = A0(1, 1) / j0;
    A0I(1, 1) = A0(0, 0) / j0;
    A0I(0, 1) = -A0(0, 1) / j0;
    A0I(1, 0) = -A0(1, 0) / j0;

    bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(functions.second, wpc); //current configuration map
    bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                                //current deformation gradient
    identity_matrix<double> I(2);                                                   //identity matrix
    bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);              //current green strain tensor

    green(0) = Ec(0, 0);
    green(1) = Ec(1, 1);
    green(3) = Ec(0, 1);

    if (ep == "EPD")
    {
        green(2) = 0.0;
    }
    else //EPT
    {
        green(2) = poisson * (green(0) + green(1)) / (poisson - 1.0);
    }

    return green;
}

// void Element::StressCalculate(const std::string &ep)
// {
// if (elementType_ == "T6")
// {
// for (int i = 0; i < 6; i++)
// {
// bounded_matrix<double, 6, 2> xsi;
// xsi(0, 0) = 1.0;
// xsi(0, 1) = 0.0;

// xsi(1, 0) = 0.0;
// xsi(1, 1) = 1.0;

// xsi(2, 0) = 0.0;
// xsi(2, 1) = 0.0;

// xsi(3, 0) = 0.5;
// xsi(3, 1) = 0.5;

// xsi(4, 0) = 0.0;
// xsi(4, 1) = 0.5;

// xsi(5, 0) = 0.5;
// xsi(5, 1) = 0.0;

// bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi(i, 0), xsi(i, 1)); //initial configuration map
// double j0 = jacobianDeterminant(A0);
// bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
// A0I(0, 0) = A0(1, 1) / j0;
// A0I(1, 1) = A0(0, 0) / j0;
// A0I(0, 1) = -A0(0, 1) / j0;
// A0I(1, 0) = -A0(1, 0) / j0;
// bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi(i, 0), xsi(i, 1));
// bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                   //current deformation gradient
// identity_matrix<double> I(2);                                      //identity matrix
// bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor
// bounded_matrix<double, 2, 2> S;                                    //second piola kirchhoff stress tensor
// double young = material_->getYoung();
// double poisson = material_->getPoisson();

// if (ep == "EPD")
// {
// S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
// S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
// S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
// S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
// }
// else
// {
// S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
// S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
// S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
// S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
// }

// bounded_matrix<double, 2, 2> sigma;
// double jac = jacobianDeterminant(Ac);
// bounded_matrix<double, 2, 2> mat1;
// mat1 = prod(Ac, S);
// sigma = (1.0 / jac) * prod(mat1, trans(Ac));

// bounded_vector<double, 3> cauchy;

// cauchy(0) = sigma(0, 0);
// cauchy(1) = sigma(1, 1);
// cauchy(2) = sigma(0, 1);
// // std::cout << thickness_ << std::endl;

// getConnection()[i]->setStressState(cauchy);
// }
// }
// else if (elementType_ == "T3")
// {
// for (int i = 0; i < 3; i++)
// {
// bounded_matrix<double, 3, 2> xsi;
// xsi(0, 0) = 1.0;
// xsi(0, 1) = 0.0;

// xsi(1, 0) = 0.0;
// xsi(1, 1) = 1.0;

// xsi(2, 0) = 0.0;
// xsi(2, 1) = 0.0;

// bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi(i, 0), xsi(i, 1)); //initial configuration map
// double j0 = jacobianDeterminant(A0);
// bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
// A0I(0, 0) = A0(1, 1) / j0;
// A0I(1, 1) = A0(0, 0) / j0;
// A0I(0, 1) = -A0(0, 1) / j0;
// A0I(1, 0) = -A0(1, 0) / j0;
// bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi(i, 0), xsi(i, 1));
// bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                   //current deformation gradient
// identity_matrix<double> I(2);                                      //identity matrix
// bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor
// bounded_matrix<double, 2, 2> S;                                    //second piola kirchhoff stress tensor
// double young = material_->getYoung();
// double poisson = material_->getPoisson();

// if (ep == "EPD")
// {
// S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
// S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
// S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
// S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
// }
// else
// {
// S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
// S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
// S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
// S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
// }

// bounded_matrix<double, 2, 2> sigma;
// double jac = jacobianDeterminant(Ac);
// bounded_matrix<double, 2, 2> mat1;
// mat1 = prod(Ac, S);
// sigma = (1.0 / jac) * prod(mat1, trans(Ac));

// bounded_vector<double, 3> cauchy;

// cauchy(0) = sigma(0, 0);
// cauchy(1) = sigma(1, 1);
// cauchy(2) = sigma(0, 1);

// getConnection()[i]->setStressState(cauchy);
// }
// }
// else if (elementType_ == "T10")
// {
// for (int i = 0; i < 10; i++)
// {
// bounded_matrix<double, 10, 2> xsi;
// xsi(0, 0) = 1.0;
// xsi(0, 1) = 0.0;

// xsi(1, 0) = 0.0;
// xsi(1, 1) = 1.0;

// xsi(2, 0) = 0.0;
// xsi(2, 1) = 0.0;

// xsi(3, 0) = 2.0 / 3.0;
// xsi(3, 1) = 1.0 / 3.0;

// xsi(4, 0) = 1.0 / 3.0;
// xsi(4, 1) = 2.0 / 3.0;

// xsi(5, 0) = 0.0;
// xsi(5, 1) = 2.0 / 3.0;

// xsi(6, 0) = 0.0;
// xsi(6, 1) = 1.0 / 3.0;

// xsi(7, 0) = 1.0 / 3.0;
// xsi(7, 1) = 0.0;

// xsi(8, 0) = 2.0 / 3.0;
// xsi(8, 1) = 0.0;

// xsi(9, 0) = 1.0 / 3.0;
// xsi(9, 1) = 1.0 / 3.0;

// bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(xsi(i, 0), xsi(i, 1)); //initial configuration map
// double j0 = jacobianDeterminant(A0);
// bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
// A0I(0, 0) = A0(1, 1) / j0;
// A0I(1, 1) = A0(0, 0) / j0;
// A0I(0, 1) = -A0(0, 1) / j0;
// A0I(1, 0) = -A0(1, 0) / j0;
// bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(xsi(i, 0), xsi(i, 1));
// bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                   //current deformation gradient
// identity_matrix<double> I(2);                                      //identity matrix
// bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I); //current green strain tensor
// bounded_matrix<double, 2, 2> S;                                    //second piola kirchhoff stress tensor
// double young = material_->getYoung();
// double poisson = material_->getPoisson();

// if (ep == "EPD")
// {
// S(0, 0) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(0, 0) + poisson * Ec(1, 1));
// S(1, 1) = (young / ((1.0 + poisson) * (1.0 - 2.0 * poisson))) * ((1.0 - poisson) * Ec(1, 1) + poisson * Ec(0, 0));
// S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
// S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
// }
// else
// {
// S(0, 0) = (young / (1.0 - poisson * poisson)) * (Ec(0, 0) + poisson * Ec(1, 1));
// S(1, 1) = (young / (1.0 - poisson * poisson)) * (Ec(1, 1) + poisson * Ec(0, 0));
// S(1, 0) = (young / (1.0 + poisson)) * Ec(1, 0);
// S(0, 1) = (young / (1.0 + poisson)) * Ec(0, 1);
// }

// bounded_matrix<double, 2, 2> sigma;
// double jac = jacobianDeterminant(Ac);
// bounded_matrix<double, 2, 2> mat1;
// mat1 = prod(Ac, S);
// sigma = (1.0 / jac) * prod(mat1, trans(Ac));

// bounded_vector<double, 3> cauchy;

// cauchy(0) = sigma(0, 0);
// cauchy(1) = sigma(1, 1);
// cauchy(2) = sigma(0, 1);

// getConnection()[i]->setStressState(cauchy);
// }
// }
// }

// matrix<double> Cell::massMatrix()
// {
//     matrix<double> mass(2 * controlPoints_.size(), 2 * controlPoints_.size(), 0.0);
//     matrix<double> domainIntegrationPoints_ = hammerQuadrature(numberOfDomainIntegrationPoints_);
//     double auxiliar = material_->getDensity() * thickness_;
//     double resul;

//     for (int ih = 0; ih < numberOfDomainIntegrationPoints_; ih++)
//     {
//         double xsi1 = domainIntegrationPoints_(ih, 0);
//         double xsi2 = domainIntegrationPoints_(ih, 1);
//         double weight = domainIntegrationPoints_(ih, 2);

//         vector<double> phi = domainShapeFunction(xsi1, xsi2);

//         for (int i = 0; i < controlPoints_.size(); i++)
//         {

//             for (int k = 0; k < controlPoints_.size(); k++)
//             {
//                 resul = auxiliar * phi(i) * phi(k);
//                 mass(2 * i, 2 * k) += resul;
//                 mass(2 * i + 1, 2 * k + 1) += resul;
//             }
//         }
//     }
//     return mass;
// }