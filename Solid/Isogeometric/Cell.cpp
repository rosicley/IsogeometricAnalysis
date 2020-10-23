#include "Cell.h"

Cell::Cell() {}

Cell::Cell(const int &index,
           Patch *patch,
           const std::vector<ControlPoint *> &controlPoints)
{
    index_ = index;
    controlPoints_ = controlPoints;
    patch_ = patch;
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
//     return maerial_;
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
        phiIso(i) = phit(i) / sum;
    }
    return phiIso;
}
std::pair<vector<double>, matrix<double>> Cell::shapeFunctionAndDerivates(const bounded_vector<double, 2> &qxsi,
                                                                          const vector<double> &wpc,
                                                                          const bounded_vector<double, 2> inc) //,
                                                                                                               //   const vector<double> &wpc,
                                                                                                               //   const bounded_vector<double, 2> inc)
{
    int npc = controlPoints_.size();

    int degm = patch_->getDegree(0);
    int degn = patch_->getDegree(1);
    int npcm = patch_->getNpc_Dir(0);
    int npcn = patch_->getNpc_Dir(1);

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
        phiIso(i) = phit(i) / sum;
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
        dphiIso(0, i) = (1.0 / (sum * sum * 2 * auxiliar(0))) * (dphit(0, i) * sum - phit(i) * sum1);
        dphiIso(1, i) = (1.0 / (sum * sum * 2 * auxiliar(1))) * (dphit(1, i) * sum - phit(i) * sum2);
    }

    return std::make_pair(phiIso, dphiIso);
};

bounded_matrix<double, 2, 2> Cell::referenceJacobianMatrix(const matrix<double> &dphi_dxsi)
{
    int npc = controlPoints_.size();

    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    for (int i = 0; i < npc; i++)
    {
        bounded_vector<double, 2> initialCoord = controlPoints_[i]->getInitialCoordinate();
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

bounded_matrix<double, 2, 2> Cell::currentJacobianMatrix(const matrix<double> &dphi_dxsi)
{
    int npc = controlPoints_.size();

    double dx1_dxsi1 = 0.0, dx1_dxsi2 = 0.0, dx2_dxsi1 = 0.0, dx2_dxsi2 = 0.0;

    for (int i = 0; i < npc; i++)
    {
        bounded_vector<double, 2> currentCoord = controlPoints_[i]->getCurrentCoordinate();
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

double Cell::jacobianDeterminant(const bounded_matrix<double, 2, 2> &jacobianMatrix)
{
    return (jacobianMatrix(0, 0) * jacobianMatrix(1, 1) - jacobianMatrix(0, 1) * jacobianMatrix(1, 0));
}

matrix<double> Cell::isoQuadrature(const int &points)
{
    matrix<double> pointCoordIso(points * points, 3); //xsi1, xsi2, weight

    if (points == 3)
    {
        pointCoordIso(0, 0) = 0.;
        pointCoordIso(0, 1) = 0.;
        pointCoordIso(0, 2) = 0.790123456790124;

        pointCoordIso(1, 0) = 0.774596669241483;
        pointCoordIso(1, 1) = 0.000000000000000;
        pointCoordIso(1, 2) = 0.493827160493828;

        pointCoordIso(2, 0) = -0.774596669241483;
        pointCoordIso(2, 1) = 0.;
        pointCoordIso(2, 2) = 0.493827160493828;

        pointCoordIso(3, 0) = 0.;
        pointCoordIso(3, 1) = 0.774596669241483;
        pointCoordIso(3, 2) = 0.493827160493828;

        pointCoordIso(4, 0) = 0.774596669241483;
        pointCoordIso(4, 1) = 0.774596669241483;
        pointCoordIso(4, 2) = 0.308641975308642;

        pointCoordIso(5, 0) = -0.774596669241483;
        pointCoordIso(5, 1) = 0.774596669241483;
        pointCoordIso(5, 2) = 0.308641975308642;

        pointCoordIso(6, 0) = 0.;
        pointCoordIso(6, 1) = -0.774596669241483;
        pointCoordIso(6, 2) = 0.493827160493828;

        pointCoordIso(7, 0) = 0.774596669241483;
        pointCoordIso(7, 1) = -0.774596669241483;
        pointCoordIso(7, 2) = 0.308641975308642;

        pointCoordIso(8, 0) = -0.774596669241483;
        pointCoordIso(8, 1) = -0.774596669241483;
        pointCoordIso(8, 2) = 0.308641975308642;
    }
    else if (points == 4)
    {
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
    }
    else if (points == 5)
    {
        pointCoordIso(0, 0) = 0.000000000000000;
        pointCoordIso(1, 0) = 0.000000000000000;
        pointCoordIso(2, 0) = 0.000000000000000;
        pointCoordIso(3, 0) = 0.000000000000000;
        pointCoordIso(4, 0) = 0.000000000000000;
        pointCoordIso(5, 0) = 0.538469310195683;
        pointCoordIso(6, 0) = 0.538469310195683;
        pointCoordIso(7, 0) = 0.538469310195683;
        pointCoordIso(8, 0) = 0.538469310195683;
        pointCoordIso(9, 0) = 0.538469310195683;
        pointCoordIso(10, 0) = -0.538469310195683;
        pointCoordIso(11, 0) = -0.538469310195683;
        pointCoordIso(12, 0) = -0.538469310195683;
        pointCoordIso(13, 0) = -0.538469310195683;
        pointCoordIso(14, 0) = -0.538469310195683;
        pointCoordIso(15, 0) = 0.906179845938664;
        pointCoordIso(16, 0) = 0.906179845938664;
        pointCoordIso(17, 0) = 0.906179845938664;
        pointCoordIso(18, 0) = 0.906179845938664;
        pointCoordIso(19, 0) = 0.906179845938664;
        pointCoordIso(20, 0) = -0.906179845938664;
        pointCoordIso(21, 0) = -0.906179845938664;
        pointCoordIso(22, 0) = -0.906179845938664;
        pointCoordIso(23, 0) = -0.906179845938664;
        pointCoordIso(24, 0) = -0.906179845938664;

        pointCoordIso(0, 1) = 0.000000000000000;
        pointCoordIso(1, 1) = 0.538469310195683;
        pointCoordIso(2, 1) = -0.538469310195683;
        pointCoordIso(3, 1) = 0.906179845938664;
        pointCoordIso(4, 1) = -0.906179845938664;
        pointCoordIso(5, 1) = 0.000000000000000;
        pointCoordIso(6, 1) = 0.538469310195683;
        pointCoordIso(7, 1) = -0.538469310195683;
        pointCoordIso(8, 1) = 0.906179845938664;
        pointCoordIso(9, 1) = -0.906179845938664;
        pointCoordIso(10, 1) = 0.000000000000000;
        pointCoordIso(11, 1) = 0.538469310195683;
        pointCoordIso(12, 1) = -0.538469310195683;
        pointCoordIso(13, 1) = 0.906179845938664;
        pointCoordIso(14, 1) = -0.906179845938664;
        pointCoordIso(15, 1) = 0.000000000000000;
        pointCoordIso(16, 1) = 0.538469310195683;
        pointCoordIso(17, 1) = -0.538469310195683;
        pointCoordIso(18, 1) = 0.906179845938664;
        pointCoordIso(19, 1) = -0.906179845938664;
        pointCoordIso(20, 1) = 0.000000000000000;
        pointCoordIso(21, 1) = 0.538469310195683;
        pointCoordIso(22, 1) = -0.538469310195683;
        pointCoordIso(23, 1) = 0.906179845938664;
        pointCoordIso(24, 1) = -0.906179845938664;

        pointCoordIso(0, 2) = 0.323634567901235;
        pointCoordIso(1, 2) = 0.27228653255075;
        pointCoordIso(2, 2) = 0.27228653255075;
        pointCoordIso(3, 2) = 0.134785072387521;
        pointCoordIso(4, 2) = 0.134785072387521;
        pointCoordIso(5, 2) = 0.272286532550750;
        pointCoordIso(6, 2) = 0.229085404223991;
        pointCoordIso(7, 2) = 0.229085404223991;
        pointCoordIso(8, 2) = 0.1134;
        pointCoordIso(9, 2) = 0.1134;
        pointCoordIso(10, 2) = 0.272286532550750;
        pointCoordIso(11, 2) = 0.229085404223991;
        pointCoordIso(12, 2) = 0.229085404223991;
        pointCoordIso(13, 2) = 0.1134;
        pointCoordIso(14, 2) = 0.1134;
        pointCoordIso(15, 2) = 0.134785072387521;
        pointCoordIso(16, 2) = 0.1134;
        pointCoordIso(17, 2) = 0.1134;
        pointCoordIso(18, 2) = 0.056134348862429;
        pointCoordIso(19, 2) = 0.056134348862429;
        pointCoordIso(20, 2) = 0.134785072387521;
        pointCoordIso(21, 2) = 0.1134;
        pointCoordIso(22, 2) = 0.1134;
        pointCoordIso(23, 2) = 0.056134348862429;
        pointCoordIso(24, 2) = 0.056134348862429;
    }
    else if (points == 6)
    {
        pointCoordIso(0, 0) = -0.9324695142;
        pointCoordIso(1, 0) = -0.9324695142;
        pointCoordIso(2, 0) = -0.9324695142;
        pointCoordIso(3, 0) = -0.9324695142;
        pointCoordIso(4, 0) = -0.9324695142;
        pointCoordIso(5, 0) = -0.9324695142;
        pointCoordIso(6, 0) = -0.6612093864;
        pointCoordIso(7, 0) = -0.6612093864;
        pointCoordIso(8, 0) = -0.6612093864;
        pointCoordIso(9, 0) = -0.6612093864;
        pointCoordIso(10, 0) = -0.6612093864;
        pointCoordIso(11, 0) = -0.6612093864;
        pointCoordIso(12, 0) = -0.238619186;
        pointCoordIso(13, 0) = -0.238619186;
        pointCoordIso(14, 0) = -0.238619186;
        pointCoordIso(15, 0) = -0.238619186;
        pointCoordIso(16, 0) = -0.238619186;
        pointCoordIso(17, 0) = -0.238619186;
        pointCoordIso(18, 0) = 0.238619186;
        pointCoordIso(19, 0) = 0.238619186;
        pointCoordIso(20, 0) = 0.238619186;
        pointCoordIso(21, 0) = 0.238619186;
        pointCoordIso(22, 0) = 0.238619186;
        pointCoordIso(23, 0) = 0.238619186;
        pointCoordIso(24, 0) = 0.6612093864;
        pointCoordIso(25, 0) = 0.6612093864;
        pointCoordIso(26, 0) = 0.6612093864;
        pointCoordIso(27, 0) = 0.6612093864;
        pointCoordIso(28, 0) = 0.6612093864;
        pointCoordIso(29, 0) = 0.6612093864;
        pointCoordIso(30, 0) = 0.9324695142;
        pointCoordIso(31, 0) = 0.9324695142;
        pointCoordIso(32, 0) = 0.9324695142;
        pointCoordIso(33, 0) = 0.9324695142;
        pointCoordIso(34, 0) = 0.9324695142;
        pointCoordIso(35, 0) = 0.9324695142;

        pointCoordIso(0, 1) = -0.9324695142;
        pointCoordIso(1, 1) = -0.6612093864;
        pointCoordIso(2, 1) = -0.238619186;
        pointCoordIso(3, 1) = 0.238619186;
        pointCoordIso(4, 1) = 0.6612093864;
        pointCoordIso(5, 1) = 0.9324695142;
        pointCoordIso(6, 1) = -0.9324695142;
        pointCoordIso(7, 1) = -0.6612093864;
        pointCoordIso(8, 1) = -0.238619186;
        pointCoordIso(9, 1) = 0.238619186;
        pointCoordIso(10, 1) = 0.6612093864;
        pointCoordIso(11, 1) = 0.9324695142;
        pointCoordIso(12, 1) = -0.9324695142;
        pointCoordIso(13, 1) = -0.6612093864;
        pointCoordIso(14, 1) = -0.238619186;
        pointCoordIso(15, 1) = 0.238619186;
        pointCoordIso(16, 1) = 0.6612093864;
        pointCoordIso(17, 1) = 0.9324695142;
        pointCoordIso(18, 1) = -0.9324695142;
        pointCoordIso(19, 1) = -0.6612093864;
        pointCoordIso(20, 1) = -0.238619186;
        pointCoordIso(21, 1) = 0.238619186;
        pointCoordIso(22, 1) = 0.6612093864;
        pointCoordIso(23, 1) = 0.9324695142;
        pointCoordIso(24, 1) = -0.9324695142;
        pointCoordIso(25, 1) = -0.6612093864;
        pointCoordIso(26, 1) = -0.238619186;
        pointCoordIso(27, 1) = 0.238619186;
        pointCoordIso(28, 1) = 0.6612093864;
        pointCoordIso(29, 1) = 0.9324695142;
        pointCoordIso(30, 1) = -0.9324695142;
        pointCoordIso(31, 1) = -0.6612093864;
        pointCoordIso(32, 1) = -0.238619186;
        pointCoordIso(33, 1) = 0.238619186;
        pointCoordIso(34, 1) = 0.6612093864;
        pointCoordIso(35, 1) = 0.9324695142;

        pointCoordIso(0, 2) = 0.0293520816618528;
        pointCoordIso(1, 2) = 0.0618072933355744;
        pointCoordIso(2, 2) = 0.0801651172683080;
        pointCoordIso(3, 2) = 0.0801651172683080;
        pointCoordIso(4, 2) = 0.0618072933355744;
        pointCoordIso(5, 2) = 0.0293520816618528;
        pointCoordIso(6, 2) = 0.0618072933355744;
        pointCoordIso(7, 2) = 0.1301489125534340;
        pointCoordIso(8, 2) = 0.1688053670388390;
        pointCoordIso(9, 2) = 0.1688053670388390;
        pointCoordIso(10, 2) = 0.1301489125534340;
        pointCoordIso(11, 2) = 0.0618072933355744;
        pointCoordIso(12, 2) = 0.0801651172683080;
        pointCoordIso(13, 2) = 0.1688053670388390;
        pointCoordIso(14, 2) = 0.2189434500992700;
        pointCoordIso(15, 2) = 0.2189434500992700;
        pointCoordIso(16, 2) = 0.1688053670388390;
        pointCoordIso(17, 2) = 0.0801651172683080;
        pointCoordIso(18, 2) = 0.0801651172683080;
        pointCoordIso(19, 2) = 0.1688053670388390;
        pointCoordIso(20, 2) = 0.2189434500992700;
        pointCoordIso(21, 2) = 0.2189434500992700;
        pointCoordIso(22, 2) = 0.1688053670388390;
        pointCoordIso(23, 2) = 0.0801651172683080;
        pointCoordIso(24, 2) = 0.0618072933355744;
        pointCoordIso(25, 2) = 0.1301489125534340;
        pointCoordIso(26, 2) = 0.1688053670388390;
        pointCoordIso(27, 2) = 0.1688053670388390;
        pointCoordIso(28, 2) = 0.1301489125534340;
        pointCoordIso(29, 2) = 0.0618072933355744;
        pointCoordIso(30, 2) = 0.0293520816618528;
        pointCoordIso(31, 2) = 0.0618072933355744;
        pointCoordIso(32, 2) = 0.0801651172683080;
        pointCoordIso(33, 2) = 0.0801651172683080;
        pointCoordIso(34, 2) = 0.0618072933355744;
        pointCoordIso(35, 2) = 0.0293520816618528;
    }
    else if (points == 7)
    {
        pointCoordIso(0, 0) = -0.94910791230;
        pointCoordIso(1, 0) = -0.94910791230;
        pointCoordIso(2, 0) = -0.94910791230;
        pointCoordIso(3, 0) = -0.94910791230;
        pointCoordIso(4, 0) = -0.94910791230;
        pointCoordIso(5, 0) = -0.94910791230;
        pointCoordIso(6, 0) = -0.94910791230;
        pointCoordIso(7, 0) = -0.74153118550;
        pointCoordIso(8, 0) = -0.74153118550;
        pointCoordIso(9, 0) = -0.74153118550;
        pointCoordIso(10, 0) = -0.74153118550;
        pointCoordIso(11, 0) = -0.74153118550;
        pointCoordIso(12, 0) = -0.74153118550;
        pointCoordIso(13, 0) = -0.74153118550;
        pointCoordIso(14, 0) = -0.40584515130;
        pointCoordIso(15, 0) = -0.40584515130;
        pointCoordIso(16, 0) = -0.40584515130;
        pointCoordIso(17, 0) = -0.40584515130;
        pointCoordIso(18, 0) = -0.40584515130;
        pointCoordIso(19, 0) = -0.40584515130;
        pointCoordIso(20, 0) = -0.40584515130;
        pointCoordIso(21, 0) = 0.00000000000;
        pointCoordIso(22, 0) = 0.00000000000;
        pointCoordIso(23, 0) = 0.00000000000;
        pointCoordIso(24, 0) = 0.00000000000;
        pointCoordIso(25, 0) = 0.00000000000;
        pointCoordIso(26, 0) = 0.00000000000;
        pointCoordIso(27, 0) = 0.00000000000;
        pointCoordIso(28, 0) = 0.40584515130;
        pointCoordIso(29, 0) = 0.40584515130;
        pointCoordIso(30, 0) = 0.40584515130;
        pointCoordIso(31, 0) = 0.40584515130;
        pointCoordIso(32, 0) = 0.40584515130;
        pointCoordIso(33, 0) = 0.40584515130;
        pointCoordIso(34, 0) = 0.40584515130;
        pointCoordIso(35, 0) = 0.74153118550;
        pointCoordIso(36, 0) = 0.74153118550;
        pointCoordIso(37, 0) = 0.74153118550;
        pointCoordIso(38, 0) = 0.74153118550;
        pointCoordIso(39, 0) = 0.74153118550;
        pointCoordIso(40, 0) = 0.74153118550;
        pointCoordIso(41, 0) = 0.74153118550;
        pointCoordIso(42, 0) = 0.94910791230;
        pointCoordIso(43, 0) = 0.94910791230;
        pointCoordIso(44, 0) = 0.94910791230;
        pointCoordIso(45, 0) = 0.94910791230;
        pointCoordIso(46, 0) = 0.94910791230;
        pointCoordIso(47, 0) = 0.94910791230;
        pointCoordIso(48, 0) = 0.94910791230;

        pointCoordIso(0, 1) = -0.94910791230;
        pointCoordIso(1, 1) = -0.74153118550;
        pointCoordIso(2, 1) = -0.40584515130;
        pointCoordIso(3, 1) = 0.00000000000;
        pointCoordIso(4, 1) = 0.40584515130;
        pointCoordIso(5, 1) = 0.74153118550;
        pointCoordIso(6, 1) = 0.94910791230;
        pointCoordIso(7, 1) = -0.94910791230;
        pointCoordIso(8, 1) = -0.74153118550;
        pointCoordIso(9, 1) = -0.40584515130;
        pointCoordIso(10, 1) = 0.00000000000;
        pointCoordIso(11, 1) = 0.40584515130;
        pointCoordIso(12, 1) = 0.74153118550;
        pointCoordIso(13, 1) = 0.94910791230;
        pointCoordIso(14, 1) = -0.94910791230;
        pointCoordIso(15, 1) = -0.74153118550;
        pointCoordIso(16, 1) = -0.40584515130;
        pointCoordIso(17, 1) = 0.00000000000;
        pointCoordIso(18, 1) = 0.40584515130;
        pointCoordIso(19, 1) = 0.74153118550;
        pointCoordIso(20, 1) = 0.94910791230;
        pointCoordIso(21, 1) = -0.94910791230;
        pointCoordIso(22, 1) = -0.74153118550;
        pointCoordIso(23, 1) = -0.40584515130;
        pointCoordIso(24, 1) = 0.00000000000;
        pointCoordIso(25, 1) = 0.40584515130;
        pointCoordIso(26, 1) = 0.74153118550;
        pointCoordIso(27, 1) = 0.94910791230;
        pointCoordIso(28, 1) = -0.94910791230;
        pointCoordIso(29, 1) = -0.74153118550;
        pointCoordIso(30, 1) = -0.40584515130;
        pointCoordIso(31, 1) = 0.00000000000;
        pointCoordIso(32, 1) = 0.40584515130;
        pointCoordIso(33, 1) = 0.74153118550;
        pointCoordIso(34, 1) = 0.94910791230;
        pointCoordIso(35, 1) = -0.94910791230;
        pointCoordIso(36, 1) = -0.74153118550;
        pointCoordIso(37, 1) = -0.40584515130;
        pointCoordIso(38, 1) = 0.00000000000;
        pointCoordIso(39, 1) = 0.40584515130;
        pointCoordIso(40, 1) = 0.74153118550;
        pointCoordIso(41, 1) = 0.94910791230;
        pointCoordIso(42, 1) = -0.94910791230;
        pointCoordIso(43, 1) = -0.74153118550;
        pointCoordIso(44, 1) = -0.40584515130;
        pointCoordIso(45, 1) = 0.00000000000;
        pointCoordIso(46, 1) = 0.40584515130;
        pointCoordIso(47, 1) = 0.74153118550;
        pointCoordIso(48, 1) = 0.94910791230;

        pointCoordIso(0, 2) = 0.0167663564459182;
        pointCoordIso(1, 2) = 0.0362176431234162;
        pointCoordIso(2, 2) = 0.0494412511449538;
        pointCoordIso(3, 2) = 0.0541194307196297;
        pointCoordIso(4, 2) = 0.0494412511449538;
        pointCoordIso(5, 2) = 0.0362176431234162;
        pointCoordIso(6, 2) = 0.0167663564459182;
        pointCoordIso(7, 2) = 0.0362176431234162;
        pointCoordIso(8, 2) = 0.0782351059782272;
        pointCoordIso(9, 2) = 0.1067999237233840;
        pointCoordIso(10, 2) = 0.1169054370380620;
        pointCoordIso(11, 2) = 0.1067999237233840;
        pointCoordIso(12, 2) = 0.0782351059782272;
        pointCoordIso(13, 2) = 0.0362176431234162;
        pointCoordIso(14, 2) = 0.0494412511449538;
        pointCoordIso(15, 2) = 0.1067999237233840;
        pointCoordIso(16, 2) = 0.1457941874648330;
        pointCoordIso(17, 2) = 0.1595893761809270;
        pointCoordIso(18, 2) = 0.1457941874648330;
        pointCoordIso(19, 2) = 0.1067999237233840;
        pointCoordIso(20, 2) = 0.0494412511449538;
        pointCoordIso(21, 2) = 0.0541194307196297;
        pointCoordIso(22, 2) = 0.1169054370380620;
        pointCoordIso(23, 2) = 0.1595893761809270;
        pointCoordIso(24, 2) = 0.1746898791555790;
        pointCoordIso(25, 2) = 0.1595893761809270;
        pointCoordIso(26, 2) = 0.1169054370380620;
        pointCoordIso(27, 2) = 0.0541194307196297;
        pointCoordIso(28, 2) = 0.0494412511449538;
        pointCoordIso(29, 2) = 0.1067999237233840;
        pointCoordIso(30, 2) = 0.1457941874648330;
        pointCoordIso(31, 2) = 0.1595893761809270;
        pointCoordIso(32, 2) = 0.1457941874648330;
        pointCoordIso(33, 2) = 0.1067999237233840;
        pointCoordIso(34, 2) = 0.0494412511449538;
        pointCoordIso(35, 2) = 0.0362176431234162;
        pointCoordIso(36, 2) = 0.0782351059782272;
        pointCoordIso(37, 2) = 0.1067999237233840;
        pointCoordIso(38, 2) = 0.1169054370380620;
        pointCoordIso(39, 2) = 0.1067999237233840;
        pointCoordIso(40, 2) = 0.0782351059782272;
        pointCoordIso(41, 2) = 0.0362176431234162;
        pointCoordIso(42, 2) = 0.0167663564459182;
        pointCoordIso(43, 2) = 0.0362176431234162;
        pointCoordIso(44, 2) = 0.0494412511449538;
        pointCoordIso(45, 2) = 0.0541194307196297;
        pointCoordIso(46, 2) = 0.0494412511449538;
        pointCoordIso(47, 2) = 0.0362176431234162;
        pointCoordIso(48, 2) = 0.0167663564459182;
    }
    else if (points == 8)
    {
        pointCoordIso(0, 0) = -0.9602898564000;
        pointCoordIso(1, 0) = -0.9602898564000;
        pointCoordIso(2, 0) = -0.9602898564000;
        pointCoordIso(3, 0) = -0.9602898564000;
        pointCoordIso(4, 0) = -0.9602898564000;
        pointCoordIso(5, 0) = -0.9602898564000;
        pointCoordIso(6, 0) = -0.9602898564000;
        pointCoordIso(7, 0) = -0.9602898564000;
        pointCoordIso(8, 0) = -0.7966664774000;
        pointCoordIso(9, 0) = -0.7966664774000;
        pointCoordIso(10, 0) = -0.7966664774000;
        pointCoordIso(11, 0) = -0.7966664774000;
        pointCoordIso(12, 0) = -0.7966664774000;
        pointCoordIso(13, 0) = -0.7966664774000;
        pointCoordIso(14, 0) = -0.7966664774000;
        pointCoordIso(15, 0) = -0.7966664774000;
        pointCoordIso(16, 0) = -0.5255324099000;
        pointCoordIso(17, 0) = -0.5255324099000;
        pointCoordIso(18, 0) = -0.5255324099000;
        pointCoordIso(19, 0) = -0.5255324099000;
        pointCoordIso(20, 0) = -0.5255324099000;
        pointCoordIso(21, 0) = -0.5255324099000;
        pointCoordIso(22, 0) = -0.5255324099000;
        pointCoordIso(23, 0) = -0.5255324099000;
        pointCoordIso(24, 0) = -0.1834346424000;
        pointCoordIso(25, 0) = -0.1834346424000;
        pointCoordIso(26, 0) = -0.1834346424000;
        pointCoordIso(27, 0) = -0.1834346424000;
        pointCoordIso(28, 0) = -0.1834346424000;
        pointCoordIso(29, 0) = -0.1834346424000;
        pointCoordIso(30, 0) = -0.1834346424000;
        pointCoordIso(31, 0) = -0.1834346424000;
        pointCoordIso(32, 0) = 0.1834346424000;
        pointCoordIso(33, 0) = 0.1834346424000;
        pointCoordIso(34, 0) = 0.1834346424000;
        pointCoordIso(35, 0) = 0.1834346424000;
        pointCoordIso(36, 0) = 0.1834346424000;
        pointCoordIso(37, 0) = 0.1834346424000;
        pointCoordIso(38, 0) = 0.1834346424000;
        pointCoordIso(39, 0) = 0.1834346424000;
        pointCoordIso(40, 0) = 0.5255324099000;
        pointCoordIso(41, 0) = 0.5255324099000;
        pointCoordIso(42, 0) = 0.5255324099000;
        pointCoordIso(43, 0) = 0.5255324099000;
        pointCoordIso(44, 0) = 0.5255324099000;
        pointCoordIso(45, 0) = 0.5255324099000;
        pointCoordIso(46, 0) = 0.5255324099000;
        pointCoordIso(47, 0) = 0.5255324099000;
        pointCoordIso(48, 0) = 0.7966664774000;
        pointCoordIso(49, 0) = 0.7966664774000;
        pointCoordIso(50, 0) = 0.7966664774000;
        pointCoordIso(51, 0) = 0.7966664774000;
        pointCoordIso(52, 0) = 0.7966664774000;
        pointCoordIso(53, 0) = 0.7966664774000;
        pointCoordIso(54, 0) = 0.7966664774000;
        pointCoordIso(55, 0) = 0.7966664774000;
        pointCoordIso(56, 0) = 0.9602898564000;
        pointCoordIso(57, 0) = 0.9602898564000;
        pointCoordIso(58, 0) = 0.9602898564000;
        pointCoordIso(59, 0) = 0.9602898564000;
        pointCoordIso(60, 0) = 0.9602898564000;
        pointCoordIso(61, 0) = 0.9602898564000;
        pointCoordIso(62, 0) = 0.9602898564000;
        pointCoordIso(63, 0) = 0.9602898564000;

        pointCoordIso(0, 1) = -0.9602898564000;
        pointCoordIso(1, 1) = -0.7966664774000;
        pointCoordIso(2, 1) = -0.5255324099000;
        pointCoordIso(3, 1) = -0.1834346424000;
        pointCoordIso(4, 1) = 0.1834346424000;
        pointCoordIso(5, 1) = 0.5255324099000;
        pointCoordIso(6, 1) = 0.7966664774000;
        pointCoordIso(7, 1) = 0.9602898564000;
        pointCoordIso(8, 1) = -0.9602898564000;
        pointCoordIso(9, 1) = -0.7966664774000;
        pointCoordIso(10, 1) = -0.5255324099000;
        pointCoordIso(11, 1) = -0.1834346424000;
        pointCoordIso(12, 1) = 0.1834346424000;
        pointCoordIso(13, 1) = 0.5255324099000;
        pointCoordIso(14, 1) = 0.7966664774000;
        pointCoordIso(15, 1) = 0.9602898564000;
        pointCoordIso(16, 1) = -0.9602898564000;
        pointCoordIso(17, 1) = -0.7966664774000;
        pointCoordIso(18, 1) = -0.5255324099000;
        pointCoordIso(19, 1) = -0.1834346424000;
        pointCoordIso(20, 1) = 0.1834346424000;
        pointCoordIso(21, 1) = 0.5255324099000;
        pointCoordIso(22, 1) = 0.7966664774000;
        pointCoordIso(23, 1) = 0.9602898564000;
        pointCoordIso(24, 1) = -0.9602898564000;
        pointCoordIso(25, 1) = -0.7966664774000;
        pointCoordIso(26, 1) = -0.5255324099000;
        pointCoordIso(27, 1) = -0.1834346424000;
        pointCoordIso(28, 1) = 0.1834346424000;
        pointCoordIso(29, 1) = 0.5255324099000;
        pointCoordIso(30, 1) = 0.7966664774000;
        pointCoordIso(31, 1) = 0.9602898564000;
        pointCoordIso(32, 1) = -0.9602898564000;
        pointCoordIso(33, 1) = -0.7966664774000;
        pointCoordIso(34, 1) = -0.5255324099000;
        pointCoordIso(35, 1) = -0.1834346424000;
        pointCoordIso(36, 1) = 0.1834346424000;
        pointCoordIso(37, 1) = 0.5255324099000;
        pointCoordIso(38, 1) = 0.7966664774000;
        pointCoordIso(39, 1) = 0.9602898564000;
        pointCoordIso(40, 1) = -0.9602898564000;
        pointCoordIso(41, 1) = -0.7966664774000;
        pointCoordIso(42, 1) = -0.5255324099000;
        pointCoordIso(43, 1) = -0.1834346424000;
        pointCoordIso(44, 1) = 0.1834346424000;
        pointCoordIso(45, 1) = 0.5255324099000;
        pointCoordIso(46, 1) = 0.7966664774000;
        pointCoordIso(47, 1) = 0.9602898564000;
        pointCoordIso(48, 1) = -0.9602898564000;
        pointCoordIso(49, 1) = -0.7966664774000;
        pointCoordIso(50, 1) = -0.5255324099000;
        pointCoordIso(51, 1) = -0.1834346424000;
        pointCoordIso(52, 1) = 0.1834346424000;
        pointCoordIso(53, 1) = 0.5255324099000;
        pointCoordIso(54, 1) = 0.7966664774000;
        pointCoordIso(55, 1) = 0.9602898564000;
        pointCoordIso(56, 1) = -0.9602898564000;
        pointCoordIso(57, 1) = -0.7966664774000;
        pointCoordIso(58, 1) = -0.5255324099000;
        pointCoordIso(59, 1) = -0.1834346424000;
        pointCoordIso(60, 1) = 0.1834346424000;
        pointCoordIso(61, 1) = 0.5255324099000;
        pointCoordIso(62, 1) = 0.7966664774000;
        pointCoordIso(63, 1) = 0.9602898564000;

        pointCoordIso(0, 2) = 0.01024721654119470;
        pointCoordIso(1, 2) = 0.02251130659095380;
        pointCoordIso(2, 2) = 0.03175606455054590;
        pointCoordIso(3, 2) = 0.03671394848693700;
        pointCoordIso(4, 2) = 0.03671394848693700;
        pointCoordIso(5, 2) = 0.03175606455054590;
        pointCoordIso(6, 2) = 0.02251130659095380;
        pointCoordIso(7, 2) = 0.01024721654119470;
        pointCoordIso(8, 2) = 0.02251130659095380;
        pointCoordIso(9, 2) = 0.04945332446081400;
        pointCoordIso(10, 2) = 0.06976240839115840;
        pointCoordIso(11, 2) = 0.08065399489035940;
        pointCoordIso(12, 2) = 0.08065399489035940;
        pointCoordIso(13, 2) = 0.06976240839115840;
        pointCoordIso(14, 2) = 0.04945332446081400;
        pointCoordIso(15, 2) = 0.02251130659095380;
        pointCoordIso(16, 2) = 0.03175606455054590;
        pointCoordIso(17, 2) = 0.06976240839115840;
        pointCoordIso(18, 2) = 0.09841185961908670;
        pointCoordIso(19, 2) = 0.11377631314509700;
        pointCoordIso(20, 2) = 0.11377631314509700;
        pointCoordIso(21, 2) = 0.09841185961908670;
        pointCoordIso(22, 2) = 0.06976240839115840;
        pointCoordIso(23, 2) = 0.03175606455054590;
        pointCoordIso(24, 2) = 0.03671394848693700;
        pointCoordIso(25, 2) = 0.08065399489035940;
        pointCoordIso(26, 2) = 0.11377631314509700;
        pointCoordIso(27, 2) = 0.13153952666880100;
        pointCoordIso(28, 2) = 0.13153952666880100;
        pointCoordIso(29, 2) = 0.11377631314509700;
        pointCoordIso(30, 2) = 0.08065399489035940;
        pointCoordIso(31, 2) = 0.03671394848693700;
        pointCoordIso(32, 2) = 0.03671394848693700;
        pointCoordIso(33, 2) = 0.08065399489035940;
        pointCoordIso(34, 2) = 0.11377631314509700;
        pointCoordIso(35, 2) = 0.13153952666880100;
        pointCoordIso(36, 2) = 0.13153952666880100;
        pointCoordIso(37, 2) = 0.11377631314509700;
        pointCoordIso(38, 2) = 0.08065399489035940;
        pointCoordIso(39, 2) = 0.03671394848693700;
        pointCoordIso(40, 2) = 0.03175606455054590;
        pointCoordIso(41, 2) = 0.06976240839115840;
        pointCoordIso(42, 2) = 0.09841185961908670;
        pointCoordIso(43, 2) = 0.11377631314509700;
        pointCoordIso(44, 2) = 0.11377631314509700;
        pointCoordIso(45, 2) = 0.09841185961908670;
        pointCoordIso(46, 2) = 0.06976240839115840;
        pointCoordIso(47, 2) = 0.03175606455054590;
        pointCoordIso(48, 2) = 0.02251130659095380;
        pointCoordIso(49, 2) = 0.04945332446081400;
        pointCoordIso(50, 2) = 0.06976240839115840;
        pointCoordIso(51, 2) = 0.08065399489035940;
        pointCoordIso(52, 2) = 0.08065399489035940;
        pointCoordIso(53, 2) = 0.06976240839115840;
        pointCoordIso(54, 2) = 0.04945332446081400;
        pointCoordIso(55, 2) = 0.02251130659095380;
        pointCoordIso(56, 2) = 0.01024721654119470;
        pointCoordIso(57, 2) = 0.02251130659095380;
        pointCoordIso(58, 2) = 0.03175606455054590;
        pointCoordIso(59, 2) = 0.03671394848693700;
        pointCoordIso(60, 2) = 0.03671394848693700;
        pointCoordIso(61, 2) = 0.03175606455054590;
        pointCoordIso(62, 2) = 0.02251130659095380;
        pointCoordIso(63, 2) = 0.01024721654119470;
    }
    return pointCoordIso;
}

std::pair<vector<double>, matrix<double>> Cell::cellContributions(const std::string &ep, const std::string &typeAnalyze,
                                                                  const int &step, const int &numberOfStep, const double &deltat,
                                                                  const double &beta, const double &gama, const int &pointsQuadrature)
{
    int npc = controlPoints_.size();

    vector<double> wpc(npc);

    for (int i = 0; i < npc; i++)
    {
        wpc(i) = controlPoints_[i]->getWeight();
    }

    bounded_vector<int, 2> inc;
    inc = controlPoints_[npc - 1]->getINC();

    vector<double> rhs(2 * npc, 0.0);
    matrix<double> tangent(2 * npc, 2 * npc, 0.0);
    matrix<double> domainIntegrationPoints_ = isoQuadrature(pointsQuadrature);

    double young, poisson, density;
    patch_->getMaterial()->setProperties(young, poisson, density);
    if (typeAnalyze == "STATIC")
    {
        density = 0.0;
    }
    double thickness = patch_->getThickness();
    if (distanceFE_.size() == 0)
    {
        for (int ih = 0; ih < domainIntegrationPoints_.size1(); ih++)
        {
            double xsi1 = domainIntegrationPoints_(ih, 0);
            double xsi2 = domainIntegrationPoints_(ih, 1);
            double weight = domainIntegrationPoints_(ih, 2);
            bounded_vector<double, 2> qxsi;
            qxsi(0) = xsi1;
            qxsi(1) = xsi2;

            std::pair<vector<double>, matrix<double>> functions;

            functions = shapeFunctionAndDerivates(qxsi, wpc, inc);

            bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(functions.second); //initial configuration map
            double j0 = jacobianDeterminant(A0);
            bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
            A0I(0, 0) = A0(1, 1) / j0;
            A0I(1, 1) = A0(0, 0) / j0;
            A0I(0, 1) = -A0(0, 1) / j0;
            A0I(1, 0) = -A0(1, 0) / j0;

            bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(functions.second); //current configuration map
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
                    //double shape = 0.0;

                    //shape = phi(i) * shape;

                    if (typeAnalyze == "DYNAMIC")
                    {
                        for (int m = 0; m < controlPoints_.size(); m++)
                        {
                            accel += functions.first(m) * controlPoints_[m]->getCurrentAcceleration()(j);
                        }
                    }

                    double m = density * functions.first(i) * accel; //inertial force

                    double shape = (shapeForce_(j) * step / numberOfStep / thickness) * functions.first(i);

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
    }
    else
    {
        for (int ih = 0; ih < domainIntegrationPoints_.size1(); ih++)
        {
            if (distanceFE_[ih] <= 0.0)
            {
                double xsi1 = domainIntegrationPoints_(ih, 0);
                double xsi2 = domainIntegrationPoints_(ih, 1);
                double weight = domainIntegrationPoints_(ih, 2);
                bounded_vector<double, 2> qxsi;
                qxsi(0) = xsi1;
                qxsi(1) = xsi2;

                std::pair<vector<double>, matrix<double>> functions;

                functions = shapeFunctionAndDerivates(qxsi, wpc, inc);

                bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(functions.second); //initial configuration map
                double j0 = jacobianDeterminant(A0);
                bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
                A0I(0, 0) = A0(1, 1) / j0;
                A0I(1, 1) = A0(0, 0) / j0;
                A0I(0, 1) = -A0(0, 1) / j0;
                A0I(1, 0) = -A0(1, 0) / j0;

                bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(functions.second); //current configuration map
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
                        //double shape = 0.0;

                        //shape = phi(i) * shape;

                        if (typeAnalyze == "DYNAMIC")
                        {
                            for (int m = 0; m < controlPoints_.size(); m++)
                            {
                                accel += functions.first(m) * controlPoints_[m]->getCurrentAcceleration()(j);
                            }
                        }

                        double m = density * functions.first(i) * accel; //inertial force

                        double shape = (shapeForce_(j) * step / numberOfStep / thickness) * functions.first(i);

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
        }
    }

    return std::make_pair(rhs, tangent);
}

bounded_vector<double, 4> Cell::getCauchStress(const bounded_vector<double, 2> &qxsi, const std::string &ep) //(SIGMA11, SIGMA22, SIGMA33, SIGMA12)
{
    bounded_vector<double, 4> cauchStress;

    int npc = controlPoints_.size();
    vector<double> wpc(npc);
    bounded_vector<int, 2> inc;

    inc = controlPoints_[npc - 1]->getINC();

    for (int i = 0; i < npc; i++)
    {
        wpc(i) = controlPoints_[i]->getWeight();
    }

    double young, poisson, density;
    patch_->getMaterial()->setProperties(young, poisson, density);

    std::pair<vector<double>, matrix<double>> functions;
    functions = shapeFunctionAndDerivates(qxsi, wpc, inc); //retorna as funções de formas e derivadas calculadas nos pontos qxsi

    bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(functions.second); //initial configuration map
    double j0 = jacobianDeterminant(A0);
    bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
    A0I(0, 0) = A0(1, 1) / j0;
    A0I(1, 1) = A0(0, 0) / j0;
    A0I(0, 1) = -A0(0, 1) / j0;
    A0I(1, 0) = -A0(1, 0) / j0;

    bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(functions.second); //current configuration map
    bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                           //current deformation gradient
    identity_matrix<double> I(2);                                              //identity matrix
    bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);         //current green strain tensor
    bounded_matrix<double, 2, 2> S;                                            //second piola kirchhoff stress tensor

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

bounded_matrix<double, 2, 2> Cell::get2PiolaStress(const bounded_vector<double, 2> &qxsi, const std::string &ep)
{
    int npc = controlPoints_.size();
    vector<double> wpc(npc);
    bounded_vector<int, 2> inc;

    inc = controlPoints_[npc - 1]->getINC();

    for (int i = 0; i < npc; i++)
    {
        wpc(i) = controlPoints_[i]->getWeight();
    }

    double young, poisson, density;
    patch_->getMaterial()->setProperties(young, poisson, density);

    std::pair<vector<double>, matrix<double>> functions;
    functions = shapeFunctionAndDerivates(qxsi, wpc, inc); //retorna as funções de formas e derivadas calculadas nos pontos qxsi

    bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(functions.second); //initial configuration map
    double j0 = jacobianDeterminant(A0);
    bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
    A0I(0, 0) = A0(1, 1) / j0;
    A0I(1, 1) = A0(0, 0) / j0;
    A0I(0, 1) = -A0(0, 1) / j0;
    A0I(1, 0) = -A0(1, 0) / j0;

    bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(functions.second); //current configuration map
    bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                           //current deformation gradient
    identity_matrix<double> I(2);                                              //identity matrix
    bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);         //current green strain tensor
    bounded_matrix<double, 2, 2> S;                                            //second piola kirchhoff stress tensor

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
    return S;
}

bounded_vector<double, 4> Cell::getGreen(const bounded_vector<double, 2> &qxsi, const std::string &ep) //(SIGMA11, SIGMA22, SIGMA33, SIGMA12)
{
    bounded_vector<double, 4> green;

    int npc = controlPoints_.size();
    vector<double> wpc(npc);
    bounded_vector<int, 2> inc;

    inc = controlPoints_[npc - 1]->getINC();

    for (int i = 0; i < npc; i++)
    {
        wpc(i) = controlPoints_[i]->getWeight();
    }

    double young, poisson, density;
    patch_->getMaterial()->setProperties(young, poisson, density);

    std::pair<vector<double>, matrix<double>> functions;
    functions = shapeFunctionAndDerivates(qxsi, wpc, inc); //retorna as funções de formas e derivadas calculadas nos pontos qxsi

    bounded_matrix<double, 2, 2> A0 = referenceJacobianMatrix(functions.second); //initial configuration map
    double j0 = jacobianDeterminant(A0);
    bounded_matrix<double, 2, 2> A0I; //inverse initial configuration map
    A0I(0, 0) = A0(1, 1) / j0;
    A0I(1, 1) = A0(0, 0) / j0;
    A0I(0, 1) = -A0(0, 1) / j0;
    A0I(1, 0) = -A0(1, 0) / j0;

    bounded_matrix<double, 2, 2> A1 = currentJacobianMatrix(functions.second); //current configuration map
    bounded_matrix<double, 2, 2> Ac = prod(A1, A0I);                           //current deformation gradient
    identity_matrix<double> I(2);                                              //identity matrix
    bounded_matrix<double, 2, 2> Ec = 0.5 * (prod(trans(Ac), Ac) - I);         //current green strain tensor

    // if (Ec(0, 0) >= 10.0e-14 or Ec(0, 0) <= -10.0e-14)
    // {
    //     green(0) = Ec(0, 0);
    // }
    // else
    // {
    //     green(0) = 0.0;
    // }

    // if (Ec(1, 1) >= 10.0e-14 or Ec(1, 1) <= -10.0e-14)
    // {
    //     green(1) = Ec(1, 1);
    // }
    // else
    // {
    //     green(1) = 0.0;
    // }

    // if (Ec(0, 1) >= 10.0e-14 or Ec(0, 1) <= -10.0e-14)
    // {
    //     green(3) = Ec(0, 1);
    // }
    // else
    // {
    //     green(3) = 0.0;
    // }

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

matrix<double> Cell::massMatrix(const int &pointsQuadrature)
{
    int number = controlPoints_.size();
    matrix<double> mass(2 * number, 2 * number, 0.0);
    matrix<double> domainIntegrationPoints_ = isoQuadrature(pointsQuadrature);
    double auxiliar = patch_->getMaterial()->getDensity() * patch_->getThickness();
    double resul;
    bounded_vector<int, 2> INC_;
    vector<double> wpc2(number);

    INC_ = controlPoints_[number - 1]->getINC();
    for (int i = 0; i < number; i++)
    {
        wpc2(i) = controlPoints_[i]->getWeight();
    }

    for (int ih = 0; ih < domainIntegrationPoints_.size1(); ih++)
    {
        if (distanceFE_[ih] <= 0.0)
        {
            bounded_vector<double, 2> qxsi;
            qxsi(0) = domainIntegrationPoints_(ih, 0);
            qxsi(1) = domainIntegrationPoints_(ih, 1);
            double weight = domainIntegrationPoints_(ih, 2);

            vector<double> phi = shapeFunction(qxsi, wpc2, INC_);

            for (int i = 0; i < controlPoints_.size(); i++)
            {
                for (int k = 0; k < controlPoints_.size(); k++)
                {
                    resul = auxiliar * phi(i) * phi(k) * weight;
                    mass(2 * i, 2 * k) += resul;
                    mass(2 * i + 1, 2 * k + 1) += resul;
                }
            }
        }
    }
    return mass;
}

std::pair<vector<double>, vector<double>> Cell::boundaryShapeFunctionAndDerivates(const double &qxsi,
                                                                                  const vector<double> &wpc,
                                                                                  const bounded_vector<double, 2> inc,
                                                                                  const int &curveNumber)
{
    int deg, npc, dim, numLocal, ind, auxiliar;
    bounded_vector<int, 2> spans = patch_->getSpanNumber();
    if (curveNumber == 0 or curveNumber == 2)
    {
        deg = patch_->getDegree(0);
        npc = patch_->getNpc_Dir(0);
        auxiliar = spans(0);
    }
    else
    {
        deg = patch_->getDegree(1);
        npc = patch_->getNpc_Dir(1);
        auxiliar = spans(1);
    }

    dim = deg + npc + 1;
    numLocal = deg + 1;

    vector<double> knot(dim);
    vector<double> left(deg + 1);
    vector<double> right(deg + 1);
    vector<double> phiL(deg + 1);
    matrix<double> aux(2, deg);
    vector<double> dphiL(deg + 1);
    vector<double> phiIso(numLocal);
    vector<double> dphiIso(numLocal);

    double saved, temp;

    if (curveNumber == 0 or curveNumber == 2)
    {
        knot = patch_->getKnotVectorU();
        ind = inc(0);
    }
    else
    {
        knot = patch_->getKnotVectorV();
        ind = inc(1);
    }

    double uv1 = knot(ind);
    double uv2 = knot(ind + 1);

    const double xsi = ((uv2 - uv1) * qxsi + (uv2 + uv1)) * 0.5;
    matrix<double> BF(deg + 1, deg + 1);

    BF(0, 0) = 1.0;

    for (int j = 1; j <= deg; j++)
    {
        left(j) = xsi - knot(ind + 1 - j);
        right(j) = knot(ind + j) - xsi;
        saved = 0.0;

        for (int r = 0; r < j; r++)
        {
            //lower triangle (Piegl and Tiller pg. 68-71)
            BF(j, r) = right(r + 1) + left(j - r);
            temp = BF(r, j - 1) / BF(j, r);

            //upper triangle
            BF(r, j) = saved + right(r + 1) * temp;
            saved = left(j - r) * temp;
        }

        BF(j, j) = saved;
    }

    for (int i = 0; i <= deg; i++)
    {
        phiL(i) = BF(i, deg);
    }

    double sum = 0.0;
    for (int i = 0; i < numLocal; i++)
    {
        sum += phiL(i) * wpc(i);
    }

    for (int i = 0; i < numLocal; i++)
    {
        phiIso(i) = phiL(i) / sum;
    }

    //DERIVATES

    int s1, s2, k, rk, pk, j1, j2;
    //int , cor;
    double d;

    for (int r = 0; r <= deg; r++)
    {
        s1 = 0;
        s2 = 1;
        aux(0, 0) = 1.0;
        k = 1;
        d = 0.0;
        rk = r - k;
        pk = deg - k;

        if (r >= k)
        {
            aux(s2, 0) = aux(s1, 0) / BF(pk + 1, rk);
            d = aux(s2, 0) * BF(rk, pk);
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
            j2 = deg - r;
        }

        for (int j = j1; j <= j2; j++)
        {
            aux(s2, j) = (aux(s1, j) - aux(s1, j - 1)) / BF(pk + 1, rk + j);
            d = d + aux(s2, j) * BF(rk + j, pk);
        }

        if (r <= pk)
        {
            aux(s2, k) = -aux(s1, k - 1) / BF(pk + 1, r);
            d = d + aux(s2, k) * BF(r, pk);
        }

        dphiL(r) = d;

        int j = s1;
        s1 = s2;
        s2 = j;
    }

    for (int i = 0; i < (deg + 1); i++)
    {
        dphiL(i) = dphiL(i) * deg;
    }

    sum = 0.0;
    double sum1 = 0.0;

    for (int i = 0; i < numLocal; i++)
    {
        sum += phiL(i) * wpc(i);
        sum1 += dphiL(i) * wpc(i);
    }

    for (int i = 0; i < numLocal; i++)
    {
        dphiIso(i) = (1.0 / (sum * sum * 2 * auxiliar)) * (dphiL(i) * sum - phiL(i) * sum1);
    }

    return std::make_pair(phiIso, dphiIso);
};

matrix<double> Cell::boundaryIsoQuadrature(const int &points)
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

    return pointCoordIso;
}

vector<double> Cell::computeDistribuitedLoads(const bounded_vector<double, 2> &value, const int &quadraturePoints, const int &curveNumber)
{
    std::vector<ControlPoint *> points = getControlPointsOnSide(curveNumber);
    int npc = points.size();
    vector<double> distribuitedLoad(2 * npc, 0.0);
    matrix<double> integrationPoints(quadraturePoints, 2);
    integrationPoints = boundaryIsoQuadrature(quadraturePoints);
    vector<double> wpc(npc);
    bounded_vector<int, 2> inc_;
    double thickness = patch_->getThickness();

    inc_ = points[npc - 1]->getINC();

    for (int i = 0; i < npc; i++)
    {
        wpc(i) = points[i]->getWeight();
    }

    for (int iq = 0; iq < quadraturePoints; iq++)
    {
        double xsi = integrationPoints(iq, 0);
        double weight = integrationPoints(iq, 1);

        std::pair<vector<double>, vector<double>> functions = boundaryShapeFunctionAndDerivates(xsi, wpc, inc_, curveNumber);
        bounded_vector<double, 2> tangent;
        tangent(0) = 0.0;
        tangent(1) = 0.0;
        //double sum = 0.0;
        for (int ih = 0; ih < npc; ih++)
        {
            tangent(0) += functions.second(ih) * points[ih]->getCurrentCoordinate()(0);
            tangent(1) += functions.second(ih) * points[ih]->getCurrentCoordinate()(1);
            //sum += functions.first(ih);
        }

        double jacobian = sqrt(pow(tangent(0), 2) + pow(tangent(1), 2));
        // bounded_vector<double, 2> normal;
        // normal(0) = tangent(1) / jacobian;
        // normal(1) = -tangent(0) / jacobian;

        for (int ih = 0; ih < npc; ih++)
        {
            distribuitedLoad(2 * ih) += value(0) * functions.first(ih) * weight * thickness * jacobian;
            distribuitedLoad(2 * ih + 1) += value(1) * functions.first(ih) * weight * thickness * jacobian;
        }
    }
    return distribuitedLoad;
}

std::vector<ControlPoint *> Cell::getControlPointsOnSide(const int &side)
{
    std::vector<ControlPoint *> points;
    if (side == 0)
    {
        int numberOfPoints = patch_->getDegree(0) + 1;
        for (int i = 0; i < numberOfPoints; i++)
        {
            points.push_back(controlPoints_[i]);
        }
    }
    else if (side == 1)
    {
        int numberOfPoints = patch_->getDegree(1) + 1;
        int aux = patch_->getDegree(0) + 1;
        int j = aux - 1;

        for (int i = 0; i < numberOfPoints; i++)
        {
            points.push_back(controlPoints_[j]);
            j = j + aux;
        }
    }
    else if (side == 2)
    {
        int numberOfPoints = patch_->getDegree(0) + 1;
        int aux = patch_->getDegree(1);
        int j = aux * numberOfPoints;

        for (int i = 0; i < numberOfPoints; i++)
        {
            points.push_back(controlPoints_[j]);
            j = j + 1;
        }
    }
    else if (side == 3)
    {
        int numberOfPoints = patch_->getDegree(1) + 1;
        int aux = patch_->getDegree(0) + 1;
        int j = 0;

        for (int i = 0; i < numberOfPoints; i++)
        {
            points.push_back(controlPoints_[j]);
            j = j + aux;
        }
    }
    return points;
}

// void Cell::computeDistanceFromFEBoundary(const int &pointsQuadrature, std::vector<BoundaryElement *> boundaryFE)
// {
//     int npc = controlPoints_.size();
//     matrix<double> domainIntegrationPoints_ = isoQuadrature(pointsQuadrature);

//     if (distanceFE_.size() > 1)
//     {
//         distanceFE_.erase(distanceFE_.begin(), distanceFE_.begin() + distanceFE_.size());
//     }

//     distanceFE_.reserve(domainIntegrationPoints_.size1());

//     for (int ip = 0; ip < domainIntegrationPoints_.size1(); ip++)
//     {
//         distanceFE_[ip] = 100000000000.0;
//         bounded_vector<double, 2> xsi, coordIP;

//         xsi(0) = domainIntegrationPoints_(ip, 0);
//         xsi(1) = domainIntegrationPoints_(ip, 1);
//         coordIP = calculateGlobalCoordinate(xsi);

//         double distance;

//         std::cout << "COORD: " << coordIP(0) << " " << coordIP(1) << std::endl;
//         for (BoundaryElement *bound : boundaryFE)
//         {
//             double xsiBoundary = 0.0; //primeira tentativa
//             double deltaxsi = 1000.0;
//             std::vector<Node *> boundaryNodes = bound->getNodes();

//             int cont = 0;
//             bounded_vector<double, 2> coordBoundary;
//             bounded_vector<double, 2> firstDerivate;
//             bounded_vector<double, 2> secondDerivate;
//             bounded_vector<double, 2> normal;
//             double aux;

//             while (fabs(deltaxsi) >= 1.0e-06 and xsiBoundary >= -1.0 and xsiBoundary <= 1.0 and cont <= 15)
//             {
//                 matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(xsiBoundary); //PHI, PHI', PHI''
//                 coordBoundary(0) = 0.0;
//                 coordBoundary(1) = 0.0;
//                 firstDerivate(0) = 0.0;
//                 firstDerivate(1) = 0.0;
//                 secondDerivate(0) = 0.0;
//                 secondDerivate(1) = 0.0;
//                 int aux = 0;
//                 bounded_vector<double, 2> coordinateNode;
//                 for (Node *node : boundaryNodes)
//                 {
//                     coordinateNode = node->getCurrentCoordinate();

//                     coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
//                     coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

//                     firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
//                     firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

//                     secondDerivate(0) += boundaryFunctions(aux, 2) * coordinateNode(0);
//                     secondDerivate(1) += boundaryFunctions(aux, 2) * coordinateNode(1);
//                     aux = aux + 1;
//                 }
//                 double h = -((coordIP(0) - coordBoundary(0)) * (-firstDerivate(0)) + (coordIP(1) - coordBoundary(1)) * (-firstDerivate(1)));
//                 double dh_dxsi = (-firstDerivate(0)) * (-firstDerivate(0)) + (-secondDerivate(0)) * (coordIP(0) - coordBoundary(0)) +
//                                  (-firstDerivate(1)) * (-firstDerivate(1)) + (-secondDerivate(1)) * (coordIP(1) - coordBoundary(1));

//                 deltaxsi = h / dh_dxsi;

//                 xsiBoundary = xsiBoundary + deltaxsi;

//                 cont++;
//             }

//             if (xsiBoundary <= -1.0)
//             {
//                 matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(-1.0); //PHI, PHI', PHI''
//                 int aux = 0;
//                 bounded_vector<double, 2> coordinateNode;
//                 coordBoundary(0) = 0.0;
//                 coordBoundary(1) = 0.0;
//                 firstDerivate(0) = 0.0;
//                 firstDerivate(1) = 0.0;
//                 for (Node *node : boundaryNodes)
//                 {
//                     coordinateNode = node->getCurrentCoordinate();
//                     coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
//                     coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

//                     firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
//                     firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);
//                     aux = aux + 1;
//                 }
//             }
//             else if (xsiBoundary >= 1.0)
//             {
//                 matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(1.0); //PHI, PHI', PHI''
//                 int aux = 0;
//                 bounded_vector<double, 2> coordinateNode;
//                 coordBoundary(0) = 0.0;
//                 coordBoundary(1) = 0.0;
//                 firstDerivate(0) = 0.0;
//                 firstDerivate(1) = 0.0;
//                 for (Node *node : boundaryNodes)
//                 {
//                     coordinateNode = node->getCurrentCoordinate();
//                     coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
//                     coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

//                     firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
//                     firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);
//                     aux = aux + 1;
//                 }
//             }

//             distance = sqrt((coordIP(0) - coordBoundary(0)) * (coordIP(0) - coordBoundary(0)) + (coordIP(1) - coordBoundary(1)) * (coordIP(1) - coordBoundary(1)));
//             normal(0) = firstDerivate(1);
//             normal(1) = -firstDerivate(0);
//             aux = (coordIP(0) - coordBoundary(0)) * normal(0) + (coordIP(1) - coordBoundary(1)) * normal(1);

//             if (fabs(distance) < 1.0e-10)
//             {
//                 distance = 0.0;
//             }

//             if (fabs(distance) < fabs(distanceFE_[ip]))
//             {
//                 if (aux > 0)
//                 {
//                     distanceFE_[ip] = -distance;
//                 }
//                 else
//                 {
//                     distanceFE_[ip] = distance;
//                 }
//             }
//         }
//         std::cout << "DISTANCE: " << distanceFE_[ip] << std::endl;
//         std::cout << std::endl;
//     }

//     // std::cout<<"DISTANCE: ";
//     // for(int i = 0; i<16; i++)
//     // {
//     //     std::cout<<distanceFE_[i]<<" ";
//     // }
//     // std::cout<<std::endl;
//     std::cout << "FOI DE NOVO... " << std::endl;
// }

bounded_vector<double, 2> Cell::calculateGlobalCoordinate(const bounded_vector<double, 2> &qxsi)
{
    int npc = controlPoints_.size();
    bounded_vector<int, 2> inc_;
    vector<double> wpc2(npc);
    vector<double> phi;
    inc_ = controlPoints_[npc - 1]->getINC();
    for (int i = 0; i < npc; i++)
    {
        wpc2(i) = controlPoints_[i]->getWeight();
    }

    phi = shapeFunction(qxsi, wpc2, inc_);

    bounded_vector<double, 2> coordIP;

    coordIP(0) = 0.0;
    coordIP(1) = 0.0;

    for (int cp = 0; cp < npc; cp++) //global coordinates of the integration point
    {
        bounded_vector<double, 2> coordinateCP = controlPoints_[cp]->getCurrentCoordinate();
        coordIP(0) += phi(cp) * coordinateCP(0);
        coordIP(1) += phi(cp) * coordinateCP(1);
    }
    return coordIP;
}

void Cell::setDistanceFromFEBoundary(const std::vector<double> &distance)
{
    distanceFE_ = distance;
}

vector<double> Cell::diagonalMass(const int &points)
{
    int number = controlPoints_.size();
    vector<double> mass(number, 0.0);
    matrix<double> domainIntegrationPoints_ = isoQuadrature(points);
    bounded_vector<int, 2> INC_;
    vector<double> wpc2(number);

    INC_ = controlPoints_[number - 1]->getINC();

    for (int i = 0; i < number; i++)
    {
        wpc2(i) = controlPoints_[i]->getWeight();
    }

    for (int ih = 0; ih < domainIntegrationPoints_.size1(); ih++)
    {
        if (distanceFE_[ih] <= 0.0)
        {
            bounded_vector<double, 2> qxsi;
            qxsi(0) = domainIntegrationPoints_(ih, 0);
            qxsi(1) = domainIntegrationPoints_(ih, 1);
            double weight = domainIntegrationPoints_(ih, 2);

            vector<double> phi = shapeFunction(qxsi, wpc2, INC_);

            for (int i = 0; i < number; i++)
            {
                mass(i) += phi(i) * phi(i) * weight;
            }
        }
    }
    return mass;
}
