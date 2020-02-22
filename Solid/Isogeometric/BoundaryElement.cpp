#include "BoundaryElement.h"

std::vector<ControlPoint *> BoundaryElement::getControlPoints()
{
    return controlPoints_;
}

ControlPoint *BoundaryElement::getControlPoint(const int &index)
{
    return controlPoints_[index];
}

std::pair<vector<double>, vector<double>> BoundaryElement::shapeFunctionAndDerivates(const double &qxsi,
                                                                                     const vector<double> &wpc,
                                                                                     const bounded_vector<double, 2> inc)
{
    int deg, npc, dim, numLocal, ind;
    if (cunverNumber_ == 0 or cunverNumber_ == 2)
    {
        deg = patch_->getDegree(0);
        npc = patch_->getNpc_Dir(0);
    }
    else
    {
        deg = patch_->getDegree(1);
        npc = patch_->getNpc_Dir(1);
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

    if (cunverNumber_ == 0 or cunverNumber_ == 2)
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
        phiIso(i) = phiL(i) * wpc(i) / sum;
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
        dphiIso(i) = (wpc(i) * dphiL(i) * sum - wpc(i) * phiL(i) * sum1) / (sum * sum);
    }

    return std::make_pair(phiIso, dphiIso);
};

matrix<double> BoundaryElement::isoQuadrature(const int &points)
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
    int npc = controlPoints_.size();
    vector<double> distribuitedLoad(2 * npc, 0.0);
    matrix<double> integrationPoints(quadraturePoints, 2);
    integrationPoints = isoQuadrature(quadraturePoints);
    vector<double> wpc(npc);
    bounded_vector<int, 2> inc_;
    double thickness = patch_->getThickness();

    inc_ = controlPoints_[npc - 1]->getINC();

    for (int i = 0; i < npc; i++)
    {
        wpc(i) = controlPoints_[i]->getWeight();
    }

    for (int iq = 0; iq < quadraturePoints; iq++)
    {
        double xsi = integrationPoints(iq, 0);
        double weight = integrationPoints(iq, 1);

        std::pair<vector<double>, vector<double>> functions = shapeFunctionAndDerivates(xsi, wpc, inc_);
        bounded_vector<double, 2> tangent; // normal;

        for (int ih = 0; ih < npc; ih++)
        {
            tangent(0) += functions.second(ih) * controlPoints_[ih]->getCurrentCoordinate()(0);
            tangent(1) += functions.second(ih) * controlPoints_[ih]->getCurrentCoordinate()(1);
        }

        // double j = norm_2(tangent);

        // normal(0) = tangent(1) / j;
        // normal(1) = -tangent(0) / j;

        double aux = 0.0;
        for (int m = 0; m < npc; m++)
        {
            aux += functions.first(m);
        }

        for (int ih = 0; ih < npc; ih++)
        {
            distribuitedLoad(2 * ih) += aux * value(0) * functions.first(ih) * tangent(1) * weight * thickness;
            distribuitedLoad(2 * ih + 1) += aux * value(1) * functions.first(ih) * (-1.0*tangent(0)) * weight * thickness;
        }
    }
    return distribuitedLoad;
}
