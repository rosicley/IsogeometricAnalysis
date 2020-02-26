#pragma once
//#include "Cell.h"
#include "../Material.h"
#include "ControlPoint.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

class Patch
{
public:
    Patch();

    Patch(const int &index, const int &npc, Material *mat, const double &thickness);

    ~Patch();

    int getDegree(const int &dir);

    int getIndex();

    bounded_vector<int, 2> getDegrees();

    int getNpc_Dir(const int &dir);

    void setDegree(const bounded_vector<int, 2> &degree);

    void setNpc_Dir(const bounded_vector<int, 2> &npc_dir);

    int getNumberOfCells();

    int getControlPointsNumber(); //number of control points

    matrix<double> isoQuadrature();

    vector<double> getKnotVectorU();

    vector<double> getKnotVectorV();

    ControlPoint *getControlPoint(int index);

    std::vector<ControlPoint *> getControlPoints();

    void addControlPoint(const int &index,
                         const bounded_vector<double, 2> &initialCoordinate,
                         const double &weight);

    void setKnotU(const vector<double> &uknot);

    void setKnotV(const vector<double> &vknot);

    void setKnotsVectors(const vector<double> &uknot, const vector<double> &vknot);

    Material* getMaterial();

    double getThickness();

    void setSpanNumber(bounded_vector<int, 2> spans);

    bounded_vector<int, 2> getSpanNumber();

    void removeControlPoints();

private:
    int index_; //ok

    bounded_vector<int, 2> degree_;

    bounded_vector<int, 2> npc_dir_;

    bounded_vector<int, 2> spans_;

    int numberOfCells_; // = numberOfSpanU * numberOfSpanV;

    int npc_; //ok

    vector<double> knotVectorU_; //ok

    vector<double> knotVectorV_; //ok

    std::vector<ControlPoint *> controlPoints_; //ok

    Material *material_; //ok

    double thickness_; //ok
};