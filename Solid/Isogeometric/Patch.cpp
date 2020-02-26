#include "Patch.h"

Patch::Patch() {}

Patch::Patch(const int &index, const int &npc, Material *mat, const double &thickness)

{
    index_ = index;
    npc_ = npc;
    material_ = mat;
    thickness_ = thickness;
}

Patch::~Patch() {}

int Patch::getNumberOfCells()
{
    return numberOfCells_;
}

int Patch::getControlPointsNumber()
{
    return npc_;
}

int Patch::getIndex()
{
    return index_;
}

vector<double> Patch::getKnotVectorU()
{
    return knotVectorU_;
}

Material *Patch::getMaterial()
{
    return material_;
}

vector<double> Patch::getKnotVectorV()
{
    return knotVectorV_;
}

ControlPoint *Patch::getControlPoint(int index)
{
    return controlPoints_[index];
}

std::vector<ControlPoint *> Patch::getControlPoints()
{
    return controlPoints_;
}

double Patch::getThickness()
{
    return thickness_;
}

bounded_vector<int, 2> Patch::getDegrees()
{
    return degree_;
}

void Patch::setSpanNumber(bounded_vector<int, 2> spans)
{
    spans_ = spans;
}

bounded_vector<int, 2> Patch::getSpanNumber()
{
    return spans_;
}

int Patch::getDegree(const int &dir)
{
    return degree_(dir);
}

int Patch::getNpc_Dir(const int &dir)
{
    return npc_dir_(dir);
}

void Patch::setDegree(const bounded_vector<int, 2> &degree)
{
    degree_ = degree;
}

void Patch::setNpc_Dir(const bounded_vector<int, 2> &npc_dir)
{
    npc_dir_ = npc_dir;
}

void Patch::setKnotU(const vector<double> &uknot)
{
    knotVectorU_ = uknot;
}

void Patch::setKnotV(const vector<double> &vknot)
{
    knotVectorV_ = vknot;
}

void Patch::addControlPoint(const int &index,
                            const bounded_vector<double, 2> &initialCoordinate,
                            const double &weight)
{
    ControlPoint *point = new ControlPoint(index, initialCoordinate, weight);
    controlPoints_.push_back(point);
}

void Patch::setKnotsVectors(const vector<double> &uknot, const vector<double> &vknot)
{
    knotVectorU_ = uknot;
    knotVectorV_ = vknot;
}

void Patch::removeControlPoints()
{
    controlPoints_.erase(controlPoints_.begin(),controlPoints_.begin()+controlPoints_.size());
}