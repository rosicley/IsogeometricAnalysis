#include "Mesh.h"

Mesh::Mesh() {}

Mesh::Mesh(const int &index, Material *mat, const double &thickness, const std::string &elementType)
{
    index_ = index;
    material_ = mat;
    thickness_ = thickness;
    elementType_ = elementType;
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
}

Mesh::~Mesh() {}

double Mesh::getThickness()
{
    return thickness_;
}

int Mesh::getOrder()
{
    return order_;
}

Material *Mesh::getMaterial()
{
    return material_;
}

int Mesh::getIndex()
{
    return index_;
}

std::string Mesh::getElementType()
{
    return elementType_;
}

