#include "Material.h"

Material::Material() {}

Material::Material(const int &index,
                   const double &young,
                   const double &poisson,
                   const double &density)
{
    index_ = index;
    young_ = young;
    poisson_ = poisson;
    density_ = density;
}

Material::~Material() {}

int Material::getIndex()
{
    return index_;
}

double Material::getYoung()
{
    return young_;
}

double Material::getPoisson()
{
    return poisson_;
}

double Material::getDensity()
{
    return density_;
}

void Material::setProperties(double &young, double &poisson, double &density)
{
    young = young_;
    poisson = poisson_;
    density = density_;
}
