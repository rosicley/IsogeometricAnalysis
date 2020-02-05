#pragma once
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

class Material
{
public:
    Material();

    Material(const int &index,
             const double &young,
             const double &poisson,
             const double &density);

    ~Material();

    int getIndex();

    double getYoung();

    double getPoisson();

    double getDensity();

    void setProperties(double &young, double &poisson, double &density);

private:
    int index_;

    double young_;

    double poisson_;

    double density_;
};
