#include "ControlPoint.h"

ControlPoint::ControlPoint() {}

ControlPoint::ControlPoint(const int &index,
                           const bounded_vector<double, 2> &initialCoordinate,
                           const double &weight)
{
    index_ = index;
    initialCoordinate_ = initialCoordinate;
    pastCoordinate_ = initialCoordinate;
    currentCoordinate_ = initialCoordinate;
    weight_ = weight;
    for (size_t i = 0; i < 2; i++)
    {
        currentVelocity_(i) = 0.0;
        currentAcceleration_(i) = 0.0;
        pastVelocity_(i) = 0.0;
        pastAcceleration_(i) = 0.0;
        // stressState_(i) = 0.0;
        // stressState_(i + 2) = 0.0;
    }
    // stressState_(2)=1.0;
}

ControlPoint::~ControlPoint() {}

int ControlPoint::getIndex()
{
    return index_;
}

double ControlPoint::getWeight()
{
    return weight_;
}

bounded_vector<double, 2> ControlPoint::getInitialCoordinate()
{
    return initialCoordinate_;
}

bounded_vector<double, 2> ControlPoint::getCurrentDisplacement()
{
    return (currentCoordinate_ - initialCoordinate_);
}

bounded_vector<double, 2> ControlPoint::getPastCoordinate()
{
    return pastCoordinate_;
}

bounded_vector<double, 2> ControlPoint::getPastVelocity()
{
    return pastVelocity_;
}

bounded_vector<double, 2> ControlPoint::getPastAcceleration()
{
    return pastAcceleration_;
}

bounded_vector<double, 2> ControlPoint::getCurrentCoordinate()
{
    return currentCoordinate_;
}

bounded_vector<double, 2> ControlPoint::getCurrentVelocity()
{
    return currentVelocity_;
}

bounded_vector<double, 2> ControlPoint::getCurrentAcceleration()
{
    return currentAcceleration_;
}

bounded_vector<int, 2> ControlPoint::getINC()
{
    return inc_;
}

// bounded_vector<double, 4> ControlPoint::getStressState()
// {
//     return stressState_;
// }

void ControlPoint::setPastCoordinate(const bounded_vector<double, 2> &pastCoordinate)
{
    pastCoordinate_ = pastCoordinate;
}

void ControlPoint::setPastVelocity(const bounded_vector<double, 2> &pastVelocity)
{
    pastVelocity_ = pastVelocity;
}

void ControlPoint::setPastAcceleration(const bounded_vector<double, 2> &pastAcceleration)
{
    pastAcceleration_ = pastAcceleration;
}

void ControlPoint::setCurrentCoordinate(const bounded_vector<double, 2> &currentCoordinate)
{
    currentCoordinate_ = currentCoordinate;
}

void ControlPoint::setCurrentVelocity(const bounded_vector<double, 2> &currentVelocity)
{
    currentVelocity_ = currentVelocity;
}

void ControlPoint::setCurrentAcceleration(const bounded_vector<double, 2> &currentAcceleration)
{
    currentAcceleration_ = currentAcceleration;
}

void ControlPoint::incrementCurrentCoordinate(const int &direction, const double &value)
{
    currentCoordinate_(direction) += value;
}

void ControlPoint::setINC(const bounded_vector<int, 2> inc)
{
    inc_ = inc;
}

void ControlPoint::setIndex(const int &newIndex)
{
    index_ = newIndex;
}


void ControlPoint::updatePastValue()
{
    pastAcceleration_ = currentAcceleration_;
    pastCoordinate_ = currentCoordinate_;
    pastVelocity_ = currentVelocity_;
}