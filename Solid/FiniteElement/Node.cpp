#include "Node.h"

Node::Node() {}

Node::Node(const int &index, const int &indexFE, const bounded_vector<double, 2> &initialCoordinate)
{
    index_ = index;
    indexFE_ = indexFE;
    initialCoordinate_ = initialCoordinate;
    pastCoordinate_ = initialCoordinate;
    currentCoordinate_ = initialCoordinate;
    for (size_t i = 0; i < 2; i++)
    {
        currentVelocity_(i) = 0.0;
        currentAcceleration_(i) = 0.0;
        pastVelocity_(i) = 0.0;
        pastAcceleration_(i) = 0.0;
    }
    stressState_(0) = 0.0;
    stressState_(1) = 0.0;
    stressState_(2) = 0.0;
    stressState_(3) = 0.0;
    stressState_(4) = 1.0;
    distanceToBoundary_ = -10000.0;
    cellIndex_ = -1;
}

Node::~Node() {}

int Node::getIndex()
{
    return index_;
}

int Node::getIndexFE()
{
    return indexFE_;
}

bounded_vector<double, 2> Node::getInitialCoordinate()
{
    return initialCoordinate_;
}

bounded_vector<double, 2> Node::getPastCoordinate()
{
    return pastCoordinate_;
}

bounded_vector<double, 2> Node::getPastVelocity()
{
    return pastVelocity_;
}

bounded_vector<double, 2> Node::getPastAcceleration()
{
    return pastAcceleration_;
}

bounded_vector<double, 2> Node::getCurrentCoordinate()
{
    return currentCoordinate_;
}

bounded_vector<double, 2> Node::getCurrentVelocity()
{
    return currentVelocity_;
}

bounded_vector<double, 2> Node::getCurrentAcceleration()
{
    return currentAcceleration_;
}

bounded_vector<double, 2> Node::getCurrentDisplacement()
{
    return (currentCoordinate_ - initialCoordinate_);
}

bounded_vector<double, 5> Node::getStressState()
{
    return stressState_;
}

double Node::getDistanceToBoundary()
{
    return distanceToBoundary_;
}

bool Node::getPoint()
{
    return point_;
}

int Node::getCellIndex()
{
    return cellIndex_;
}

bounded_vector<double, 2> Node::getXsisGlobal()
{
    return xsisglobal_;
}

bounded_vector<double, 2> Node::getValuesOfBlendingFunction()
{
    return blendingFunction_;
}

void Node::setPastCoordinate(const bounded_vector<double, 2> &pastCoordinate)
{
    pastCoordinate_ = pastCoordinate;
}

void Node::setPastVelocity(const bounded_vector<double, 2> &pastVelocity)
{
    pastVelocity_ = pastVelocity;
}

void Node::setPastAcceleration(const bounded_vector<double, 2> &pastAcceleration)
{
    pastAcceleration_ = pastAcceleration;
}

void Node::setCurrentCoordinate(const bounded_vector<double, 2> &currentCoordinate)
{
    currentCoordinate_ = currentCoordinate;
}

void Node::setCurrentVelocity(const bounded_vector<double, 2> &currentVelocity)
{
    currentVelocity_ = currentVelocity;
}

void Node::setCurrentAcceleration(const bounded_vector<double, 2> &currentAcceleration)
{
    currentAcceleration_ = currentAcceleration;
}

void Node::setInitialCoordinate(const bounded_vector<double, 2> &initialCoordinate)
{
    initialCoordinate_ = initialCoordinate;
    currentCoordinate_ = initialCoordinate;
}

void Node::setStressState(const bounded_vector<double, 4> &stressState)
{
    stressState_(0) = stressState_(0) + stressState(0);
    stressState_(1) = stressState_(1) + stressState(1);
    stressState_(2) = stressState_(2) + stressState(2);
    stressState_(3) = stressState_(3) + stressState(3);
    stressState_(4) = stressState_(4) + 1.0;
}

void Node::setZeroStressState()
{
    stressState_(0) = 0.0;
    stressState_(1) = 0.0;
    stressState_(2) = 0.0;
    stressState_(3) = 0.0;
    stressState_(4) = 0.0;
}

void Node::incrementCurrentCoordinate(const int &direction, const double &value)
{
    currentCoordinate_(direction) += value;
}

void Node::setIndex(const int &index)
{
    index_ = index;
}

void Node::updatePastValue()
{
    pastAcceleration_ = currentAcceleration_;
    pastCoordinate_ = currentCoordinate_;
    pastVelocity_ = currentVelocity_;
}

void Node::setDistanceToBoundary(const double &distance)
{
    distanceToBoundary_ = distance;
}

void Node::setCellIndex(const int &cellIndex)
{
    cellIndex_ = cellIndex;
}

void Node::setXsisGlobal(const bounded_vector<double, 2> &xsis)
{
    xsisglobal_ = xsis;
}

void Node::setValuesOfBlendingFunction(const bounded_vector<double, 2> &bvalues)
{
    blendingFunction_ = bvalues;
}

void Node::setPoint()
{
    point_ = true;
}
