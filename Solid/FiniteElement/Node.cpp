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
        stressState_(i) = 0.0;
        stressState_(i + 2) = 0.0;
    }
    stressState_(3) = 1.0;
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

bounded_vector<double, 4> Node::getStressState()
{
    return stressState_;
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

void Node::setStressState(const bounded_vector<double, 3> &stressState)
{
    stressState_(0) = stressState_(0) + stressState(0);
    stressState_(1) = stressState_(1) + stressState(1);
    stressState_(2) = stressState_(2) + stressState(2);
    stressState_(3) = stressState_(3) + 1.0;
}

void Node::setZeroStressState()
{
    stressState_(0) = 0.0;
    stressState_(1) = 0.0;
    stressState_(2) = 0.0;
    stressState_(3) = 0.0;
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

// void Node::setIndexFE(const int &index)
// {
//     indexFE_ = index;
// }
