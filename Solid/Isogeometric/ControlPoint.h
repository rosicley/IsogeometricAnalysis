#pragma once

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>

using namespace boost::numeric::ublas;

class ControlPoint
{
public:
	ControlPoint();

	ControlPoint(const int &index,
				 const bounded_vector<double, 2> &initialCoordinate,
				 const double &weight);

	~ControlPoint();

	int getIndex();

	double getWeight();

	bounded_vector<double, 2> getInitialCoordinate();

	bounded_vector<double, 2> getPastCoordinate();

	bounded_vector<double, 2> getPastVelocity();

	bounded_vector<double, 2> getPastAcceleration();

	bounded_vector<double, 2> getCurrentCoordinate();

	bounded_vector<double, 2> getCurrentVelocity();

	bounded_vector<double, 2> getCurrentAcceleration();

	bounded_vector<int, 2> getINC();

	bounded_vector<double, 2> getCurrentDisplacement();

	// bounded_vector<double, 4> getStressState();2

	void setPastCoordinate(const bounded_vector<double, 2> &pastCoordinate);

	void setPastVelocity(const bounded_vector<double, 2> &pastVelocity);

	void setPastAcceleration(const bounded_vector<double, 2> &pastAcceleration);

	void setCurrentCoordinate(const bounded_vector<double, 2> &currentCoordinate);

	void setCurrentVelocity(const bounded_vector<double, 2> &currentVelocity);

	void setCurrentAcceleration(const bounded_vector<double, 2> &currentAcceleration);

	void incrementCurrentCoordinate(const int &direction, const double &value);

	void setINC(const bounded_vector<int, 2> inc);
	// void setStressState(const bounded_vector<double, 3> &stressState);

	// void setZeroStressState();

private:
	int index_;

	double weight_;

	bounded_vector<double, 2> initialCoordinate_;

	bounded_vector<double, 2> pastCoordinate_;

	bounded_vector<double, 2> pastVelocity_;

	bounded_vector<double, 2> pastAcceleration_;

	bounded_vector<double, 2> currentCoordinate_;

	bounded_vector<double, 2> currentVelocity_;

	bounded_vector<double, 2> currentAcceleration_;

	bounded_vector<int, 2> inc_;

	// bounded_vector<double, 4> stressState_; //{SigmaX1, SigmaX2, TalX1X2, contador}
};