#pragma once
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>

using namespace boost::numeric::ublas;

class Node
{
public:
	Node();

	Node(const int &index, const int &indexFE,
		 const bounded_vector<double, 2> &initialCoordinate); //(nº do nó, {x1, x2})

	~Node();

	int getIndex();

	int getIndexFE();

	bounded_vector<double, 2> getInitialCoordinate();

	bounded_vector<double, 2> getPastCoordinate();

	bounded_vector<double, 2> getPastVelocity();

	bounded_vector<double, 2> getPastAcceleration();

	bounded_vector<double, 2> getCurrentCoordinate();

	bounded_vector<double, 2> getCurrentVelocity();

	bounded_vector<double, 2> getCurrentAcceleration();

	bounded_vector<double, 2> getCurrentDisplacement();

	bounded_vector<double, 5> getStressState();

	double getDistanceToBoundary();

	int getCellIndex();

	bounded_vector<double, 2> getXsisGlobal();

	bounded_vector<double, 2> getValuesOfBlendingFunction();

	void setPastCoordinate(const bounded_vector<double, 2> &pastCoordinate);

	void setPastVelocity(const bounded_vector<double, 2> &pastVelocity);

	void setPastAcceleration(const bounded_vector<double, 2> &pastAcceleration);

	void setCurrentCoordinate(const bounded_vector<double, 2> &currentCoordinate);

	void setCurrentVelocity(const bounded_vector<double, 2> &currentVelocity);

	void setCurrentAcceleration(const bounded_vector<double, 2> &currentAcceleration);

	void setStressState(const bounded_vector<double, 4> &stressState);

	void setZeroStressState();

	void incrementCurrentCoordinate(const int &direction, const double &value);

	void setIndex(const int &index);

	void updatePastValue();

	void setDistanceToBoundary(const double &distance);

	void setCellIndex(const int &cellIndex);

	void setXsisGlobal(const bounded_vector<double, 2> &xsis);

	void setValuesOfBlendingFunction(const bounded_vector<double, 2> &bvalues);

private:
	int index_;

	int indexFE_;

	bounded_vector<double, 2> initialCoordinate_;

	bounded_vector<double, 2> pastCoordinate_;

	bounded_vector<double, 2> pastVelocity_;

	bounded_vector<double, 2> pastAcceleration_;

	bounded_vector<double, 2> currentCoordinate_;

	bounded_vector<double, 2> currentVelocity_;

	bounded_vector<double, 2> currentAcceleration_;

	bounded_vector<double, 5> stressState_; //{SigmaX1, SigmaX2, TalX1X2, contador}

	double distanceToBoundary_;

	int cellIndex_; 

	bounded_vector<double, 2> xsisglobal_;

	bounded_vector<double, 2> blendingFunction_;
	
};