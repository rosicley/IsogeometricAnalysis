#include "Point.h"

Point::Point() {}

Point::Point(const std::string &name,
			 const bounded_vector<double, 2> &coordinates,
			 const double &lcar,
			 const bool &discretization)
{
	name_ = name;
	coordinates_ = coordinates;
	lcar_ = lcar;
	discretization_ = discretization;
}
Point::~Point() {}

std::string Point::getName()
{
	return name_;
}

std::string Point::getGmshCode()
{
	std::stringstream text;
	text << std::fixed;
	if (discretization_)
	{
		text << name_ << " = newp; Point(" << name_ << ") = {" << coordinates_[0] << ", " << coordinates_[1] << ", 0.0, " << lcar_
			 << "}; Physical Point('" << name_ << "') = {" << name_ << "};\n//\n";
	}
	else
	{
		text << name_ << " = newp; Point(" << name_ << ") = {" << coordinates_[0] << ", " << coordinates_[1] << ", 0.0, " << lcar_ << "};\n//\n";
	}
	return text.str();
}

void Point::addNodeToPoint(Node *node)
{
	pointNode_ = node;
	pointNode_->setPoint();
}

Node *Point::getPointNode()
{
	return pointNode_;
}

bounded_vector<double, 2> Point::getCoordinates()
{
	return coordinates_;
}

double Point::getlcar()
{
	return lcar_;
}

void Point::setlcar(const double &newlcar)
{
	lcar_ = newlcar;
}

// void Point::setCrackPoint()
// {
// 	crackPoint_ = true;
// }

// bool Point::getCrackPoint()
// {
// 	return crackPoint_;
// }

void Point::setCoordinates(const bounded_vector<double, 2> &newCoordinates)
{
	coordinates_ = newCoordinates;
}
