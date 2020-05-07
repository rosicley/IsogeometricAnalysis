#include "Point.h"

Point::Point() {}

Point::Point(const std::string &name,
			 const std::vector<double> &coordinates,
			 const double &lcar,
			 const bool &discretization)
{
	name_ = name;
	coordinates_.reserve(2);
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
		return text.str();
	}
	else
	{
		text << name_ << " = newp; Point(" << name_ << ") = {" << coordinates_[0] << ", " << coordinates_[1] << ", 0.0, " << lcar_ << "};\n//\n";
		return text.str();
	}
}

void Point::addNodeToPoint(Node *node)
{
	pointNode_ = node;
}

Node *Point::getPointNode()
{
	return pointNode_;
}