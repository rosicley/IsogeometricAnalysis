#include "Crack.h"

Crack::Crack() {}

Crack::Crack(const std::string &name, std::vector<Point *> points, PlaneSurface *surface, const std::string &openBoundary)
{
    name_ = name;
    points_ = points;
    surface_ = surface;
    openBoundary_ = openBoundary;
}

Crack::~Crack() {}

std::string Crack::getName()
{
    return name_;
}

std::vector<Point *> Crack::getPoints()
{
    return points_;
}

PlaneSurface *Crack::getPlaneSurface()
{
    return surface_;
}

Point *Crack::getLastPoint()
{
    return points_[points_.size() - 1];
}

Point *Crack::getFirstPoint()
{
    return points_[0];
}

std::string Crack::getGmshCodeCrackPlugin()
{
    std::stringstream text;
    text << "Plugin(Crack).Dimension = 1;\n";
    text << "Plugin(Crack).PhysicalGroup = " << name_ << ";\n";
    if (openBoundary_ == "first")
    {
        text << "Plugin(Crack).OpenBoundaryPhysicalGroup = " << points_[0]->getName() << ";\n";
    }
    else if (openBoundary_ == "second")
    {
        text << "Plugin(Crack).OpenBoundaryPhysicalGroup = " << points_[points_.size() - 1]->getName() << ";\n";
    }
    text << "Plugin(Crack).Run; \n//\n";

    return text.str();
}

// std::string Crack::getGmshCode()
// {
// }

// std::string Crack::getGmshCodeEmbedded()
// {
//     return embeddedCode_;
// }

void Crack::addLastPoint(Point *newPoint)
{
    points_.push_back(newPoint);
    //points_[points_.size() - 2]->setlcar(points_[points_.size() - 3]->getlcar());
    points_[points_.size() - 2]->setlcar(1.0);
}

double Crack::getLastAngle()
{
    bounded_vector<double, 2> inclination = points_[points_.size() - 1]->getCoordinates() - points_[points_.size() - 2]->getCoordinates();
    double angle = atan2(inclination(1), inclination(0));
    return angle;
}

std::pair<std::vector<Point *>, std::vector<Point *>> Crack::getPointsOfLocalGeometry()
{
    return std::make_pair(geometryBoundary1_, geometryBoundary2_);
}

void Crack::addLastPointsOfLocalGeometry(Point *boundary1, Point *boundary2)
{
    geometryBoundary1_.push_back(boundary1);
    geometryBoundary2_.push_back(boundary2);
}

void Crack::setSurface(PlaneSurface *newSurface)
{
    surface_ = newSurface;
}

std::string Crack::getOpenBoundary()
{
    return openBoundary_;
}

void Crack::setAuxPoint(Point *auxPoint)
{
    auxPoint_ = auxPoint;
}

Point *Crack::getAuxPoint()
{
    return auxPoint_;
}
