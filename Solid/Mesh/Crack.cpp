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

std::string Crack::getGmshCode()
{
}

std::string Crack::getGmshCodeEmbedded()
{
    return embeddedCode_;
}