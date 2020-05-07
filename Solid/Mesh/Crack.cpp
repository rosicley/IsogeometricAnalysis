#include "Crack.h"

Crack::Crack() {}

Crack::Crack(const std::string &name, Line *line, PlaneSurface *surface, const std::string &openBoundary)
{
    name_ = name;
    line_ = line;
    surface_ = surface;
    openBoundary_ = openBoundary;
}

Crack::~Crack() {}

std::string Crack::getName()
{
    return name_;
}

Line *Crack::getLine()
{
    return line_;
}

PlaneSurface *Crack::getPlaneSurface()
{
    return surface_;
}

std::string Crack::getGmshCode()
{
    std::stringstream text;
    text << "Plugin(Crack).Dimension = 1;\n";
    text << "Plugin(Crack).PhysicalGroup = " << line_->getName() << ";\n";
    if (openBoundary_ == "first")
    {
        text << "Plugin(Crack).OpenBoundaryPhysicalGroup = " << line_->getInitialPoint()->getName() << ";\n";
    }
    else if (openBoundary_ == "second")
    {
        text << "Plugin(Crack).OpenBoundaryPhysicalGroup = " << line_->getEndPoint()->getName() << ";\n";
    }
    text << "Plugin(Crack).Run; \n//\n";

    return text.str();
}