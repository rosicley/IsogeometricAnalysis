#include "PlaneSurface.h"

PlaneSurface::PlaneSurface() {}

PlaneSurface::PlaneSurface(const std::string &name, std::vector<LineLoop *> lineLoop, const int &indexMaterial, const double &thickness)
{
    name_ = name;
    lineLoop_.reserve(lineLoop.size());
    lineLoop_ = lineLoop;
    indexMaterial_ = indexMaterial;
    thickness_ = thickness;
}

PlaneSurface::~PlaneSurface() {}

std::string PlaneSurface::getName()
{
    return name_;
}

double PlaneSurface::getThickness()
{
    return thickness_;
}

int PlaneSurface::getIndexMaterial()
{
    return indexMaterial_;
}

std::string PlaneSurface::getGmshCode()
{
    std::stringstream text;
    text << name_ << " = news; Plane Surface(" << name_ << ") = {";
    for (int i = 0; i < lineLoop_.size(); i++)
    {
        text << lineLoop_[i]->getName();
        if (lineLoop_.size() - i != 1)
        {
            text << ", ";
        }
    }

    text << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";
    return text.str();
}

LineLoop *PlaneSurface::getLineLoop(const int &index)
{
    return lineLoop_[index];
}