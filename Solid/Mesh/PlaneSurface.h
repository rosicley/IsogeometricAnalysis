#pragma once

#include "LineLoop.h"
#include "../Material.h"

class PlaneSurface
{
public:
    PlaneSurface();

    PlaneSurface(const std::string &name, std::vector<LineLoop *>lineLoop, const int &indexMaterial, const double &thickness = 1.0);

    ~PlaneSurface();

    std::string getName();

    double getThickness();

    int getIndexMaterial();

    std::string getGmshCode();

    LineLoop *getLineLoop(const int &index);

private:
    std::string name_;
    double thickness_;
    std::vector<LineLoop *> lineLoop_;
    int indexMaterial_;
    //std::vector<Element *> elements_;
};
