#pragma once

#include "Line.h"
#include "PlaneSurface.h"

class Crack
{
public:
    Crack();

    Crack(const std::string& name, Line* line, PlaneSurface *surface, const std::string &openBoundary);

    ~Crack();

    std::string getName();

    Line *getLine();

    PlaneSurface *getPlaneSurface();

    std::string getGmshCode();


private:
    std::string name_;
    std::string openBoundary_;
    Line *line_;
    PlaneSurface *surface_;

};