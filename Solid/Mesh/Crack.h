#pragma once

#include "Line.h"
#include "PlaneSurface.h"

class Crack
{
public:
    Crack();

    Crack(const std::string &name, std::vector<Point *> points, PlaneSurface *surface, const std::string &openBoundary);

    ~Crack();

    std::string getName();

    std::vector<Point *> getPoints();

    PlaneSurface *getPlaneSurface();

    std::string getGmshCodeCrackPlugin();

    std::string getGmshCode();

    std::string getGmshCodeEmbedded();

    Point *getLastPoint();

    Point *getFirstPoint();



private:
    std::string name_;
    std::string openBoundary_;
    std::vector<Point *> points_;
    //std::vector<Line *> lines_;
    //std::vector<Line *> linesOfJ_; 
    // Line *line_;
    PlaneSurface *surface_;
  
    std::string embeddedCode_;
};