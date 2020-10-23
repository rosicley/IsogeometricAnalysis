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

  //  std::string getGmshCode();

    // std::string getGmshCodeEmbedded();

    Point *getLastPoint();

    Point *getFirstPoint();

    void addLastPoint(Point *newPoint);

    double getLastAngle();

    std::pair<std::vector<Point *>, std::vector<Point *>> getPointsOfLocalGeometry();

    void addLastPointsOfLocalGeometry(Point *boundary1, Point *boundary2);

    void setSurface(PlaneSurface *newSurface);

    void setAuxPoint(Point *auxPoint);

    Point *getAuxPoint();

    std::string getOpenBoundary();

private:
    std::string name_;
    std::string openBoundary_;
    std::vector<Point *> points_;
    //std::vector<Line *> lines_;
    //std::vector<Line *> linesOfJ_;
    // Line *line_;
    PlaneSurface *surface_;
    std::vector<Point *> geometryBoundary1_; //only local
    std::vector<Point *> geometryBoundary2_; //only local
    Point *auxPoint_; //transition point between geometryBoundary1 and geometryBoundary2

    //std::string embeddedCode_;
};