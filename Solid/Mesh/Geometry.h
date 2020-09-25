#pragma once

#include "PlaneSurface.h"
#include "Crack.h"
#include <math.h>
//#include "GeometricBoundaryCondition.h"
#include <unordered_map>

class Geometry
{
public:
    Geometry();

    Geometry(const std::string &domain);

    ~Geometry();

    int getNumberOfPoints();

    int getNumberOfLines();

    int getNumberOfLineLoops();

    int getNumberOfPlaneSurfaces();

    int getNumberOfCrackes();

    int getNumberOfBoundaryConditions(const std::string &type);

    Point *getPoint(const std::string &name);

    Line *getLine(const std::string &name);

    LineLoop *getLineLoop(const std::string &name);

    PlaneSurface *getPlaneSurface(const std::string &name);

    std::unordered_map<std::string, PlaneSurface *> getPlaneSurfaces();

    std::unordered_map<std::string, Crack *> getCrackes();

    std::unordered_map<std::string, bounded_vector<double, 2>> getDirichletCondition();

    std::unordered_map<std::string, bounded_vector<double, 2>> getNeumannCondition();

    //std::vector<GeometricBoundaryCondition*> getBoundaryConditions(const std::string& type);

    std::string getGmshCode();

    std::string createGmshCode();

    Point *addPoint(std::vector<double> coordinates, const double &lcar = 1.0, const bool &discretization = true);

    Line *addLine(std::vector<Point *> points, const bool& discretization = true);

    Line *addCircle(std::vector<Point *> points, const bool& discretization = true); //{initial point, center point, end point}

    Line *addEllipse(std::vector<Point *> points, const bool& discretization = true); //{initial point, center point, a point in major axis, end point}

    // Line *addCrackLine(std::vector<Point *> points);

    Line *addSpline(std::vector<Point *> points, const bool& discretization = true);

    Line *addBSpline(std::vector<Point *> points, const bool& discretization = true);

    LineLoop *addLineLoop(std::vector<Line *> lines, const bool &verify = true);

    PlaneSurface *addPlaneSurface(std::vector<LineLoop *> lineLoop, const int &indexMaterial = 0, const double &thickness = 1.0);

    Crack *addCrack(std::vector<Point *> points, PlaneSurface *surface, const std::string &openBoundary, const double &jradius, const double &lcarOfJintegral, const double &lengthMesh = 0.0);

    // PlaneSurface *addPlaneSurface(std::vector<Line *> lines, double thickness = 1.0);

    void appendGmshCode(std::string text);

    void transfiniteLine(std::vector<Line *> lines, const int &divisions, const double &progression = 1);

    void transfiniteSurface(std::vector<PlaneSurface *> surfaces, std::string oientation = "Left", std::vector<Point *> points = std::vector<Point *>());

    void addNeumannCondition(Line *line, const std::vector<double> &directionX1, const std::vector<double> &directionX2);

    void addNeumannCondition(Point *point, const std::vector<double> &directionX1, const std::vector<double> &directionX2);

    void addDirichletCondition(Line *line, const std::vector<double> &directionX1, const std::vector<double> &directionX2);

    void addDirichletCondition(Point *point, const std::vector<double> &directionX1, const std::vector<double> &directionX2);

    Crack *getCrack(const std::string &name);

    double getRadiusJintegral();


    // void addBoundaryCondition(const std::string &type, Point *point, const std::vector<double> &componentX = std::vector<double>(), const std::vector<double> &componentY = std::vector<double>(), const std::string &referenceSystem = "GLOBAL", const std::string &method = "STRONG", const double &penaltyParameter = 1.0e6);

    // void addBoundaryCondition(const std::string &type, Line *line, const std::vector<double> &componentX = std::vector<double>(), const std::vector<double> &componentY = std::vector<double>(), const std::string &referenceSystem = "GLOBAL", const std::string &method = "STRONG", const double &penaltyParameter = 1.0e6);

private:
    std::unordered_map<std::string, Point *> points_;
    std::unordered_map<std::string, Line *> lines_;
    std::unordered_map<std::string, LineLoop *> lineLoops_;
    std::unordered_map<std::string, PlaneSurface *> planeSurfaces_;
    std::unordered_map<std::string, Crack *> crackes_;

    std::unordered_map<std::string, bounded_vector<double, 2>> neumannConditions_;
    std::unordered_map<std::string, bounded_vector<double, 2>> diricheletConditions_;

    //std::unordered_map<std::string, std::vector<GeometricBoundaryCondition*>> boundaryConditions_;
    std::string gmshCode_;
    std::string domain_; //local ou global
    double Jradius_;
    double lengthMesh_;
    double lcarJ_;
};
