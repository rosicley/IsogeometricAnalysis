#include "Geometry.h"

Geometry::Geometry() {}

Geometry::~Geometry() {}

// int Geometry::getNumberOfBoundaryConditions(const std::string &type)
// {
//     return boundaryConditions_[type].size();
// }

Point *Geometry::getPoint(const std::string &name)
{
    return points_[name];
}

Line *Geometry::getLine(const std::string &name)
{
    return lines_[name];
}

LineLoop *Geometry::getLineLoop(const std::string &name)
{
    return lineLoops_[name];
}

PlaneSurface *Geometry::getPlaneSurface(const std::string &name)
{
    return planeSurfaces_[name];
}

std::unordered_map<std::string, PlaneSurface *> Geometry::getPlaneSurfaces()
{
    return planeSurfaces_;
}

std::unordered_map<std::string, Crack *> Geometry::getCrackes()
{
    return crackes_;
}

std::unordered_map<std::string, bounded_vector<double, 2>> Geometry::getDirichletCondition()
{
    return diricheletConditions_;
}

std::unordered_map<std::string, bounded_vector<double, 2>> Geometry::getNeumannCondition()
{
    return neumannConditions_;
}

// std::vector<GeometricBoundaryCondition *> Geometry::getBoundaryConditions(const std::string &type)
// {
//     return boundaryConditions_[type];
// }

std::string Geometry::getGmshCode()
{
    return gmshCode_;
}

Point *Geometry::addPoint(std::vector<double> coordinates, const double &lcar, const bool &discretization)
{
    int index = points_.size();
    std::stringstream name;
    name << "p" << index;
    Point *p = new Point(name.str(), coordinates, lcar, discretization);
    points_[name.str()] = p;
    gmshCode_ += p->getGmshCode();
    return p;
}

Line *Geometry::addLine(std::vector<Point *> points)
{
    int index = lines_.size();
    std::stringstream name;
    name << "l" << index;
    Line *l = new Line(name.str(), points);
    lines_[name.str()] = l;
    gmshCode_ += l->getGmshCodeLine();
    return l;
}

Line *Geometry::addCircle(std::vector<Point *> points)
{
    int index = lines_.size();
    std::stringstream name;
    name << "l" << index;
    Line *l = new Line(name.str(), points);
    lines_[name.str()] = l;
    gmshCode_ += l->getGmshCodeCircle();
    return l;
}

Line *Geometry::addCrackLine(std::vector<Point *> points)
{
    std::stringstream name;
    name << "c" << crack_++;
    Line *l = new Line(name.str(), points);
    lines_[name.str()] = l;
    gmshCode_ += l->getGmshCodeLine();
    return l;
}

LineLoop *Geometry::addLineLoop(std::vector<Line *> lines)
{
    int index = lineLoops_.size();
    std::stringstream name;
    name << "ll" << index;
    LineLoop *ll = new LineLoop(name.str(), lines);
    ll->verification();
    lineLoops_[name.str()] = ll;
    gmshCode_ += ll->getGmshCode();
    return ll;
}

PlaneSurface *Geometry::addPlaneSurface(std::vector<LineLoop *> lineLoop, const int &indexMaterial, const double &thickness)
{
    int index = planeSurfaces_.size();
    std::stringstream name;
    name << "s" << index;
    PlaneSurface *s = new PlaneSurface(name.str(), lineLoop, indexMaterial, thickness);
    planeSurfaces_[name.str()] = s;
    gmshCode_ += s->getGmshCode();
    return s;
}

Crack *Geometry::addCrack(Line *line, PlaneSurface *surface, const std::string &openBoundary)
{
    int index = crackes_.size();
    std::stringstream name;
    name << "c" << index;

    Crack *c = new Crack(name.str(), line, surface, openBoundary);
    crackes_[name.str()] = c;

    std::stringstream text;
    text << "Line{" << name.str() << "} In Surface{" << surface->getName() << "};\n//\n";
    gmshCode_ += text.str();
}

// PlaneSurface *Geometry::addPlaneSurface(std::vector<Line *> lines, double thickness)
// {
//     int index = getNumberOfPlaneSurfaces();
//     std::stringstream name;
//     name << "s" << index;
//     LineLoop *ll = addLineLoop(lines);
//     std::vector<LineLoop *> lineloops;
//     lineloops.push_back(ll);
//     PlaneSurface *s = new PlaneSurface(index, name.str(), lineloops, thickness);
//     planeSurfaces_[s->getName()] = s;
//     gmshCode_ += s->getGmshCode();
//     return s;
// }

void Geometry::appendGmshCode(std::string text)
{
    gmshCode_ += text;
}

void Geometry::transfiniteLine(std::vector<Line *> lines, const int &divisions, const double &progression)
{
    std::stringstream text;
    text << "Transfinite Line {";
    for (size_t i = 0; i < lines.size(); i++)
    {
        text << lines[i]->getName();
        if (i != (lines.size() - 1))
            text << ", ";
    }
    text << "} = " << divisions << " Using Progression " << progression << ";\n//\n";
    gmshCode_ += text.str();
}

void Geometry::transfiniteSurface(std::vector<PlaneSurface *> planeSurfaces, std::string orientation, std::vector<Point *> points)
{
    std::stringstream text;
    text << "Transfinite Surface {";
    for (size_t i = 0; i < planeSurfaces.size(); i++)
    {
        text << planeSurfaces[i]->getName();
        if (i != (planeSurfaces.size() - 1))
            text << ", ";
    }
    text << "} ";
    if (points.size() != 0)
    {
        text << "= {";
        for (size_t i = 0; i < points.size(); i++)
        {
            text << points[i]->getName();
            if (i != (points.size() - 1))
                text << ", ";
        }
        text << "} " << orientation << ";\n//\n";
    }
    else
    {
        text << orientation << ";\n//\n";
    }
    gmshCode_ += text.str();
}

void Geometry::addNeumannCondition(Line *line, const std::vector<double> &directionX1, const std::vector<double> &directionX2)
{
    bounded_vector<double, 2> force;
    if (directionX1.size() > 0)
    {
        force(0) = directionX1[0];
    }
    else
    {
        force(0) = 0.0;
    }

    if (directionX2.size() > 0)
    {
        force(1) = directionX2[0];
    }
    else
    {
        force(1) = 0.0;
    }

    neumannConditions_[line->getName()] = force;
}

void Geometry::addNeumannCondition(Point *point, const std::vector<double> &directionX1, const std::vector<double> &directionX2)
{
    bounded_vector<double, 2> force;
    if (directionX1.size() > 0)
    {
        force(0) = directionX1[0];
    }
    else
    {
        force(0) = 0.0;
    }

    if (directionX2.size() > 0)
    {
        force(1) = directionX2[0];
    }
    else
    {
        force(1) = 0.0;
    }

    neumannConditions_[point->getName()] = force;
}

void Geometry::addDirichletCondition(Line *line, const std::vector<double> &directionX1, const std::vector<double> &directionX2)
{
    bounded_vector<double, 2> desloc;
    if (directionX1.size() > 0)
    {
        desloc(0) = directionX1[0];
    }
    else
    {
        desloc(0) = 1.0e-240;
    }

    if (directionX2.size() > 0)
    {
        desloc(1) = directionX2[0];
    }
    else
    {
        desloc(1) = 1.0e-240;
    }

    diricheletConditions_[line->getName()] = desloc;
}

void Geometry::addDirichletCondition(Point *point, const std::vector<double> &directionX1, const std::vector<double> &directionX2)
{
    bounded_vector<double, 2> desloc;
    if (directionX1.size() > 0)
    {
        desloc(0) = directionX1[0];
    }
    else
    {
        desloc(0) = 1.0e-240;
    }

    if (directionX2.size() > 0)
    {
        desloc(1) = directionX2[0];
    }
    else
    {
        desloc(1) = 1.0e-240;
    }

    diricheletConditions_[point->getName()] = desloc;
}