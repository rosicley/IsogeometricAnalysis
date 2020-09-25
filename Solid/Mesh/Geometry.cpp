#include "Geometry.h"

Geometry::Geometry() {}

Geometry::Geometry(const std::string &domain)
{
    domain_ = domain;
}

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

std::string Geometry::getGmshCode()
{
    return gmshCode_;
}

Crack *Geometry::getCrack(const std::string &name)
{
    return crackes_[name];
}

double Geometry::getRadiusJintegral()
{
    return Jradius_;
}

Point *Geometry::addPoint(std::vector<double> coordinates, const double &lcar, const bool &discretization)
{
    int index = points_.size();
    std::stringstream name;
    name << "p" << index;
    Point *p = new Point(name.str(), coordinates, lcar, discretization);
    points_[name.str()] = p;
    return p;
}

Line *Geometry::addLine(std::vector<Point *> points, const bool &discretization)
{
    int index = lines_.size();
    std::stringstream name;
    name << "l" << index;
    Line *l = new Line(name.str(), points, "line", discretization);
    lines_[name.str()] = l;
    return l;
}

Line *Geometry::addCircle(std::vector<Point *> points, const bool &discretization)
{
    int index = lines_.size();
    std::stringstream name;
    name << "l" << index;
    Line *l = new Line(name.str(), points, "circle", discretization);
    lines_[name.str()] = l;
    return l;
}

Line *Geometry::addSpline(std::vector<Point *> points, const bool &discretization)
{
    int index = lines_.size();
    std::stringstream name;
    name << "l" << index;
    Line *l = new Line(name.str(), points, "spline", discretization);
    lines_[name.str()] = l;
    return l;
}

Line *Geometry::addBSpline(std::vector<Point *> points, const bool &discretization)
{
    int index = lines_.size();
    std::stringstream name;
    name << "l" << index;
    Line *l = new Line(name.str(), points, "bspline", discretization);
    lines_[name.str()] = l;
    return l;
}

LineLoop *Geometry::addLineLoop(std::vector<Line *> lines, const bool &verify)
{
    int index = lineLoops_.size();
    std::stringstream name;
    name << "ll" << index;
    LineLoop *ll = new LineLoop(name.str(), lines);
    if (verify)
    {
        ll->verification();
    }

    lineLoops_[name.str()] = ll;
    return ll;
}

PlaneSurface *Geometry::addPlaneSurface(std::vector<LineLoop *> lineLoop, const int &indexMaterial, const double &thickness)
{
    int index = planeSurfaces_.size();
    std::stringstream name;
    name << "s" << index;
    PlaneSurface *s = new PlaneSurface(name.str(), lineLoop, indexMaterial, thickness);
    planeSurfaces_[name.str()] = s;
    return s;
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

Crack *Geometry::addCrack(std::vector<Point *> points, PlaneSurface *surface, const std::string &openBoundary, const double &jradius, const double &lcarJ, const double &lengthMesh)
{
    int index = crackes_.size();
    std::stringstream name;
    name << "c" << index;

    Crack *c = new Crack(name.str(), points, surface, openBoundary);
    crackes_[name.str()] = c;

    Jradius_ = jradius;
    lengthMesh_ = lengthMesh;
    lcarJ_ = lcarJ;
    return c;
}

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

std::string Geometry::createGmshCode()
{
    std::string gmshCode;
    for (int in = 0; in < points_.size(); in++)
    {
        std::string name = "p" + std::to_string(in);
        gmshCode += points_[name]->getGmshCode();
    }

    for (int il = 0; il < lines_.size(); il++)
    {
        std::string name = "l" + std::to_string(il);
        gmshCode += lines_[name]->getGmshCode();
    }

    int auxC = 0;
    int npoints = points_.size();

    for (int ic = 0; ic < crackes_.size(); ic++)
    {
        std::string nameOfCrack = "c" + std::to_string(ic);

        std::vector<Point *> pointsCrack = crackes_[nameOfCrack]->getPoints();

        bounded_vector<double, 2> coordLastPoint = pointsCrack[pointsCrack.size() - 1]->getCoordinates();

        bounded_vector<double, 2> auxvec = coordLastPoint - pointsCrack[pointsCrack.size() - 2]->getCoordinates();
        bounded_vector<double, 2> versor = (1.0 / norm_2(auxvec)) * auxvec;
        bounded_vector<double, 2> versorNormal;
        versorNormal(0) = -versor(1);
        versorNormal(1) = versor(0);

        bounded_vector<double, 2> coordaux = coordLastPoint - versor * Jradius_;
        Point *p1aux = new Point("p" + std::to_string(npoints++), {coordaux(0), coordaux(1)}, lcarJ_, true);
        gmshCode += p1aux->getGmshCode();

        coordaux = coordLastPoint - versorNormal * Jradius_;
        Point *p2aux = new Point("p" + std::to_string(npoints++), {coordaux(0), coordaux(1)}, lcarJ_, true);
        gmshCode += p2aux->getGmshCode();

        coordaux = coordLastPoint + versor * Jradius_;
        Point *p3aux = new Point("p" + std::to_string(npoints++), {coordaux(0), coordaux(1)}, lcarJ_, true);
        gmshCode += p3aux->getGmshCode();

        coordaux = coordLastPoint + versorNormal * Jradius_;
        Point *p4aux = new Point("p" + std::to_string(npoints++), {coordaux(0), coordaux(1)}, lcarJ_, true);
        gmshCode += p4aux->getGmshCode();

        gmshCode += nameOfCrack + " = newl;\n //\n";

        if (pointsCrack.size() > 2)
        {
            for (int i = 0; i < pointsCrack.size() - 1; i++)
            {
                Line *l = new Line(nameOfCrack + std::to_string(auxC++), {pointsCrack[i], pointsCrack[i + 1]}, "line", false);
                gmshCode += l->getGmshCode();
                delete l;
            }
        }

        Line *c00 = new Line(nameOfCrack + std::to_string(auxC++), {pointsCrack[pointsCrack.size() - 2], p1aux}, "line", false);
        gmshCode += c00->getGmshCode();
        delete c00;

        Line *c01 = new Line(nameOfCrack + std::to_string(auxC++), {p1aux, pointsCrack[pointsCrack.size() - 1]}, "line", false);
        gmshCode += c00->getGmshCode();
        delete c01;

        gmshCode += "Physical Line('" + nameOfCrack + "') = {";
        for (int i = 0; i < auxC; i++)
        {
            gmshCode += nameOfCrack + std::to_string(i);
            if (i != auxC - 1)
            {
                gmshCode += ", ";
            }
        }
        gmshCode += "}; \n//\n";

        Line *j1 = new Line("j" + std::to_string(ic) + "0", {p1aux, pointsCrack[pointsCrack.size() - 1], p2aux}, "circle", false);
        gmshCode += j1->getGmshCode();

        Line *j2 = new Line("j" + std::to_string(ic) + "1", {p2aux, pointsCrack[pointsCrack.size() - 1], p3aux}, "circle", false);
        gmshCode += j2->getGmshCode();

        Line *j3 = new Line("j" + std::to_string(ic) + "2", {p3aux, pointsCrack[pointsCrack.size() - 1], p4aux}, "circle", false);
        gmshCode += j3->getGmshCode();

        Line *j4 = new Line("j" + std::to_string(ic) + "3", {p4aux, pointsCrack[pointsCrack.size() - 1], p1aux}, "circle", false);
        gmshCode += j4->getGmshCode();

        gmshCode += "Physical Line('j" + std::to_string(ic) + "') = {j" + std::to_string(ic) + "0, j" + std::to_string(ic) + "1, j" + std::to_string(ic) + "2, j" + std::to_string(ic) + "3};\n//\n";
    }
    for (int ill = 0; ill < lineLoops_.size(); ill++)
    {
        std::string name = "ll" + std::to_string(ill);
        gmshCode += lineLoops_[name]->getGmshCode();
    }
    for (int is = 0; is < planeSurfaces_.size(); is++)
    {
        std::string name = "s" + std::to_string(is);
        gmshCode += planeSurfaces_[name]->getGmshCode();
    }
    std::stringstream text;
    for (int ic = 0; ic < crackes_.size(); ic++)
    {
        std::string nameOfCrack = "c" + std::to_string(ic);
        for (int i = 0; i < auxC; i++)
        {
            gmshCode += "Line{" + nameOfCrack + std::to_string(i) + "} In Surface{" + crackes_[nameOfCrack]->getPlaneSurface()->getName() + "};\n//\n";
        }
        for (int i = 0; i < 4; i++)
        {
            gmshCode += "Line{j" + std::to_string(ic) + std::to_string(i) + "} In Surface{" + crackes_[nameOfCrack]->getPlaneSurface()->getName() + "};\n//\n";
        }
    }

    return gmshCode;
}

int Geometry::getNumberOfPoints()
{
    return points_.size();
}

int Geometry::getNumberOfCrackes()
{
    return crackes_.size();
}
