#include "Geometry.h"

Geometry::Geometry() {}

Geometry::Geometry(const std::string &domain)
{
    std::string aux = domain;
    while (aux != "GLOBAL" and aux != "LOCAL")
    {
        std::cout << "Is the created geometry GLOBAL or LOCAL?" << std::endl;
        std::cin >> aux;
    }
    domain_ = aux;
    remeshNumber_ = 0;
}

Geometry::~Geometry() {}

Point *Geometry::getPoint(const std::string &name)
{
    return points_[name];
}

std::unordered_map<std::string, Point *> Geometry::getPoints()
{
    return points_;
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

Crack *Geometry::getCrack(const std::string &name)
{
    return crackes_[name];
}

double Geometry::getRadiusJintegral()
{
    return Jradius_;
}

std::string Geometry::getTypeOfDomain()
{
    return domain_;
}

Point *Geometry::addPoint(std::vector<double> coordinates, const double &lcar, const bool &discretization)
{
    int index = points_.size();
    std::stringstream name;
    name << "p" << index;
    bounded_vector<double, 2> coord;
    coord(0) = coordinates[0];
    coord(1) = coordinates[1];
    Point *p = new Point(name.str(), coord, lcar, discretization);
    points_[name.str()] = p;
    return p;
}

void Geometry::addCrackPoint(const std::string &crackName, const bounded_vector<double, 2> &coordinates)
{
    if (domain_ == "GLOBAL")
    {
        int index = points_.size();
        std::stringstream name;
        name << "p" << index;
        Point *p = new Point(name.str(), coordinates, crackes_[crackName]->getLastPoint()->getlcar(), true);
        points_[name.str()] = p;

        crackes_[crackName]->addLastPoint(p);
    }
    else //LOCAL
    {
        std::stringstream name;
        name << "p" << points_.size();
        Point *p = new Point(name.str(), coordinates, crackes_[crackName]->getLastPoint()->getlcar(), true);
        points_[name.str()] = p;

        crackes_[crackName]->addLastPoint(p);
    }
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

void Geometry::addCrack(std::vector<Point *> points, PlaneSurface *surface, const std::string &openBoundary, const double &jradius, const double &lcarJ)
{
    int index = crackes_.size();
    std::stringstream name;
    name << "c" << index;

    Crack *c = new Crack(name.str(), points, surface, openBoundary);
    crackes_[name.str()] = c;

    Jradius_ = jradius;
    lcarJ_ = lcarJ;
}

void Geometry::addCrackOnGlobal(std::vector<Point *> points, const std::string &openBoundary, const double &jradius, const double &lcarJ, const double &lengthOffset, const double &lcarOfLocalBoundary)
{
    int index = crackes_.size();
    std::stringstream name;
    name << "c" << index;

    PlaneSurface *surface;

    Crack *c = new Crack(name.str(), points, surface, openBoundary);
    crackes_[name.str()] = c;

    Jradius_ = jradius;
    lcarJ_ = lcarJ;
    lengthOffset_ = lengthOffset;
    lcarOfLocalBoundary_ = lcarOfLocalBoundary;
}

// void Geometry::transfiniteLine(std::vector<Line *> lines, const int &divisions, const double &progression)
// {
//     std::stringstream text;
//     text << "Transfinite Line {";
//     for (size_t i = 0; i < lines.size(); i++)
//     {
//         text << lines[i]->getName();
//         if (i != (lines.size() - 1))
//             text << ", ";
//     }
//     text << "} = " << divisions << " Using Progression " << progression << ";\n//\n";
//     gmshCode_ += text.str();
// }

// void Geometry::transfiniteSurface(std::vector<PlaneSurface *> planeSurfaces, std::string orientation, std::vector<Point *> points)
// {
//     std::stringstream text;
//     text << "Transfinite Surface {";
//     for (size_t i = 0; i < planeSurfaces.size(); i++)
//     {
//         text << planeSurfaces[i]->getName();
//         if (i != (planeSurfaces.size() - 1))
//             text << ", ";
//     }
//     text << "} ";
//     if (points.size() != 0)
//     {
//         text << "= {";
//         for (size_t i = 0; i < points.size(); i++)
//         {
//             text << points[i]->getName();
//             if (i != (points.size() - 1))
//                 text << ", ";
//         }
//         text << "} " << orientation << ";\n//\n";
//     }
//     else
//     {
//         text << orientation << ";\n//\n";
//     }
//     gmshCode_ += text.str();
// }

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
        Point *p1aux = new Point("p" + std::to_string(npoints++), coordaux, lcarJ_, true);
        gmshCode += p1aux->getGmshCode();

        coordaux = coordLastPoint - versorNormal * Jradius_;
        Point *p2aux = new Point("p" + std::to_string(npoints++), coordaux, lcarJ_, true);
        gmshCode += p2aux->getGmshCode();

        coordaux = coordLastPoint + versor * Jradius_;
        Point *p3aux = new Point("p" + std::to_string(npoints++), coordaux, lcarJ_, true);
        gmshCode += p3aux->getGmshCode();

        coordaux = coordLastPoint + versorNormal * Jradius_;
        Point *p4aux = new Point("p" + std::to_string(npoints++), coordaux, lcarJ_, true);
        gmshCode += p4aux->getGmshCode();

        gmshCode += nameOfCrack + " = newl;\n //\n";

        if (pointsCrack.size() > 2)
        {
            for (int i = 0; i < pointsCrack.size() - 2; i++)
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
    remeshNumber_ += 1;
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

bounded_matrix<double, 2, 2> Geometry::inverseMatrix(const bounded_matrix<double, 2, 2> &matrix)
{
    bounded_matrix<double, 2, 2> inverse;
    double detinv = 1.0 / (matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0));

    inverse(0, 0) = detinv * matrix(1, 1);
    inverse(1, 0) = -detinv * matrix(1, 0);
    inverse(0, 1) = -detinv * matrix(0, 1);
    inverse(1, 1) = detinv * matrix(0, 0);

    return inverse;
}

void Geometry::createGeometryFromCrack()
{
    if (remeshNumber_ == 0) //Primeira malha criada
    {
        for (std::unordered_map<std::string, Crack *>::const_iterator crack = crackes_.begin(); crack != crackes_.end(); crack++)
        {
            std::vector<Point *> crackPoints = crack->second->getPoints();
            std::string open = crack->second->getOpenBoundary();

            if (crackPoints.size() > 2)
            {
                std::cout << "NÃO FOI IMPLEMENTADO FISSURAS INICIAIS COM MAIS DE DOIS PONTOS." << std::endl;
            }
            double initialAngle = crack->second->getLastAngle();
            const double pi = 3.14159265359;

            bounded_vector<double, 2> coord1, coord2, coordInitial, coordLast;

            coordInitial = crackPoints[0]->getCoordinates();
            coordLast = crackPoints[1]->getCoordinates();

            if (open == "first")
            {
                coord1 = coordInitial;
                coord2 = coordInitial;

                coord1(0) += cos(initialAngle - 0.5 * pi) * lengthOffset_;
                coord1(1) += sin(initialAngle - 0.5 * pi) * lengthOffset_;

                coord2(0) += cos(initialAngle + 0.5 * pi) * lengthOffset_;
                coord2(1) += sin(initialAngle + 0.5 * pi) * lengthOffset_;
            }
            else
            {
                coord1 = coordInitial;
                coord2 = coordInitial;

                double length = lengthOffset_ / sin(0.25 * pi);

                coord1(0) += cos(initialAngle - 0.75 * pi) * length;
                coord1(1) += sin(initialAngle - 0.75 * pi) * length;

                coord2(0) += cos(initialAngle + 0.75 * pi) * length;
                coord2(1) += sin(initialAngle + 0.75 * pi) * length;
            }

            Point *paux1_1 = addPoint({coord1(0), coord1(1)}, lcarOfLocalBoundary_, false);
            Point *paux1_2 = addPoint({coord2(0), coord2(1)}, lcarOfLocalBoundary_, false);

            crack->second->addLastPointsOfLocalGeometry(paux1_1, paux1_2);

            if (open == "second")
            {
                coord1 = coordLast;
                coord2 = coordLast;

                coord1(0) += cos(initialAngle - 0.5 * pi) * lengthOffset_;
                coord1(1) += sin(initialAngle - 0.5 * pi) * lengthOffset_;

                coord2(0) += cos(initialAngle + 0.5 * pi) * lengthOffset_;
                coord2(1) += sin(initialAngle + 0.5 * pi) * lengthOffset_;
            }
            else
            {
                coord1 = coordLast;
                coord2 = coordLast;

                double length = lengthOffset_ / sin(0.25 * pi);

                coord1(0) += cos(initialAngle - 0.25 * pi) * length;
                coord1(1) += sin(initialAngle - 0.25 * pi) * length;

                coord2(0) += cos(initialAngle + 0.25 * pi) * length;
                coord2(1) += sin(initialAngle + 0.25 * pi) * length;
            }
            Point *paux1 = addPoint({coord1(0), coord1(1)}, lcarOfLocalBoundary_, false);
            Point *paux2 = addPoint({coord2(0), coord2(1)}, lcarOfLocalBoundary_, false);

            crack->second->addLastPointsOfLocalGeometry(paux1, paux2);

            std::pair<std::vector<Point *>, std::vector<Point *>> boundaryPoints = crack->second->getPointsOfLocalGeometry();

            std::vector<Point *> conec;

            int aux = boundaryPoints.first.size();
            if (open == "first")
            {
                bounded_vector<double, 2> auxCoord = 0.5 * (boundaryPoints.first[aux - 1]->getCoordinates() + boundaryPoints.second[aux - 1]->getCoordinates());
                Point *paux = addPoint({auxCoord(0), auxCoord(1)}, lcarOfLocalBoundary_, false);
                crack->second->setAuxPoint(paux);
                for (int i = 0; i < aux; i++)
                {
                    conec.push_back(boundaryPoints.first[i]);
                }
                conec.push_back(paux);
                for (int i = 0; i < aux; i++)
                {
                    conec.push_back(boundaryPoints.second[aux - 1 - i]);
                }
            }
            else if (open == "second")
            {
                bounded_vector<double, 2> auxCoord = 0.5 * (boundaryPoints.first[0]->getCoordinates() + boundaryPoints.second[0]->getCoordinates());
                Point *paux = addPoint({auxCoord(0), auxCoord(1)}, lcarOfLocalBoundary_, false);
                crack->second->setAuxPoint(paux);

                conec.push_back(crackPoints[1]);
                for (int i = 0; i < aux; i++)
                {
                    conec.push_back(boundaryPoints.second[aux - 1 - i]);
                }
                conec.push_back(paux);
                for (int i = 0; i < aux; i++)
                {
                    conec.push_back(boundaryPoints.first[i]);
                }
                conec.push_back(crackPoints[1]);
            }
            else
            {
                for (int i = 0; i < aux; i++)
                {
                    conec.push_back(boundaryPoints.first[i]);
                }
                for (int i = 0; i < aux; i++)
                {
                    conec.push_back(boundaryPoints.second[aux - 1 - i]);
                }
                conec.push_back(boundaryPoints.first[0]);
            }

            Line *l0 = addBSpline(conec, true);
            LineLoop *ll;
            if (open == "first")
            {
                Line *l1 = addLine({boundaryPoints.second[0], crackPoints[0]});
                Line *l2 = addLine({crackPoints[0], boundaryPoints.first[0]});
                ll = addLineLoop({l0, l1, l2});
            }
            else if (open == "second")
            {
                Line *l1 = addLine({boundaryPoints.first[aux - 1], crackPoints[aux - 1]});
                Line *l2 = addLine({crackPoints[aux - 1], boundaryPoints.second[aux - 1]});
                ll = addLineLoop({l0, l1, l2});
            }
            else
            {
                ll = addLineLoop({l0});
            }

            PlaneSurface *s = addPlaneSurface({ll}, indexMaterial_, thickness_);
            crack->second->setSurface(s);
        }
        remeshNumber_ += 1;
    }
    else //AQUI PRECISAMOS ATUALIZAR A POSIÇÃO DO ÚLTIMO BOUNDARYPOINTS E ADICIONAR MAIS UM
    {
        for (std::unordered_map<std::string, Crack *>::const_iterator crack = crackes_.begin(); crack != crackes_.end(); crack++)
        {
            std::vector<Point *> crackPoints = crack->second->getPoints();
            std::string open = crack->second->getOpenBoundary();

            std::pair<std::vector<Point *>, std::vector<Point *>> boundaryPoints = crack->second->getPointsOfLocalGeometry();

            if (crackPoints.size() > boundaryPoints.first.size() and crackPoints.size() > boundaryPoints.second.size()) //FOI ADICIONADO OUTRO PONTO NA FISSURA, LOGO PRECISAMOS RECONSTRUIR O CONTORNO
            {
                const double pi = 3.14159265359;
                bounded_vector<double, 2> coord0, coord1, coord2;
                int aux = crackPoints.size();
                coord0 = crackPoints[aux - 3]->getCoordinates();
                coord1 = crackPoints[aux - 2]->getCoordinates();
                coord2 = crackPoints[aux - 1]->getCoordinates();

                double teta1, teta0;

                teta0 = atan2(coord1(1) - coord0(1), coord1(0) - coord0(0));
                teta1 = atan2(coord2(1) - coord1(1), coord2(0) - coord1(0));

                bounded_matrix<double, 2, 2> mat;
                mat(0, 0) = cos(teta0);
                mat(1, 0) = sin(teta0);
                mat(0, 1) = -cos(teta1);
                mat(1, 1) = -sin(teta1);

                //MODIFICANDO ÚLTIMO DO CONTORNO 1
                bounded_vector<double, 2> boundCoord0, auxCoord1;
                boundCoord0 = boundaryPoints.first[aux - 3]->getCoordinates();
                auxCoord1(0) = coord1(0) + cos(teta1 - 0.5 * pi) * lengthOffset_;
                auxCoord1(1) = coord1(1) + sin(teta1 - 0.5 * pi) * lengthOffset_;

                bounded_vector<double, 2> auxDelta = auxCoord1 - boundCoord0;
                bounded_vector<double, 2> values = prod(inverseMatrix(mat), auxDelta);

                bounded_vector<double, 2> newCoord;
                newCoord(0) = boundCoord0(0) + cos(teta0) * values(0);
                newCoord(1) = boundCoord0(1) + sin(teta0) * values(0);
                boundaryPoints.first[aux - 2]->setCoordinates(newCoord);

                //MODIFICANDO ÚLTIMO DO CONTORNO 2
                boundCoord0 = boundaryPoints.second[aux - 3]->getCoordinates();
                auxCoord1(0) = coord1(0) + cos(teta1 + 0.5 * pi) * lengthOffset_;
                auxCoord1(1) = coord1(1) + sin(teta1 + 0.5 * pi) * lengthOffset_;

                auxDelta = auxCoord1 - boundCoord0;
                values = prod(inverseMatrix(mat), auxDelta);

                newCoord(0) = boundCoord0(0) + cos(teta0) * values(0);
                newCoord(1) = boundCoord0(1) + sin(teta0) * values(0);
                boundaryPoints.second[aux - 2]->setCoordinates(newCoord);

                //CRIANDO DOIS ÚLTIMOS PONTOS DO CONTORNO
                bounded_vector<double, 2> coord2_1, coord2_2;
                coord2_1 = coord2;
                coord2_2 = coord2;

                double length = lengthOffset_ / sin(0.25 * pi);

                coord2_1(0) += cos(teta1 - 0.25 * pi) * length;
                coord2_1(1) += sin(teta1 - 0.25 * pi) * length;

                coord2_2(0) += cos(teta1 + 0.25 * pi) * length;
                coord2_2(1) += sin(teta1 + 0.25 * pi) * length;

                Point *paux1 = addPoint({coord2_1(0), coord2_1(1)}, lcarOfLocalBoundary_, false);
                Point *paux2 = addPoint({coord2_2(0), coord2_2(1)}, lcarOfLocalBoundary_, false);
                crack->second->addLastPointsOfLocalGeometry(paux1, paux2);

                //EDITANDO BSPLINE DO CONTORNO LOCAL

                std::vector<Point *> conec;

                boundaryPoints = crack->second->getPointsOfLocalGeometry();

                if (open == "first") //IMPLEMENTADO SOMENTE PARA ESSE TIPO DE PROPAGAÇÃO!
                {
                    bounded_vector<double, 2> auxCoord = 0.5 * (boundaryPoints.first[aux - 1]->getCoordinates() + boundaryPoints.second[aux - 1]->getCoordinates());
                    crack->second->getAuxPoint()->setCoordinates(auxCoord);
                    for (int i = 0; i < aux; i++)
                    {
                        conec.push_back(boundaryPoints.first[i]);
                    }
                    conec.push_back(crack->second->getAuxPoint());
                    for (int i = 0; i < aux; i++)
                    {
                        conec.push_back(boundaryPoints.second[aux - 1 - i]);
                    }
                }
                else if (open == "second")
                {
                    conec.push_back(crackPoints[1]);
                    for (int i = 0; i < aux; i++)
                    {
                        conec.push_back(boundaryPoints.second[aux - 1 - i]);
                    }
                    for (int i = 0; i < aux; i++)
                    {
                        conec.push_back(boundaryPoints.first[i]);
                    }
                    conec.push_back(crackPoints[1]);
                }
                else
                {
                    for (int i = 0; i < aux; i++)
                    {
                        conec.push_back(boundaryPoints.first[i]);
                    }
                    for (int i = 0; i < aux; i++)
                    {
                        conec.push_back(boundaryPoints.second[aux - 1 - i]);
                    }
                    conec.push_back(boundaryPoints.first[0]);
                }

                crack->second->getPlaneSurface()->getLineLoop(0)->getLine(0)->setPoints(conec);
            }
            else
            {
                break;
            }
        }
        remeshNumber_ += 1;
    }
}