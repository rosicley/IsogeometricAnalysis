#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include "Isogeometric/Patch.h"
#include "Isogeometric/Cell.h"
#include "BoundaryConditions/DirichletCondition.h"
#include "BoundaryConditions/NeumannCondition.h"
#include "Material.h"
#include <boost/timer.hpp>
#include <boost/thread.hpp>
#include <metis.h>
#include <petscksp.h>

using namespace boost::numeric::ublas;

class GlobalSolid
{
public:
    GlobalSolid();

    GlobalSolid(const std::string &planeState);

    ~GlobalSolid();

    void addMaterial(const int &index,
                     const double &young,
                     const double &poisson,
                     const double &density);

    void addDirichletCondition(const int &index,
                               const int &patchIndex,
                               const bounded_vector<int, 2> &free,
                               const bounded_vector<double, 2> &values);

    void addNeumannCondition(const int &index,
                             const int &patchIndex,
                             const bounded_vector<double, 2> &values);

    void addPatch(const int &index, 
                  const int &npc, 
                  const int &indexMaterial, 
                  const double &thickness);

    int solveStaticProblem();

    int solveDynamicProblem();

    int firstAccelerationCalculation();

    void exportToParaview(const int &loadstep);

    void dataReading(const std::string &inputParameters, 
                     const std::string &inputProperties, 
                     const std::string &inputMeshIso,
                     const bool &rhino);

    Material *getMaterial(const int &index);

    Patch *getPatch(const int &index);

    std::string getAnalysisType();

    std::string getPlaneState();

    matrix<double> coordinatesForInterpolation(const int &orderElemement);

    void ISOdomainDecompositionMETIS();


private:
    std::vector<Patch *> patches_;

    std::vector<DirichletCondition *> dirichletConditions_;

    std::vector<NeumannCondition *> neumannConditions_;

    std::vector<Material *> materials_;

    //std::string elementType_;

    std::string planeState_;

    std::string problemType_;

    double deltat_;

    double gamma_;

    double beta_;

    bounded_vector<double, 2> shapeForces_;

    //int numberOfHammer_;

    //int order_;

    idx_t *cellPartition_;

    idx_t *pointsPartition_;

    int numberOfSteps_;

    int maximumOfIteration_;

    double tolerance_;

    std::vector<Cell *> cells_;

    std::vector<Cell *> cells_part;

    std::vector<ControlPoint *> controlPoints_; 

    int orderParaview_;

    int quadrature_;

    int cpnumber_;
};