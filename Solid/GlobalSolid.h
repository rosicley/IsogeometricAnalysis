#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <unordered_set>
#include "Isogeometric/Patch.h"
#include "Isogeometric/Cell.h"
#include "BoundaryConditions/DirichletCondition.h"
#include "BoundaryConditions/NeumannCondition.h"
#include "Material.h"
#include "FiniteElement/Element.h"
#include "FiniteElement/Node.h"
#include "FiniteElement/Mesh.h"
#include "FiniteElement/BoundaryElement.h"
#include "FiniteElement/Inactive.h"
#include "Mesh/Geometry.h"
#include <boost/timer.hpp>
#include <boost/thread.hpp>
#include <metis.h>
#include <petscksp.h>
#include <petscvec.h>

using namespace boost::numeric::ublas;
#ifdef _WIN32
#include <direct.h>
#define getCurrentDir _getcwd
#define remove "del "
#else
#include <unistd.h>
#define getCurrentDir getcwd
#define remove "rm "
#endif

class GlobalSolid
{
public:
    GlobalSolid();

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

    void addDirichletConditionFE(const int &index,
                                 const bounded_vector<int, 2> &free,
                                 const bounded_vector<double, 2> &values);

    void addNeumannConditionFE(const int &index,
                               const bounded_vector<double, 2> &values);

    void addPatch(const int &index,
                  const int &npc,
                  const int &indexMaterial,
                  const double &thickness);

    void addMesh(const int &index, const int &indexMaterial, const double &thickness,
                 const std::string &elementType);

    void addNode(const int &index, const int &indexFE,
                 const bounded_vector<double, 2> &initialCoordinate);

    // void addElement(const int &index,
    // const std::vector<int> &nodesIndex,
    // const int &materialIndex,
    // const double &thickness,
    // const std::string &elementType);

    void addBoundaryFiniteElement(const int &lineNumber, const std::vector<int> &nodesIndex);

    int solveStaticProblem();

    int solveDynamicProblem();

    int firstAccelerationCalculation();

    void exportToParaviewISO(const int &loadstep);

    void exportToParaviewFEM(const int &number);

    //void exportToParaviewFEM_Blended(const int &number);

    void exportMirror();

    void dataReading(const std::string &inputParameters,
                     const std::string &inputProperties,
                     const std::string &inputMeshIso,
                     const bool &rhino);

    void dataFromGmsh(const std::string &inputGmesh, const std::string &elementType, Geometry *geometry);

    Material *getMaterial(const int &index);

    Patch *getPatch(const int &index);

    std::string getAnalysisType();

    std::string getPlaneState();

    bounded_vector<double, 2> blendFunction(const double &DA);

    matrix<double> coordinatesForInterpolation(const int &orderElemement);

    matrix<double> coordinatesFEM(const int &orderElement);

    bounded_matrix<double, 2, 2> inverseMatrix(const bounded_matrix<double, 2, 2> &matrix);

    void ISOdomainDecompositionMETIS();

    void domainDecompositionMETIS(const std::string &elementType);

    void incidenceLocalxGlobal();

    void computeDistanceFromFEBoundary();

    void checkInactivesCPandNode();

    void exportToParaviewHammerPoints();

    void teste();

    vector<double> diagonalMassMatrix();

    void shareDataBetweenRanks();

    void stressCalculateFEM();

    void generateMesh(Geometry *geometry, const std::string &elementType = "T3", const std::string &algorithm = "AUTO", std::string geofile = std::string(),
                      const bool &plotMesh = true, const bool &showInfo = false);

    void applyBoundaryConditions(Geometry *geometry);

    void addBlending(std::vector<Line *> lines, const double &thickness);

private:
    std::vector<Patch *> patches_;

    //std::vector<Mesh *> meshes_;
    std::unordered_map<std::string, Mesh *> meshes_;

    std::vector<DirichletCondition *> dirichletConditions_;

    std::vector<DirichletConditionFE *> dirichletConditionsFE_;

    std::vector<int> inactiveCPandNode_;

    std::vector<NeumannCondition *> neumannConditions_;

    std::vector<NeumannConditionFE *> neumannConditionsFE_;

    std::vector<InactiveCP *> inactiveCP_;

    std::vector<InactiveNode *> inactiveNode_;

    std::vector<Material *> materials_;

    std::string planeState_;

    std::string problemType_;

    double deltat_;

    double gamma_;

    double beta_;

    double blendZoneThickness_;

    bounded_vector<double, 2> shapeForces_;

    idx_t *cellPartition_;

    idx_t *pointsPartition_;

    idx_t *elementPartition_;

    idx_t *nodePartition_;

    int numberOfSteps_;

    int maximumOfIteration_;

    double tolerance_;

    std::vector<Cell *> cells_;

    std::vector<Element *> elements_;

    std::vector<Cell *> cells_part;

    std::vector<Element *> elements_part;

    std::vector<ControlPoint *> controlPoints_;

    std::vector<Node *> nodes_;

    //std::vector<BoundaryElement *> boundaryFE_;

    int orderParaview_;

    int quadrature_;

    int cpnumber_ = 0;

    int cpaux_ = 0;

    int hammerPoints_; //number of hammer points for elements outside the blend zone

    int hammerPointsBlendZone_; //number of hammer points for elements inside the blend zone

    //std::vector<std::vector<BoundaryElement *>> boundary_;

    std::unordered_map<std::string, std::vector<BoundaryElement *>> finiteElementBoundary_;

    std::vector<BoundaryElement *> blendingBoundary_;

    //Geometry *geometry_;

    std::string finiteElementType_;

    std::string current_working_dir_;
};