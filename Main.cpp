static char help[] = "Code to solve static and dynamic nonlinear analysis of 2D elastic solids";

#include "Solid/GlobalSolid.h"

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char *)0, help);

    boost::posix_time::ptime initial =
        boost::posix_time::microsec_clock::local_time();

    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    std::string file;

    clock_t t11 = clock();

    Geometry *geo = new Geometry("GLOBAL"); // LOCAL OR GLOBAL

    double L = 1.0, d = 1.0; //, L2 = 10;
    Point *p0 = geo->addPoint({0.0, 0.0, 0.0});
    Point *p1 = geo->addPoint({L, 0.0, 0.0});
    Point *p2 = geo->addPoint({L, d, 0.0});
    Point *p3 = geo->addPoint({0.0, d, 0.0});
    Line *l0 = geo->addLine({p0, p1});
    Line *l1 = geo->addLine({p1, p2});
    Line *l2 = geo->addLine({p2, p3});
    Line *l3 = geo->addLine({p3, p0});
    LineLoop *ll0 = geo->addLineLoop({l0, l1, l2, l3});
    PlaneSurface *s0 = geo->addPlaneSurface({ll0}, 0, 1.0);

    geo->addDirichletCondition(l3, {0.0}, {});
    geo->addDirichletCondition(p3, {}, {0.0});

    geo->addNeumannCondition(l1, {1.0}, {});

    GlobalSolid *problem = new GlobalSolid;

    problem->dataReading("parameters.txt", "properties.txt", "naoexisteste.txt", true);

    problem->generateMesh(geo, "T10", "AUTO", "initial", true, false);

    if (rank == 0)
    {
        problem->exportMirror();
    }

    std::string planeState = problem->getPlaneState(); // EPD OR EPT
    std::string type = problem->getAnalysisType();     // DYNAMIC OR STATIC

    if (planeState == "EPD" or planeState == "EPT")
    {
        if (type == "STATIC")
        {
            problem->solveStaticProblem();
        }
        else if (type == "DYNAMIC")
        {
            problem->solveDynamicProblem();
        }
        else
        {
            std::cout << "PLEASE SELECT THE APPROPRIATE TYPE OF ANALYSIS (STATIC OR DYNAMIC)" << std::endl;
        }
    }
    else
    {
        std::cout << "PLEASE SELECT THE APPROPRIATE PLANE STATE (EPD OR EPT)" << std::endl;
    }

    problem->error();

    MPI_Barrier(PETSC_COMM_WORLD);

    delete geo;
    delete problem;

    boost::posix_time::ptime end =
        boost::posix_time::microsec_clock::local_time();

    PetscFinalize();

    return 0;
}
