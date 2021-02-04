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

    Geometry *solid = new Geometry("GLOBAL"); //LOCAL OR GLOBAL

    Point *p0 = solid->addPoint({10.0, 0.125}, 1.0);
    Point *p1 = solid->addPoint({10.0, 0.0}, 1.0);
    Point *p2 = solid->addPoint({20.0, 0.0}, 1.0);
    Point *p3 = solid->addPoint({20.0, 0.125}, 1.0);

    Line *l0 = solid->addLine({p0, p1});
    Line *l1 = solid->addLine({p1, p2});
    Line *l2 = solid->addLine({p2, p3});
    Line *l3 = solid->addLine({p3, p0});

    LineLoop *ll0 = solid->addLineLoop({l0, l1, l2, l3});

    PlaneSurface *s0 = solid->addPlaneSurface({ll0});

    // Point *p0 = solid->addPoint({10.0, 0.125}, 1.0);
    // Point *p1 = solid->addPoint({0.0, 0.125}, 1.0);
    // Point *p2 = solid->addPoint({0.0, 0.0}, 1.0);
    // Point *p3 = solid->addPoint({20.0, 0.0}, 1.0);
    // Point *p4 = solid->addPoint({20.0, 0.125}, 1.0);

    // Line *l0 = solid->addLine({p0, p1});
    // Line *l1 = solid->addLine({p1, p2});
    // Line *l2 = solid->addLine({p2, p3});
    // Line *l3 = solid->addLine({p3, p4});
    // Line *l4 = solid->addLine({p4, p0});

    // LineLoop *ll0 = solid->addLineLoop({l0, l1, l2, l3, l4});

    // PlaneSurface *s0 = solid->addPlaneSurface({ll0});

    //solid->addCrack({p0, p1}, s0, "first", 0.1, 0.0375);

    double sigma = 500.0;
    // solid->addDirichletCondition(l0, {}, {0.0});
    solid->addDirichletCondition(l2, {0.0}, {0.0});

    // solid->addDirichletCondition(l1, {0.0}, {0.0});
    // solid->addDirichletCondition(l3, {0.0}, {0.0});

    //solid->addDirichletCondition(l3, {0.0}, {});
    //solid->addDirichletCondition(l4, {0.0}, {});

    solid->addNeumannCondition(p0, {}, {-640.0});

    GlobalSolid *problem = new GlobalSolid;

    //problem->useQuarterPointElements();

    problem->dataReading("parameters.txt", "properties.txt", "MondkarISO.txt", true);

    problem->generateMesh(solid, "T10", "AUTO", "initial", true, false);

    //problem->setParametersOfCrackPropagation(0.025, 1.0e-30, true, false, true);

    problem->addBlendingZone({l0, l1, l3}, 0.025, 4, true);

    //problem->addBlendingZoneInLocalMesh(0.15, 4, false);

    if (rank == 0)
    {
        problem->exportMirror();
    }

    std::string planeState = problem->getPlaneState(); //EPD OR EPT
    std::string type = problem->getAnalysisType();     //DYNAMIC OR STATIC

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

    bounded_vector<double, 2> sifs;
    if (rank == 0)
    {
        sifs = problem->getSIFs() / sigma;

        std::cout << sifs(0) << " " << sifs(1) << std::endl;

        file += std::to_string(sifs(0)) + " " + std::to_string(sifs(1)) + "\n";
    }

    delete solid;
    delete problem;

    boost::posix_time::ptime end =
        boost::posix_time::microsec_clock::local_time();

    if (rank == 0)
    {
        boost::posix_time::time_duration diff = end - initial;
        std::cout << "  TOTAL (s) = " << std::fixed
                  << diff.total_milliseconds() / 1000. << std::endl;

        std::cout << file << std::endl;
    }

    PetscFinalize();

    return 0;
}
