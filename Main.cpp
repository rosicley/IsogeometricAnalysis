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

    double lcar1 = 0.15;
    double lcar0 = 0.15;
    double radius = 0.3 / 2.0;
    Point *p0 = solid->addPoint({0.0, 0.0}, lcar0);
    Point *p1 = solid->addPoint({2.0, 0.0}, lcar0);
    Point *p1c = solid->addPoint({2.0, 1.0}, lcar0);
    Point *p2 = solid->addPoint({2.0, 2.0}, lcar0);
    Point *p3 = solid->addPoint({0.0, 2.0}, lcar0);

    Point *p4 = solid->addPoint({0.5, 0.0}, lcar0 / 1.5);
    Point *p5 = solid->addPoint({0.5, 0.25}, 0.0005);

    Point *pc0 = solid->addPoint({0.5, 1.5}, lcar1, false);
    Point *pc1 = solid->addPoint({0.5 + radius, 1.5}, lcar1);
    Point *pc2 = solid->addPoint({0.5, 1.5 + radius}, lcar1);
    Point *pc3 = solid->addPoint({0.5 - radius, 1.5}, lcar1);
    Point *pc4 = solid->addPoint({0.5, 1.5 - radius}, lcar1);

    Point *pd0 = solid->addPoint({1.0, 1.0}, lcar1, false);
    Point *pd2 = solid->addPoint({1.0, 1.0 + radius}, lcar1);
    Point *pd3 = solid->addPoint({1.0 - radius, 1.0}, lcar1);
    Point *pd1 = solid->addPoint({1.0 + radius, 1.0}, lcar1);
    Point *pd4 = solid->addPoint({1.0, 1.0 - radius}, lcar1);

    Point *pe0 = solid->addPoint({1.5, 0.5}, lcar1, false);
    Point *pe2 = solid->addPoint({1.5, 0.5 + radius}, lcar1);
    Point *pe3 = solid->addPoint({1.5 - radius, 0.5}, lcar1);
    Point *pe1 = solid->addPoint({1.5 + radius, 0.5}, lcar1);
    Point *pe4 = solid->addPoint({1.5, 0.5 - radius}, lcar1);

    Line *l0 = solid->addLine({p0, p4});
    Line *l0a = solid->addLine({p4, p1});
    Line *l1a = solid->addLine({p1, p1c});
    Line *l1 = solid->addLine({p1c, p2});
    Line *l2 = solid->addLine({p2, p3});
    Line *l3 = solid->addLine({p3, p0});

    Line *lc0 = solid->addCircle({pc1, pc0, pc2});
    Line *lc1 = solid->addCircle({pc2, pc0, pc3});
    Line *lc2 = solid->addCircle({pc3, pc0, pc4});
    Line *lc3 = solid->addCircle({pc4, pc0, pc1});

    Line *ld0 = solid->addCircle({pd1, pd0, pd2});
    Line *ld1 = solid->addCircle({pd2, pd0, pd3});
    Line *ld2 = solid->addCircle({pd3, pd0, pd4});
    Line *ld3 = solid->addCircle({pd4, pd0, pd1});

    Line *le0 = solid->addCircle({pe1, pe0, pe2});
    Line *le1 = solid->addCircle({pe2, pe0, pe3});
    Line *le2 = solid->addCircle({pe3, pe0, pe4});
    Line *le3 = solid->addCircle({pe4, pe0, pe1});

    LineLoop *ll0 = solid->addLineLoop({l0, l0a, l1a, l1, l2, l3, lc0, lc1, lc2, lc3, ld0, ld1, ld2, ld3, le0, le1, le2, le3});

    PlaneSurface *s0 = solid->addPlaneSurface({ll0});

    solid->addCrack({p4, p5}, s0, "first", 0.035, 0.1);

    solid->addDirichletCondition(l3, {0.0}, {0.0});

    solid->addNeumannCondition(p1c, {10.0}, {});

    // solid->addNeumannCondition(l2, {}, {20.0e3});

    GlobalSolid *problem = new GlobalSolid;

    problem->dataReading("parameters.txt", "properties.txt", "nnnnn.txt", true);

    // if (rank == 0)
    //     problem->teste();

    problem->generateMesh(solid, "T10", "AUTO", "initial", true, false);

    problem->setParametersOfCrackPropagation(0.05, 1.0e-30, true, false, true);

    // problem->addBlendingZone({l1, l2, l4}, 7.5 / 100.0 * offset, 5, true);

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

    // // bounded_vector<double, 2> sifs;

    // // sifs = problem->getSIFs();

    // // std::cout << sifs(0) << " " << sifs(1) << std::endl;

    // // file += std::to_string(sifs(0)) + " " + std::to_string(sifs(1)) + "\n";

    // boost::posix_time::ptime end =
    //     boost::posix_time::microsec_clock::local_time();

    // if (rank == 0)
    // {
    //     boost::posix_time::time_duration diff = end - initial;
    //     std::cout << "  TOTAL (s) = " << std::fixed
    //               << diff.total_milliseconds() / 1000. << std::endl;

    //     std::cout << file << std::endl;
    // }

    PetscFinalize();

    return 0;
}
