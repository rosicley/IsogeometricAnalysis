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

    Geometry *solid = new Geometry("LOCAL"); //LOCAL OR GLOBAL

    // double angle = 0.5235987756;
    // double length = 1.0;

    // Point *p0 = solid->addPoint({5.0 - cos(angle) * length * 0.5, 5.0 - sin(angle) * length * 0.5}, 0.05);
    // Point *p1 = solid->addPoint({5.0 + cos(angle) * length * 0.5, 5.0 + sin(angle) * length * 0.5}, 0.05);

    Point *p0 = solid->addPoint({1.10, 0.0}, 0.025);
    Point *p1 = solid->addPoint({1.10, 0.05}, 0.005);

    solid->addCrackOnGlobal({p0, p1}, "first", 0.02, 0.01, 0.08, 0.025);

    GlobalSolid *problem = new GlobalSolid;

    problem->useQuarterPointElements();

    problem->dataReading("parameters.txt", "properties.txt", "heiderme2.txt", true);
    problem->generateMesh(solid, "T6", "AUTO", "initial", true, false);

    problem->setParametersOfCrackPropagation(0.025, 1.0e-30, true, false, true);

    problem->addBlendingZoneInLocalMesh(0.01, false);

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

    sifs = problem->getSIFs();

    std::cout << sifs(0) / 1.0e05 << " " << sifs(1) / 1.0e05 << std::endl;

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
