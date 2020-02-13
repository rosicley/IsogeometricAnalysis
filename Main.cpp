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

    GlobalSolid *problem = new GlobalSolid;

    // MPI_Barrier(PETSC_COMM_WORLD);

    // boost::posix_time::ptime initialRead =
    //     boost::posix_time::microsec_clock::local_time();

    problem->dataReading("parameters.txt", "properties.txt", "circunferencia_2_4.txt",true);

    // boost::posix_time::ptime endRead =
    //     boost::posix_time::microsec_clock::local_time();

    //     if (rank == 0)
    // {
    //     boost::posix_time::time_duration diff = endRead - initialRead;
    //     std::cout << "   " << std::fixed
    //               << diff.total_milliseconds() / 1000. << std::endl;
    // }

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

    boost::posix_time::ptime end =
        boost::posix_time::microsec_clock::local_time();

    if (rank == 0)
    {
        boost::posix_time::time_duration diff = end - initial;
        std::cout << "  TOTAL (s) = " << std::fixed
                  << diff.total_milliseconds() / 1000. << std::endl;
    }

    PetscFinalize();

    return 0;
}
