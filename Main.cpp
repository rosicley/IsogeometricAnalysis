static char help[] = "Code to solve static and dynamic nonlinear analysis of 2D elastic solids";

#include "Solid/GlobalSolid.h"

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char *)0, help);

    std::stringstream name;
    name << "SIFs";
    std::ofstream file(name.str());

    for (int i = 0; i <= 90)
    {

        boost::posix_time::ptime initial =
            boost::posix_time::microsec_clock::local_time();

        int rank, size;

        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);

        Geometry *solid = new Geometry("GLOBAL"); //LOCAL OR GLOBAL

        Point *p0 = solid->addPoint({0.0, 0.0}, 1.0);
        Point *p1 = solid->addPoint({10.0, 0.0}, 1.0);
        Point *p2 = solid->addPoint({10.0, 10.0}, 1.0);
        Point *p3 = solid->addPoint({0.0, 10.0}, 1.0);

        const double pi = 3.14159265359;
        Point *p4 = solid->addPoint({4.566987298, 4.75}, 0.05);
        Point *p5 = solid->addPoint({5.433012702, 5.25}, 0.05);

        Point *p6 = solid->addPoint({0.0, 5.0}, 1.0);

        Line *l0 = solid->addLine({p0, p1});
        Line *l1 = solid->addLine({p1, p2});
        Line *l2 = solid->addLine({p2, p3});
        Line *l3 = solid->addLine({p3, p6});
        Line *l4 = solid->addLine({p6, p0});

        LineLoop *ll0 = solid->addLineLoop({l0, l1, l2, l3, l4});

        PlaneSurface *s0 = solid->addPlaneSurface({ll0});

        Crack *c = solid->addCrack({p4, p5}, s0, "none", 0.25, 0.05);

        // solid->transfiniteLine({l0, l2, l3, l4}, 10);
        // solid->transfiniteLine({l1}, 20);
        // solid->transfiniteLine({c0}, 10);
        //solid->transfiniteSurface({s0}, "Alternated");

        solid->addDirichletCondition(p0, {0.0}, {});
        solid->addDirichletCondition(l0, {}, {0.0});
        solid->addNeumannCondition(l2, {}, {1.0});

        GlobalSolid *problem = new GlobalSolid;
        problem->dataReading("parameters.txt", "properties.txt", "none.txt", true);

        problem->generateMesh(solid, "T6", "AUTO", "initial", true, false);

        //problem->addBlending({l0, l2, l3}, 0.025);

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

        if (rank == 0)
        {
            bounded_vector<double, 2> teste = problem->getSIFs();
        }

        boost::posix_time::ptime end =
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0)
        {
            boost::posix_time::time_duration diff = end - initial;
            std::cout << "  TOTAL (s) = " << std::fixed
                      << diff.total_milliseconds() / 1000. << std::endl;
        }

        i += 5;
    }
    PetscFinalize();

    return 0;
}
