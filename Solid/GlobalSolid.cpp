#include "GlobalSolid.h"

GlobalSolid::GlobalSolid() {}

GlobalSolid::GlobalSolid(const std::string &planeState)
{
    planeState_ = planeState;
    deltat_ = 1.0;
    gamma_ = 0.5;
    beta_ = 0.25;
    shapeForces_(0) = 0.0;
    shapeForces_(1) = 0.0;
}

GlobalSolid::~GlobalSolid() {}

Material *GlobalSolid::getMaterial(const int &index)
{
    return materials_[index];
}

Patch *GlobalSolid::getPatch(const int &index)
{
    return patches_[index];
}

void GlobalSolid::addMaterial(const int &index, const double &young, const double &poisson, const double &density)
{
    Material *mat = new Material(index, young, poisson, density);
    materials_.push_back(mat);
}

void GlobalSolid::addDirichletCondition(const int &index, const int &patchIndex, const bounded_vector<int, 2> &free, const bounded_vector<double, 2> &values)
{
    if (free(0) == 1)
    {
        DirichletCondition *cond = new DirichletCondition(patches_[patchIndex]->getControlPoint(index), 0, values(0));
        dirichletConditions_.push_back(cond);
    }
    if (free(1) == 1)
    {
        DirichletCondition *cond = new DirichletCondition(patches_[patchIndex]->getControlPoint(index), 1, values(1));
        dirichletConditions_.push_back(cond);
    }
}

void GlobalSolid::addNeumannCondition(const int &index, const int &patchIndex, const bounded_vector<double, 2> &values)
{
    if (values(0) != 0.0)
    {
        NeumannCondition *cond = new NeumannCondition(patches_[patchIndex]->getControlPoint(index), 0, values(0));
        neumannConditions_.push_back(cond);
    }
    if (values(1) != 0.0)
    {
        NeumannCondition *cond = new NeumannCondition(patches_[patchIndex]->getControlPoint(index), 1, values(1));
        neumannConditions_.push_back(cond);
    }
}

void GlobalSolid::addPatch(const int &index, const int &npc, const int &indexMaterial, const double &thickness)
{
    Patch *patch = new Patch(index, npc, materials_[indexMaterial], thickness);
    patches_.push_back(patch);
}

void GlobalSolid::dataReading(const std::string &inputParameters, const std::string &inputProperties, const std::string &inputMeshIso, const bool &rhino)
{
    std::ifstream parameters(inputParameters);
    std::ifstream properties(inputProperties);
    std::ifstream isogeometric(inputMeshIso);
    std::string line;

    std::stringstream text1;
    text1 << "mirror.txt";
    std::ofstream mirror(text1.str());

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    //READING ANALYSIS PARAMETERS
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    parameters >> problemType_;
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    parameters >> planeState_;
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    parameters >> numberOfSteps_ >> maximumOfIteration_ >> tolerance_;

    if (problemType_ == "DYNAMIC")
    {
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        parameters >> deltat_ >> beta_ >> gamma_;
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        parameters >> quadrature_;
    }
    else
    {
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        std::getline(parameters, line);
        parameters >> quadrature_;
    }

    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    parameters >> orderParaview_;

    mirror << "TYPE OF ANALYSIS: " << problemType_ << std::endl;
    mirror << "PLANE STATE: " << planeState_ << std::endl;
    mirror << "NUMBER OF LOAD/TIME STEP: " << numberOfSteps_ << std::endl;
    mirror << "TOLERANCE: " << tolerance_ << std::endl;
    if (problemType_ == "DYNAMIC")
    {
        mirror << "DELTAT: " << deltat_ << std::endl;
        mirror << "BETA: " << beta_ << std::endl;
        mirror << "GAMMA: " << gamma_ << std::endl;
    }
    mirror << "NUMBER OF INTEGRATION POINTS: " << quadrature_ * quadrature_ << std::endl;
    mirror << "DEGREE(PARAVIEW): " << orderParaview_ << std::endl;
    mirror << std::endl;

    //READING SOLID PROPERTIES
    int nmaterial;
    double young, poisson, density;

    std::getline(properties, line);
    std::getline(properties, line);
    std::getline(properties, line);
    properties >> nmaterial;

    // mirror << "NUMBER OF MATERIALS: " << nmaterial << std::endl;

    for (int i = 0; i < nmaterial; i++)
    {
        std::getline(properties, line);
        std::getline(properties, line);
        std::getline(properties, line);
        std::getline(properties, line);
        std::getline(properties, line);
        properties >> young >> poisson >> density;
        addMaterial(i, young, poisson, density);
        // mirror << "MATERIAL " << i << std::endl;
        // mirror << "YOUNG: " << young << std::endl;
        // mirror << "POISSON: " << poisson << std::endl;
        // mirror << "DENSITY: " << density << std::endl;
    }

    //READING ISOGEOMETRIC MESH

    int npatch;       //nº de patchs
    int cellCont = 0; //contador de células
    std::getline(isogeometric, line);
    std::getline(isogeometric, line);
    std::getline(isogeometric, line);
    isogeometric >> npatch;

    mirror << "NUMBER OF PATCHES: " << npatch << std::endl;

    for (int ipatch = 0; ipatch < npatch; ipatch++)
    {
        int npc, imaterial, dim_u, dim_v;
        double thickness;
        bounded_vector<int, 2> npc_dir, degree; //NUMBER OF CONTROL POINTS(U, V);ORDER(U,V)

        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        isogeometric >> npc;
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        isogeometric >> imaterial;
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        isogeometric >> thickness;
        addPatch(ipatch, npc, imaterial, thickness);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);

        mirror << "PATCH " << ipatch << std::endl;
        mirror << "THICKNESS: " << thickness << std::endl;
        mirror << "YOUNG: " << materials_[imaterial]->getYoung() << std::endl;
        mirror << "POISSON: " << materials_[imaterial]->getPoisson() << std::endl;
        mirror << "DENSITY: " << materials_[imaterial]->getDensity() << std::endl;

        mirror << "CONTROL POINTS (" << npc << ")" << std::endl;
        mirror << "INDEX  COORDINATES (X1, X2)  WEIGHT" << std::endl;

        //CREATE CONTROL POINTS IN PATCH
        matrix<double> auxil(npc, 3);
        if (rhino == true)
        {
            mirror << "///CONTROL POINTS FROM RHINO///" << std::endl;
            for (int ipc = 0; ipc < npc; ipc++)
            {
                double aux;
                isogeometric >> auxil(ipc, 0) >> auxil(ipc, 1) >> aux >> auxil(ipc, 2);

                std::getline(isogeometric, line);
            }
        }
        else
        {
            for (int ipc = 0; ipc < npc; ipc++)
            {
                bounded_vector<double, 2> coord;
                double aux, weight;

                isogeometric >> coord(0) >> coord(1) >> aux >> weight;

                patches_[ipatch]->addControlPoint(ipc, coord, weight);

                mirror << ipc << " " << coord(0) << " " << coord(1) << " " << weight << std::endl;

                // if (coord(0) / weight <= -20.0)
                // {
                //     bounded_vector<int, 2> free;
                //     bounded_vector<double, 2> value;
                //     free(0) = 1;
                //     free(1) = 1;
                //     value(0) = 0.0;
                //     value(1) = 0.0;
                //     addDirichletCondition(ipc, ipatch, free, value);
                // }
                // if (coord(0) / weight >= 20.0)
                // {
                //     bounded_vector<double, 2> value;
                //     value(0) = 200.0;
                //     value(1) = 0.0;
                //     addNeumannCondition(ipc, ipatch, value);
                // }

                std::getline(isogeometric, line);
            }
        }

        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        isogeometric >> npc_dir(0) >> npc_dir(1);

        if (rhino == true)
        {
            for (int j = 0; j < npc_dir(1); j++)
            {
                for (int i = 0; i < npc_dir(0); i++)
                {
                    bounded_vector<double, 2> coord;
                    double aux, weight;

                    coord(0) = auxil(i * npc_dir(1) + j, 0);
                    coord(1) = auxil(i * npc_dir(1) + j, 1);
                    weight = auxil(i * npc_dir(1) + j, 2);

                    patches_[ipatch]->addControlPoint(j * npc_dir(0) + i, coord, weight);

                    mirror << j * npc_dir(0) + i << " " << coord(0) << " " << coord(1) << " " << weight << std::endl;

                    if (coord(0) / weight == 0.0)
                    {
                        bounded_vector<int, 2> free;
                        bounded_vector<double, 2> value;
                        free(0) = 1;
                        free(1) = 1;
                        value(0) = 0.0;
                        value(1) = 0.0;
                        addDirichletCondition(j * npc_dir(0) + i, ipatch, free, value);
                    }
                    if (coord(0) / weight == 500.0)
                    {
                        bounded_vector<double, 2> value;
                        value(0) = 0.0;
                        value(1) = -4000.0;
                        addNeumannCondition(j * npc_dir(0) + i, ipatch, value);
                    }
                }
            }
        }

        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        isogeometric >> degree(0) >> degree(1);

        patches_[ipatch]->setDegree(degree);
        patches_[ipatch]->setNpc_Dir(npc_dir);
        dim_u = npc_dir(0) + degree(0) + 1;
        dim_v = npc_dir(1) + degree(1) + 1;
        vector<double> uknot(dim_u);
        vector<double> vknot(dim_v);

        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);

        for (int iu = 0; iu < dim_u; iu++)
        {
            isogeometric >> uknot(iu);
            std::getline(isogeometric, line);
        }

        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);

        for (int iv = 0; iv < dim_v; iv++)
        {
            isogeometric >> vknot(iv);
            std::getline(isogeometric, line);
        }

        patches_[ipatch]->setKnotsVectors(uknot, vknot);

        //CONDIÇÕES DE CONTORNO RELACIONADAS COM O PATCH 0

        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);

        int nneuman, ndirichlet, ipc;
        bounded_vector<double, 2> value;
        bounded_vector<int, 2> free;

        isogeometric >> nneuman;

        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);

        for (int ineuman = 0; ineuman < nneuman; ineuman++)
        {
            isogeometric >> ipc >> value(0) >> value(1);
            addNeumannCondition(ipc, ipatch, value);
            std::getline(isogeometric, line);
        }

        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);

        isogeometric >> ndirichlet;

        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);
        std::getline(isogeometric, line);

        for (int idirichlet = 0; idirichlet < ndirichlet; idirichlet++)
        {
            isogeometric >> ipc >> free(0) >> free(1) >> value(0) >> value(1);
            addDirichletCondition(ipc, ipatch, free, value);
            std::getline(isogeometric, line);
        }

        /////////////////////////////////FIM DE LEITURA DO PATCH//////////////////////////////////
        /////////////////////////////////CRIANDO CÉLULAS DO PATCH/////////////////////////////////

        int ngp = 0;
        bounded_vector<int, 2> INC, span;
        span(0) = 0;
        span(1) = 0;
        std::vector<ControlPoint *> pointsPatch;
        pointsPatch = patches_[ipatch]->getControlPoints();
        ngp = 0;
        int auxU = (degree(0) + 1);
        int auxV = (degree(1) + 1);
        int auxUV = auxU * auxV;

        for (int i = 0; i < npc_dir(0); i++)
        {
            if ((uknot(i)) != (uknot(i + 1)))
            {
                span(0) = span(0) + 1;
            }
        }
        for (int j = 0; j < npc_dir(1); j++)
        {
            if ((vknot(j)) != (vknot(j + 1)))
            {
                span(1) = span(1) + 1;
            }
        }
        patches_[ipatch]->setSpanNumber(span);

        mirror << "CELLS (" << span(0) * span(1) << ")" << std::endl;
        mirror << "INDEX   {CONNECTION}" << std::endl;

        for (int j = 0; j < npc_dir(1); j++)
        {
            for (int i = 0; i < npc_dir(0); i++)
            {
                INC(0) = i;
                INC(1) = j;
                pointsPatch[ngp]->setINC(INC);
                ngp++;

                // finding the elements
                if ((i >= degree(0)) && (j >= degree(1)))
                {
                    if (((uknot(i)) != (uknot(i + 1))) && ((vknot(j)) != (vknot(j + 1))))
                    {

                        vector<double> connect(auxUV);
                        std::vector<ControlPoint *> connection(auxUV);

                        for (int jloc = 0; jloc < auxV; jloc++)
                        {
                            for (int iloc = 0; iloc < auxU; iloc++)
                            {
                                // global function number
                                int ngf = ngp - jloc * npc_dir(0) - iloc - 1; //

                                // local number function
                                int nlf = (auxUV - 1) - jloc * auxU - iloc;

                                connect(nlf) = ngf;
                            }
                        }

                        mirror << cellCont << " { ";

                        for (int z = 0; z < auxUV; z++)
                        {
                            connection[z] = pointsPatch[connect(z)];
                            mirror << connect(z) << " ";
                        }
                        mirror << "}" << std::endl;

                        Cell *cell = new Cell(cellCont++, patches_[ipatch], connection);
                        cells_.push_back(cell);

                        // for (int k = 0; k < 9; k++)
                        // {

                        //     nodes_[connect(k)]->pushInverseIncidence(index);
                        // };
                    }
                }
            }
        }
    }

    mirror << "DIRICHLET CONDITION (" << dirichletConditions_.size() << ")" << std::endl;
    mirror << "INDEX POINT   DIRECTION   VALUE" << std::endl;
    for (DirichletCondition *dir : dirichletConditions_)
    {
        mirror << dir->getControlPoint()->getIndex() << " " << dir->getDirection() << " " << dir->getValue() << std::endl;
    }

    mirror << "NEUMAN CONDITION (" << neumannConditions_.size() << ")" << std::endl;
    mirror << "INDEX POINT   DIRECTION   VALUE" << std::endl;
    for (NeumannCondition *dir : neumannConditions_)
    {
        mirror << dir->getControlPoint()->getIndex() << " " << dir->getDirection() << " " << dir->getValue() << std::endl;
    }

    parameters.close();
    properties.close();
    isogeometric.close();

    ///CRIAR LOOPING PARA ORGANIZAR OS PONTOS DE CONTROLE DOS PATCHES NO GLOBAL
    controlPoints_ = patches_[0]->getControlPoints();

    for (Patch *pat : patches_)
    {
        pat->removeControlPoints();
    }

    ISOdomainDecompositionMETIS();

    for (int i = 0; i < cells_.size(); i++)
    {
        if (rank == cellPartition_[i])
        {
            cells_part.push_back(cells_[i]);
        }
        else if (rank != 0)
        {

            delete cells_[i];
        }
    }

    if (rank != 0)
    {
        cells_.erase(cells_.begin(), cells_.begin() + cells_.size());
    }
}

std::string GlobalSolid::getAnalysisType()
{
    return problemType_;
}

std::string GlobalSolid::getPlaneState()
{
    return planeState_;
}

int GlobalSolid::solveStaticProblem()
{
    Mat A;
    Vec b, x, All;
    PetscErrorCode ierr;
    PetscInt Istart, Iend, Idof, Ione, iterations, *dof;
    KSP ksp;
    PC pc;
    VecScatter ctx;
    PetscScalar val, value;

    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    //std::stringstream text1;

    if (rank == 0)
    {
        exportToParaview(0);
        // text1 << "ForcasInternas.txt";
    }
    //std::ofstream file1(text1.str());

    double initialNorm = 0.0;
    for (ControlPoint *node : controlPoints_)
    {
        double x1 = node->getInitialCoordinate()(0);
        double x2 = node->getInitialCoordinate()(1);
        initialNorm += x1 * x1 + x2 * x2;
    }

    PetscMalloc1(dirichletConditions_.size(), &dof);
    int idir = 0;
    for (DirichletCondition *dir : dirichletConditions_)
    {
        int indexNode = dir->getControlPoint()->getIndex();
        int direction = dir->getDirection();
        dof[idir++] = (2 * indexNode + direction);
    }

    int nnnumber = controlPoints_.size();

    for (int loadStep = 1; loadStep <= numberOfSteps_; loadStep++)
    {
        boost::posix_time::ptime t1 =
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0)
        {
            std::cout << "------------------------- LOAD STEP = "
                      << loadStep << " -------------------------\n";
        }
        //double norm = 100.0;
        for (int iteration = 0; iteration < maximumOfIteration_; iteration++) //definir o máximo de interações por passo de carga
        {
            //Create PETSc sparse parallel matrix
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                2 * nnnumber, 2 * nnnumber,
                                100, NULL, 300, NULL, &A);
            CHKERRQ(ierr);

            ierr = MatGetOwnershipRange(A, &Istart, &Iend);
            CHKERRQ(ierr);

            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b);
            CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, 2 * nnnumber);
            CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);
            CHKERRQ(ierr);
            ierr = VecDuplicate(b, &x);
            CHKERRQ(ierr);
            ierr = VecDuplicate(b, &All);
            CHKERRQ(ierr);

            if (rank == 0)
            {
                for (NeumannCondition *con : neumannConditions_)
                {
                    int ind = con->getControlPoint()->getIndex();
                    int dir = con->getDirection();
                    double val1 = con->getValue() * (1.0 * loadStep / (1.0 * numberOfSteps_));
                    int dof = 2 * ind + dir;
                    ierr = VecSetValues(b, 1, &dof, &val1, ADD_VALUES);
                }
            }

            if (iteration == 0)
            {
                for (DirichletCondition *con : dirichletConditions_)
                {
                    ControlPoint *cp = con->getControlPoint();
                    int dir = con->getDirection();
                    double val1 = (con->getValue()) / (1.0 * numberOfSteps_);

                    cp->incrementCurrentCoordinate(dir, val1);
                }
            }

            // if (rank == 0)
            // {
            for (Cell *el : cells_part)
            {
                std::pair<vector<double>, matrix<double>> elementMatrices;
                elementMatrices = el->cellContributions(planeState_, "STATIC", loadStep, numberOfSteps_, 1.0, 0.25, 0.5, quadrature_);
                int num = el->getControlPoints().size();

                for (size_t i = 0; i < num; i++)
                {
                    if (fabs(elementMatrices.first(2 * i)) >= 1.0e-11)
                    {
                        int dof = 2 * el->getControlPoints()[i]->getIndex();
                        ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
                    }
                    if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-11)
                    {
                        int dof = 2 * el->getControlPoints()[i]->getIndex() + 1;
                        ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
                    }

                    for (size_t j = 0; j < num; j++)
                    {
                        if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-11)
                        {
                            int dof1 = 2 * el->getControlPoints()[i]->getIndex();
                            int dof2 = 2 * el->getControlPoints()[j]->getIndex();
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-11)
                        {
                            int dof1 = 2 * el->getControlPoints()[i]->getIndex() + 1;
                            int dof2 = 2 * el->getControlPoints()[j]->getIndex();
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-11)
                        {
                            int dof1 = 2 * el->getControlPoints()[i]->getIndex();
                            int dof2 = 2 * el->getControlPoints()[j]->getIndex() + 1;
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-11)
                        {
                            int dof1 = 2 * el->getControlPoints()[i]->getIndex() + 1;
                            int dof2 = 2 * el->getControlPoints()[j]->getIndex() + 1;
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
                        }

                        //std::cout<<elementMatrices.second(2 * i, 2 * j)<<" "<<elementMatrices.second(2 * i + 1, 2 * j)<<" "<<elementMatrices.second(2 * i, 2 * j +1)<< " "<< elementMatrices.second(2 * i+1, 2 * j +1)  <<std::endl;
                    }
                }

                // }
            }

            //MPI_Barrier(PETSC_COMM_WORLD);

            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);
            ierr = VecAssemblyBegin(b);
            CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);
            CHKERRQ(ierr);

            // std::cout<<"MATRIZ A"<<std::endl;
            // MPI_Barrier(PETSC_COMM_WORLD);
            //             std::cout<<" "<<std::endl;

            // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            // std::cout<<" "<<std::endl;

            //             std::cout<<" "<<std::endl;

            // //             MPI_Barrier(PETSC_COMM_WORLD);

            // VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

            //             std::cout<<" "<<std::endl;

            // std::cout<<"VECTOR X"<<std::endl;
            // std::cout << std::endl;
            // VecView(b, PETSC_VIEWER_STDOUT_WORLD);
            // CHKERRQ(ierr);
            // std::cout << std::endl;
            //   std::cout<<std::endl;

            // VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

            MatZeroRowsColumns(A, dirichletConditions_.size(), dof, 1.0, x, b);

            // std::cout<<"MATRIZ A"<<std::endl;

            // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            // std::cout<<std::endl;
            // for (DirichletCondition *dir : dirichletConditions_)
            // {
            //     PetscInt teste;
            //     teste = 2 * dir->getControlPoint()->getIndex() + dir->getDirection();
            //     MatZeroRowsColumns(A, 1, dof, 1.0, x, b);
            // }

            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
            CHKERRQ(ierr);
            ierr = KSPSetOperators(ksp, A, A);
            CHKERRQ(ierr);

            //Solve using MUMPS
#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp, KSPPREONLY);
            ierr = KSPGetPC(ksp, &pc);
            ierr = PCSetType(pc, PCLU);
#endif
            ierr = KSPSetFromOptions(ksp);
            CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);

            //Solve linear system
            ierr = KSPSolve(ksp, b, x);
            CHKERRQ(ierr);
            ierr = KSPGetTotalIterations(ksp, &iterations);

            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(x, &ctx, &All);
            CHKERRQ(ierr);
            ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            ierr = VecScatterDestroy(&ctx);
            CHKERRQ(ierr);

            //Updates nodal variables
            double norm = 0.0;
            Ione = 1;

            for (size_t i = 0; i < nnnumber; i++)
            {
                Idof = 2 * i;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                controlPoints_[i]->incrementCurrentCoordinate(0, val);

                Idof = 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                controlPoints_[i]->incrementCurrentCoordinate(1, val);
            }

            boost::posix_time::ptime t2 =
                boost::posix_time::microsec_clock::local_time();

            if (rank == 0)
            {
                boost::posix_time::time_duration diff = t2 - t1;
                std::cout << "Iteration = " << iteration
                          << " (" << loadStep << ")"
                          << "   x Norm = " << std::scientific << sqrt(norm / initialNorm)
                          << "  Time (s) = " << std::fixed
                          << diff.total_milliseconds() / 1000. << std::endl;
            }

            ierr = KSPDestroy(&ksp);
            CHKERRQ(ierr);
            ierr = VecDestroy(&b);
            CHKERRQ(ierr);
            ierr = VecDestroy(&x);
            CHKERRQ(ierr);
            ierr = VecDestroy(&All);
            CHKERRQ(ierr);
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);

            if (sqrt(norm / initialNorm) <= tolerance_)
            {
                break;
            }
        }

        if (rank == 0)
        {
            // for (Node *n : nodes_)
            // {
            // n->setZeroStressState();
            // }

            // for (int i = 0; i < elements_.size(); i++)
            // {
            // elements_[i]->StressCalculate(planeState_);
            // }
            exportToParaview(loadStep);
            //file1 << (1.0 * loadStep) / (1.0 * numberOfSteps) << " " << -nodes_[0]->getCurrentCoordinate()[1] + nodes_[0]->getInitialCoordinate()[1] << std::endl;
        }
    }
    PetscFree(dof);
    return 0;
}

void GlobalSolid::exportToParaview(const int &loadstep)
{
    matrix<double> qxsi2 = coordinatesForInterpolation(orderParaview_);

    int auxxx = (orderParaview_ + 1) * (orderParaview_ + 1);
    std::stringstream text;
    text << "output" << loadstep << ".vtu";
    std::ofstream file(text.str());

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << auxxx * cells_.size()
         << "\"  NumberOfCells=\"" << cells_.size()
         << "\">"
         << "\n";
    //nodal coordinates
    file << "    <Points>"
         << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">"
         << "\n";

    for (Cell *cell : cells_)
    {

        bounded_vector<double, 2> x, qxsi;

        matrix<double> xInt(auxxx, 2, 0.0);
        bounded_vector<int, 2> INC_;
        std::vector<ControlPoint *> connection = cell->getControlPoints();
        vector<double> wpc2(connection.size());
        vector<double> phi2_;

        INC_ = connection[connection.size() - 1]->getINC();

        for (int i = 0; i < connection.size(); i++)
        {
            wpc2(i) = connection[i]->getWeight();
        }

        for (int j = 0; j < auxxx; j++)
        {
            qxsi(0) = qxsi2(j, 0);
            qxsi(1) = qxsi2(j, 1);
            for (int i = 0; i < connection.size(); i++)
            {
                x = connection[i]->getCurrentCoordinate();

                phi2_ = cell->shapeFunction(qxsi, wpc2, INC_);

                xInt(j, 0) += phi2_(i) * x(0) / wpc2(i);
                xInt(j, 1) += phi2_(i) * x(1) / wpc2(i);
            }
        }

        for (int i = 0; i < auxxx; i++)
        {
            file << xInt(i, 0) << " " << xInt(i, 1) << " " << 0.0 << std::endl;
        }
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Points>"
         << "\n";

    //element connectivity
    file << "    <Cells>"
         << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">"
         << "\n";

    int Ncp = 0;
    for (Cell *cell : cells_)
    {
        for (int i = 0; i < auxxx; i++)
        {
            file << Ncp++ << " ";
        }
        file << std::endl;
    }

    file << "      </DataArray>"
         << "\n";
    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";
    int aux = 0;
    for (Cell *cell : cells_)
    {
        file << aux + auxxx << std::endl;
        aux += auxxx;
    }

    file << "      </DataArray>"
         << "\n";
    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";

    for (Cell *cell : cells_)
    {
        file << 70 << "\n";
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Cells>"
         << "\n";

    //nodal results
    file << "    <PointData>"
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"Displacement\" format=\"ascii\">"
         << "\n";

    for (Cell *cell : cells_)
    {
        //matrix<double> qxsi2 = coordinatesForInterpolation(orderParaview_);
        bounded_vector<double, 2> x, qxsi;
        matrix<double> xInt(auxxx, 2, 0.0);
        bounded_vector<int, 2> INC_;
        std::vector<ControlPoint *> connection = cell->getControlPoints();
        vector<double> wpc2(connection.size());
        vector<double> phi2_;

        INC_ = connection[connection.size() - 1]->getINC();

        for (int i = 0; i < connection.size(); i++)
        {
            wpc2(i) = connection[i]->getWeight();
        }

        for (int j = 0; j < auxxx; j++)
        {
            qxsi(0) = qxsi2(j, 0);
            qxsi(1) = qxsi2(j, 1);
            for (int i = 0; i < connection.size(); i++)
            {
                x = connection[i]->getCurrentDisplacement();

                phi2_ = cell->shapeFunction(qxsi, wpc2, INC_);

                xInt(j, 0) += phi2_(i) * x(0);
                xInt(j, 1) += phi2_(i) * x(1);
            }
        }

        for (int i = 0; i < auxxx; i++)
        {
            file << xInt(i, 0) << " " << xInt(i, 1) << std::endl;
        }
    }

    file << "      </DataArray> "
         << "\n";

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
    //      << "Name=\"Velocity\" format=\"ascii\">"
    //      << "\n";
    // for (Node *n : nodes_)
    // {
    //     bounded_vector<double, 2> currentVelocity = n->getCurrentVelocity();
    //     file << currentVelocity(0) << " " << currentVelocity(1) << "\n";
    // }
    // file << "      </DataArray> "
    //      << "\n";

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
    //      << "Name=\"Acceleration\" format=\"ascii\">"
    //      << "\n";
    // for (Node *n : nodes_)
    // {
    //     file << n->getCurrentAcceleration()(0) << " " << n->getCurrentAcceleration()(1) << "\n";
    // }
    // file << "      </DataArray> "
    //      << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"4\" "
         << "Name=\"CauchyStress\" format=\"ascii\">"
         << "\n";
    for (Cell *cell : cells_)
    {
        bounded_vector<double, 2> qxsi;
        bounded_vector<double, 4> cauchy;

        for (int j = 0; j < auxxx; j++)
        {
            qxsi(0) = qxsi2(j, 0);
            qxsi(1) = qxsi2(j, 1);

            cauchy = cell->getCauchStress(qxsi, planeState_);

            file << cauchy(0) << " " << cauchy(1) << " " << cauchy(2) << " " << cauchy(3) << std::endl;
        }
    }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"4\" "
         << "Name=\"GreenStrain\" format=\"ascii\">"
         << "\n";
    for (Cell *cell : cells_)
    {
        bounded_vector<double, 2> qxsi;
        bounded_vector<double, 4> green;

        for (int j = 0; j < auxxx; j++)
        {
            qxsi(0) = qxsi2(j, 0);
            qxsi(1) = qxsi2(j, 1);

            green = cell->getGreen(qxsi, planeState_);

            file << green(0) << " " << green(1) << " " << green(2) << " " << green(3) << std::endl;
        }
    }
    file << "      </DataArray> "
         << "\n";

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //      << "Name=\"CauchyShearStress\" format=\"ascii\">"
    //      << "\n";
    // for (Node *n : nodes_)
    // {
    //     double cont = n->getStressState()(3);
    //     double aux3 = n->getStressState()(2);
    //     file << aux3 / cont << "\n";
    // }
    // file << "      </DataArray> "
    //      << "\n";

    file << "    </PointData>"
         << "\n";
    //elemental results
    file << "    <CellData>"
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"Process\" format=\"ascii\">" << std::endl;
    for (Cell *el : cells_)
    {
        file << cellPartition_[el->getIndex()] << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //      << "Name=\"ElementType\" format=\"ascii\">" << std::endl;
    // for (Element *el : elements_)
    // {
    //     file << 0 << "\n";
    // }
    // file << "      </DataArray> "
    //      << "\n";

    file << "    </CellData>"
         << "\n";
    //footnote
    file << "  </Piece>"
         << "\n"
         << "  </UnstructuredGrid>"
         << "\n"
         << "</VTKFile>"
         << "\n";
}

matrix<double> GlobalSolid::coordinatesForInterpolation(const int &orderElemement)
{
    int n;
    n = (orderElemement + 1) * (orderElemement + 1);
    matrix<double> coord(n, 2, 0.0);

    if (orderElemement == 1)
    {
        coord(0, 0) = -1.0;
        coord(0, 1) = -1.0;

        coord(1, 0) = 1.0;
        coord(1, 1) = -1.0;

        coord(2, 0) = 1.0;
        coord(2, 1) = 1.0;

        coord(3, 0) = -1.0;
        coord(3, 1) = 1.0;
    }

    if (orderElemement == 2)
    {
        coord(0, 0) = -1.0;
        coord(0, 1) = -1.0;

        coord(1, 0) = 1.0;
        coord(1, 1) = -1.0;

        coord(2, 0) = 1.0;
        coord(2, 1) = 1.0;

        coord(3, 0) = -1.0;
        coord(3, 1) = 1.0;

        coord(4, 0) = 0.0;
        coord(4, 1) = -1.0;

        coord(5, 0) = 1.0;
        coord(5, 1) = 0.0;

        coord(6, 0) = 0.0;
        coord(6, 1) = 1.0;

        coord(7, 0) = -1.0;
        coord(7, 1) = 0.0;

        coord(8, 0) = 0.0;
        coord(8, 1) = 0.0;
    }

    if (orderElemement == 3)
    {
        coord(0, 0) = -1.0;
        coord(0, 1) = -1.0;

        coord(1, 0) = 1.0;
        coord(1, 1) = -1.0;

        coord(2, 0) = 1.0;
        coord(2, 1) = 1.0;

        coord(3, 0) = -1.0;
        coord(3, 1) = 1.0;

        coord(4, 0) = -0.33333333;
        coord(4, 1) = -1.0;

        coord(5, 0) = 0.33333333;
        coord(5, 1) = -1.0;

        coord(6, 0) = 1.0;
        coord(6, 1) = -0.33333333;

        coord(7, 0) = 1.0;
        coord(7, 1) = 0.33333333;

        coord(9, 0) = 0.33333333;
        coord(9, 1) = 1.0;

        coord(8, 0) = -0.33333333;
        coord(8, 1) = 1.0;

        coord(11, 0) = -1.0;
        coord(11, 1) = 0.33333333;

        coord(10, 0) = -1.0;
        coord(10, 1) = -0.33333333;

        coord(12, 0) = -0.33333333;
        coord(12, 1) = -0.33333333;

        coord(13, 0) = 0.33333333;
        coord(13, 1) = -0.33333333;

        coord(15, 0) = 0.33333333;
        coord(15, 1) = 0.33333333;

        coord(14, 0) = -0.33333333;
        coord(14, 1) = 0.33333333;
    }

    return coord;
}

void GlobalSolid::ISOdomainDecompositionMETIS()
{
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());

    int size;
    bounded_vector<int, 2> aux = patches_[0]->getDegrees();

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = cells_.size();
    idx_t numNd = controlPoints_.size();
    idx_t ssize = size;
    idx_t one = 1;
    idx_t n;
    n = (aux(0) + 1) * (aux(1) + 1);
    idx_t elem_start[numEl + 1], elem_connec[n * numEl];

    cellPartition_ = new idx_t[numEl];
    pointsPartition_ = new idx_t[numNd];

    for (idx_t i = 0; i < numEl + 1; i++)
    {
        elem_start[i] = n * i;
    }

    for (idx_t jel = 0; jel < numEl; jel++)
    {

        for (idx_t i = 0; i < cells_[jel]->getControlPoints().size(); i++)
        {
            int nodeIndex = cells_[jel]->getControlPoint(i)->getIndex();
            elem_connec[n * jel + i] = nodeIndex;
        }
    }

    //Performs the domain decomposition
    METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec,
                       NULL, NULL, &one, &ssize, NULL, NULL,
                       &objval, cellPartition_, pointsPartition_);

    mirrorData << std::endl
               << "DOMAIN DECOMPOSITION - CELLS" << std::endl;
    for (int i = 0; i < cells_.size(); i++)
    {
        mirrorData << "process = " << cellPartition_[i]
                   << ", cell = " << i << std::endl;
    }

    mirrorData << std::endl
               << "DOMAIN DECOMPOSITION - CONTROL POINTS" << std::endl;
    for (int i = 0; i < controlPoints_.size(); i++)
    {
        mirrorData << "process = " << pointsPartition_[i]
                   << ", control point = " << i << std::endl;
    }
}

int GlobalSolid::solveDynamicProblem()
{
    firstAccelerationCalculation();

    Mat A;
    Vec b, x, All;
    PetscErrorCode ierr;
    PetscInt Istart, Iend, Idof, Ione, iterations, *dof;
    KSP ksp;
    PC pc;
    VecScatter ctx;
    PetscScalar val, value;

    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    //int n = (order_ + 1) * (order_ + 2) / 2.0;
    //std::stringstream text1;

    if (rank == 0)
    {
        exportToParaview(0);
        //text1 << "DeslocamentoxTempo.txt";
    }
    //std::ofstream file1(text1.str());

    double initialNorm = 0.0;
    for (ControlPoint *node : controlPoints_)
    {
        double x1 = node->getInitialCoordinate()(0);
        double x2 = node->getInitialCoordinate()(1);
        initialNorm += x1 * x1 + x2 * x2;
    }

    PetscMalloc1(dirichletConditions_.size(), &dof);
    for (size_t i = 0; i < dirichletConditions_.size(); i++)
    {
        int indexNode = dirichletConditions_[i]->getControlPoint()->getIndex();
        int direction = dirichletConditions_[i]->getDirection();
        dof[i] = (2 * indexNode + direction);
    }

    int nnnumber = controlPoints_.size();

    for (int timeStep = 1; timeStep <= numberOfSteps_; timeStep++)
    {
        boost::posix_time::ptime t1 =
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0)
        {
            std::cout << "------------------------- TIME STEP = "
                      << timeStep << " -------------------------\n";
        }

        //double norm = 100.0;

        for (int iteration = 0; iteration < maximumOfIteration_; iteration++) //definir o máximo de interações por passo de carga
        {
            //Create PETSc sparse parallel matrix
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                2 * nnnumber, 2 * nnnumber,
                                100, NULL, 300, NULL, &A);
            CHKERRQ(ierr);

            ierr = MatGetOwnershipRange(A, &Istart, &Iend);
            CHKERRQ(ierr);

            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b);
            CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, 2 * nnnumber);
            CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);
            CHKERRQ(ierr);
            ierr = VecDuplicate(b, &x);
            CHKERRQ(ierr);
            ierr = VecDuplicate(b, &All);
            CHKERRQ(ierr);

            if (rank == 0)
            {
                for (NeumannCondition *con : neumannConditions_)
                {
                    int ind = con->getControlPoint()->getIndex();
                    int dir = con->getDirection();
                    double val1 = con->getValue(); //AS FORÇAS APLICADAS PODEM VARIAR AO LONGO DO TEMPO
                    int dof = 2 * ind + dir;
                    ierr = VecSetValues(b, 1, &dof, &val1, ADD_VALUES);
                }
            }

            if (iteration == 0)
            {
                for (DirichletCondition *con : dirichletConditions_)
                {
                    ControlPoint *cp = con->getControlPoint();
                    int dir = con->getDirection();
                    double val1 = (con->getValue()) / (1.0 * numberOfSteps_);

                    cp->incrementCurrentCoordinate(dir, val1);
                }
            }

            for (Cell *el : cells_part)
            {
                std::pair<vector<double>, matrix<double>> elementMatrices;
                elementMatrices = el->cellContributions(planeState_, "DYNAMIC", 1, 1, deltat_, beta_, gamma_, quadrature_); //COM 1 E 1 AS FORÇAS DE DOMINIO PERMANCEM CONSTATEM AO LONGO DO TEMPO
                int num = el->getControlPoints().size();

                for (size_t i = 0; i < num; i++)
                {
                    if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
                    {
                        int dof = 2 * el->getControlPoint(i)->getIndex();
                        ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
                    }
                    if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
                    {
                        int dof = 2 * el->getControlPoint(i)->getIndex() + 1;
                        ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
                    }

                    for (size_t j = 0; j < num; j++)
                    {
                        if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-15)
                        {
                            int dof1 = 2 * el->getControlPoint(i)->getIndex();
                            int dof2 = 2 * el->getControlPoint(j)->getIndex();
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-15)
                        {
                            int dof1 = 2 * el->getControlPoint(i)->getIndex() + 1;
                            int dof2 = 2 * el->getControlPoint(j)->getIndex();
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-15)
                        {
                            int dof1 = 2 * el->getControlPoint(i)->getIndex();
                            int dof2 = 2 * el->getControlPoint(j)->getIndex() + 1;
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-15)
                        {
                            int dof1 = 2 * el->getControlPoint(i)->getIndex() + 1;
                            int dof2 = 2 * el->getControlPoint(j)->getIndex() + 1;
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
                        }
                    }
                }
            }

            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
            CHKERRQ(ierr);

            ierr = VecAssemblyBegin(b);
            CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);
            CHKERRQ(ierr);

            MatZeroRowsColumns(A, dirichletConditions_.size(), dof, 1.0, x, b);

            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
            CHKERRQ(ierr);
            ierr = KSPSetOperators(ksp, A, A);
            CHKERRQ(ierr);

            //Solve using MUMPS
#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp, KSPPREONLY);
            ierr = KSPGetPC(ksp, &pc);
            ierr = PCSetType(pc, PCLU);
#endif
            ierr = KSPSetFromOptions(ksp);
            CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);

            //Solve linear system
            ierr = KSPSolve(ksp, b, x);
            CHKERRQ(ierr);
            ierr = KSPGetTotalIterations(ksp, &iterations);

            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(x, &ctx, &All);
            CHKERRQ(ierr);
            ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            ierr = VecScatterDestroy(&ctx);
            CHKERRQ(ierr);

            // VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

            //Updates nodal variables
            double norm = 0.0;
            Ione = 1;

            for (ControlPoint *cp : controlPoints_)
            {
                int i = cp->getIndex();

                Idof = 2 * i;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                cp->incrementCurrentCoordinate(0, val);

                Idof = 2 * i + 1;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                cp->incrementCurrentCoordinate(1, val);

                bounded_vector<double, 2> vel, accel;
                accel = cp->getCurrentCoordinate() / (beta_ * deltat_ * deltat_) - cp->getPastCoordinate() / (beta_ * deltat_ * deltat_) -
                        cp->getPastVelocity() / (beta_ * deltat_) - cp->getPastAcceleration() * (0.5 / beta_ - 1.0);
                cp->setCurrentAcceleration(accel);

                vel = gamma_ * deltat_ * cp->getCurrentAcceleration() + cp->getPastVelocity() + deltat_ * (1.0 - gamma_) * cp->getPastAcceleration();
                cp->setCurrentVelocity(vel);
            }

            boost::posix_time::ptime t2 =
                boost::posix_time::microsec_clock::local_time();

            if (rank == 0)
            {
                boost::posix_time::time_duration diff = t2 - t1;
                std::cout << "Iteration = " << iteration
                          << " (" << timeStep << ")"
                          << "   x Norm = " << std::scientific << sqrt(norm / initialNorm)
                          << "  Time (s) = " << std::fixed
                          << diff.total_milliseconds() / 1000. << std::endl;
            }

            ierr = KSPDestroy(&ksp);
            CHKERRQ(ierr);
            ierr = VecDestroy(&b);
            CHKERRQ(ierr);
            ierr = VecDestroy(&x);
            CHKERRQ(ierr);
            ierr = VecDestroy(&All);
            CHKERRQ(ierr);
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);

            if (sqrt(norm / initialNorm) <= tolerance_)
            {
                break;
            }
        }

        if (rank == 0)
        {
            exportToParaview(timeStep);
            // file1 << timeStep * deltat_ << " " << nodes_[1]->getCurrentCoordinate()[1] - nodes_[1]->getInitialCoordinate()[1] << std::endl;
        }

        for (ControlPoint *cp : controlPoints_)
        {
            bounded_vector<double, 2> coordinate = cp->getCurrentCoordinate();
            cp->setPastCoordinate(coordinate);
            bounded_vector<double, 2> vel = cp->getCurrentVelocity();
            cp->setPastVelocity(vel);
            bounded_vector<double, 2> accel = cp->getCurrentAcceleration();
            cp->setPastAcceleration(accel);
        }

        for (ControlPoint *cp : controlPoints_)
        {
            bounded_vector<double, 2> vel, accel;
            accel = cp->getCurrentCoordinate() / (beta_ * deltat_ * deltat_) - cp->getPastCoordinate() / (beta_ * deltat_ * deltat_) -
                    cp->getPastVelocity() / (beta_ * deltat_) - cp->getPastAcceleration() * (0.5 / beta_ - 1.0);
            cp->setCurrentAcceleration(accel);
            vel = gamma_ * deltat_ * cp->getCurrentAcceleration() + cp->getPastVelocity() + deltat_ * (1.0 - gamma_) * cp->getPastAcceleration();
            cp->setCurrentVelocity(vel);
        }
    }
    PetscFree(dof);
    return 0;
}

int GlobalSolid::firstAccelerationCalculation()
{
    Mat A;
    Vec b, x, All;
    PetscErrorCode ierr;
    PetscInt Istart, Iend, Idof, Ione, iterations, *dof;
    KSP ksp;
    PC pc;
    VecScatter ctx;
    PetscScalar val, value;

    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int nnnumber = controlPoints_.size();

    //int n = (order_ + 1) * (order_ + 2) / 2.0;

    //Create PETSc sparse parallel matrix
    PetscMalloc1(dirichletConditions_.size(), &dof);
    for (size_t i = 0; i < dirichletConditions_.size(); i++)
    {
        int indexNode = dirichletConditions_[i]->getControlPoint()->getIndex();
        int direction = dirichletConditions_[i]->getDirection();
        dof[i] = (2 * indexNode + direction);
    }
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        2 * nnnumber, 2 * nnnumber,
                        100, NULL, 300, NULL, &A);
    CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, 2 * nnnumber);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &All);
    CHKERRQ(ierr);

    if (rank == 0)
    {
        for (NeumannCondition *con : neumannConditions_)
        {
            int ind = con->getControlPoint()->getIndex();
            int dir = con->getDirection();
            double val1 = con->getValue();
            int dof = 2 * ind + dir;
            ierr = VecSetValues(b, 1, &dof, &val1, ADD_VALUES);
        }
    }

    for (DirichletCondition *con : dirichletConditions_)
    {
        int dir = con->getDirection();
        double val1 = con->getValue();

        con->getControlPoint()->incrementCurrentCoordinate(dir, val1);
    }

    for (Cell *el : cells_part)
    {
        std::pair<vector<double>, matrix<double>> elementMatrices;
        elementMatrices = el->cellContributions(planeState_, "STATIC", 1, 1, deltat_, beta_, gamma_, quadrature_);
        int num = el->getControlPoints().size();
        matrix<double> massLocal(2 * num, 2 * num, 0.0);
        massLocal = el->massMatrix(quadrature_);

        for (size_t i = 0; i < num; i++)
        {
            if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
            {
                int dof = 2 * el->getControlPoint(i)->getIndex();
                ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
            }
            if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
            {
                int dof = 2 * el->getControlPoint(i)->getIndex() + 1;
                ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
            }

            for (size_t j = 0; j < num; j++)
            {
                if (fabs(massLocal(2 * i, 2 * j)) >= 1.e-15)
                {
                    int dof1 = 2 * el->getControlPoint(i)->getIndex();
                    int dof2 = 2 * el->getControlPoint(j)->getIndex();
                    ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j), ADD_VALUES);
                }
                if (fabs(massLocal(2 * i + 1, 2 * j)) >= 1.e-15)
                {
                    int dof1 = 2 * el->getControlPoint(i)->getIndex() + 1;
                    int dof2 = 2 * el->getControlPoint(j)->getIndex();
                    ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i + 1, 2 * j), ADD_VALUES);
                }
                if (fabs(massLocal(2 * i, 2 * j + 1)) >= 1.e-15)
                {
                    int dof1 = 2 * el->getControlPoint(i)->getIndex();
                    int dof2 = 2 * el->getControlPoint(j)->getIndex() + 1;
                    ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j + 1), ADD_VALUES);
                }
                if (fabs(massLocal(2 * i + 1, 2 * j + 1)) >= 1.e-15)
                {
                    int dof1 = 2 * el->getControlPoint(i)->getIndex() + 1;
                    int dof2 = 2 * el->getControlPoint(j)->getIndex() + 1;
                    ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i + 1, 2 * j + 1), ADD_VALUES);
                }
            }
        }
    }

    //Assemble matrices and vectors
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ;

    ierr = VecAssemblyBegin(b);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);
    CHKERRQ(ierr);

    MatZeroRowsColumns(A, dirichletConditions_.size(), dof, 1.0, x, b);

    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);

    //Solve using MUMPS
#if defined(PETSC_HAVE_MUMPS)
    ierr = KSPSetType(ksp, KSPPREONLY);
    ierr = KSPGetPC(ksp, &pc);
    ierr = PCSetType(pc, PCLU);
#endif
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);

    //Solve linear system
    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);
    ierr = KSPGetTotalIterations(ksp, &iterations);

    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(x, &ctx, &All);
    CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);
    CHKERRQ(ierr);

    Ione = 1;

    for (ControlPoint *cp : controlPoints_)
    {
        bounded_vector<double, 2> firstAccel;
        int i = cp->getIndex();
        Idof = 2 * i;
        ierr = VecGetValues(All, Ione, &Idof, &val);
        CHKERRQ(ierr);
        firstAccel(0) = val;

        Idof = 2 * i + 1;
        ierr = VecGetValues(All, Ione, &Idof, &val);
        CHKERRQ(ierr);
        firstAccel(1) = val;
        cp->setCurrentAcceleration(firstAccel);
        cp->setPastAcceleration(firstAccel);
    }

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&All);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    PetscFree(dof);
    if (rank == 0)
    {
        std::cout << "ACELERAÇÕES NO PRIMEIRO PASSO DE TEMPO CALCULADAS." << std::endl;
    }
    return 0;
}
