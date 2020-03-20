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

void GlobalSolid::addDirichletConditionFE(const int &index, const bounded_vector<int, 2> &free, const bounded_vector<double, 2> &values)
{
    if (free(0) == 1)
    {
        DirichletConditionFE *cond = new DirichletConditionFE(nodes_[index], 0, values(0));
        dirichletConditionsFE_.push_back(cond);
    }
    if (free(1) == 1)
    {
        DirichletConditionFE *cond = new DirichletConditionFE(nodes_[index], 1, values(1));
        dirichletConditionsFE_.push_back(cond);
    }
}

void GlobalSolid::addNeumannConditionFE(const int &index, const bounded_vector<double, 2> &values)
{
    if (values(0) != 0.0)
    {
        NeumannConditionFE *cond = new NeumannConditionFE(nodes_[index], 0, values(0));
        neumannConditionsFE_.push_back(cond);
    }
    if (values(1) != 0.0)
    {
        NeumannConditionFE *cond = new NeumannConditionFE(nodes_[index], 1, values(1));
        neumannConditionsFE_.push_back(cond);
    }
}

void GlobalSolid::addPatch(const int &index, const int &npc, const int &indexMaterial, const double &thickness)
{
    Patch *patch = new Patch(index, npc, materials_[indexMaterial], thickness);
    patches_.push_back(patch);
}

void GlobalSolid::addMesh(const int &index, const int &indexMaterial, const double &thickness, const std::string &elementType)
{
    Mesh *mesh = new Mesh(index, materials_[indexMaterial], thickness, elementType);
    meshes_.push_back(mesh);
}

void GlobalSolid::addNode(const int &index, const int &indexFE,
                          const bounded_vector<double, 2> &initialCoordinate)
{
    Node *node = new Node(index, indexFE, initialCoordinate);
    nodes_.push_back(node);
}

void GlobalSolid::dataReading(const std::string &inputParameters, const std::string &inputProperties, const std::string &inputMeshIso, const std::string &inputMeshFE, const bool &rhino)
{
    std::ifstream parameters(inputParameters);
    std::ifstream properties(inputProperties);
    std::ifstream isogeometric(inputMeshIso);
    std::string line;

    //std::stringstream text1;
    //text1 << "ISOmirror.txt";
    //std::ofstream mirror(text1.str());

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
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    parameters >> blendZoneThickness_;

    //READING SOLID PROPERTIES
    int nmaterial;
    double young, poisson, density;
    std::getline(properties, line);
    std::getline(properties, line);
    std::getline(properties, line);
    properties >> nmaterial;
    for (int i = 0; i < nmaterial; i++)
    {
        std::getline(properties, line);
        std::getline(properties, line);
        std::getline(properties, line);
        std::getline(properties, line);
        std::getline(properties, line);
        properties >> young >> poisson >> density;
        addMaterial(i, young, poisson, density);
    }

    //READING ISOGEOMETRIC MESH
    int npatch;       //nº de patchs
    int cellCont = 0; //contador de células
    int cellpatches = 0;
    std::getline(isogeometric, line);
    std::getline(isogeometric, line);
    std::getline(isogeometric, line);
    isogeometric >> npatch;
    std::getline(isogeometric, line);
    for (int ipatch = 0; ipatch < npatch; ipatch++)
    {
        int npc, imaterial, dim_u, dim_v;
        double thickness;
        bounded_vector<int, 2> npc_dir, degree; //NUMBER OF CONTROL POINTS(U, V);ORDER(U,V)
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

        //CREATE CONTROL POINTS IN PATCH
        matrix<double> auxil(npc, 3);
        if (rhino == true)
        {
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

        //CONDIÇÕES DE CONTORNO RELACIONADAS COM O PATCH ipatch
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

        std::string typecond;
        //std::vector<int> neumanCurve;
        //std::vector<double> neumanValue0;
        //std::vector<double> neumanValue1;

        std::vector<bounded_vector<double, 3>> neummanCurve;

        for (int ineuman = 0; ineuman < nneuman; ineuman++)
        {
            isogeometric >> typecond;
            if (typecond == "POINT")
            {
                isogeometric >> ipc >> value(0) >> value(1);
                addNeumannCondition(ipc, ipatch, value);
                std::getline(isogeometric, line);
            }
            else if (typecond == "CURVE")
            {
                bounded_vector<double, 3> aux;

                isogeometric >> ipc >> aux(1) >> aux(2);
                aux(0) = static_cast<double>(ipc);

                neummanCurve.push_back(aux);
                std::getline(isogeometric, line);
            }
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
            isogeometric >> typecond;
            if (typecond == "POINT")
            {
                isogeometric >> ipc >> free(0) >> free(1) >> value(0) >> value(1);
                addDirichletCondition(ipc, ipatch, free, value);
                std::getline(isogeometric, line);
            }
            else if (typecond == "CURVE")
            {
                isogeometric >> ipc >> free(0) >> free(1) >> value(0) >> value(1);
                if (ipc == 0)
                {
                    for (int cp = 0; cp < npc_dir(0); cp++)
                    {
                        addDirichletCondition(cp, ipatch, free, value);
                    }
                }
                else if (ipc == 1)
                {
                    int aux = npc_dir(0) - 1;
                    for (int ih = 0; ih < npc_dir(1); ih++)
                    {
                        addDirichletCondition(aux, ipatch, free, value);
                        aux = aux + npc_dir(0);
                    }
                }
                else if (ipc == 2)
                {
                    int aux = (npc_dir(1) - 1) * npc_dir(0);
                    for (int ih = 0; ih < npc_dir(0); ih++)
                    {
                        addDirichletCondition(aux, ipatch, free, value);
                        aux = aux + 1;
                    }
                }
                else if (ipc == 3)
                {
                    int aux = 0;
                    for (int ih = 0; ih < npc_dir(1); ih++)
                    {
                        addDirichletCondition(aux, ipatch, free, value);
                        aux = aux + npc_dir(0);
                    }
                }
                std::getline(isogeometric, line);
            }
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

                        for (int z = 0; z < auxUV; z++)
                        {
                            connection[z] = pointsPatch[connect(z)];
                        }

                        Cell *cell = new Cell(cellCont++, patches_[ipatch], connection);
                        cells_.push_back(cell);
                    }
                }
            }
        }
        if (rank == 0)
        {
            for (int i = 0; i < neummanCurve.size(); i++)
            {
                bounded_vector<double, 3> neu = neummanCurve[i];
                value(0) = neu(1);
                value(1) = neu(2);
                int curve = static_cast<int>(neu(0));

                if (curve == 0)
                {
                    int j = cellpatches;
                    vector<double> force(2 * npc_dir(0), 0.0);
                    for (int c = 0; c < span(0); c++)
                    {
                        vector<double> contribuition = cells_[j]->computeDistribuitedLoads(value, quadrature_, 0);
                        j = j + 1;

                        for (int in = 0; in < contribuition.size(); in++)
                        {
                            force(2 * c + in) += contribuition(in);
                        }
                    }
                    for (int ih = 0; ih < npc_dir(0); ih++)
                    {
                        value(0) = force(2 * ih);
                        value(1) = force(2 * ih + 1);
                        addNeumannCondition(ih, ipatch, value);
                    }
                }
                else if (curve == 1)
                {
                    int j = cellpatches + span(0) - 1;
                    vector<double> force(2 * npc_dir(1), 0.0);

                    for (int c = 0; c < span(1); c++)
                    {
                        vector<double> contribuition = cells_[j]->computeDistribuitedLoads(value, quadrature_, 1);
                        j = j + span(0);

                        for (int in = 0; in < contribuition.size(); in++)
                        {
                            force(2 * c + in) += contribuition(in);
                        }
                    }

                    int aux = npc_dir(0) - 1;
                    for (int ih = 0; ih < npc_dir(1); ih++)
                    {
                        value(0) = force(2 * ih);
                        value(1) = force(2 * ih + 1);
                        addNeumannCondition(aux, ipatch, value);
                        aux = aux + npc_dir(0);
                    }
                }
                else if (curve == 2)
                {
                    int j = cellpatches + (span(1) - 1) * span(0);
                    vector<double> force(2 * npc_dir(0), 0.0);
                    for (int c = 0; c < span(0); c++)
                    {
                        vector<double> contribuition = cells_[j]->computeDistribuitedLoads(value, quadrature_, 2);
                        std::vector<ControlPoint *> points = cells_[j]->getControlPointsOnSide(2);
                        j = j + 1;

                        for (int in = 0; in < contribuition.size(); in++)
                        {
                            force(2 * c + in) += contribuition(in);
                        }
                    }
                    int aux = (npc_dir(1) - 1) * npc_dir(0);
                    for (int ih = 0; ih < npc_dir(0); ih++)
                    {
                        value(0) = force(2 * ih);
                        value(1) = force(2 * ih + 1);
                        addNeumannCondition(aux, ipatch, value);
                        aux = aux + 1;
                    }
                }
                else if (curve == 3)
                {
                    int j = cellpatches;
                    vector<double> force(2 * npc_dir(1), 0.0);
                    for (int c = 0; c < span(1); c++)
                    {
                        vector<double> contribuition = cells_[j]->computeDistribuitedLoads(value, quadrature_, 3);
                        j = j + span(0);

                        for (int in = 0; in < contribuition.size(); in++)
                        {
                            force(2 * c + in) += contribuition(in);
                        }
                    }
                    int aux = 0;
                    for (int ih = 0; ih < npc_dir(1); ih++)
                    {
                        value(0) = force(2 * ih);
                        value(1) = force(2 * ih + 1);
                        addNeumannCondition(aux, ipatch, value);
                        aux = aux + npc_dir(0);
                    }
                }
            }
        }
        cellpatches += span(0) * span(1);
    }

    int cpindex = 0;
    double dist = 1.0e-6;
    if (cells_.size() > 0)
    {
        parameters.close();
        properties.close();
        isogeometric.close();

        ///LOOPING PARA ORGANIZAR OS PONTOS DE CONTROLE DOS PATCHES NO GLOBAL
        for (Patch *patch : patches_)
        {
            std::vector<ControlPoint *> cpaux = patch->getControlPoints();

            for (ControlPoint *cp_patch : cpaux)
            {
                bounded_vector<double, 2> coord_patch = cp_patch->getInitialCoordinate();
                int ver = -1;

                for (ControlPoint *cp_global : controlPoints_)
                {
                    bounded_vector<double, 2> coord_global = cp_global->getInitialCoordinate();
                    int index = cp_global->getIndex();

                    bounded_vector<double, 2> aux = coord_global - coord_patch;
                    double dist2 = aux(0) * aux(0) + aux(1) * aux(1);

                    if (dist2 <= dist)
                    {
                        ver = index;
                        break;
                    }
                }
                if (ver != -1)
                {
                    cp_patch->setIndex(ver);
                    controlPoints_.push_back(cp_patch);
                }
                else
                {
                    cp_patch->setIndex(cpindex);
                    controlPoints_.push_back(cp_patch);
                    cpindex = cpindex + 1;
                }
            }
        }

        cpnumber_ = cpindex;
        cpaux_ = cpnumber_;

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
            // else if (rank != 0)
            // {
            //     delete cells_[i];
            // }
        }
        // if (rank != 0)
        // {
        //     cells_.erase(cells_.begin(), cells_.begin() + cells_.size());
        // }
    }

    dataFromGmsh(inputMeshFE);
}

void GlobalSolid::dataFromGmsh(const std::string &inputGmesh)
{
    std::ifstream feMesh(inputGmesh);

    // std::stringstream text1;
    // text1 << "mirror.txt";
    // std::ofstream mirror(text1.str());

    std::string line;
    std::string elementType;
    int nmesh, nnode, nelem, auxEl, auxBo;

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    feMesh >> elementType;
    if (elementType == "T3")
    {
        auxEl = 3;
        auxBo = 2;
    }
    else if (elementType == "T6")
    {
        auxEl = 6;
        auxBo = 3;
    }
    else if (elementType == "T10")
    {
        auxEl = 10;
        auxBo = 4;
    }
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    feMesh >> nmesh;
    std::getline(feMesh, line);
    for (int i = 0; i < nmesh; i++)
    {
        int indexMaterial;
        double thickness;
        std::getline(feMesh, line);
        std::getline(feMesh, line);
        std::getline(feMesh, line);
        std::getline(feMesh, line);
        feMesh >> indexMaterial >> thickness;
        addMesh(i, indexMaterial, thickness, elementType);
        std::getline(feMesh, line);
    }
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    feMesh >> nnode;
    std::getline(feMesh, line);
    int trash;
    bounded_vector<double, 2> coord;
    for (int inode = 0; inode < nnode; inode++)
    {
        feMesh >> trash >> coord(0) >> coord(1);
        addNode(cpnumber_, inode, coord);
        cpnumber_ = cpnumber_ + 1;
        std::getline(feMesh, line);
    }
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    feMesh >> nelem;
    std::getline(feMesh, line);

    int type, num, auxconec;
    int elementCont = 0;
    for (int ielem = 0; ielem < nelem; ielem++)
    {
        feMesh >> trash >> type >> trash >> trash >> num;
        if (type == 1 or type == 8 or type == 26)
        {
            std::vector<Node *> conec;

            for (int in = 0; in < auxBo; in++)
            {
                feMesh >> auxconec;
                conec.push_back(nodes_[auxconec - 1]);
            }

            BoundaryElement *bfe = new BoundaryElement(num, conec);
            boundaryFE_.push_back(bfe);
        }
        else if (type == 2 or type == 9 or type == 21)
        {
            std::vector<Node *> conec;

            for (int in = 0; in < auxEl; in++)
            {
                feMesh >> auxconec;
                conec.push_back(nodes_[auxconec - 1]);
            }

            Element *el = new Element(elementCont++, meshes_[num - 1], conec);
            elements_.push_back(el);
        }
        std::getline(feMesh, line);
    }
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);

    int nneuman, ndirichlet, inode;
    bounded_vector<double, 2> value;
    bounded_vector<int, 2> free;
    std::string typecond;

    feMesh >> nneuman;
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);

    for (int in = 0; in < nneuman; in++)
    {
        feMesh >> typecond;
        if (typecond == "NODE")
        {
            feMesh >> inode >> value(0) >> value(1);
            if (fabs(value(0)) >= 1.0e-10)
            {
                NeumannConditionFE *neu0 = new NeumannConditionFE(nodes_[inode - 1], 0, value(0));
                neumannConditionsFE_.push_back(neu0);
            }
            if (fabs(value(1)) >= 1.0e-10)
            {
                NeumannConditionFE *neu1 = new NeumannConditionFE(nodes_[inode - 1], 1, value(1));
                neumannConditionsFE_.push_back(neu1);
            }
            std::getline(feMesh, line);
        }
        else if (typecond == "LINE")
        {
            feMesh >> inode >> value(0) >> value(1);

            for (BoundaryElement *bound : boundaryFE_)
            {
                if (inode == bound->getBoundaryIndex())
                {
                    vector<double> force = bound->computeDistribuitedLoads(value, quadrature_);
                    std::vector<Node *> conec = bound->getNodes();
                    double value0, value1;
                    for (int i = 0; i < conec.size(); i++)
                    {
                        value0 = force(2 * i);
                        value1 = force(2 * i + 1);
                        if (fabs(value(0)) >= 1.0e-10)
                        {
                            NeumannConditionFE *neu0 = new NeumannConditionFE(conec[i], 0, value0);
                            neumannConditionsFE_.push_back(neu0);
                        }
                        if (fabs(value(1)) >= 1.0e-10)
                        {
                            NeumannConditionFE *neu1 = new NeumannConditionFE(conec[i], 1, value1);
                            neumannConditionsFE_.push_back(neu1);
                        }
                    }
                }
            }
            std::getline(feMesh, line);
        }
    }
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    feMesh >> ndirichlet;

    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);

    for (int id = 0; id < ndirichlet; id++)
    {
        feMesh >> typecond;
        if (typecond == "NODE")
        {
            feMesh >> inode >> free(0) >> free(1) >> value(0) >> value(1);
            if (free(0) != 0)
            {
                DirichletConditionFE *dir = new DirichletConditionFE(nodes_[inode - 1], 0, value(0));
                dirichletConditionsFE_.push_back(dir);
            }
            if (free(1) != 0)
            {
                DirichletConditionFE *dir = new DirichletConditionFE(nodes_[inode - 1], 1, value(1));
                dirichletConditionsFE_.push_back(dir);
            }
            std::getline(feMesh, line);
        }
        else if (typecond == "LINE")
        {
            feMesh >> inode >> free(0) >> free(1) >> value(0) >> value(1);
            std::vector<int> nodexIndex;
            for (BoundaryElement *bound : boundaryFE_)
            {
                if (inode == bound->getBoundaryIndex())
                {
                    std::vector<Node *> conec = bound->getNodes();

                    for (Node *no : conec)
                    {
                        int existe = 0;
                        for (int tes = 0; tes < nodexIndex.size(); tes++)
                        {
                            if (no->getIndexFE() == nodexIndex[tes])
                            {
                                existe = -1;
                            }
                        }
                        if (existe == 0)
                        {
                            nodexIndex.push_back(no->getIndexFE());
                        }
                    }
                }
            }
            for (int i = 0; i < nodexIndex.size(); i++)
            {
                if (free(0) != 0)
                {
                    DirichletConditionFE *dir = new DirichletConditionFE(nodes_[nodexIndex[i]], 0, value(0));
                    dirichletConditionsFE_.push_back(dir);
                }
                if (free(1) != 0)
                {
                    DirichletConditionFE *dir = new DirichletConditionFE(nodes_[nodexIndex[i]], 1, value(1));
                    dirichletConditionsFE_.push_back(dir);
                }
            }
            std::getline(feMesh, line);
        }
    }

    if (rank == 0)
    {
        exportMirror();
    }

    domainDecompositionMETIS(elementType);

    for (int i = 0; i < elements_.size(); i++)
    {
        if (rank == elementPartition_[i])
        {
            elements_part.push_back(elements_[i]);
        }
        // else if (rank != 1) //DEPOIS COLOCAR PARA O RANK 1 TER TODOS!!!
        // {
        //     delete elements_[i];
        // }
    } //DESCOMENTAR AQUI
    // if (rank != 1)
    // {
    //     elements_.erase(elements_.begin(), elements_.begin() + elements_.size());
    // }
}

void GlobalSolid::exportMirror()
{
    std::stringstream name;
    name << "mirror.txt";
    std::ofstream mirror(name.str());

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

    mirror << std::endl;
    mirror << "NUMBER OF FREEDOM DEGREE: " << cpnumber_ << std::endl;
    mirror << "CONTROL POINTS: " << cpaux_ << std::endl;
    mirror << "NODES: " << nodes_.size() << std::endl;
    mirror << std::endl;

    mirror << "NUMBER OF PATCHES: " << patches_.size() << std::endl;
    for (Patch *pat : patches_)
    {
        double young, poisson, densisity;
        pat->getMaterial()->setProperties(young, poisson, densisity);
        mirror << "PATCH: " << pat->getIndex() << std::endl;
        mirror << "THICKNESS: " << pat->getThickness() << std::endl;
        mirror << "YOUNG: " << young << std::endl;
        mirror << "POISSON: " << poisson << std::endl;
        mirror << "DENSITY: " << densisity << std::endl;
        mirror << "NUMBER OF CONTROL POINTS (U,V): " << pat->getNpc_Dir(0) << " " << pat->getNpc_Dir(1) << " (" << pat->getControlPointsNumber() << ")" << std::endl;
        mirror << "NUMBER OF CELLS: " << pat->getNumberOfCells() << std::endl;
        mirror << "DEGREE (U,V): " << pat->getDegree(0) << " " << pat->getDegree(1) << std::endl;
        mirror << std::endl;
    }

    mirror << "NUMBER OF CONTROL POINTS: " << controlPoints_.size() << std::endl;
    mirror << "GLOBAL INDEX   COORDINATES   WEIGHT" << std::endl;
    for (ControlPoint *con : controlPoints_)
    {
        mirror << con->getIndex() << "   " << con->getInitialCoordinate()(0) << " " << con->getInitialCoordinate()(1) << "    " << con->getWeight() << std::endl;
    }
    mirror << std::endl;

    mirror << "NUMBER OF CELLS: " << cells_.size() << std::endl;
    mirror << "INDEX   PATCH   { CONNECTION }" << std::endl;
    for (Cell *cell : cells_)
    {
        std::vector<ControlPoint *> conec = cell->getControlPoints();
        mirror << cell->getIndex() << "    " << cell->getPatch()->getIndex() << "    { ";
        for (ControlPoint *con : conec)
        {
            mirror << con->getIndex() << " ";
        }
        mirror << "}" << std::endl;
    }
    mirror << std::endl;

    mirror << "NUMBER OF DIRICHLET CONDITION(ISO): " << dirichletConditions_.size() << std::endl;
    mirror << "INDEX POINT   DIRECTION   VALUE" << std::endl;
    for (DirichletCondition *dir : dirichletConditions_)
    {
        mirror << dir->getControlPoint()->getIndex() << "    " << dir->getDirection() << "    " << dir->getValue() << std::endl;
    }
    mirror << std::endl;

    mirror << "NUMBER OF NEUMANN CONDITION(ISO): " << neumannConditions_.size() << std::endl;
    mirror << "INDEX POINT   DIRECTION   VALUE" << std::endl;
    for (NeumannCondition *dir : neumannConditions_)
    {
        mirror << dir->getControlPoint()->getIndex() << "    " << dir->getDirection() << "    " << dir->getValue() << std::endl;
    }
    mirror << std::endl;

    //FEM
    mirror << "NUMBER OF SURFACE: " << patches_.size() << std::endl;
    for (Mesh *pat : meshes_)
    {
        double young, poisson, densisity;
        pat->getMaterial()->setProperties(young, poisson, densisity);
        mirror << "SURFACE: " << pat->getIndex() << std::endl;
        mirror << "THICKNESS: " << pat->getThickness() << std::endl;
        mirror << "YOUNG: " << young << std::endl;
        mirror << "POISSON: " << poisson << std::endl;
        mirror << "DENSITY: " << densisity << std::endl;
        mirror << "ELEMENT TYPE: " << pat->getElementType() << std::endl;
        mirror << std::endl;
    }

    mirror << "NUMBER OF NODES: " << nodes_.size() << std::endl;
    mirror << "GLOBAL INDEX    FEM INDEX   COORDINATES" << std::endl;
    for (Node *con : nodes_)
    {
        mirror << con->getIndex() << "   " << con->getIndexFE() << "   " << con->getInitialCoordinate()(0) << " " << con->getInitialCoordinate()(1) << std::endl;
    }
    mirror << std::endl;

    mirror << "NUMBER OF ELEMENTS: " << elements_.size() << std::endl;
    mirror << "INDEX   SURFACE   { CONNECTION }" << std::endl;
    for (Element *el : elements_)
    {
        std::vector<Node *> conec = el->getConnection();
        mirror << el->getIndex() << "    " << el->getMesh()->getIndex() << "    { ";
        for (Node *con : conec)
        {
            mirror << con->getIndex() << " ";
        }
        mirror << "}" << std::endl;
    }
    mirror << std::endl;

    mirror << "NUMBER OF DIRICHLET CONDITION(FEM): " << dirichletConditionsFE_.size() << std::endl;
    mirror << "INDEX NODE   DIRECTION   VALUE" << std::endl;
    for (DirichletConditionFE *dir : dirichletConditionsFE_)
    {
        mirror << dir->getNode()->getIndex() << "    " << dir->getDirection() << "    " << dir->getValue() << std::endl;
    }
    mirror << std::endl;

    mirror << "NUMBER OF NEUMANN CONDITION(FEM): " << neumannConditionsFE_.size() << std::endl;
    mirror << "INDEX NODE   DIRECTION   VALUE" << std::endl;
    for (NeumannConditionFE *dir : neumannConditionsFE_)
    {
        mirror << dir->getNode()->getIndex() << "    " << dir->getDirection() << "    " << dir->getValue() << std::endl;
    }
    mirror << std::endl;
}

void GlobalSolid::teste()
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    // std::cout << "RANK: " << rank << " NODES: " << nodes_.size() << std::endl;
    // std::cout << "RANK: " << rank << " Elements: " << elements_.size() << std::endl;
    // std::cout << "RANK: " << rank << " BoundaryElements: " << boundaryFE_.size() << std::endl;
    // std::cout<<"TESTANDO "<<elements_.size()<<std::endl;

    // std::cout<<"TESTANDO 1 "<<std::endl;
    //computeDistanceFromFEBoundary();
    // std::cout << "teste... " << std::endl;
    // incidenceLocalxGlobal();

    if (rank == 0)
    {
        computeDistanceFromFEBoundary();
        incidenceLocalxGlobal();

        vector<double> mass = diagonalMassMatrix();
        std::cout << "GLOBAL: " << std::endl;
        for (int i = 0; i < controlPoints_.size(); i++)
        {
            std::cout << mass(i) << std::endl;
        }
        std::cout << "LOCAL: " << std::endl;
        for (int i = controlPoints_.size(); i < cpnumber_; i++)
        {
            std::cout << mass(i) << std::endl;
        }
    }
    // else if (rank == 1)
    // {
    //     exportToParaviewFEM(0);
    //     testeParaviewFE();
    // }

    // for (Cell *cell : cells_part)
    // {
    //     cell->computeDistanceFromFEBoundary(quadrature_, boundaryFE_);
    // }
    // for (Cell *cell : cells_part)
    // {
    //     cell->computeDistanceFromFEBoundary(quadrature_, boundaryFE_);
    // }
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
        exportToParaviewISO(0);
    }
    else if (rank == 1)
    {
        exportToParaviewFEM(0);
    }
    //std::ofstream file1(text1.str());

    double initialNorm = 0.0;
    for (ControlPoint *node : controlPoints_)
    {
        double x1 = node->getInitialCoordinate()(0);
        double x2 = node->getInitialCoordinate()(1);
        initialNorm += x1 * x1 + x2 * x2;
    }
    for (Node *node : nodes_)
    {
        double x1 = node->getInitialCoordinate()(0);
        double x2 = node->getInitialCoordinate()(1);
        initialNorm += x1 * x1 + x2 * x2;
    }

    //////REVER QUANDO RECONSTRUIR MALHA!
    computeDistanceFromFEBoundary();
    incidenceLocalxGlobal();
    checkInactivesCPs();
    // checkControlPointsOutsideTheDomain();
    if (rank == 1)
    {
        testeParaviewFE();
    }
    /////

    int ndir = dirichletConditions_.size() + dirichletConditionsFE_.size() + cpOutsideDomain_.size();
    PetscMalloc1(ndir, &dof);
    int idir = 0;
    for (DirichletCondition *dir : dirichletConditions_)
    {
        int indexNode = dir->getControlPoint()->getIndex();
        int direction = dir->getDirection();
        dof[idir++] = (2 * indexNode + direction);
    }
    for (DirichletConditionFE *dir : dirichletConditionsFE_)
    {
        int indexNode = dir->getNode()->getIndex();
        int direction = dir->getDirection();
        dof[idir++] = (2 * indexNode + direction);
    }
    for (DirichletCondition *dir : cpOutsideDomain_)
    {
        int indexNode = dir->getControlPoint()->getIndex();
        int direction = dir->getDirection();
        // std::cout << "CP: " << indexNode << " DIRECTION: " << direction << std::endl;
        dof[idir++] = (2 * indexNode + direction);
    }

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
                                2 * cpnumber_, 2 * cpnumber_,
                                100, NULL, 300, NULL, &A);
            CHKERRQ(ierr);

            ierr = MatGetOwnershipRange(A, &Istart, &Iend);
            CHKERRQ(ierr);

            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b);
            CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, 2 * cpnumber_);
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
                for (NeumannConditionFE *con : neumannConditionsFE_)
                {
                    int ind = con->getNode()->getIndex();
                    int dir = con->getDirection();
                    double val1 = con->getValue() * (1.0 * loadStep / (1.0 * numberOfSteps_));
                    int dof = 2 * ind + dir;
                    ierr = VecSetValues(b, 1, &dof, &val1, ADD_VALUES);
                    // std::cout<<"NODE: "<<ind<<" DIRECTION: "<<dir<<" VALUE: "<<val1<<std::endl;
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
                for (DirichletConditionFE *con : dirichletConditionsFE_)
                {
                    Node *no = con->getNode();
                    int dir = con->getDirection();
                    double val1 = (con->getValue()) / (1.0 * numberOfSteps_);

                    no->incrementCurrentCoordinate(dir, val1);
                }
            }
            // computeDistanceFromFEBoundary();
            // incidenceLocalxGlobal();

            for (Cell *el : cells_part)
            {
                std::pair<vector<double>, matrix<double>> elementMatrices;
                elementMatrices = el->cellContributions(planeState_, "STATIC", loadStep, numberOfSteps_, 1.0, 0.25, 0.5, quadrature_);
                int num = el->getControlPoints().size();

                for (size_t i = 0; i < num; i++)
                {

                    int dof = 2 * el->getControlPoint(i)->getIndex();
                    ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);

                    dof = 2 * el->getControlPoint(i)->getIndex() + 1;
                    ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);

                    for (size_t j = 0; j < num; j++)
                    {

                        int dof1 = 2 * el->getControlPoint(i)->getIndex();
                        int dof2 = 2 * el->getControlPoint(j)->getIndex();
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);

                        dof1 = 2 * el->getControlPoint(i)->getIndex() + 1;
                        dof2 = 2 * el->getControlPoint(j)->getIndex();
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);

                        dof1 = 2 * el->getControlPoint(i)->getIndex();
                        dof2 = 2 * el->getControlPoint(j)->getIndex() + 1;
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);

                        dof1 = 2 * el->getControlPoint(i)->getIndex() + 1;
                        dof2 = 2 * el->getControlPoint(j)->getIndex() + 1;
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
                    }
                }
            }

            for (Element *el : elements_part)
            {
                std::pair<vector<double>, matrix<double>> elementMatrices;
                elementMatrices = el->elementContributions(planeState_, "STATIC", loadStep, numberOfSteps_, 1.0, 0.25, 0.5);
                std::vector<int> freedom = el->getFreedomDegree();
                int num = freedom.size();

                for (size_t i = 0; i < num; i++)
                {

                    int dof = 2 * freedom[i];
                    ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);

                    dof = 2 * freedom[i] + 1;
                    ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);

                    for (size_t j = 0; j < num; j++)
                    {

                        int dof1 = 2 * freedom[i];
                        int dof2 = 2 * freedom[j];
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);

                        dof1 = 2 * freedom[i] + 1;
                        dof2 = 2 * freedom[j];
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);

                        dof1 = 2 * freedom[i];
                        dof2 = 2 * freedom[j] + 1;
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);

                        dof1 = 2 * freedom[i] + 1;
                        dof2 = 2 * freedom[j] + 1;
                        ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j + 1), ADD_VALUES);
                    }
                }
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

            // MatView(A, PETSC_VIEWER_STDOUT_WORLD);
            // CHKERRQ(ierr);

            // VecView(b, PETSC_VIEWER_STDOUT_WORLD);
            // CHKERRQ(ierr);

            //Zerando linhas e colunas
            MatZeroRowsColumns(A, ndir, dof, 1.0, x, b);

            //MatView(A, PETSC_VIEWER_STDOUT_WORLD);
            // CHKERRQ(ierr);

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

            for (ControlPoint *cp : controlPoints_)
            {
                int newIndex = cp->getIndex();
                Idof = 2 * newIndex;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                cp->incrementCurrentCoordinate(0, val);

                Idof = 2 * newIndex + 1;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                cp->incrementCurrentCoordinate(1, val);
            }
            for (Node *node : nodes_)
            {
                int newIndex = node->getIndex();
                Idof = 2 * newIndex;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                node->incrementCurrentCoordinate(0, val);

                Idof = 2 * newIndex + 1;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                node->incrementCurrentCoordinate(1, val);
            }
            for (CP_BlendingZone *ble : cpsInsideBledingZone_)
            {
                ble->interpolateGlobalCoordinate();
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
            exportToParaviewISO(loadStep);
        }
        else if (rank == 1)
        {
            for (Node *n : nodes_)
            {
                n->setZeroStressState();
            }
            for (Element *el : elements_)
            {
                el->StressCalculate(planeState_);
            }
            exportToParaviewFEM(loadStep);
        }
    }
    PetscFree(dof);
    return 0;
}

void GlobalSolid::exportToParaviewISO(const int &loadstep)
{
    matrix<double> qxsi2 = coordinatesForInterpolation(orderParaview_);

    int auxxx = (orderParaview_ + 1) * (orderParaview_ + 1);
    int nIsoPoints = auxxx * cells_.size();
    std::stringstream text;
    text << "outputISO" << loadstep << ".vtu";
    std::ofstream file(text.str());

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << nIsoPoints
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
        bounded_vector<double, 2> coord;
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
            phi2_ = cell->shapeFunction(qxsi, wpc2, INC_);
            coord(0) = 0.0;
            coord(1) = 0.0;
            for (int i = 0; i < connection.size(); i++)
            {
                x = connection[i]->getCurrentCoordinate();
                coord(0) += phi2_(i) * x(0) / wpc2(i);
                coord(1) += phi2_(i) * x(1) / wpc2(i);
            }
            file << coord(0) << " " << coord(1) << " " << 0.0 << std::endl;
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
    int aux = auxxx;
    for (Cell *cell : cells_)
    {
        file << aux << std::endl;
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
        bounded_vector<double, 2> x, qxsi;
        bounded_vector<double, 2> disp;
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
            phi2_ = cell->shapeFunction(qxsi, wpc2, INC_);
            disp(0) = 0.0;
            disp(1) = 0.0;
            for (int i = 0; i < connection.size(); i++)
            {
                x = connection[i]->getCurrentDisplacement();
                disp(0) += phi2_(i) * x(0);
                disp(1) += phi2_(i) * x(1);
            }
            file << disp(0) << " " << disp(1) << std::endl;
        }
    }
    file << "      </DataArray> "
         << "\n";

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

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"4\" "
    //      << "Name=\"GreenStrain\" format=\"ascii\">"
    //      << "\n";
    // for (Cell *cell : cells_)
    // {
    //     bounded_vector<double, 2> qxsi;
    //     bounded_vector<double, 4> green;
    //     for (int j = 0; j < auxxx; j++)
    //     {
    //         qxsi(0) = qxsi2(j, 0);
    //         qxsi(1) = qxsi2(j, 1);
    //         green = cell->getGreen(qxsi, planeState_);
    //         file << green(0) << " " << green(1) << " " << green(2) << " " << green(3) << std::endl;
    //     }
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

void GlobalSolid::exportToParaviewFEM(const int &num)
{
    std::stringstream name;
    name << "outputFEM" << num << ".vtu";
    std::ofstream file(name.str());

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << nodes_.size()
         << "\"  NumberOfCells=\"" << elements_.size()
         << "\">"
         << "\n";

    //nodal coordinates
    file << "    <Points>"
         << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">"
         << "\n";
    for (Node *no : nodes_)
    {
        bounded_vector<double, 2> coord = no->getCurrentCoordinate();
        file << coord(0) << " " << coord(1) << " " << 0.0 << "\n";
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
    for (Element *el : elements_)
    {
        std::vector<Node *> conec = el->getConnection();
        for (Node *no : conec)
        {
            file << no->getIndexFE() << " ";
        }
        file << "\n";
    }
    file << "      </DataArray>"
         << "\n";

    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";
    int aux = 0;
    for (Element *el : elements_)
    {
        int n = el->getConnection().size();
        aux += n;
        file << aux << "\n";
    }
    file << "      </DataArray>"
         << "\n";

    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";
    for (Element *e : elements_)
    {
        file << 69 << "\n";
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
    for (Node *n : nodes_)
    {
        bounded_vector<double, 2> disp = n->getCurrentDisplacement();

        file << disp(0) << " " << disp(1) << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"4\" "
         << "Name=\"CauchyStress\" format=\"ascii\">"
         << "\n";
    for (Node *n : nodes_)
    {
        bounded_vector<double, 5> stress = n->getStressState(); //(x1, x2, x12, contador)
        file << stress(0) / stress(4) << " " << stress(1) / stress(4) << " " << stress(2) / stress(4) << " " << stress(3) / stress(4) << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"4\" "
    //      << "Name=\"GreenStrain\" format=\"ascii\">"
    //      << "\n";

    // for (Node *n : nodes_)
    // {
    //     file << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
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

    for (Element *el : elements_)
    {
        file << elementPartition_[el->getIndex()] << "\n";
    }
    file << "      </DataArray> "
         << "\n";

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
    mirror2 = "ISOdomain_decomposition.txt";
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

void GlobalSolid::domainDecompositionMETIS(const std::string &elementType)
{
    std::string mirror2;
    mirror2 = "FEdomain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());

    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = elements_.size();
    idx_t numNd = nodes_.size();
    idx_t ssize = size;
    idx_t one = 1;
    idx_t n;
    if (elementType == "T3")
        n = 3;
    else if (elementType == "T6")
        n = 6;
    else if (elementType == "T10")
        n = 10;
    idx_t elem_start[numEl + 1], elem_connec[n * numEl];
    elementPartition_ = new idx_t[numEl];
    nodePartition_ = new idx_t[numNd];
    for (idx_t i = 0; i < numEl + 1; i++)
    {
        elem_start[i] = n * i;
    }
    for (idx_t jel = 0; jel < numEl; jel++)
    {
        for (idx_t i = 0; i < n; i++)
        {
            int nodeIndex = elements_[jel]->getNode(i)->getIndexFE();
            elem_connec[n * jel + i] = nodeIndex;
        }
    }
    //Performs the domain decomposition
    METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec,
                       NULL, NULL, &one, &ssize, NULL, NULL,
                       &objval, elementPartition_, nodePartition_);

    mirrorData << std::endl
               << "DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
    for (int i = 0; i < elements_.size(); i++)
    {
        mirrorData << "process = " << elementPartition_[i]
                   << ", element = " << i << std::endl;
    }

    mirrorData << std::endl
               << "DOMAIN DECOMPOSITION - NODES" << std::endl;
    for (int i = 0; i < nodes_.size(); i++)
    {
        mirrorData << "process = " << nodePartition_[i]
                   << ", node = " << i << std::endl;
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
        exportToParaviewISO(0);
        //text1 << "DeslocamentoxTempo.txt";
    }
    else if (rank == 1)
    {
        exportToParaviewFEM(0);
    }
    // std::ofstream file1(text1.str());

    double initialNorm = 0.0;
    for (ControlPoint *node : controlPoints_)
    {
        double x1 = node->getInitialCoordinate()(0);
        double x2 = node->getInitialCoordinate()(1);
        initialNorm += x1 * x1 + x2 * x2;
    }
    for (Node *node : nodes_)
    {
        double x1 = node->getInitialCoordinate()(0);
        double x2 = node->getInitialCoordinate()(1);
        initialNorm += x1 * x1 + x2 * x2;
    }

    int ndir = dirichletConditions_.size() + dirichletConditionsFE_.size();
    PetscMalloc1(ndir, &dof);
    int idir = 0;
    for (DirichletCondition *dir : dirichletConditions_)
    {
        int indexNode = dir->getControlPoint()->getIndex();
        int direction = dir->getDirection();
        dof[idir++] = (2 * indexNode + direction);
    }
    for (DirichletConditionFE *dir : dirichletConditionsFE_)
    {
        int indexNode = dir->getNode()->getIndex();
        int direction = dir->getDirection();
        dof[idir++] = (2 * indexNode + direction);
    }

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
                                2 * cpnumber_, 2 * cpnumber_,
                                100, NULL, 300, NULL, &A);
            CHKERRQ(ierr);

            ierr = MatGetOwnershipRange(A, &Istart, &Iend);
            CHKERRQ(ierr);

            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD, &b);
            CHKERRQ(ierr);
            ierr = VecSetSizes(b, PETSC_DECIDE, 2 * cpnumber_);
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
                for (NeumannConditionFE *con : neumannConditionsFE_)
                {
                    int ind = con->getNode()->getIndex();
                    int dir = con->getDirection();
                    double val1 = con->getValue();
                    int dof = 2 * ind + dir;
                    ierr = VecSetValues(b, 1, &dof, &val1, ADD_VALUES);
                    // std::cout<<"NODE: "<<ind<<" DIRECTION: "<<dir<<" VALUE: "<<val1<<std::endl;
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
                for (DirichletConditionFE *con : dirichletConditionsFE_)
                {
                    Node *no = con->getNode();
                    int dir = con->getDirection();
                    double val1 = (con->getValue()) / (1.0 * numberOfSteps_);

                    no->incrementCurrentCoordinate(dir, val1);
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

            for (Element *el : elements_part)
            {
                std::pair<vector<double>, matrix<double>> elementMatrices;
                elementMatrices = el->elementContributions(planeState_, "DYNAMIC", 1, 1, deltat_, beta_, gamma_);
                int num = el->getConnection().size();

                for (size_t i = 0; i < num; i++)
                {
                    if (fabs(elementMatrices.first(2 * i)) >= 1.0e-11)
                    {
                        int dof = 2 * el->getNode(i)->getIndex();
                        ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
                    }
                    if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-11)
                    {
                        int dof = 2 * el->getNode(i)->getIndex() + 1;
                        ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
                    }

                    for (size_t j = 0; j < num; j++)
                    {
                        if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-11)
                        {
                            int dof1 = 2 * el->getNode(i)->getIndex();
                            int dof2 = 2 * el->getNode(j)->getIndex();
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-11)
                        {
                            int dof1 = 2 * el->getNode(i)->getIndex() + 1;
                            int dof2 = 2 * el->getNode(j)->getIndex();
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i + 1, 2 * j), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-11)
                        {
                            int dof1 = 2 * el->getNode(i)->getIndex();
                            int dof2 = 2 * el->getNode(j)->getIndex() + 1;
                            ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &elementMatrices.second(2 * i, 2 * j + 1), ADD_VALUES);
                        }
                        if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-11)
                        {
                            int dof1 = 2 * el->getNode(i)->getIndex() + 1;
                            int dof2 = 2 * el->getNode(j)->getIndex() + 1;
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

            MatZeroRowsColumns(A, ndir, dof, 1.0, x, b);

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
            for (Node *node : nodes_)
            {
                int newIndex = node->getIndex();
                Idof = 2 * newIndex;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                node->incrementCurrentCoordinate(0, val);

                Idof = 2 * newIndex + 1;
                ierr = VecGetValues(All, Ione, &Idof, &val);
                CHKERRQ(ierr);
                norm += val * val;
                node->incrementCurrentCoordinate(1, val);

                bounded_vector<double, 2> vel, accel;
                accel = node->getCurrentCoordinate() / (beta_ * deltat_ * deltat_) - node->getPastCoordinate() / (beta_ * deltat_ * deltat_) -
                        node->getPastVelocity() / (beta_ * deltat_) - node->getPastAcceleration() * (0.5 / beta_ - 1.0);
                node->setCurrentAcceleration(accel);

                vel = gamma_ * deltat_ * node->getCurrentAcceleration() + node->getPastVelocity() + deltat_ * (1.0 - gamma_) * node->getPastAcceleration();
                node->setCurrentVelocity(vel);
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
            exportToParaviewISO(timeStep);
            // file1 << timeStep * deltat_ << " " << controlPoints_[1999]->getCurrentCoordinate()[1] - controlPoints_[1999]->getInitialCoordinate()[1] << " " << controlPoints_[1999]->getCurrentCoordinate()[0] - controlPoints_[1999]->getInitialCoordinate()[0] << std::endl;
        }
        else if (rank == 1)
        {
            for (Node *n : nodes_)
            {
                n->setZeroStressState();
            }
            for (Element *el : elements_)
            {
                el->StressCalculate(planeState_);
            }
            exportToParaviewFEM(timeStep);
        }

        // for (ControlPoint *cp : controlPoints_)
        // {
        //     bounded_vector<double, 2> coordinate = cp->getCurrentCoordinate();
        //     cp->setPastCoordinate(coordinate);
        //     bounded_vector<double, 2> vel = cp->getCurrentVelocity();
        //     cp->setPastVelocity(vel);
        //     bounded_vector<double, 2> accel = cp->getCurrentAcceleration();
        //     cp->setPastAcceleration(accel);
        // }

        for (ControlPoint *cp : controlPoints_)
        {
            cp->updatePastValue();
            bounded_vector<double, 2> vel, accel;
            accel = cp->getCurrentCoordinate() / (beta_ * deltat_ * deltat_) - cp->getPastCoordinate() / (beta_ * deltat_ * deltat_) -
                    cp->getPastVelocity() / (beta_ * deltat_) - cp->getPastAcceleration() * (0.5 / beta_ - 1.0);
            cp->setCurrentAcceleration(accel);
            vel = gamma_ * deltat_ * cp->getCurrentAcceleration() + cp->getPastVelocity() + deltat_ * (1.0 - gamma_) * cp->getPastAcceleration();
            cp->setCurrentVelocity(vel);
        }
        for (Node *node : nodes_)
        {
            node->updatePastValue();
            bounded_vector<double, 2> vel, accel;
            accel = node->getCurrentCoordinate() / (beta_ * deltat_ * deltat_) - node->getPastCoordinate() / (beta_ * deltat_ * deltat_) -
                    node->getPastVelocity() / (beta_ * deltat_) - node->getPastAcceleration() * (0.5 / beta_ - 1.0);
            node->setCurrentAcceleration(accel);
            vel = gamma_ * deltat_ * node->getCurrentAcceleration() + node->getPastVelocity() + deltat_ * (1.0 - gamma_) * node->getPastAcceleration();
            node->setCurrentVelocity(vel);
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

    int ndir = dirichletConditions_.size() + dirichletConditionsFE_.size();
    PetscMalloc1(ndir, &dof);
    int idir = 0;
    for (DirichletCondition *dir : dirichletConditions_)
    {
        int indexNode = dir->getControlPoint()->getIndex();
        int direction = dir->getDirection();
        dof[idir++] = (2 * indexNode + direction);
    }
    for (DirichletConditionFE *dir : dirichletConditionsFE_)
    {
        int indexNode = dir->getNode()->getIndex();
        int direction = dir->getDirection();
        dof[idir++] = (2 * indexNode + direction);
    }

    //Create PETSc sparse parallel matrix
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        2 * cpnumber_, 2 * cpnumber_,
                        100, NULL, 300, NULL, &A);
    CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, 2 * cpnumber_);
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
        for (NeumannConditionFE *con : neumannConditionsFE_)
        {
            int ind = con->getNode()->getIndex();
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
    for (DirichletConditionFE *con : dirichletConditionsFE_)
    {
        Node *no = con->getNode();
        int dir = con->getDirection();
        double val1 = (con->getValue()) / (1.0 * numberOfSteps_);

        no->incrementCurrentCoordinate(dir, val1);
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
    for (Element *el : elements_part)
    {
        std::pair<vector<double>, matrix<double>> elementMatrices;
        elementMatrices = el->elementContributions(planeState_, "STATIC", 1, 1, deltat_, beta_, gamma_);
        int n = elementMatrices.first.size();
        matrix<double> massLocal;
        massLocal = el->massMatrix();
        int num = el->getConnection().size();

        for (size_t i = 0; i < num; i++)
        {
            if (fabs(elementMatrices.first(2 * i)) >= 1.0e-15)
            {
                int dof = 2 * el->getNode(i)->getIndex();
                ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i), ADD_VALUES);
            }
            if (fabs(elementMatrices.first(2 * i + 1)) >= 1.0e-15)
            {
                int dof = 2 * el->getNode(i)->getIndex() + 1;
                ierr = VecSetValues(b, 1, &dof, &elementMatrices.first(2 * i + 1), ADD_VALUES);
            }

            for (size_t j = 0; j < num; j++)
            {
                if (fabs(elementMatrices.second(2 * i, 2 * j)) >= 1.e-15)
                {
                    int dof1 = 2 * el->getNode(i)->getIndex();
                    int dof2 = 2 * el->getNode(j)->getIndex();
                    ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j), ADD_VALUES);
                }
                if (fabs(elementMatrices.second(2 * i + 1, 2 * j)) >= 1.e-15)
                {
                    int dof1 = 2 * el->getNode(i)->getIndex() + 1;
                    int dof2 = 2 * el->getNode(j)->getIndex();
                    ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i + 1, 2 * j), ADD_VALUES);
                }
                if (fabs(elementMatrices.second(2 * i, 2 * j + 1)) >= 1.e-15)
                {
                    int dof1 = 2 * el->getNode(i)->getIndex();
                    int dof2 = 2 * el->getNode(j)->getIndex() + 1;
                    ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j + 1), ADD_VALUES);
                }
                if (fabs(elementMatrices.second(2 * i + 1, 2 * j + 1)) >= 1.e-15)
                {
                    int dof1 = 2 * el->getNode(i)->getIndex() + 1;
                    int dof2 = 2 * el->getNode(j)->getIndex() + 1;
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

    // VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);
    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);
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
    for (Node *node : nodes_)
    {
        bounded_vector<double, 2> firstAccel;
        int i = node->getIndex();
        Idof = 2 * i;
        ierr = VecGetValues(All, Ione, &Idof, &val);
        CHKERRQ(ierr);
        firstAccel(0) = val;

        Idof = 2 * i + 1;
        ierr = VecGetValues(All, Ione, &Idof, &val);
        CHKERRQ(ierr);
        firstAccel(1) = val;
        node->setCurrentAcceleration(firstAccel);
        node->setPastAcceleration(firstAccel);
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

void GlobalSolid::incidenceLocalxGlobal()
{
    // int rank;
    // MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /*Looping para:
        -Encontrar pontos de hammer dentro da blending zone;
            -Célula em que o ponto de hammer está inserido;
            -Xsis correspondentes ao ponto de hammer na célula;
            -Computar valor de b e db_dxsiLocal
    */
    for (Element *el : elements_) ///Looping para encontrar pontos de hammer dentro da blending zone;
    {
        std::vector<Node *> conec = el->getConnection();
        matrix<double> integrationPoints = el->hammerQuadrature();

        //Variáveis para serem exportadas para os elementos
        std::vector<bool> insideBlendZone;                       //verdadeiro se o ponto de hammer está na blend zone e encontrou uma célula;
        std::vector<bounded_vector<double, 2>> xsiIncidenceCell; //coordenadas adimensionais do global que corresponde ao ponto de hammer
        std::vector<Cell *> incidenceCell;                       //célula que houve incidência do ponto de hammer
        std::vector<double> bValue;                              //value of assigned distance
        std::vector<bounded_vector<double, 2>> db_dxsiValue;     //derivada de b em relação a xsi

        for (int ih = 0; ih < integrationPoints.size1(); ih++)
        {
            double DAhammer = 0.0;             //distância interpolada do ponto de hammer ao contorno de elementos finitos
            bounded_vector<double, 2> coordFE; //coordenada global do ponto de hammer
            coordFE(0) = 0.0;
            coordFE(1) = 0.0;
            vector<double> phi = el->domainShapeFunction(integrationPoints(ih, 0), integrationPoints(ih, 1));

            for (int i = 0; i < conec.size(); i++)
            {
                DAhammer += phi(i) * conec[i]->getDistanceToBoundary();
                coordFE += phi(i) * conec[i]->getCurrentCoordinate();
            }

            if (DAhammer - blendZoneThickness_ <= 1.0e-10) //inside of blend zone
            {
                insideBlendZone.push_back(false);

                for (Cell *cell : cells_)
                {
                    bounded_vector<double, 2> xsiISO, coordISO, deltaxsi;
                    deltaxsi(0) = 100000.0;
                    deltaxsi(1) = 100000.0;
                    xsiISO(0) = 0.0; //first attempt xs1
                    xsiISO(1) = 0.0; //first attempt xs2

                    std::vector<ControlPoint *> conPoints = cell->getControlPoints();
                    vector<double> wpc2(conPoints.size());
                    for (int i = 0; i < conPoints.size(); i++)
                    {
                        wpc2(i) = conPoints[i]->getWeight();
                    }

                    int cont = 0;
                    while (deltaxsi(0) >= 1.0e-06 and deltaxsi(1) >= 1.0e-06 and //newton-raphson para encontrar os xsis
                           xsiISO(0) >= -1.0 and xsiISO(0) <= 1.0 and
                           xsiISO(1) >= -1.0 and xsiISO(1) <= 1.0 and cont <= 15)
                    {
                        std::pair<vector<double>, matrix<double>> functions;
                        functions = cell->shapeFunctionAndDerivates(xsiISO);

                        coordISO(0) = 0.0;
                        coordISO(1) = 0.0;

                        for (int cp = 0; cp < conPoints.size(); cp++)
                        {
                            bounded_vector<double, 2> coordinateCP = conPoints[cp]->getCurrentCoordinate();
                            coordISO(0) += functions.first(cp) * coordinateCP(0) / wpc2(cp);
                            coordISO(1) += functions.first(cp) * coordinateCP(1) / wpc2(cp);
                        }

                        bounded_matrix<double, 2, 2> jacobian = cell->referenceJacobianMatrix(functions.second);
                        bounded_matrix<double, 2, 2> inverse = inverseMatrix(jacobian);
                        deltaxsi = prod(inverse, coordFE - coordISO);

                        xsiISO = xsiISO + deltaxsi;

                        cont++;
                    }

                    if (xsiISO(0) >= -1.0 and xsiISO(0) <= 1.0 and //ponto de hammer encontrou uma célula
                        xsiISO(1) >= -1.0 and xsiISO(1) <= 1.0)
                    {
                        insideBlendZone[ih] = true;
                        xsiIncidenceCell.push_back(xsiISO);
                        incidenceCell.push_back(cell);

                        bounded_vector<double, 2> blend = blendFunction(DAhammer); //(b, db_dDAhammer)
                        bValue.push_back(blend(0));

                        bounded_vector<double, 2> dDA_dxsi;
                        matrix<double> dphi_dxsi = el->domainDerivativeShapeFunction(integrationPoints(ih, 0), integrationPoints(ih, 1));
                        dDA_dxsi(0) = 0.0;
                        dDA_dxsi(1) = 0.0;
                        for (int i = 0; i < conec.size(); i++)
                        {
                            double DAnode = conec[i]->getDistanceToBoundary();
                            dDA_dxsi(0) = dphi_dxsi(0, i) * DAnode; //dDAhammer_dxsiLocal1
                            dDA_dxsi(1) = dphi_dxsi(1, i) * DAnode; //dDAhammer_dxsiLocal1
                        }
                        bounded_vector<double, 2> db_dxsi; 
                        db_dxsi(0) = blend(1) * dDA_dxsi(0); //(db_dxsiLocal1)
                        db_dxsi(1) = blend(1) * dDA_dxsi(1); //(db_dxsiLocal2) 
                        db_dxsiValue.push_back(db_dxsi);

                        break;
                    }
                }
            }
            else
            {
                insideBlendZone.push_back(false);
            }
        }
        el->setIncidenceOnGlobal(insideBlendZone, incidenceCell, xsiIncidenceCell, bValue, db_dxsiValue);
    }

    // /*Looping para:
    //     -Encontrar pontos de controle dentro da malha local;
    //         -Elemento em que o ponto de controle está inserido;
    //         -Xsis correspondentes ao ponto de controle no elemento;
    //     -Encontrar pontos de controle fora do domínio global;
    //         -Adicionar condição de contorno;
    // */
    // if (cpOutsideDomain_.size() > 1)
    // {
    //     cpOutsideDomain_.erase(cpOutsideDomain_.begin(), cpOutsideDomain_.begin() + cpOutsideDomain_.size());
    // }
    // for (ControlPoint *cp : controlPoints_)
    // {
    //     bounded_vector<double, 2> coord = cp->getCurrentCoordinate() / cp->getWeight();
    //     double distanceCP = 1000000.0;
    //     double distance;

    //     for (BoundaryElement *bound : boundaryFE_)
    //     {
    //         double xsiBoundary = 0.0; //primeira tentativa
    //         double deltaxsi = 100000.0;
    //         int cont = 0;
    //         std::vector<Node *> boundaryConec = bound->getNodes();
    //         bounded_vector<double, 2> coordBoundary, firstDerivate, secondDerivate, normal;

    //         while (fabs(deltaxsi) >= 1.0e-06 and xsiBoundary >= -1.0 and xsiBoundary <= 1.0 and cont <= 15)
    //         {
    //             matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(xsiBoundary); //PHI, PHI', PHI''
    //             coordBoundary(0) = 0.0;
    //             coordBoundary(1) = 0.0;
    //             firstDerivate(0) = 0.0;
    //             firstDerivate(1) = 0.0;
    //             secondDerivate(0) = 0.0;
    //             secondDerivate(1) = 0.0;
    //             int aux = 0;
    //             for (Node *node : boundaryConec)
    //             {
    //                 bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();

    //                 coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
    //                 coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

    //                 firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
    //                 firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

    //                 secondDerivate(0) += boundaryFunctions(aux, 2) * coordinateNode(0);
    //                 secondDerivate(1) += boundaryFunctions(aux, 2) * coordinateNode(1);
    //                 aux = aux + 1;
    //             }
    //             double h = -((coord(0) - coordBoundary(0)) * (-firstDerivate(0)) + (coord(1) - coordBoundary(1)) * (-firstDerivate(1)));
    //             double dh_dxsi = (-firstDerivate(0)) * (-firstDerivate(0)) + (-secondDerivate(0)) * (coord(0) - coordBoundary(0)) +
    //                              (-firstDerivate(1)) * (-firstDerivate(1)) + (-secondDerivate(1)) * (coord(1) - coordBoundary(1));

    //             deltaxsi = h / dh_dxsi;

    //             xsiBoundary = xsiBoundary + deltaxsi;

    //             cont++;
    //         }

    //         if (xsiBoundary <= -1.0)
    //         {
    //             matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(-1.0); //PHI, PHI', PHI''
    //             int aux = 0;
    //             coordBoundary(0) = 0.0;
    //             coordBoundary(1) = 0.0;
    //             firstDerivate(0) = 0.0;
    //             firstDerivate(1) = 0.0;
    //             for (Node *node : boundaryConec)
    //             {
    //                 bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
    //                 coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
    //                 coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

    //                 firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
    //                 firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

    //                 aux = aux + 1;
    //             }
    //         }
    //         else if (xsiBoundary >= 1.0)
    //         {
    //             matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(1.0); //PHI, PHI', PHI''
    //             int aux = 0;
    //             coordBoundary(0) = 0.0;
    //             coordBoundary(1) = 0.0;
    //             firstDerivate(0) = 0.0;
    //             firstDerivate(1) = 0.0;
    //             for (Node *node : boundaryConec)
    //             {
    //                 bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
    //                 coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
    //                 coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

    //                 firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
    //                 firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

    //                 aux = aux + 1;
    //             }
    //         }

    //         distance = sqrt((coord(0) - coordBoundary(0)) * (coord(0) - coordBoundary(0)) + (coord(1) - coordBoundary(1)) * (coord(1) - coordBoundary(1)));
    //         normal(0) = firstDerivate(1);
    //         normal(1) = -firstDerivate(0);
    //         double aux = (coord(0) - coordBoundary(0)) * normal(0) + (coord(1) - coordBoundary(1)) * normal(1);

    //         if (distance < 1.0e-10)
    //         {
    //             distance = 0.0;
    //         }

    //         if (distance < fabs(distanceCP))
    //         {
    //             if (aux > 0.0) //fora do domínio local
    //             {
    //                 distanceCP = -distance;
    //             }
    //             else if (aux < 0.0) //dentro do domínio local
    //             {
    //                 distanceCP = distance;
    //             }
    //             else if (aux == 0.0 and distance == 0.0) //caso único: exatamente sobre um elemento do contorno
    //             {
    //                 distanceCP = distance;
    //             }
    //         }
    //     }
    //     if (distanceCP - blendZoneThickness_ > 1.0e-10) //dentro do local, mas fora da blend zone
    //     {
    //         std::cout << cp->getIndex() << " " << distanceCP << " " << distanceCP - blendZoneThickness_ << std::endl;

    //         DirichletCondition *dir0 = new DirichletCondition(cp, 0, 0.0);
    //         cpOutsideDomain_.push_back(dir0);

    //         DirichletCondition *dir1 = new DirichletCondition(cp, 1, 0.0);
    //         cpOutsideDomain_.push_back(dir1);

    //         for (Element *el : elements_)
    //         {
    //             bounded_vector<double, 2> xsiFE, coordFE, deltaxsi;
    //             deltaxsi(0) = 100000.0;
    //             deltaxsi(1) = 100000.0;
    //             xsiFE(0) = 0.5; //first attempt xs1
    //             xsiFE(1) = 0.5; //first attempt xs2

    //             std::vector<Node *> conec = el->getConnection();

    //             int cont = 0;
    //             while (deltaxsi(0) >= 1.0e-06 and deltaxsi(1) >= 1.0e-06 and
    //                    xsiFE(0) >= 0.0 and xsiFE(0) <= 1.0 and
    //                    xsiFE(1) >= 0.0 and xsiFE(1) <= 1.0 and cont <= 15)
    //             {

    //                 coordFE = el->calculateGlobalCoordinate(xsiFE);

    //                 bounded_matrix<double, 2, 2> jacobian = el->referenceJacobianMatrix(xsiFE(0), xsiFE(1));
    //                 bounded_matrix<double, 2, 2> inverse = inverseMatrix(jacobian);
    //                 deltaxsi = prod(inverse, coord - coordFE);

    //                 xsiFE = xsiFE + deltaxsi;

    //                 cont++;
    //             }

    //             if (xsiFE(0) >= 0.0 and xsiFE(0) <= 1.0 and //ponto de hammer encontrou uma célula
    //                 xsiFE(1) >= 0.0 and xsiFE(1) <= 1.0 and
    //                 xsiFE(0) + xsiFE(1) <= 1.0) //só vai funcionar para triângulos não muito deformados...
    //             {
    //                 CP_BlendingZone *inside = new CP_BlendingZone(cp, el, xsiFE);
    //                 cpsInsideBledingZone_.push_back(inside);

    //                 break;
    //             }
    //         }
    //     }
    // }
    // if (distanceCP >= 0.0) //teste
    // {
    //     for (Element *el : elements_)
    //     {
    //         bounded_vector<double, 2> xsiFE, coordFE, deltaxsi;
    //         deltaxsi(0) = 100000.0;
    //         deltaxsi(1) = 100000.0;
    //         xsiFE(0) = 0.5; //first attempt xs1
    //         xsiFE(1) = 0.5; //first attempt xs2

    //         std::vector<Node *> conec = el->getConnection();

    //         int cont = 0;
    //         while (deltaxsi(0) >= 1.0e-06 and deltaxsi(1) >= 1.0e-06 and
    //                xsiFE(0) >= 0.0 and xsiFE(0) <= 1.0 and
    //                xsiFE(1) >= 0.0 and xsiFE(1) <= 1.0 and cont <= 15)
    //         {

    //             coordFE = el->calculateGlobalCoordinate(xsiFE);

    //             bounded_matrix<double, 2, 2> jacobian = el->referenceJacobianMatrix(xsiFE(0), xsiFE(1));
    //             bounded_matrix<double, 2, 2> inverse = inverseMatrix(jacobian);
    //             deltaxsi = prod(inverse, coord - coordFE);

    //             xsiFE = xsiFE + deltaxsi;

    //             cont++;
    //         }

    //         if (xsiFE(0) >= 0.0 and xsiFE(0) <= 1.0 and //ponto de hammer encontrou uma célula
    //             xsiFE(1) >= 0.0 and xsiFE(1) <= 1.0 and
    //             xsiFE(0) + xsiFE(1) <= 1.0) //só vai funcionar para triângulos não muito deformados...
    //         {
    //             CP_BlendingZone *inside = new CP_BlendingZone(cp, el, xsiFE);
    //             cpsInsideBledingZone_.push_back(inside);

    //             break;
    //         }
    //     }
    // }

    // else
    // {
    //     cp->setInsideLocal(false);
    // }
}
// for (Element *el : elements_) ///COLOCAR ELEMENTS_PART
// {
//     // std::cout << "rank: " << rank << " index: " << el->getIndex() << std::endl;
//     matrix<double> integrationPoints = el->hammerQuadrature();
//     std::vector<double> distanceFE = el->getDistanceFromFEBoundary();
//     std::vector<bool> insideBlendZone;
//     std::vector<bounded_vector<double, 2>> xsiIncidenceCell;
//     std::vector<Cell *> incidenceCell;

//     for (int ih = 0; ih < integrationPoints.size1(); ih++)
//     {
//         if (distanceFE[ih] - blendZoneThickness_ <= 0.0) //inside of blend zone
//         {
//             bounded_vector<double, 2> xsiFE, coordFE;
//             xsiFE(0) = integrationPoints(ih, 0);
//             xsiFE(1) = integrationPoints(ih, 1);

//             coordFE = el->calculateGlobalCoordinate(xsiFE); //coordenates of hammer points
//             insideBlendZone.push_back(false);

//             for (Cell *cell : cells_)
//             {
//                 bounded_vector<double, 2> xsiISO, coordISO, deltaxsi;
//                 deltaxsi(0) = 100000.0;
//                 deltaxsi(1) = 100000.0;
//                 xsiISO(0) = 0.0; //first attempt xs1
//                 xsiISO(1) = 0.0; //first attempt xs2

//                 std::vector<ControlPoint *> conPoints = cell->getControlPoints();
//                 vector<double> wpc2(conPoints.size());
//                 for (int i = 0; i < conPoints.size(); i++)
//                 {
//                     wpc2(i) = conPoints[i]->getWeight();
//                 }

//                 int cont = 0;
//                 while (deltaxsi(0) >= 1.0e-06 and deltaxsi(1) >= 1.0e-06 and
//                        xsiISO(0) >= -1.0 and xsiISO(0) <= 1.0 and
//                        xsiISO(1) >= -1.0 and xsiISO(1) <= 1.0 and cont <= 15)
//                 {
//                     std::pair<vector<double>, matrix<double>> functions;
//                     functions = cell->shapeFunctionAndDerivates(xsiISO);

//                     coordISO(0) = 0.0;
//                     coordISO(1) = 0.0;

//                     for (int cp = 0; cp < conPoints.size(); cp++)
//                     {
//                         bounded_vector<double, 2> coordinateCP = conPoints[cp]->getCurrentCoordinate();
//                         coordISO(0) += functions.first(cp) * coordinateCP(0) / wpc2(cp);
//                         coordISO(1) += functions.first(cp) * coordinateCP(1) / wpc2(cp);
//                     }

//                     bounded_matrix<double, 2, 2> jacobian = cell->referenceJacobianMatrix(functions.second);
//                     bounded_matrix<double, 2, 2> inverse = inverseMatrix(jacobian);
//                     deltaxsi = prod(inverse, coordFE - coordISO);

//                     xsiISO = xsiISO + deltaxsi;

//                     cont++;
//                 }

//                 if (xsiISO(0) >= -1.0 and xsiISO(0) <= 1.0 and
//                     xsiISO(1) >= -1.0 and xsiISO(1) <= 1.0)
//                 {
//                     insideBlendZone[ih] = true;
//                     xsiIncidenceCell.push_back(xsiISO);
//                     incidenceCell.push_back(cell);

//                     // std::cout << "INDEX: " << el->getIndex() << " HAMMER: " << ih << " ISO: " << cell->calculateGlobalCoordinate(xsiISO)(0) << " " << cell->calculateGlobalCoordinate(xsiISO)(1) << std::endl;
//                     // std::cout << "INDEX: " << el->getIndex() << " HAMMER: " << ih << " fe: " << el->calculateGlobalCoordinate(xsiFE)(0) << " " << el->calculateGlobalCoordinate(xsiFE)(1) << std::endl;

//                     break;
//                 }
//             }
//         }
//         else
//         {
//             insideBlendZone.push_back(false);
//         }
//     }
//     el->setIncidenceOnGlobal(insideBlendZone, incidenceCell, xsiIncidenceCell);
// }

bounded_matrix<double, 2, 2> GlobalSolid::inverseMatrix(const bounded_matrix<double, 2, 2> &matrix)
{
    bounded_matrix<double, 2, 2> inverse;
    double detinv = 1.0 / (matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0));

    inverse(0, 0) = detinv * matrix(1, 1);
    inverse(1, 0) = -detinv * matrix(1, 0);
    inverse(0, 1) = -detinv * matrix(0, 1);
    inverse(1, 1) = detinv * matrix(0, 0);

    return inverse;
}

void GlobalSolid::computeDistanceFromFEBoundary()
{
    for (Node *no : nodes_) //paralelizar utilizando elements_part
    {
        double distanceNode = 1000000.0;
        bounded_vector<double, 2> coord = no->getCurrentCoordinate();

        for (BoundaryElement *bound : boundaryFE_)
        {
            double distance;
            double xsiBoundary = 0.0; //primeira tentativa
            double deltaxsi = 100000.0;
            int cont = 0;
            std::vector<Node *> boundaryConec = bound->getNodes();
            bounded_vector<double, 2> coordBoundary, firstDerivate, secondDerivate;

            while (fabs(deltaxsi) >= 1.0e-06 and xsiBoundary >= -1.0 and xsiBoundary <= 1.0 and cont <= 15)
            {
                matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(xsiBoundary); //PHI, PHI', PHI''
                coordBoundary(0) = 0.0;
                coordBoundary(1) = 0.0;
                firstDerivate(0) = 0.0;
                firstDerivate(1) = 0.0;
                secondDerivate(0) = 0.0;
                secondDerivate(1) = 0.0;
                int aux = 0;
                for (Node *node : boundaryConec)
                {
                    bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();

                    coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
                    coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

                    firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
                    firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

                    secondDerivate(0) += boundaryFunctions(aux, 2) * coordinateNode(0);
                    secondDerivate(1) += boundaryFunctions(aux, 2) * coordinateNode(1);
                    aux = aux + 1;
                }
                double h = -((coord(0) - coordBoundary(0)) * (-firstDerivate(0)) + (coord(1) - coordBoundary(1)) * (-firstDerivate(1)));
                double dh_dxsi = (-firstDerivate(0)) * (-firstDerivate(0)) + (-secondDerivate(0)) * (coord(0) - coordBoundary(0)) +
                                 (-firstDerivate(1)) * (-firstDerivate(1)) + (-secondDerivate(1)) * (coord(1) - coordBoundary(1));

                deltaxsi = h / dh_dxsi;

                xsiBoundary = xsiBoundary + deltaxsi;

                cont++;
            }

            if (xsiBoundary <= -1.0)
            {
                matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(-1.0); //PHI, PHI', PHI''
                int aux = 0;
                coordBoundary(0) = 0.0;
                coordBoundary(1) = 0.0;
                for (Node *node : boundaryConec)
                {
                    bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
                    coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
                    coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

                    aux = aux + 1;
                }
            }
            else if (xsiBoundary >= 1.0)
            {
                matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(1.0); //PHI, PHI', PHI''
                int aux = 0;
                coordBoundary(0) = 0.0;
                coordBoundary(1) = 0.0;
                for (Node *node : boundaryConec)
                {
                    bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
                    coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
                    coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

                    aux = aux + 1;
                }
            }

            distance = sqrt((coord(0) - coordBoundary(0)) * (coord(0) - coordBoundary(0)) + (coord(1) - coordBoundary(1)) * (coord(1) - coordBoundary(1)));

            if (distance < distanceNode)
            {
                distanceNode = distance;
            }
        }
        no->setDistanceToBoundary(distanceNode);
        // std::cout<<"NODE: "<<no->getIndex() <<" DISTANCE: "<<distanceNode<<std::endl;
    }

    // for (Element *el : elements_) //voltar para elements_part
    // {
    //     matrix<double> integrationPoints = el->hammerQuadrature();
    //     std::vector<double> distanceFE;

    //     for (int ip = 0; ip < integrationPoints.size1(); ip++)
    //     {
    //         distanceFE.push_back(100000000.0);
    //         bounded_vector<double, 2> xsi, coordIP; //coordinate of hammer point
    //         xsi(0) = integrationPoints(ip, 0);
    //         xsi(1) = integrationPoints(ip, 1);
    //         coordIP = el->calculateGlobalCoordinate(xsi);
    //         double distance;

    //         for (BoundaryElement *bound : boundaryFE_)
    //         {
    //             double xsiBoundary = 0.0; //primeira tentativa
    //             double deltaxsi = 100000.0;
    //             int cont = 0;
    //             std::vector<Node *> boundaryConec = bound->getNodes();
    //             bounded_vector<double, 2> coordBoundary, firstDerivate, secondDerivate, normal;

    //             while (fabs(deltaxsi) >= 1.0e-06 and xsiBoundary >= -1.0 and xsiBoundary <= 1.0 and cont <= 15)
    //             {
    //                 matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(xsiBoundary); //PHI, PHI', PHI''
    //                 coordBoundary(0) = 0.0;
    //                 coordBoundary(1) = 0.0;
    //                 firstDerivate(0) = 0.0;
    //                 firstDerivate(1) = 0.0;
    //                 secondDerivate(0) = 0.0;
    //                 secondDerivate(1) = 0.0;
    //                 int aux = 0;
    //                 for (Node *node : boundaryConec)
    //                 {
    //                     bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();

    //                     coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
    //                     coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

    //                     firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
    //                     firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

    //                     secondDerivate(0) += boundaryFunctions(aux, 2) * coordinateNode(0);
    //                     secondDerivate(1) += boundaryFunctions(aux, 2) * coordinateNode(1);
    //                     aux = aux + 1;
    //                 }
    //                 double h = -((coordIP(0) - coordBoundary(0)) * (-firstDerivate(0)) + (coordIP(1) - coordBoundary(1)) * (-firstDerivate(1)));
    //                 double dh_dxsi = (-firstDerivate(0)) * (-firstDerivate(0)) + (-secondDerivate(0)) * (coordIP(0) - coordBoundary(0)) +
    //                                  (-firstDerivate(1)) * (-firstDerivate(1)) + (-secondDerivate(1)) * (coordIP(1) - coordBoundary(1));

    //                 deltaxsi = h / dh_dxsi;

    //                 xsiBoundary = xsiBoundary + deltaxsi;

    //                 cont++;
    //             }

    //             if (xsiBoundary <= -1.0)
    //             {
    //                 matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(-1.0); //PHI, PHI', PHI''
    //                 int aux = 0;
    //                 coordBoundary(0) = 0.0;
    //                 coordBoundary(1) = 0.0;
    //                 for (Node *node : boundaryConec)
    //                 {
    //                     bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
    //                     coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
    //                     coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

    //                     aux = aux + 1;
    //                 }
    //             }
    //             else if (xsiBoundary >= 1.0)
    //             {
    //                 matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(1.0); //PHI, PHI', PHI''
    //                 int aux = 0;
    //                 coordBoundary(0) = 0.0;
    //                 coordBoundary(1) = 0.0;
    //                 for (Node *node : boundaryConec)
    //                 {
    //                     bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
    //                     coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
    //                     coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

    //                     aux = aux + 1;
    //                 }
    //             }

    //             distance = sqrt((coordIP(0) - coordBoundary(0)) * (coordIP(0) - coordBoundary(0)) + (coordIP(1) - coordBoundary(1)) * (coordIP(1) - coordBoundary(1)));

    //             if (distance < distanceFE[ip])
    //             {
    //                 distanceFE[ip] = distance;
    //             }
    //         }
    //     }
    //     el->setDistanceFromFEBoundary(distanceFE);
    // }

    for (Cell *cell : cells_) //voltar para cells_part
    {
        matrix<double> integrationPoints = cell->isoQuadrature(quadrature_);
        std::vector<double> distanceFE;

        for (int ip = 0; ip < integrationPoints.size1(); ip++)
        {
            distanceFE.push_back(100000000.0);
            bounded_vector<double, 2> xsi, coordIP; //coordinate of hammer point
            xsi(0) = integrationPoints(ip, 0);
            xsi(1) = integrationPoints(ip, 1);
            coordIP = cell->calculateGlobalCoordinate(xsi);
            double distance;

            for (BoundaryElement *bound : boundaryFE_)
            {
                double xsiBoundary = 0.0; //primeira tentativa
                double deltaxsi = 100000.0;
                int cont = 0;
                std::vector<Node *> boundaryConec = bound->getNodes();
                bounded_vector<double, 2> coordBoundary, firstDerivate, secondDerivate, normal;

                while (fabs(deltaxsi) >= 1.0e-06 and xsiBoundary >= -1.0 and xsiBoundary <= 1.0 and cont <= 15)
                {
                    matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(xsiBoundary); //PHI, PHI', PHI''
                    coordBoundary(0) = 0.0;
                    coordBoundary(1) = 0.0;
                    firstDerivate(0) = 0.0;
                    firstDerivate(1) = 0.0;
                    secondDerivate(0) = 0.0;
                    secondDerivate(1) = 0.0;
                    int aux = 0;
                    for (Node *node : boundaryConec)
                    {
                        bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();

                        coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
                        coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

                        firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
                        firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

                        secondDerivate(0) += boundaryFunctions(aux, 2) * coordinateNode(0);
                        secondDerivate(1) += boundaryFunctions(aux, 2) * coordinateNode(1);
                        aux = aux + 1;
                    }
                    double h = -((coordIP(0) - coordBoundary(0)) * (-firstDerivate(0)) + (coordIP(1) - coordBoundary(1)) * (-firstDerivate(1)));
                    double dh_dxsi = (-firstDerivate(0)) * (-firstDerivate(0)) + (-secondDerivate(0)) * (coordIP(0) - coordBoundary(0)) +
                                     (-firstDerivate(1)) * (-firstDerivate(1)) + (-secondDerivate(1)) * (coordIP(1) - coordBoundary(1));

                    deltaxsi = h / dh_dxsi;

                    xsiBoundary = xsiBoundary + deltaxsi;

                    cont++;
                }

                if (xsiBoundary <= -1.0)
                {
                    matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(-1.0); //PHI, PHI', PHI''
                    int aux = 0;
                    coordBoundary(0) = 0.0;
                    coordBoundary(1) = 0.0;
                    firstDerivate(0) = 0.0;
                    firstDerivate(1) = 0.0;
                    for (Node *node : boundaryConec)
                    {
                        bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
                        coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
                        coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

                        firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
                        firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

                        aux = aux + 1;
                    }
                }
                else if (xsiBoundary >= 1.0)
                {
                    matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(1.0); //PHI, PHI', PHI''
                    int aux = 0;
                    coordBoundary(0) = 0.0;
                    coordBoundary(1) = 0.0;
                    firstDerivate(0) = 0.0;
                    firstDerivate(1) = 0.0;
                    for (Node *node : boundaryConec)
                    {
                        bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
                        coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
                        coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

                        firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
                        firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

                        aux = aux + 1;
                    }
                }

                distance = sqrt((coordIP(0) - coordBoundary(0)) * (coordIP(0) - coordBoundary(0)) + (coordIP(1) - coordBoundary(1)) * (coordIP(1) - coordBoundary(1)));
                normal(0) = firstDerivate(1);
                normal(1) = -firstDerivate(0);
                double aux = (coordIP(0) - coordBoundary(0)) * normal(0) + (coordIP(1) - coordBoundary(1)) * normal(1);

                if (distance < 1.0e-10)
                {
                    distance = 0.0;
                }

                if (distance < fabs(distanceFE[ip]))
                {
                    if (aux > 0.0) //fora do domínio local
                    {
                        distanceFE[ip] = -distance;
                    }
                    else if (aux < 0.0) //dentro do domínio local
                    {
                        distanceFE[ip] = distance;
                    }
                    else if (aux == 0.0 and distance == 0.0) //caso único: exatamente sobre um elemento do contorno
                    {
                        distanceFE[ip] = distance;
                    }
                }
            }
        }
        cell->setDistanceFromFEBoundary(distanceFE);
        // for(int i=0; i<distanceFE.size();i++)
        // {
        //     std::cout<<"CELL: "<<cell->getIndex()<<" DENTRO? "<<distanceFE[i]<<std::endl;
        // }
    }
}

void GlobalSolid::testeParaviewFE()
{
    std::stringstream name;
    name << "hammerPointsFE.vtu";
    std::ofstream file(name.str());

    matrix<double> hammer = elements_[0]->hammerQuadrature();
    int nn = hammer.size1();

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << nn * elements_.size()
         << "\"  NumberOfCells=\"" << nn * elements_.size()
         << "\">"
         << "\n";

    //nodal coordinates
    file << "    <Points>"
         << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">"
         << "\n";
    for (Element *el : elements_)
    {
        matrix<double> hammer = el->hammerQuadrature();
        for (int i = 0; i < hammer.size1(); i++)
        {
            bounded_vector<double, 2> xsi;
            xsi(0) = hammer(i, 0);
            xsi(1) = hammer(i, 1);
            bounded_vector<double, 2> coord = el->calculateGlobalCoordinate(xsi);
            file << coord(0) << " " << coord(1) << " " << 0.0 << "\n";
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
    int aux = 0;
    for (Element *el : elements_)
    {
        for (int i = 0; i < nn; i++)
        {
            file << aux++ << std::endl;
        }
    }
    file << "      </DataArray>"
         << "\n";

    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";
    aux = 0;
    for (Element *el : elements_)
    {
        for (int i = 0; i < nn; i++)
        {
            file << aux++ << std::endl;
        }
    }
    file << "      </DataArray>"
         << "\n";

    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";
    for (Element *el : elements_)
    {
        for (int i = 0; i < nn; i++)
        {
            file << 1 << std::endl;
        }
    }
    file << "      </DataArray>"
         << "\n"
         << "    </Cells>"
         << "\n";

    //nodal results
    file << "    <PointData>"
         << "\n";
    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         << "Name=\"Blend\" format=\"ascii\">"
         << "\n";
    for (Element *el : elements_)
    {
        std::vector<bool> blend = el->ipInsideBlendZone();
        std::vector<double> b = el->getBlendValue();
        int au = 0;
        for (int i = 0; i < blend.size(); i++)
        {
            double tt = 0.0;
            if (blend[i] == true)
            {
                tt = b[au++];
            }
            file << blend[i] << " " << tt << " " << 1.0 - tt << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"6\" "
    //      << "Name=\"Local\" format=\"ascii\">"
    //      << "\n";
    // for (Element *el : elements_)
    // {
    //     matrix<double> hammer = el->hammerQuadrature();
    //     std::vector<bool> blend = el->ipInsideBlendZone();
    //     std::vector<double> dist = el->getDistanceFromFEBoundary();
    //     vector<double> phi;

    //     for (int i = 0; i < nn; i++)
    //     {
    //         phi = el->domainShapeFunction(hammer(i, 0), hammer(i, 1));
    //         double tt = 0.0;
    //         if (blend[i] == true)
    //         {
    //             tt = el->blendFunction(dist[i], blendZoneThickness_);
    //         }
    //         for (int j = 0; j < 6; j++)
    //         {
    //             file << (1.0 - tt) * phi(j) << " ";
    //         }
    //         file << "\n";
    //     }
    // }
    // file << "      </DataArray> "
    //      << "\n";

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"9\" "
    //      << "Name=\"Global\" format=\"ascii\">"
    //      << "\n";
    // for (Element *el : elements_)
    // {
    //     matrix<double> hammer = el->hammerQuadrature();
    //     std::vector<bool> blend = el->ipInsideBlendZone();
    //     std::vector<double> dist = el->getDistanceFromFEBoundary();
    //     std::vector<Cell *> cell = el->incidenceCell();
    //     std::vector<bounded_vector<double, 2>> xsiGlobal = el->xsiIncidenceCell();

    //     //vector<double> phi;
    //     int aa = 0;
    //     for (int i = 0; i < hammer.size1(); i++)
    //     {
    //         double tt = 0.0;
    //         if (blend[i] == true)
    //         {
    //             tt = el->blendFunction(dist[i], blendZoneThickness_);
    //             std::pair<vector<double>, matrix<double>> phi = cell[aa]->shapeFunctionAndDerivates(xsiGlobal[aa]); //PHI, PHI'
    //             phi.first = tt * phi.first;

    //             for (int j = 0; j < phi.first.size(); j++)
    //             {
    //                 file << phi.first(j) << " ";
    //             }

    //             aa++;
    //         }
    //         else
    //         {
    //             for (int j = 0; j < 9; j++)
    //             {
    //                 file << 0.0 << " ";
    //             }
    //         }
    //         file << "\n";
    //     }
    // }
    // file << "      </DataArray> "
    //      << "\n";

    file << "    </PointData>"
         << "\n";
    //elemental results
    file << "    <CellData>"
         << "\n";

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //      << "Name=\"Process\" format=\"ascii\">" << std::endl;

    // for (Element *el : elements_)
    // {
    //     file << elementPartition_[el->getIndex()] << "\n";
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

bounded_vector<double, 2> GlobalSolid::blendFunction(const double &DA)
{
    bounded_vector<double, 2> blend; //(b, db_dDA)
    double aux = DA / blendZoneThickness_;
    blend(0) = 2.0 * aux * aux * aux - 3.0 * aux * aux + 1.0;
    blend(1) = (6.0 / blendZoneThickness_) * aux * aux - (6.0 / blendZoneThickness_) * aux;

    return blend;
}

vector<double> GlobalSolid::diagonalMassMatrix()
{
    vector<double> diagonal(cpnumber_, 0.0); //ver para paralelizar com vetor do PETSC;s
    for (Cell *cell : cells_)
    {
        vector<double> mass = cell->diagonalMass(quadrature_);
        std::vector<ControlPoint *> conec = cell->getControlPoints();
        for (int i = 0; i < conec.size(); i++)
        {
            diagonal(conec[i]->getIndex()) += mass(i);
        }
    }
    for (Element *el : elements_)
    {
        vector<double> mass = el->diagonalMass();
        std::vector<int> freedom = el->getFreedomDegree();
        for (int i = 0; i < freedom.size(); i++)
        {
            diagonal(freedom[i]) += mass(i);
        }
    }

    // for (int i = 0; i < controlPoints_.size(); i++)
    // {
    //     if (diagonal(i) < 5.0e-03)
    //     {
    //         DirichletCondition *dir0 = new DirichletCondition(controlPoints_[i], 0, 0.0);
    //         cpOutsideDomain_.push_back(dir0);

    //         DirichletCondition *dir1 = new DirichletCondition(controlPoints_[i], 1, 0.0);
    //         cpOutsideDomain_.push_back(dir1);
    //     }
    // }
    return diagonal;
}

void GlobalSolid::checkInactivesCPs()
{
    vector<double> pseudoDiagonal = diagonalMassMatrix();

    for (int i = 0; i < controlPoints_.size(); i++)
    {
        if (pseudoDiagonal(i) < 5.0e-06)
        {
            ControlPoint *cp=controlPoints_[i];
            bounded_vector<double, 2> coord = cp->getCurrentCoordinate();
            DirichletCondition *dir0 = new DirichletCondition(cp, 0, 0.0);
            cpOutsideDomain_.push_back(dir0);

            DirichletCondition *dir1 = new DirichletCondition(cp, 1, 0.0);
            cpOutsideDomain_.push_back(dir1);
            std::cout<<"CP "<<cp->getIndex()<<" foi desativado..."<<std::endl;

            for (Element *el : elements_)
            {
                bounded_vector<double, 2> xsiFE, coordFE, deltaxsi;
                deltaxsi(0) = 100000.0;
                deltaxsi(1) = 100000.0;
                xsiFE(0) = 0.5; //first attempt xs1
                xsiFE(1) = 0.5; //first attempt xs2

                std::vector<Node *> conec = el->getConnection();

                int cont = 0;
                while (deltaxsi(0) >= 1.0e-06 and deltaxsi(1) >= 1.0e-06 and
                       xsiFE(0) >= 0.0 and xsiFE(0) <= 1.0 and
                       xsiFE(1) >= 0.0 and xsiFE(1) <= 1.0 and cont <= 15)
                {

                    coordFE = el->calculateGlobalCoordinate(xsiFE);

                    bounded_matrix<double, 2, 2> jacobian = el->referenceJacobianMatrix(xsiFE(0), xsiFE(1));
                    bounded_matrix<double, 2, 2> inverse = inverseMatrix(jacobian);
                    deltaxsi = prod(inverse, coord - coordFE);

                    xsiFE = xsiFE + deltaxsi;

                    cont++;
                }

                if (xsiFE(0) >= 0.0 and xsiFE(0) <= 1.0 and //ponto de hammer encontrou uma célula
                    xsiFE(1) >= 0.0 and xsiFE(1) <= 1.0 and
                    xsiFE(0) + xsiFE(1) <= 1.0) //só vai funcionar para triângulos não muito deformados...
                {
                    CP_BlendingZone *inside = new CP_BlendingZone(cp, el, xsiFE);
                    cpsInsideBledingZone_.push_back(inside);

                    break;
                }
            }
        }
    }

    /*Looping para:
        -Encontrar pontos de controle dentro da malha local;
            -Elemento em que o ponto de controle está inserido;
            -Xsis correspondentes ao ponto de controle no elemento;
        -Encontrar pontos de controle fora do domínio global;
            -Adicionar condição de contorno;
    */
    // if (cpOutsideDomain_.size() > 1)
    // {
    //     cpOutsideDomain_.erase(cpOutsideDomain_.begin(), cpOutsideDomain_.begin() + cpOutsideDomain_.size());
    // }
    // for (ControlPoint *cp : controlPoints_)
    // {
    //     bounded_vector<double, 2> coord = cp->getCurrentCoordinate() / cp->getWeight();
    //     double distanceCP = 1000000.0;
    //     double distance;

    //     for (BoundaryElement *bound : boundaryFE_)
    //     {
    //         double xsiBoundary = 0.0; //primeira tentativa
    //         double deltaxsi = 100000.0;
    //         int cont = 0;
    //         std::vector<Node *> boundaryConec = bound->getNodes();
    //         bounded_vector<double, 2> coordBoundary, firstDerivate, secondDerivate, normal;

    //         while (fabs(deltaxsi) >= 1.0e-06 and xsiBoundary >= -1.0 and xsiBoundary <= 1.0 and cont <= 15)
    //         {
    //             matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(xsiBoundary); //PHI, PHI', PHI''
    //             coordBoundary(0) = 0.0;
    //             coordBoundary(1) = 0.0;
    //             firstDerivate(0) = 0.0;
    //             firstDerivate(1) = 0.0;
    //             secondDerivate(0) = 0.0;
    //             secondDerivate(1) = 0.0;
    //             int aux = 0;
    //             for (Node *node : boundaryConec)
    //             {
    //                 bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();

    //                 coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
    //                 coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

    //                 firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
    //                 firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

    //                 secondDerivate(0) += boundaryFunctions(aux, 2) * coordinateNode(0);
    //                 secondDerivate(1) += boundaryFunctions(aux, 2) * coordinateNode(1);
    //                 aux = aux + 1;
    //             }
    //             double h = -((coord(0) - coordBoundary(0)) * (-firstDerivate(0)) + (coord(1) - coordBoundary(1)) * (-firstDerivate(1)));
    //             double dh_dxsi = (-firstDerivate(0)) * (-firstDerivate(0)) + (-secondDerivate(0)) * (coord(0) - coordBoundary(0)) +
    //                              (-firstDerivate(1)) * (-firstDerivate(1)) + (-secondDerivate(1)) * (coord(1) - coordBoundary(1));

    //             deltaxsi = h / dh_dxsi;

    //             xsiBoundary = xsiBoundary + deltaxsi;

    //             cont++;
    //         }

    //         if (xsiBoundary <= -1.0)
    //         {
    //             matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(-1.0); //PHI, PHI', PHI''
    //             int aux = 0;
    //             coordBoundary(0) = 0.0;
    //             coordBoundary(1) = 0.0;
    //             firstDerivate(0) = 0.0;
    //             firstDerivate(1) = 0.0;
    //             for (Node *node : boundaryConec)
    //             {
    //                 bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
    //                 coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
    //                 coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

    //                 firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
    //                 firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

    //                 aux = aux + 1;
    //             }
    //         }
    //         else if (xsiBoundary >= 1.0)
    //         {
    //             matrix<double> boundaryFunctions = bound->shapeFunctionsAndDerivates(1.0); //PHI, PHI', PHI''
    //             int aux = 0;
    //             coordBoundary(0) = 0.0;
    //             coordBoundary(1) = 0.0;
    //             firstDerivate(0) = 0.0;
    //             firstDerivate(1) = 0.0;
    //             for (Node *node : boundaryConec)
    //             {
    //                 bounded_vector<double, 2> coordinateNode = node->getCurrentCoordinate();
    //                 coordBoundary(0) += boundaryFunctions(aux, 0) * coordinateNode(0);
    //                 coordBoundary(1) += boundaryFunctions(aux, 0) * coordinateNode(1);

    //                 firstDerivate(0) += boundaryFunctions(aux, 1) * coordinateNode(0);
    //                 firstDerivate(1) += boundaryFunctions(aux, 1) * coordinateNode(1);

    //                 aux = aux + 1;
    //             }
    //         }

    //         distance = sqrt((coord(0) - coordBoundary(0)) * (coord(0) - coordBoundary(0)) + (coord(1) - coordBoundary(1)) * (coord(1) - coordBoundary(1)));
    //         normal(0) = firstDerivate(1);
    //         normal(1) = -firstDerivate(0);
    //         double aux = (coord(0) - coordBoundary(0)) * normal(0) + (coord(1) - coordBoundary(1)) * normal(1);

    //         if (distance < 1.0e-10)
    //         {
    //             distance = 0.0;
    //         }

    //         if (distance < fabs(distanceCP))
    //         {
    //             if (aux > 0.0) //fora do domínio local
    //             {
    //                 distanceCP = -distance;
    //             }
    //             else if (aux < 0.0) //dentro do domínio local
    //             {
    //                 distanceCP = distance;
    //             }
    //             else if (aux == 0.0 and distance == 0.0) //caso único: exatamente sobre um elemento do contorno
    //             {
    //                 distanceCP = distance;
    //             }
    //         }
    //     }
        
    // }
}