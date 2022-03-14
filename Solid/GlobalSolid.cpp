#include "GlobalSolid.h"

GlobalSolid::GlobalSolid()
{
    char buff[FILENAME_MAX];
    getCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    current_working_dir_ = current_working_dir;
    overlappingAnalysis_ = false;
    plotHammerPointsInBlendingZone_ = false;
    analysisOfCrackPropagation_ = false;
    plotSIFs_ = true;
    viewNewMesh_ = false;
    exportUndeformedMesh_ = false;
    localOverlapping_ = false;
    quarterPointElement_ = false;
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

int GlobalSolid::getNodesNumber()
{
    return nodes_.size();
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

// void GlobalSolid::addMesh(const int &index, const int &indexMaterial, const double &thickness, const std::string &elementType)
// {
//     Mesh *mesh = new Mesh(index, materials_[indexMaterial], thickness, elementType);
//     meshes_.push_back(mesh);
// }

void GlobalSolid::addNode(const int &index, const int &indexFE,
                          const bounded_vector<double, 2> &initialCoordinate)
{
    Node *node = new Node(index, indexFE, initialCoordinate);
    nodes_.push_back(node);
}

void GlobalSolid::dataReading(const std::string &inputParameters, const std::string &inputProperties, const std::string &inputMeshIso, const bool &rhino)
{
    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    std::string line;

    //Opening parameters file
    std::ifstream parameters(inputParameters);
    //Reading analysis parameters
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
        parameters >> quadrature_ >> quadratureBlendingZone_;
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
        parameters >> quadrature_ >> quadratureBlendingZone_;
    }
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    parameters >> hammerPoints_ >> hammerPointsBlendZone_;
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    std::getline(parameters, line);
    parameters >> orderParaview_;
    //Closing parameters file
    parameters.close();

    //Opening properties file
    std::ifstream properties(inputProperties);
    //Reading solid properties
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
    //Closing properties file
    properties.close();

    FILE *iso;
    iso = fopen(inputMeshIso.c_str(), "r");
    if (iso != NULL) //the file exists
    {
        fclose(iso);

        //Opening the file
        std::ifstream isogeometric(inputMeshIso);
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
                            vector<double> contribuition = cells_[j]->computeDistribuitedLoads(value, 10, 0);
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
                            vector<double> contribuition = cells_[j]->computeDistribuitedLoads(value, 10, 1);
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
                            vector<double> contribuition = cells_[j]->computeDistribuitedLoads(value, 10, 2);
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
                            vector<double> contribuition = cells_[j]->computeDistribuitedLoads(value, 10, 3);
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

        //Closing the file
        isogeometric.close();

        int cpindex = 0;
        double dist = 1.0e-06;

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

                    //bounded_vector<double, 2> aux = coord_global - coord_patch;
                    //double dist2 = aux(0) * aux(0) + aux(1) * aux(1);

                    if (abs(coord_global(0) - coord_patch(0)) <= dist and abs(coord_global(1) - coord_patch(1)) <= dist)
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

        if (size > 1)
        {
            ISOdomainDecompositionMETIS();
            for (int i = 0; i < cells_.size(); i++)
            {
                if (rank == cellPartition_[i])
                {
                    cells_part.push_back(cells_[i]);
                }
            }
        }
        else
        {
            cellPartition_ = new idx_t[cells_.size()];
            pointsPartition_ = new idx_t[controlPoints_.size()];
            for (int i = 0; i < cells_.size(); i++)
            {
                cells_part.push_back(cells_[i]);
            }
        }
    }
}

std::vector<std::string> split(std::string str, std::string delim)
{
    std::istringstream is(str);
    std::vector<std::string> values;
    std::string token;
    while (getline(is, token, ' '))
        values.push_back(token);
    return values;
}

void GlobalSolid::dataFromGmsh(const std::string &inputGmesh, const std::string &elementType, Geometry *geometry)
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    int nentities, nmesh, nnode, nelem, auxEl, auxBo;

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

    //opening .msh file
    std::ifstream feMesh(inputGmesh);
    std::string line;
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);
    std::getline(feMesh, line);

    //reading physical entities
    feMesh >> nentities;
    std::getline(feMesh, line);
    std::unordered_map<int, std::string> physicalEntities;
    physicalEntities.reserve(nentities);
    for (int i = 0; i < nentities; i++)
    {
        std::getline(feMesh, line);
        std::vector<std::string> tokens = split(line, " ");
        int index;
        std::istringstream(tokens[1]) >> index;
        physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
    }

    std::unordered_map<std::string, PlaneSurface *> planeSurfaces = geometry->getPlaneSurfaces();

    int cont = 0;
    for (std::unordered_map<int, std::string>::const_iterator entities = physicalEntities.begin(); entities != physicalEntities.end(); entities++)
    {
        std::string name = entities->second;
        if (name[0] == 's') //surface
        {
            Mesh *mesh = new Mesh(cont++, materials_[planeSurfaces[entities->second]->getIndexMaterial()], planeSurfaces[entities->second]->getThickness(), elementType);
            meshes_[entities->second] = mesh;
        }
        else if (name[0] == 'l' or name[0] == 'j') //line
        {
            std::vector<BoundaryElement *> aux;
            finiteElementBoundary_.insert(std::unordered_map<std::string, std::vector<BoundaryElement *>>::value_type(name, aux));
        }
        else if (name[0] == 'c') //crack
        {
            std::vector<BoundaryElement *> aux;
            finiteElementBoundary_.insert(std::unordered_map<std::string, std::vector<BoundaryElement *>>::value_type(name + "0", aux));
            finiteElementBoundary_.insert(std::unordered_map<std::string, std::vector<BoundaryElement *>>::value_type(name + "1", aux));
        }
    }

    std::getline(feMesh, line);
    std::getline(feMesh, line);

    //reading nodes
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

    //reading elements
    feMesh >> nelem;
    cont = 0;
    int auxconec, entitie, auxentitie;
    int elementCont = 0;
    int num = 0;
    int auxn = 0;
    for (int ielem = 0; ielem < nelem; ielem++)
    {
        feMesh >> trash >> trash >> trash >> entitie >> auxentitie;
        std::string name = physicalEntities[entitie];
        if (name[0] == 's') //2D elements
        {
            std::vector<Node *> conec;

            for (int in = 0; in < auxEl; in++)
            {
                feMesh >> auxconec;
                conec.push_back(nodes_[auxconec - 1]);
            }

            Element *el = new Element(elementCont++, meshes_[name], conec);
            elements_.push_back(el);
        }
        else if (name[0] == 'l' or name[0] == 'j')
        {
            std::vector<Node *> conec;

            for (int in = 0; in < auxBo; in++)
            {
                feMesh >> auxconec;
                conec.push_back(nodes_[auxconec - 1]);
            }

            if (name[0] == 'j')
            {
                auxn = 1;
            }

            BoundaryElement *bfe = new BoundaryElement(num++, conec);
            finiteElementBoundary_[name].push_back(bfe);
        }
        else if (name[0] == 'c')
        {
            std::vector<Node *> conec;

            for (int in = 0; in < auxBo; in++)
            {
                feMesh >> auxconec;
                conec.push_back(nodes_[auxconec - 1]);
            }

            BoundaryElement *bfe = new BoundaryElement(num++, conec);

            if (finiteElementBoundary_[name + "0"].size() == 0 or auxn == 0)
            {
                finiteElementBoundary_[name + "0"].push_back(bfe);
                auxn == 0;
            }
            else
            {
                finiteElementBoundary_[name + "1"].push_back(bfe);
            }
        }
        else if (name[0] == 'p' and geometry->getPoints().count(name) >= 1)
        {
            feMesh >> auxconec;
            geometry->getPoint(name)->addNodeToPoint(nodes_[auxconec - 1]);
        }
        std::getline(feMesh, line);
    }

    for (int i = 0; i < geometry->getNumberOfCrackes(); i++)
    {
        std::vector<BoundaryElement *> integralj = finiteElementBoundary_["j0"];
        double radius = geometry_->getRadiusJintegral();
        bounded_vector<double, 2> tip = geometry_->getCrack("c" + std::to_string(i))->getLastPoint()->getCoordinates();

        for (Element *el : elements_)
        {
            int teste = -1;
            for (BoundaryElement *bound : integralj)
            {
                std::vector<Node *> conecbound = bound->getNodes();
                for (int i = 0; i < 3; i++)
                {
                    if (conecbound == el->getSideNodes(i))
                    {
                        teste = i;
                        JintegralElements_.push_back(el);
                        break;
                    }
                }
            }

            // if (problemType_ == "DYNAMIC")
            // {
            //     bounded_vector<double, 2> coord0, coord1, coord2;
            //     coord0 = el->getConnection()[0]->getInitialCoordinate();
            //     coord1 = el->getConnection()[1]->getInitialCoordinate();
            //     coord2 = el->getConnection()[2]->getInitialCoordinate();

            //     if ((norm_2(coord0 - tip) < radius + 1.0e-04) and (norm_2(coord1 - tip) < radius + 1.0e-04) and (norm_2(coord2 - tip) < radius + 1.0e-04))
            //     {
            //         dynamicJintegralElements_.push_back(el);
            //         teste += 11;
            //     }
            // }
            elementsSideJintegral_.push_back(teste);
        }
    }

    //Closing the .msh file
    feMesh.close();
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
    mirror << "NUMBER OF INTEGRATION POINTS IN BLENDINGZONE: " << quadratureBlendingZone_ * quadratureBlendingZone_ << std::endl;
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
    for (auto me : meshes_)
    {
        double young, poisson, densisity;
        me.second->getMaterial()->setProperties(young, poisson, densisity);
        mirror << "SURFACE: " << me.second->getIndex() << std::endl;
        mirror << "THICKNESS: " << me.second->getThickness() << std::endl;
        mirror << "YOUNG: " << young << std::endl;
        mirror << "POISSON: " << poisson << std::endl;
        mirror << "DENSITY: " << densisity << std::endl;
        mirror << "ELEMENT TYPE: " << me.second->getElementType() << std::endl;
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
    Cell *el = cells_[48];
    std::pair<vector<double>, matrix<double>> elementMatrices;
    elementMatrices = el->cellContributions(planeState_, "STATIC", numberOfSteps_, numberOfSteps_, 1.0, 0.25, 0.5, quadrature_);
    // std::vector<int> freedom = el->getFreedomDegree();
    // int num = freedom.size();

    for (size_t i = 0; i < el->getControlPoints().size(); i++)
    {
        for (size_t j = 0; j < el->getControlPoints().size() * 2; j++)
        {
            // std::cout << elementMatrices.second(i, j) << " ";
        }
        // std::cout << std::endl;
        // std::cout << elementMatrices.first(i) << std::endl;
        std::cout << el->getControlPoint(i)->getIndex() * 2 << std::endl;
        std::cout << el->getControlPoint(i)->getIndex() * 2 + 1 << std::endl;
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
    PetscInt Istart, Iend, Idof, Idof2, Ione, iterations, *dof;
    KSP ksp;
    PC pc;
    VecScatter ctx;
    PetscScalar val, value;

    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    std::stringstream fileName;
    if (rank == 0 and plotSIFs_ == true)
    {
        fileName << "tangentMatrix"
                 << ".txt";
    }
    std::ofstream sifFile(fileName.str());

    if (rank == 0 and cells_.size() > 0)
    {
        if (analysisOfCrackPropagation_ == true)
        {
            exportToParaviewISO(std::to_string(0) + "_" + std::to_string(0));
        }
        else
        {
            exportToParaviewISO(std::to_string(0));
        }
    }

    if (rank == 0 and elements_.size() > 0)
    {
        if (analysisOfCrackPropagation_ == true)
        {
            exportToParaviewFEM(std::to_string(0) + "_" + std::to_string(0));
        }
        else
        {
            exportToParaviewFEM(std::to_string(0));
        }
    }

    if (overlappingAnalysis_ == true)
    {
        if (localOverlapping_ == true)
        {
            setBlendingZoneLines();
        }
        computeDistanceFromFEBoundary();
        incidenceLocalxGlobal();
        checkInactivesCPandNode();
        if (size > 1)
            shareDataBetweenRanks();

        if (plotHammerPointsInBlendingZone_ == true)
        {
            exportToParaviewHammerPoints();
        }
    }

    double length = 0.0;
    int auxprop = 0;

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

    for (int loadStep = 1; loadStep <= numberOfSteps_; loadStep++)
    {
        boost::posix_time::ptime t1 =
            boost::posix_time::microsec_clock::local_time();
        if (rank == 0 and analysisOfCrackPropagation_ == false)
        {
            std::cout << "------------------------- LOAD STEP = "
                      << loadStep << " -------------------------\n";
        }
        int propagation = 0;
        bool testPropagation = true;

        double valueForce;

        while (testPropagation == true)
        {
            if (rank == 0 and analysisOfCrackPropagation_ == true)
            {
                std::cout << "------------ LOAD STEP = "
                          << loadStep << " ------ PROPAGATION = " << propagation << " ------------\n";

                Crack *c = geometry_->getCrackes()["c0"];
                std::vector<Point *> points = c->getPoints();

                for (Point *p : points)
                {
                    bounded_vector<double, 2> coord = p->getCoordinates();

                    std::cout << coord(0) << " " << coord(1) << std::endl;
                }
            }

            int ndir = dirichletConditions_.size() + dirichletConditionsFE_.size() + 2 * inactiveCPandNode_.size();
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
                // std::cout<<"INDEX NODE "<<indexNode<<"; DIRECTION: "<<direction<<std::endl;
            }
            for (int index : inactiveCPandNode_)
            {
                dof[idir++] = 2 * index;
                dof[idir++] = 2 * index + 1;
            }

            for (int iteration = 0; iteration < maximumOfIteration_; iteration++) //definir o máximo de interações por passo de carga
            {
                MPI_Barrier(PETSC_COMM_WORLD);

                //Create PETSc sparse parallel matrix
                ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                    2 * cpnumber_, 2 * cpnumber_,
                                    1800, NULL, 3000, NULL, &A);
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
                    elementMatrices = el->elementContributions(planeState_, "STATIC", loadStep, numberOfSteps_, 1.0, 0.25, 0.5, hammerPoints_, hammerPointsBlendZone_);
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

                double menorDP = 1.0e240;
                if (size == 1)
                {
                    for (int i = 0; i < 2 * cpnumber_; i++)
                    {
                        Ione = 1;
                        Idof = i;
                        ierr = MatGetValues(A, Ione, &Idof, Ione, &Idof, &val);
                        CHKERRQ(ierr);
                        if (val < menorDP)
                        {
                            menorDP = val;
                        }
                    }
                }

                // menorDP = 1.0;

                if (size == 1)
                {

                    // sifFile.precision(5);
                    // sifFile << std::scientific;

                    sifFile << "mat = SparseArray[{";
                    // sifFile << "mat = {";
                    Ione = 1;
                    for (int i = 0; i < 2 * cpnumber_; i++)
                    {
                        Idof = i;
                        // sifFile << "{";
                        for (int j = 0; j < 2 * cpnumber_; j++)
                        {
                            Idof2 = j;
                            ierr = MatGetValues(A, Ione, &Idof, Ione, &Idof2, &val);
                            CHKERRQ(ierr);
                            // sifFile << val;
                            // if (j != 2 * cpnumber_ - 1)
                            // {
                            //     sifFile << ", ";
                            // }
                            // else
                            // {
                            //     if (i != 2 * cpnumber_ - 1)
                            //     {
                            //         sifFile << "},\n";
                            //     }
                            //     else
                            //     {
                            //         sifFile << "}";
                            //     }
                            // }
                            if (fabs(val) > 1.0e-25)
                            {
                                sifFile << "{" << i + 1 << ", " << j + 1 << "} -> " << val / menorDP;
                                if (i == j and i == 2 * cpnumber_ - 1)
                                {
                                    sifFile << "}];\nSingularValueList[mat]";
                                }
                                else
                                {
                                    sifFile << ",";
                                }
                            }
                            else if (fabs(val) > 0.0)
                            {
                                std::cout << "VALORES DESCARTADOS... VERIFICAR! " << val << "\n";
                            }

                            // if(i==j)
                            // {
                            //     std::cout<<val<<std::endl;
                            // }
                        }
                        sifFile << "\n";
                    }
                }

                // MatView(A, PETSC_VIEWER_STDOUT_WORLD);
                // CHKERRQ(ierr);

                //Create KSP context to solve the linear system
                ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
                CHKERRQ(ierr);
                ierr = KSPSetOperators(ksp, A, A);
                CHKERRQ(ierr);

// ////Solve using GMRES
// ierr = KSPSetTolerances(ksp, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, 5000);
// CHKERRQ(ierr);
// ierr = KSPGetPC(ksp, &pc);
// ierr = PCSetType(pc, PCNONE);
// ierr = KSPSetType(ksp, KSPDGMRES);
// CHKERRQ(ierr);
// ierr = KSPGMRESSetRestart(ksp, 1000);
// CHKERRQ(ierr);

//Solve using MUMPS
#if defined(PETSC_HAVE_MUMPS)
                ierr = KSPSetType(ksp, KSPPREONLY);
                ierr = KSPGetPC(ksp, &pc);
                ierr = PCSetType(pc, PCLU);
#endif
                ierr = KSPSetFromOptions(ksp);
                CHKERRQ(ierr);
                ierr = KSPSetUp(ksp);

                //Solve linear systemll0
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

                // for (InactiveNode *no : inactiveNode_)
                // {
                //     no->interpolateGlobalCoordinate();
                // }
                // for (InactiveCP *cp : inactiveCP_)
                // {
                //     cp->interpolateGlobalCoordinate();
                // }

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

            if (analysisOfCrackPropagation_ == true)
            {
                if (rank == 0 and cells_.size() > 0)
                {
                    exportToParaviewISO(std::to_string(loadStep) + "_" + std::to_string(propagation));
                    // exportToParaviewAuxiliarMesh(std::to_string(loadStep) + "_" + std::to_string(propagation));
                }
                else if (rank == 0 and elements_.size() > 0)
                {
                    stressCalculateFEM();
                    exportToParaviewFEM(std::to_string(loadStep) + "_" + std::to_string(propagation));
                }

                if (rank == 0 and plotSIFs_ == true)
                {
                    bounded_vector<double, 2> sifs = getSIFs();
                    sifFile << sifs(0) << " " << sifs(1) << std::endl;
                    length += crackLengthPropagation_;
                }

                verifyCrackPropagation(testPropagation);

                if (testPropagation == true)
                {
                    if (rank == 1 and exportUndeformedMesh_ == true)
                    {
                        exportToParaviewFEM(std::to_string(0) + "_" + std::to_string(++auxprop));
                    }

                    for (ControlPoint *cp : controlPoints_)
                    {
                        cp->setCurrentCoordinate(cp->getInitialCoordinate());
                    }

                    if (overlappingAnalysis_ == true)
                    {
                        if (localOverlapping_ == true)
                        {
                            setBlendingZoneLines();
                        }
                        else
                        {
                            addBlendingZone(linesOfBlendingZone_, blendingZoneThickness_, blendingFunctionOrder_, plotHammerPointsInBlendingZone_);
                        }

                        computeDistanceFromFEBoundary();
                        incidenceLocalxGlobal();
                        checkInactivesCPandNode();
                        shareDataBetweenRanks();

                        if (plotHammerPointsInBlendingZone_ == true)
                        {
                            exportToParaviewHammerPoints();
                        }
                    }

                    initialNorm = 0.0;
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
                }

                propagation += 1;
            }
            else
            {
                // for (Node *n : nodes_)
                // {
                //     double sigmaA = -800.0, sigmaB = 600.0;
                //     double young, poisson, density;
                //     materials_[0]->setProperties(young, poisson, density);
                //     double K1, K2;
                //     if (planeState_ == "EPD")
                //     {
                //         K1 = (1 - poisson * poisson) / young;
                //         K2 = poisson / (1 - poisson);
                //     }
                //     else
                //     {
                //     }
                //     double a = 1.0;
                //     double b = 1.5;
                //     double C = (a * a * b * b * (sigmaA - sigmaB)) / (b * b - a * a);
                //     double b2 = (sigmaB * b * b - sigmaA * a * a) / (b * b - a * a);

                //     bounded_vector<double, 2> coord = n->getInitialCoordinate();
                //     double tetaa = atan2(coord(1), coord(0));
                //     double radius = norm_2(coord);
                //     double urr = K1 * (b2 * (1 - K2) * radius - C / radius * (1 + K2));
                //     n->incrementCurrentCoordinate(0, cos(tetaa) * urr);
                //     n->incrementCurrentCoordinate(1, sin(tetaa) * urr);
                // }
                // for (Node *n : nodes_)
                // {
                //     double a = 0.25;
                //     double S = 1000.0;
                //     double young, poisson, density;
                //     materials_[0]->setProperties(young, poisson, density);
                //     double k = (3.0 - poisson) / (1.0 + poisson);
                //     double mi = young / (2.0 * (1.0 + poisson));

                //     bounded_vector<double, 2> coord = n->getInitialCoordinate();
                //     double theta = atan2(coord(1), coord(0));
                //     double radius = norm_2(coord);
                //     bounded_vector<double, 2> teoDisp, dif;
                //     teoDisp(0) = S * a / (8.0 * mi) * (radius / a * (k + 1.0) * cos(theta) + 2.0 * a / radius * ((1.0 + k) * cos(theta) + cos(3.0 * theta)) - 2.0 * a * a * a / (radius * radius * radius) * cos(3.0 * theta));
                //     teoDisp(1) = S * a / (8.0 * mi) * (radius / a * (k - 3.0) * sin(theta) + 2.0 * a / radius * ((1.0 - k) * sin(theta) + sin(3.0 * theta)) - 2.0 * a * a * a / (radius * radius * radius) * sin(3.0 * theta));
                //     dif = (n->getCurrentDisplacement() - (teoDisp)) + coord;
                //     n->setCurrentCoordinate(dif);
                // }
                // std::cout << "PASSOU\n";
                if (rank == 0 and cells_.size() > 0)
                {
                    exportToParaviewISO(std::to_string(loadStep));
                }
                if (rank == 0 and elements_.size() > 0)
                {
                    stressCalculateFEM();
                    exportToParaviewFEM(std::to_string(loadStep));
                }

                // if (rank == 0)
                // {
                //     // sifFile << valueForce << std::en
                //     sifFile << controlPoints_[0]->getCurrentCoordinate()(0) << " " << controlPoints_[0]->getCurrentCoordinate()(1) << std::endl;
                // }

                // if (rank == 0 and plotSIFs_ == true)
                // {
                //     bounded_vector<double, 2> sifs = getSIFs();
                //     sifFile << sifs(0) << " " << sifs(1) << " " << nodes_[2]->getCurrentCoordinate()(1) << " " << nodes_[3]->getCurrentCoordinate()(1) << std::endl;
                //     length += crackLengthPropagation_;
                // }
                testPropagation = false;
            }
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    // remesh(teste);

    // if (rank == 1 and elements_.size() > 0)
    // {
    //     exportToParaviewFEM(10000);
    // }

    PetscFree(dof);
    return 0;
}

void GlobalSolid::exportToParaviewISO(const std::string &nameF)
{
    matrix<double> qxsi2 = coordinatesForInterpolation(orderParaview_);

    int auxxx = (orderParaview_ + 1) * (orderParaview_ + 1);
    int nIsoPoints = auxxx * cells_.size();
    std::stringstream text;
    text << "outputISO" << nameF << ".vtu";
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
                coord(0) += phi2_(i) * x(0);
                coord(1) += phi2_(i) * x(1);
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

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"TheoreticalDisplacement\" format=\"ascii\">"
         << "\n";
    for (Cell *cell : cells_)
    {
        double sigmaA = -8.0, sigmaB = 6.0;
        double young, poisson, density;
        materials_[0]->setProperties(young, poisson, density);
        double K1, K2;
        if (planeState_ == "EPD")
        {
            K1 = (1.0 - poisson * poisson) / young;
            K2 = poisson / (1.0 - poisson);
        }
        else
        {
        }
        double a = 1.0;
        double b = 2.0;
        double C = (a * a * b * b * (sigmaA - sigmaB)) / (b * b - a * a);
        double b2 = (sigmaB * b * b - sigmaA * a * a) / (b * b - a * a);

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
                disp += phi2_(i) * connection[i]->getInitialCoordinate();
            }
            double tetaa = atan2(disp(1), disp(0));
            double radius = norm_2(disp);
            double urr = K1 * (b2 * (1.0 - K2) * radius - C / radius * (1.0 + K2));
            file << cos(tetaa) * urr << " " << sin(tetaa) * urr << std::endl;
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

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"Index\" format=\"ascii\">" << std::endl;
    for (Cell *el : cells_)
    {
        file << el->getIndex() << "\n";
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

void GlobalSolid::exportToParaviewFEM(const std::string &nameF)
{
    std::stringstream name;
    name << "outputFEM" << nameF << ".vtu";
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
        bounded_vector<double, 2> coordNode = no->getCurrentCoordinate();
        int cellIndex = no->getCellIndex();
        if (cellIndex == -1)
        {
            file << coordNode(0) << " " << coordNode(1) << " " << 0.1 << "\n";
        }
        else
        {
            Cell *cell = cells_[cellIndex];
            bounded_vector<double, 2> xsiglobal = no->getXsisGlobal();
            bounded_vector<double, 2> blendF = blendFunction(no->getDistanceToBoundary());

            std::vector<ControlPoint *> connection = cell->getControlPoints();
            vector<double> wpc2(connection.size());
            bounded_vector<int, 2> INC_;
            INC_ = connection[connection.size() - 1]->getINC();
            for (int i = 0; i < connection.size(); i++)
            {
                wpc2(i) = connection[i]->getWeight();
            }
            vector<double> phi = cell->shapeFunction(xsiglobal, wpc2, INC_);

            bounded_vector<double, 2> coordInter;
            coordInter = (1.0 - blendF(0)) * coordNode;

            for (int i = 0; i < connection.size(); i++)
            {
                bounded_vector<double, 2> coordCP = connection[i]->getCurrentCoordinate();
                coordInter += blendF(0) * phi(i) * coordCP;
            }

            file << coordInter(0) << " " << coordInter(1) << " " << 0.1 << "\n";
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
        bounded_vector<double, 2> dispNode = n->getCurrentDisplacement();

        int cellIndex = n->getCellIndex();
        if (cellIndex == -1)
        {
            file << dispNode(0) << " " << dispNode(1) << "\n";
        }
        else
        {
            Cell *cell = cells_[cellIndex];
            bounded_vector<double, 2> xsiglobal = n->getXsisGlobal();
            bounded_vector<double, 2> blendF = blendFunction(n->getDistanceToBoundary());

            std::vector<ControlPoint *> connection = cell->getControlPoints();
            vector<double> wpc2(connection.size());
            bounded_vector<int, 2> INC_;
            INC_ = connection[connection.size() - 1]->getINC();
            for (int i = 0; i < connection.size(); i++)
            {
                wpc2(i) = connection[i]->getWeight();
            }

            vector<double> phi = cell->shapeFunction(xsiglobal, wpc2, INC_);

            bounded_vector<double, 2> dispInter;
            dispInter = (1.0 - blendF(0)) * dispNode;

            for (int i = 0; i < connection.size(); i++)
            {
                bounded_vector<double, 2> dispCP = connection[i]->getCurrentDisplacement();
                dispInter += blendF(0) * phi(i) * dispCP;
            }

            file << dispInter(0) << " " << dispInter(1) << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";
    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"TheoreticalDisplacement\" format=\"ascii\">"
         << "\n";
    for (Node *n : nodes_)
    {
        // double a = 0.25;
        // double S = 1000.0;
        // double young, poisson, density;
        // materials_[0]->setProperties(young, poisson, density);
        // double k = (3.0 - poisson) / (1.0 + poisson);
        // double mi = young / (2.0 * (1.0 + poisson));

        // bounded_vector<double, 2> coord = n->getInitialCoordinate();
        // double theta = atan2(coord(1), coord(0));
        // double radius = norm_2(coord);
        // bounded_vector<double, 2> teoDisp, dif;
        // teoDisp(0) = S * a / (8.0 * mi) * (radius / a * (k + 1.0) * cos(theta) + 2.0 * a / radius * ((1.0 + k) * cos(theta) + cos(3.0 * theta)) - 2.0 * a * a * a / (radius * radius * radius) * cos(3.0 * theta));
        // teoDisp(1) = S * a / (8.0 * mi) * (radius / a * (k - 3.0) * sin(theta) + 2.0 * a / radius * ((1.0 - k) * sin(theta) + sin(3.0 * theta)) - 2.0 * a * a * a / (radius * radius * radius) * sin(3.0 * theta));

        // file << teoDisp(0) << " " << teoDisp(1) << std::endl;

        ///ANEL
        double sigmaA = -8.0, sigmaB = 6.0;
        double young, poisson, density;
        materials_[0]->setProperties(young, poisson, density);
        double K1, K2;
        if (planeState_ == "EPD")
        {
            K1 = (1.0 - poisson * poisson) / young;
            K2 = poisson / (1.0 - poisson);
        }
        else
        {
        }
        double a = 2.0;
        double b = 6.0;
        double C = (a * a * b * b * (sigmaA - sigmaB)) / (b * b - a * a);
        double b2 = (sigmaB * b * b - sigmaA * a * a) / (b * b - a * a);

        bounded_vector<double, 2> coord = n->getInitialCoordinate();
        double tetaa = atan2(coord(1), coord(0));
        double radius = norm_2(coord);
        double urr = K1 * (b2 * (1.0 - K2) * radius - C / radius * (1.0 + K2));
        file << cos(tetaa) * urr << " " << sin(tetaa) * urr << std::endl;

        //BARRINHA FORÇA COSSENO
        // double young, poisson, density;
        // materials_[0]->setProperties(young, poisson, density);
        // double l = 10.0;
        // double S = 1.0;
        // double pi = 3.1415926535897932384626433;

        // double ux = -1.0 / (young * S) * (-l * l / (pi * pi) * cos(pi * coord(0) / l) - 2.0 * l * coord(0) / (pi * pi) + l * l / (pi * pi));

        // file << ux << " " << 0.0 << std::endl;
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

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"blendingZone\" format=\"ascii\">"
         << "\n";
    for (Node *n : nodes_)
    {
        int cellIndex = n->getCellIndex();
        if (cellIndex == -1)
        {
            file << 0 << "\n";
        }
        else
        {
            file << 1 << "\n";
        }
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

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"Index\" format=\"ascii\">" << std::endl;

    for (Element *el : elements_)
    {
        file << el->getIndex() << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"J-Integral\" format=\"ascii\">" << std::endl;
    if (elementsSideJintegral_.size() >= 1)
    {
        for (Element *el : elements_)
        {
            file << elementsSideJintegral_[el->getIndex()] << "\n";
        }
    }
    else
    {
        for (Element *el : elements_)
        {
            file << -1 << "\n";
        }
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

void GlobalSolid::exportToParaviewAuxiliarMesh(const std::string &name)
{
    std::vector<std::string> lines;
    lines.push_back("l0");
    lines.push_back("l1");
    lines.push_back("l2");

    std::string gmshCode;

    int point = 0;
    int numberbound = 0;
    for (std::string line : lines)
    {
        std::vector<BoundaryElement *> boundL = finiteElementBoundary_[line];
        int aux = 0;
        for (int i = 0; i < boundL.size(); i++)
        {
            if (aux == 0)
            {
                bounded_vector<double, 2> coord = boundL[i]->getNode(0)->getCurrentCoordinate();
                point++;
                gmshCode += "Point(" + std::to_string(point) + ") = {" + std::to_string(coord(0)) + ", " + std::to_string(coord(1)) + ", 0.05, 1.0};\n//\n";

                point++;
                gmshCode += "Point(" + std::to_string(point) + ") = {" + std::to_string(coord(0)) + ", " + std::to_string(coord(1)) + ", 0.05, 1.0};\n//\n";

                coord = boundL[i]->getNode(1)->getCurrentCoordinate();
                point++;
                gmshCode += "Point(" + std::to_string(point) + ") = {" + std::to_string(coord(0)) + ", " + std::to_string(coord(1)) + ", 0.05, 1.0};\n//\n";

                aux++;
            }
            else
            {
                bounded_vector<double, 2> coord = boundL[i]->getNode(1)->getCurrentCoordinate();
                point++;
                gmshCode += "Point(" + std::to_string(point) + ") = {" + std::to_string(coord(0)) + ", " + std::to_string(coord(1)) + ", 0.05, 1.0};\n//\n";

                if (i == boundL.size() - 1)
                {
                    point++;
                    gmshCode += "Point(" + std::to_string(point) + ") = {" + std::to_string(coord(0)) + ", " + std::to_string(coord(1)) + ", 0.05, 1.0};\n//\n";
                }
            }
            numberbound++;
        }
    }

    gmshCode += "BSpline(1) = {";
    for (int i = 1; i <= point; i++)
    {
        gmshCode += std::to_string(i);
        if (i != point)
        {
            gmshCode += ", ";
        }
        else
        {
            gmshCode += ", 1};\n//\n";
        }
    }

    gmshCode += "Line Loop(1) = {1}; \n//\n";
    gmshCode += "Plane Surface(1) = {1};\n//\n";
    gmshCode += "Transfinite Line {1} = " + std::to_string(numberbound) + " Using Progression 1;\n//\n";
    gmshCode += "Mesh.ElementOrder = 1; \n//\n";
    gmshCode += "Mesh 2; \n//\n";
    gmshCode += "Save \"outputAUX" + name + ".vtk\";\n//\n";

    std::ofstream file("temp.geo");
    file << gmshCode;
    file.close();

    std::string cmd = current_working_dir_ + "/Solid/Mesh/gmsh ";
    cmd += "temp.geo -";
    cmd += " -v 0";

    system(cmd.c_str());

    std::string aux = "temp*";
    system((remove + aux).c_str());
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
    else if (orderElemement == 2)
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
    else if (orderElemement == 3)
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
    else if (orderParaview_ == 19)
    {
        coord(0, 0) = -1.0000000000000000;
        coord(1, 0) = 1.0000000000000000;
        coord(2, 0) = 1.0000000000000000;
        coord(3, 0) = -1.0000000000000000;
        coord(4, 0) = -0.8947368421052630;
        coord(5, 0) = -0.7894736842105260;
        coord(6, 0) = -0.6842105263157900;
        coord(7, 0) = -0.5789473684210530;
        coord(8, 0) = -0.4736842105263160;
        coord(9, 0) = -0.3684210526315790;
        coord(10, 0) = -0.2631578947368420;
        coord(11, 0) = -0.1578947368421050;
        coord(12, 0) = -0.0526315789473685;
        coord(13, 0) = 0.0526315789473684;
        coord(14, 0) = 0.1578947368421050;
        coord(15, 0) = 0.2631578947368420;
        coord(16, 0) = 0.3684210526315790;
        coord(17, 0) = 0.4736842105263160;
        coord(18, 0) = 0.5789473684210530;
        coord(19, 0) = 0.6842105263157890;
        coord(20, 0) = 0.7894736842105260;
        coord(21, 0) = 0.8947368421052630;
        coord(22, 0) = 1.0000000000000000;
        coord(23, 0) = 1.0000000000000000;
        coord(24, 0) = 1.0000000000000000;
        coord(25, 0) = 1.0000000000000000;
        coord(26, 0) = 1.0000000000000000;
        coord(27, 0) = 1.0000000000000000;
        coord(28, 0) = 1.0000000000000000;
        coord(29, 0) = 1.0000000000000000;
        coord(30, 0) = 1.0000000000000000;
        coord(31, 0) = 1.0000000000000000;
        coord(32, 0) = 1.0000000000000000;
        coord(33, 0) = 1.0000000000000000;
        coord(34, 0) = 1.0000000000000000;
        coord(35, 0) = 1.0000000000000000;
        coord(36, 0) = 1.0000000000000000;
        coord(37, 0) = 1.0000000000000000;
        coord(38, 0) = 1.0000000000000000;
        coord(39, 0) = 1.0000000000000000;
        coord(40, 0) = -0.8947368421052630;
        coord(41, 0) = -0.7894736842105260;
        coord(42, 0) = -0.6842105263157900;
        coord(43, 0) = -0.5789473684210530;
        coord(44, 0) = -0.4736842105263160;
        coord(45, 0) = -0.3684210526315790;
        coord(46, 0) = -0.2631578947368420;
        coord(47, 0) = -0.1578947368421050;
        coord(48, 0) = -0.0526315789473685;
        coord(49, 0) = 0.0526315789473684;
        coord(50, 0) = 0.1578947368421050;
        coord(51, 0) = 0.2631578947368420;
        coord(52, 0) = 0.3684210526315790;
        coord(53, 0) = 0.4736842105263160;
        coord(54, 0) = 0.5789473684210530;
        coord(55, 0) = 0.6842105263157890;
        coord(56, 0) = 0.7894736842105260;
        coord(57, 0) = 0.8947368421052630;
        coord(58, 0) = -1.0000000000000000;
        coord(59, 0) = -1.0000000000000000;
        coord(60, 0) = -1.0000000000000000;
        coord(61, 0) = -1.0000000000000000;
        coord(62, 0) = -1.0000000000000000;
        coord(63, 0) = -1.0000000000000000;
        coord(64, 0) = -1.0000000000000000;
        coord(65, 0) = -1.0000000000000000;
        coord(66, 0) = -1.0000000000000000;
        coord(67, 0) = -1.0000000000000000;
        coord(68, 0) = -1.0000000000000000;
        coord(69, 0) = -1.0000000000000000;
        coord(70, 0) = -1.0000000000000000;
        coord(71, 0) = -1.0000000000000000;
        coord(72, 0) = -1.0000000000000000;
        coord(73, 0) = -1.0000000000000000;
        coord(74, 0) = -1.0000000000000000;
        coord(75, 0) = -1.0000000000000000;
        coord(76, 0) = -0.8947368421052630;
        coord(77, 0) = -0.7894736842105260;
        coord(78, 0) = -0.6842105263157900;
        coord(79, 0) = -0.5789473684210530;
        coord(80, 0) = -0.4736842105263160;
        coord(81, 0) = -0.3684210526315790;
        coord(82, 0) = -0.2631578947368420;
        coord(83, 0) = -0.1578947368421050;
        coord(84, 0) = -0.0526315789473685;
        coord(85, 0) = 0.0526315789473684;
        coord(86, 0) = 0.1578947368421050;
        coord(87, 0) = 0.2631578947368420;
        coord(88, 0) = 0.3684210526315790;
        coord(89, 0) = 0.4736842105263160;
        coord(90, 0) = 0.5789473684210530;
        coord(91, 0) = 0.6842105263157890;
        coord(92, 0) = 0.7894736842105260;
        coord(93, 0) = 0.8947368421052630;
        coord(94, 0) = -0.8947368421052630;
        coord(95, 0) = -0.7894736842105260;
        coord(96, 0) = -0.6842105263157900;
        coord(97, 0) = -0.5789473684210530;
        coord(98, 0) = -0.4736842105263160;
        coord(99, 0) = -0.3684210526315790;
        coord(100, 0) = -0.2631578947368420;
        coord(101, 0) = -0.1578947368421050;
        coord(102, 0) = -0.0526315789473685;
        coord(103, 0) = 0.0526315789473684;
        coord(104, 0) = 0.1578947368421050;
        coord(105, 0) = 0.2631578947368420;
        coord(106, 0) = 0.3684210526315790;
        coord(107, 0) = 0.4736842105263160;
        coord(108, 0) = 0.5789473684210530;
        coord(109, 0) = 0.6842105263157890;
        coord(110, 0) = 0.7894736842105260;
        coord(111, 0) = 0.8947368421052630;
        coord(112, 0) = -0.8947368421052630;
        coord(113, 0) = -0.7894736842105260;
        coord(114, 0) = -0.6842105263157900;
        coord(115, 0) = -0.5789473684210530;
        coord(116, 0) = -0.4736842105263160;
        coord(117, 0) = -0.3684210526315790;
        coord(118, 0) = -0.2631578947368420;
        coord(119, 0) = -0.1578947368421050;
        coord(120, 0) = -0.0526315789473685;
        coord(121, 0) = 0.0526315789473684;
        coord(122, 0) = 0.1578947368421050;
        coord(123, 0) = 0.2631578947368420;
        coord(124, 0) = 0.3684210526315790;
        coord(125, 0) = 0.4736842105263160;
        coord(126, 0) = 0.5789473684210530;
        coord(127, 0) = 0.6842105263157890;
        coord(128, 0) = 0.7894736842105260;
        coord(129, 0) = 0.8947368421052630;
        coord(130, 0) = -0.8947368421052630;
        coord(131, 0) = -0.7894736842105260;
        coord(132, 0) = -0.6842105263157900;
        coord(133, 0) = -0.5789473684210530;
        coord(134, 0) = -0.4736842105263160;
        coord(135, 0) = -0.3684210526315790;
        coord(136, 0) = -0.2631578947368420;
        coord(137, 0) = -0.1578947368421050;
        coord(138, 0) = -0.0526315789473685;
        coord(139, 0) = 0.0526315789473684;
        coord(140, 0) = 0.1578947368421050;
        coord(141, 0) = 0.2631578947368420;
        coord(142, 0) = 0.3684210526315790;
        coord(143, 0) = 0.4736842105263160;
        coord(144, 0) = 0.5789473684210530;
        coord(145, 0) = 0.6842105263157890;
        coord(146, 0) = 0.7894736842105260;
        coord(147, 0) = 0.8947368421052630;
        coord(148, 0) = -0.8947368421052630;
        coord(149, 0) = -0.7894736842105260;
        coord(150, 0) = -0.6842105263157900;
        coord(151, 0) = -0.5789473684210530;
        coord(152, 0) = -0.4736842105263160;
        coord(153, 0) = -0.3684210526315790;
        coord(154, 0) = -0.2631578947368420;
        coord(155, 0) = -0.1578947368421050;
        coord(156, 0) = -0.0526315789473685;
        coord(157, 0) = 0.0526315789473684;
        coord(158, 0) = 0.1578947368421050;
        coord(159, 0) = 0.2631578947368420;
        coord(160, 0) = 0.3684210526315790;
        coord(161, 0) = 0.4736842105263160;
        coord(162, 0) = 0.5789473684210530;
        coord(163, 0) = 0.6842105263157890;
        coord(164, 0) = 0.7894736842105260;
        coord(165, 0) = 0.8947368421052630;
        coord(166, 0) = -0.8947368421052630;
        coord(167, 0) = -0.7894736842105260;
        coord(168, 0) = -0.6842105263157900;
        coord(169, 0) = -0.5789473684210530;
        coord(170, 0) = -0.4736842105263160;
        coord(171, 0) = -0.3684210526315790;
        coord(172, 0) = -0.2631578947368420;
        coord(173, 0) = -0.1578947368421050;
        coord(174, 0) = -0.0526315789473685;
        coord(175, 0) = 0.0526315789473684;
        coord(176, 0) = 0.1578947368421050;
        coord(177, 0) = 0.2631578947368420;
        coord(178, 0) = 0.3684210526315790;
        coord(179, 0) = 0.4736842105263160;
        coord(180, 0) = 0.5789473684210530;
        coord(181, 0) = 0.6842105263157890;
        coord(182, 0) = 0.7894736842105260;
        coord(183, 0) = 0.8947368421052630;
        coord(184, 0) = -0.8947368421052630;
        coord(185, 0) = -0.7894736842105260;
        coord(186, 0) = -0.6842105263157900;
        coord(187, 0) = -0.5789473684210530;
        coord(188, 0) = -0.4736842105263160;
        coord(189, 0) = -0.3684210526315790;
        coord(190, 0) = -0.2631578947368420;
        coord(191, 0) = -0.1578947368421050;
        coord(192, 0) = -0.0526315789473685;
        coord(193, 0) = 0.0526315789473684;
        coord(194, 0) = 0.1578947368421050;
        coord(195, 0) = 0.2631578947368420;
        coord(196, 0) = 0.3684210526315790;
        coord(197, 0) = 0.4736842105263160;
        coord(198, 0) = 0.5789473684210530;
        coord(199, 0) = 0.6842105263157890;
        coord(200, 0) = 0.7894736842105260;
        coord(201, 0) = 0.8947368421052630;
        coord(202, 0) = -0.8947368421052630;
        coord(203, 0) = -0.7894736842105260;
        coord(204, 0) = -0.6842105263157900;
        coord(205, 0) = -0.5789473684210530;
        coord(206, 0) = -0.4736842105263160;
        coord(207, 0) = -0.3684210526315790;
        coord(208, 0) = -0.2631578947368420;
        coord(209, 0) = -0.1578947368421050;
        coord(210, 0) = -0.0526315789473685;
        coord(211, 0) = 0.0526315789473684;
        coord(212, 0) = 0.1578947368421050;
        coord(213, 0) = 0.2631578947368420;
        coord(214, 0) = 0.3684210526315790;
        coord(215, 0) = 0.4736842105263160;
        coord(216, 0) = 0.5789473684210530;
        coord(217, 0) = 0.6842105263157890;
        coord(218, 0) = 0.7894736842105260;
        coord(219, 0) = 0.8947368421052630;
        coord(220, 0) = -0.8947368421052630;
        coord(221, 0) = -0.7894736842105260;
        coord(222, 0) = -0.6842105263157900;
        coord(223, 0) = -0.5789473684210530;
        coord(224, 0) = -0.4736842105263160;
        coord(225, 0) = -0.3684210526315790;
        coord(226, 0) = -0.2631578947368420;
        coord(227, 0) = -0.1578947368421050;
        coord(228, 0) = -0.0526315789473685;
        coord(229, 0) = 0.0526315789473684;
        coord(230, 0) = 0.1578947368421050;
        coord(231, 0) = 0.2631578947368420;
        coord(232, 0) = 0.3684210526315790;
        coord(233, 0) = 0.4736842105263160;
        coord(234, 0) = 0.5789473684210530;
        coord(235, 0) = 0.6842105263157890;
        coord(236, 0) = 0.7894736842105260;
        coord(237, 0) = 0.8947368421052630;
        coord(238, 0) = -0.8947368421052630;
        coord(239, 0) = -0.7894736842105260;
        coord(240, 0) = -0.6842105263157900;
        coord(241, 0) = -0.5789473684210530;
        coord(242, 0) = -0.4736842105263160;
        coord(243, 0) = -0.3684210526315790;
        coord(244, 0) = -0.2631578947368420;
        coord(245, 0) = -0.1578947368421050;
        coord(246, 0) = -0.0526315789473685;
        coord(247, 0) = 0.0526315789473684;
        coord(248, 0) = 0.1578947368421050;
        coord(249, 0) = 0.2631578947368420;
        coord(250, 0) = 0.3684210526315790;
        coord(251, 0) = 0.4736842105263160;
        coord(252, 0) = 0.5789473684210530;
        coord(253, 0) = 0.6842105263157890;
        coord(254, 0) = 0.7894736842105260;
        coord(255, 0) = 0.8947368421052630;
        coord(256, 0) = -0.8947368421052630;
        coord(257, 0) = -0.7894736842105260;
        coord(258, 0) = -0.6842105263157900;
        coord(259, 0) = -0.5789473684210530;
        coord(260, 0) = -0.4736842105263160;
        coord(261, 0) = -0.3684210526315790;
        coord(262, 0) = -0.2631578947368420;
        coord(263, 0) = -0.1578947368421050;
        coord(264, 0) = -0.0526315789473685;
        coord(265, 0) = 0.0526315789473684;
        coord(266, 0) = 0.1578947368421050;
        coord(267, 0) = 0.2631578947368420;
        coord(268, 0) = 0.3684210526315790;
        coord(269, 0) = 0.4736842105263160;
        coord(270, 0) = 0.5789473684210530;
        coord(271, 0) = 0.6842105263157890;
        coord(272, 0) = 0.7894736842105260;
        coord(273, 0) = 0.8947368421052630;
        coord(274, 0) = -0.8947368421052630;
        coord(275, 0) = -0.7894736842105260;
        coord(276, 0) = -0.6842105263157900;
        coord(277, 0) = -0.5789473684210530;
        coord(278, 0) = -0.4736842105263160;
        coord(279, 0) = -0.3684210526315790;
        coord(280, 0) = -0.2631578947368420;
        coord(281, 0) = -0.1578947368421050;
        coord(282, 0) = -0.0526315789473685;
        coord(283, 0) = 0.0526315789473684;
        coord(284, 0) = 0.1578947368421050;
        coord(285, 0) = 0.2631578947368420;
        coord(286, 0) = 0.3684210526315790;
        coord(287, 0) = 0.4736842105263160;
        coord(288, 0) = 0.5789473684210530;
        coord(289, 0) = 0.6842105263157890;
        coord(290, 0) = 0.7894736842105260;
        coord(291, 0) = 0.8947368421052630;
        coord(292, 0) = -0.8947368421052630;
        coord(293, 0) = -0.7894736842105260;
        coord(294, 0) = -0.6842105263157900;
        coord(295, 0) = -0.5789473684210530;
        coord(296, 0) = -0.4736842105263160;
        coord(297, 0) = -0.3684210526315790;
        coord(298, 0) = -0.2631578947368420;
        coord(299, 0) = -0.1578947368421050;
        coord(300, 0) = -0.0526315789473685;
        coord(301, 0) = 0.0526315789473684;
        coord(302, 0) = 0.1578947368421050;
        coord(303, 0) = 0.2631578947368420;
        coord(304, 0) = 0.3684210526315790;
        coord(305, 0) = 0.4736842105263160;
        coord(306, 0) = 0.5789473684210530;
        coord(307, 0) = 0.6842105263157890;
        coord(308, 0) = 0.7894736842105260;
        coord(309, 0) = 0.8947368421052630;
        coord(310, 0) = -0.8947368421052630;
        coord(311, 0) = -0.7894736842105260;
        coord(312, 0) = -0.6842105263157900;
        coord(313, 0) = -0.5789473684210530;
        coord(314, 0) = -0.4736842105263160;
        coord(315, 0) = -0.3684210526315790;
        coord(316, 0) = -0.2631578947368420;
        coord(317, 0) = -0.1578947368421050;
        coord(318, 0) = -0.0526315789473685;
        coord(319, 0) = 0.0526315789473684;
        coord(320, 0) = 0.1578947368421050;
        coord(321, 0) = 0.2631578947368420;
        coord(322, 0) = 0.3684210526315790;
        coord(323, 0) = 0.4736842105263160;
        coord(324, 0) = 0.5789473684210530;
        coord(325, 0) = 0.6842105263157890;
        coord(326, 0) = 0.7894736842105260;
        coord(327, 0) = 0.8947368421052630;
        coord(328, 0) = -0.8947368421052630;
        coord(329, 0) = -0.7894736842105260;
        coord(330, 0) = -0.6842105263157900;
        coord(331, 0) = -0.5789473684210530;
        coord(332, 0) = -0.4736842105263160;
        coord(333, 0) = -0.3684210526315790;
        coord(334, 0) = -0.2631578947368420;
        coord(335, 0) = -0.1578947368421050;
        coord(336, 0) = -0.0526315789473685;
        coord(337, 0) = 0.0526315789473684;
        coord(338, 0) = 0.1578947368421050;
        coord(339, 0) = 0.2631578947368420;
        coord(340, 0) = 0.3684210526315790;
        coord(341, 0) = 0.4736842105263160;
        coord(342, 0) = 0.5789473684210530;
        coord(343, 0) = 0.6842105263157890;
        coord(344, 0) = 0.7894736842105260;
        coord(345, 0) = 0.8947368421052630;
        coord(346, 0) = -0.8947368421052630;
        coord(347, 0) = -0.7894736842105260;
        coord(348, 0) = -0.6842105263157900;
        coord(349, 0) = -0.5789473684210530;
        coord(350, 0) = -0.4736842105263160;
        coord(351, 0) = -0.3684210526315790;
        coord(352, 0) = -0.2631578947368420;
        coord(353, 0) = -0.1578947368421050;
        coord(354, 0) = -0.0526315789473685;
        coord(355, 0) = 0.0526315789473684;
        coord(356, 0) = 0.1578947368421050;
        coord(357, 0) = 0.2631578947368420;
        coord(358, 0) = 0.3684210526315790;
        coord(359, 0) = 0.4736842105263160;
        coord(360, 0) = 0.5789473684210530;
        coord(361, 0) = 0.6842105263157890;
        coord(362, 0) = 0.7894736842105260;
        coord(363, 0) = 0.8947368421052630;
        coord(364, 0) = -0.8947368421052630;
        coord(365, 0) = -0.7894736842105260;
        coord(366, 0) = -0.6842105263157900;
        coord(367, 0) = -0.5789473684210530;
        coord(368, 0) = -0.4736842105263160;
        coord(369, 0) = -0.3684210526315790;
        coord(370, 0) = -0.2631578947368420;
        coord(371, 0) = -0.1578947368421050;
        coord(372, 0) = -0.0526315789473685;
        coord(373, 0) = 0.0526315789473684;
        coord(374, 0) = 0.1578947368421050;
        coord(375, 0) = 0.2631578947368420;
        coord(376, 0) = 0.3684210526315790;
        coord(377, 0) = 0.4736842105263160;
        coord(378, 0) = 0.5789473684210530;
        coord(379, 0) = 0.6842105263157890;
        coord(380, 0) = 0.7894736842105260;
        coord(381, 0) = 0.8947368421052630;
        coord(382, 0) = -0.8947368421052630;
        coord(383, 0) = -0.7894736842105260;
        coord(384, 0) = -0.6842105263157900;
        coord(385, 0) = -0.5789473684210530;
        coord(386, 0) = -0.4736842105263160;
        coord(387, 0) = -0.3684210526315790;
        coord(388, 0) = -0.2631578947368420;
        coord(389, 0) = -0.1578947368421050;
        coord(390, 0) = -0.0526315789473685;
        coord(391, 0) = 0.0526315789473684;
        coord(392, 0) = 0.1578947368421050;
        coord(393, 0) = 0.2631578947368420;
        coord(394, 0) = 0.3684210526315790;
        coord(395, 0) = 0.4736842105263160;
        coord(396, 0) = 0.5789473684210530;
        coord(397, 0) = 0.6842105263157890;
        coord(398, 0) = 0.7894736842105260;
        coord(399, 0) = 0.8947368421052630;
        coord(0, 1) = -1.0000000000000000;
        coord(1, 1) = -1.0000000000000000;
        coord(2, 1) = 1.0000000000000000;
        coord(3, 1) = 1.0000000000000000;
        coord(4, 1) = -1.0000000000000000;
        coord(5, 1) = -1.0000000000000000;
        coord(6, 1) = -1.0000000000000000;
        coord(7, 1) = -1.0000000000000000;
        coord(8, 1) = -1.0000000000000000;
        coord(9, 1) = -1.0000000000000000;
        coord(10, 1) = -1.0000000000000000;
        coord(11, 1) = -1.0000000000000000;
        coord(12, 1) = -1.0000000000000000;
        coord(13, 1) = -1.0000000000000000;
        coord(14, 1) = -1.0000000000000000;
        coord(15, 1) = -1.0000000000000000;
        coord(16, 1) = -1.0000000000000000;
        coord(17, 1) = -1.0000000000000000;
        coord(18, 1) = -1.0000000000000000;
        coord(19, 1) = -1.0000000000000000;
        coord(20, 1) = -1.0000000000000000;
        coord(21, 1) = -1.0000000000000000;
        coord(22, 1) = -0.8947368421052630;
        coord(23, 1) = -0.7894736842105260;
        coord(24, 1) = -0.6842105263157900;
        coord(25, 1) = -0.5789473684210530;
        coord(26, 1) = -0.4736842105263160;
        coord(27, 1) = -0.3684210526315790;
        coord(28, 1) = -0.2631578947368420;
        coord(29, 1) = -0.1578947368421050;
        coord(30, 1) = -0.0526315789473685;
        coord(31, 1) = 0.0526315789473684;
        coord(32, 1) = 0.1578947368421050;
        coord(33, 1) = 0.2631578947368420;
        coord(34, 1) = 0.3684210526315790;
        coord(35, 1) = 0.4736842105263160;
        coord(36, 1) = 0.5789473684210530;
        coord(37, 1) = 0.6842105263157890;
        coord(38, 1) = 0.7894736842105260;
        coord(39, 1) = 0.8947368421052630;
        coord(40, 1) = 1.0000000000000000;
        coord(41, 1) = 1.0000000000000000;
        coord(42, 1) = 1.0000000000000000;
        coord(43, 1) = 1.0000000000000000;
        coord(44, 1) = 1.0000000000000000;
        coord(45, 1) = 1.0000000000000000;
        coord(46, 1) = 1.0000000000000000;
        coord(47, 1) = 1.0000000000000000;
        coord(48, 1) = 1.0000000000000000;
        coord(49, 1) = 1.0000000000000000;
        coord(50, 1) = 1.0000000000000000;
        coord(51, 1) = 1.0000000000000000;
        coord(52, 1) = 1.0000000000000000;
        coord(53, 1) = 1.0000000000000000;
        coord(54, 1) = 1.0000000000000000;
        coord(55, 1) = 1.0000000000000000;
        coord(56, 1) = 1.0000000000000000;
        coord(57, 1) = 1.0000000000000000;
        coord(58, 1) = -0.8947368421052630;
        coord(59, 1) = -0.7894736842105260;
        coord(60, 1) = -0.6842105263157900;
        coord(61, 1) = -0.5789473684210530;
        coord(62, 1) = -0.4736842105263160;
        coord(63, 1) = -0.3684210526315790;
        coord(64, 1) = -0.2631578947368420;
        coord(65, 1) = -0.1578947368421050;
        coord(66, 1) = -0.0526315789473685;
        coord(67, 1) = 0.0526315789473684;
        coord(68, 1) = 0.1578947368421050;
        coord(69, 1) = 0.2631578947368420;
        coord(70, 1) = 0.3684210526315790;
        coord(71, 1) = 0.4736842105263160;
        coord(72, 1) = 0.5789473684210530;
        coord(73, 1) = 0.6842105263157890;
        coord(74, 1) = 0.7894736842105260;
        coord(75, 1) = 0.8947368421052630;
        coord(76, 1) = -0.8947368421052630;
        coord(77, 1) = -0.8947368421052630;
        coord(78, 1) = -0.8947368421052630;
        coord(79, 1) = -0.8947368421052630;
        coord(80, 1) = -0.8947368421052630;
        coord(81, 1) = -0.8947368421052630;
        coord(82, 1) = -0.8947368421052630;
        coord(83, 1) = -0.8947368421052630;
        coord(84, 1) = -0.8947368421052630;
        coord(85, 1) = -0.8947368421052630;
        coord(86, 1) = -0.8947368421052630;
        coord(87, 1) = -0.8947368421052630;
        coord(88, 1) = -0.8947368421052630;
        coord(89, 1) = -0.8947368421052630;
        coord(90, 1) = -0.8947368421052630;
        coord(91, 1) = -0.8947368421052630;
        coord(92, 1) = -0.8947368421052630;
        coord(93, 1) = -0.8947368421052630;
        coord(94, 1) = -0.7894736842105260;
        coord(95, 1) = -0.7894736842105260;
        coord(96, 1) = -0.7894736842105260;
        coord(97, 1) = -0.7894736842105260;
        coord(98, 1) = -0.7894736842105260;
        coord(99, 1) = -0.7894736842105260;
        coord(100, 1) = -0.7894736842105260;
        coord(101, 1) = -0.7894736842105260;
        coord(102, 1) = -0.7894736842105260;
        coord(103, 1) = -0.7894736842105260;
        coord(104, 1) = -0.7894736842105260;
        coord(105, 1) = -0.7894736842105260;
        coord(106, 1) = -0.7894736842105260;
        coord(107, 1) = -0.7894736842105260;
        coord(108, 1) = -0.7894736842105260;
        coord(109, 1) = -0.7894736842105260;
        coord(110, 1) = -0.7894736842105260;
        coord(111, 1) = -0.7894736842105260;
        coord(112, 1) = -0.6842105263157900;
        coord(113, 1) = -0.6842105263157900;
        coord(114, 1) = -0.6842105263157900;
        coord(115, 1) = -0.6842105263157900;
        coord(116, 1) = -0.6842105263157900;
        coord(117, 1) = -0.6842105263157900;
        coord(118, 1) = -0.6842105263157900;
        coord(119, 1) = -0.6842105263157900;
        coord(120, 1) = -0.6842105263157900;
        coord(121, 1) = -0.6842105263157900;
        coord(122, 1) = -0.6842105263157900;
        coord(123, 1) = -0.6842105263157900;
        coord(124, 1) = -0.6842105263157900;
        coord(125, 1) = -0.6842105263157900;
        coord(126, 1) = -0.6842105263157900;
        coord(127, 1) = -0.6842105263157900;
        coord(128, 1) = -0.6842105263157900;
        coord(129, 1) = -0.6842105263157900;
        coord(130, 1) = -0.5789473684210530;
        coord(131, 1) = -0.5789473684210530;
        coord(132, 1) = -0.5789473684210530;
        coord(133, 1) = -0.5789473684210530;
        coord(134, 1) = -0.5789473684210530;
        coord(135, 1) = -0.5789473684210530;
        coord(136, 1) = -0.5789473684210530;
        coord(137, 1) = -0.5789473684210530;
        coord(138, 1) = -0.5789473684210530;
        coord(139, 1) = -0.5789473684210530;
        coord(140, 1) = -0.5789473684210530;
        coord(141, 1) = -0.5789473684210530;
        coord(142, 1) = -0.5789473684210530;
        coord(143, 1) = -0.5789473684210530;
        coord(144, 1) = -0.5789473684210530;
        coord(145, 1) = -0.5789473684210530;
        coord(146, 1) = -0.5789473684210530;
        coord(147, 1) = -0.5789473684210530;
        coord(148, 1) = -0.4736842105263160;
        coord(149, 1) = -0.4736842105263160;
        coord(150, 1) = -0.4736842105263160;
        coord(151, 1) = -0.4736842105263160;
        coord(152, 1) = -0.4736842105263160;
        coord(153, 1) = -0.4736842105263160;
        coord(154, 1) = -0.4736842105263160;
        coord(155, 1) = -0.4736842105263160;
        coord(156, 1) = -0.4736842105263160;
        coord(157, 1) = -0.4736842105263160;
        coord(158, 1) = -0.4736842105263160;
        coord(159, 1) = -0.4736842105263160;
        coord(160, 1) = -0.4736842105263160;
        coord(161, 1) = -0.4736842105263160;
        coord(162, 1) = -0.4736842105263160;
        coord(163, 1) = -0.4736842105263160;
        coord(164, 1) = -0.4736842105263160;
        coord(165, 1) = -0.4736842105263160;
        coord(166, 1) = -0.3684210526315790;
        coord(167, 1) = -0.3684210526315790;
        coord(168, 1) = -0.3684210526315790;
        coord(169, 1) = -0.3684210526315790;
        coord(170, 1) = -0.3684210526315790;
        coord(171, 1) = -0.3684210526315790;
        coord(172, 1) = -0.3684210526315790;
        coord(173, 1) = -0.3684210526315790;
        coord(174, 1) = -0.3684210526315790;
        coord(175, 1) = -0.3684210526315790;
        coord(176, 1) = -0.3684210526315790;
        coord(177, 1) = -0.3684210526315790;
        coord(178, 1) = -0.3684210526315790;
        coord(179, 1) = -0.3684210526315790;
        coord(180, 1) = -0.3684210526315790;
        coord(181, 1) = -0.3684210526315790;
        coord(182, 1) = -0.3684210526315790;
        coord(183, 1) = -0.3684210526315790;
        coord(184, 1) = -0.2631578947368420;
        coord(185, 1) = -0.2631578947368420;
        coord(186, 1) = -0.2631578947368420;
        coord(187, 1) = -0.2631578947368420;
        coord(188, 1) = -0.2631578947368420;
        coord(189, 1) = -0.2631578947368420;
        coord(190, 1) = -0.2631578947368420;
        coord(191, 1) = -0.2631578947368420;
        coord(192, 1) = -0.2631578947368420;
        coord(193, 1) = -0.2631578947368420;
        coord(194, 1) = -0.2631578947368420;
        coord(195, 1) = -0.2631578947368420;
        coord(196, 1) = -0.2631578947368420;
        coord(197, 1) = -0.2631578947368420;
        coord(198, 1) = -0.2631578947368420;
        coord(199, 1) = -0.2631578947368420;
        coord(200, 1) = -0.2631578947368420;
        coord(201, 1) = -0.2631578947368420;
        coord(202, 1) = -0.1578947368421050;
        coord(203, 1) = -0.1578947368421050;
        coord(204, 1) = -0.1578947368421050;
        coord(205, 1) = -0.1578947368421050;
        coord(206, 1) = -0.1578947368421050;
        coord(207, 1) = -0.1578947368421050;
        coord(208, 1) = -0.1578947368421050;
        coord(209, 1) = -0.1578947368421050;
        coord(210, 1) = -0.1578947368421050;
        coord(211, 1) = -0.1578947368421050;
        coord(212, 1) = -0.1578947368421050;
        coord(213, 1) = -0.1578947368421050;
        coord(214, 1) = -0.1578947368421050;
        coord(215, 1) = -0.1578947368421050;
        coord(216, 1) = -0.1578947368421050;
        coord(217, 1) = -0.1578947368421050;
        coord(218, 1) = -0.1578947368421050;
        coord(219, 1) = -0.1578947368421050;
        coord(220, 1) = -0.0526315789473685;
        coord(221, 1) = -0.0526315789473685;
        coord(222, 1) = -0.0526315789473685;
        coord(223, 1) = -0.0526315789473685;
        coord(224, 1) = -0.0526315789473685;
        coord(225, 1) = -0.0526315789473685;
        coord(226, 1) = -0.0526315789473685;
        coord(227, 1) = -0.0526315789473685;
        coord(228, 1) = -0.0526315789473685;
        coord(229, 1) = -0.0526315789473685;
        coord(230, 1) = -0.0526315789473685;
        coord(231, 1) = -0.0526315789473685;
        coord(232, 1) = -0.0526315789473685;
        coord(233, 1) = -0.0526315789473685;
        coord(234, 1) = -0.0526315789473685;
        coord(235, 1) = -0.0526315789473685;
        coord(236, 1) = -0.0526315789473685;
        coord(237, 1) = -0.0526315789473685;
        coord(238, 1) = 0.0526315789473684;
        coord(239, 1) = 0.0526315789473684;
        coord(240, 1) = 0.0526315789473684;
        coord(241, 1) = 0.0526315789473684;
        coord(242, 1) = 0.0526315789473684;
        coord(243, 1) = 0.0526315789473684;
        coord(244, 1) = 0.0526315789473684;
        coord(245, 1) = 0.0526315789473684;
        coord(246, 1) = 0.0526315789473684;
        coord(247, 1) = 0.0526315789473684;
        coord(248, 1) = 0.0526315789473684;
        coord(249, 1) = 0.0526315789473684;
        coord(250, 1) = 0.0526315789473684;
        coord(251, 1) = 0.0526315789473684;
        coord(252, 1) = 0.0526315789473684;
        coord(253, 1) = 0.0526315789473684;
        coord(254, 1) = 0.0526315789473684;
        coord(255, 1) = 0.0526315789473684;
        coord(256, 1) = 0.1578947368421050;
        coord(257, 1) = 0.1578947368421050;
        coord(258, 1) = 0.1578947368421050;
        coord(259, 1) = 0.1578947368421060;
        coord(260, 1) = 0.1578947368421060;
        coord(261, 1) = 0.1578947368421060;
        coord(262, 1) = 0.1578947368421060;
        coord(263, 1) = 0.1578947368421060;
        coord(264, 1) = 0.1578947368421070;
        coord(265, 1) = 0.1578947368421070;
        coord(266, 1) = 0.1578947368421070;
        coord(267, 1) = 0.1578947368421070;
        coord(268, 1) = 0.1578947368421070;
        coord(269, 1) = 0.1578947368421080;
        coord(270, 1) = 0.1578947368421080;
        coord(271, 1) = 0.1578947368421080;
        coord(272, 1) = 0.1578947368421080;
        coord(273, 1) = 0.1578947368421080;
        coord(274, 1) = 0.2631578947368420;
        coord(275, 1) = 0.2631578947368420;
        coord(276, 1) = 0.2631578947368420;
        coord(277, 1) = 0.2631578947368420;
        coord(278, 1) = 0.2631578947368420;
        coord(279, 1) = 0.2631578947368420;
        coord(280, 1) = 0.2631578947368420;
        coord(281, 1) = 0.2631578947368420;
        coord(282, 1) = 0.2631578947368420;
        coord(283, 1) = 0.2631578947368420;
        coord(284, 1) = 0.2631578947368420;
        coord(285, 1) = 0.2631578947368420;
        coord(286, 1) = 0.2631578947368420;
        coord(287, 1) = 0.2631578947368420;
        coord(288, 1) = 0.2631578947368420;
        coord(289, 1) = 0.2631578947368420;
        coord(290, 1) = 0.2631578947368420;
        coord(291, 1) = 0.2631578947368420;
        coord(292, 1) = 0.3684210526315790;
        coord(293, 1) = 0.3684210526315790;
        coord(294, 1) = 0.3684210526315790;
        coord(295, 1) = 0.3684210526315790;
        coord(296, 1) = 0.3684210526315790;
        coord(297, 1) = 0.3684210526315790;
        coord(298, 1) = 0.3684210526315790;
        coord(299, 1) = 0.3684210526315790;
        coord(300, 1) = 0.3684210526315790;
        coord(301, 1) = 0.3684210526315790;
        coord(302, 1) = 0.3684210526315790;
        coord(303, 1) = 0.3684210526315790;
        coord(304, 1) = 0.3684210526315790;
        coord(305, 1) = 0.3684210526315790;
        coord(306, 1) = 0.3684210526315790;
        coord(307, 1) = 0.3684210526315790;
        coord(308, 1) = 0.3684210526315790;
        coord(309, 1) = 0.3684210526315790;
        coord(310, 1) = 0.4736842105263160;
        coord(311, 1) = 0.4736842105263160;
        coord(312, 1) = 0.4736842105263160;
        coord(313, 1) = 0.4736842105263160;
        coord(314, 1) = 0.4736842105263160;
        coord(315, 1) = 0.4736842105263160;
        coord(316, 1) = 0.4736842105263160;
        coord(317, 1) = 0.4736842105263160;
        coord(318, 1) = 0.4736842105263160;
        coord(319, 1) = 0.4736842105263160;
        coord(320, 1) = 0.4736842105263160;
        coord(321, 1) = 0.4736842105263160;
        coord(322, 1) = 0.4736842105263160;
        coord(323, 1) = 0.4736842105263160;
        coord(324, 1) = 0.4736842105263160;
        coord(325, 1) = 0.4736842105263160;
        coord(326, 1) = 0.4736842105263160;
        coord(327, 1) = 0.4736842105263160;
        coord(328, 1) = 0.5789473684210530;
        coord(329, 1) = 0.5789473684210530;
        coord(330, 1) = 0.5789473684210530;
        coord(331, 1) = 0.5789473684210530;
        coord(332, 1) = 0.5789473684210530;
        coord(333, 1) = 0.5789473684210530;
        coord(334, 1) = 0.5789473684210530;
        coord(335, 1) = 0.5789473684210530;
        coord(336, 1) = 0.5789473684210530;
        coord(337, 1) = 0.5789473684210530;
        coord(338, 1) = 0.5789473684210530;
        coord(339, 1) = 0.5789473684210530;
        coord(340, 1) = 0.5789473684210530;
        coord(341, 1) = 0.5789473684210530;
        coord(342, 1) = 0.5789473684210530;
        coord(343, 1) = 0.5789473684210530;
        coord(344, 1) = 0.5789473684210530;
        coord(345, 1) = 0.5789473684210530;
        coord(346, 1) = 0.6842105263157890;
        coord(347, 1) = 0.6842105263157890;
        coord(348, 1) = 0.6842105263157890;
        coord(349, 1) = 0.6842105263157890;
        coord(350, 1) = 0.6842105263157890;
        coord(351, 1) = 0.6842105263157890;
        coord(352, 1) = 0.6842105263157890;
        coord(353, 1) = 0.6842105263157890;
        coord(354, 1) = 0.6842105263157890;
        coord(355, 1) = 0.6842105263157890;
        coord(356, 1) = 0.6842105263157890;
        coord(357, 1) = 0.6842105263157890;
        coord(358, 1) = 0.6842105263157890;
        coord(359, 1) = 0.6842105263157890;
        coord(360, 1) = 0.6842105263157890;
        coord(361, 1) = 0.6842105263157890;
        coord(362, 1) = 0.6842105263157890;
        coord(363, 1) = 0.6842105263157890;
        coord(364, 1) = 0.7894736842105260;
        coord(365, 1) = 0.7894736842105260;
        coord(366, 1) = 0.7894736842105260;
        coord(367, 1) = 0.7894736842105260;
        coord(368, 1) = 0.7894736842105260;
        coord(369, 1) = 0.7894736842105260;
        coord(370, 1) = 0.7894736842105260;
        coord(371, 1) = 0.7894736842105260;
        coord(372, 1) = 0.7894736842105260;
        coord(373, 1) = 0.7894736842105260;
        coord(374, 1) = 0.7894736842105260;
        coord(375, 1) = 0.7894736842105260;
        coord(376, 1) = 0.7894736842105260;
        coord(377, 1) = 0.7894736842105260;
        coord(378, 1) = 0.7894736842105260;
        coord(379, 1) = 0.7894736842105260;
        coord(380, 1) = 0.7894736842105260;
        coord(381, 1) = 0.7894736842105260;
        coord(382, 1) = 0.8947368421052630;
        coord(383, 1) = 0.8947368421052630;
        coord(384, 1) = 0.8947368421052630;
        coord(385, 1) = 0.8947368421052630;
        coord(386, 1) = 0.8947368421052630;
        coord(387, 1) = 0.8947368421052630;
        coord(388, 1) = 0.8947368421052630;
        coord(389, 1) = 0.8947368421052630;
        coord(390, 1) = 0.8947368421052630;
        coord(391, 1) = 0.8947368421052630;
        coord(392, 1) = 0.8947368421052630;
        coord(393, 1) = 0.8947368421052630;
        coord(394, 1) = 0.8947368421052630;
        coord(395, 1) = 0.8947368421052630;
        coord(396, 1) = 0.8947368421052630;
        coord(397, 1) = 0.8947368421052630;
        coord(398, 1) = 0.8947368421052630;
        coord(399, 1) = 0.8947368421052630;
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
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (rank == 0 and cells_.size() > 0)
    {
        if (analysisOfCrackPropagation_ == true)
        {
            exportToParaviewISO(std::to_string(0) + "_" + std::to_string(0));
        }
        else
        {
            exportToParaviewISO(std::to_string(-1));
        }
    }
    else if (rank == 1 and elements_.size() > 0)
    {
        if (analysisOfCrackPropagation_ == true)
        {
            exportToParaviewFEM(std::to_string(0) + "_" + std::to_string(0));
        }
        else
        {
            exportToParaviewFEM(std::to_string(-1));
        }
    }

    if (overlappingAnalysis_ == true)
    {
        if (localOverlapping_ == true)
        {
            setBlendingZoneLines();
        }
        computeDistanceFromFEBoundary();
        incidenceLocalxGlobal();
        checkInactivesCPandNode();
        shareDataBetweenRanks();

        if (plotHammerPointsInBlendingZone_ == true)
        {
            exportToParaviewHammerPoints();
        }
    }

    firstAccelerationCalculation();

    Mat A;
    Vec b, x, All;
    PetscErrorCode ierr;
    PetscInt Istart, Iend, Idof, Ione, iterations, *dof;
    KSP ksp;
    PC pc;
    VecScatter ctx;
    PetscScalar val, value;

    std::stringstream fileName;

    if (rank == 1)
    {
        fileName << "SIFs"
                 << ".txt";
    }

    std::ofstream sifFile(fileName.str());

    double length = 0.0;
    int auxprop = 0;

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

    bool testPropagation = true;
    int propagation = 0;
    double time = 0.0;

    for (int timeStep = 1; timeStep <= numberOfSteps_; timeStep++)
    {
        boost::posix_time::ptime t1 =
            boost::posix_time::microsec_clock::local_time();

        if (rank == 0)
        {
            std::cout << "------------------------- TIME STEP = "
                      << timeStep << " -------------------------\n";
        }

        int ndir = dirichletConditions_.size() + dirichletConditionsFE_.size() + 2 * inactiveCPandNode_.size();
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
        for (int index : inactiveCPandNode_)
        {
            dof[idir++] = 2 * index;
            dof[idir++] = 2 * index + 1;
        }

        for (int iteration = 0; iteration < maximumOfIteration_; iteration++) //definir o máximo de interações por passo de carga
        {
            //Create PETSc sparse parallel matrix
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                2 * cpnumber_, 2 * cpnumber_,
                                1800, NULL, 3000, NULL, &A);
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
                }
            }

            if (iteration == 0)
            {
                for (DirichletCondition *con : dirichletConditions_)
                {
                    ControlPoint *cp = con->getControlPoint();
                    int dir = con->getDirection();
                    double val1 = (con->getValue()); // / (1.0 * numberOfSteps_);

                    cp->incrementCurrentCoordinate(dir, val1);
                }
                for (DirichletConditionFE *con : dirichletConditionsFE_)
                {
                    Node *no = con->getNode();
                    int dir = con->getDirection();
                    double val1 = (con->getValue()); // / (1.0 * numberOfSteps_);

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
                elementMatrices = el->elementContributions(planeState_, "DYNAMIC", 1, 1, deltat_, beta_, gamma_, hammerPoints_, hammerPointsBlendZone_);
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
                // if (node->getIndex() == 2 and rank == 0)
                // {
                //     // std::cout << node->getCurrentCoordinate()(0) << " " << node->getCurrentCoordinate()(1) << std::endl;
                //     // std::cout << node->getPastCoordinate()(0) << " " << node->getPastCoordinate()(1) << std::endl;
                //     std::cout << accel(0) << " " << accel(1) << std::endl;
                //     std::cout << vel(0) << " " << vel(1) << std::endl;
                // }
            }

            for (InactiveNode *no : inactiveNode_)
            {
                no->interpolateGlobalCoordinate();
            }
            for (InactiveCP *cp : inactiveCP_)
            {
                cp->interpolateGlobalCoordinate();
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

        if (rank == 0 and cells_.size() > 0)
        {
            exportToParaviewISO(std::to_string(timeStep));
        }
        else if (rank == 1 and elements_.size() > 0)
        {
            stressCalculateFEM();
            exportToParaviewFEM(std::to_string(timeStep));
        }

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

        if (rank == 1)
        {
            // bounded_vector<double, 2> sifs = getSIFs();
            // sifFile << time << " " << sifs(0) << " " << sifs(1) << std::endl;

            // sifFile << time << " " << (controlPoints_[0]->getCurrentDisplacement()(1) + controlPoints_[4]->getCurrentDisplacement()(0)) * 0.5
            //         << " " << (controlPoints_[20]->getCurrentDisplacement()(1) + controlPoints_[24]->getCurrentDisplacement()(0)) * 0.5
            //         << " " << (nodes_[0]->getCurrentDisplacement()(0) + nodes_[3]->getCurrentDisplacement()(1)) * 0.5
            //         << " " << (nodes_[1]->getCurrentDisplacement()(0) + nodes_[2]->getCurrentDisplacement()(1)) * 0.5
            //         << std::endl;

            sifFile << controlPoints_[0]->getCurrentCoordinate()(0) << " " << controlPoints_[0]->getCurrentCoordinate()(1) << std::endl;

            time += deltat_;

            // Node *no = neumannConditionsFE_[0]->getNode();
            // bounded_vector<double, 2> coordNode = no->getCurrentCoordinate();
            // int cellIndex = no->getCellIndex();
            // if (cellIndex == -1)
            // {
            //     std::cout << "O NÓ NÃO ENCONTROU UMA CÉLULA.."
            //               << "\n";
            // }
            // else
            // {
            //     Cell *cell = cells_[cellIndex];
            //     bounded_vector<double, 2> xsiglobal = no->getXsisGlobal();
            //     bounded_vector<double, 2> blendF = blendFunction(no->getDistanceToBoundary());

            //     std::vector<ControlPoint *> connection = cell->getControlPoints();
            //     vector<double> wpc2(connection.size());
            //     bounded_vector<int, 2> INC_;
            //     INC_ = connection[connection.size() - 1]->getINC();
            //     for (int i = 0; i < connection.size(); i++)
            //     {
            //         wpc2(i) = connection[i]->getWeight();
            //     }
            //     vector<double> phi = cell->shapeFunction(xsiglobal, wpc2, INC_);

            //     bounded_vector<double, 2> coordInter;
            //     coordInter = (1.0 - blendF(0)) * coordNode;

            //     for (int i = 0; i < connection.size(); i++)
            //     {
            //         bounded_vector<double, 2> coordCP = connection[i]->getCurrentCoordinate();
            //         coordInter += blendF(0) * phi(i) * coordCP;
            //     }

            //     coordNode = coordInter;
            // }
            // sifFile << -coordNode(1) + no->getInitialCoordinate()(1) << " " << -no->getCurrentCoordinate()(1) + no->getInitialCoordinate()(1) << std::endl;
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

    int ndir = dirichletConditions_.size() + dirichletConditionsFE_.size() + 2 * inactiveCPandNode_.size();
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
    for (int index : inactiveCPandNode_)
    {
        dof[idir++] = 2 * index;
        dof[idir++] = 2 * index + 1;
    }

    //Create PETSc sparse parallel matrix
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        2 * cpnumber_, 2 * cpnumber_,
                        1800, NULL, 3000, NULL, &A);
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
        matrix<double> massLocal = el->massMatrix(quadrature_);

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
                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j), ADD_VALUES);

                dof1 = 2 * el->getControlPoint(i)->getIndex() + 1;
                dof2 = 2 * el->getControlPoint(j)->getIndex();
                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i + 1, 2 * j), ADD_VALUES);

                dof1 = 2 * el->getControlPoint(i)->getIndex();
                dof2 = 2 * el->getControlPoint(j)->getIndex() + 1;
                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j + 1), ADD_VALUES);

                dof1 = 2 * el->getControlPoint(i)->getIndex() + 1;
                dof2 = 2 * el->getControlPoint(j)->getIndex() + 1;
                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i + 1, 2 * j + 1), ADD_VALUES);
            }
        }
    }
    for (Element *el : elements_part)
    {
        std::pair<vector<double>, matrix<double>> elementMatrices;
        elementMatrices = el->elementContributions(planeState_, "STATIC", 1, 1, deltat_, beta_, gamma_, hammerPoints_, hammerPointsBlendZone_);
        matrix<double> massLocal = el->massMatrix(hammerPoints_, hammerPointsBlendZone_);
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
                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j), ADD_VALUES);

                dof1 = 2 * freedom[i] + 1;
                dof2 = 2 * freedom[j];
                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i + 1, 2 * j), ADD_VALUES);

                dof1 = 2 * freedom[i];
                dof2 = 2 * freedom[j] + 1;
                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i, 2 * j + 1), ADD_VALUES);

                dof1 = 2 * freedom[i] + 1;
                dof2 = 2 * freedom[j] + 1;
                ierr = MatSetValues(A, 1, &dof1, 1, &dof2, &massLocal(2 * i + 1, 2 * j + 1), ADD_VALUES);
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
    double norm = 0.0;
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

        bounded_vector<double, 2> inicialCoord = cp->getInitialCoordinate();
        cp->setCurrentCoordinate(inicialCoord);
    }
    for (Node *node : nodes_)
    {
        bounded_vector<double, 2> firstAccel;
        int i = node->getIndex();
        Idof = 2 * i;
        ierr = VecGetValues(All, Ione, &Idof, &val);
        CHKERRQ(ierr);
        firstAccel(0) = val;
        norm += val;

        Idof = 2 * i + 1;
        ierr = VecGetValues(All, Ione, &Idof, &val);
        CHKERRQ(ierr);
        firstAccel(1) = val;
        norm += val;

        node->setCurrentAcceleration(firstAccel);
        node->setPastAcceleration(firstAccel);

        bounded_vector<double, 2> inicialCoord = node->getInitialCoordinate();
        node->setCurrentCoordinate(inicialCoord);
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
    /*Looping para:
        -Encontrar pontos de hammer dentro da blending zone;
            -Célula em que o ponto de hammer está inserido;
            -Xsis correspondentes ao ponto de hammer na célula;
            -Computar valor de b e db_dxsiLocal
    */
    for (Element *el : elements_part) ///Looping para encontrar pontos de hammer dentro da blending zone;
    {
        std::vector<Node *> conec = el->getConnection();
        matrix<double> integrationPoints = el->hammerQuadrature(hammerPointsBlendZone_);
        bool blendZone = false; //test to check if the finite element is inside blend zone;

        //Variáveis para serem exportadas para os elementos
        std::vector<bool> insideBlendZone;                       //verdadeiro se o ponto de hammer está na blend zone e encontrou uma célula;
        std::vector<bounded_vector<double, 2>> xsiIncidenceCell; //coordenadas adimensionais do global que corresponde ao ponto de hammer
        std::vector<Cell *> incidenceCell;                       //célula que houve incidência do ponto de hammer
        std::vector<double> bValue;                              //value of assigned distance
        std::vector<bounded_vector<double, 2>> db_dxsiValue;     //derivada de b em relação a xsi

        for (int ih = 0; ih < hammerPointsBlendZone_; ih++)
        {
            double DAhammer = 0.0;             //distância interpolada do ponto de hammer ao contorno de elementos finitos
            bounded_vector<double, 2> coordFE; //coordenada global do ponto de hammer
            coordFE(0) = 0.0;
            coordFE(1) = 0.0;
            vector<double> phi = el->domainShapeFunction(integrationPoints(ih, 0), integrationPoints(ih, 1));

            for (int i = 0; i < conec.size(); i++)
            {
                DAhammer += phi(i) * conec[i]->getDistanceToBoundary();
                coordFE += phi(i) * conec[i]->getInitialCoordinate();
            }

            if (DAhammer <= blendingZoneThickness_) //inside of blend zone
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
                    int npc = conPoints.size();
                    vector<double> wpc(npc);
                    for (int i = 0; i < npc; i++)
                    {
                        wpc(i) = conPoints[i]->getWeight();
                    }
                    bounded_vector<int, 2> inc;
                    inc = conPoints[npc - 1]->getINC();

                    int cont = 0;
                    double erro = 1.0;
                    while (erro >= 1.0e-08 and cont <= 20)
                    {
                        std::pair<vector<double>, matrix<double>> functions;
                        functions = cell->shapeFunctionAndDerivates(xsiISO, wpc, inc);

                        coordISO(0) = 0.0;
                        coordISO(1) = 0.0;

                        for (int cp = 0; cp < npc; cp++)
                        {
                            bounded_vector<double, 2> coordinateCP = conPoints[cp]->getInitialCoordinate();
                            coordISO(0) += functions.first(cp) * coordinateCP(0);
                            coordISO(1) += functions.first(cp) * coordinateCP(1);
                        }

                        bounded_matrix<double, 2, 2> jacobian = cell->referenceJacobianMatrix(functions.second);
                        bounded_matrix<double, 2, 2> inverse = inverseMatrix(jacobian);
                        deltaxsi = prod(inverse, coordFE - coordISO);

                        xsiISO = xsiISO + deltaxsi;

                        erro = norm_2(deltaxsi);

                        cont++;
                    }

                    if (xsiISO(0) >= -1.0 and xsiISO(0) <= 1.0 and //ponto de hammer encontrou uma célula
                        xsiISO(1) >= -1.0 and xsiISO(1) <= 1.0 and erro <= 1.0e-08)
                    {
                        insideBlendZone[ih] = true;
                        blendZone = true;
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
                            dDA_dxsi(0) += dphi_dxsi(0, i) * DAnode; //dDAhammer_dxsiLocal1
                            dDA_dxsi(1) += dphi_dxsi(1, i) * DAnode; //dDAhammer_dxsiLocal1
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
        if (blendZone == true)
        {
            el->setIncidenceOnGlobal(insideBlendZone, incidenceCell, xsiIncidenceCell, bValue, db_dxsiValue);
        }
    }
}

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
    for (Element *el : elements_part)
    {
        std::vector<Node *> conec = el->getConnection();
        for (Node *no : conec)
        {
            if (no->getDistanceToBoundary() == -10000.0) //-10000.0 é o valor que a variável foi iniciada na criação do nó
            {
                double distanceNode = 1000000.0;
                bounded_vector<double, 2> coord = no->getInitialCoordinate();

                for (BoundaryElement *bound : blendingBoundary_)
                {
                    double distance;
                    double xsiBoundary = 0.0; //primeira tentativa
                    double deltaxsi = 100000.0;
                    int cont = 0;
                    std::vector<Node *> boundaryConec = bound->getNodes();
                    bounded_vector<double, 2> coordBoundary, firstDerivate, secondDerivate;

                    while ((fabs(deltaxsi) >= 1.0e-08) and (cont <= 20))
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
                            bounded_vector<double, 2> coordinateNode = node->getInitialCoordinate();

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
                            bounded_vector<double, 2> coordinateNode = node->getInitialCoordinate();
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
                            bounded_vector<double, 2> coordinateNode = node->getInitialCoordinate();
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

                bounded_vector<double, 2> coordFE = no->getInitialCoordinate();
                if (distanceNode <= blendingZoneThickness_)
                {
                    for (Cell *cell : cells_)
                    {
                        bounded_vector<double, 2> xsiISO, coordISO, deltaxsi;
                        deltaxsi(0) = 100000.0;
                        deltaxsi(1) = 100000.0;
                        xsiISO(0) = 0.0; //first attempt xs1
                        xsiISO(1) = 0.0; //first attempt xs2

                        std::vector<ControlPoint *> conPoints = cell->getControlPoints();
                        int npc = conPoints.size();
                        vector<double> wpc(npc);
                        for (int i = 0; i < npc; i++)
                        {
                            wpc(i) = conPoints[i]->getWeight();
                        }
                        bounded_vector<int, 2> inc;
                        inc = conPoints[npc - 1]->getINC();

                        int cont = 0;
                        double erro = 1.0;

                        while ((erro >= 1.0e-08) and cont <= 20)
                        {
                            std::pair<vector<double>, matrix<double>> functions;
                            functions = cell->shapeFunctionAndDerivates(xsiISO, wpc, inc);

                            coordISO(0) = 0.0;
                            coordISO(1) = 0.0;

                            for (int cp = 0; cp < npc; cp++)
                            {
                                bounded_vector<double, 2> coordinateCP = conPoints[cp]->getInitialCoordinate();
                                coordISO(0) += functions.first(cp) * coordinateCP(0);
                                coordISO(1) += functions.first(cp) * coordinateCP(1);
                            }

                            bounded_matrix<double, 2, 2> jacobian = cell->currentJacobianMatrix(functions.second);
                            bounded_matrix<double, 2, 2> inverse = inverseMatrix(jacobian);
                            deltaxsi = prod(inverse, coordFE - coordISO);

                            xsiISO = xsiISO + deltaxsi;

                            erro = norm_2(deltaxsi);

                            cont++;
                        }

                        if (xsiISO(0) >= -1.0 - tolerance_ and xsiISO(0) <= 1.0 + tolerance_ and //o nó encontrou uma célula
                            xsiISO(1) >= -1.0 - tolerance_ and xsiISO(1) <= 1.0 + tolerance_ and erro <= 1.0e-08)
                        {
                            no->setCellIndex(cell->getIndex());
                            no->setXsisGlobal(xsiISO);

                            break;
                        }
                    }
                }
            }
        }
    }

    for (Cell *cell : cells_part)
    {
        matrix<double> integrationPoints = cell->isoQuadrature(quadratureBlendingZone_);
        std::vector<double> distanceFE;
        bool in = false;

        for (int ip = 0; ip < integrationPoints.size1(); ip++)
        {
            distanceFE.push_back(-100000000.0);
            bounded_vector<double, 2> xsi, coordIP; //coordinate of hammer point
            xsi(0) = integrationPoints(ip, 0);
            xsi(1) = integrationPoints(ip, 1);
            coordIP = cell->calculateGlobalCoordinate(xsi);
            double distance;

            for (BoundaryElement *bound : blendingBoundary_)
            {
                double xsiBoundary = 0.0; //primeira tentativa
                double deltaxsi = 100000.0;
                int cont = 0;
                std::vector<Node *> boundaryConec = bound->getNodes();
                bounded_vector<double, 2> coordBoundary, firstDerivate, secondDerivate, normal;

                while (fabs(deltaxsi) >= 1.0e-08 and cont <= 20)
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

                if (fabs(distance) < 1.0e-10)
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

        for (int ipp = 0; ipp < distanceFE.size(); ipp++)
        {
            if (distanceFE[ipp] > 0.0)
            {
                in = true;
                break;
            }
        }

        if (in == true) //elemento tem pelo menos um ponto dentro da malha local
        {
            cell->setDistanceFromFEBoundary(distanceFE);
        }
        else
        {
            distanceFE.clear();
            cell->setDistanceFromFEBoundary(distanceFE);
        }
    }
}

void GlobalSolid::exportToParaviewHammerPoints()
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Barrier(PETSC_COMM_WORLD);

    std::stringstream name;
    if (rank == 0)
    {
        name << "hammerPointsFE.vtu";
    }
    std::ofstream file(name.str());

    int nn = 0;
    for (int i = 0; i < elements_.size(); i++)
    {
        int teste;
        if (elementPartition_[i] == rank)
        {
            if (elements_[i]->ipInsideBlendZone().size() > 0)
            {
                teste = 1;
            }
            else
            {
                teste = 0;
            }
        }
        MPI_Bcast(&teste, 1, MPI_INT, elementPartition_[i], PETSC_COMM_WORLD);
        MPI_Barrier(PETSC_COMM_WORLD);
        if (teste == 0)
        {
            nn += hammerPoints_;
        }
        else
        {
            nn += hammerPointsBlendZone_;
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }

    //matrix<double> hammer = elements_[0]->hammerQuadrature();
    //int nn = hammer.size1();

    //header
    if (rank == 0)
    {
        file << "<?xml version=\"1.0\"?>"
             << "\n"
             << "<VTKFile type=\"UnstructuredGrid\">"
             << "\n"
             << "  <UnstructuredGrid>"
             << "\n"
             << "  <Piece NumberOfPoints=\"" << nn
             << "\"  NumberOfCells=\"" << nn
             << "\">"
             << "\n";

        //nodal coordinates
        file << "    <Points>"
             << "\n"
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">"
             << "\n";
    }

    for (int i = 0; i < elements_.size(); i++)
    {
        matrix<double> hammer;
        int teste;
        if (elementPartition_[i] == rank)
        {
            if (elements_[i]->ipInsideBlendZone().size() > 0)
            {
                teste = 1;
            }
            else
            {
                teste = 0;
            }
        }
        MPI_Bcast(&teste, 1, MPI_INT, elementPartition_[i], PETSC_COMM_WORLD);
        MPI_Barrier(PETSC_COMM_WORLD);
        if (teste == 0)
        {
            hammer = elements_[i]->hammerQuadrature(hammerPoints_);
        }
        else
        {
            hammer = elements_[i]->hammerQuadrature(hammerPointsBlendZone_);
        }
        if (rank == 0)
        {
            for (int j = 0; j < hammer.size1(); j++)
            {
                bounded_vector<double, 2> xsi;
                xsi(0) = hammer(j, 0);
                xsi(1) = hammer(j, 1);
                bounded_vector<double, 2> coord = elements_[i]->calculateGlobalCoordinate(xsi);
                file << coord(0) << " " << coord(1) << " " << 1.0 << "\n";
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }
    if (rank == 0)
    {
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

        for (int i = 0; i < nn; i++)
        {
            file << aux++ << std::endl;
        }

        file << "      </DataArray>"
             << "\n";

        //offsets
        file << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">"
             << "\n";

        aux = 0;

        for (int i = 0; i < nn; i++)
        {
            file << aux++ << std::endl;
        }

        file << "      </DataArray>"
             << "\n";

        //elements type
        file << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">"
             << "\n";

        for (int i = 0; i < nn; i++)
        {
            file << 1 << std::endl;
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
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    for (int i = 0; i < elements_.size(); i++)
    {
        int teste;
        std::vector<bool> blend;
        std::vector<double> b;
        double bvalue;
        int inside;

        if (elementPartition_[i] == rank)
        {
            if (elements_[i]->ipInsideBlendZone().size() > 0)
            {
                teste = 1;
                blend = elements_[i]->ipInsideBlendZone();
                b = elements_[i]->getBlendValue();
            }
            else
            {
                teste = 0;
            }
        }
        MPI_Bcast(&teste, 1, MPI_INT, elementPartition_[i], PETSC_COMM_WORLD);
        MPI_Barrier(PETSC_COMM_WORLD);

        if (teste == 0 and rank == 0)
        {
            for (int j = 0; j < hammerPoints_; j++)
            {
                file << 0 << " " << 0.0 << " " << 1.0 << "\n";
            }
        }

        if (teste == 1)
        {
            int aux = 0;
            for (int h = 0; h < hammerPointsBlendZone_; h++)
            {
                if (elementPartition_[i] == rank)
                {
                    inside = blend[h];
                    if (inside == 1)
                    {
                        bvalue = b[aux++];
                    }
                    else
                    {
                        bvalue = 0;
                    }
                }
                MPI_Bcast(&bvalue, 1, MPI_DOUBLE, elementPartition_[i], PETSC_COMM_WORLD);
                MPI_Bcast(&inside, 1, MPI_INT, elementPartition_[i], PETSC_COMM_WORLD);
                MPI_Barrier(PETSC_COMM_WORLD);
                if (rank == 0)
                {
                    file << inside << " " << bvalue << " " << 1.0 - bvalue << "\n";
                }
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }

    if (rank == 0)
    {
        file << "      </DataArray> "
             << "\n";

        file << "    </PointData>"
             << "\n";
        //elemental results
        file << "    <CellData>"
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
}

bounded_vector<double, 2> GlobalSolid::blendFunction(const double &DA)
{
    bounded_vector<double, 2> blend; //(b, db_dDA)
    double aux = DA / blendingZoneThickness_;

    if (blendingFunctionOrder_ == 3)
    {
        blend(0) = 2.0 * aux * aux * aux - 3.0 * aux * aux + 1.0;
        blend(1) = (6.0 / blendingZoneThickness_) * aux * aux - (6.0 / blendingZoneThickness_) * aux;
    }
    else if (blendingFunctionOrder_ == 5)
    {
        blend(0) = -6.0 * pow(aux, 5) + 15.0 * pow(aux, 4) - 10.0 * pow(aux, 3) + 1.0;
        blend(1) = -30.0 / blendingZoneThickness_ * pow(aux, 4) + 60 / blendingZoneThickness_ * pow(aux, 3) - 30.0 / blendingZoneThickness_ * pow(aux, 2);

        // blend(0) = pow(aux, 5.0) - 2.5 * pow(aux, 4.0) + 4.0 * pow(aux, 3.0) - 3.5 * pow(aux, 2.0) + 1.0;
        // blend(1) = 5.0 / blendingZoneThickness_ * pow(aux, 4.0) - 10.0 / blendingZoneThickness_ * pow(aux, 3.0) + 12.0 / blendingZoneThickness_ * pow(aux, 2.0) - 7.0 / blendingZoneThickness_ * pow(aux, 1.0);

        // blend(0) = aux * aux * aux * aux * aux - 2.5 * aux * aux * aux * aux + 4.0 * aux * aux * aux - 3.5 * aux * aux + 1.0;
        // blend(1) = (5.0 / blendingZoneThickness_) * aux * aux * aux * aux - (10.0 / blendingZoneThickness_) * aux * aux * aux + (12.0 / blendingZoneThickness_) * aux * aux - (7.0 / blendingZoneThickness_) * aux;
    }
    else if (blendingFunctionOrder_ == 4)
    {
        if (DA <= blendingZoneThickness_ * 0.5)
        {
            blend(0) = 8.0 * aux * aux * aux * aux - 8.0 * aux * aux * aux + 1.0;
            blend(1) = (32.0 / blendingZoneThickness_) * aux * aux * aux - (24.0 / blendingZoneThickness_) * aux * aux;
        }
        else
        {
            double dif = 1.0 - aux;
            blend(0) = -8.0 * dif * dif * dif * dif + 8.0 * dif * dif * dif;
            blend(1) = (32.0 / blendingZoneThickness_) * dif * dif * dif - (24.0 / blendingZoneThickness_) * dif * dif;
        }
        // blend(0) = aux * aux * aux * aux - 2.0 * aux * aux + 1.0;
        // blend(1) = (4.0 / blendingZoneThickness_) * aux * aux * aux - (4.0 / blendingZoneThickness_) * aux;

        // blend(0) = 2.0 * aux * aux * aux * aux - 2.0 * aux * aux * aux - aux * aux + 1.0;
        // blend(1) = (8.0 / blendingZoneThickness_) * aux * aux * aux - (6.0 / blendingZoneThickness_) * aux * aux - (2.0 / blendingZoneThickness_) * aux;

        // blend(0) = 3.0 * aux * aux * aux * aux - 4.0 * aux * aux * aux + 1.0;
        // blend(1) = (12.0 / blendingZoneThickness_) * aux * aux * aux - (12.0 / blendingZoneThickness_) * aux * aux;
    }
    else if (blendingFunctionOrder_ == 1)
    {
        blend(0) = -aux + 1.0;
        blend(1) = -1.0 / blendingZoneThickness_;
    }
    else if (blendingFunctionOrder_ == 2)
    {
        blend(0) = aux * aux - 2.0 * aux + 1.0;
        blend(1) = (2.0 / blendingZoneThickness_) * aux - (2.0 / blendingZoneThickness_);
    }
    else
    {
        std::cout << "ERRO NA ORDEM DA BLENDING FUNCTION!!!" << std::endl;
    }

    return blend;
}

void GlobalSolid::checkInactivesCPandNode()
{
    Vec b, All;
    PetscInt Idof, Ione = 1;
    PetscScalar val;
    VecScatter ctx;

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    VecCreate(PETSC_COMM_WORLD, &b);
    VecSetSizes(b, PETSC_DECIDE, cpnumber_);
    VecSetFromOptions(b);

    for (Cell *cell : cells_part)
    {
        vector<double> mass = cell->diagonalMass(quadrature_);
        std::vector<ControlPoint *> conec = cell->getControlPoints();
        for (int i = 0; i < conec.size(); i++)
        {
            int cp = conec[i]->getIndex();
            VecSetValues(b, 1, &cp, &mass(i), ADD_VALUES);
        }
    }
    for (Element *el : elements_part)
    {
        vector<double> mass = el->diagonalMass(hammerPoints_, hammerPointsBlendZone_);
        std::vector<int> freedom = el->getFreedomDegree();
        for (int i = 0; i < freedom.size(); i++)
        {
            int dof = freedom[i];
            VecSetValues(b, 1, &dof, &mass(i), ADD_VALUES);
        }
    }

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);

    VecScatterCreateToAll(b, &ctx, &All);
    VecScatterBegin(ctx, b, All, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, b, All, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&ctx);

    std::unordered_set<ControlPoint *> inactiveCP;
    std::unordered_set<int> inactiveCPaux;

    std::unordered_set<int> inactiveNode;

    double min = 1.0e200;
    double max = 1.0e-200;

    double limN = 1.0e-15;

    for (Node *node : nodes_)
    {
        Idof = node->getIndex();
        VecGetValues(All, Ione, &Idof, &val);

        if (val < min)
        {
            min = val;
        }
        if (val > max)
        {
            max = val;
        }

        if (val < limN)
        {
            inactiveNode.insert(node->getIndexFE());
            inactiveCPandNode_.insert(Idof);
        }
    }

    double limCP = min / 10000.0;

    if (rank == 0)
        std::cout << "phi*phi\nNode - max: " << max << " // min: " << min << std::endl;

    min = 1.0e200;
    max = 1.0e-200;
    for (ControlPoint *cp : controlPoints_)
    {
        Idof = cp->getIndex();
        VecGetValues(All, Ione, &Idof, &val);
        if (val < min and val != 0.0)
        {
            min = val;
        }
        if (val > max)
        {
            max = val;
        }
        if (val <= limCP)
        {
            if (inactiveCPaux.count(Idof) == 0)
            {
                inactiveCP.insert(cp);
                inactiveCPaux.insert(Idof);
                inactiveCPandNode_.insert(Idof);
            }
        }
    }

    if (rank == 0)
        std::cout << "CP - max: " << max << " // min: " << min << "\n\n";

    // for (Node *node : nodes_)
    // {
    //     Idof = node->getIndex();
    //     VecGetValues(All, Ione, &Idof, &val);

    //     if (val < tolerance_)
    //     {
    //         inactiveNode.insert(node->getIndexFE());
    //         inactiveCPandNode_.insert(Idof);
    //         // std::cout << "The node " << node->getIndex() << " was desactived." << std::endl;
    //     }
    // }

    VecDestroy(&b);

    for (ControlPoint *cp : inactiveCP)
    {
        if (rank == 0)
        {
            Idof = cp->getIndex();
            VecGetValues(All, Ione, &Idof, &val);
            std::cout << "The CP " << cp->getIndex() << "(" << cp->getInitialCoordinate()(0) / cp->getWeight() << " " << cp->getInitialCoordinate()(1) / cp->getWeight() << ") "
                      << " was desactived ( " << val << " )." << std::endl;
        }
        bounded_vector<double, 2> coord = cp->getInitialCoordinate();

        for (Element *el : elements_)
        {
            bounded_vector<double, 2> xsiFE, coordFE, deltaxsi;
            deltaxsi(0) = 100000.0;
            deltaxsi(1) = 100000.0;
            xsiFE(0) = 0.5; //first attempt xs1
            xsiFE(1) = 0.5; //first attempt xs2

            std::vector<Node *> conec = el->getConnection();

            int cont = 0;
            double error = 1.0;
            while (error >= 1.0e-08 and cont <= 15)
            {

                coordFE = el->calculateGlobalCoordinate(xsiFE);

                bounded_matrix<double, 2, 2> jacobian = el->referenceJacobianMatrix(xsiFE(0), xsiFE(1));
                bounded_matrix<double, 2, 2> inverse = inverseMatrix(jacobian);
                deltaxsi = prod(inverse, coord / cp->getWeight() - coordFE);

                xsiFE = xsiFE + deltaxsi;

                error = norm_2(deltaxsi);

                cont++;
            }

            if (xsiFE(0) >= -tolerance_ and xsiFE(0) <= 1.0 + tolerance_ and //ponto de hammer encontrou uma célula
                xsiFE(1) >= -tolerance_ and xsiFE(1) <= 1.0 + tolerance_ and
                xsiFE(0) + xsiFE(1) <= 1.0 + tolerance_ and error <= 1.0e-08) //só vai funcionar para triângulos não muito deformados...
            {
                InactiveCP *inactive = new InactiveCP(cp, el, xsiFE);
                inactiveCP_.push_back(inactive);

                if (rank == 0)
                {
                    std::cout << "CP " << cp->getIndex() << " found the element " << el->getIndex() << std::endl;
                }

                break;
            }
        }
    }

    for (int nn : inactiveNode)
    {
        Node *no = nodes_[nn];
        Idof = no->getIndex();
        VecGetValues(All, Ione, &Idof, &val);
        if (rank == 0)
            std::cout << "NODE " << no->getIndex() << "(" << no->getInitialCoordinate()(0) << " " << no->getInitialCoordinate()(1) << ") "
                      << " was desactived ( " << val << " )." << std::endl;
        bounded_vector<double, 2> coordFE = no->getInitialCoordinate();

        for (Cell *cell : cells_)
        {
            bounded_vector<double, 2> xsiISO, coordISO, deltaxsi;
            deltaxsi(0) = 100000.0;
            deltaxsi(1) = 100000.0;
            xsiISO(0) = 0.0; //first attempt xs1
            xsiISO(1) = 0.0; //first attempt xs2

            std::vector<ControlPoint *> conPoints = cell->getControlPoints();
            int npc = conPoints.size();
            vector<double> wpc(npc);
            for (int i = 0; i < npc; i++)
            {
                wpc(i) = conPoints[i]->getWeight();
            }
            bounded_vector<int, 2> inc;
            inc = conPoints[npc - 1]->getINC();

            int cont = 0;
            double error = 1.0;
            while (error >= 1.0e-08 and cont <= 15)
            {
                std::pair<vector<double>, matrix<double>> functions;
                functions = cell->shapeFunctionAndDerivates(xsiISO, wpc, inc);

                coordISO(0) = 0.0;
                coordISO(1) = 0.0;

                for (int cp = 0; cp < npc; cp++)
                {
                    bounded_vector<double, 2> coordinateCP = conPoints[cp]->getInitialCoordinate();
                    coordISO(0) += functions.first(cp) * coordinateCP(0);
                    coordISO(1) += functions.first(cp) * coordinateCP(1);
                }

                bounded_matrix<double, 2, 2> jacobian = cell->referenceJacobianMatrix(functions.second);
                bounded_matrix<double, 2, 2> inverse = inverseMatrix(jacobian);
                deltaxsi = prod(inverse, coordFE - coordISO);

                xsiISO = xsiISO + deltaxsi;

                error = norm_2(deltaxsi);

                cont++;
            }

            if (xsiISO(0) >= -1.0 - 1.0e-08 and xsiISO(0) <= 1.0 + 1.0e-08 and //ponto de hammer encontrou uma célula
                xsiISO(1) >= -1.0 - 1.0e-08 and xsiISO(1) <= 1.0 + 1.0e-08 and error <= 1.0e-08)
            {

                InactiveNode *inactive = new InactiveNode(no, cell, xsiISO);
                inactiveNode_.push_back(inactive);

                if (rank == 0)
                {
                    std::cout << "NODE " << no->getIndex() << " found the cell " << cell->getIndex() << std::endl;
                }

                break;
            }
        }
    }
}

void GlobalSolid::shareDataBetweenRanks()
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    for (int el = 0; el < elements_.size(); el++) //compartilhar distância dos nós ao contorno, célula em que o nó está inserido e xsisGlobal correspondentes com o RANK 1, que irá necessitar de todos quando for exportar os dados
    {
        if (elementPartition_[el] != 1)
        {
            std::vector<Node *> conec = elements_[el]->getConnection();

            vector<int> cellIndex(conec.size(), -1);
            vector<double> distanceBoundary(conec.size(), 10000.0);
            vector<double> xsisGlobal1(conec.size(), 10000.0);
            vector<double> xsisGlobal2(conec.size(), 10000.0);

            bounded_vector<double, 2> xsi;

            if (elementPartition_[el] == rank)
            {
                for (int i = 0; i < conec.size(); i++)
                {
                    xsi = conec[i]->getXsisGlobal();
                    cellIndex(i) = conec[i]->getCellIndex();
                    distanceBoundary(i) = conec[i]->getDistanceToBoundary();
                    xsisGlobal1(i) = xsi(0);
                    xsisGlobal2(i) = xsi(1);
                }
            }
            //std::cout<<"TESETE 1   "<< rank<<std::endl;
            MPI_Bcast(&cellIndex[0], conec.size(), MPI_INT, elementPartition_[el], PETSC_COMM_WORLD);
            MPI_Bcast(&distanceBoundary[0], conec.size(), MPI_DOUBLE, elementPartition_[el], PETSC_COMM_WORLD);
            MPI_Bcast(&xsisGlobal1[0], conec.size(), MPI_DOUBLE, elementPartition_[el], PETSC_COMM_WORLD);
            MPI_Bcast(&xsisGlobal2[0], conec.size(), MPI_DOUBLE, elementPartition_[el], PETSC_COMM_WORLD);
            //std::cout<<"TESETE 2   "<<rank<< std::endl;
            MPI_Barrier(PETSC_COMM_WORLD);

            if (rank == 1)
            {
                for (int i = 0; i < conec.size(); i++)
                {
                    if (conec[i]->getDistanceToBoundary() == -10000.0 and distanceBoundary(i) != -10000.0)
                    {
                        conec[i]->setCellIndex(cellIndex(i));
                        conec[i]->setDistanceToBoundary(distanceBoundary(i));
                        xsi(0) = xsisGlobal1(i);
                        xsi(1) = xsisGlobal2(i);
                        conec[i]->setXsisGlobal(xsi);
                    }
                }
            }
        }
    }
    if (rank == 1)
    {
        for (Node *no : nodes_)
        {
            no->setValuesOfBlendingFunction(blendFunction(no->getDistanceToBoundary()));
        }
    }
}

matrix<double> GlobalSolid::coordinatesFEM(const int &order)
{
    int n = (order + 1) * (order + 2) / 2.0;
    matrix<double> xsi(n, 2);

    if (order == 1)
    {
        xsi(0, 0) = 1.0;
        xsi(0, 1) = 0.0;

        xsi(1, 0) = 0.0;
        xsi(1, 1) = 1.0;

        xsi(2, 0) = 0.0;
        xsi(2, 1) = 0.0;
    }
    else if (order == 2)
    {
        xsi(0, 0) = 1.0;
        xsi(0, 1) = 0.0;

        xsi(1, 0) = 0.0;
        xsi(1, 1) = 1.0;

        xsi(2, 0) = 0.0;
        xsi(2, 1) = 0.0;

        xsi(3, 0) = 0.5;
        xsi(3, 1) = 0.5;

        xsi(4, 0) = 0.0;
        xsi(4, 1) = 0.5;

        xsi(5, 0) = 0.5;
        xsi(5, 1) = 0.0;
    }
    else
    {
        xsi(0, 0) = 1.0;
        xsi(0, 1) = 0.0;

        xsi(1, 0) = 0.0;
        xsi(1, 1) = 1.0;

        xsi(2, 0) = 0.0;
        xsi(2, 1) = 0.0;

        xsi(3, 0) = 0.666666666666667;
        xsi(3, 1) = 0.333333333333333;

        xsi(4, 0) = 0.333333333333333;
        xsi(4, 1) = 0.666666666666667;

        xsi(5, 0) = 0.0;
        xsi(5, 1) = 0.666666666666667;

        xsi(6, 0) = 0.0;
        xsi(6, 1) = 0.333333333333333;

        xsi(7, 0) = 0.333333333333333;
        xsi(7, 1) = 0.0;

        xsi(8, 0) = 0.666666666666667;
        xsi(8, 1) = 0.0;

        xsi(9, 0) = 0.333333333333333;
        xsi(9, 1) = 0.333333333333333;
    }
    return xsi;
}

void GlobalSolid::stressCalculateFEM()
{
    for (Node *n : nodes_)
    {
        n->setZeroStressState();
    }

    for (Element *el : elements_)
    {
        std::vector<bool> insideBlendZone = el->ipInsideBlendZone();
        std::vector<Node *> connection = el->getConnection();
        int nnodes = connection.size();
        if (insideBlendZone.size() == 0) //o elemento não tem nenhum ponto de hammer na blending zone
        {
            int nh;

            if (nnodes == 6)
            {
                nh = 7;
            }
            else if (nnodes == 10)
            {
                nh = 12;
            }

            matrix<double> pointCoord(nh, 3);

            pointCoord = el->hammerQuadrature(nh);

            matrix<double, column_major> phiMatrix(nnodes, nnodes, 0.0);
            matrix<double, column_major> cauchyStress(nnodes, 4, 0.0);

            for (int i = 0; i < nnodes; i++)
            {
                bounded_vector<double, 4> stress;

                bounded_vector<double, 2> qxsi;
                qxsi(0) = pointCoord(i, 0);
                qxsi(1) = pointCoord(i, 1);

                vector<double> phi = el->domainShapeFunction(qxsi(0), qxsi(1));

                for (int j = 0; j < nnodes; j++)
                {
                    phiMatrix(i, j) = phi(j);
                }

                stress = el->getCauchyStress(qxsi, planeState_);

                cauchyStress(i, 0) = stress(0);
                cauchyStress(i, 1) = stress(1);
                cauchyStress(i, 2) = stress(2);
                cauchyStress(i, 3) = stress(3);
            }

            vector<int> c(nnodes, 0);

            boost::numeric::bindings::lapack::gesv(phiMatrix, c, cauchyStress);

            for (int i = 0; i < nnodes; i++)
            {
                bounded_vector<double, 4> stress;
                stress(0) = cauchyStress(i, 0);
                stress(1) = cauchyStress(i, 1);
                stress(2) = cauchyStress(i, 2);
                stress(3) = cauchyStress(i, 3);
                connection[i]->setStressState(stress);
            }
        }
        else
        {
            int order;
            if (nnodes == 3)
            {
                order = 1;
            }
            else if (nnodes == 6)
            {
                order = 2;
            }
            else if (nnodes == 10)
            {
                order = 3;
            }

            matrix<double> xsiLocal = coordinatesFEM(order);
            for (int i = 0; i < nnodes; i++)
            {
                bounded_vector<double, 2> qxsi;
                qxsi(0) = xsiLocal(i, 0);
                qxsi(1) = xsiLocal(i, 1);

                int cellIndex = connection[i]->getCellIndex();

                if (cellIndex == -1)
                {
                    connection[i]->setStressState(el->getCauchyStress(qxsi, planeState_));
                }
                else
                {
                    Cell *cell = cells_[cellIndex];
                    connection[i]->setStressState(el->getCauchyStressBlendZone(qxsi, planeState_, cell, connection[i]->getXsisGlobal(), connection[i]->getValuesOfBlendingFunction()));
                }
            }
        }
    }
}

void GlobalSolid::generateMesh(Geometry *geometry, const std::string &elementType, const std::string &algorithm, std::string geofile,
                               const bool &plotMesh, const bool &showInfo)
{
    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    MPI_Barrier(PETSC_COMM_WORLD);

    boost::posix_time::ptime geoIni =
        boost::posix_time::microsec_clock::local_time();

    finiteElementType_ = elementType;
    geometry_ = geometry;
    typeAlgorithm_ = algorithm;

    if (geometry->getTypeOfDomain() == "LOCAL")
    {
        geometry->createGeometryFromCrack();
    }

    std::string mshfile;
    bool deleteFiles = false;

    if (geofile.empty())
    {
        mshfile = "temp.msh";
        geofile = "temp.geo";
        deleteFiles = true;
    }
    else
    {
        mshfile = geofile + ".msh";
        geofile = geofile + ".geo";
    }

    if (rank == 0)
    {
        std::string gmshCode = geometry->createGmshCode();

        // gmshCode += "Transfinite Curve {l0, l2} = 7 Using Progression 1;\nTransfinite Curve {l1, l3} = 2 Using Progression 1;\nTransfinite Surface {s0} Left;\n";

        if (elementType == "T10")
        {
            gmshCode += "Mesh.ElementOrder = 3;\n//\n";
        }
        else if (elementType == "T6")
        {
            gmshCode += "Mesh.ElementOrder = 2;\n//\n";
        }
        else if (elementType == "T3")
        {
            gmshCode += "Mesh.ElementOrder = 1;\n//\n";
        }
        else
        {
            std::cout << elementType << " is not supported. Please select another type of finite element. \n ";
            exit(EXIT_FAILURE);
        }

        if (algorithm == "AUTO")
        {
            gmshCode += "Mesh.Algorithm = 2;\n//\n";
        }
        else if (algorithm == "DELAUNAY")
        {
            gmshCode += "Mesh.Algorithm = 5;\n//\n";
        }

        else if (algorithm == "FRONT")
        {
            gmshCode += "Mesh.Algorithm = 6;\n//\n";
        }
        else if (algorithm == "ADAPT")
        {
            gmshCode += "Mesh.Algorithm = 1;\n//\n";
        }
        else if (algorithm == "BAMG")
        {
            gmshCode += "Mesh.Algorithm = 7;\n//\n";
        }
        else
        {
            std::cout << elementType << " is not supported. Please select another type of algorithm. \n ";
            exit(EXIT_FAILURE);
        }

        gmshCode += "Mesh 2;\n//\n";

        std::unordered_map<std::string, Crack *> crackes = geometry->getCrackes();

        for (std::unordered_map<std::string, Crack *>::const_iterator c = crackes.begin(); c != crackes.end(); c++)
        {
            gmshCode += c->second->getGmshCodeCrackPlugin();
        }

        gmshCode += "Mesh.MshFileVersion = 2.2;\n";
        gmshCode += "Save \"" + mshfile + "\";\n";

        std::ofstream file(geofile);
        file << gmshCode;
        file.close();

        std::string cmd = current_working_dir_ + "/Solid/Mesh/gmsh ";
        if (plotMesh)
        {
            cmd += geofile; // + " &";
        }
        else
        {
            cmd += geofile + " -";
        }

        if (!showInfo and !plotMesh)
        {
            cmd += " -v 0";
        }

        system(cmd.c_str());
    }

    MPI_Barrier(PETSC_COMM_WORLD);

    boost::posix_time::ptime geoFin =
        boost::posix_time::microsec_clock::local_time();

    dataFromGmsh(mshfile, elementType, geometry);

    boost::posix_time::ptime dataFin =
        boost::posix_time::microsec_clock::local_time();

    applyBoundaryConditions(geometry);

    MPI_Barrier(PETSC_COMM_WORLD);

    if (rank == 0 and deleteFiles)
    {
        std::string aux = "temp*";
        system((remove + aux).c_str());
    }

    boost::posix_time::ptime metisIni =
        boost::posix_time::microsec_clock::local_time();

    if (size > 1)
        domainDecompositionMETIS(elementType);

    boost::posix_time::ptime metisFin =
        boost::posix_time::microsec_clock::local_time();

    boost::posix_time::time_duration diff = geoFin - geoIni;
    boost::posix_time::time_duration diff1 = dataFin - geoFin;
    boost::posix_time::time_duration diff2 = metisFin - metisIni;

    time__ << diff.total_milliseconds() / 1000. << " " << diff1.total_milliseconds() / 1000. << " " << diff2.total_milliseconds() / 1000. << "\n";

    if (rank == 0)
        std::cout << time__.str() << std::endl;

    if (quarterPointElement_ == true)
    {
        quarterPointElements();
    }

    if (size > 1)
    {
        for (int i = 0; i < elements_.size(); i++)
        {
            if (rank == elementPartition_[i])
            {
                elements_part.push_back(elements_[i]);
            }
        }
    }
    else
    {
        elementPartition_ = new idx_t[elements_.size()];
        nodePartition_ = new idx_t[nodes_.size()];
        for (int i = 0; i < elements_.size(); i++)
        {
            elements_part.push_back(elements_[i]);
            elementPartition_[i] = 0;
        }
    }

    MPI_Barrier(PETSC_COMM_WORLD);

    // for (Node *n : nodes_)
    // {
    //     bounded_vector<double, 2> initial = n->getInitialCoordinate();
    //     double angle = atan2(initial(1), initial(0));
    //     double radius = norm_2(initial);

    //     double curX, curY;

    //     double pi = 3.14159265358979323846;

    //     curX = cos(angle - 0.5*pi) * radius;
    //     curY = sin(angle - 0.5*pi) * radius;

    //     initial(0) = curX;
    //     initial(1) = curY;

    //    n->setCurrentCoordinate(initial);
    // }
}

void GlobalSolid::applyBoundaryConditions(Geometry *geometry)
{
    std::unordered_map<std::string, bounded_vector<double, 2>> condition = geometry->getNeumannCondition();

    for (std::unordered_map<std::string, bounded_vector<double, 2>>::const_iterator neumman = condition.begin(); neumman != condition.end(); neumman++)
    {
        std::string name = neumman->first;
        bounded_vector<double, 2> force = neumman->second;
        if (name[0] == 'p')
        {
            Node *no = geometry->getPoint(name)->getPointNode();
            if (force[0] != 0.0)
            {
                NeumannConditionFE *neu0 = new NeumannConditionFE(no, 0, force(0));
                neumannConditionsFE_.push_back(neu0);
            }
            if (force[1] != 0.0)
            {
                NeumannConditionFE *neu1 = new NeumannConditionFE(no, 1, force(1));
                neumannConditionsFE_.push_back(neu1);
            }
        }
        else
        {
            std::vector<BoundaryElement *> boundLine = finiteElementBoundary_[name];
            for (BoundaryElement *bound : boundLine)
            {
                vector<double> forcevec = bound->computeDistribuitedLoads(force, 5);
                std::vector<Node *> conec = bound->getNodes();
                double value0, value1;

                for (int i = 0; i < conec.size(); i++)
                {
                    value0 = forcevec(2 * i);
                    value1 = forcevec(2 * i + 1);

                    if (fabs(value0) >= 1.0e-24)
                    {
                        bool exist = false;
                        for (NeumannConditionFE *neu : neumannConditionsFE_)
                        {
                            if (neu->getNode() == conec[i] and neu->getDirection() == 0)
                            {
                                neu->incrementValue(value0);
                                exist = true;
                            }
                        }

                        if (exist == false)
                        {
                            NeumannConditionFE *neu0 = new NeumannConditionFE(conec[i], 0, value0);
                            neumannConditionsFE_.push_back(neu0);
                        }
                    }
                    if (fabs(value1) >= 1.0e-24)
                    {
                        bool exist = false;
                        for (NeumannConditionFE *neu : neumannConditionsFE_)
                        {
                            if (neu->getNode() == conec[i] and neu->getDirection() == 1)
                            {
                                neu->incrementValue(value1);
                                exist = true;
                            }
                        }
                        if (exist == false)
                        {
                            NeumannConditionFE *neu1 = new NeumannConditionFE(conec[i], 1, value1);
                            neumannConditionsFE_.push_back(neu1);
                        }
                    }
                }
            }
        }
    }

    condition = geometry->getDirichletCondition();

    for (std::unordered_map<std::string, bounded_vector<double, 2>>::const_iterator dirichlet = condition.begin(); dirichlet != condition.end(); dirichlet++)
    {
        std::string name = dirichlet->first;
        bounded_vector<double, 2> desloc = dirichlet->second;

        if (name[0] == 'p')
        {
            Node *no = geometry->getPoint(name)->getPointNode();
            if (desloc(0) != 1.0e-240)
            {
                DirichletConditionFE *dir = new DirichletConditionFE(no, 0, desloc(0));
                dirichletConditionsFE_.push_back(dir);
            }

            if (desloc(1) != 1.0e-240)
            {
                DirichletConditionFE *dir = new DirichletConditionFE(no, 1, desloc(1));
                dirichletConditionsFE_.push_back(dir);
            }
        }
        else
        {
            std::vector<BoundaryElement *> boundLine = finiteElementBoundary_[name];
            std::unordered_set<Node *> nodes;
            for (BoundaryElement *bound : boundLine)
            {
                std::vector<Node *> conec = bound->getNodes();
                for (Node *no : conec)
                {
                    nodes.insert(no);
                }
            }
            for (Node *no : nodes)
            {
                if (desloc(0) != 1.0e-240)
                {
                    DirichletConditionFE *dir = new DirichletConditionFE(no, 0, desloc(0));
                    dirichletConditionsFE_.push_back(dir);
                }
                if (desloc(1) != 1.0e-240)
                {
                    DirichletConditionFE *dir = new DirichletConditionFE(no, 1, desloc(1));
                    dirichletConditionsFE_.push_back(dir);
                }
            }
        }
    }
}

void GlobalSolid::addBlendingZone(std::vector<Line *> lines, const double &thickness, const int &blendingFunctionOrder, const bool &plotHammerPoints)
{
    for (Line *l : lines)
    {
        std::vector<BoundaryElement *> boundary = finiteElementBoundary_[l->getName()];
        for (BoundaryElement *b : boundary)
        {
            blendingBoundary_.push_back(b);
        }
    }
    linesOfBlendingZone_ = lines;
    blendingZoneThickness_ = thickness;
    overlappingAnalysis_ = true;
    plotHammerPointsInBlendingZone_ = plotHammerPoints;
    blendingFunctionOrder_ = blendingFunctionOrder;
}

void GlobalSolid::addBlendingZoneInLocalMesh(const double &thickness, const int &blendingFunctionOrder, const bool &plotHammerPoints)
{
    blendingZoneThickness_ = thickness;
    plotHammerPointsInBlendingZone_ = plotHammerPoints;
    overlappingAnalysis_ = true;
    localOverlapping_ = true;
    blendingFunctionOrder_ = blendingFunctionOrder;
}

void GlobalSolid::setBlendingZoneLines()
{
    std::unordered_map<std::string, Crack *> crackes = geometry_->getCrackes();
    for (std::unordered_map<std::string, Crack *>::const_iterator crack = crackes.begin(); crack != crackes.end(); crack++)
    {
        std::vector<BoundaryElement *> boundary = finiteElementBoundary_[crack->second->getPlaneSurface()->getLineLoop(0)->getLine(0)->getName()];
        for (BoundaryElement *b : boundary)
        {
            blendingBoundary_.push_back(b);
        }
    }
}

void GlobalSolid::setParametersOfCrackPropagation(const double &length, const double &fractureThougness, const bool &plotSIFs, const bool &viewNewMesh, const bool &exportUndeformedMesh)
{
    crackLengthPropagation_ = length;
    fractureThougness_ = fractureThougness;
    viewNewMesh_ = viewNewMesh;
    exportUndeformedMesh_ = exportUndeformedMesh;
    analysisOfCrackPropagation_ = true;
    plotSIFs_ = plotSIFs;
}

bounded_vector<double, 2> GlobalSolid::getSIFs()
{
    bounded_vector<double, 2> sif;
    sif(0) = 0.0;
    sif(1) = 0.0;

    Node *node = geometry_->getCrackes()["c0"]->getLastPoint()->getPointNode();

    double angle = geometry_->getCrackes()["c0"]->getLastAngle();

    for (Element *el : JintegralElements_)
    {
        sif += el->contributionJ_IntegralInitial(quadrature_, elementsSideJintegral_[el->getIndex()], planeState_, angle, node->getInitialCoordinate());
        // sif += el->contributionJ_IntegralFromRice(quadrature_, elementsSideJintegral_[el->getIndex()], planeState_, angle, node);
    }

    // for (Element *el : dynamicJintegralElements_)
    // {
    //      // sif += el->domainJ_IntegralInitial2(hammerPoints_, planeState_, problemType_, geometry_->getRadiusJintegral(), angle, node->getInitialCoordinate());
    //      //sif += el->domainJ_IntegralInitial(hammerPoints_, planeState_, "STATIC", 0.2, angle, node->getInitialCoordinate());
    // }

    ///CORRELAÇÃO DE DESLOCAMENTO///

    // std::vector<BoundaryElement *> c00 = finiteElementBoundary_["c00"];
    // std::vector<BoundaryElement *> c01 = finiteElementBoundary_["c01"];

    // double parametro = geometry_->getRadiusJintegral();

    // double KIsum = 0.0;
    // double KIIsum = 0.0;
    // double sum = 0.0;

    // bounded_vector<double, 2> coordTip = node->getInitialCoordinate();

    // double mi, pi, k;
    // mi = 199.992e06 / (2.0 * (1.0 + 0.3));
    // k = 27.0 / 13.0;
    // pi = 3.14159265358979323846;

    // for (int ie = 0; ie < c00.size(); ie++)
    // {
    //     std::vector<Node *> nodc00 = c00[ie]->getNodes();
    //     std::vector<Node *> nodc01 = c01[ie]->getNodes();

    //     if (norm_2(coordTip - nodc00[0]->getInitialCoordinate()) <= parametro and norm_2(coordTip - nodc00[1]->getInitialCoordinate()) <= parametro)
    //     {
    //         for (int in = 0; in < nodc00.size(); in++)
    //         {
    //             bounded_vector<double, 2> initial = nodc00[in]->getInitialCoordinate();

    //             if (initial(0) == coordTip(0) and initial(1) == coordTip(1))
    //             {
    //             }
    //             else
    //             {
    //                 bounded_vector<double, 2> delta = nodc00[in]->getCurrentCoordinate() - nodc01[in]->getCurrentCoordinate();
    //                 double radius = norm_2(initial - coordTip);

    //                 double COD = delta(0) * (-sin(angle)) + delta(1) * cos(angle);
    //                 double CSD = delta(0) * cos(angle) + delta(1) * sin(angle);

    //                 KIsum += (mi / (k + 1.0)) * sqrt(2.0 * pi / radius) * COD;
    //                 KIIsum += (mi / (k + 1.0)) * sqrt(2.0 * pi / radius) * CSD;
    //                 sum += 1.0;
    //             }
    //         }
    //     }
    // }

    // std::cout << sum << std::endl;

    // sif(0) = KIsum / sum;
    // sif(1) = KIIsum / sum;

    return sif;
}

void GlobalSolid::verifyCrackPropagation(bool &testPropagation)
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    bounded_vector<double, 2> SIFs = getSIFs();

    double tetap;

    if (SIFs(1) == 0)
    {
        tetap = 0.0;
    }
    else if (SIFs(1) > 0)
    {
        double aux = SIFs(0) / SIFs(1);
        tetap = 2.0 * atan(0.25 * (aux - sqrt(aux * aux + 8.0)));
    }
    else
    {
        double aux = SIFs(0) / SIFs(1);
        tetap = 2.0 * atan(0.25 * (aux + sqrt(aux * aux + 8.0)));
    }

    double SIFeq = cos(0.5 * tetap) * (SIFs(0) * cos(0.5 * tetap) * cos(0.5 * tetap) - 1.5 * SIFs(1) * sin(tetap));

    bounded_vector<double, 2> tipCrack = geometry_->getCrackes()["c0"]->getLastPoint()->getCoordinates();
    double alfa = geometry_->getCrackes()["c0"]->getLastAngle();

    bounded_vector<double, 2> newTipCrack;
    if (SIFeq >= fractureThougness_)
    {
        newTipCrack(0) = tipCrack(0) + crackLengthPropagation_ * cos(tetap + alfa);
        newTipCrack(1) = tipCrack(1) + crackLengthPropagation_ * sin(tetap + alfa);
        geometry_->addCrackPoint("c0", newTipCrack);

        for (int i = 0; i < elements_.size(); i++)
        {
            delete elements_[i];
        }
        for (int i = 0; i < dirichletConditionsFE_.size(); i++)
        {
            delete dirichletConditionsFE_[i];
        }
        for (int i = 0; i < neumannConditionsFE_.size(); i++)
        {
            delete neumannConditionsFE_[i];
        }
        for (int i = 0; i < nodes_.size(); i++)
        {
            delete nodes_[i];
        }
        for (std::unordered_map<std::string, std::vector<BoundaryElement *>>::const_iterator bound = finiteElementBoundary_.begin(); bound != finiteElementBoundary_.end(); bound++)
        {
            std::vector<BoundaryElement *> boundElements = bound->second;
            for (int i = 0; i < boundElements.size(); i++)
            {
                delete boundElements[i];
            }
        }
        for (std::unordered_map<std::string, Mesh *>::const_iterator m = meshes_.begin(); m != meshes_.end(); m++)
        {
            delete m->second;
        }
        for (int i = 0; i < inactiveCP_.size(); i++)
        {
            delete inactiveCP_[i];
        }
        for (int i = 0; i < inactiveNode_.size(); i++)
        {
            delete inactiveNode_[i];
        }
        elements_.clear();
        elements_part.clear();
        nodes_.clear();
        finiteElementBoundary_.clear();
        meshes_.clear();
        dirichletConditionsFE_.clear();
        neumannConditionsFE_.clear();
        JintegralElements_.clear();
        elementsSideJintegral_.clear();
        inactiveCP_.clear();
        inactiveNode_.clear();
        inactiveCPandNode_.clear();
        blendingBoundary_.clear();

        cpnumber_ = cpaux_;

        testPropagation = true;

        generateMesh(geometry_, finiteElementType_, typeAlgorithm_, "newMesh", viewNewMesh_, false);

        if (rank == 0)
        {
            std::cout << "The crack tip propagated and a new mesh have been done." << std::endl;
        }
    }
    else
    {
        testPropagation = false;
    }
}

void GlobalSolid::quarterPointElements()
{
    std::unordered_map<std::string, Crack *> crackes = geometry_->getCrackes();
    for (std::unordered_map<std::string, Crack *>::const_iterator crack = crackes.begin(); crack != crackes.end(); crack++)
    {
        std::string open = crack->second->getOpenBoundary();
        if (open != "second")
        {
            Node *lastnode = crack->second->getLastPoint()->getPointNode();
            if (finiteElementType_ == "T6")
            {
                for (Element *el : elements_)
                {
                    std::vector<Node *> conection = el->getConnection();
                    bounded_vector<double, 2> coord0, coord1, coord2;
                    coord0 = conection[0]->getInitialCoordinate();
                    coord1 = conection[1]->getInitialCoordinate();
                    coord2 = conection[2]->getInitialCoordinate();

                    if (conection[0] == lastnode)
                    {
                        bounded_vector<double, 2> coord3, coord5;
                        coord3 = 0.25 * (3.0 * coord0 + coord1);
                        coord5 = 0.25 * (3.0 * coord0 + coord2);
                        conection[3]->setInitialCoordinate(coord3);
                        conection[5]->setInitialCoordinate(coord5);
                    }
                    else if (conection[1] == lastnode)
                    {

                        bounded_vector<double, 2> coord3, coord4;
                        coord3 = 0.25 * (coord0 + 3.0 * coord1);
                        coord4 = 0.25 * (3.0 * coord1 + coord2);
                        conection[3]->setInitialCoordinate(coord3);
                        conection[4]->setInitialCoordinate(coord4);
                    }
                    else if (conection[2] == lastnode)
                    {
                        bounded_vector<double, 2> coord5, coord4;
                        coord5 = 0.25 * (coord0 + 3.0 * coord2);
                        coord4 = 0.25 * (coord1 + 3.0 * coord2);
                        conection[5]->setInitialCoordinate(coord5);
                        conection[4]->setInitialCoordinate(coord4);
                    }
                }
            }
            else if (finiteElementType_ == "T10")
            {
                for (Element *el : elements_)
                {
                    std::vector<Node *> conection = el->getConnection();
                    bounded_vector<double, 2> coord0, coord1, coord2;
                    coord0 = conection[0]->getInitialCoordinate();
                    coord1 = conection[1]->getInitialCoordinate();
                    coord2 = conection[2]->getInitialCoordinate();

                    if (conection[0] == lastnode)
                    {
                        bounded_vector<double, 2> coord3, coord8, coord4, coord7;
                        // coord3 = (7.0 * coord0 + 2.0 * coord1) / 9.0;
                        // coord8 = (7.0 * coord0 + 2.0 * coord2) / 9.0;
                        // conection[3]->setInitialCoordinate(coord3);
                        // conection[8]->setInitialCoordinate(coord8);

                        // coord3 = (31.0 * coord0 + 5.0 * coord1) / 36.0;
                        // coord8 = (31.0 * coord0 + 5.0 * coord2) / 36.0;
                        // conection[3]->setInitialCoordinate(coord3);
                        // conection[8]->setInitialCoordinate(coord8);

                        // coord4 = (coord0 + coord1) / 2.0;
                        // coord7 = (coord0 + coord2) / 2.0;
                        // conection[4]->setInitialCoordinate(coord4);
                        // conection[7]->setInitialCoordinate(coord7);

                        coord3 = (17.0 * coord0 + coord1) / 18.0;
                        coord8 = (17.0 * coord0 + coord2) / 18.0;
                        conection[3]->setInitialCoordinate(coord3);
                        conection[8]->setInitialCoordinate(coord8);

                        coord4 = (2.0 * coord0 + coord1) / 3.0;
                        coord7 = (2.0 * coord0 + coord2) / 3.0;
                        conection[4]->setInitialCoordinate(coord4);
                        conection[7]->setInitialCoordinate(coord7);

                        conection[9]->setInitialCoordinate((coord4 + coord7) * 0.5);
                    }
                    else if (conection[1] == lastnode)
                    {
                        bounded_vector<double, 2> coord4, coord5, coord3, coord6;
                        // coord4 = (2.0 * coord0 + 7.0 * coord1) / 9.0;
                        // coord5 = (7.0 * coord1 + 2.0 * coord2) / 9.0;
                        // conection[4]->setInitialCoordinate(coord4);
                        // conection[5]->setInitialCoordinate(coord5);

                        // coord4 = (5.0 * coord0 + 31.0 * coord1) / 36.0;
                        // coord5 = (31.0 * coord1 + 5.0 * coord2) / 36.0;
                        // conection[4]->setInitialCoordinate(coord4);
                        // conection[5]->setInitialCoordinate(coord5);

                        // coord3 = (coord0 + coord1) / 2.0;
                        // coord6 = (coord1 + coord2) / 2.0;
                        // conection[3]->setInitialCoordinate(coord3);
                        // conection[6]->setInitialCoordinate(coord6);

                        coord4 = (coord0 + 17.0 * coord1) / 18.0;
                        coord5 = (17.0 * coord1 + coord2) / 18.0;
                        conection[4]->setInitialCoordinate(coord4);
                        conection[5]->setInitialCoordinate(coord5);

                        coord3 = (coord0 + 2.0 * coord1) / 3.0;
                        coord6 = (2.0 * coord1 + coord2) / 3.0;
                        conection[3]->setInitialCoordinate(coord3);
                        conection[6]->setInitialCoordinate(coord6);

                        conection[9]->setInitialCoordinate((coord6 + coord3) * 0.5);
                    }
                    else if (conection[2] == lastnode)
                    {
                        bounded_vector<double, 2> coord6, coord7, coord5, coord8;
                        // coord6 = (2.0 * coord0 + 7.0 * coord2) / 9.0;
                        // coord7 = (2.0 * coord1 + 7.0 * coord2) / 9.0;
                        // conection[6]->setInitialCoordinate(coord6);
                        // conection[7]->setInitialCoordinate(coord7);

                        // coord6 = (5.0 * coord0 + 31.0 * coord2) / 36.0;
                        // coord7 = (5.0 * coord1 + 31.0 * coord2) / 36.0;
                        // conection[6]->setInitialCoordinate(coord6);
                        // conection[7]->setInitialCoordinate(coord7);

                        // coord5 = (coord0 + coord2) / 2.0;
                        // coord8 = (coord1 + coord2) / 2.0;
                        // conection[5]->setInitialCoordinate(coord5);
                        // conection[8]->setInitialCoordinate(coord8);

                        coord6 = (coord0 + 17.0 * coord2) / 18.0;
                        coord7 = (coord1 + 17.0 * coord2) / 18.0;
                        conection[6]->setInitialCoordinate(coord6);
                        conection[7]->setInitialCoordinate(coord7);

                        coord5 = (coord0 + 2.0 * coord2) / 3.0;
                        coord8 = (coord1 + 2.0 * coord2) / 3.0;
                        conection[5]->setInitialCoordinate(coord5);
                        conection[8]->setInitialCoordinate(coord8);

                        conection[9]->setInitialCoordinate((coord5 + coord8) * 0.5);
                    }
                }
            }
        }
        if (open != "first")
        {
            Node *initialnode = crack->second->getFirstPoint()->getPointNode();
            if (finiteElementType_ == "T6")
            {
                for (Element *el : elements_)
                {
                    std::vector<Node *> conection = el->getConnection();
                    bounded_vector<double, 2> coord0, coord1, coord2;
                    coord0 = conection[0]->getInitialCoordinate();
                    coord1 = conection[1]->getInitialCoordinate();
                    coord2 = conection[2]->getInitialCoordinate();

                    if (conection[0] == initialnode)
                    {
                        bounded_vector<double, 2> coord3, coord5;
                        coord3 = 0.25 * (3.0 * coord0 + coord1);
                        coord5 = 0.25 * (3.0 * coord0 + coord2);
                        conection[3]->setInitialCoordinate(coord3);
                        conection[5]->setInitialCoordinate(coord5);
                    }
                    else if (conection[1] == initialnode)
                    {

                        bounded_vector<double, 2> coord3, coord4;
                        coord3 = 0.25 * (coord0 + 3.0 * coord1);
                        coord4 = 0.25 * (3.0 * coord1 + coord2);
                        conection[3]->setInitialCoordinate(coord3);
                        conection[4]->setInitialCoordinate(coord4);
                    }
                    else if (conection[2] == initialnode)
                    {
                        bounded_vector<double, 2> coord5, coord4;
                        coord5 = 0.25 * (coord0 + 3.0 * coord2);
                        coord4 = 0.25 * (coord1 + 3.0 * coord2);
                        conection[5]->setInitialCoordinate(coord5);
                        conection[4]->setInitialCoordinate(coord4);
                    }
                }
            }
            else if (finiteElementType_ == "T10")
            {
                for (Element *el : elements_)
                {
                    std::vector<Node *> conection = el->getConnection();
                    bounded_vector<double, 2> coord0, coord1, coord2;
                    coord0 = conection[0]->getInitialCoordinate();
                    coord1 = conection[1]->getInitialCoordinate();
                    coord2 = conection[2]->getInitialCoordinate();

                    if (conection[0] == initialnode)
                    {
                        bounded_vector<double, 2> coord3, coord8, coord4, coord7;
                        // coord3 = (7.0 * coord0 + 2.0 * coord1) / 9.0;
                        // coord8 = (7.0 * coord0 + 2.0 * coord2) / 9.0;
                        // conection[3]->setInitialCoordinate(coord3);
                        // conection[8]->setInitialCoordinate(coord8);

                        // coord3 = (31.0 * coord0 + 5.0 * coord1) / 36.0;
                        // coord8 = (31.0 * coord0 + 5.0 * coord2) / 36.0;
                        // conection[3]->setInitialCoordinate(coord3);
                        // conection[8]->setInitialCoordinate(coord8);

                        // coord4 = (coord0 + coord1) / 2.0;
                        // coord7 = (coord0 + coord2) / 2.0;
                        // conection[4]->setInitialCoordinate(coord4);
                        // conection[7]->setInitialCoordinate(coord7);

                        coord3 = (17.0 * coord0 + coord1) / 18.0;
                        coord8 = (17.0 * coord0 + coord2) / 18.0;
                        conection[3]->setInitialCoordinate(coord3);
                        conection[8]->setInitialCoordinate(coord8);

                        coord4 = (2.0 * coord0 + coord1) / 3.0;
                        coord7 = (2.0 * coord0 + coord2) / 3.0;
                        conection[4]->setInitialCoordinate(coord4);
                        conection[7]->setInitialCoordinate(coord7);

                        conection[9]->setInitialCoordinate((coord4 + coord7) * 0.5);
                    }
                    else if (conection[1] == initialnode)
                    {
                        bounded_vector<double, 2> coord4, coord5, coord3, coord6;
                        // coord4 = (2.0 * coord0 + 7.0 * coord1) / 9.0;
                        // coord5 = (7.0 * coord1 + 2.0 * coord2) / 9.0;
                        // conection[4]->setInitialCoordinate(coord4);
                        // conection[5]->setInitialCoordinate(coord5);

                        // coord4 = (5.0 * coord0 + 31.0 * coord1) / 36.0;
                        // coord5 = (31.0 * coord1 + 5.0 * coord2) / 36.0;
                        // conection[4]->setInitialCoordinate(coord4);
                        // conection[5]->setInitialCoordinate(coord5);

                        // coord3 = (coord0 + coord1) / 2.0;
                        // coord6 = (coord1 + coord2) / 2.0;
                        // conection[3]->setInitialCoordinate(coord3);
                        // conection[6]->setInitialCoordinate(coord6);

                        coord4 = (coord0 + 17.0 * coord1) / 18.0;
                        coord5 = (17.0 * coord1 + coord2) / 18.0;
                        conection[4]->setInitialCoordinate(coord4);
                        conection[5]->setInitialCoordinate(coord5);

                        coord3 = (coord0 + 2.0 * coord1) / 3.0;
                        coord6 = (2.0 * coord1 + coord2) / 3.0;
                        conection[3]->setInitialCoordinate(coord3);
                        conection[6]->setInitialCoordinate(coord6);

                        conection[9]->setInitialCoordinate((coord6 + coord3) * 0.5);
                    }
                    else if (conection[2] == initialnode)
                    {
                        bounded_vector<double, 2> coord6, coord7, coord5, coord8;
                        // coord6 = (2.0 * coord0 + 7.0 * coord2) / 9.0;
                        // coord7 = (2.0 * coord1 + 7.0 * coord2) / 9.0;
                        // conection[6]->setInitialCoordinate(coord6);
                        // conection[7]->setInitialCoordinate(coord7);

                        // coord6 = (5.0 * coord0 + 31.0 * coord2) / 36.0;
                        // coord7 = (5.0 * coord1 + 31.0 * coord2) / 36.0;
                        // conection[6]->setInitialCoordinate(coord6);
                        // conection[7]->setInitialCoordinate(coord7);

                        // coord5 = (coord0 + coord2) / 2.0;
                        // coord8 = (coord1 + coord2) / 2.0;
                        // conection[5]->setInitialCoordinate(coord5);
                        // conection[8]->setInitialCoordinate(coord8);

                        coord6 = (coord0 + 17.0 * coord2) / 18.0;
                        coord7 = (coord1 + 17.0 * coord2) / 18.0;
                        conection[6]->setInitialCoordinate(coord6);
                        conection[7]->setInitialCoordinate(coord7);

                        coord5 = (coord0 + 2.0 * coord2) / 3.0;
                        coord8 = (coord1 + 2.0 * coord2) / 3.0;
                        conection[5]->setInitialCoordinate(coord5);
                        conection[8]->setInitialCoordinate(coord8);

                        conection[9]->setInitialCoordinate((coord5 + coord8) * 0.5);
                    }
                }
            }
        }
    }
}

void GlobalSolid::useQuarterPointElements()
{
    quarterPointElement_ = true;
}

void GlobalSolid::error()
{
    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    bounded_vector<double, 2> normsError;
    normsError(0) = 0.0;
    normsError(1) = 0.0;
    for (Element *el : elements_part)
    {
        normsError += el->errorL2(hammerPoints_, hammerPointsBlendZone_);
    }
    std::cout << rank << std::scientific << " L2: " << normsError(0) << " - H1: " << normsError(1) << std::endl;
    normsError(0) = 0.0;
    normsError(1) = 0.0;
    for (Cell *el : cells_part)
    {
        normsError += el->errorL2(quadrature_);
    }
    std::cout << rank << std::scientific << " L2: " << normsError(0) << " - H1: " << normsError(1) << std::endl;
}