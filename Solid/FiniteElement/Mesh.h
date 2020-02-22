#pragma once

#include "../Material.h"
#include "Node.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <string>

class Mesh
{
public:
    Mesh();

    Mesh(const int &index, const int &nnodes, const int &nelem, Material *mat, const double &thickness, const std::string &elementType);

    ~Mesh();

    double getThickness();

    Material *getMaterial();

    void addNode(const int &index, const int &indexFE, const bounded_vector<double, 2> &initialCoordinate);

    Node *getNode(int index);

    void removeNodes();

    std::string getElementType();

    std::vector<Node *> getNodes();

private:
    int index_;

    int nnodes_;
    
    int nelem_;

    Material *material_;

    double thickness_;

    std::vector<Node *> nodes_;

    std::string elementType_;
};