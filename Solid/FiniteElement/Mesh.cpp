#include "Mesh.h"

Mesh::Mesh() {}

Mesh::Mesh(const int &index, Material *mat, const double &thickness, const std::string &elementType)
{
    index_ = index;
    //nnodes_ = nnodes;
    //nelem_ = nelem;
    material_ = mat;
    thickness_ = thickness;
    elementType_ = elementType;
}

Mesh::~Mesh() {}

double Mesh::getThickness()
{
    return thickness_;
}

Material *Mesh::getMaterial()
{
    return material_;
}

int Mesh::getIndex()
{
    return index_;
}

// void Mesh::addNode(const int &index, const int &indexFE, const bounded_vector<double, 2> &initialCoordinate)
// {
//     Node *node = new Node(index, indexFE, initialCoordinate);
//     nodes_.push_back(node);
// }

// void Mesh::removeNodes()
// {
//     nodes_.erase(nodes_.begin(), nodes_.begin() + nodes_.size());
// }

std::string Mesh::getElementType()
{
    return elementType_;
}

// Node *Mesh::getNode(int index)
// {
//     return nodes_[index];
// }

// std::vector<Node *> Mesh::getNodes()
// {
//     return nodes_;
// }

//void Mesh::addBoundaryElement(const int &boundaryIndex, const std::vector<int> &connection);

