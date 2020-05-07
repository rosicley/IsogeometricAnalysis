#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "../FiniteElement/Node.h"

class Point
{
public:
    Point();
    Point(const std::string &name,
          const std::vector<double> &coordinates,
          const double &lcar,
          const bool &discretization);
    ~Point();

    std::string getName();

    std::string getGmshCode();

    void addNodeToPoint(Node *node);

    Node *getPointNode();

private:
    std::string name_;                // Gmsh Physical entity name
    std::vector<double> coordinates_; // Coordinates vector (x,y)
    double lcar_;                     // Characteristic length Gmsh parameter
    bool discretization_;             // Choose to discretize a point with a mesh node
    Node *pointNode_;                 // Defines the mesh node discretizing the point
};