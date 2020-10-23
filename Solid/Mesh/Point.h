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
          const bounded_vector<double, 2> &coordinates,
          const double &lcar,
          const bool &discretization);
    ~Point();

    std::string getName();

    std::string getGmshCode();

    void addNodeToPoint(Node *node);

    bounded_vector<double, 2> getCoordinates();

    void setCoordinates(const bounded_vector<double, 2> &newCoordinates);

    double getlcar();

    Node *getPointNode();

    void setlcar(const double &newlcar);

    // void setCrackPoint();

    // bool getCrackPoint();

private:
    std::string name_;                      // Gmsh Physical entity name
    bounded_vector<double, 2> coordinates_; // Coordinates vector (x,y)
    double lcar_;                           // Characteristic length Gmsh parameter
    bool discretization_;                   // Choose to discretize a point with a mesh node
    Node *pointNode_;                       // Defines the mesh node discretizing the point
  //  bool crackPoint_ = false;
};