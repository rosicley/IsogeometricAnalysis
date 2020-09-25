#pragma once

#include "Point.h"
#include "../FiniteElement/Node.h"

class Line
{
public:
    Line();

    Line(const std::string &name, std::vector<Point *> points, const std::string &type, const bool& discretization); //, const bool &discretization = true);

    ~Line();

    Line *operator-();

    std::string getName();

    Point *getInitialPoint();

    Point *getEndPoint();

    std::string getGmshCode();

    // std::string getGmshCodeCircle();

    // std::string getGmshCodeEllipse();

    // std::string getGmshCodeSpline();

    void appendNodes(std::vector<Node *> nodes);

    void setName(const std::string &name);

private:
    std::string name_;
    std::vector<Point *> points_;
    bool discretization_;
    std::vector<Node *> lineNodes_;
    std::string typeLine_; //typeOfLine line, circle, bpsline, spline 
};
