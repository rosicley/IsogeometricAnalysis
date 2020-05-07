#pragma once

#include "Point.h"
#include "../FiniteElement/Node.h"

class Line
{
public:
    Line();

    Line(const std::string &name, std::vector<Point *> points); //, const bool &discretization = true);

    ~Line();

    Line *operator-();

    std::string getName();

    Point *getInitialPoint();

    Point *getEndPoint();

    std::string getGmshCodeLine();

    std::string getGmshCodeCircle();

    std::string getGmshCodeEllipse();

    void appendNodes(std::vector<Node *> nodes);

    void setName(const std::string &name);

private:
    std::string name_;
    std::vector<Point *> points_;
    //bool discretization_;
    std::vector<Node *> lineNodes_;
};
