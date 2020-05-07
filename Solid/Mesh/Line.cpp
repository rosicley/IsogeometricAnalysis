#include "Line.h"

Line::Line() {}

Line::Line(const std::string &name, std::vector<Point *> points) //, const bool &discretization)
{
    name_ = name;
    //discretization_ = discretization;
    points_.reserve(points.size());

    points_ = points;
    // for (Point *point : points)
    //     points_.push_back(point);
}

Line::~Line() {}

Line *Line::operator-()
{
    Line *copy = new Line(name_, {points_});
    copy->setName("-" + name_);
    return copy;
}

std::string Line::getName()
{
    return name_;
}
Point *Line::getInitialPoint()
{
    return points_[0];
}
Point *Line::getEndPoint()
{
    return points_[points_.size()-1];
}

std::string Line::getGmshCodeLine()
{
    std::stringstream text;
    // if (discretization_)
    // {
    text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName()
         << "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
    return text.str();
    //}
    // else
    // {
    //     text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName()
    //          << "};\n//\n";
    //     return text.str();
    // }
}

std::string Line::getGmshCodeCircle()
{
    std::stringstream text;
    // if (discretization_)
    // {
    text << name_ << " = newl; Circle(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName()
         << ", " << points_[2]->getName() << "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
    return text.str();
    // }
    // else
    // {
    text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName()
         << ", " << points_[2]->getName() << "};\n//\n";
    return text.str();
    //}
}

std::string Line::getGmshCodeEllipse()
{
    std::stringstream text;
    // if (discretization_)
    // {
    text << name_ << " = newl; Circle(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName()
         << ", " << points_[2]->getName() << ", " << points_[3]->getName() << "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
    return text.str();
    // }
    // else
    // {
    //     text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName()
    //          << ", " << points_[2]->getName() << ", " << points_[3]->getName() << "};\n//\n";
    //     return text.str();
    // }
}

void Line::setName(const std::string &name)
{
    name_ = name;
}

void Line::appendNodes(std::vector<Node *> nodes)
{
    for (Node *node1 : nodes)
    {
        bool notDuplicate = true;
        for (Node *node2 : lineNodes_)
        {
            if (node1->getIndex() == node2->getIndex())
            {
                notDuplicate = false;
                break;
            }
        }
        if (notDuplicate)
            lineNodes_.push_back(node1);
    }
}
