#include "LineLoop.h"

LineLoop::LineLoop() {}

LineLoop::LineLoop(const std::string &name, std::vector<Line *> lines)
{
    name_ = name;
    lines_.reserve(lines.size());
    for (Line *line : lines)
    {
        lines_.push_back(line);
    }
}

LineLoop::~LineLoop() {}

std::string LineLoop::getName()
{
    return name_;
}

// Line *LineLoop::getLine(const int &index)
// {
//     return lines_[index];
// }

std::vector<Line *> LineLoop::getLines()
{
    return lines_;
}

std::string LineLoop::getGmshCode()
{
    std::stringstream text;
    text << name_ << " = newll; Line Loop(" << name_ << ") = {";
    for (size_t i = 0; i < lines_.size(); i++)
    {
        text << lines_[i]->getName();
        if (i != (lines_.size() - 1))
            text << ", ";
    }
    text << "};\n//\n";
    return text.str();
}

void LineLoop::verification()
{
    for (size_t i = 0; i < lines_.size(); i++)
    {
        std::string name = lines_[i]->getName();
        int index = i;
        std::string previous_name = (i == 0) ? lines_[lines_.size() - 1]->getName() : lines_[i - 1]->getName();
        int previous_index = (i == 0) ? lines_.size() - 1 : i - 1;
        Point *initial_point = (name[0] == '-') ? lines_[index]->getEndPoint() : lines_[index]->getInitialPoint();
        Point *end_point = (previous_name[0] == '-') ? lines_[previous_index]->getInitialPoint() : lines_[previous_index]->getEndPoint();
        if (initial_point->getName() != end_point->getName())
        {
            std::stringstream text;
            text << "The lines ";
            for (size_t j = 0; j < lines_.size(); j++)
            {
                text << lines_[j]->getName();
                if (j != (lines_.size() - 1))
                    text << ", ";
            }
            text << " do not form a closed loop." << std::endl;
            std::cout << text.str();
            //break;
        }
    }
}