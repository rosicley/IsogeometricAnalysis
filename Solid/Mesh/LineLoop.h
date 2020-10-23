#pragma once

#include "Line.h"

class LineLoop
{
public:
    LineLoop();

    LineLoop(const std::string &name, std::vector<Line *> lines);

    ~LineLoop();

    std::string getName();

    Line *getLine(const int &index);

    std::vector<Line *> getLines();

    std::string getGmshCode();

    void verification();

private:
    std::string name_;
    std::vector<Line *> lines_;
};
