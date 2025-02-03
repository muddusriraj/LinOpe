#ifndef GRAPH_H
#define GRAPH_H

#include <GL/freeglut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <vector>
#include "operators.h"


void drawArrow(double x, double y, double angle);

void drawText(const char* text, double x, double y);

void drawVector(const std::vector<double>& vec, const std::string color);

void plot(const std::vector<darray>& vectors, const std::vector<std::string> colors = { "blue" });


#endif
