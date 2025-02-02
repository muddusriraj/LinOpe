#include "graph.h"

std::vector<std::vector<double>> g_vectors;
std::vector<std::string> g_colors;

void drawArrow(double x, double y, double angle) {
    double arrowSize = 0.5f;
    double angleRad = angle * 3.14159f / 180.0f;

    double x1 = x - arrowSize * cos(angleRad - 0.3f);
    double y1 = y - arrowSize * sin(angleRad - 0.3f);

    double x2 = x - arrowSize * cos(angleRad + 0.3f); 
    double y2 = y - arrowSize * sin(angleRad + 0.3f);

    glBegin(GL_LINES);
        glVertex2f(x, y);
        glVertex2f(x - arrowSize * cos(angleRad), y - arrowSize * sin(angleRad));
        glVertex2f(x, y);
        glVertex2f(x1, y1);
        glVertex2f(x, y);
        glVertex2f(x2, y2);
    glEnd();
}

void drawText(const char* text, double x, double y) {
    glRasterPos2f(x, y);
    for (size_t i = 0; i < strlen(text); i++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
    }
}

void drawVector(const std::vector<double>& vec, const std::string color) {
    double vx = vec[0];
    double vy = vec[1];
    double c1=0.0;
    double c2=0.0;
    double c3=0.0;
    if(color == "red"){
        c1=1.0;
    } else if (color == "green"){
        c2=1.0;
    } else if (color == "blue"){
        c3=1.0;
    }
    glColor3f(c1, c2, c3);
    glBegin(GL_LINES);
        glVertex2d(0.0, 0.0);
        glVertex2d(vx, vy);
    glEnd();

    double angle = atan2(vy, vx) * 180.0 / 3.14159;
    drawArrow(vx, vy, angle);
}

void plot(const std::vector<darray>& vectors, const std::vector<std::string> colors) {
    std::vector<std::vector<double>> tempVecs(vectors.size());
    for(size_t i=0; i<vectors.size(); i++){
        tempVecs[i] = vectors[i].getVector();
    }
    g_vectors = tempVecs;
    if(colors.size() != tempVecs.size()) {
        g_colors.resize(tempVecs.size());
        for(size_t i=0; i<g_colors.size(); i++) {
            g_colors[i] = "blue";
        }
    } else {
        g_colors = colors;
    }

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 600);
    glutCreateWindow("2D Grid Plane");

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); 
    gluOrtho2D(-10.0, 10.0, -10.0, 10.0);  

    glutDisplayFunc([]() {
        glClear(GL_COLOR_BUFFER_BIT); 

        glColor3f(0.8f, 0.8f, 0.8f);  
        for (float i = -9.0f; i <= 9.0f; i++) {
            glBegin(GL_LINES);
                glVertex2f(i, -10.0f); glVertex2f(i, 10.0f); 
                glVertex2f(-10.0f, i); glVertex2f(10.0f, i); 
            glEnd();
        }

        glColor3f(0.0f, 0.0f, 0.0f);
        glBegin(GL_LINES);
            glVertex2f(-10.0f, 0.0f);
            glVertex2f(10.0f, 0.0f); 
            glVertex2f(0.0f, -10.0f); 
            glVertex2f(0.0f, 10.0f);   
        glEnd();

        drawArrow(10.0f, 0.0f, 0.0f); 
        drawArrow(-10.0f, 0.0f, 180.0f);
        drawArrow(0.0f, 10.0f, 90.0f); 
        drawArrow(0.0f, -10.0f, 270.0f);

        for (size_t i = 0; i < g_vectors.size(); i++) {
            drawVector(g_vectors[i], g_colors[i]);
        }

        glColor3f(0.0f, 0.0f, 0.0f); 
        drawText("X", 10.5f, 0.0f);
        drawText("Y", 0.0f, 10.5f); 

        glutSwapBuffers(); 
    });

    glutMainLoop(); 
}

