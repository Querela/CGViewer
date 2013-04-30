#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <QtGui>
#include <GL/glew.h>
#include <GL/glext.h>
#include <QGLWidget>
#include <vector>

#include "CGMath.h"

#define MAX_DEPTH 4

class Raytracer : public QGLWidget
{
    Q_OBJECT
        public:
            Raytracer( QString file );
            void genImage();

        protected:

            void paintGL();
            void resizeGL( int width, int height );
            void initializeGL();
            void keyPressEvent(QKeyEvent *event);

        private:
            void init();
            bool initShaders();
            QColor raytrace(Vector r, Vector e, int depth);
            QImage *image;
            QImage finalImage;
            std::vector<Triangle> triangles;

            //Scene setup
            Vector camera, center, upVector;
            float focalLength;
            QColor backgroundColor, ambientLight;
            
            std::vector<Lightsource> lights;

            float superSamplingRate;

            //acceleration
            char** renderedImage;
            unsigned int *idx;
            GLuint displayList;
            GLuint screenTexID;
            GLuint fboId;
            
};

#endif
