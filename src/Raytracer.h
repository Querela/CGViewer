#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <QtGui>
#include <GL/glew.h>
#include <GL/glext.h>
#include <QGLWidget>
#include <vector>
#include <cfloat>

#include "CGMath.h"
#include "Octree.h"

#define MAX_DEPTH 4
#define BOX_MARGIN 0.1f

// undefine/comment out if you want no textures
#define TEXTURE_ON
// undefine/comment out if you want no reflections (faster raytracing w/o)
#define REFLECTION_ON
// undefine/comment out if you want no refractions (faster raytracing w/o)
#define REFRACTION_ON
// undefine/comment out if you want no shadows
#define SHADOWS_ON
// uncomment if you want only the shadow part of the image
//#define INVERT_SHADOWS

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
            void close();
            QColor raytrace(Vector r, Vector e, int depth);
            QColor raytrace(Vector start, Vector dir, int depth, float density);

            QImage *image;
            QImage finalImage;

            // scene data
            std::vector<Triangle> triangles;
            std::vector<Lightsource> lights;
            Octree *octree;

            //Scene setup
            Vector camera, center, upVector;
            QColor backgroundColor, ambientLight;

            float focalLength;
            float superSamplingRate;

            //acceleration for background
            char** renderedImage;
            unsigned int *idx;
            GLuint displayList;
            GLuint screenTexID;
            GLuint fboId;
};

#endif
