TEMPLATE        = app
TARGET          = CGViewer

CONFIG          += qt warn_on debug
QT              += opengl

CONFIG          += console

QMAKE_CXXFLAGS  += -fopenmp

LIBS            += -lGLEW -lGLU -fopenmp

OBJECTS_DIR     = ./obj 
MOC_DIR         = ./moc
VPATH           += ./src

HEADERS         = main.h GUI.h CGMath.h Matrix.h Model.h Scene.h EditWidgets.h Light.h SaveSceneDialog.h Raytracer.h ShaderStuff.h Voxel.h Octree.h

SOURCES         = main.cpp GUI.cpp Matrix.cpp Model.cpp Scene.cpp EditWidgets.cpp Light.cpp SaveSceneDialog.cpp Raytracer.cpp ShaderStuff.cpp Voxel.cpp Octree.cpp

 
