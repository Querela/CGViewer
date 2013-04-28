#include <QApplication>
#include "main.h"

int main(int argc, char *argv[])
{
     QApplication app(argc, argv);

     GUI window;
     window.resize(1024,768);

     window.show();
     return app.exec();
} 
