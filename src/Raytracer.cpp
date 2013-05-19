#include "Raytracer.h"
#include <fstream>
#include <omp.h>

#include <iostream>
using namespace std;

Raytracer::Raytracer( QString path )
{
    superSamplingRate = 1.1; //TODO: setting this to 1 and enabling on the fly output will cause a crash
    //read file and init everything

    vector<Vector> vectors;
    vector<Vector> normals;
    vector<Material> materials;
    vector<Vector> tex_coords;

    int nLights = 0, nMat =0, nVert=0, nNorm=0, nTexCoords=0, nFaces=0;

    //QString path("scene.tri");
    //create File
    ifstream file(path.toUtf8().data());

    file>>nLights>>nMat>>nVert>>nNorm>>nTexCoords>>nFaces;    //read Header

    //0.5 read Raytracer-Properties
    int width,height;
    float tx,ty,tz;    
    file>>width;file>>height;
    file>>tx; file>>ty; file>>tz;
    camera = Vector(tx,ty,tz);
    file>>tx; file>>ty;file>>tz;
    center = Vector(tx,ty,tz);
    file>>focalLength;
    float ux, uy, uz;
    file>>ux>>uy>>uz;
    upVector = Vector(ux,uy,uz);
    //file>>rotationAngle;
    float back_r,back_g,back_b;
    file>>back_r>>back_g>>back_b;
    backgroundColor.setRgb(back_r*255,back_g*255,back_b*255);
    float amb_r,amb_g,amb_b;
    file>>amb_r>>amb_g>>amb_b;
    ambientLight.setRgb(amb_r*255,amb_g*255,amb_b*255);

    image = new QImage(width*superSamplingRate, height*superSamplingRate, QImage::Format_RGB32);
    image->fill(qRgb(255,255,255));    

    //0.75 read the lights
    for (int i=0;i<nLights;++i)
    {
        Lightsource l;        
        file>>tx>>ty>>tz;
        l.position = Vector(tx,ty,tz);
        file>>l.ambient[0]>>l.ambient[1]>>l.ambient[2];
        file>>l.diffuse[0]>>l.diffuse[1]>>l.diffuse[2];
        file>>l.specular[0]>>l.specular[1]>>l.specular[2];
        file>>l.constAtt>>l.linAtt>>l.quadAtt;        
        lights.push_back(l);        
    }

    //1. Materials->Meshes
    for (int i=0;i<nMat;++i)
    {
        Material mat;
        float r,g,b,shininess,alpha,sharpness,density;
        int isTexture;
        file>>r>>g>>b;
        mat.ambient[0]=r; mat.ambient[1]=g; mat.ambient[2]=b; mat.ambient[3]=1.0;
        file>>r>>g>>b;
        mat.diffuse[0]=r; mat.diffuse[1]=g; mat.diffuse[2]=b; mat.diffuse[3]=1.0;
        file>>r>>g>>b;
        mat.specular[0]=r; mat.specular[1]=g; mat.specular[2]=b; mat.specular[3]=1.0;
        file>>shininess;
        mat.shininess = shininess;
        file>>alpha;
        mat.alpha=alpha;
        file>>sharpness;
        mat.sharpness=sharpness;
        file>>density;
        mat.density=density;
        mat.isTexture=false;
        mat.hasNormalMap = false;
        file>>isTexture;        
        if (isTexture > 0) //load texture
        {
            string texName;
            file>>texName;
            cout<<"Texture: "<<texName<<endl;
            QString filepath; //path in were the
            filepath = path;
            filepath.chop(path.length()-path.lastIndexOf("/"));
            filepath.append(QString("/%1").arg(QString(texName.c_str())));
            QImage texture(filepath);
            if(!texture.isNull())
            {
                mat.texture = texture;
                mat.isTexture = true;
            }
            else
                QMessageBox::warning(this, tr("Texture Loader"), tr("Texture couldn't be loaded!"),QMessageBox::Ok);            
        }
        if (isTexture == 2) //load normal map (heavy code duplicates TODO: change that)
        {
            string bumpName;
            file>>bumpName;
            cout<<"Normal Map: "<<bumpName<<endl;
            QString filepath; //path in were the
            filepath = path;
            filepath.chop(path.length()-path.lastIndexOf("/"));
            filepath.append(QString("/%1").arg(QString(bumpName.c_str())));
            QImage bumpMap(filepath);
            if(!bumpMap.isNull())
            {
                mat.hasNormalMap = true;
                mat.normalMap = bumpMap;
            }
            else
                QMessageBox::warning(this, tr("Texture Loader"), tr("Normal map couldn't be loaded!"),QMessageBox::Ok);
        }
        materials.push_back(mat);
    }
    //2.Fill vertices
    for (int i=0;i<nVert;++i)
    {
        float x,y,z;
        file>>x>>y>>z;
        vectors.push_back(Vector(x,y,z));
    }
    //3.Fill normals
    for (int i=0;i<nNorm;++i)
    {
        float x,y,z;
        file>>x>>y>>z;
        Vector n(x,y,z);
        n.normalize();
        normals.push_back(n);
    }
    //3.Fill texCoords
    for (int i=0;i<nTexCoords;++i)
    {
        float u,v;
        file>>u>>v;
        tex_coords.push_back(Vector(u,v,1.0));
    }
    //4.Fill faces /Triangles
    for (int i=0;i<nFaces;++i)
    {
        int matNr,vertNr1,vertNr2,vertNr3,NormNr1,NormNr2,NormNr3,TexCoordsNr1,TexCoordsNr2,TexCoordsNr3;
        file>>matNr>>vertNr1>>vertNr2>>vertNr3>>NormNr1>>NormNr2>>NormNr3>>TexCoordsNr1>>TexCoordsNr2>>TexCoordsNr3;
        Triangle t;
        //cout<<"Size: "<<materials.size()<<", Nr: "<<matNr<<endl;cout.flush();
        t.material = materials[matNr];
        t.vertices[0]=vectors[vertNr1];
        t.vertices[1]=vectors[vertNr2];
        t.vertices[2]=vectors[vertNr3];
        t.normals[0]=normals[NormNr1];
        t.normals[1]=normals[NormNr2];
        t.normals[2]=normals[NormNr3];
        if (materials[matNr].isTexture)
        {
            t.texCoords[0]=tex_coords[TexCoordsNr1];
            t.texCoords[1]=tex_coords[TexCoordsNr2];
            t.texCoords[2]=tex_coords[TexCoordsNr3];
        }
        // computing planeNormals
        t.planeNormal = crossProduct(t.vertices[1] - t.vertices[0],
                                      t.vertices[2] - t.vertices[0]);
        triangles.push_back(t);
    }
    cout<<"Got "<<triangles.size()<<" Triangles\n";

}

void Raytracer::init()
{
    //draw the geometry into a texture and save the primary-Rays
    GLubyte *renderedImage = new GLubyte[image->height()*image->width()*4];
    
    gluLookAt(camera[0], camera[1], camera[2], center[0], center[1], center[2], upVector[0], upVector[1], upVector[2]);
    setGeometry(0,0,1024,768);
    makeCurrent();
    GLenum err = glewInit();
    if (GLEW_OK != err)
        cout << "Error: glewInit() failed\n";
    else
        cout << "Succesfully initiated GLEW\n";

    
    //create DisplayList
    displayList = glGenLists(1);
    glNewList(displayList, GL_COMPILE);
    for (unsigned int i=0; i<triangles.size(); ++i)
    {
        float r = (float)(i % 255);
        float g = (float)((i/255) % 255);
        float b = (float)(((i/255)/255) % 255);
        float a = (float)(((((i/255)/255))/255) % 255);    
        glColor4f(r/255.0,g/255.0,b/255.0,a/255.0);
        glBegin(GL_TRIANGLES);
        glVertex3fv( triangles[i].vertices[0].getValues() );
        glVertex3fv( triangles[i].vertices[1].getValues() );
        glVertex3fv( triangles[i].vertices[2].getValues() );
        glEnd();
    }
    glEndList();

    
    glGenTextures(1, &screenTexID);
    glBindTexture(GL_TEXTURE_2D,screenTexID);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, 4, image->width(), image->height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

    // create a renderbuffer object to store depth info
    GLuint rboId;
    glGenRenderbuffersEXT(1, &rboId);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rboId);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, image->width(), image->height());
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
    
    //GLuint fboId;
    glGenFramebuffersEXT(1, &fboId);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fboId);    

    // attach the texture to FBO color attachment point
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, screenTexID, 0);
    // attach the renderbuffer to depth attachment point
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, rboId);
    
    resizeGL(image->width(), image->height());
    
    paintGL();    
    
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    
    glBindTexture(GL_TEXTURE_2D, screenTexID);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, (void*)renderedImage);

    
    idx = new unsigned int[image->height()*image->width()];

    for (int i=0; i<image->height()*image->width()*4; i+=4)
        idx[i/4] = (int)renderedImage[i] + (int)(renderedImage[i+1])*255 + (int)(renderedImage[i+2])*255*255 + (int)(renderedImage[i+3])*255*255*255;

    delete renderedImage;

}

void Raytracer::resizeGL(int w, int h)
{
    // Reset the viewport
    glViewport(0, 0, w, h);
    // Reset the projection and modelview matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (GLfloat)w/(GLfloat)h, 1.0f, 200.0f);
    glMatrixMode(GL_MODELVIEW);
}

void Raytracer::initializeGL()
{
    glClearColor((float)backgroundColor.red()/255.0, (float)backgroundColor.green()/255.0, (float)backgroundColor.blue()/255.0, 1.0);
    //glClearColor(1.0, 1.0, 1.0, 0.0);
    cout<<"Back_Color: "<<backgroundColor.red()/255.0<<", "<<(float)backgroundColor.green()/255.0<<", "<<(float)backgroundColor.blue()/255.0<<endl;
    
    glDepthFunc(GL_LEQUAL);                            // Type Of Depth Testing
    glEnable(GL_DEPTH_TEST);                        // Enable Depth Testing    
}


void Raytracer::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);    
    glDisable(GL_CULL_FACE);
    glEnable(GL_COLOR_MATERIAL);
    glCallList(displayList);
}


void Raytracer::generateVoxels()
{
    QTime t;
    t.start();

    // get some initial values
    float minX = triangles[0].vertices[0][0];
    float minY = triangles[0].vertices[0][1];
    float minZ = triangles[0].vertices[0][2];
    float maxX = triangles[0].vertices[0][0];
    float maxY = triangles[0].vertices[0][1];
    float maxZ = triangles[0].vertices[0][2];

    // get min and max values -> bounding box of entire scene
    for (unsigned int i = 0; i < triangles.size(); i ++)
    {
        Triangle t = triangles[i];
        for (unsigned int j = 0; j < 3; j ++)
        {
           if (minX > t.vertices[j][0]) minX = t.vertices[j][0];
           if (maxX < t.vertices[j][0]) maxX = t.vertices[j][0];
           if (minY > t.vertices[j][1]) minY = t.vertices[j][1];
           if (maxY < t.vertices[j][1]) maxY = t.vertices[j][1];
           if (minZ > t.vertices[j][2]) minZ = t.vertices[j][2];
           if (maxZ < t.vertices[j][2]) maxZ = t.vertices[j][2];
        }
    }

    cout << "Time (bounding box): " << t.restart() << " msec" << endl;
    // DEBUG: Output bounding box dimensions/coordinates
    cout << "BoundingBox X[" << minX << "/" << maxX << "], Y[" << minY << "/"
         << maxY << "], Z[" << minZ << "/" << maxZ << "]" << endl;

    // ------------------------------------------------------------------------

    // Set origin and dimension of voxels
    voxels.origin = Vector(minX, minY, minZ);
    voxels.size = Vector((maxX - minX), (maxY - minY), (maxZ - minZ));

    // DEBUG: Output voxelGroup
    cout << "VoxelGroup origin: (" << minX << "/" << minY << "/" << minZ << ")" << endl;
    cout << "VoxelGroup size  : (" << voxels.size[0] << "/" << voxels.size[1] << "/"
         << voxels.size[2] << ")" << endl;

    // TEST: compute grid resolution
    float cubeRoot = powf(VOXEL_LAMBDA * triangles.size() /
                     ((maxX - minX) * (maxY - minY) * (maxZ - minZ)), 1 / 3.f);
    unsigned int voxelNumX = std::max(1, std::min((int) ((maxX - minX) * cubeRoot), 128));
    unsigned int voxelNumY = std::max(1, std::min((int) ((maxY - minY) * cubeRoot), 128));
    unsigned int voxelNumZ = std::max(1, std::min((int) ((maxZ - minZ) * cubeRoot), 128));

    // DEBUG: Output voxel grid size
    cout << "Voxel grid size (lambda = " << VOXEL_LAMBDA
                            << "): X = " << voxelNumX
                             << ", Y = " << voxelNumY
                             << ", Z = " << voxelNumZ << endl;

    // Set to const. value
    voxelNumX = VOXEL_NUM_PER_DIM;
    voxelNumY = VOXEL_NUM_PER_DIM;
    voxelNumZ = VOXEL_NUM_PER_DIM;

    // store for later use
    voxels.resolution[0] = voxelNumX;
    voxels.resolution[1] = voxelNumY;
    voxels.resolution[2] = voxelNumZ;

    voxels.voxelSize = Vector((maxX - minX) / voxelNumX,
                              (maxY - minY) / voxelNumY,
                              (maxZ - minZ) / voxelNumZ);
    // DEBUG: Output single voxel size
    //cout << "Voxel-Vector: (" << voxels.voxelSize[0] << "/"
    //                          << voxels.voxelSize[1] << "/"
    //                          << voxels.voxelSize[2] << ")" << endl;

    // generate voxels (assign triangles later)
    for (unsigned int z = 0; z < voxelNumZ; z ++)
    {
        for (unsigned int y = 0; y < voxelNumY; y ++)
        {
            for (unsigned int x = 0; x < voxelNumX; x ++)
            {
                Voxel v;
                v.pos = Vector(minX + x * voxels.voxelSize[0],
                               minY + y * voxels.voxelSize[1],
                               minZ + z * voxels.voxelSize[2]);
                voxels.voxels.push_back(v);
                // DEBUG: Output generated voxels (nr) ?
                //cout << (z * voxelNumX * voxelNumY + y * voxelNumX + x)
                //     << ": (" << v.pos[0] << "/" << v.pos[1] << "/" << v.pos[2] << ")" << endl;
            }
        }
    }

    // insert triangles into voxels
    for (unsigned int i = 0; i < triangles.size(); i ++)
    {
        Triangle t = triangles[i];

        // get bounding box of triangle
        float tMinX = maxX;
        float tMaxX = minX;
        float tMinY = maxY;
        float tMaxY = minY;
        float tMinZ = maxZ;
        float tMaxZ = minZ;

        for (unsigned int j = 0; j < 3; j ++)
        {
            if (t.vertices[j][0] < tMinX) tMinX = t.vertices[j][0];
            if (t.vertices[j][0] > tMaxX) tMaxX = t.vertices[j][0];
            if (t.vertices[j][1] < tMinY) tMinY = t.vertices[j][1];
            if (t.vertices[j][1] > tMaxY) tMaxY = t.vertices[j][1];
            if (t.vertices[j][2] < tMinZ) tMinZ = t.vertices[j][2];
            if (t.vertices[j][2] > tMaxZ) tMaxZ = t.vertices[j][2];
        }

        // get voxel coordinates
        unsigned int vMinX = std::max((unsigned int) 0, std::min((unsigned int)
                                      ((tMinX - minX) / voxelNumX), voxelNumX - 1));
        unsigned int vMaxX = std::max((unsigned int) 0, std::min((unsigned int)
                                      ((tMaxX - minX) / voxelNumX), voxelNumX - 1));
        unsigned int vMinY = std::max((unsigned int) 0, std::min((unsigned int)
                                      ((tMinY - minY) / voxelNumY), voxelNumY - 1));
        unsigned int vMaxY = std::max((unsigned int) 0, std::min((unsigned int)
                                      ((tMaxY - minY) / voxelNumY), voxelNumY - 1));
        unsigned int vMinZ = std::max((unsigned int) 0, std::min((unsigned int)
                                      ((tMinZ - minZ) / voxelNumZ), voxelNumZ - 1));
        unsigned int vMaxZ = std::max((unsigned int) 0, std::min((unsigned int)
                                      ((tMaxZ - minZ) / voxelNumZ), voxelNumZ - 1));

        // DEBUG: Output min/max voxel coordinates
        //cout << "Voxel coordinates: X[" << vMinX << "-" << vMaxX << "], Y["
        //                                << vMinY << "-" << vMaxY << "], Z["
        //                                << vMinZ << "-" << vMaxZ << "]" << endl;

        // insert vertices into voxels
        // take all voxels which intersect with the bounding box of the triangle
        // TODO: optimize
        for (unsigned int z = vMinZ; z <= vMaxZ; z ++)
        {
            for (unsigned int y = vMinY; y <= vMaxY; y ++)
            {
                for (unsigned int x = vMinX; x <= vMaxX; x ++)
                {
                    voxels.voxels[z * voxelNumX * voxelNumY +
                                  y * voxelNumX + x].vertices.push_back(i);
//                    Voxel v = voxels.voxels[z * voxelNumX * voxelNumY +
//                                            y * voxelNumX + x];
//                    Triangle t = triangles[i];
//                    float vminx = v.pos[0], vminy = v.pos[1], vminz = v.pos[2];
//                    float vmaxx = vminx + voxels.voxelSize[0], vmaxy = vminy +
//                    voxels.voxelSize[1], vmaxz = vminz + voxels.voxelSize[2];
                }
            }
        }
        // TEST: assign all triangles to all voxels
        //for (unsigned int j = 0; j < voxels.voxels.size(); j ++)
        //{
        //    voxels.voxels[j].vertices.push_back(i);
        //}
    }

    // DEBUG: Output voxels with vertices
    //for (unsigned int i = 0; i < voxels.voxels.size(); i ++)
    //{
    //    cout << "Voxel " << i << " (" << voxels.voxels[i].vertices.size()
    //         << " Tri) (" << voxels.voxels[i].pos[0] << "/"
    //                      << voxels.voxels[i].pos[1] << "/"
    //                      << voxels.voxels[i].pos[2] << "):";
    //    for (unsigned int j = 0; j < voxels.voxels[i].vertices.size(); j ++)
    //    {
    //        cout << " " << voxels.voxels[i].vertices[j];
    //    }
    //    cout << endl;
    //}

    cout << "Time (" << triangles.size() << " triangles -> "
         << voxels.voxels.size() << " voxels): " << t.elapsed()  << " msec" << endl;
}


void Raytracer::genImage()
{
    QTime t;
    t.start();
    
    init();

    generateVoxels();

    int count = 0;

    unsigned int backgroundCode = backgroundColor.red() + backgroundColor.green()*255 + backgroundColor.blue()*255*255 + (unsigned int)((unsigned int)255*(unsigned int)255*(unsigned int)255*(unsigned int)255); //the alpha value of background is always set  to 1.0 (255)
    cout<<"Back_Code: "<<backgroundCode<<endl;
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fboId);    
    paintGL();
    GLdouble mvMatrix[16];
    GLdouble projMatrix[16];
    int viewPort[4];
    float* zValues = new float[image->height()*image->width()];
    glReadPixels( 0, 0, image->width(), image->height(), GL_DEPTH_COMPONENT, GL_FLOAT, zValues );
    glGetIntegerv(GL_VIEWPORT, viewPort);
    glGetDoublev(GL_MODELVIEW_MATRIX,mvMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);

    QWidget *w = new QWidget(NULL);
    QLabel *label = new QLabel(w);
    label->setGeometry(0,0,image->width()/superSamplingRate, image->height()/superSamplingRate);
    w->resize(label->size());
    w->setWindowTitle(QString("Rendered image"));
    w->show();
    

#pragma omp parallel for schedule(dynamic, 1)
  for (int j=0; j<image->height(); j+=1)
    {    
        for (int i=0; i<image->width(); i+=1)
        {
            unsigned int index = idx[j*image->width()+i];
            QColor c;
            if (index != backgroundCode)
            {
                double dX=0.0, dY=0.0, dZ=0.0;
                gluUnProject( (float)i, (float)j, zValues[j*image->width()+i], mvMatrix, projMatrix, viewPort, &dX, &dY, &dZ );
                Vector intersection(dX, dY, dZ);
                Vector dir = intersection - camera;
                //dir.normalize();
                c = raytrace(camera, dir, MAX_DEPTH);
                                                
            }
            else
            {
                c = backgroundColor;
            }

            image->setPixel(i,image->height()-(j+1), QRgb(c.rgb()));

        }

        #pragma omp critical
        {
            ++count;
            cout<<"\r----"<<(float)count/(float)image->height()*100.0<<"----";
        }
       
        // get the thread number (0 == master thread) 
        const int threadID = omp_get_thread_num();
        #pragma omp critical
        {
            if ((!(superSamplingRate <= 1.0) || omp_get_num_threads()==1) && (threadID == 0))
            {
                label->setPixmap(QPixmap::fromImage(image->scaled ( image->width()/superSamplingRate, image->height()/superSamplingRate, Qt::IgnoreAspectRatio, Qt::FastTransformation )));
                label->repaint();
            }
        }
    }
    delete zValues;
    cout<<endl;

    label->setPixmap(QPixmap::fromImage(image->scaled ( image->width()/superSamplingRate, image->height()/superSamplingRate, Qt::IgnoreAspectRatio, Qt::SmoothTransformation )));
    label->repaint();
    cout<<"Time: "<<t.elapsed()/1000.0<<endl;

    updateGL();
}


QColor Raytracer::raytrace(Vector start, Vector dir, int depth)
{
    QColor color = backgroundColor;

    if (depth <= 0)
    {
        return color;
    }

    // DEBUG:
    //cout << "start (" << start[0] << "/" << start[1] << "/" << start[2]
    //     << "), dir (" << dir[0] << "/" << dir[1] << "/" << dir[2] << ")" << endl;

    float minT = INFINITY, maxT = -INFINITY;
    float tempT;

    // ------------------------------------------------------------------------
    // get nearest and farthest t for ray

    // for all dimensions -> 3 d with 2 planes each
    for (unsigned int i = 0; i < 3; i ++)
    {
        // != division by 0
        if (dir[i] != 0)
        {
            tempT = (voxels.origin[i] - start[i]) / dir[i];
            // tempT > 0 for object behind us
            if (tempT > 0)
            {
                Vector p = start + dir * tempT;
                // check if point of square and not outside
                if ((p[(i + 1) % 3] >= voxels.origin[(i + 1) % 3]) &&
                    (p[(i + 1) % 3] <= voxels.origin[(i + 1) % 3] + voxels.size[(i + 1) % 3]) &&
                    (p[(i + 2) % 3] >= voxels.origin[(i + 2) % 3]) &&
                    (p[(i + 2) % 3] <= voxels.origin[(i + 2) % 3] + voxels.size[(i + 2) % 3]))
                {
                    if (tempT < minT) minT = tempT;
                    if (tempT > maxT) maxT = tempT;
                }
            }

            tempT = (voxels.origin[i] + voxels.size[i] - start[i]) / dir[i];
            if (tempT > 0)
            {
                Vector p = start + dir * tempT;
                if ((p[(i + 1) % 3] >= voxels.origin[(i + 1) % 3]) &&
                    (p[(i + 1) % 3] <= voxels.origin[(i + 1) % 3] + voxels.size[(i + 1) % 3]) &&
                    (p[(i + 2) % 3] >= voxels.origin[(i + 2) % 3]) &&
                    (p[(i + 2) % 3] <= voxels.origin[(i + 2) % 3] + voxels.size[(i + 2) % 3]))
                {
                    if (tempT < minT) minT = tempT;
                    if (tempT > maxT) maxT = tempT;
                }
            }
        }
    }

    // no intersection with box
    if ((minT == INFINITY) || (maxT == -INFINITY))
    {
        //cout << "outside" << endl;
        color.setRgb(0, 255, 0);
        return color;
    }

    // ------------------------------------------------------------------------
 
    int voxel[3], lastVoxel[3], step[3];
    Vector deltaT, nextCrossingT;

//    Vector rayI = start + dir * minT;
//    cout << "Inter min: (" << rayI[0] << "/" << rayI[1] << "/" << rayI[2] << ")" << endl;
//    rayI = start + dir * maxT;
//    cout << "Inter max: (" << rayI[0] << "/" << rayI[1] << "/" << rayI[2] << ")" << endl;


    // need voxel coordinates and step
    for (unsigned int i = 0; i < 3; i ++)
    {
        float rayOrigVoxel = ((start[0] + (dir[0] * minT)) - voxels.origin[i]);
        // DEBUG:
        //cout << "rayOrigVoxel = (" << rayOrigVoxel << "), voxel[" << i << "] = ("
        //     << rayOrigVoxel / voxels.voxelSize[i] << ")" << endl;
        voxel[i] = std::max(0, std::min((int) (rayOrigVoxel / voxels.voxelSize[i]),
                                        (int) voxels.resolution[i] - 1));
        if (dir[i] < 0)
        {
            deltaT[i] = - voxels.voxelSize[i] * 1 / dir[i]; // ?
            nextCrossingT[i] = minT + (voxel[i] * voxels.voxelSize[i] - rayOrigVoxel) * 1 / dir[i]; // ?

            lastVoxel[i] = -1;
            step[i] = -1;
        }
        else
        {
            deltaT[i] = voxels.voxelSize[i] * 1 / dir[i];
            nextCrossingT[i] = minT + ((voxel[i] + 1) * voxels.voxelSize[i] -
                                rayOrigVoxel) * 1 / dir[i];

            lastVoxel[i] = voxels.resolution[i];
            step[i] = 1;
        }
    }
    // DEBUG:
//    cout << "vStart = (" << voxel[0] << "/" << voxel[1] << "/" << voxel[2]
//         << "), vStop = (" << lastVoxel[0] << "/" << lastVoxel[1] << "/" << lastVoxel[2] 
//         << "), dT = (" << deltaT[0] << "/" << deltaT[1] << "/" << deltaT[2]
//         << "), nCT = (" << nextCrossingT[0] << "/" << nextCrossingT[1] << "/"
//                         << nextCrossingT[2] << ")" << endl;

    // ------------------------------------------------------------------------

    // raytrace voxels

    float lastT = INFINITY;
    float bn;
    int intersect = -1;
    const static unsigned int jumpMapAxis[8] = {2, 1, 2, 1, 2, 2, 0, 0};
    Vector p;

    while(true)
    {
        unsigned int v = voxel[2] * voxels.resolution[0] * voxels.resolution[1] +
                         voxel[1] * voxels.resolution[0] + voxel[0];
        Voxel vox = voxels.voxels[v];
//        cout << "-> [" << v << "] ";

        // compute triangle intersections
        for (unsigned int i = 0; i < vox.vertices.size(); i ++)
        {
            Triangle tri = triangles[vox.vertices[i]];
            bn = scalarProduct(dir, tri.planeNormal);

            // check if parallel
            if (bn == 0)
            {
                continue;
            }

            // get t =  ((p - e) * n) / (b * n)
            float t = scalarProduct((tri.vertices[0] - start), tri.planeNormal) / bn;
            if (t < 0)
            {
                continue;
            }

            // ray(t) = e + (b * t)
            p = start + (dir * t);
            //cout << "Intersection at (" << p[0] << "/" << p[1] << "/" << p[2] << ")" << endl;

            //float alpha0, alpha1, alpha2, area;
            Vector area0v, area1v, area2v;

            // area(p0, p1, p2) = 0.5 * || (p1 - p0) x (p2 - p0) ||
            //area = crossProduct(tri.vertices[1] - tri.vertices[0],
            //                    tri.vertices[2] - tri.vertices[0]).norm() * 0.5;

            // check if point in triangle
            area0v = crossProduct(tri.vertices[1] - p,
                                  tri.vertices[2] - p);
            if (scalarProduct(area0v, tri.planeNormal) < 0) continue;
            area1v = crossProduct(tri.vertices[2] - p,
                                  tri.vertices[0] - p);
            if (scalarProduct(area1v, tri.planeNormal) < 0) continue;
            area2v = crossProduct(tri.vertices[0] - p,
                                  tri.vertices[1] - p);
            if (scalarProduct(area2v, tri.planeNormal) < 0) continue;

            //alpha0 = (area0v.norm() * 0.5) / area;
            //alpha1 = (area1v.norm() * 0.5) / area;
            //alpha2 = (area2v.norm() * 0.5) / area;
            //float alpha = alpha0 + alpha1 + alpha2;
            //if (alpha < 0.99999 || alpha > 1.00001) continue;

            //Vector p = alpha0 * tri.vertices[0] +
            //           alpha1 * tri.vertices[1] +
            //           alpha2 * tri.vertices[2];

            if (lastT > t)
            {
                lastT = t;
                intersect = i;
            }
        }

//        cout << "dir = (" << dir[0] << "/" << dir[1] << "/" << dir[2] << "), t = (" << lastT
//             << "), p = (" << p[0] << "/" << p[1] << "/" << p[2] << "), tri = (" << intersect
//             << "), vox = (" << v  << ")" << endl;


        if (intersect != -1)
        {
            // compute color
            float col = 255 * ((intersect + 1) / (float) triangles.size());
            color.setRgb(col, col, col);

            for (unsigned int i = 0; i < lights.size(); i ++)
            {
                // compute ray to light source
                // ...
            }

            // no other voxel to test
            break;
        }

        unsigned int k = ((nextCrossingT[0] < nextCrossingT[1]) << 2) +
                         ((nextCrossingT[1] < nextCrossingT[2]) << 1) +
                         ((nextCrossingT[2] < nextCrossingT[0]));
        unsigned int axis = jumpMapAxis[k];

        //cout << "maxT = " << maxT << ", nextCrossingT[" << axis << "] = " << nextCrossingT[axis] << endl;
        if (maxT < nextCrossingT[axis]) break;
        voxel[axis] += step[axis];
        //cout << "voxel[" << axis << "] = " << voxel[axis]
        //     << ", lastVoxel[" << axis << "] = " << lastVoxel[axis] << endl;
        if (voxel[axis] == lastVoxel[axis]) break;
        nextCrossingT[axis] += deltaT[axis];
    }
//    cout << endl;

//    cout << "dir = (" << dir[0] << "/" << dir[1] << "/" << dir[2] << "), lastT = (" << lastT
//         << "), p = (" << p[0] << "/" << p[1] << "/" << p[2] << "), v = (" << intersect
//         << "), voxel = (" << v  << ")" << endl;


    // no voxel with no triangle for intersection ...
    if (intersect == -1)
    {
        color.setRgb(255, 0, 0);
    }

    // ------------------------------------------------------------------------

    return color;
}


void Raytracer::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_S)
    {
        QString path = QFileDialog::getSaveFileName ( NULL, QString("Save Image"));
        finalImage.save ( path, "PNG", -1 );
    }
}

