#include "Raytracer.h"
#include <fstream>
#include <omp.h>

#include <iostream>
using namespace std;

Raytracer::Raytracer( QString path )
{
// display dis/enabled constants used for rendering
#ifndef TEXTURE_ON
    cout << "- Textures disabled!" << endl;
#else
    cout << "- Using textures." << endl;
#endif

#ifndef REFLECTION_ON
    cout << "- Reflection rays disabled! (faster rendering)" << endl;
#else
    cout << "- Reflection enabled." << endl;
#endif

#ifndef REFRACTION_ON
    cout << "- Refraction rays disabled! (faster rendering)" << endl;
#else
    cout << "- Refraction enabled." << endl;
#endif

#ifndef SHADOWS_ON
    cout << "- Shadows disabled!" << endl;
#else
    cout << "- Shadows enabled." << endl;
#endif

#ifdef INVERT_SHADOWS
    cout << "- Showing only shadow part! (Debug?)" << endl;
#endif

    QString basepath = path;
    basepath.chop(path.length() - path.lastIndexOf("/"));
    cout << "Base path: " << basepath.toStdString() << endl;

    superSamplingRate = 1.1; //TODO: setting this to 1 and enabling on the fly output will cause a crash
    //read file and init everything

    vector<Vector> vectors;
    vector<Vector> normals;
    vector<Material> materials;
    vector<Vector> tex_coords;

    int nLights = 0, nMat =0, nVert=0, nNorm=0, nTexCoords=0, nFaces=0;

    // take loading time
    QTime t;
    t.start();

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
            QString filepath; //path in were the
            filepath = basepath;
            filepath.append(QString("/%1").arg(QString(texName.c_str())));
            QImage texture(filepath);
            if(!texture.isNull())
            {
                cout << "Texture: " << texName << endl;
                mat.texture = texture;
                mat.isTexture = true;
            }
            else
            {
                //QMessageBox::warning(this, tr("Texture Loader"),
                //                     tr("Texture couldn't be loaded!"), QMessageBox::Ok);
                cout << "Texture (" << texName << ") couldn't be loaded!" << endl;
            }
        }
        if (isTexture == 2) //load normal map (heavy code duplicates TODO: change that)
        {
            string bumpName;
            file>>bumpName;
            QString filepath;
            filepath = basepath;
            filepath.append(QString("/%1").arg(QString(bumpName.c_str())));
            QImage bumpMap(filepath);
            if(!bumpMap.isNull())
            {
                cout << "Normal Map: " << bumpName << endl;
                mat.hasNormalMap = true;
                mat.normalMap = bumpMap;
            }
            else
            {
                //QMessageBox::warning(this, tr("Texture Loader"),
                //                     tr("Normal map couldn't be loaded!"), QMessageBox::Ok);
                cout << "Normal Map (" << bumpName << ") couldn't be loaded!" << endl;
            }
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
    float minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX, maxx = FLT_MIN, maxy = FLT_MIN, maxz = FLT_MIN;
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

        // computing planeNormals
        t.planeNormal = crossProduct(t.vertices[0] - t.vertices[1],
                                     t.vertices[0] - t.vertices[2]);
        t.planeNormal.normalize();

        // acceleration with precomputation
        Vector b = t.vertices[1]-t.vertices[0];
        Vector c = t.vertices[2]-t.vertices[0];
        // scalarProduct(v,v) == v.normSquare();
        float d = (scalarProduct(c,c) * scalarProduct(b,b)) - (scalarProduct(c,b) * scalarProduct(b,c));

        t.ubeta = (b * (scalarProduct(c, c) / d)) - (c * (scalarProduct(c, b) / d));
        t.ugamma = (c * (scalarProduct(b, b) / d)) - (b * (scalarProduct(b, c) / d));
        t.kbeta = scalarProduct(-t.vertices[0], t.ubeta);
        t.kgamma = scalarProduct(-t.vertices[0], t.ugamma);

        // get bounding box for octree
        // vector 1
        if ( minx > t.vertices[0][0] )
            minx = t.vertices[0][0];
        else if ( maxx < t.vertices[0][0] )
            maxx = t.vertices[0][0];
        if ( miny > t.vertices[0][1] )
            miny = t.vertices[0][1];
        else if ( maxy < t.vertices[0][1] )
            maxy = t.vertices[0][1];
        if ( minz > t.vertices[0][2] )
            minz = t.vertices[0][2];
        else if ( maxz < t.vertices[0][2] )
            maxz = t.vertices[0][2];

        // vector 2
        if ( minx > t.vertices[1][0] )
            minx = t.vertices[1][0];
        else if ( maxx < t.vertices[1][0] )
            maxx = t.vertices[1][0];
        if ( miny > t.vertices[1][1] )
            miny = t.vertices[1][1];
        else if ( maxy < t.vertices[1][1] )
            maxy = t.vertices[1][1];
        if ( minz > t.vertices[1][2] )
            minz = t.vertices[1][2];
        else if ( maxz < t.vertices[1][2] )
            maxz = t.vertices[1][2];

        // vector 3
        if ( minx > t.vertices[2][0] )
            minx = t.vertices[2][0];
        else if ( maxx < t.vertices[2][0] )
            maxx = t.vertices[2][0];
        if ( miny > t.vertices[2][1] )
            miny = t.vertices[2][1];
        else if ( maxy < t.vertices[2][1] )
            maxy = t.vertices[2][1];
        if ( minz > t.vertices[2][2] )
            minz = t.vertices[2][2];
        else if ( maxz < t.vertices[2][2] )
            maxz = t.vertices[2][2];

        if (materials[matNr].isTexture)
        {
            t.texCoords[0]=tex_coords[TexCoordsNr1];
            t.texCoords[1]=tex_coords[TexCoordsNr2];
            t.texCoords[2]=tex_coords[TexCoordsNr3];
        }
        triangles.push_back(t);
    }
    cout << "Loading data: " << t.restart()/1000.0f << " sec" << endl;
    cout << "Got " << triangles.size() << " triangles" << endl;

    // building Octree
    minx -= BOX_MARGIN;
    maxx += BOX_MARGIN;
    miny -= BOX_MARGIN;
    maxy += BOX_MARGIN;
    minz -= BOX_MARGIN;
    maxz += BOX_MARGIN;

    octree = new Octree();
    octree->build(&triangles, minx, miny, minz, maxx, maxy, maxz);

    cout << "Building linear Octree (" << triangles.size() << " triangles ~~> "
         << octree->size() << " voxels): " << t.elapsed()/1000.0f << " sec" << endl;
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

void Raytracer::close()
{
    setVisible(false);
}

void Raytracer::genImage()
{
    QTime t;
    t.start();
    
    init();

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
            cout<<"\r----"<<(float)count/(float)image->height()*100.0f<<"----";
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
    cout << "Time rendering: " << t.elapsed()/1000.0f << endl;

    updateGL();
}



QColor Raytracer::raytrace(Vector start, Vector dir, int depth)
{
    return raytrace(start, dir, depth, 1.0f);
} 

QColor Raytracer::raytrace(Vector start, Vector dir, int depth, float density)
{
    // abort if too much recursions
    if (depth <= 0)
    {
        return backgroundColor;
    }
    depth --;

    float dis = FLT_MAX;
    Triangle triangle;
    Vector p;

    for (int i = 0; i < octree->size(); i ++)
    {
        // cut each voxel
        if (octree->cutVoxel(i, &start, &dir, dis))
        {
            Triangle tr;
            Vector q;

            // cut each triangle in voxel
            float tmp = octree->cutTriangles(i, &start, &dir, NULL, &tr, &q);
            if (tmp < dis)
            {
                dis = tmp;
                triangle = tr;
                p = q;
            }
        }
    }

    if (dis != FLT_MAX)
    {
        float alpha[3];

        // alpha values for interpolation
        alpha[1] = dot(p, triangle.ubeta) + triangle.kbeta;
        alpha[2] = dot(p, triangle.ugamma) + triangle.kgamma;
        alpha[0] = 1.0f - alpha[1] - alpha[2];

        // normal on intersection point
        // n = alpha0 * n0 + alpha1 * n1 + alpha2 * n2
        Vector pn = triangle.normals[0] * alpha[0] +
                    triangle.normals[1] * alpha[1] +
                    triangle.normals[2] * alpha[2];
        pn.normalize(); // might be faster ??
                        // you wouldn't think so but with normalization
                        // it is slightly faster ...

        // inverse view direction
        // a = - b (== dir)
        // or: invDir = (p - start) ?
        Vector invDir = dir * (-1.f);
        // should also always be normalized
        invDir.normalize();

        // really important !
        // check if angle greater 90° -> invert normal vector at intersection
        float inv = scalarProduct(invDir, pn);
        if (inv < 0) pn = pn * (-1.f);


        // 1st compute color

        Material mat = triangle.material;
        float r = 0.0f, g = 0.0f, b = 0.0f;
        //float r = 1.0f, g = 1.0f, b = 1.0f;
        // use 1.0f (and comment out next) if you want a scene with too much light ...
        r = ambientLight.redF()   * mat.ambient[0];
        g = ambientLight.greenF() * mat.ambient[1];
        b = ambientLight.blueF()  * mat.ambient[2];


        // --------------------------------------------------------------------
        // for all lightsources
        for (unsigned int i = 0; i < lights.size(); i ++)
        {
            Lightsource lig = lights[i];
            // l - vector to lightsource
            Vector l = lig.position - p;
            // d - distance to lightsource
            float d = l.norm();
            // great effect if not normalized ...
            l.normalize();
            
            // reflection ray
            // r = 2 * (n x l) * n - l
            Vector refl = (pn * scalarProduct(pn, l) * 2) - l;
            //refl.normalize();

            // check visibility ?
            float pnl = scalarProduct(pn, l);
            // n * l < 90° --> pnl > 0 (cos)
            if (pnl < 0) continue;

            // 0.0f .. 1.0f refraction factor
            float alphaFactor = 1.f;

#ifdef SHADOWS_ON
            // check shadows
            // check for intersection with triangle before light source
            bool shadow = false;
            for (int j = 0; j < octree->size(); j ++)
            {
                if (octree->cutVoxel(j, &p, &l, d))
                {
                    shadow = octree->cutTriangles(j, &p, &l, &triangle, d, &alphaFactor);
                    if (shadow) break;
                }
            }
#endif

#ifdef INVERT_SHADOWS
            // invert (! shadow) if you want only the shadow part ...
            shadow = ! shadow;
#endif

#ifdef SHADOWS_ON
            if (! shadow)
            {
#endif
                float fatt = 1.f;
                if ((lig.constAtt != 0) && (lig.linAtt != 0) && (lig.quadAtt != 0))
                {
                    fatt = 1.f / (lig.constAtt + d * lig.linAtt + d * d * lig.quadAtt);
                }

                // f_att(d) * (I_d * O_d * (n x l) + I_s * O_s * (max(0, (r x a)))^s
                float ridp = pow(std::max(0.f, scalarProduct(refl, invDir)), mat.shininess);
                r += fatt * (lig.diffuse[0] * mat.diffuse[0] * pnl +
                             lig.specular[0] * mat.specular[0] * ridp);// * alphaFactor;
                g += fatt * (lig.diffuse[1] * mat.diffuse[1] * pnl +
                             lig.specular[1] * mat.specular[1] * ridp);// * alphaFactor;
                b += fatt * (lig.diffuse[2] * mat.diffuse[2] * pnl +
                             lig.specular[2] * mat.specular[2] * ridp);// * alphaFactor;
#ifdef SHADOWS_ON
            }
#endif
        }

        // 2nd compute (recursive) effects (reflection/refraction)

#ifdef REFLECTION_ON
        // --------------------------------------------------------------------
        // reflection ray (recursive)
        if (mat.sharpness != 0.f)
        {
            Vector refl = (pn * scalarProduct(pn, invDir) * 2) - invDir;
            QColor c_refl = raytrace(p, refl, depth, density);

            r = (1.f - mat.sharpness) * r + mat.sharpness * c_refl.redF();
            g = (1.f - mat.sharpness) * g + mat.sharpness * c_refl.greenF();
            b = (1.f - mat.sharpness) * b + mat.sharpness * c_refl.blueF();
        }
#endif

#ifdef REFRACTION_ON
        // --------------------------------------------------------------------
        // refraction ray (e. g. glass) (recursive)
        // TODO: check code - not sure if properly working ... :(
        if (mat.alpha < 1.0f)
        {
            Vector refrac;
            float c1 = scalarProduct(invDir, triangle.planeNormal);
            float ns = density / mat.density;
            float c2_nw = ns * ns * (1.f - c1 * c1);

            if (c2_nw > 1.f)
            {
                refrac = dir - (triangle.planeNormal * c1 * 2.f);
            }
            else
            {
                //refrac = (dir * ns) + (triangle.planeNormal * (ns * (-c1) - sqrt(1.f - c2_nw)));
                refrac = (dir * ns) - (triangle.planeNormal * (ns - sqrt(1.f - c2_nw)));
                density = mat.density;
            }
            //refrac = (dir * ns) - (triangle.planeNormal * (ns - sqrt(1.f - c2_nw)));
            //refrac = (dir * ns) + (triangle.planeNormal * (ns * (-c1) - sqrt(1.f - c2_nw)));

            QColor c_refrac = raytrace(p, refrac, depth, density);

            r = mat.alpha * c_refrac.redF()   + (1.f - mat.alpha) * r;
            g = mat.alpha * c_refrac.greenF() + (1.f - mat.alpha) * g;
            b = mat.alpha * c_refrac.blueF()  + (1.f - mat.alpha) * b;
        }
#endif

        // 3rd texture

#ifdef TEXTURE_ON
        // --------------------------------------------------------------------
        // if material of triangle is a texture get mapping
        if (triangle.material.isTexture)
        {
            // get texture coordinate with interpolation
            Vector tex = triangle.texCoords[0] * alpha[0] +
                         triangle.texCoords[1] * alpha[1] +
                         triangle.texCoords[2] * alpha[2];
            int x = std::abs((int)(tex[0] * (triangle.material.texture.width() - 1)))
                        % triangle.material.texture.width();
            int y = std::abs((int)((1 - tex[1]) * (triangle.material.texture.height() - 1)))
                        % triangle.material.texture.height();

            // get color of point on texture
            QColor tc = triangle.material.texture.pixel(x, y);

            r *= tc.redF();
            g *= tc.greenF();
            b *= tc.blueF();

        }
#endif

        return QColor(max(min((int)(r * 255.f), 255), 0),
                      max(min((int)(g * 255.f), 255), 0),
                      max(min((int)(b * 255.f), 255), 0)); 
    }

    // default background
    return backgroundColor;
}

void Raytracer::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_S)
    {
        QString path = QFileDialog::getSaveFileName ( NULL, QString("Save Image"));
        finalImage.save ( path, "PNG", -1 );
    }
}

