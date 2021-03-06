#ifndef CGMATH_H
#define CGMATH_H

//DataTypes, useful Functions and Math Stuff

#include <vector>
#include <math.h>
#include <string>
#include <cstring>
#include <QImage>

#define EPSILON 0.0001f

#define dot(v1, v2) scalarProduct(v1, v2)
#define cross(v1, v2) crossProduct(v1, v2)

class Vector
{
    public:
        Vector()
        {
            values[0] = values[1] = values[2] = 0.0;
              values[3] = 1.0;
        }

        Vector(float x, float y, float z=0, float w=1.0)
        {
            values[0] = x;
            values[1] = y;
            values[2] = z;
            values[3] = w;
        }

        
        float& operator () (size_t i)
        {
            return values[i];
        }

        float operator () (size_t i) const
        {
            return values[i];
        }

        float& operator [] (size_t i)
        {
            return values[i];
        }

        float operator [] (size_t i) const
        {
            return values[i];
        }

        void setValue(int index, float value)
        {
            values[index] = value;
        }

        Vector& operator = ( const Vector& v ) 
        {
            memcpy(values, v.values, 4*sizeof(float));
              return *this;
        }

        bool operator == ( const Vector& v )
        {
            return (values[0] == v[0]) && (values[1] == v[1]) && (values[2] == v[2]);
        }

        Vector operator + ( const Vector& v ) const
        {
            return Vector( values[0] + v[0], values[1] + v[1], values[2] + v[2] );
        }

        Vector operator - ( const Vector& v ) const
        {
            return Vector( values[0] - v[0], values[1] - v[1], values[2] - v[2] );
        }

        Vector operator - ( ) const
        {
            return Vector( -values[0], -values[1], -values[2] );
        }

        Vector operator * ( const float& scalar ) const
        {
            return Vector( values[0]*scalar, values[1]*scalar, values[2]*scalar );
        }

        Vector operator / ( const float& scalar ) const
        {
            return Vector( values[0]/scalar, values[1]/scalar, values[2]/scalar );
        }

        double normSquare ()
        {
            return values[0]*values[0] + values[1]*values[1] + values[2]*values[2];
        }

        double norm ()
        {
            return sqrt( values[0]*values[0] + values[1]*values[1] + values[2]*values[2] );
        }

        void normalize()
        {
            double n = norm();
            values[0] /= n;
            values[1] /= n;
            values[2] /= n;
        }

        void invert()
        {
             values[0] *= -1.0f;
             values[1] *= -1.0f;
             values[2] *= -1.0f;
        }

        float* getValues ()
        {
            return values;
        }

    private:
        float values[4];
};

inline float scalarProduct( Vector v1, Vector v2 )
{
      return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
};

inline Vector crossProduct( Vector v1, Vector v2 )
{
    return Vector(v1[1]*v2[2] - v2[1]*v1[2], v1[2]*v2[0] - v2[2]*v1[0], v1[0]*v2[1] - v2[0]*v1[1]);
};


//since this programm works only with triangles every vector in the Face-struct has dimension 3
struct Face
{
    std::vector<unsigned int> vectors;
    std::vector<unsigned int> normals;
    std::vector<unsigned int> tex_coords;
};

struct Mesh
{
    std::vector<unsigned int> faces;
    unsigned int material;
};


//--------------esp. for the raytracer-----------------------


struct Lightsource
{
    Vector position;
    float ambient[3];
    float diffuse[3];
    float specular[3];
    float constAtt;
    float linAtt;
    float quadAtt;    
};

struct Material
{
    std::string name;
    float ambient[4];
    float diffuse[4];
    float specular[4];
    float shininess;
    float alpha; //how light is going through an object (1.0 no transparency, 0.0 total transparency)
    float sharpness; //how is does the material reflect (1.0 -> perfect mirror, 0.0 -> no reflections)
    float density; //index of refraction i.e. how does the light bend when going through an object (1.0 means no bending)
    //diffuse (ambient) map
    bool isTexture;
    unsigned int tex_id;
    std::string texName;
    QImage texture;
    
    //normal map -- not important to the raytracer
    bool hasNormalMap;
    unsigned int normalMap_id;
    std::string normalMapName;
    QImage normalMap;
};


struct Triangle //just to have all the informations of Material and 
{
    Material material;
    Vector vertices[3];
    Vector normals[3];
    Vector texCoords[3]; //only [0] and [1] are used to store texInformations
    Vector planeNormal;

    Vector ubeta;
    Vector ugamma;
    float kbeta;
    float kgamma;

    bool operator == ( const Triangle& t )
    {
        return (vertices[0] == t.vertices[0]) && (vertices[1] == t.vertices[1]) &&
               (vertices[2] == t.vertices[2]) && (normals[0] == t.normals[0]) &&
               (normals[1] == t.normals[1]) && (normals[2] == t.normals[2]) &&
               (texCoords[0] == t.texCoords[0]) && (texCoords[1] == t.texCoords[1]) &&
               (texCoords[2] == t.texCoords[2]) && (planeNormal == t.planeNormal) &&
               (ubeta == t.ubeta) && (ugamma == t.ugamma) && (kbeta == t.kbeta) && (kgamma == t.kgamma);
    }
};

/**
 * Checks whether the triangle is intersected by
 * the ray. Returns true if it is so and the
 * parameter t of the ray function.
 */
inline bool cut(Vector *start, Vector *dir, Triangle *triangle, float *t)
{
    float d;
    if ( (d = dot(*dir, (*triangle).planeNormal)) != 0  )
    {
        float t_temp = dot(((*triangle).vertices[0] - *start), (*triangle).planeNormal) / d;

        if ( t_temp <= EPSILON )
        {
            return false;
        }

        Vector p_temp = *start + *dir * t_temp;
        float beta = dot(p_temp, (*triangle).ubeta) + (*triangle).kbeta;

        if ( beta < 0 )
        {
            return false;
        }

        float gamma = dot(p_temp, (*triangle).ugamma) + (*triangle).kgamma;

        if ( gamma < 0 )
        {
            return false;
        }

        if ( (beta + gamma) > 1 )
        {
            return false;
        }

        *t = t_temp;
        return true;
    }

    return false;
};

/**
 * Checks whether the triangle is intersected by
 * the ray. Returns true if it is so and the
 * intersection point.
 */
inline bool cut(Vector *start, Vector *dir, Triangle *triangle, Vector *p)
{
    float t;
    bool r;
    if ( (r = cut(start, dir, triangle, &t)) )
    {
            *p = *start + *dir * t;
    }

    return r;
};

#endif
 
