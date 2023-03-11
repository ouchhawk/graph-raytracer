#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

namespace parser {
    struct Vec3f {
        float x, y, z;
        Vec3f() : x(0), y(0), z(0) {}
        Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
        Vec3f operator + (Vec3f v) { return Vec3f(x+v.x, y+v.y, z+v.z); }
        Vec3f operator + (float f) { return Vec3f(x+f, y+f, z+f); }
        Vec3f operator - (Vec3f v) { return Vec3f(x-v.x, y-v.y, z-v.z); } // distance
        Vec3f operator - (float f) { return Vec3f(x-f, y-f, z-f); }
        Vec3f operator * (float a) { return Vec3f(x*a, y*a, z*a); }
        Vec3f operator / (float a) { return Vec3f(x/a, y/a, z/a); }
        bool operator < (Vec3f v) { return  (x<v.x && y<v.y && z<v.z); }
        bool operator > (Vec3f v) { return  (x>v.x && y>v.y && z>v.z); }
        float operator * (Vec3f v) { return x*v.x + y*v.y + z*v.z; } // dot
        Vec3f operator ^ (Vec3f v) { // cross
            return Vec3f(y*v.z - v.y*z, v.x*z - x*v.z, x*v.y - v.x*y); 
        }
        double length() { return sqrt(x*x + y*y + z*z); }
        Vec3f normalize() { return Vec3f(x/length(), y/length(), z/length()); }
        void print() const { cout << x << " " << y << " " << z << endl; }
    };

    struct Vec3i {
        int x, y, z;
        void print() const { cout << x << " " << y << " " << z << endl; }
    };

    struct Vec4f {
        float x, y, z, w;
        void print() const { cout << x << " " << y << " " << z << " " << w << endl; }
    };

    struct Matrix {
        float ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ox, oy, oz;
    };

    struct Camera {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        float near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight {
        Vec3f position;
        Vec3f intensity;
    };

    struct Material {
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
    };

    struct Face {
        int v0_id;
        int v1_id;
        int v2_id;

        Face() {}
        Face(int v0, int v1, int v2) : v0_id(v0), v1_id(v1), v2_id(v2){}
    };

    struct Mesh {
        int material_id;
        std::vector<Face> faces;
    };

    struct Triangle {
        int material_id;
        Face indices;
    };

    struct Sphere {
        int material_id;
        int center_vertex_id;
        float radius;
    };

    struct Check {
        bool hit = false;
        Vec3f hitPos;   // pos
        Vec3f n;        // normal
        int m=0;        // material_id
        string type;    // m - mesh, s - sphere, t - tri  for shadow
        int sira;       // for shadow
        float t;
        Check(bool hit) : hit(hit) {}
        Check(bool hit, Vec3f hitPos) : hit(hit), hitPos(hitPos) {}
        Check(bool hit, Vec3f hitPos, float t) : hit(hit), hitPos(hitPos), t(t) {}

    };

    struct Scene {
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;

        void loadFromXml(const std::string& filepath);
    };
}

#endif