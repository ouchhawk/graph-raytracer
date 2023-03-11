#include <iostream>
#include "parser.h"
#include "ppm.h"

using namespace std;
using namespace parser;

//GROUP40

Vec3f applyEffects(Scene &scene, Vec3f pos, Check &result, Vec3f v, Vec3f wi, int light_id, vector<vector<Vec3f>> &normals, Vec3f amb, Vec3f d, int recursion);
Vec3f findColor(Vec3f pos, Check &result, Vec3f amb, Scene &scene, Vec3f d, vector<vector<Vec3f>> &normals, int recursion);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - FUNCTIONS
float determinant(float a, float d, float g, float b, float e, float h, float c, float f, float i) {
    return a*(e*i - h*f) + b*(g*f - d*i) + c*(d*h - e*g);
}

float beta(Matrix a) {
    float acx = a.ax - a.cx;
    float acy = a.ay - a.cy;
    float acz = a.az - a.cz;
    float x = determinant(a.ax - a.ox, acx, a.dx, a.ay - a.oy, acy, a.dy, a.az - a.oz, acz, a.dz);
    float y = determinant(a.ax - a.bx, acx, a.dx, a.ay - a.by, acy, a.dy, a.az - a.bz, acz, a.dz);
    if(y)
        return x/y;
    
    return 0;
}

float gamma(Matrix a) {
    float abx = a.ax - a.bx;
    float aby = a.ay - a.by;
    float abz = a.az - a.bz;
    float x = determinant(abx, a.ax - a.ox, a.dx, aby, a.ay - a.oy, a.dy, abz, a.az - a.oz, a.dz);
    float y = determinant(abx, a.ax - a.cx, a.dx, aby, a.ay - a.cy, a.dy, abz, a.az - a.cz, a.dz);
    if(y)
        return x/y;
    
    return 0;
}

float t(Matrix a) {
    float abx = a.ax - a.bx;
    float acx = a.ax - a.cx;
    float aby = a.ay - a.by;
    float acy = a.ay - a.cy;
    float abz = a.az - a.bz;
    float acz = a.az - a.cz;
    float x = determinant(abx, acx, a.ax - a.ox, aby, acy, a.ay - a.oy, abz, acz, a.az - a.oz);
    float y = determinant(abx, acx, a.dx, aby, acy, a.dy, abz, acz, a.dz);
    if(y)
        return x/y;
    
    return 0;
}

Vec3f triangleNormal(vector<Vec3f> data, Face t) {
    Vec3f p0 = data[t.v0_id - 1];
    Vec3f p1 = data[t.v1_id - 1];
    Vec3f p2 = data[t.v2_id - 1];
    Vec3f crossV = (p1 - p0)^(p2 - p0);
    Vec3f normal = crossV.normalize();
    return normal;
}

Check rayTri(vector<Vec3f> &data, Vec3f e, Vec3f d, Face c, int mat) { // e = camPos
    Matrix A;
    A.ax = data[c.v0_id - 1].x;
    A.ay = data[c.v0_id - 1].y;
    A.az = data[c.v0_id - 1].z;
    A.bx = data[c.v1_id - 1].x;
    A.by = data[c.v1_id - 1].y;
    A.bz = data[c.v1_id - 1].z;
    A.cx = data[c.v2_id - 1].x;
    A.cy = data[c.v2_id - 1].y;
    A.cz = data[c.v2_id - 1].z;
    A.dx = d.x;
    A.dy = d.y;
    A.dz = d.z;
    A.ox = e.x;
    A.oy = e.y;
    A.oz = e.z;

    float tCramer = t(A);
    float betaCramer = beta(A);
    float gammaCramer = gamma(A);
    if(betaCramer + gammaCramer <= 1 && betaCramer >= 0 && gammaCramer >= 0 && tCramer >= 0) {
        Check res = Check(true, e + d*tCramer, tCramer);
        res.n = triangleNormal(data, c);
        res.m = mat;
        res.t = tCramer;
        return res;
    }
    return Check(false);
}
//wi - towards light
bool shadowTriangle(vector<Vec3f> data, Vec3f point, Vec3f wi, Face c, float epsilon, Vec3f wiNotNormalized) {
    Vec3f x = point + wi*epsilon;
    
    Matrix A;
    A.ax = data[c.v0_id - 1].x;
    A.ay = data[c.v0_id - 1].y;
    A.az = data[c.v0_id - 1].z;
    A.bx = data[c.v1_id - 1].x;
    A.by = data[c.v1_id - 1].y;
    A.bz = data[c.v1_id - 1].z;
    A.cx = data[c.v2_id - 1].x;
    A.cy = data[c.v2_id - 1].y;
    A.cz = data[c.v2_id - 1].z;
    A.dx = wi.x;
    A.dy = wi.y;
    A.dz = wi.z;
    A.ox = point.x;
    A.oy = point.y;
    A.oz = point.z;

    float tCramer = t(A);
    float betaCramer = beta(A);
    float gammaCramer = gamma(A);
    if(betaCramer + gammaCramer < 1 && betaCramer > 0 && gammaCramer > 0 && tCramer > 0) {
        if(tCramer/(wiNotNormalized.length()/wi.length()) < 1)
            return true;
    }
    return false;
}

Check raySphere(vector<Vec3f> &data, Vec3f e, Vec3f d, Sphere sphere, int mat) {
    Vec3f c = data[sphere.center_vertex_id - 1];
    Vec3f ec = e - c;
    float dd = d*d;

    float dis = pow((d*ec), 2) - dd*(ec*ec - pow(sphere.radius, 2));
    if(dis >= 0) {
        float t = -(d*ec + sqrt(dis))/dd;
        if(t > 0) {
            Check res = Check(true, e + d*t, t);
            res.n = ((e + d*t - c) / sphere.radius).normalize();
            res.m = mat;
            res.t = t;
            return res;
        }
    }
    return Check(false);
}

bool shadowSphere(Vec3f point, Vec3f wi, Vec3f c, float r, float epsilon, Vec3f wiNotNormalized) {
    Vec3f x = point + wi*epsilon;
    Vec3f xc = x - c;
    float dis = pow((wi*xc), 2) - wi*wi*(xc*xc - pow(r, 2));
    if(dis >= 0) {
        if(-(wi*(point - c) + sqrt(dis))/(wi*wi) > 0) { //isikla noktanin arasinda cisim var mi diye check et, objenin arkasinda kalan objeleri degil
            if((-(wi*xc + sqrt(dis))/(wi*wi))/(wiNotNormalized.length()/wi.length()) < 1)
                return true;
        }
    }
    return false;
}

Check rayMesh(vector<Vec3f> &data, Vec3f e, Vec3f d, Mesh mesh, int mat, vector<vector<Vec3f>> &normals, int j) {
    Check min = Check(false, Vec3f(999, 999, 999));
    for (int i = 0; i < mesh.faces.size(); ++i) {
        if( normals[j][i]*d <= 0 ) { // kameranin gormedigi arkada kalan yarisini hesaplamiyoruz
            Check c = rayTri(data, e, d, mesh.faces[i], mat);
            if(c.hit && (c.hitPos.length() < min.hitPos.length())) { //&& c.t > 0
                min = c;
                min.hit=true;
                // min.t = c.t; 
            }
        }
    }
    if(min.hit) 
        return min;

    return Check(false);
}

bool shadowMesh(vector<Vec3f> data, Vec3f point, Vec3f wi, vector<Face> faces, float epsilon, Vec3f wiNotNormalized, vector<vector<Vec3f>> &normals, int j3) {
    for (int i = 0; i < faces.size(); ++i) {
        if( normals[j3][i]*wi <= 0 ) {
            if (shadowTriangle(data, point, wi, faces[i], epsilon, wiNotNormalized)) 
                return true;
        }
    }
    return false;
}

Check findClosest(vector<Vec3f> &data, Vec3f e, Vec3f d, vector<Triangle> &triangles, vector<Sphere> &spheres, vector<Mesh> &meshes, vector<vector<Vec3f>> &normals) {
    Check tmin = Check(false, Vec3f(999, 999, 999));

    for (int x = 0; x < triangles.size(); ++x) {
        Check t = rayTri(data, e, d, triangles[x].indices, triangles[x].material_id);
        if(tmin.hitPos.length() > t.hitPos.length() && t.hitPos.length() > 0 && t.t>0) {
            tmin = t;
            tmin.type = "t";
            tmin.sira = x;
        }
    }
    for (int x = 0; x < spheres.size(); ++x) {
        Check s = raySphere(data, e, d, spheres[x], spheres[x].material_id);
        if(tmin.hitPos.length() > s.hitPos.length() && s.hitPos.length() > 0 && s.t>0) {
            tmin = s;
            tmin.type = "s";
            tmin.sira = x;
        }
    }
    for (int x = 0; x < meshes.size(); ++x) {
        Check m = rayMesh(data, e, d, meshes[x], meshes[x].material_id, normals, x);
        if(tmin.hitPos.length() > m.hitPos.length() && m.hitPos.length() > 0 && m.t>0) {
            tmin = m;
            tmin.type = "m";
            tmin.sira = x;
        }
    }
    return tmin;
}

Vec3f ambient(Vec3f ambient_light, vector<Material> materials, int m) {
    Vec3f I, ka;
    I = ambient_light;
    ka = materials[m].ambient; // !!!
    return Vec3f(I.x*ka.x, I.y*ka.y, I.z*ka.z);
}

Vec3f applyEffects(Scene &scene, Vec3f pos, Check &result, Vec3f v, Vec3f wi, int light_id, vector<vector<Vec3f>> &normals, Vec3f amb, Vec3f d, int recursion) { 
    int phongExp;
    Vec3f Lpos, I, h, kd, ks, km, w0, wr;
    I = scene.point_lights[light_id].intensity;
    Lpos = scene.point_lights[light_id].position;
    phongExp = scene.materials[result.m - 1].phong_exponent;
    kd = scene.materials[result.m - 1].diffuse;
    ks = scene.materials[result.m - 1].specular;
    km = scene.materials[result.m - 1].mirror;
    Vec3f mirrorColor = Vec3f(0, 0, 0);

    h = (v + wi).normalize();
    float Ir = pow((Lpos - result.hitPos).length(), 2);
    float calc = fmax(0, result.n*wi)/Ir;
    float calc2 = pow(fmax(0, result.n*h), phongExp)/Ir;

    if(recursion > 0) { //recursive reflectance
        w0 = pos - result.hitPos;
        wr = w0*-1 + result.n*(result.n*w0)*2;

        Check temp = findClosest(scene.vertex_data, result.hitPos, wr, scene.triangles, scene.spheres, scene.meshes, normals);
        Vec3f temp2 = findColor(result.hitPos, temp, amb, scene, wr, normals, --recursion);
        mirrorColor.x = mirrorColor.x + km.x*temp2.x;
        mirrorColor.y = mirrorColor.y + km.y*temp2.y;
        mirrorColor.z = mirrorColor.z + km.z*temp2.z;
    } 

    float x = kd.x*I.x*calc + ks.x*I.x*calc2 + mirrorColor.x;
    float y = kd.y*I.y*calc + ks.y*I.y*calc2 + mirrorColor.y;
    float z = kd.z*I.z*calc + ks.z*I.z*calc2 + mirrorColor.z;
    mirrorColor = Vec3f(0, 0, 0);
    return Vec3f(x, y, z);
}

Vec3f findColor(Vec3f pos, Check &result, Vec3f amb, Scene &scene, Vec3f d, vector<vector<Vec3f>> &normals, int recursion) { //For Reflectance
    Vec3f color = Vec3f(0, 0, 0);

    if (result.hit) {
        Vec3f v = (pos - result.hitPos).normalize();
        for (int i = 0; i < scene.point_lights.size(); ++i) {
            bool gotShadow = false;
            Vec3f wi = (scene.point_lights[i].position - result.hitPos).normalize();
            Vec3f wiNotNormalized = scene.point_lights[i].position - result.hitPos;
            for (int j = 0; j < scene.spheres.size(); ++j) {
                Sphere sphere = scene.spheres[j];
                if(!(result.type == "s" && result.sira == j) && result.t>0) {
                    if(shadowSphere(result.hitPos, wi, scene.vertex_data[sphere.center_vertex_id - 1], sphere.radius, scene.shadow_ray_epsilon, wiNotNormalized))
                        gotShadow = true;
                }
            }
            if(gotShadow == false) {
                for (int j2 = 0; j2 < scene.triangles.size(); ++j2) {
                    if(!(result.type == "t" && result.sira == j2) && result.t>0) {
                        if(shadowTriangle(scene.vertex_data, result.hitPos, wi, scene.triangles[j2].indices, scene.shadow_ray_epsilon, wiNotNormalized))
                            gotShadow = true;
                    }
                }
            }
            // if(gotShadow == false) {
            //     for (int j3 = 0; j3 < scene.meshes.size(); ++j3) {
            //         gotShadow = shadowMesh(scene.vertex_data, result.hitPos, wi, scene.meshes[j3].faces, scene.shadow_ray_epsilon, wiNotNormalized, normals, j3);
            //     }
            // }
            if(!gotShadow)
                color = color + applyEffects(scene, pos, result, v, wi, i, normals, amb, d, recursion);
        }
        color = color + amb;
    } else {
        color = color + Vec3f(scene.background_color.x, scene.background_color.y, scene.background_color.z);
    }
    return color;
}

int clamp(float f) {
    if (f > 255)
        return 255;
    else if (f < 0)
        return 0;
    return f;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - MAIN
int main(int argc, char* argv[]) {
    Scene scene;
    scene.loadFromXml(argv[1]);

    for (int z = 0; z < scene.cameras.size(); ++z) {
        int count = 0;
        Camera cam = scene.cameras[z];
        int width = cam.image_width;
        int height = cam.image_height;
        Vec3f w = cam.gaze.normalize()*-1;
        Vec3f u = cam.up.normalize()^w;
        cam.up = w^u;

        unsigned char* image = new unsigned char [width*height*3];
        vector<vector<Vec3f>> normals;
         // Calculate normals of meshes beforehand in order to not calculate anything for the triangles we cant see
        for (int i = 0; i < scene.meshes.size(); ++i) {
            vector<Vec3f> temp;
            normals.push_back(temp);
            for (int j = 0; j < scene.meshes[i].faces.size(); ++j) {
                Vec3f p1 = scene.vertex_data[scene.meshes[i].faces[j].v0_id - 1];
                Vec3f p2 = scene.vertex_data[scene.meshes[i].faces[j].v1_id - 1];
                Vec3f p3 = scene.vertex_data[scene.meshes[i].faces[j].v2_id - 1];
                Vec3f normal = ((p2 - p1)^(p3 - p1)).normalize();
                normals[i].push_back(normal);
            }
        }

        // Render for all pixels
        for (int i = 0; i < cam.image_height; ++i) {
            for (int j = 0; j < cam.image_width; ++j) {
                float left = cam.near_plane.x;
                float right = cam.near_plane.y;
                float bot = cam.near_plane.z;
                float top = cam.near_plane.w;

                float su = (right - left)*(j + 0.5)/cam.image_width;
                float sv = (top - bot)*(i + 0.5)/cam.image_height;
                Vec3f m = cam.position + cam.gaze*cam.near_distance;
                Vec3f q = m + u*left + cam.up*top;
                Vec3f s = q + u*su - cam.up*sv;
                Vec3f d = s - cam.position;

                Check dir = findClosest(scene.vertex_data, cam.position, d, scene.triangles, scene.spheres, scene.meshes, normals);

                Vec3f amb = ambient(scene.ambient_light, scene.materials, dir.m-1);
                
                Vec3f color = findColor(cam.position, dir, amb, scene, d, normals, scene.max_recursion_depth);

                for (int k = 0; k < 3; ++k) {
                    if(k==0)
                        image[count++] = clamp(color.x);
                    else if(k==1)
                        image[count++] = clamp(color.y);
                    else
                        image[count++] = clamp(color.z);
                }
            }
        }
        write_ppm(cam.image_name.c_str(), image, width, height);
    }
}