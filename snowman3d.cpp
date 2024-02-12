#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "geometry.h"

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;

struct Light {
    Light(const Vec3f &p, const float i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Material {
    float (*displacement)(Vec3f&);

    Material(const float r, const Vec4f &a, const Vec3f &color, const float spec, float (*d)(Vec3f&)) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec), displacement(d) {}
    Material() : refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};

struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float r, const Material &m) : center(c), radius(r), material(m) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t, const float eta_i=1.f) { // Snell's law
    float cosi = - std::max(-1.f, std::min(1.f, I*N));
    if (cosi<0) return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k<0 ? Vec3f(1,0,0) : I*eta + N*(eta*cosi - sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

float signed_distance(const Vec3f &p, Sphere sphere) { // this function defines the implicit surface we render
    Vec3f s = Vec3f(p).normalize(sphere.radius);
    return std::sqrt(p.x*p.x+p.y*p.y+p.z*p.z) - (sphere.radius + sphere.material.displacement(s));
}

Vec3f distance_field_normal(const Vec3f &pos, Sphere s) { // simple finite differences, very sensitive to the choice of the eps constant
    const float eps = 0.1;
    float d = signed_distance(pos, s);
    float nx = signed_distance(pos + Vec3f(eps, 0, 0), s) - d;
    float ny = signed_distance(pos + Vec3f(0, eps, 0), s) - d;
    float nz = signed_distance(pos + Vec3f(0, 0, eps), s) - d;
    return Vec3f(nx, ny, nz).normalize();
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = distance_field_normal(hit - spheres[i].center, spheres[i]);
            material = spheres[i].material;
        }
    }

    float checkerboard_dist = std::numeric_limits<float>::max();
    if(fabs(dir.y)>1e-3){
        float d = -(orig.y+8)/dir.y;
        Vec3f pt = orig + dir*d;
        if (d>0 && pt.x < 14 &&  pt.x > -15 && pt.z>8 && pt.z<32 && d<spheres_dist) {
            Vec3f pt = orig + dir*d;
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0,1,0);
            material.diffuse_color = Vec3f(.3, .3, .3);
        }
    }

    return std::min(spheres_dist, checkerboard_dist)<1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth=0) {
    Vec3f point, N;
    Material material;

    if (depth>4 || !scene_intersect(orig, dir, spheres, point, N, material)) {
        
        float phi = atan2(dir.y, sqrtf(dir.x * dir.x + dir.z * dir.z));
        float theta = atan2(dir.z, dir.x);

        int row = ( 1 - (phi / (M_PI/2)) ) * (envmap_height/2);
        int col = ( 1 - (-theta / (M_PI)) ) * (envmap_width/2);

        return envmap[col + row * envmap_width];
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // offset the original point to avoid occlusion by the object itself
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, depth + 1);
    Vec3f refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflect_color*material.albedo[2] + refract_color*material.albedo[3];
}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    const int   width    = 1024;
    const int   height   = 768;
    const float fov      = M_PI/3.;
    std::vector<Vec3f> framebuffer(width*height);

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) { // actual rendering loop
        for (size_t i = 0; i<width; i++) {
            float dir_x =  (i + 0.5) -  3*width/5.;
            float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
            float dir_z = height/(2.*tan(fov/2.));
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), Vec3f(dir_x+50, dir_y-200, dir_z).normalize(), spheres, lights);
        }
    }

    std::vector<unsigned char> pixmap(width*height*3);
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            pixmap[i*3+j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    stbi_write_jpg("out.jpg", width, height, 3, pixmap.data(), 100);
}

float displacementSnow(Vec3f &s){
    return sin(16*s.x)*sin(16*s.y)*sin(16*s.z)*0.005;
}

float displacementEye(Vec3f &s){
    return sin(16*s.x)*sin(16*s.y)*sin(16*s.z)*0.005;
}

float displacementCarrot(Vec3f &s){
    return Vec2f(sin(75.0), cos(75.0)) * Vec2f(Vec2f(s.x,s.y).norm() , s.z);
}

int main() {
    int n = -1;
    unsigned char *pixmap = stbi_load("../snow.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    envmap = std::vector<Vec3f>(envmap_width*envmap_height);
    for (int j = envmap_height-1; j>=0 ; j--) {
        for (int i = 0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);

    Material      carrot(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(1.0, 127.0/255.0, 0.01), 50., displacementCarrot);
    Material     snow(1.0, Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(1.0, 1.0, 1.0), 10., displacementSnow);
    Material     eye(1.0, Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(0.0, 0.0, 0.0), 10., displacementEye);

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-1, -5, 12), 2, snow));
    spheres.push_back(Sphere(Vec3f(-1, -4, 12), 1.8, snow));
    spheres.push_back(Sphere(Vec3f(-1, -3, 12), 1.6, snow));
    spheres.push_back(Sphere(Vec3f(-1, -1, 12), 1, snow));

    spheres.push_back(Sphere(Vec3f(-1, -0.7, 11), 0.165, carrot));

    spheres.push_back(Sphere(Vec3f(-1.40, -0.41, 11.3), 0.115, eye));
    spheres.push_back(Sphere(Vec3f(-0.6, -0.41, 11.3), 0.115, eye));

    std::vector<Light>  lights;
    lights.push_back(Light(Vec3f(50, 60, -120), 0.5));
    lights.push_back(Light(Vec3f(0, 240, -120), 0.5));
    lights.push_back(Light(Vec3f(-50, 60, -120), 0.5));

    render(spheres, lights);

    return 0;
}

