#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include "geometry.h"

#define SPHERE 0
#define SEGMENT 1

double dot(Vec2f u, Vec2f v) { return u.x * v.x + u.y * v.y; }
double dot(Vec3f u, Vec3f v) { return u.x * v.x + u.y * v.y + u.z * v.z; }

double dot2(Vec2f v) { return dot(v,v); }
double dot2(Vec3f v) { return dot(v,v); }

double clamp(double v, double min, double max) { return v < min ? min : (v > max ? max : v); }
template <typename T> inline T lerp(const T &v0, const T &v1, float t) {
    return v0 + (v1-v0)*std::max(0.f, std::min(1.f, t));
}

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

struct Primitive{ 

    Vec3f color;
    Primitive(const Vec3f &pColor) : color(pColor) {};
    virtual float signed_distance(const Vec3f &p) const{};
};

/*
SDF functions are from : https://iquilezles.org/articles/distfunctions/
*/
struct Sphere : Primitive {
    Vec3f center;
    float radius;

    Sphere(const Vec3f &pColor, const Vec3f &pCenter, const float &pRadius) : Primitive(pColor), 
        center(pCenter), radius(pRadius) {};
    float signed_distance(const Vec3f &p) const {
        return (p - center).norm() - radius;
    }    
};

struct Segment : Primitive{

    Vec3f a, b;
    float thickness;

    Segment(const Vec3f &pColor, const Vec3f &pA, const Vec3f &pB, const float &pThickness) : Primitive(pColor), 
        a(pA), b(pB), thickness(pThickness) {};
    float signed_distance(const Vec3f &p) const {
        Vec3f pa = p - a, ba = b - a;
        float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
        return ( pa - ba*h ).norm() - thickness;
    }
};

struct RoundCone : Primitive{

    Vec3f a, b;
    float r1, r2;

    RoundCone(const Vec3f &pColor, const Vec3f &pA, const Vec3f &pB, const float &pR1, const float &pR2) : Primitive(pColor), 
        a(pA), b(pB), r1(pR1), r2(pR2) {};
    float signed_distance(const Vec3f &p) const {
        Vec3f  ba = b - a;
        float l2 = dot(ba,ba);
        float rr = r1 - r2;
        float a2 = l2 - rr*rr;
        float il2 = 1.0/l2;

        Vec3f pa = p - a;
        float y = dot(pa,ba);
        float z = y - l2;
        float x2 = dot2( pa*l2 - ba*y );
        float y2 = y*y*l2;
        float z2 = z*z*l2;

        float k = sign(rr)*rr*rr*x2;
        if( sign(z)*a2*z2>k )
            return  sqrt(x2 + z2) * il2 - r2;
        if( sign(y)*a2*y2<k )
            return  sqrt(x2 + y2) * il2 - r1;

        return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
    }
};

struct CappedCone : Primitive{

    Vec3f pos;
    double h;
    double ra, rb;

    CappedCone(const Vec3f &pColor, const Vec3f &pPos, const double &pH, const double &pRa, const double &pRb) : Primitive(pColor), 
        pos(pPos), h(pH), ra(pRa), rb(pRb) {};

    float signed_distance(const Vec3f &p) const {
        
        Vec2f base = Vec2f(p.x - pos.x, p.z - pos.z);

        Vec2f q = Vec2f( base.norm(), p.y - pos.y );
        Vec2f k1 = Vec2f(rb, h);
        Vec2f k2 = Vec2f(rb-ra,2.0*h);
        Vec2f ca = Vec2f(q.x-std::min((double)q.x,(q.y<0.0)?ra:rb), abs(q.y)-h);
        Vec2f cb = q - k1 + k2*clamp( dot(k1-q,k2)/dot2(k2), 0.0, 1.0 );
        float s = (cb.x<0.0 && ca.y<0.0) ? -1.0 : 1.0;
        return s*sqrt( std::min(dot2(ca),dot2(cb)) );
    }
};

struct CappedCone2 : Primitive{

    Vec3f a, b;
    double ra, rb;

    CappedCone2(const Vec3f &pColor, const Vec3f &pA, const Vec3f &pB, const double &pRa, const double &pRb) : Primitive(pColor), 
        a(pA), b(pB), ra(pRa), rb(pRb) {};

    float signed_distance(const Vec3f &p) const {
        
        float rba  = rb-ra;
        float baba = dot(b-a,b-a);
        float papa = dot(p-a,p-a);
        float paba = dot(p-a,b-a)/baba;
        float x = sqrt( papa - paba*paba*baba );
        float cax = std::max(0.0,x-((paba<0.5)?ra:rb));
        float cay = abs(paba-0.5)-0.5;
        float k = rba*rba + baba;
        float f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );
        float cbx = x-ra - f*rba;
        float cby = paba - f;
        float s = (cbx<0.0 && cay<0.0) ? -1.0 : 1.0;
        return s*sqrt( std::min(cax*cax + cay*cay*baba,
                            cbx*cbx + cby*cby*baba) );
    }
};

struct Pyramid : Primitive{

    Vec3f pos;
    double h;

    Pyramid(const Vec3f &pColor, const Vec3f &pPos, const double &pH) : Primitive(pColor), 
        pos(pPos), h(pH) {};

    float signed_distance(const Vec3f &p) const {
        float m2 = h*h + 0.25;

        Vec3f base = p - pos;
        Vec3f c = Vec3f(abs(base.x), base.y, abs(base.z));
        Vec3f c2 = Vec3f(((c.z > c.x) ? c.z : c.x) - 0.5, c.y, ((c.z > c.x) ? c.x : c.z) - 0.5);

        Vec3f q = Vec3f( c2.z, h*c2.y - 0.5*c2.x, h*c2.x + 0.5*c2.y);
        
        float s = std::max((double)-q.x, 0.0);
        float t = clamp( (q.y-0.5*q.x)/(m2+0.25), 0.0, 1.0 );
            
        float a = m2*(q.x+s)*(q.x+s) + q.y*q.y;
        float b = m2*(q.x+0.5*t)*(q.x+0.5*t) + (q.y-m2*t)*(q.y-m2*t);
            
        float d2 = std::max((double)-q.y, q.x*m2+q.y*0.5) < 0.0 ? 0.0 : std::min(a,b);
            
        return sqrt( (d2+q.z*q.z)/m2 ) * sign(std::max(q.z,-c2.y));
    }
};

bool trace(const Primitive *pr, const Vec3f &orig, const Vec3f &dir, Vec3f &pos) {
    pos = orig;
    for (size_t i=0; i<256; i++) {
        float d = pr->signed_distance(pos);
        if (d < 0) return true;
        pos = pos + dir*std::max(d*0.1f, .01f);
    }
    return false;
}

float elevation(const float &x, const float &z){

    return 0.04 * (sin(2 * x) * sin (2* z)) + 0.03*(cos(2.5 * x) + cos(2.5*z)) -0.9;
}

bool intersect_terrain(const Vec3f &orig, const Vec3f &dir, Vec3f &pos){

    float dt = 0.01f;
    float mint = 0.001f;
    float maxt = 30.0f;

    for (float t=mint; t<maxt; t+=dt) {

        Vec3f p = orig + dir * t;
        if(p.y < elevation(p.x,p.z)){
            pos = orig + dir * (t - 0.5* dt);
            return true;
        }
    }
    return false;
}

Vec3f terrain_normal( const Vec3f & p )
{
    const float eps = 0.1;
    return Vec3f( elevation(p.x-eps,p.z) - elevation(p.x+eps,p.z), 2.0f*eps, elevation(p.x,p.z-eps) - elevation(p.x,p.z+eps) ).normalize();
}

Vec3f distance_field_normal(const Primitive *pr, const Vec3f &pos) {
    const float eps = 0.1;
    float d = pr->signed_distance(pos);
    float nx = pr->signed_distance(pos + Vec3f(eps, 0, 0)) - d;
    float ny = pr->signed_distance(pos + Vec3f(0, eps, 0)) - d;
    float nz = pr->signed_distance(pos + Vec3f(0, 0, eps)) - d;
    return Vec3f(nx, ny, nz).normalize();
}

float shadow(const Primitive *pr, const Vec3f &orig, const Vec3f &dir, const float k){

    float res = 1.0;
    float t = 0;
    for(int i=0; i < 256 && t < 50; i++){

        float h = pr->signed_distance(orig + dir * t);
        if(h < 0.001){
            return 0;
        }
        res = std::min(res, k*h/t);
        t += h;
    }
    return res;
}

int main() {
    const int   width    = 800;
    const int   height   = 600;
    const float fov      = M_PI/3.;
    std::vector<Vec3f> framebuffer(width*height);
    Vec3f camPos = Vec3f(4, 4.5, 3);
    Vec3f sun_light_dir = Vec3f(1,0.9,0.5).normalize();
    Vec3f ind_light_dir = Vec3f(-sun_light_dir.x, 0, -sun_light_dir.z);
    float a = -M_PI/7.;
    float b = M_PI/5.;
    float c = 0;

    Vec3f snow_color = Vec3f(1,1,1);
    Vec3f sky_color = Vec3f(0.4,0.57,0.69);
    Vec3f leaves_color = Vec3f(0.033,0.084,0.041);
    Vec3f wood_color = Vec3f(0.125,0.068,0.032);
    Vec3f stone_color = Vec3f(0.1,0.1,0.1);
    Vec3f carrot_color = Vec3f(1,0.5,0);
    Vec3f hat_color = Vec3f(0,0,0);

    // Pre-compute to be use camera rotation matrix
    float cosA = cos(a);
    float sinA = sin(a);
    float cosB = cos(b);
    float sinB = sin(b);
    float cosC = cos(c);
    float sinC = sin(c);

    std::srand(1000);

    std::vector<Primitive*> scene = {
        
        // Body
        new Sphere(snow_color, {0,0,-3}, 1.5), 
        new Sphere(snow_color, {0,2.05,-3}, 1.1),
        new Sphere(snow_color, {0,3.4,-3}, 0.75),

        // Buttons
        new Sphere(stone_color, {0,1.95,-1.93}, 0.1),
        new Sphere(stone_color, {0,2.3,-1.92}, 0.1),
        new Sphere(stone_color, {0,2.65,-2.03}, 0.1),
        
        // Carrot
        new RoundCone(carrot_color, {0,3.45,-3}, {0,3.45,-2}, 0.2, 0.08),

        // Eyes
        new Sphere(stone_color, {-0.3,3.6,-2.35}, 0.1),
        new Sphere(stone_color, {0.3,3.6,-2.35}, 0.1),

        // Hat
        new CappedCone2(hat_color, {0,4.2,-3}, {0,4.6,-3}, 0.6, 0.65),
        new CappedCone2(hat_color, {0,4,-3}, {0,4.03,-3}, 0.75, 0.75),

        // Arms
        new Segment(wood_color, {1,2.05,-3}, {2.4,3.15,-3}, 0.07),
        new Segment(wood_color, {2.05,2.9,-3}, {2.07,3.3,-3}, 0.06),

        new Segment(wood_color, {-1,2.05,-3}, {-2.4,3.15,-3}, 0.07),
        new Segment(wood_color, {-2.05,2.9,-3}, {-2.07,3.3,-3}, 0.06),

        // Trees
        new CappedCone(wood_color, {-7,0,-4}, 3, 0.85, 0.4),
        new CappedCone(leaves_color, {-7,4,-4}, 1.5, 2, 1),
        new CappedCone(leaves_color, {-7,6,-4}, 1.5, 1.8, 0.8),
        new CappedCone(leaves_color, {-7,8,-4}, 1.5, 1.6, 0),

        new CappedCone(wood_color, {-5,0,-7}, 3, 0.8, 0.4),
        new CappedCone(leaves_color, {-5,4,-7}, 1.5, 2.2, 1),
        new CappedCone(leaves_color, {-5,6,-7}, 1.5, 2, 0.8),
        new CappedCone(leaves_color, {-5,8,-7}, 1.5, 1.8, 0),

        new CappedCone(wood_color, {-9,0,0}, 3, 0.8, 0.4),
        new CappedCone(leaves_color, {-9,4,0}, 1.5, 2.2, 1),
        new CappedCone(leaves_color, {-9,6,0}, 1.5, 2, 0.8),
        new CappedCone(leaves_color, {-9,8,0}, 1.5, 1.8, 0),
        new CappedCone(leaves_color, {-9,9.3,0}, 1.2, 1.6, 0),

        new CappedCone(wood_color, {0,0,-11}, 3, 0.8, 0.4),
        new CappedCone(leaves_color, {0,4,-11}, 1.5, 2.2, 1),
        new CappedCone(leaves_color, {0,6,-11}, 1.5, 2, 0.8),
        new CappedCone(leaves_color, {0,8,-11}, 1.5, 1.8, 0),
        new CappedCone(leaves_color, {0,9.3,-11}, 1.2, 1.6, 0),

        new CappedCone(wood_color, {4,0,-8}, 3, 0.85, 0.4),
        new CappedCone(leaves_color, {4,4,-8}, 1.5, 2.2, 1),
        new CappedCone(leaves_color, {4,6,-8}, 1.5, 2, 0.8),
        new CappedCone(leaves_color, {4,8,-8}, 1.5, 1.8, 0),

        new CappedCone(wood_color, {-13,0,3}, 3, 0.85, 0.4),
        new CappedCone(leaves_color, {-13,4,3}, 1.5, 2.2, 1),
        new CappedCone(leaves_color, {-13,6,3}, 1.5, 2, 0.8),
        new CappedCone(leaves_color, {-13,8,3}, 1.5, 1.8, 0),

        new CappedCone(wood_color, {8.5,0,-7}, 3, 0.8, 0.4),
        new CappedCone(leaves_color, {8.5,4,-7}, 1.5, 2.2, 1),
        new CappedCone(leaves_color, {8.5,6,-7}, 1.5, 2, 0.8),
        new CappedCone(leaves_color, {8.5,8,-7}, 1.5, 1.8, 0),
        new CappedCone(leaves_color, {8.5,9.3,-7}, 1.2, 1.6, 0),

        new Segment(wood_color, {-1.5,-0.8,0.5}, {-5,-0.8,1.8}, 0.09),
        new Segment(wood_color, {-2.5,-0.8,1}, {-1.3,-0.8,1.7}, 0.07),

        // Rocks
        new Sphere(stone_color, {5.1,-0.8,-2.4}, 0.4),
        new Sphere(stone_color, {5.6,-0.8,-2.5}, 0.2),
    };

#pragma omp parallel for
    for (size_t j = 0; j<height; j++) { // actual rendering loop
        for (size_t i = 0; i<width; i++) {
            float dir_x = (i + 0.5) -  width/2.;
            float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
            float dir_z = -0.5*height/(2.*tan(fov/2.));

            // Camera rotation
            Vec3f xRot = Vec3f(dir_x, dir_y * cosA - dir_z * sinA, dir_y * sinA + dir_z * cosA);
            Vec3f yRot = Vec3f(xRot.x * cosB + xRot.z * sinB, xRot.y, -xRot.x * sinB + xRot.z * cosB);
            Vec3f dir = Vec3f(yRot.x * cosC - yRot.y * sinC, yRot.x * sinC + yRot.y * cosC, yRot.z).normalize();

            Vec3f hit;
            bool intersect = false;
            float nearest_dist = 1e10;
            float shadow_value = 1;
            float sun_intensity = 1;
            float sky_intensity = 1;
            float indirect_intensity = 1;

            // Terrain raymarching
            if(intersect_terrain(camPos, dir, hit)){
                framebuffer[i+j*width] = snow_color * 0.5;

                Vec3f normal = terrain_normal(hit);
                sun_intensity = clamp(sun_light_dir*normal, 0, 1);
                sky_intensity = clamp(0.5 + normal.y * 0.5, 0, 1);
                indirect_intensity = clamp(normal * ind_light_dir.normalize(), 0.0, 1.0 );
                
                // Compute shadows
                for(const Primitive *p : scene){
                    shadow_value *= shadow(p, hit, sun_light_dir, 16);
                    if(shadow_value < 0.001){
                        break;
                    }
                }
                nearest_dist = (hit - camPos).norm();
                intersect = true;
            }

            // Primitives raymarching
            for(const Primitive *p : scene){
                if (trace(p, camPos, dir, hit)) {
                    float dist = (hit - camPos).norm();
                    if(dist < nearest_dist){
                        Vec3f normal = distance_field_normal(p, hit);
                        sun_intensity = clamp(sun_light_dir*normal, 0, 1);
                        sky_intensity = clamp(0.5 + normal.y * 0.5, 0, 1);
                        indirect_intensity = clamp(normal * ind_light_dir.normalize(), 0.0, 1.0 );
                        
                        // Compute shadows
                        for(const Primitive *ps : scene){
                            if(ps != p){
                                shadow_value = shadow(ps, hit, sun_light_dir, 16);
                                if(shadow_value < 0.1)
                                    break;
                            }
                        }
                        framebuffer[i+j*width] = p->color * 0.3;
                                                        
                        nearest_dist = dist;
                    }
                    intersect = true;
                }
            }
            if(!intersect){
                framebuffer[i+j*width] = sky_color;
            }else{

                // Compute lighting

                // Shadow colorization
                Vec3f cShadow = Vec3f(shadow_value, pow(shadow_value, 1.2), pow(shadow_value, 1.5));

                // Sun light
                Vec3f val = Vec3f(framebuffer[i+j*width].x * sun_intensity * 1.64 * cShadow.x,
                                    framebuffer[i+j*width].y * sun_intensity * 1.27 * cShadow.y,
                                    framebuffer[i+j*width].z * sun_intensity * 0.99 * cShadow.z) * 1.5;
                
                // Sky light
                val = val + Vec3f(0.16,0.20,0.28) * sky_intensity * 0.2;

                // Indirect light
                val = val + Vec3f(0.40,0.28,0.20) * indirect_intensity * 0.3;

                framebuffer[i+j*width] = val;
            }

            // gamma correction
            framebuffer[i+j*width] = Vec3f(
                                        pow(framebuffer[i+j*width].x, 1.0/2.2), 
                                        pow(framebuffer[i+j*width].y, 1.0/2.2),
                                        pow(framebuffer[i+j*width].z, 1.0/2.2));
            
        }
    }

    std::ofstream ofs("./out.ppm", std::ios::binary); // save the framebuffer to file
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(std::max(0, std::min(255, static_cast<int>(255*framebuffer[i][j]))));
        }
    }
    ofs.close();

    return 0;
}