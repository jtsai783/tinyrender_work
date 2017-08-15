#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include <iostream>
#include "model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(p0.x-p1.x)<std::abs(p0.y-p1.y)) {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    if (p0.x>p1.x) {
        std::swap(p0, p1);
    }

    for (int x=p0.x; x<=p1.x; x++) {
        float t = (x-p0.x)/(float)(p1.x-p0.x);
        int y = p0.y*(1.-t) + p1.y*t;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

void bary(int x, int y, Vec3f *pts, float *bc){
    Vec3f v0 = pts[1] - pts[0];
    Vec3f v1 = pts[2] - pts[0];
    Vec3f P = Vec3f(x,y, 0);
    Vec3f v2 = P - pts[0];

    float dot00 = v0.x * v0.x + v0.y * v0.y;
    float dot01 = v0.x * v1.x + v0.y * v1.y;
    float dot02 = v0.x * v2.x + v0.y * v2.y;
    float dot11 = v1.x * v1.x + v1.y * v1.y;
    float dot12 = v1.x * v2.x + v1.y * v2.y;

    float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    bc[1] = u;
    bc[2] = v;
    bc[0] = 1 - (u + v);
    // return (u >= 0) && (v >= 0) && (u + v <= 1); 
}

void triangle(float *zbuffer, Vec3f *pts, TGAImage &image, TGAColor color) {
    //find bounding box
    int y_max = pts[0].y;
    int y_min = pts[0].y;
    int x_max = pts[0].x;
    int x_min = pts[0].x;
    for(int i = 0; i < 3; i++){
        if(y_max < pts[i].y){
            y_max = pts[i].y;
        }
        if(y_min > pts[i].y){
            y_min = pts[i].y;
        }
        if(x_max < pts[i].x){
            x_max = pts[i].x;
        }
        if(x_min > pts[i].x){
            x_min = pts[i].x;
        }
    }


    //clip bounding box with screen rectangle
    if(y_max < 0){ y_max = 0;};
    if(y_min < 0){ y_min = 0;};
    if(x_max < 0){ x_max = 0;};
    if(x_min < 0){ x_min = 0;};

    if(y_max > image.get_height() - 1){ y_max = image.get_height() - 1;};
    if(y_min > image.get_height() - 1){ y_min = image.get_height() - 1;};
    if(x_max > image.get_width() - 1){ x_max = image.get_width() - 1;};
    if(x_min > image.get_width() - 1){ x_min = image.get_width() - 1;};

    //for each point in the bounding box, test to see if its in the triangle
    for(int x = x_min; x <= x_max; x++){
        for(int y = y_min; y <= y_max; y++){
            float bc[3];
            bary(x, y, pts, bc);
            //calculate z-value of the point


            if((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)){
                // std::cout << pts[0];
                float z = bc[0] * (float)pts[0].z + bc[1] * (float)pts[1].z + bc[2] * (float)pts[2].z;
                    // std::cout << bc[0];
                    // std::cout << " ";
                    // std::cout << bc[1];
                    // std::cout << " ";
                    // std::cout << bc[2];
                    // std::cout << "\n";
                    // std::cout << z;
                    // std::cout << "\n";
                if(zbuffer[x + y * width] < z){
                    zbuffer[x + y * width] = z;
                    image.set(x,y,color);    
                }
                
            }
        }
    }




    // line(pts[0], pts[1], image, green);
    // line(pts[1], pts[2], image, green);
    // line(pts[2], pts[0], image, green);
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);

    // Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)};
    // Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)};
    // Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)};

    // triangle(t0, image, red);
    // triangle(t1, image, white);
    // triangle(t2, image, green);

    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    Vec3f light_dir(0,0,-1);
    for (int i=0; i<model->nfaces(); i++) { 
        std::vector<int> face = model->face(i); 
        Vec3f screen_coords[3]; 
        Vec3f world_coords[3]; 
        for (int j=0; j<3; j++) { 
            Vec3f v = model->vert(face[j]); 
            screen_coords[j] = Vec3f((v.x+1.)*width/2., (v.y+1.)*height/2., v.z); 
            world_coords[j]  = v; 
        } 
        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]); 
        n.normalize(); 
        float intensity = n*light_dir; 
        // if (intensity>0) { 
            triangle(zbuffer, screen_coords, image, TGAColor(intensity*255, intensity*255, intensity*255, 255)); 
        // } 
    }


    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}

