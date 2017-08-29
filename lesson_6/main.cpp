#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include <iostream>
#include "model.h"
#include "matrix.h"
#include "our_gl.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
TGAImage *texture_img = NULL;
const int width  = 800;
const int height = 800;

Vec3f cross(Vec3f a, Vec3f b){
    float x = a.y * b.z - a.z * b.y;
    float y = a.z * b.x - a.x * b.z;
    float z = a.x * b.y - a.y * b.x;
    return Vec3f(x,y,z);
}

Vec3f normalize(Vec3f a){
    double mag = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    return Vec3f(a.x / mag, a.y / mag, a.z / mag);
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    TGAImage texture_img(width, height, TGAImage::RGB);

    texture_img.read_tga_file("african_head_diffuse.tga");
    texture_img.flip_vertically();

    TGAImage image(width, height, TGAImage::RGB);

    Vec3f camera(1, 1, 3);
    Vec3f center(0, 0, 0);
    Vec3f up(0, 1, 0);

    Matrix view_m = view(center, camera, up);

    Matrix model_m = Matrix::createIdentity(4);

    Matrix proj_m = Matrix::createIdentity(4);
    proj_m(3, 2) = -1.0/(camera - center).norm();

    Matrix clip_m = clip(width/8, height/8, width*3/4, height*3/4);

    Matrix cob = clip_m * proj_m * view_m * model_m;
    Matrix norm_cob = cob.inverse().transpose();

    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    Vec3f light_dir(0,2,-1);

    for (int i= 0; i<model->nfaces(); i++) { 
        std::vector<int> face = model->face(i);
        std::vector<int> texture_face = model->texture_face(i);
        Vec3f screen_coords[3]; 
        Vec3f texture_coords[3];
        Vec3f normal_coords[3];
        for (int j=0; j<3; j++) { 
            Vec3f v = model->vert(face[j]);
            Matrix v_m = Matrix(4,1);
            v_m(0, 0) = v.x;
            v_m(1, 0) = v.y;
            v_m(2, 0) = v.z;
            v_m(3, 0) = 1;
            Matrix screen =  cob * v_m;
            screen /= screen(3,0);
            Vec3f vt = model->texture_vert(texture_face[j]);
            Vec3f vn = model->normal_vert(face[j]);
            screen_coords[j] = Vec3f(screen(0,0), screen(1,0), screen(2,0));
            Matrix vn_m = Matrix(4,1);
            vn_m(0, 0) = vn.x;
            vn_m(1, 0) = vn.y;
            vn_m(2, 0) = vn.z;
            vn_m(3, 0) = 0;
            vn_m = norm_cob * vn_m;
            normal_coords[j] = Vec3f(vn_m(0,0),vn_m(1,0),vn_m(2,0));
            texture_coords[j] = vt;
        }
        triangle(texture_coords, zbuffer, screen_coords, image, texture_img, light_dir, normal_coords); 
    }


    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    delete model;
    return 0;
}

