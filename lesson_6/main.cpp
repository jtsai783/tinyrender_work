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

Vec3f camera(0, -1, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);
Vec3f light_dir(0,0,-1);

struct GouraudShader : public IShader {
    Vec3f varying_intensity;

    virtual Matrix vertex(Vec3f v, int nthvert){
        varying_intensity[nthvert] = 1;
        Matrix v_m = Matrix(4,1);
        v_m(0, 0) = v.x;
        v_m(1, 0) = v.y;
        v_m(2, 0) = v.z;
        v_m(3, 0) = 1;
        v_m = Viewport * Projection * ModelView * v_m;
        // v_m /= v_m(3,0);
        // std::cout << v_m;
        return v_m;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color){

        return false;
    }
};

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

    view(center, camera, up);
    proj(-1.0/(camera - center).norm());
    clip(width/8, height/8, width*3/4, height*3/4);
    light_dir.normalize();


    int *zbuffer = new int[width*height];
    for (int i=width*height; i--; zbuffer[i] = 0);


    GouraudShader shader;
    for (int i= 0; i < model->nfaces(); i++) { 
        std::vector<int> face = model->face(i);
        // std::vector<int> texture_face = model->texture_face(i);
        Matrix screen_coords[3]; 
        // Vec3f texture_coords[3];
        // Vec3f normal_coords[3];
        for (int j=0; j<3; j++) { 

            screen_coords[j] = shader.vertex(model->vert(face[j]), j);

            // Vec3f v = model->vert(face[j]); 




            // Matrix screen =  cob * v_m;
            // screen /= screen(3,0);
            // Vec3f vt = model->texture_vert(texture_face[j]);
            // Vec3f vn = model->normal_vert(face[j]);
            // screen_coords[j] = Vec3f(screen(0,0), screen(1,0), screen(2,0));
            // Matrix vn_m = Matrix(4,1);
            // vn_m(0, 0) = vn.x;
            // vn_m(1, 0) = vn.y;
            // vn_m(2, 0) = vn.z;
            // vn_m(3, 0) = 0;
            // vn_m = norm_cob * vn_m;
            // normal_coords[j] = Vec3f(vn_m(0,0),vn_m(1,0),vn_m(2,0));
            // texture_coords[j] = vt;
        }
        triangle(zbuffer, screen_coords, image, shader, width); 
    }


    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    delete model;
    return 0;
}

