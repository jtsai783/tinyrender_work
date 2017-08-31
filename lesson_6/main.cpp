#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include <iostream>
#include "model.h"
#include "matrix.h"
#include "our_gl.h"
// #include <algorithm>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
TGAImage *texture_img = NULL;
const int width  = 800;
const int height = 800;

Vec3f camera(1, 1, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);
Vec3f light_dir(1,1,1);

struct Shader : public IShader {
    // Vec3f varying_intensity;
    int thisFace;
    Matrix cob;
    Matrix cob_IT;

    virtual Vec3f vertex(int nface, int nthvert){
        int normal_index = (model->face(nface))[nthvert * 3 + 2];
        int vert_index = (model->face(nface))[nthvert * 3];
        Vec3f v = model->vert(vert_index);
        // varying_intensity[nthvert] = std::max(0.f, model->normal_vert(normal_index) * light_dir);
        Matrix v_m = Matrix(4,1);
        v_m(0, 0) = v.x;
        v_m(1, 0) = v.y;
        v_m(2, 0) = v.z;
        v_m(3, 0) = 1;
        v_m = cob * v_m;
        v_m /= v_m(3,0);
        thisFace = nface;
        return Vec3f(v_m(0,0),v_m(1,0),v_m(2,0));
    }

    virtual bool fragment(float *bc, TGAColor &color){

        std::vector<int> face = model->face(thisFace);
        Vec3f texture_coords_0 = model->texture_vert(face[1]);
        Vec3f texture_coords_1 = model->texture_vert(face[4]);
        Vec3f texture_coords_2 = model->texture_vert(face[7]);

        float texture_x = bc[0] * texture_coords_0.x + bc[1] * texture_coords_1.x + bc[2] * texture_coords_2.x;
        float texture_y = bc[0] * texture_coords_0.y + bc[1] * texture_coords_1.y + bc[2] * texture_coords_2.y;

        // float intensity = bc[0] * varying_intensity[0] + bc[1] * varying_intensity[1] + bc[2] * varying_intensity[2];

        // intensity = 1;

        // if( intensity >= 0 && intensity < 0.25){
        //     intensity = 0.25;
        // } else if ( intensity >= 0.25 && intensity < 0.5){
        //     intensity = 0.5;
        // } else if( intensity >= 0.5 && intensity < 0.75){
        //     intensity = 0.75;
        // } else if( intensity >= 0.75 && intensity < 1){
        //     intensity = 1;
        // }

        texture_x = texture_x * (model->diffuse).get_width();
        texture_y = texture_y * (model->diffuse).get_height();


        // //get the color on normal map
        TGAColor normal = (model->normal_map).get(roundf(texture_x), roundf(texture_y));
        Vec3f normal_vec = Vec3f((float)normal.r, (float)normal.g, (float)normal.b);

        Matrix vn_m = Matrix(4,1);
        vn_m(0,0) = normal_vec.x;
        vn_m(1,0) = normal_vec.y;
        vn_m(2,0) = normal_vec.z;
        vn_m(3,0) = 1;
        vn_m = cob_IT * vn_m;
        vn_m = vn_m / vn_m(3,0);
        normal_vec = Vec3f(vn_m(0,0),vn_m(1,0),vn_m(2,0));
        // normal_vec.normalize();

        Matrix light_m = Matrix(4,1);
        light_m(0,0) = light_dir.x;
        light_m(1,0) = light_dir.y;
        light_m(2,0) = light_dir.z;
        light_m(3,0) = 1;
        light_m = cob * light_m;
        light_m = light_m / light_m(3,0);
        Vec3f light_vec = Vec3f(light_m(0,0),light_m(1,0),light_m(2,0));
        // light_vec.normalize();

        float intensity = light_vec.normalize() * normal_vec.normalize();
        intensity = std::max(0.f, -intensity);
       

        color = (model->diffuse).get(roundf(texture_x), roundf(texture_y));
        color = TGAColor((float)color.r * intensity, (float)color.g * intensity, (float)color.b * intensity, 255);

        return false;
    }
};

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    // TGAImage texture_img(width, height, TGAImage::RGB);

    // texture_img.read_tga_file("african_head_diffuse.tga");
    // texture_img.flip_vertically();

    TGAImage image(width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    view(center, camera, up);
    proj(-1.0/(camera - center).norm());
    clip(width/8, height/8, width*3/4, height*3/4);
    light_dir.normalize();


    // int *zbuffer = new int[width*height];
    // for (int i=width*height; i--; zbuffer[i] = 0);


    Shader shader;
    shader.cob = Viewport * Projection * ModelView;
    shader.cob_IT = shader.cob.inverse().transpose();
    for (int i= 0; i < model->nfaces(); i++) { 
        std::vector<int> face = model->face(i);
        Vec3f screen_coords[3]; 
        for (int j=0; j<3; j++) { 
            screen_coords[j] = shader.vertex(i, j);
        }
        triangle(zbuffer, screen_coords, image, shader); 
    }


    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    zbuffer.flip_vertically();
    image.write_tga_file("output.tga");
    zbuffer.write_tga_file("zbuffer.tga");

    delete model;
    return 0;
}

