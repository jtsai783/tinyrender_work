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

const Vec3f camera(1, 1, 3);
const Vec3f center(0, 0, 0);
const Vec3f up(0, 1, 0);
const Vec3f light_dir(1,1,1);

Vec3f cross2(Vec3f a, Vec3f b){
    float x = a.y * b.z - a.z * b.y;
    float y = a.z * b.x - a.x * b.z;
    float z = a.x * b.y - a.y * b.x;
    return Vec3f(x,y,z);
}

Vec3f transform_v(Matrix m, Vec3f v){
    Matrix v_m = Matrix(4,1);
    v_m(0,0) = v.x;
    v_m(1,0) = v.y;
    v_m(2,0) = v.z;
    v_m(3,0) = 0.0;
    v_m = m * v_m;
    return Vec3f(v_m(0,0),v_m(1,0),v_m(2,0));
}

Vec3f transform_p(Matrix m, Vec3f v){
    Matrix v_m = Matrix(4,1);
    v_m(0,0) = v.x;
    v_m(1,0) = v.y;
    v_m(2,0) = v.z;
    v_m(3,0) = 1.0;
    v_m = m * v_m;
    v_m = v_m / v_m(3,0);
    return Vec3f(v_m(0,0),v_m(1,0),v_m(2,0));
}

struct Shader : public IShader {
    Vec3f texture_coords[3];
    // Vec3f normal_vec[3];
    Matrix cob;
    Matrix cob_IT;

    Matrix M;
    Matrix M_IT;

    Vec3f triangle_normal;
    Vec3f verts[3];

    virtual Vec3f vertex(int nface, int nthvert){
        std::vector<int> face = model->face(nface);
        Vec3f v = model->vert(face[nthvert * 3]);
        Vec3f thisVert = transform_p(cob, v);
        texture_coords[nthvert] = model->texture_vert(face[nthvert * 3 + 1]);
        verts[nthvert] = v;
        if(nthvert == 2){
           triangle_normal = cross2((verts[0] - verts[1]), (verts[0] - verts[2]));
        }
        // normal_vec[nthvert] = model->normal_vert(face[nthvert * 3 + 2]);
        // normal_vec[nthvert] = transform_v(M_IT, normal_vec[nthvert]);
        return thisVert;
    }

    virtual bool fragment(float *bc, TGAColor &color){

        //get color
        float texture_x = bc[0] * texture_coords[0].x + bc[1] * texture_coords[1].x + bc[2] * texture_coords[2].x;
        float texture_y = bc[0] * texture_coords[0].y + bc[1] * texture_coords[1].y + bc[2] * texture_coords[2].y;

        texture_x = texture_x * (model->diffuse).get_width();
        texture_y = texture_y * (model->diffuse).get_height();
        color = (model->diffuse).get(roundf(texture_x), roundf(texture_y));

        //interpolate the vectors
        // float normal_x = bc[0] * normal_vec[0].x + bc[1] * normal_vec[1].x + bc[2] * normal_vec[2].x;
        // float normal_y = bc[0] * normal_vec[0].y + bc[1] * normal_vec[1].y + bc[2] * normal_vec[2].y;
        // float normal_z = bc[0] * normal_vec[0].z + bc[1] * normal_vec[1].z + bc[2] * normal_vec[2].z;
        // Vec3f normal_vec = Vec3f(normal_x, normal_y, normal_z);

        // normal_x = ((normal_vec.normalize() + Vec3f(1.0,1.0,1.0)) * 0.5 * 255).x;
        // normal_y = ((normal_vec.normalize() + Vec3f(1.0,1.0,1.0)) * 0.5 * 255).y;
        // normal_z = ((normal_vec.normalize() + Vec3f(1.0,1.0,1.0)) * 0.5 * 255).z;

        // TGAColor normal_color = TGAColor(roundf(normal_x), roundf(normal_y), roundf(normal_z), 255);

        // color = normal_color;

        TGAColor normal_color = (model->normal_map).get(roundf(texture_x), roundf(texture_y));
        Vec3f normal_vec = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        normal_vec = normal_vec * (1.0/255.0) * 2.0 - Vec3f(1.0,1.0,1.0);

        Vec3f temp_light = transform_v(M, light_dir);
        Vec3f temp_norm = transform_v(M_IT, normal_vec);

        // normal_vec = transform_p((Viewport * ModelView).inverse().transpose(), normal_vec);

        // light_dir = transform_v(cob, light_dir);
        // normal_vec = transform_v(cob_IT, normal_vec);

        // float intensity = std::max(0.f,  light_dir.normalize() * normal_vec.normalize());

        // triangle_normal = transform_v(ModelView.inverse().transpose(), triangle_normal);
        // light_dir = transform_v(ModelView, light_dir);

        float intensity = std::max(0.f,  temp_norm.normalize() * temp_light.normalize());

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
    // clip(0, 0, width, height);
    // light_dir.normalize();


    // int *zbuffer = new int[width*height];
    // for (int i=width*height; i--; zbuffer[i] = 0);


    Shader shader;
    shader.cob = Viewport * Projection * ModelView;
    shader.cob_IT = shader.cob.inverse().transpose();

    shader.M = Projection * ModelView;
    shader.M_IT = shader.M.inverse().transpose();

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

