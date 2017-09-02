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

const int width  = 2000;
const int height = 2000;

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
    Matrix cob;
    Matrix cob_IT;

    Matrix M;
    Matrix M_IT;

    Vec3f triangle_normal;
    Vec3f vertex_normal[3];
    Vec3f verts[3];

    Matrix A = Matrix(3,3);

    Matrix iA = Matrix(3,1);
    Matrix jA = Matrix(3,1);

    virtual Vec3f vertex(int nface, int nthvert, TGAImage &image){
        std::vector<int> face = model->face(nface);
        Vec3f v = model->vert(face[nthvert * 3]);
        Vec3f thisVert = transform_p(cob, v);
        texture_coords[nthvert] = model->texture_vert(face[nthvert * 3 + 1]);
        Vec3f normal = model->normal_vert(face[nthvert * 3 + 2]);
        vertex_normal[nthvert] = transform_v(M_IT, normal);
        verts[nthvert] = v;
        if(nthvert == 2){
           //calculate the A matrix
            Vec3f p0p1 = transform_p(M,verts[1]) - transform_p(M,verts[0]);
            Vec3f p0p2 = transform_p(M,verts[2]) - transform_p(M,verts[0]);
            // Vec3f p0p1 = verts[1] - verts[0];
            // Vec3f p0p2 = verts[2] - verts[0];
            A(0,0) = p0p1.x;
            A(0,1) = p0p1.y;
            A(0,2) = p0p1.z;

            A(1,0) = p0p2.x;
            A(1,1) = p0p2.y;
            A(1,2) = p0p2.z;

            A(2,0) = vertex_normal[2].x;
            A(2,1) = vertex_normal[2].y;
            A(2,2) = vertex_normal[2].z;

            iA(0,0) = texture_coords[1].x - texture_coords[0].x;
            iA(1,0) = texture_coords[2].x - texture_coords[0].x;
            iA(2,0) = 0;

            jA(0,0) = texture_coords[1].y - texture_coords[0].y;
            jA(1,0) = texture_coords[2].y - texture_coords[0].y;
            jA(2,0) = 0;

            Matrix i = A.inverse() * iA;
            Matrix j = A.inverse() * jA;

            Vec3f i_vec(i(0,0), i(0,1), i(0,2));
            i_vec = i_vec.normalize() * 30;

            Vec2i p0 = Vec2i((int)thisVert.x, (int)thisVert.y);
            Vec2i p1 = p0 + Vec2i(i_vec.x, i_vec.y);

            line(p0, p1, image, green);


        }
        return thisVert;
    }

    virtual bool fragment(float *bc, TGAColor &color){

        //get color
        float texture_x = bc[0] * texture_coords[0].x + bc[1] * texture_coords[1].x + bc[2] * texture_coords[2].x;
        float texture_y = bc[0] * texture_coords[0].y + bc[1] * texture_coords[1].y + bc[2] * texture_coords[2].y;

        texture_x = texture_x * (model->diffuse).get_width();
        texture_y = texture_y * (model->diffuse).get_height();
        color = (model->diffuse).get(roundf(texture_x), roundf(texture_y));



        //setup A
        Vec3f n = vertex_normal[0] * bc[0] + vertex_normal[1] * bc[1] + vertex_normal[2] * bc[2];
        // n = n.normalize();

        // A(2,0) = n.x;
        // A(2,1) = n.y;
        // A(2,2) = n.z;

        // Matrix i = A.inverse() * iA; 
        // Matrix j = A.inverse() * jA;

        // //read tangent space normal from map
        // TGAColor normal_color = (model->normal_map).get(roundf(texture_x), roundf(texture_y));
        // Vec3f normal_vec = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        // normal_vec = normal_vec * (1.0/255.0) * 2.0 - Vec3f(1.0,1.0,1.0);

        // Matrix Tan = Matrix(3,3);
        // Tan(0,0) = i(0,0);
        // Tan(0,1) = i(0,1);
        // Tan(0,2) = i(0,2);

        // Tan(1,0) = j(0,0);
        // Tan(1,1) = j(0,1);
        // Tan(1,2) = j(0,2);

        // Tan(2,0) = n.x;
        // Tan(2,1) = n.y;
        // Tan(2,2) = n.z;

        // n = transform_v(Tan.transpose(), normal_vec);

        // n = (n + Vec3f(1.0,1.0,1.0)) * 0.5 * 255.0;

        // std::cout << n << std::endl;

        // color = TGAColor(n.x, n.y, n.z, 255);

        Vec3f temp_light = transform_v(M, light_dir);
        // Vec3f temp_norm = transform_v(M_IT, normal_vec);

        // float intensity = std::max(0.f,  temp_norm.normalize() * temp_light.normalize());


        float intensity = std::max(0.f,  n.normalize() * temp_light.normalize());

        // float intensity = 1;

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


    TGAImage image(width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    // image = TGAImage(width, height, TGAImage::RGB);

    view(center, camera, up);
    proj(-1.0/(camera - center).norm());
    clip(width/8, height/8, width*3/4, height*3/4);


    Shader shader;
    shader.cob = Viewport * Projection * ModelView;
    shader.cob_IT = shader.cob.inverse().transpose();

    shader.M = Projection * ModelView;
    shader.M_IT = shader.M.inverse().transpose();

    for (int i= 0; i < model->nfaces(); i++) { 
        std::vector<int> face = model->face(i);
        Vec3f screen_coords[3]; 
        for (int j=0; j<3; j++) { 
            screen_coords[j] = shader.vertex(i, j, image);
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

