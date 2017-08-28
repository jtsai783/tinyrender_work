#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include <iostream>
#include "model.h"
#include "matrix.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
TGAImage *texture_img = NULL;
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

Matrix view(Vec3f center, Vec3f eye, Vec3f up){
    // Matrix m = Matrix::createIdentity(4);
    Vec3f z = (eye-center);
    z = normalize(z);
    Vec3f x = cross(up,z);
    x = normalize(x);
    Vec3f y = cross(z,x);
    y = normalize(y);
    Matrix Minv = Matrix::createIdentity(4);
    Matrix Tr   = Matrix::createIdentity(4);
    for (int i=0; i<3; i++) {
        Minv(0,i) = x[i];
        Minv(1,i) = y[i];
        Minv(2,i) = z[i];
        Tr(i,3) = -center[i];
    }
    Matrix ModelView = Minv*Tr;

    return ModelView;
}

Matrix clip(int x, int y, int width, int height){
    Matrix m = Matrix::createIdentity(4);
    m(0,3) = x + width/2.0;
    m(1,3) = y+ height/2.0;
    m(2,3) = 255/2.0;

    m(0,0) = width / 2.0;
    m(1,1) = height/2.0;
    m(2,2) = 255/2.0;
    return m;
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


void triangle(Vec3f *texture_coords, float *zbuffer, Vec3f *pts, TGAImage &image, TGAImage &texture_img, Vec3f light_dir, Vec3f *normal_coords) {

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

    for(int i = 0; i < 3 ; i++){
        texture_coords[i].x *= texture_img.get_width();
        texture_coords[i].y *= texture_img.get_height();
    }

    // std::cout << texture_coords[0] << " " << texture_coords[1] << " " << texture_coords[2];
    //for each point in the bounding box, test to see if its in the triangle
    for(int x = x_min; x <= x_max; x++){
        for(int y = y_min; y <= y_max; y++){
            float bc[3];
            bary(x, y, pts, bc);
            //calculate z-value of the point


            if((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)){
                float z = bc[0] * pts[0].z + bc[1] * pts[1].z + bc[2] * pts[2].z;
                // std::cout << z << std::endl;


                if(zbuffer[x + y * width] < z){
                    zbuffer[x + y * width] = z;

                    float normal_x = bc[0] * normal_coords[0].x + bc[1] * normal_coords[1].x + bc[2] * normal_coords[2].x;
                    float normal_y = bc[0] * normal_coords[0].y + bc[1] * normal_coords[1].y + bc[2] * normal_coords[2].y;
                    float normal_z = bc[0] * normal_coords[0].z + bc[1] * normal_coords[1].z + bc[2] * normal_coords[2].z;

                    Vec3f normal = Vec3f(normal_x, normal_y, normal_z);
                    normal.normalize();
                    float intensity = normal * light_dir;
                    // std::cout << intensity << std::endl;
                    intensity = -intensity;
                    if(intensity < 0){
                        intensity = 0;
                    }

                    float texture_x = bc[0] * texture_coords[0].x + bc[1] * texture_coords[1].x + bc[2] * texture_coords[2].x;
                    float texture_y = bc[0] * texture_coords[0].y + bc[1] * texture_coords[1].y + bc[2] * texture_coords[2].y;
                    TGAColor text_color = texture_img.get(roundf(texture_x), roundf(texture_y));
                    TGAColor adjusted_color = TGAColor((float)text_color.r * intensity, (float)text_color.g * intensity, (float)text_color.g * intensity, 255);
                    image.set(x,y,adjusted_color);    
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

    // std::cout << cob;

    // Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)};
    // Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)};
    // Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)};

    // triangle(t0, image, red);
    // triangle(t1, image, white);
    // triangle(t2, image, green);

    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    Vec3f light_dir(0,2,-1);



    for (int i= 0; i<model->nfaces(); i++) { 
        std::vector<int> face = model->face(i);
        std::vector<int> texture_face = model->texture_face(i);
        Vec3f screen_coords[3]; 
        // Vec3f world_coords[3];
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
            // screen = clip_m * screen;
            // std::cout << texture_face[j];
            Vec3f vt = model->texture_vert(texture_face[j]);
            Vec3f vn = model->normal_vert(face[j]);
            // float coeff = 1.0 - (v.z / 2.0);
            // screen_coords[j] = Vec3f((v.x/coeff+1.)*800/2. + 100, (v.y/coeff+1.)*800/2. + 100, v.z/coeff);
            screen_coords[j] = Vec3f(screen(0,0), screen(1,0), screen(2,0));
            // std::cout << v;
            // std::cout << cob;
            // std::cout << screen_coords[j];
            Matrix vn_m = Matrix(4,1);
            vn_m(0, 0) = vn.x;
            vn_m(1, 0) = vn.y;
            vn_m(2, 0) = vn.z;
            vn_m(3, 0) = 0;
            vn_m = norm_cob * vn_m;
            // normal_coords[j] = Vec3f(vn.x/coeff, vn.y/coeff, vn.z/coeff);
            normal_coords[j] = Vec3f(vn_m(0,0),vn_m(1,0),vn_m(2,0));
            // world_coords[j]  = v;
            texture_coords[j] = vt;
            // std::cout << vt;
        }

        // Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]); 
        // n.normalize(); 
        // float intensity = n*light_dir; 
        // if (intensity>0) { 
            // std::cout << screen_coords[2] << std::endl;
            triangle(texture_coords, zbuffer, screen_coords, image, texture_img, light_dir, normal_coords);
        // } 
    }


    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    // texture_img.write_tga_file("text_text.tga");
    delete model;
    return 0;
}

