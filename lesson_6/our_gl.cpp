#include "our_gl.h"
#include "matrix.h"

Matrix Viewport;
Matrix ModelView;
Matrix Projection;

IShader::~IShader() {}



Vec3f normalize(Vec3f a){
    double mag = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    return Vec3f(a.x / mag, a.y / mag, a.z / mag);
}

Vec3f cross(Vec3f a, Vec3f b){
    float x = a.y * b.z - a.z * b.y;
    float y = a.z * b.x - a.x * b.z;
    float z = a.x * b.y - a.y * b.x;
    return Vec3f(x,y,z);
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
}

void view(Vec3f center, Vec3f eye, Vec3f up){
    ModelView = Matrix::createIdentity(4);
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
    ModelView = Minv*Tr;
    // return ModelView;
}

void clip(int x, int y, int width, int height){
    Viewport = Matrix::createIdentity(4);
    Viewport(0,3) = x + width/2.0;
    Viewport(1,3) = y+ height/2.0;
    Viewport(2,3) = 255/2.0;

    Viewport(0,0) = width / 2.0;
    Viewport(1,1) = height/2.0;
    Viewport(2,2) = 255/2.0;
    // return Viewport;
}

void proj(float coeff){
    Projection = Matrix::createIdentity(4);
    Projection(3,2) = coeff;
}

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

void triangle(int *zbuffer, Matrix *pts_m, TGAImage &image, IShader &shader, int width) {

    Vec3f pts[3];
    pts[0] = Vec3f(pts_m[0](0,0), pts_m[0](0,1), pts_m[0](0,2));
    pts[1] = Vec3f(pts_m[1](0,0), pts_m[1](0,1), pts_m[1](0,2));


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

    if(y_max < 0){ y_max = 0;};
    if(y_min < 0){ y_min = 0;};
    if(x_max < 0){ x_max = 0;};
    if(x_min < 0){ x_min = 0;};

    if(y_max > image.get_height() - 1){ y_max = image.get_height() - 1;};
    if(y_min > image.get_height() - 1){ y_min = image.get_height() - 1;};
    if(x_max > image.get_width() - 1){ x_max = image.get_width() - 1;};
    if(x_min > image.get_width() - 1){ x_min = image.get_width() - 1;};

    // for(int i = 0; i < 3 ; i++){
    //     texture_coords[i].x *= texture_img.get_width();
    //     texture_coords[i].y *= texture_img.get_height();
    // }

    for(int x = x_min; x <= x_max; x++){
        for(int y = y_min; y <= y_max; y++){
            float bc[3];
            bary(x, y, pts, bc);

            if((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)){
                float z = bc[0] * pts[0].z + bc[1] * pts[1].z + bc[2] * pts[2].z;


                if(zbuffer[x + y * width] < z){
            //         zbuffer[x + y * width] = z;

            //         float normal_x = bc[0] * normal_coords[0].x + bc[1] * normal_coords[1].x + bc[2] * normal_coords[2].x;
            //         float normal_y = bc[0] * normal_coords[0].y + bc[1] * normal_coords[1].y + bc[2] * normal_coords[2].y;
            //         float normal_z = bc[0] * normal_coords[0].z + bc[1] * normal_coords[1].z + bc[2] * normal_coords[2].z;

            //         Vec3f normal = Vec3f(normal_x, normal_y, normal_z);
            //         normal.normalize();
            //         float intensity = normal * light_dir;
            //         intensity = -intensity;
            //         if(intensity < 0){
            //             intensity = 0;
            //         }

            //         float texture_x = bc[0] * texture_coords[0].x + bc[1] * texture_coords[1].x + bc[2] * texture_coords[2].x;
            //         float texture_y = bc[0] * texture_coords[0].y + bc[1] * texture_coords[1].y + bc[2] * texture_coords[2].y;
            //         TGAColor text_color = texture_img.get(roundf(texture_x), roundf(texture_y));
            //         TGAColor adjusted_color = TGAColor((float)text_color.r * intensity, (float)text_color.g * intensity, (float)text_color.g * intensity, 255);
            //         image.set(x,y,adjusted_color);    
                }
                
            }
        }
    }

}