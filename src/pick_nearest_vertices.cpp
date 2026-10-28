#include <pick_nearest_vertices.h>
#include <iostream>
#include <igl/Hit.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/unproject.h>
#include <algorithm>


//NOTE: You can use the libigl igl::unproject function to implement this method

//Input:
//  win - window coordinate of mouse click (x_window, y_window, 0)
//  view - view transformation matrix
//  proj - projection matrix
//  viewport - viewport coordinates .
//  V - 3xn dense matrix of mesh vertices, each row of the matrix is a single vertex.
//  radius - selection radius for vertex picking
//Output:
//  verts - vertex ids (rows in V) of selected vertices

bool pick_nearest_vertices(std::vector<unsigned int> &verts, Eigen::Ref<const Eigen::Vector3d> win, 
                           Eigen::Ref<const Eigen::Matrix44f> view, Eigen::Ref<const Eigen::Matrix44f> proj, Eigen::Vector4f viewport,
                           Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double radius) {

    verts.clear();
   Eigen::Vector3d win_near = win;
    Eigen::Vector3d win_far = win;

    win_near(2) = 0.0;  // 近平面
    win_far(2) = 1.0;   // 远平面

    Eigen::Vector3d ray_origin;
    Eigen::Vector3d ray_far;

    igl::unproject(win_near, view, proj, viewport, ray_origin);
    igl::unproject(win_far,  view, proj, viewport, ray_far);

    Eigen::Vector3d ray_direction =
        (ray_far - ray_origin).normalized();

        std::vector<igl::Hit> hits;

    bool did_hit = igl::ray_mesh_intersect(
        ray_origin.transpose(),
        ray_direction.transpose(),
        V,
        F,
        hits
    );
    

    if (!did_hit || hits.empty()) {
        return false;
    }

    auto closest_it = std::min_element(
        hits.begin(), hits.end(),
        [](const igl::Hit &a, const igl::Hit &b) {
            return a.t < b.t;
        }
    );

    const igl::Hit &closest_hit = *closest_it;

    Eigen::RowVector3d hit_position =
    ray_origin.transpose()
    + closest_hit.t * ray_direction.transpose();

    for (int i = 0; i < V.rows(); i ++) {
        double distance = (V.row(i) - hit_position).norm();

        if (distance < radius) {
            verts.push_back(i);
        }
    }


    return !verts.empty();;
}