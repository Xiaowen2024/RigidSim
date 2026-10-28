#include <rigid_body_jacobian.h>

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> x) {

    Eigen::Matrix3d x_cross;
    x_cross <<
        0.0, -x(2),  x(1),
        x(2),  0.0, -x(0),
        -x(1),  x(0),  0.0;

    J.block<3, 3>(0, 0) = R * x_cross.transpose() * R.transpose();
    J.block<3, 3>(0, 3).setIdentity();             
}

