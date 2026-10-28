#include <rodrigues.h>
#include <cmath>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {

    double theta = omega.norm();
    if (theta < 1e-12) {
        R.setIdentity();
        return;
    }

    Eigen::Vector3d k = omega / theta;
    Eigen::Matrix3d K;
    K << 0.0,   -k.z(),  k.y(),
         k.z(),  0.0,   -k.x(),
        -k.y(),  k.x(),  0.0;

    R = Eigen::Matrix3d::Identity()
        + std::sin(theta) * K
        + (1.0 - std::cos(theta)) * K * K;
}