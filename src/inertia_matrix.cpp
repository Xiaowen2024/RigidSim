#include <inertia_matrix.h>
#include <cassert>

//Input:
//  V - the nx3 matrix of vertices.
//  F - the mx3 matrix of triangle vertex indices.
//  density - the material density.
//Output:
//  I - the 3x3 angular inertia matrix
//  center - the center of mass of the object
//  mass - the total mass of the object

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {

    auto volume_monomial = [&](int px, int py, int pz) -> double {
        double result = 0.0;

        for (int e = 0; e < F.rows(); ++e) {
            Eigen::Vector3d a = V.row(F(e, 0)).transpose();
            Eigen::Vector3d b = V.row(F(e, 1)).transpose();
            Eigen::Vector3d c = V.row(F(e, 2)).transpose();

            Eigen::Vector3d cross = (b - a).cross(c - a);
            double area = 0.5 * cross.norm();

            // 必须是非退化三角形
            if (area == 0.0) {
                continue;
            }

            Eigen::Vector3d normal = cross.normalized();

            // 三角形上的 4 个积分点
            Eigen::Vector3d p0 = (a + b + c) / 3.0;
            Eigen::Vector3d p1 = 0.6 * a + 0.2 * b + 0.2 * c;
            Eigen::Vector3d p2 = 0.2 * a + 0.6 * b + 0.2 * c;
            Eigen::Vector3d p3 = 0.2 * a + 0.2 * b + 0.6 * c;

            // g = [x^(px+1)y^py z^pz / (px+1), 0, 0]
            // div(g) = x^px y^py z^pz
            auto g_x = [&](const Eigen::Vector3d &p) {
                return std::pow(p.x(), px + 1)
                    * std::pow(p.y(), py)
                    * std::pow(p.z(), pz)
                    / static_cast<double>(px + 1);
            };

            // 精确积分所有三次及以下多项式
            double face_integral = area * (
                -27.0 / 48.0 * g_x(p0)
                +25.0 / 48.0 * (g_x(p1) + g_x(p2) + g_x(p3))
            );

            // ∫surface g·n dA；g 只有 x 分量
            result += normal.x() * face_integral;
        }

        return result;
    };

    double volume = volume_monomial(0, 0, 0);

    Eigen::Vector3d first_moment;
    first_moment <<
        volume_monomial(1, 0, 0),  // ∫x dV
        volume_monomial(0, 1, 0),  // ∫y dV
        volume_monomial(0, 0, 1);  // ∫z dV

    center = first_moment / volume;
    mass = density * volume;

    Eigen::Matrix3d S;

    S(0, 0) = volume_monomial(2, 0, 0);  // ∫x² dV
    S(1, 1) = volume_monomial(0, 2, 0);  // ∫y² dV
    S(2, 2) = volume_monomial(0, 0, 2);  // ∫z² dV

    S(0, 1) = S(1, 0) = volume_monomial(1, 1, 0); // ∫xy dV
    S(0, 2) = S(2, 0) = volume_monomial(1, 0, 1); // ∫xz dV
    S(1, 2) = S(2, 1) = volume_monomial(0, 1, 1); // ∫yz dV

    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();

    Eigen::Matrix3d I_origin =
        density * (S.trace() * identity - S);

    I = I_origin
        - mass * (
            center.squaredNorm() * identity
            - center * center.transpose()
        );

}