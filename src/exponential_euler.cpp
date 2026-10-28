#include <exponential_euler.h>
#include <rodrigues.h>
#include <cassert>

//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities 
void exponential_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces) {
    
    int num_rigid_bodies = masses.size();

    for (int irb = 0; irb < num_rigid_bodies; irb ++) {
        Eigen::Matrix3d I_body = masses[irb].block<3, 3>(0, 0);
        Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(
            q.segment<9>(12 * irb).data()
        );

        Eigen::Vector3d angular_v = qdot.segment<3>(6 * irb);
        Eigen::Vector3d linear_v = qdot.segment<3>(6 * irb + 3);

        Eigen::Vector3d com = q.segment<3>(12 * irb + 9);
   
        Eigen::Vector3d force = forces.segment<3>(6 * irb + 3);

        Eigen::Vector3d torque = forces.segment<3>(6 * irb);

        double mass = masses[irb](3, 3);

        Eigen::Matrix3d I_world = R * I_body * R.transpose();
        Eigen::Vector3d linear_a = force / mass;
        Eigen::Vector3d angular_a = I_world.ldlt().solve(
            angular_v.cross(I_world * angular_v) + torque
        );

        angular_v += angular_a * dt;
        linear_v += linear_a * dt;

        Eigen::Matrix3d delta_R;
        rodrigues(delta_R, dt * angular_v);

        com += linear_v * dt;
        R = delta_R * R;
        Eigen::Map<Eigen::Matrix3d> R_in_q(
            q.segment<9>(12 * irb).data()
        );
        R_in_q = R;

        q.segment<3>(12 * irb + 9) = com;
        qdot.segment<3>(6 * irb) = angular_v;
        qdot.segment<3>(6 * irb + 3) = linear_v;
    }
}