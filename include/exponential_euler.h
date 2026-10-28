#pragma once

#include <Eigen/Dense>
#include <EigenTypes.h>
#include <vector>

void exponential_euler(
    Eigen::VectorXd &q,
    Eigen::VectorXd &qdot,
    double dt,
    std::vector<Eigen::Matrix66d> &masses,
    Eigen::Ref<const Eigen::VectorXd> forces
);