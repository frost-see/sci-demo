#pragma once

#include <functional>
#include <vector>

namespace sci {

class Derivative {
public:
    static double centralDifference(std::function<double(double)> f, double x, double h = 1e-5);
    static double forwardDifference(std::function<double(double)> f, double x, double h = 1e-5);
    static double backwardDifference(std::function<double(double)> f, double x, double h = 1e-5);
    static double secondDerivative(std::function<double(double)> f, double x, double h = 1e-5);
    static std::vector<double> derivatives(
        std::function<double(double)> f,
        const std::vector<double>& points,
        double h = 1e-5);
};

} // namespace sci
