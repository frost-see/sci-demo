#include "derivative.h"

namespace sci {

double Derivative::centralDifference(std::function<double(double)> f, double x, double h) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

double Derivative::forwardDifference(std::function<double(double)> f, double x, double h) {
    return (f(x + h) - f(x)) / h;
}

double Derivative::backwardDifference(std::function<double(double)> f, double x, double h) {
    return (f(x) - f(x - h)) / h;
}

double Derivative::secondDerivative(std::function<double(double)> f, double x, double h) {
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
}

std::vector<double> Derivative::derivatives(
    std::function<double(double)> f,
    const std::vector<double>& points,
    double h) {
    
    std::vector<double> result;
    result.reserve(points.size());
    
    for (double point : points) {
        result.push_back(centralDifference(f, point, h));
    }
    
    return result;
}

} // namespace sci
