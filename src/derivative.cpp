#include "derivative.h"

#include <cmath>
#include <stdexcept>

namespace {

void validateInputs(const std::function<double(double)>& f, double x, double h) {
    if (!f) {
        throw std::invalid_argument("function must be callable");
    }
    if (!std::isfinite(x)) {
        throw std::invalid_argument("x must be finite");
    }
    if (!std::isfinite(h) || h <= 0.0) {
        throw std::invalid_argument("h must be finite and positive");
    }
}

double evaluateFinite(const std::function<double(double)>& f, double point) {
    if (!std::isfinite(point)) {
        throw std::domain_error("sample point is out of range");
    }

    const double value = f(point);
    if (!std::isfinite(value)) {
        throw std::domain_error("function evaluation produced non-finite value");
    }

    return value;
}

} // namespace

namespace sci {

double Derivative::centralDifference(std::function<double(double)> f, double x, double h) {
    validateInputs(f, x, h);
    return (evaluateFinite(f, x + h) - evaluateFinite(f, x - h)) / (2.0 * h);
}

double Derivative::forwardDifference(std::function<double(double)> f, double x, double h) {
    validateInputs(f, x, h);
    return (evaluateFinite(f, x + h) - evaluateFinite(f, x)) / h;
}

double Derivative::backwardDifference(std::function<double(double)> f, double x, double h) {
    validateInputs(f, x, h);
    return (evaluateFinite(f, x) - evaluateFinite(f, x - h)) / h;
}

double Derivative::secondDerivative(std::function<double(double)> f, double x, double h) {
    validateInputs(f, x, h);
    return (evaluateFinite(f, x + h) - 2.0 * evaluateFinite(f, x) + evaluateFinite(f, x - h)) / (h * h);
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
