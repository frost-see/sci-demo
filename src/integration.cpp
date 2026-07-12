#include "integration.h"

#include <cmath>
#include <stdexcept>

namespace {

void validateInputs(const std::function<double(double)>& f, double a, double b, int intervals) {
    if (!f) {
        throw std::invalid_argument("function must be callable");
    }
    if (!std::isfinite(a) || !std::isfinite(b)) {
        throw std::invalid_argument("integration bounds must be finite");
    }
    if (intervals <= 0) {
        throw std::invalid_argument("intervals must be positive");
    }
}

double evaluateFinite(const std::function<double(double)>& f, double point) {
    const double value = f(point);
    if (!std::isfinite(value)) {
        throw std::domain_error("function evaluation produced non-finite value");
    }
    return value;
}

} // namespace

namespace sci {

double Integration::trapezoidal(const std::function<double(double)>& f, double a, double b, int intervals) {
    validateInputs(f, a, b, intervals);
    if (a == b) {
        return 0.0;
    }

    const double start = a < b ? a : b;
    const double end = a < b ? b : a;
    const double sign = a < b ? 1.0 : -1.0;
    const double h = (end - start) / static_cast<double>(intervals);

    double sum = 0.5 * (evaluateFinite(f, start) + evaluateFinite(f, end));
    for (int i = 1; i < intervals; ++i) {
        sum += evaluateFinite(f, start + static_cast<double>(i) * h);
    }

    return sign * h * sum;
}

double Integration::simpson(const std::function<double(double)>& f, double a, double b, int intervals) {
    validateInputs(f, a, b, intervals);
    if (intervals % 2 != 0) {
        throw std::invalid_argument("simpson rule requires an even number of intervals");
    }
    if (a == b) {
        return 0.0;
    }

    const double start = a < b ? a : b;
    const double end = a < b ? b : a;
    const double sign = a < b ? 1.0 : -1.0;
    const double h = (end - start) / static_cast<double>(intervals);

    double sum = evaluateFinite(f, start) + evaluateFinite(f, end);
    for (int i = 1; i < intervals; ++i) {
        const double coefficient = (i % 2 == 0) ? 2.0 : 4.0;
        sum += coefficient * evaluateFinite(f, start + static_cast<double>(i) * h);
    }

    return sign * h * sum / 3.0;
}

} // namespace sci
