#pragma once

#include <cstddef>
#include <vector>

namespace sci {

struct PolynomialFitResult {
    std::vector<double> coefficients;
    double rSquared = 0.0;

    double evaluate(double x) const;
};

class CurveFitting {
public:
    static PolynomialFitResult polynomialFit(
        const std::vector<double>& x,
        const std::vector<double>& y,
        std::size_t degree);
};

} // namespace sci
