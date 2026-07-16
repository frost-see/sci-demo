#pragma once

#include <vector>

namespace sci {

// CurveFit provides least-squares polynomial fitting for degree 2–5.
//
// fit() returns coefficients [a0, a1, ..., an] for the polynomial:
//   p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n
//
// evaluate() computes p(x) for a given coefficient vector.
class CurveFit {
public:
    // Fit a polynomial of the given degree to the data points (x, y).
    //
    // Constraints:
    //   - degree must be in [2, 5]
    //   - x and y must have equal, non-zero size
    //   - number of data points must exceed degree
    //   - all x and y values must be finite
    //
    // Throws std::invalid_argument for invalid inputs.
    // Throws std::domain_error if the data produces a singular system.
    static std::vector<double> fit(
        const std::vector<double>& x,
        const std::vector<double>& y,
        int degree);

    // Evaluate the polynomial defined by coefficients at point x.
    //
    // coefficients[i] is the coefficient of x^i.
    //
    // Throws std::invalid_argument if coefficients is empty or x is non-finite.
    static double evaluate(const std::vector<double>& coefficients, double x);
};

} // namespace sci
