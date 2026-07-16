#include "curve_fit.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

// Solve the linear system A*x = b via Gaussian elimination with partial
// pivoting.  Both A and b are modified in-place; the solution is returned.
// Throws std::domain_error when the matrix is (near-)singular.
std::vector<double> solveLinearSystem(
    std::vector<std::vector<double>>& A,
    std::vector<double>& b) {

    const int n = static_cast<int>(b.size());

    for (int col = 0; col < n; ++col) {
        // Partial pivoting: find the row with the largest absolute value in
        // the current column.
        int pivotRow = col;
        double maxVal = std::abs(A[col][col]);
        for (int row = col + 1; row < n; ++row) {
            if (std::abs(A[row][col]) > maxVal) {
                maxVal = std::abs(A[row][col]);
                pivotRow = row;
            }
        }

        if (maxVal < 1e-14) {
            throw std::domain_error(
                "singular or near-singular matrix; data may be degenerate");
        }

        if (pivotRow != col) {
            std::swap(A[col], A[pivotRow]);
            std::swap(b[col], b[pivotRow]);
        }

        // Forward elimination.
        for (int row = col + 1; row < n; ++row) {
            const double factor = A[row][col] / A[col][col];
            for (int j = col; j < n; ++j) {
                A[row][j] -= factor * A[col][j];
            }
            b[row] -= factor * b[col];
        }
    }

    // Back substitution.
    std::vector<double> solution(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        solution[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            solution[i] -= A[i][j] * solution[j];
        }
        solution[i] /= A[i][i];
    }

    return solution;
}

} // namespace

namespace sci {

std::vector<double> CurveFit::fit(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int degree) {

    if (degree < 2 || degree > 5) {
        throw std::invalid_argument("degree must be between 2 and 5");
    }
    if (x.empty() || y.empty()) {
        throw std::invalid_argument("x and y must not be empty");
    }
    if (x.size() != y.size()) {
        throw std::invalid_argument("x and y must have the same size");
    }
    const int n = static_cast<int>(x.size());
    if (n <= degree) {
        throw std::invalid_argument(
            "number of data points must exceed polynomial degree");
    }

    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) {
            throw std::invalid_argument(
                "all x and y values must be finite");
        }
    }

    // Number of coefficients.
    const int m = degree + 1;

    // Build the normal equations (A^T A) * c = (A^T y) where A is the
    // Vandermonde matrix A[i][j] = x[i]^j.
    // (A^T A)[p][q] = sum_i  x[i]^(p+q)
    // (A^T y)[p]    = sum_i  x[i]^p * y[i]
    std::vector<std::vector<double>> AtA(m, std::vector<double>(m, 0.0));
    std::vector<double> Aty(m, 0.0);

    for (int i = 0; i < n; ++i) {
        // Precompute powers x[i]^0 .. x[i]^(2*degree).
        std::vector<double> xpow(2 * m - 1, 1.0);
        for (int k = 1; k < 2 * m - 1; ++k) {
            xpow[k] = xpow[k - 1] * x[i];
        }

        for (int p = 0; p < m; ++p) {
            for (int q = 0; q < m; ++q) {
                AtA[p][q] += xpow[p + q];
            }
            Aty[p] += xpow[p] * y[i];
        }
    }

    return solveLinearSystem(AtA, Aty);
}

double CurveFit::evaluate(const std::vector<double>& coefficients, double x) {
    if (coefficients.empty()) {
        throw std::invalid_argument("coefficients must not be empty");
    }
    if (!std::isfinite(x)) {
        throw std::invalid_argument("x must be finite");
    }

    // Horner's method for numerical stability.
    double result = 0.0;
    for (int i = static_cast<int>(coefficients.size()) - 1; i >= 0; --i) {
        result = result * x + coefficients[i];
    }
    return result;
}

} // namespace sci
