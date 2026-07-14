#include "curve_fitting.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

void validateFitInputs(const std::vector<double>& x, const std::vector<double>& y, std::size_t degree) {
    if (x.empty() || y.empty()) {
        throw std::invalid_argument("x and y must not be empty");
    }
    if (x.size() != y.size()) {
        throw std::invalid_argument("x and y must have the same size");
    }
    if (x.size() <= degree) {
        throw std::invalid_argument("number of points must exceed polynomial degree");
    }

    for (double value : x) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("x values must be finite");
        }
    }
    for (double value : y) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("y values must be finite");
        }
    }
}

double evaluatePolynomial(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
        result = result * x + *it;
    }
    return result;
}

std::vector<double> solveAugmentedSystem(std::vector<std::vector<double>> matrix) {
    const std::size_t size = matrix.size();

    for (std::size_t pivot = 0; pivot < size; ++pivot) {
        std::size_t bestRow = pivot;
        for (std::size_t row = pivot + 1; row < size; ++row) {
            if (std::abs(matrix[row][pivot]) > std::abs(matrix[bestRow][pivot])) {
                bestRow = row;
            }
        }

        if (std::abs(matrix[bestRow][pivot]) <= std::numeric_limits<double>::epsilon()) {
            throw std::domain_error("fit system is singular");
        }

        if (bestRow != pivot) {
            std::swap(matrix[pivot], matrix[bestRow]);
        }

        const double pivotValue = matrix[pivot][pivot];
        for (std::size_t column = pivot; column <= size; ++column) {
            matrix[pivot][column] /= pivotValue;
        }

        for (std::size_t row = 0; row < size; ++row) {
            if (row == pivot) {
                continue;
            }

            const double factor = matrix[row][pivot];
            for (std::size_t column = pivot; column <= size; ++column) {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
        }
    }

    std::vector<double> solution(size);
    for (std::size_t row = 0; row < size; ++row) {
        solution[row] = matrix[row][size];
    }
    return solution;
}

} // namespace

namespace sci {

double PolynomialFitResult::evaluate(double x) const {
    if (coefficients.empty()) {
        throw std::logic_error("coefficients must not be empty");
    }
    if (!std::isfinite(x)) {
        throw std::invalid_argument("x must be finite");
    }
    return evaluatePolynomial(coefficients, x);
}

PolynomialFitResult CurveFitting::polynomialFit(
    const std::vector<double>& x,
    const std::vector<double>& y,
    std::size_t degree) {
    validateFitInputs(x, y, degree);

    const std::size_t order = degree + 1;
    std::vector<std::vector<double>> augmented(order, std::vector<double>(order + 1, 0.0));
    std::vector<double> powerSums(2 * degree + 1, 0.0);
    std::vector<double> rhsSums(order, 0.0);

    for (std::size_t index = 0; index < x.size(); ++index) {
        const double xValue = x[index];
        const double yValue = y[index];
        double power = 1.0;
        for (std::size_t exponent = 0; exponent < powerSums.size(); ++exponent) {
            powerSums[exponent] += power;
            if (exponent < order) {
                rhsSums[exponent] += yValue * power;
            }
            power *= xValue;
        }
    }

    for (std::size_t row = 0; row < order; ++row) {
        for (std::size_t column = 0; column < order; ++column) {
            augmented[row][column] = powerSums[row + column];
        }
        augmented[row][order] = rhsSums[row];
    }

    PolynomialFitResult result;
    result.coefficients = solveAugmentedSystem(std::move(augmented));

    double yMean = 0.0;
    for (double value : y) {
        yMean += value;
    }
    yMean /= static_cast<double>(y.size());

    double residualSumSquares = 0.0;
    double totalSumSquares = 0.0;
    for (std::size_t index = 0; index < x.size(); ++index) {
        const double fitted = evaluatePolynomial(result.coefficients, x[index]);
        const double residual = y[index] - fitted;
        residualSumSquares += residual * residual;

        const double centered = y[index] - yMean;
        totalSumSquares += centered * centered;
    }

    if (totalSumSquares <= std::numeric_limits<double>::epsilon()) {
        result.rSquared = residualSumSquares <= std::numeric_limits<double>::epsilon() ? 1.0 : 0.0;
    } else {
        result.rSquared = 1.0 - residualSumSquares / totalSumSquares;
    }

    return result;
}

} // namespace sci
