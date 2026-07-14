#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "../src/curve_fitting.h"

namespace {

using sci::CurveFitting;

constexpr double kTolerance = 1e-9;

TEST(CurveFittingTest, FitsLinearDataExactly) {
    const std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    const std::vector<double> y = {1.0, 3.0, 5.0, 7.0};

    const auto fit = CurveFitting::polynomialFit(x, y, 1);

    ASSERT_EQ(fit.coefficients.size(), 2U);
    EXPECT_NEAR(fit.coefficients[0], 1.0, kTolerance);
    EXPECT_NEAR(fit.coefficients[1], 2.0, kTolerance);
    EXPECT_NEAR(fit.evaluate(4.0), 9.0, kTolerance);
    EXPECT_NEAR(fit.rSquared, 1.0, kTolerance);
}

TEST(CurveFittingTest, FitsQuadraticDataExactly) {
    const std::vector<double> x = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<double> y;
    y.reserve(x.size());
    for (double value : x) {
        y.push_back(1.0 + 2.0 * value + 3.0 * value * value);
    }

    const auto fit = CurveFitting::polynomialFit(x, y, 2);

    ASSERT_EQ(fit.coefficients.size(), 3U);
    EXPECT_NEAR(fit.coefficients[0], 1.0, kTolerance);
    EXPECT_NEAR(fit.coefficients[1], 2.0, kTolerance);
    EXPECT_NEAR(fit.coefficients[2], 3.0, kTolerance);
    EXPECT_NEAR(fit.evaluate(3.0), 34.0, kTolerance);
    EXPECT_NEAR(fit.rSquared, 1.0, kTolerance);
}

TEST(CurveFittingTest, RejectsInvalidInputsAndSingularData) {
    EXPECT_THROW(CurveFitting::polynomialFit({}, {}, 1), std::invalid_argument);
    EXPECT_THROW(CurveFitting::polynomialFit({0.0, 1.0}, {1.0}, 1), std::invalid_argument);
    EXPECT_THROW(
        CurveFitting::polynomialFit({0.0, std::numeric_limits<double>::infinity()}, {1.0, 2.0}, 1),
        std::invalid_argument);
    EXPECT_THROW(
        CurveFitting::polynomialFit({0.0, 1.0}, {1.0, std::numeric_limits<double>::quiet_NaN()}, 1),
        std::invalid_argument);
    EXPECT_THROW(CurveFitting::polynomialFit({0.0, 1.0}, {1.0, 2.0}, 2), std::invalid_argument);
    EXPECT_THROW(
        CurveFitting::polynomialFit({1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}, 1),
        std::domain_error);
}

TEST(CurveFittingTest, RejectsInvalidEvaluationPoints) {
    const sci::PolynomialFitResult emptyResult;
    EXPECT_THROW(emptyResult.evaluate(0.0), std::logic_error);

    const auto fit = CurveFitting::polynomialFit({0.0, 1.0, 2.0}, {0.0, 1.0, 4.0}, 2);
    EXPECT_THROW(fit.evaluate(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
}

} // namespace
