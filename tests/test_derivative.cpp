#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "../src/derivative.h"

namespace {

using sci::Derivative;

constexpr double kTightTolerance = 1e-6;
constexpr double kLooseTolerance = 1e-3;

double linear(double x) {
    return 2.0 * x + 3.0;
}

double quadratic(double x) {
    return x * x;
}

double cubic(double x) {
    return x * x * x;
}

double sine(double x) {
    return std::sin(x);
}

double exponential(double x) {
    return std::exp(x);
}

TEST(DerivativeMathTest, CoversBasicFunctionsAndArithmeticRules) {
    EXPECT_NEAR(Derivative::centralDifference(linear, 1.0), 2.0, kTightTolerance);
    EXPECT_NEAR(Derivative::centralDifference(quadratic, 2.0), 4.0, kTightTolerance);
    EXPECT_NEAR(Derivative::centralDifference(sine, 0.0), 1.0, kTightTolerance);
    EXPECT_NEAR(Derivative::centralDifference(exponential, 0.0), 1.0, kTightTolerance);

    const auto arithmetic = [](double x) {
        return (x * x + 3.0 * x) / (x + 2.0);
    };
    const double x = 1.5;
    const double expected = (x * x + 4.0 * x + 6.0) / ((x + 2.0) * (x + 2.0));
    EXPECT_NEAR(Derivative::forwardDifference(arithmetic, x), expected, kLooseTolerance);
    EXPECT_NEAR(Derivative::backwardDifference(arithmetic, x), expected, kLooseTolerance);
}

TEST(DerivativeMathTest, CoversChainHigherOrderAndPartialDerivatives) {
    const auto chained = [](double x) {
        return std::sin(x * x);
    };
    const double chainPoint = 0.7;
    EXPECT_NEAR(
        Derivative::centralDifference(chained, chainPoint),
        2.0 * chainPoint * std::cos(chainPoint * chainPoint),
        kLooseTolerance);

    EXPECT_NEAR(Derivative::secondDerivative(quadratic, 1.0), 2.0, kLooseTolerance);
    EXPECT_NEAR(Derivative::secondDerivative(cubic, 1.0), 6.0, 1e-2);

    const auto surface = [](double x, double y) {
        return x * x * y + std::sin(y);
    };
    const double px = 2.0;
    const double py = 0.4;
    EXPECT_NEAR(Derivative::centralDifference([&](double x) { return surface(x, py); }, px), 2.0 * px * py, kLooseTolerance);
    EXPECT_NEAR(Derivative::centralDifference([&](double y) { return surface(px, y); }, py), px * px + std::cos(py), kLooseTolerance);
}

TEST(DerivativeMathTest, ComputesBatchDerivativesAcrossMultiplePoints) {
    const std::vector<double> points = {0.0, 1.0, 2.0};
    const std::vector<double> derivs = Derivative::derivatives(quadratic, points);

    ASSERT_EQ(derivs.size(), points.size());
    EXPECT_NEAR(derivs[0], 0.0, kLooseTolerance);
    EXPECT_NEAR(derivs[1], 2.0, kLooseTolerance);
    EXPECT_NEAR(derivs[2], 4.0, kLooseTolerance);
}

TEST(DerivativeSpecialPointTest, HandlesStationaryAndNonDifferentiablePoints) {
    EXPECT_NEAR(Derivative::centralDifference(quadratic, 0.0), 0.0, kTightTolerance);

    const auto absoluteValue = [](double x) {
        return std::abs(x);
    };
    EXPECT_NEAR(Derivative::centralDifference(absoluteValue, 0.0), 0.0, kTightTolerance);
    EXPECT_NEAR(Derivative::forwardDifference(absoluteValue, 0.0), 1.0, kLooseTolerance);
    EXPECT_NEAR(Derivative::backwardDifference(absoluteValue, 0.0), -1.0, kLooseTolerance);

    const auto relu = [](double x) {
        return x < 0.0 ? 0.0 : x;
    };
    EXPECT_NEAR(Derivative::centralDifference(relu, 0.0), 0.5, kTightTolerance);
}

TEST(DerivativeNumericTest, RespectsToleranceAndStepSizeChoices) {
    const double x = 0.3;
    const double expected = std::cos(x);
    const double coarse = Derivative::centralDifference(sine, x, 1e-2);
    const double fine = Derivative::centralDifference(sine, x, 1e-5);

    EXPECT_NEAR(coarse, expected, 1e-4);
    EXPECT_NEAR(fine, expected, kTightTolerance);
    EXPECT_LT(std::abs(fine - expected), std::abs(coarse - expected));
}

TEST(DerivativeNumericTest, SupportsVerySmallFiniteResults) {
    const auto tinyScale = [](double x) {
        return 1e-150 * x * x;
    };
    EXPECT_NEAR(Derivative::centralDifference(tinyScale, 1.0), 2e-150, 1e-156);
}

TEST(DerivativeRobustnessTest, RejectsInvalidInputsAndSpecialFloatingPoints) {
    const std::function<double(double)> emptyFunction;
    EXPECT_THROW(Derivative::centralDifference(emptyFunction, 0.0), std::invalid_argument);
    EXPECT_THROW(Derivative::centralDifference(linear, 0.0, 0.0), std::invalid_argument);
    EXPECT_THROW(Derivative::forwardDifference(linear, 0.0, -1e-3), std::invalid_argument);
    EXPECT_THROW(Derivative::backwardDifference(linear, std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    EXPECT_THROW(Derivative::secondDerivative(linear, std::numeric_limits<double>::infinity()), std::invalid_argument);
    EXPECT_THROW(Derivative::centralDifference(linear, 0.0, std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    EXPECT_THROW(Derivative::derivatives(linear, {0.0, std::numeric_limits<double>::infinity()}), std::invalid_argument);
}

TEST(DerivativeRobustnessTest, ReportsDomainErrorsWithoutCrashing) {
    const auto logarithm = [](double x) {
        return std::log(x);
    };
    EXPECT_THROW(Derivative::centralDifference(logarithm, 1e-6), std::domain_error);

    const auto explosive = [](double x) {
        return std::exp(x);
    };
    EXPECT_THROW(Derivative::forwardDifference(explosive, 710.0), std::domain_error);

    const auto syntaxLikeFailure = [](double) -> double {
        throw std::runtime_error("syntax error");
    };
    EXPECT_THROW(Derivative::centralDifference(syntaxLikeFailure, 0.0), std::runtime_error);
}

} // namespace
