#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

#include <gtest/gtest.h>

#include "../src/integration.h"

namespace {

using sci::Integration;

constexpr double kTrapezoidalTolerance = 1e-4;
constexpr double kSimpsonTolerance = 1e-7;
constexpr double kPi = 3.14159265358979323846;

double quadratic(double x) {
    return x * x;
}

double sine(double x) {
    return std::sin(x);
}

TEST(IntegrationMathTest, IntegratesBasicFunctions) {
    EXPECT_NEAR(Integration::trapezoidal(quadratic, 0.0, 1.0), 1.0 / 3.0, kTrapezoidalTolerance);
    EXPECT_NEAR(Integration::simpson(quadratic, 0.0, 1.0), 1.0 / 3.0, kSimpsonTolerance);

    EXPECT_NEAR(Integration::trapezoidal(sine, 0.0, kPi), 2.0, kTrapezoidalTolerance);
    EXPECT_NEAR(Integration::simpson(sine, 0.0, kPi), 2.0, kSimpsonTolerance);
}

TEST(IntegrationMathTest, HandlesReversedAndDegenerateBounds) {
    EXPECT_NEAR(Integration::trapezoidal(quadratic, 1.0, 0.0), -1.0 / 3.0, kTrapezoidalTolerance);
    EXPECT_NEAR(Integration::simpson(quadratic, 1.0, 0.0), -1.0 / 3.0, kSimpsonTolerance);

    EXPECT_DOUBLE_EQ(Integration::trapezoidal(quadratic, 2.0, 2.0), 0.0);
    EXPECT_DOUBLE_EQ(Integration::simpson(quadratic, 2.0, 2.0), 0.0);
}

TEST(IntegrationRobustnessTest, RejectsInvalidInputs) {
    const std::function<double(double)> emptyFunction;
    EXPECT_THROW(Integration::trapezoidal(emptyFunction, 0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(Integration::trapezoidal(quadratic, std::numeric_limits<double>::quiet_NaN(), 1.0), std::invalid_argument);
    EXPECT_THROW(Integration::trapezoidal(quadratic, 0.0, 1.0, 0), std::invalid_argument);

    EXPECT_THROW(Integration::simpson(quadratic, 0.0, 1.0, 3), std::invalid_argument);
    EXPECT_THROW(Integration::simpson(quadratic, 0.0, std::numeric_limits<double>::infinity()), std::invalid_argument);
}

TEST(IntegrationRobustnessTest, ReportsDomainErrors) {
    const auto logarithm = [](double x) {
        return std::log(x);
    };
    EXPECT_THROW(Integration::trapezoidal(logarithm, -1.0, 1.0), std::domain_error);
    EXPECT_THROW(Integration::simpson(logarithm, -1.0, 1.0), std::domain_error);
}

} // namespace
