// Unit tests for sci::CurveFit
//
// Test design methods used:
//   - Equivalence class partitioning (valid and invalid input classes)
//   - Boundary value analysis (degree limits, minimum data points)
//   - Error guessing (degenerate data, edge numerics)

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "../src/curve_fit.h"

namespace {

using sci::CurveFit;

constexpr double kTightTolerance = 1e-8;
constexpr double kLooseTolerance = 1e-6;

// ---------------------------------------------------------------------------
// Helpers: generate (x, y) samples from a known polynomial
// ---------------------------------------------------------------------------

// Evaluate polynomial given as coefficient vector via Horner's method.
double evalPoly(const std::vector<double>& c, double x) {
    double r = 0.0;
    for (int i = static_cast<int>(c.size()) - 1; i >= 0; --i) {
        r = r * x + c[i];
    }
    return r;
}

// Build n equally-spaced samples over [xMin, xMax] from polynomial c.
std::pair<std::vector<double>, std::vector<double>> makeSamples(
    const std::vector<double>& c,
    double xMin,
    double xMax,
    int n) {
    std::vector<double> xs(n), ys(n);
    for (int i = 0; i < n; ++i) {
        xs[i] = xMin + (xMax - xMin) * i / (n - 1);
        ys[i] = evalPoly(c, xs[i]);
    }
    return {xs, ys};
}

// ---------------------------------------------------------------------------
// CurveFitMathTest – equivalence class: valid inputs, correct coefficients
// ---------------------------------------------------------------------------

TEST(CurveFitMathTest, FitsQuadraticExactly) {
    // p(x) = 1 + 2x + 3x^2
    const std::vector<double> expected = {1.0, 2.0, 3.0};
    auto [xs, ys] = makeSamples(expected, -2.0, 2.0, 10);
    const auto coeffs = CurveFit::fit(xs, ys, 2);
    ASSERT_EQ(coeffs.size(), 3u);
    EXPECT_NEAR(coeffs[0], 1.0, kTightTolerance);
    EXPECT_NEAR(coeffs[1], 2.0, kTightTolerance);
    EXPECT_NEAR(coeffs[2], 3.0, kTightTolerance);
}

TEST(CurveFitMathTest, FitsCubicExactly) {
    // p(x) = -1 + 0.5x - 2x^2 + x^3
    const std::vector<double> expected = {-1.0, 0.5, -2.0, 1.0};
    auto [xs, ys] = makeSamples(expected, -3.0, 3.0, 15);
    const auto coeffs = CurveFit::fit(xs, ys, 3);
    ASSERT_EQ(coeffs.size(), 4u);
    EXPECT_NEAR(coeffs[0], -1.0, kTightTolerance);
    EXPECT_NEAR(coeffs[1],  0.5, kTightTolerance);
    EXPECT_NEAR(coeffs[2], -2.0, kTightTolerance);
    EXPECT_NEAR(coeffs[3],  1.0, kTightTolerance);
}

TEST(CurveFitMathTest, FitsQuarticExactly) {
    // p(x) = 0 + x + x^2 - x^3 + 0.5*x^4
    const std::vector<double> expected = {0.0, 1.0, 1.0, -1.0, 0.5};
    auto [xs, ys] = makeSamples(expected, -2.0, 2.0, 20);
    const auto coeffs = CurveFit::fit(xs, ys, 4);
    ASSERT_EQ(coeffs.size(), 5u);
    for (int i = 0; i < 5; ++i) {
        EXPECT_NEAR(coeffs[i], expected[i], kLooseTolerance);
    }
}

TEST(CurveFitMathTest, FitsQuinticExactly) {
    // p(x) = 1 - x + x^2 + x^3 - x^4 + 0.1*x^5
    const std::vector<double> expected = {1.0, -1.0, 1.0, 1.0, -1.0, 0.1};
    auto [xs, ys] = makeSamples(expected, -1.0, 1.0, 25);
    const auto coeffs = CurveFit::fit(xs, ys, 5);
    ASSERT_EQ(coeffs.size(), 6u);
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(coeffs[i], expected[i], kLooseTolerance);
    }
}

TEST(CurveFitMathTest, OverdeterminedSystemMinimisesResiduals) {
    // Data lies exactly on y = 2 + 3x^2; fitting should recover coefficients.
    const std::vector<double> xs = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
    std::vector<double> ys;
    for (double xi : xs) {
        ys.push_back(2.0 + 3.0 * xi * xi);
    }
    const auto coeffs = CurveFit::fit(xs, ys, 2);
    EXPECT_NEAR(coeffs[0], 2.0, kTightTolerance);
    EXPECT_NEAR(coeffs[1], 0.0, kTightTolerance);
    EXPECT_NEAR(coeffs[2], 3.0, kTightTolerance);
}

TEST(CurveFitMathTest, FitWithNegativeCoefficients) {
    // p(x) = -5 - 3x + 2x^2
    const std::vector<double> expected = {-5.0, -3.0, 2.0};
    auto [xs, ys] = makeSamples(expected, 0.0, 4.0, 12);
    const auto coeffs = CurveFit::fit(xs, ys, 2);
    EXPECT_NEAR(coeffs[0], -5.0, kTightTolerance);
    EXPECT_NEAR(coeffs[1], -3.0, kTightTolerance);
    EXPECT_NEAR(coeffs[2],  2.0, kTightTolerance);
}

TEST(CurveFitMathTest, FitWithZeroConstantTerm) {
    // p(x) = x + x^2 (constant term is zero)
    const std::vector<double> expected = {0.0, 1.0, 1.0};
    auto [xs, ys] = makeSamples(expected, -5.0, 5.0, 20);
    const auto coeffs = CurveFit::fit(xs, ys, 2);
    EXPECT_NEAR(coeffs[0], 0.0, kTightTolerance);
    EXPECT_NEAR(coeffs[1], 1.0, kTightTolerance);
    EXPECT_NEAR(coeffs[2], 1.0, kTightTolerance);
}

// ---------------------------------------------------------------------------
// CurveFitBoundaryTest – boundary value analysis
// ---------------------------------------------------------------------------

TEST(CurveFitBoundaryTest, MinimumDegree2WithMinimumPoints) {
    // Degree 2 requires at least 3 data points (degree + 1).
    const std::vector<double> xs = {0.0, 1.0, 2.0};
    const std::vector<double> ys = {0.0, 1.0, 4.0};  // y = x^2
    const auto coeffs = CurveFit::fit(xs, ys, 2);
    ASSERT_EQ(coeffs.size(), 3u);
    EXPECT_NEAR(coeffs[0], 0.0, kLooseTolerance);
    EXPECT_NEAR(coeffs[1], 0.0, kLooseTolerance);
    EXPECT_NEAR(coeffs[2], 1.0, kLooseTolerance);
}

TEST(CurveFitBoundaryTest, MaximumDegree5WithMinimumPoints) {
    // Degree 5 requires at least 6 data points.
    const std::vector<double> expected = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    auto [xs, ys] = makeSamples(expected, -2.0, 2.0, 6);
    const auto coeffs = CurveFit::fit(xs, ys, 5);
    ASSERT_EQ(coeffs.size(), 6u);
    for (int i = 0; i < 6; ++i) {
        EXPECT_NEAR(coeffs[i], expected[i], kLooseTolerance);
    }
}

TEST(CurveFitBoundaryTest, MinimumValidDegreeIs2) {
    const std::vector<double> xs = {0.0, 1.0, 2.0, 3.0};
    const std::vector<double> ys = {1.0, 2.0, 5.0, 10.0};
    EXPECT_NO_THROW(CurveFit::fit(xs, ys, 2));
}

TEST(CurveFitBoundaryTest, MaximumValidDegreeIs5) {
    const std::vector<double> expected = {0.0, 1.0, -1.0, 0.5, -0.1, 0.01};
    auto [xs, ys] = makeSamples(expected, -1.0, 1.0, 10);
    EXPECT_NO_THROW(CurveFit::fit(xs, ys, 5));
}

TEST(CurveFitBoundaryTest, ExactlyDegree1AboveMinimumPoints) {
    // degree=3 with exactly 4 points (degree+1): just enough
    const std::vector<double> expected = {1.0, 2.0, 3.0, 4.0};
    auto [xs, ys] = makeSamples(expected, 0.0, 3.0, 4);
    const auto coeffs = CurveFit::fit(xs, ys, 3);
    ASSERT_EQ(coeffs.size(), 4u);
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(coeffs[i], expected[i], kLooseTolerance);
    }
}

// ---------------------------------------------------------------------------
// CurveFitEvaluateTest – evaluate() correctness
// ---------------------------------------------------------------------------

TEST(CurveFitEvaluateTest, EvaluatesConstantPolynomial) {
    const std::vector<double> coeffs = {7.0};
    EXPECT_NEAR(CurveFit::evaluate(coeffs, 0.0), 7.0, kTightTolerance);
    EXPECT_NEAR(CurveFit::evaluate(coeffs, 100.0), 7.0, kTightTolerance);
    EXPECT_NEAR(CurveFit::evaluate(coeffs, -50.0), 7.0, kTightTolerance);
}

TEST(CurveFitEvaluateTest, EvaluatesLinearPolynomial) {
    // p(x) = 2 + 3x
    const std::vector<double> coeffs = {2.0, 3.0};
    EXPECT_NEAR(CurveFit::evaluate(coeffs, 0.0), 2.0, kTightTolerance);
    EXPECT_NEAR(CurveFit::evaluate(coeffs, 1.0), 5.0, kTightTolerance);
    EXPECT_NEAR(CurveFit::evaluate(coeffs, -1.0), -1.0, kTightTolerance);
}

TEST(CurveFitEvaluateTest, EvaluatesQuadraticPolynomial) {
    // p(x) = 1 + 2x + 3x^2
    const std::vector<double> coeffs = {1.0, 2.0, 3.0};
    EXPECT_NEAR(CurveFit::evaluate(coeffs, 0.0), 1.0, kTightTolerance);
    EXPECT_NEAR(CurveFit::evaluate(coeffs, 1.0), 6.0, kTightTolerance);
    EXPECT_NEAR(CurveFit::evaluate(coeffs, -1.0), 2.0, kTightTolerance);
    EXPECT_NEAR(CurveFit::evaluate(coeffs, 2.0), 17.0, kTightTolerance);
}

TEST(CurveFitEvaluateTest, EvaluatesMatchesRoundTrip) {
    // Fit then evaluate should reproduce y values closely.
    const std::vector<double> expected = {1.0, -2.0, 3.0};
    auto [xs, ys] = makeSamples(expected, -2.0, 2.0, 10);
    const auto coeffs = CurveFit::fit(xs, ys, 2);
    for (std::size_t i = 0; i < xs.size(); ++i) {
        EXPECT_NEAR(CurveFit::evaluate(coeffs, xs[i]), ys[i], kLooseTolerance);
    }
}

TEST(CurveFitEvaluateTest, EvaluatesAtNegativeX) {
    // p(x) = x^2: p(-3) = 9
    const std::vector<double> coeffs = {0.0, 0.0, 1.0};
    EXPECT_NEAR(CurveFit::evaluate(coeffs, -3.0), 9.0, kTightTolerance);
}

TEST(CurveFitEvaluateTest, EvaluatesAtZero) {
    // p(0) = a0 regardless of other terms.
    const std::vector<double> coeffs = {5.0, -2.0, 1.0, 3.0};
    EXPECT_NEAR(CurveFit::evaluate(coeffs, 0.0), 5.0, kTightTolerance);
}

// ---------------------------------------------------------------------------
// CurveFitInvalidInputTest – equivalence class: invalid inputs
// ---------------------------------------------------------------------------

TEST(CurveFitInvalidInputTest, RejectsDegreeZero) {
    const std::vector<double> xs = {0.0, 1.0, 2.0};
    const std::vector<double> ys = {1.0, 2.0, 3.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 0), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsDegreeOne) {
    const std::vector<double> xs = {0.0, 1.0, 2.0};
    const std::vector<double> ys = {1.0, 2.0, 3.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 1), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsDegreeAboveFive) {
    const std::vector<double> xs = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    const std::vector<double> ys = {0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 6), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsNegativeDegree) {
    const std::vector<double> xs = {0.0, 1.0, 2.0};
    const std::vector<double> ys = {1.0, 2.0, 3.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, -1), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsEmptyXVector) {
    const std::vector<double> xs;
    const std::vector<double> ys;
    EXPECT_THROW(CurveFit::fit(xs, ys, 2), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsMismatchedVectorSizes) {
    const std::vector<double> xs = {0.0, 1.0, 2.0};
    const std::vector<double> ys = {1.0, 2.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 2), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsInsufficientDataPointsForDegree2) {
    // Degree 2 requires more than 2 points; exactly 2 should fail.
    const std::vector<double> xs = {0.0, 1.0};
    const std::vector<double> ys = {0.0, 1.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 2), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsInsufficientDataPointsForDegree5) {
    // Degree 5 requires more than 5 points; exactly 5 should fail.
    const std::vector<double> xs = {0.0, 1.0, 2.0, 3.0, 4.0};
    const std::vector<double> ys = {0.0, 1.0, 4.0, 9.0, 16.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 5), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsNaNInX) {
    const std::vector<double> xs = {0.0, std::numeric_limits<double>::quiet_NaN(), 2.0};
    const std::vector<double> ys = {0.0, 1.0, 4.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 2), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsInfinityInX) {
    const std::vector<double> xs = {0.0, std::numeric_limits<double>::infinity(), 2.0};
    const std::vector<double> ys = {0.0, 1.0, 4.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 2), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsNaNInY) {
    const std::vector<double> xs = {0.0, 1.0, 2.0};
    const std::vector<double> ys = {0.0, std::numeric_limits<double>::quiet_NaN(), 4.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 2), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsNegativeInfinityInY) {
    const std::vector<double> xs = {0.0, 1.0, 2.0};
    const std::vector<double> ys = {0.0, -std::numeric_limits<double>::infinity(), 4.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 2), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, RejectsAllIdenticalXValues) {
    // All x equal → Vandermonde matrix is singular.
    const std::vector<double> xs = {1.0, 1.0, 1.0, 1.0};
    const std::vector<double> ys = {2.0, 3.0, 4.0, 5.0};
    EXPECT_THROW(CurveFit::fit(xs, ys, 2), std::domain_error);
}

TEST(CurveFitInvalidInputTest, EvaluateRejectsEmptyCoefficients) {
    const std::vector<double> coeffs;
    EXPECT_THROW(CurveFit::evaluate(coeffs, 1.0), std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, EvaluateRejectsNaNX) {
    const std::vector<double> coeffs = {1.0, 2.0, 3.0};
    EXPECT_THROW(
        CurveFit::evaluate(coeffs, std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
}

TEST(CurveFitInvalidInputTest, EvaluateRejectsInfiniteX) {
    const std::vector<double> coeffs = {1.0, 2.0, 3.0};
    EXPECT_THROW(
        CurveFit::evaluate(coeffs, std::numeric_limits<double>::infinity()),
        std::invalid_argument);
}

// ---------------------------------------------------------------------------
// CurveFitAccuracyTest – numerical accuracy with noisy-adjacent exact data
// ---------------------------------------------------------------------------

TEST(CurveFitAccuracyTest, PredictedValuesCloseToTrueForAllDegrees) {
    for (int degree = 2; degree <= 5; ++degree) {
        std::vector<double> trueCoeffs(degree + 1);
        for (int i = 0; i <= degree; ++i) {
            trueCoeffs[i] = (i % 2 == 0) ? 1.0 : -0.5;
        }
        auto [xs, ys] = makeSamples(trueCoeffs, -1.0, 1.0, 30);
        const auto coeffs = CurveFit::fit(xs, ys, degree);
        // Verify fitted polynomial reproduces training data to high accuracy.
        for (std::size_t i = 0; i < xs.size(); ++i) {
            EXPECT_NEAR(CurveFit::evaluate(coeffs, xs[i]), ys[i], kLooseTolerance)
                << "degree=" << degree << " point index=" << i;
        }
    }
}

TEST(CurveFitAccuracyTest, CoefficientsRecoveredForAllValidDegrees) {
    for (int degree = 2; degree <= 5; ++degree) {
        std::vector<double> trueCoeffs(degree + 1, 1.0);
        auto [xs, ys] = makeSamples(trueCoeffs, -1.0, 1.0, 30);
        const auto fitted = CurveFit::fit(xs, ys, degree);
        ASSERT_EQ(fitted.size(), static_cast<std::size_t>(degree + 1));
        for (int i = 0; i <= degree; ++i) {
            EXPECT_NEAR(fitted[i], trueCoeffs[i], kLooseTolerance)
                << "degree=" << degree << " coeff index=" << i;
        }
    }
}

TEST(CurveFitAccuracyTest, FittedPolynomialHasSmallResiduals) {
    // Generate data from a known cubic and verify sum-of-squared residuals is tiny.
    const std::vector<double> trueCoeffs = {2.0, -1.0, 0.5, -0.25};
    auto [xs, ys] = makeSamples(trueCoeffs, -2.0, 2.0, 20);
    const auto coeffs = CurveFit::fit(xs, ys, 3);

    double ssr = 0.0;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        const double residual = CurveFit::evaluate(coeffs, xs[i]) - ys[i];
        ssr += residual * residual;
    }
    EXPECT_NEAR(ssr, 0.0, 1e-20);
}

} // namespace
