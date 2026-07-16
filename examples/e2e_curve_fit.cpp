// End-to-end validation for sci::CurveFit
//
// This program exercises the full curve fitting pipeline using realistic
// physical scenarios for each supported polynomial degree (2–5).  It prints
// a human-readable report and exits with a non-zero status if any scenario
// fails its acceptance criterion, making it suitable as a CI smoke test.
//
// Scenarios
//   1. Quadratic  – projectile height vs. time  (h = h0 + v0*t - 0.5*g*t^2)
//   2. Cubic      – thermistor resistance vs. temperature
//   3. Quartic    – deflection of a simply-supported beam under distributed load
//   4. Quintic    – displacement of a damped oscillator at discrete time steps
//
// For every scenario the program:
//   a) generates noise-free sample points from a known polynomial,
//   b) fits a polynomial of the matching degree,
//   c) computes R² and max absolute residual,
//   d) prints the report, and
//   e) asserts that R² > 0.9999 and max_residual < 1e-6.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "curve_fit.h"

namespace {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Evaluate a polynomial (Horner) given coefficient vector c[0..n].
double evalPoly(const std::vector<double>& c, double x) {
    double r = 0.0;
    for (int i = static_cast<int>(c.size()) - 1; i >= 0; --i) {
        r = r * x + c[i];
    }
    return r;
}

// Build n equally-spaced points over [xMin, xMax] from polynomial trueCoeffs.
std::pair<std::vector<double>, std::vector<double>> makeSamples(
    const std::vector<double>& trueCoeffs,
    double xMin, double xMax, int n) {

    std::vector<double> xs(n), ys(n);
    for (int i = 0; i < n; ++i) {
        xs[i] = xMin + (xMax - xMin) * i / (n - 1);
        ys[i] = evalPoly(trueCoeffs, xs[i]);
    }
    return {xs, ys};
}

// Compute coefficient of determination R².
double computeR2(const std::vector<double>& yTrue,
                 const std::vector<double>& yPred) {
    const int n = static_cast<int>(yTrue.size());
    double mean = 0.0;
    for (double v : yTrue) mean += v;
    mean /= n;

    double ssTot = 0.0;
    double ssRes = 0.0;
    for (int i = 0; i < n; ++i) {
        ssTot += (yTrue[i] - mean) * (yTrue[i] - mean);
        ssRes += (yTrue[i] - yPred[i]) * (yTrue[i] - yPred[i]);
    }
    if (ssTot < 1e-300) return 1.0;  // constant function
    return 1.0 - ssRes / ssTot;
}

// ---------------------------------------------------------------------------
// Scenario runner
// ---------------------------------------------------------------------------

struct ScenarioResult {
    std::string name;
    int degree;
    std::vector<double> trueCoeffs;
    std::vector<double> fittedCoeffs;
    double r2;
    double maxResidual;
    bool passed;
};

ScenarioResult runScenario(
    const std::string& name,
    int degree,
    const std::vector<double>& trueCoeffs,
    double xMin, double xMax, int nPoints) {

    auto [xs, ys] = makeSamples(trueCoeffs, xMin, xMax, nPoints);
    const auto fitted = sci::CurveFit::fit(xs, ys, degree);

    std::vector<double> predicted(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        predicted[i] = sci::CurveFit::evaluate(fitted, xs[i]);
    }

    const double r2 = computeR2(ys, predicted);

    double maxRes = 0.0;
    for (int i = 0; i < nPoints; ++i) {
        maxRes = std::max(maxRes, std::abs(ys[i] - predicted[i]));
    }

    const bool passed = (r2 > 0.9999) && (maxRes < 1e-6);
    return {name, degree, trueCoeffs, fitted, r2, maxRes, passed};
}

// ---------------------------------------------------------------------------
// Report printer
// ---------------------------------------------------------------------------

void printReport(const ScenarioResult& r) {
    std::cout << "==========================================================\n";
    std::cout << "Scenario : " << r.name << "\n";
    std::cout << "Degree   : " << r.degree << "\n";

    std::cout << std::fixed << std::setprecision(8);

    std::cout << "True coefficients   : [";
    for (std::size_t i = 0; i < r.trueCoeffs.size(); ++i) {
        std::cout << r.trueCoeffs[i];
        if (i + 1 < r.trueCoeffs.size()) std::cout << ", ";
    }
    std::cout << "]\n";

    std::cout << "Fitted coefficients : [";
    for (std::size_t i = 0; i < r.fittedCoeffs.size(); ++i) {
        std::cout << r.fittedCoeffs[i];
        if (i + 1 < r.fittedCoeffs.size()) std::cout << ", ";
    }
    std::cout << "]\n";

    std::cout << "R²               : " << r.r2 << "\n";
    std::cout << "Max residual     : " << r.maxResidual << "\n";
    std::cout << "Result           : " << (r.passed ? "PASS" : "FAIL") << "\n";
}

} // namespace

int main() {
    // -----------------------------------------------------------------------
    // Scenario 1 – Quadratic: projectile height
    //   h(t) = 10 + 20t - 4.905t²
    //   (h0=10 m, v0=20 m/s, g=9.81 m/s²)
    // -----------------------------------------------------------------------
    const auto s1 = runScenario(
        "Projectile height (quadratic, degree=2)",
        2,
        {10.0, 20.0, -4.905},
        0.0, 4.0, 20);

    // -----------------------------------------------------------------------
    // Scenario 2 – Cubic: thermistor resistance approximation
    //   R(T) = 100 - 0.5T + 0.002T² - 0.000003T³
    // -----------------------------------------------------------------------
    const auto s2 = runScenario(
        "Thermistor resistance (cubic, degree=3)",
        3,
        {100.0, -0.5, 0.002, -0.000003},
        0.0, 150.0, 25);

    // -----------------------------------------------------------------------
    // Scenario 3 – Quartic: beam deflection profile
    //   y(x) = 0.001*(x^4 - 2*L*x^3 + L^2*x^2)  with L=1
    //        = 0.001x^4 - 0.002x^3 + 0.001x^2
    // -----------------------------------------------------------------------
    const auto s3 = runScenario(
        "Beam deflection (quartic, degree=4)",
        4,
        {0.0, 0.0, 0.001, -0.002, 0.001},
        0.0, 1.0, 20);

    // -----------------------------------------------------------------------
    // Scenario 4 – Quintic: damped oscillator displacement polynomial
    //   x(t) = 1 - t + 0.5t² - (1/6)t³ + (1/24)t^4 - (1/120)t^5
    //   (first six terms of e^{-t} Taylor series)
    // -----------------------------------------------------------------------
    const auto s4 = runScenario(
        "Damped oscillator (quintic, degree=5)",
        5,
        {1.0, -1.0, 0.5, -1.0 / 6.0, 1.0 / 24.0, -1.0 / 120.0},
        0.0, 2.0, 30);

    // -----------------------------------------------------------------------
    // Print reports
    // -----------------------------------------------------------------------
    const std::vector<ScenarioResult> results = {s1, s2, s3, s4};
    for (const auto& r : results) {
        printReport(r);
    }

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------
    int passed = 0;
    int failed = 0;
    std::cout << "\n==================== SUMMARY ====================\n";
    for (const auto& r : results) {
        std::cout << (r.passed ? "[PASS]" : "[FAIL]")
                  << " " << r.name << "\n";
        if (r.passed) ++passed; else ++failed;
    }
    std::cout << passed << "/" << (passed + failed)
              << " scenarios passed.\n";

    return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
