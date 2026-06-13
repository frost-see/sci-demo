#include <cassert>
#include <cmath>
#include <iostream>
#include "../src/derivative.h"

using namespace sci;

// Test functions
double linearFunc(double x) {
    return 2.0 * x + 3.0;  // derivative = 2.0
}

double quadraticFunc(double x) {
    return x * x;  // derivative = 2x
}

double cubicFunc(double x) {
    return x * x * x;  // derivative = 3x^2
}

double sinFunc(double x) {
    return std::sin(x);  // derivative = cos(x)
}

double expFunc(double x) {
    return std::exp(x);  // derivative = exp(x)
}

void testCentralDifference() {
    std::cout << "Testing central difference method..." << std::endl;
    
    // Test linear function (derivative should be 2.0)
    double deriv = Derivative::centralDifference(linearFunc, 1.0);
    assert(std::abs(deriv - 2.0) < 1e-4);
    std::cout << "  Linear function: PASSED" << std::endl;
    
    // Test quadratic function at x=2 (derivative should be 4.0)
    deriv = Derivative::centralDifference(quadraticFunc, 2.0);
    assert(std::abs(deriv - 4.0) < 1e-4);
    std::cout << "  Quadratic function: PASSED" << std::endl;
    
    // Test sine function at x=0 (derivative should be 1.0)
    deriv = Derivative::centralDifference(sinFunc, 0.0);
    assert(std::abs(deriv - 1.0) < 1e-4);
    std::cout << "  Sine function: PASSED" << std::endl;
}

void testForwardDifference() {
    std::cout << "Testing forward difference method..." << std::endl;
    
    // Test linear function (derivative should be 2.0)
    double deriv = Derivative::forwardDifference(linearFunc, 1.0);
    assert(std::abs(deriv - 2.0) < 1e-3);
    std::cout << "  Linear function: PASSED" << std::endl;
}

void testBackwardDifference() {
    std::cout << "Testing backward difference method..." << std::endl;
    
    // Test linear function (derivative should be 2.0)
    double deriv = Derivative::backwardDifference(linearFunc, 1.0);
    assert(std::abs(deriv - 2.0) < 1e-3);
    std::cout << "  Linear function: PASSED" << std::endl;
}

void testSecondDerivative() {
    std::cout << "Testing second derivative..." << std::endl;
    
    // Test quadratic function (second derivative should be 2.0)
    double deriv = Derivative::secondDerivative(quadraticFunc, 1.0);
    assert(std::abs(deriv - 2.0) < 1e-3);
    std::cout << "  Quadratic function: PASSED" << std::endl;
    
    // Test cubic function at x=1 (second derivative should be 6.0)
    deriv = Derivative::secondDerivative(cubicFunc, 1.0);
    assert(std::abs(deriv - 6.0) < 1e-2);
    std::cout << "  Cubic function: PASSED" << std::endl;
}

void testDerivatives() {
    std::cout << "Testing derivatives at multiple points..." << std::endl;
    
    std::vector<double> points = {0.0, 1.0, 2.0};
    std::vector<double> derivs = Derivative::derivatives(quadraticFunc, points);
    
    assert(derivs.size() == 3);
    assert(std::abs(derivs[0] - 0.0) < 1e-3);   // 2*0 = 0
    assert(std::abs(derivs[1] - 2.0) < 1e-3);   // 2*1 = 2
    assert(std::abs(derivs[2] - 4.0) < 1e-3);   // 2*2 = 4
    std::cout << "  Multiple points: PASSED" << std::endl;
}

int main() {
    std::cout << "Running derivative tests...\n" << std::endl;
    
    try {
        testCentralDifference();
        testForwardDifference();
        testBackwardDifference();
        testSecondDerivative();
        testDerivatives();
        
        std::cout << "\nAll tests passed!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cout << "Test failed: " << e.what() << std::endl;
        return 1;
    }
}
