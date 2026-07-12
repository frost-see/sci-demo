#pragma once

#include <functional>

namespace sci {

class Integration {
public:
    static double trapezoidal(std::function<double(double)> f, double a, double b, int intervals = 1000);
    static double simpson(std::function<double(double)> f, double a, double b, int intervals = 1000);
};

} // namespace sci
