#pragma once

#include <functional>

namespace sci {

class Integration {
public:
    static double trapezoidal(const std::function<double(double)>& f, double a, double b, int intervals = 1000);
    static double simpson(const std::function<double(double)>& f, double a, double b, int intervals = 1000);
};

} // namespace sci
