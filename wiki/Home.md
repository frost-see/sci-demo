# Sci Demo Wiki

## 求导函数功能帮助

项目提供 `sci::Derivative` 类用于数值求导，头文件路径为 `src/derivative.h`。

### 可用接口

- `centralDifference(f, x, h = 1e-5)`：中心差分，精度通常最好。  
  公式：`(f(x+h) - f(x-h)) / (2h)`
- `forwardDifference(f, x, h = 1e-5)`：前向差分。  
  公式：`(f(x+h) - f(x)) / h`
- `backwardDifference(f, x, h = 1e-5)`：后向差分。  
  公式：`(f(x) - f(x-h)) / h`
- `secondDerivative(f, x, h = 1e-5)`：二阶导数中心差分。  
  公式：`(f(x+h) - 2f(x) + f(x-h)) / (h^2)`
- `derivatives(f, points, h = 1e-5)`：批量计算多个点的一阶导（内部使用中心差分）。

### 参数说明

- `f`：`std::function<double(double)>`，待求导函数。
- `x`：求导点。
- `h`：步长，必须是有限值（非 `NaN`、非 `±inf`）且必须大于 0；默认 `1e-5`。

### 异常与边界行为

- `f` 为空、`x` 非有限值、`h` 小于等于 0 或非有限值时，抛出 `std::invalid_argument`。
- 若采样点（如 `x+h`、`x-h`）非有限，或 `f` 在对应采样点上不可计算/返回非有限值，抛出 `std::domain_error`。
- 若 `f` 本身抛出异常，异常会向上传递。

### 使用示例

```cpp
#include "derivative.h"
#include <cmath>
#include <vector>

double f(double x) { return std::sin(x); }

double d1 = sci::Derivative::centralDifference(f, 0.3);
double d2 = sci::Derivative::secondDerivative(f, 0.3);
std::vector<double> ds = sci::Derivative::derivatives(f, {0.0, 0.5, 1.0});
```

## 曲线拟合功能帮助

项目提供 `sci::CurveFitting` 类用于多项式最小二乘拟合，头文件路径为 `src/curve_fitting.h`。

### 可用接口

- `polynomialFit(x, y, degree)`：对离散点进行 `degree` 次多项式最小二乘拟合。
- `PolynomialFitResult::evaluate(x)`：使用拟合结果预测指定点的函数值。

### 返回结果

- `coefficients`：按常数项到高次项顺序存放的多项式系数。
- `rSquared`：拟合优度，范围通常为 `(-∞, 1]`，越接近 `1` 表示拟合越好。

### 参数说明

- `x`、`y`：输入样本点坐标，长度必须相同且不能为空。
- `degree`：拟合多项式次数，要求样本点数量严格大于 `degree`。

### 异常与边界行为

- `x` / `y` 为空、长度不一致、包含非有限值，或样本点数量不足时，抛出 `std::invalid_argument`。
- 若法方程不可解（例如所有 `x` 完全相同导致矩阵奇异），抛出 `std::domain_error`。
- 对空拟合结果调用 `evaluate` 会抛出 `std::logic_error`；传入非有限 `x` 会抛出 `std::invalid_argument`。

### 使用示例

```cpp
#include "curve_fitting.h"
#include <vector>

std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
std::vector<double> y = {1.0, 3.0, 5.0, 7.0};

auto fit = sci::CurveFitting::polynomialFit(x, y, 1);
double predicted = fit.evaluate(4.0);  // 9.0
```
