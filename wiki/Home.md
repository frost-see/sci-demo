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

---

## 曲线拟合功能帮助

项目提供 `sci::CurveFit` 类用于最小二乘多项式拟合，头文件路径为 `src/curve_fit.h`。

### 可用接口

- `fit(x, y, degree)`：对数据点 `(x, y)` 进行 `degree` 次多项式拟合。  
  返回系数向量 `[a0, a1, ..., an]`，对应多项式 `p(x) = a0 + a1*x + ... + an*x^n`。
- `evaluate(coefficients, x)`：用给定系数向量在点 `x` 处求多项式值（Horner 法）。

### 参数说明

- `x`、`y`：`std::vector<double>`，数据点横纵坐标，长度必须相等且大于 `degree`。
- `degree`：多项式次数，必须在 `[2, 5]` 之间（含边界）。
- `coefficients`：`fit()` 返回的系数向量，或用户自定义的系数。

### 异常与边界行为

- `degree` 不在 `[2, 5]` 范围内时，抛出 `std::invalid_argument`。
- `x` 或 `y` 为空，或长度不匹配，或数据点数不超过 `degree` 时，抛出 `std::invalid_argument`。
- `x` 或 `y` 中存在 `NaN`、`±inf` 等非有限值时，抛出 `std::invalid_argument`。
- 数据退化（如所有 `x` 值相同）导致法方程奇异时，抛出 `std::domain_error`。
- `evaluate()` 中 `coefficients` 为空或 `x` 非有限时，抛出 `std::invalid_argument`。

### 使用示例

```cpp
#include "curve_fit.h"
#include <vector>

// 数据点
std::vector<double> xs = {-2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
std::vector<double> ys = { 7.0,  2.0, 1.0, 4.0, 11.0, 22.0};

// 拟合二次多项式
std::vector<double> coeffs = sci::CurveFit::fit(xs, ys, 2);
// coeffs[0] = a0, coeffs[1] = a1, coeffs[2] = a2

// 在新点上求值
double predicted = sci::CurveFit::evaluate(coeffs, 1.5);
```

### 端到端验证示例

`examples/e2e_curve_fit.cpp` 提供了完整的端到端验证程序，包含 4 个真实物理场景：

| 场景 | 多项式次数 | 描述 |
|---|---|---|
| 抛体高度 | 2次 | `h(t) = h0 + v0·t - 0.5·g·t²` |
| 热敏电阻阻值 | 3次 | 温度与电阻的三次多项式近似 |
| 简支梁挠度 | 4次 | 分布载荷下梁的四次挠度曲线 |
| 阻尼振荡器位移 | 5次 | `e^{-t}` 泰勒展开前六项 |

程序对每个场景计算 R² 和最大残差，当 R² > 0.9999 且最大残差 < 1e-6 时判定为通过，否则以非零退出码结束（适合作为 CI 验收门禁）。

编译与运行：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target e2e_curve_fit
./build/examples/e2e_curve_fit
```
