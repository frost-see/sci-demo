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
- `h`：步长，必须是有限值且必须大于 0；默认 `1e-5`。

### 异常与边界行为

- `f` 为空、`x` 非有限值、`h` 小于等于 0 或非有限值时，抛出 `std::invalid_argument`。
- 若采样点越界或函数返回非有限值，抛出 `std::domain_error`。
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
