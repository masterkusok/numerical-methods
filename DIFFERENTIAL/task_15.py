import math
import matplotlib.pyplot as plt


def explicit_rk(f, x0, y0, x_end, h, butcher_table):
    a, b, c = butcher_table
    s = len(b)
    
    x_vals, y_vals = [x0], [y0]
    x, y = x0, y0
    
    while x < x_end - h/2:
        k = [0] * s
        for i in range(s):
            sum_a = sum(a[i][j] * k[j] for j in range(i))
            k[i] = f(x + c[i] * h, y + h * sum_a)
        
        y += h * sum(b[i] * k[i] for i in range(s))
        x += h
        x_vals.append(x)
        y_vals.append(y)
    
    return x_vals, y_vals


def implicit_rk(f, x0, y0, x_end, h, butcher_table):
    a, b, c = butcher_table
    s = len(b)
    
    x_vals, y_vals = [x0], [y0]
    x, y = x0, y0
    
    while x < x_end - h/2:
        k = [f(x, y)] * s
        
        for _ in range(100):
            k_new = [0] * s
            for i in range(s):
                sum_a = sum(a[i][j] * k[j] for j in range(s))
                k_new[i] = f(x + c[i] * h, y + h * sum_a)
            
            if all(abs(k_new[i] - k[i]) < 1e-10 for i in range(s)):
                break
            k = k_new
        
        y += h * sum(b[i] * k[i] for i in range(s))
        x += h
        x_vals.append(x)
        y_vals.append(y)
    
    return x_vals, y_vals


def ode_function(x, y):
    return (x**4 * y**2 - x**2 * y - 20) / x**3

def analytical_solution(x):
    return (4-5*(x**9))/((x**2)+(x**11))

EULER_METHOD = ([[0]], [1], [0])
EULER_METHOD_CAUCHY = ([[0, 0], [1, 0]], [0.5, 0.5], [0, 1])
RC_4 = ([[0, 0, 0, 0], [0.5, 0, 0, 0], [0, 0.5, 0, 0], [0, 0, 1, 0]], [1/6, 1/3, 1/3, 1/6], [0, 0.5, 0.5, 1])

RADO_IIA = ([[1]], [1], [1])
MIDPOINT = ([[0.5]], [1], [0.5])
GAUSS_4_ORDER = ([[0.25, 0.25 - math.sqrt(3)/6], [0.25 + math.sqrt(3)/6, 0.25]], [0.5, 0.5], [0.5 - math.sqrt(3)/6, 0.5 + math.sqrt(3)/6])


def main():
    x0, y0, x_end, h = 1.0, -0.5, 3.0, 0.1
    
    print("РЕШЕНИЕ ЗАДАЧИ КОШИ")
    print(f"Уравнение: x³y' - x⁴y² + x²y + 20 = 0")
    print(f"Начальные условия: y({x0}) = {y0}")
    print(f"Интервал: [{x0}, {x_end}], шаг h = {h}")
    print()
    
    x_analyt = [x0 + i * h for i in range(int((x_end - x0) / h) + 1)]
    y_analyt = [analytical_solution(x) for x in x_analyt]
    
    print("\nАНАЛИТИЧЕСКОЕ РЕШЕНИЕ:")
    print(f"y({x_end:.1f}) = {y_analyt[-1]:.6f}")
    
    print("\nЯВНЫЕ МЕТОДЫ:")
    x_euler, y_euler = explicit_rk(ode_function, x0, y0, x_end, h, EULER_METHOD)
    print(f"Метод Эйлера: y({x_end:.1f}) = {y_euler[-1]:.6f}, погрешность = {abs(y_euler[-1] - y_analyt[-1]):.6f}")
    
    x_ek, y_ek = explicit_rk(ode_function, x0, y0, x_end, h, EULER_METHOD_CAUCHY)
    print(f"Метод Эйлера-Коши: y({x_end:.1f}) = {y_ek[-1]:.6f}, погрешность = {abs(y_ek[-1] - y_analyt[-1]):.6f}")
    
    x_kl, y_kl = explicit_rk(ode_function, x0, y0, x_end, h, RC_4)
    print(f"Классический метод: y({x_end:.1f}) = {y_kl[-1]:.6f}, погрешность = {abs(y_kl[-1] - y_analyt[-1]):.6f}")
    
    print("\nНЕЯВНЫЕ МЕТОДЫ:")
    x_rado, y_rado = implicit_rk(ode_function, x0, y0, x_end, h, RADO_IIA)
    print(f"Метод Радо IIA: y({x_end:.1f}) = {y_rado[-1]:.6f}, погрешность = {abs(y_rado[-1] - y_analyt[-1]):.6f}")
    
    x_pst, y_pst = implicit_rk(ode_function, x0, y0, x_end, h, MIDPOINT)
    print(f"Правило средней точки: y({x_end:.1f}) = {y_pst[-1]:.6f}, погрешность = {abs(y_pst[-1] - y_analyt[-1]):.6f}")
    
    x_g4, y_g4 = implicit_rk(ode_function, x0, y0, x_end, h, GAUSS_4_ORDER)
    print(f"Метод Гаусса 4 порядка: y({x_end:.1f}) = {y_g4[-1]:.6f}, погрешность = {abs(y_g4[-1] - y_analyt[-1]):.6f}")
    
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 1, 1)
    plt.plot(x_analyt, y_analyt, 'k-', label='Аналитическое', linewidth=2)
    plt.plot(x_euler, y_euler, 'o-', label='Метод Эйлера', markersize=3)
    plt.plot(x_ek, y_ek, 's-', label='Метод Эйлера-Коши', markersize=3)
    plt.plot(x_kl, y_kl, '^-', label='Классический метод', markersize=3)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Явные методы Рунге-Кутты')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(2, 1, 2)
    plt.plot(x_analyt, y_analyt, 'k-', label='Аналитическое', linewidth=2)
    plt.plot(x_rado, y_rado, 'o-', label='Метод Радо IIA', markersize=3)
    plt.plot(x_pst, y_pst, 's-', label='Правило средней точки', markersize=3)
    plt.plot(x_g4, y_g4, '^-', label='Метод Гаусса 4 порядка', markersize=3)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Неявные методы Рунге-Кутты')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()