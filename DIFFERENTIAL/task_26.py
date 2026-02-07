import matplotlib.pyplot as plt


def system(x, y):
    """Система ОДУ: y'' = f(x, y, y')
    y[0] = y, y[1] = y'
    y'[0] = y[1]
    y'[1] = (2x² + y - xy') / x²
    """
    y0, y1 = y
    return [y1, (2*x**2 + y0 - x*y1) / x**2]


def analytical(x):
    """Аналитическое решение"""
    return 2/3 * x**2 + x + 1/x


def runge_kutta_4(x0, y0, s, x_end, h):
    """Метод Рунге-Кутты 4 порядка для системы ОДУ"""
    x_vals, y_vals = [x0], [[y0, s]]
    x, y = x0, [y0, s]
    
    while x < x_end - h/2:
        k1 = [h * fi for fi in system(x, y)]
        k2 = [h * fi for fi in system(x + h/2, [y[i] + k1[i]/2 for i in range(2)])]
        k3 = [h * fi for fi in system(x + h/2, [y[i] + k2[i]/2 for i in range(2)])]
        k4 = [h * fi for fi in system(x + h, [y[i] + k3[i] for i in range(2)])]
        
        y = [y[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6 for i in range(2)]
        x += h
        x_vals.append(x)
        y_vals.append(y[:])
    
    return x_vals, y_vals


def shooting_dichotomy(x0, y0, x_end, y_end, h, tol=1e-6):
    s_low, s_high = -100, 100
    history = []
    
    for i in range(100):
        s_mid = (s_low + s_high) / 2
        x_vals, y_vals = runge_kutta_4(x0, y0, s_mid, x_end, h)
        residual = y_vals[-1][0] - y_end
        
        if i % 3 == 0:
            history.append((x_vals, y_vals))
        
        if abs(residual) < tol:
            return s_mid, runge_kutta_4(x0, y0, s_mid, x_end, h), history
        
        _, y_low = runge_kutta_4(x0, y0, s_low, x_end, h)
        if (y_low[-1][0] - y_end) * residual < 0:
            s_high = s_mid
        else:
            s_low = s_mid
    
    return s_mid, runge_kutta_4(x0, y0, s_mid, x_end, h), history


def main():
    x0, y0 = 0.5, 2.666
    x_end, y_end = 3.5, 11.952
    h = 0.15
    
    print("РЕШЕНИЕ КРАЕВОЙ ЗАДАЧИ МЕТОДОМ СТРЕЛЬБЫ")
    print(f"Уравнение: x²y'' + xy' - y = 2x²")
    print(f"Граничные условия: y({x0}) = {y0}, y({x_end}) = {y_end}")
    print(f"Шаг: h = {h}\n")
    
    s, (x_vals, y_vals), history = shooting_dichotomy(x0, y0, x_end, y_end, h)
    print(f"Метод дихотомии: y'({x0}) = {s:.6f}")
    print(f"Количество сохраненных траекторий: {len(history)}\n")
    
    x_analyt = [x0 + i*h for i in range(int((x_end - x0)/h) + 1)]
    y_analyt = [analytical(x) for x in x_analyt]
    y_numeric = [y[0] for y in y_vals]
    
    errors = [abs(y_numeric[i] - y_analyt[i]) for i in range(len(x_analyt))]
    print(f"Максимальная погрешность: {max(errors):.6e}")
    
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    linestyles = ['--', '-.', ':']
    for idx, (x_h, y_h) in enumerate(history):
        y_h_vals = [y[0] for y in y_h]
        plt.plot(x_h, y_h_vals, linestyles[idx % 3], color='gray', linewidth=1)
    plt.plot([], [], '--', color='gray', label='Промежуточные попытки')
    plt.plot(x_analyt, y_analyt, 'k-', label='Аналитическое', linewidth=2)
    plt.plot(x_vals, y_numeric, 'o-', label='Численное решение', markersize=4)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Решение краевой задачи')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    from matplotlib.ticker import LogLocator, NullFormatter
    ax = plt.gca()
    plt.plot(x_analyt, errors, 'o-', markersize=4)
    plt.xlabel('x')
    plt.ylabel('Погрешность')
    plt.title('Погрешность численного решения')
    plt.yscale('log')
    ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs='all'))
    plt.grid(True, which='both', alpha=0.3)
    
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
