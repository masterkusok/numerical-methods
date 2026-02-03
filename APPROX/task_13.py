import math
import matplotlib.pyplot as plt

def f(x):
    """Исходная функция"""
    return math.log(2*x + 3)/math.sqrt(x + 2) + (x + 2)**x/(x**2 + 1)

def df(x):
    """Первая производная"""
    return (-(2*x*(x + 2)**x)/(x**2 + 1)**2 + 
            ((x + 2)**x * (x/(x + 2) + math.log(x + 2)))/(x**2 + 1) + 
            2/((2*x + 3) * math.sqrt(x + 2)) - 
            math.log(2*x + 3)/(2 * (x + 2)**(3/2)))

def df2(x):
    """Вторая производная"""
    return (((8*x**2)/(x**2 + 1)**3 - 2/(x**2 + 1)**2) * (x + 2)**x - 
            (4*x*(x + 2)**x * (x/(x + 2) + math.log(x + 2)))/(x**2 + 1)**2 + 
            ((2/(x + 2) - x/(x + 2)**2) * (x + 2)**x + (x + 2)**x * (x/(x + 2) + math.log(x + 2))**2)/(x**2 + 1) - 
            4/((2*x + 3)**2 * math.sqrt(x + 2)) - 
            2/((2*x + 3) * (x + 2)**(3/2)) + 
            (3 * math.log(2*x + 3))/(4 * (x + 2)**(5/2)))


def first_deriv_right_diff(fn, x, h, a=-1.0, b=2.0):
    if x + h <= b + 1e-10:
        return (fn(x + h) - fn(x)) / h
    return float('nan')


def first_deriv_left_diff(fn, x, h, a=-1.0, b=2.0):
    if x - h >= a - 1e-10:
        return (fn(x) - fn(x - h)) / h
    return float('nan')


def first_deriv_central(fn, x, h, a=-1.0, b=2.0):
    if x - h >= a - 1e-10 and x + h <= b + 1e-10:
        return (fn(x + h) - fn(x - h)) / (2 * h)
    return float('nan')


def first_deriv_three_point(fn, x, h, a=-1.0, b=2.0):
    if x + 2*h <= b + 1e-10:
        return (-3 * fn(x) + 4 * fn(x + h) - fn(x + 2 * h)) / (2 * h)
    return float('nan')


def first_deriv_four_point(fn, x, h, a=-1.0, b=2.0):
    if (x - h >= a - 1e-10) and (x + 2*h <= b + 1e-10):
        return (-2 * fn(x - h) - 3 * fn(x) + 6 * fn(x + h) - fn(x + 2 * h)) / (6 * h)
    return float('nan')


def first_deriv_five_point(fn, x, h, a=-1.0, b=2.0):
    if x - 2*h >= a - 1e-10 and x + 2*h <= b + 1e-10:
        return (fn(x - 2 * h) - 8 * fn(x - h) + 8 * fn(x + h) - fn(x + 2 * h)) / (12 * h)
    return float('nan')


# Схемы численного дифференцирования второго порядка
def second_deriv_three_point(fn, x, h, a=-1.0, b=2.0):
    if x - h >= a - 1e-10 and x + h <= b + 1e-10:
        return (fn(x - h) - 2 * fn(x) + fn(x + h)) / (h * h)
    return float('nan')


def second_deriv_four_point_forward(fn, x, h, a=-1.0, b=2.0):
    if x + 3*h <= b + 1e-10:
        return (2 * fn(x) - 5 * fn(x + h) + 4 * fn(x + 2 * h) - fn(x + 3 * h)) / (h * h)
    return float('nan')


def second_deriv_four_point_backward(fn, x, h, a=-1.0, b=2.0):
    if x - 3*h >= a - 1e-10:
        return (2 * fn(x) - 5 * fn(x - h) + 4 * fn(x - 2 * h) - fn(x - 3 * h)) / (h * h)
    return float('nan')


def second_deriv_five_point(fn, x, h, a=-1.0, b=2.0):
    if x - 2*h >= a - 1e-10 and x + 2*h <= b + 1e-10:
        return (-fn(x - 2 * h) + 16 * fn(x - h) - 30 * fn(x) + 16 * fn(x + h) - fn(x + 2 * h)) / (12 * h * h)
    return float('nan')

if __name__ == "__main__":
    h = 0.15
    h1 = h
    h2 = h / 2
    
    first_schemes = [
        ("Правая разность", first_deriv_right_diff, '-'),
        ("Левая разность", first_deriv_left_diff, '--'),
        ("Центральная", first_deriv_central, '-.'),
        ("3-точечная", first_deriv_three_point, ':'),
        ("4-точечная", first_deriv_four_point, (0, (3, 1, 1, 1))),
        ("5-точечная", first_deriv_five_point, (0, (5, 2, 1, 2)))
    ]
    
    second_schemes = [
        ("3-точечная", second_deriv_three_point, '-'),
        ("4-точечная вперед", second_deriv_four_point_forward, '--'),
        ("4-точечная назад", second_deriv_four_point_backward, '-.'),
        ("5-точечная", second_deriv_five_point, ':')
    ]
    
    x_range = [-1.0 + i * 0.05 for i in range(61)]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    ax = axes[0, 0]
    ax.plot(x_range, [df(x) for x in x_range], 'k-', linewidth=3, label='Аналитическая')
    for name, scheme, style in first_schemes:
        ax.plot(x_range, [scheme(f, x, h1) for x in x_range], linestyle=style, linewidth=1.5, label=name, color='gray')
    ax.set_title(f"Первая производная, h = {h1}")
    ax.set_xlim([-1.0, 2.0])
    ax.legend()
    ax.grid(True)
    
    ax = axes[0, 1]
    ax.plot(x_range, [df(x) for x in x_range], 'k-', linewidth=3, label='Аналитическая')
    for name, scheme, style in first_schemes:
        ax.plot(x_range, [scheme(f, x, h2) for x in x_range], linestyle=style, linewidth=1.5, label=name, color='gray')
    ax.set_title(f"Первая производная, h = {h2}")
    ax.set_xlim([-1.0, 2.0])
    ax.legend()
    ax.grid(True)
    
    ax = axes[1, 0]
    ax.plot(x_range, [df2(x) for x in x_range], 'k-', linewidth=3, label='Аналитическая')
    for name, scheme, style in second_schemes:
        ax.plot(x_range, [scheme(f, x, h1) for x in x_range], linestyle=style, linewidth=1.5, label=name, color='gray')
    ax.set_title(f"Вторая производная, h = {h1}")
    ax.set_xlim([-1.0, 2.0])
    ax.legend()
    ax.grid(True)
    
    ax = axes[1, 1]
    ax.plot(x_range, [df2(x) for x in x_range], 'k-', linewidth=3, label='Аналитическая')
    for name, scheme, style in second_schemes:
        ax.plot(x_range, [scheme(f, x, h2) for x in x_range], linestyle=style, linewidth=1.5, label=name, color='gray')
    ax.set_title(f"Вторая производная, h = {h2}")
    ax.set_xlim([-1.0, 2.0])
    ax.legend()
    ax.grid(True)
    
    plt.tight_layout()
    plt.show()
