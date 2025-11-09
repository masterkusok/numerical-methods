import matplotlib.pyplot as plt
import math

def phi2(x1, x2):
    denom = 2*x1 + 3*x2
    if abs(denom) < 1e-10:
        return None
    return (3*x1**2 - 5*x1 + 3) / denom

# Якобиан для метода простой итерации
def jacobian_simple(x1, x2):
    val = 8 - x2**2
    if val <= 0:
        return None

    dphi1_dx1 = 0
    dphi1_dx2 = -x2 / math.sqrt(val)

    denom = 2*x1 + 3*x2
    if abs(denom) < 1e-10:
        denom = 1e-10

    num = 3*x1**2 - 5*x1 + 3
    dphi2_dx1 = ((6*x1 - 5)*(2*x1 + 3*x2) - 2*num) / (denom**2)
    dphi2_dx2 = -3*num / (denom**2)

    return [[dphi1_dx1, dphi1_dx2],
            [dphi2_dx1, dphi2_dx2]]

# Якобиан для метода Зейделя
def jacobian_seidel(x1, x2):
    val = 8 - x2**2
    if val <= 0:
        return None

    # φ1 такая же
    dphi1_dx1 = 0
    dphi1_dx2 = -x2 / math.sqrt(val)

    # Для Зейделя φ2 зависит от φ1(x1, x2)
    phi1_val = phi1(x1, x2)
    if phi1_val is None:
        return None

    denom = 2*phi1_val + 3*x2
    if abs(denom) < 1e-10:
        denom = 1e-10

    num = 3*phi1_val**2 - 5*phi1_val + 3
    dphi2_dphi1 = ((6*phi1_val - 5)*(2*phi1_val + 3*x2) - 2*num) / (denom**2)
    dphi2_dx2 = -3*num / (denom**2)

    # dφ2/dx1 = dφ2/dφ1 * dφ1/dx1  (а dφ1/dx1 = 0)
    dphi2_dx1 = dphi2_dphi1 * 0
    # dφ2/dx2 = dφ2/dφ1 * dφ1/dx2 + ∂φ2/∂x2
    dphi2_dx2 = dphi2_dphi1 * dphi1_dx2 + dphi2_dx2

    return [[dphi1_dx1, dphi1_dx2],
            [dphi2_dx1, dphi2_dx2]]

# Общая функция проверки
def check_convergence(method, x1, x2):
    if method == "simple":
        J = jacobian_simple(x1, x2)
        title = "Метод простой итерации"
    elif method == "seidel":
        J = jacobian_seidel(x1, x2)
        title = "Метод Зейделя"
    else:
        raise ValueError("method должен быть 'simple' или 'seidel'")

    if J is None:
        print(f"{title}: невозможно вычислить Якобиан в точке ({x1}, {x2})")
        return

    norm_inf = max(abs(J[0][0]) + abs(J[0][1]),
                   abs(J[1][0]) + abs(J[1][1]))

    print(f"\n{title} — проверка сходимости в точке ({x1:.2f}, {x2:.2f})")
    print("Якобиан Φ(x):")
    for row in J:
        print(f"  {row}")
    print(f"‖JΦ(x)‖∞ = {norm_inf:.4f}")

    if norm_inf < 1:
        print("✅ Условие сходимости выполнено\n")
    else:
        print("❌ Условие сходимости не выполнено\n")


# Система уравнений:
# x1^2 + x2^2 = 8
# 3x1^2 - 2x1x2 - 3x2^2 - 5x1 + 3 = 0

def f1(x1, x2):
    return x1**2 + x2**2 - 8

def f2(x1, x2):
    return 3*x1**2 - 2*x1*x2 - 3*x2**2 - 5*x1 + 3

# Эквивалентные функции для графиков
def x2_from_f1_pos(x1):
    val = 8 - x1**2
    return val**0.5 if val >= 0 else None

def x2_from_f1_neg(x1):
    val = 8 - x1**2
    return -val**0.5 if val >= 0 else None

def x2_from_f2(x1):
    # 3x1^2 - 2x1x2 - 3x2^2 - 5x1 + 3 = 0
    # -3x2^2 - 2x1x2 + 3x1^2 - 5x1 + 3 = 0
    # 3x2^2 + 2x1x2 - 3x1^2 + 5x1 - 3 = 0
    a = 3
    b = 2*x1
    c = -3*x1**2 + 5*x1 - 3
    d = b**2 - 4*a*c
    if d < 0:
        return None, None
    return (-b + d**0.5)/(2*a), (-b - d**0.5)/(2*a)

# Построение графиков
x1_vals = [i*0.01 for i in range(-300, 301)]
x2_pos = [x2_from_f1_pos(x1) for x1 in x1_vals]
x2_neg = [x2_from_f1_neg(x1) for x1 in x1_vals]
x2_f2_pos = []
x2_f2_neg = []
for x1 in x1_vals:
    y1, y2 = x2_from_f2(x1)
    x2_f2_pos.append(y1)
    x2_f2_neg.append(y2)

plt.figure(figsize=(10, 8))
plt.plot(x1_vals, x2_pos, 'b-', label='x1² + x2² = 8 (верх)')
plt.plot(x1_vals, x2_neg, 'b--', label='x1² + x2² = 8 (низ)')
plt.plot(x1_vals, x2_f2_pos, 'r-', label='f2 = 0 (ветвь 1)')
plt.plot(x1_vals, x2_f2_neg, 'r--', label='f2 = 0 (ветвь 2)')
plt.grid(True)
plt.xlabel('x1')
plt.ylabel('x2')
plt.legend()
plt.title('Графики эквивалентных функций')
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.show()

# Начальные приближения (из графика)
initial_guesses = [(2.5, 1.0), (-0.5, 2.8), (-0.5, -2.8), (2.5, -1.0)]

eps = 0.0001

# Метод Ньютона
def newton_method(x1_0, x2_0, eps):
    x1, x2 = x1_0, x2_0
    iterations = 0
    print(f"\nМетод Ньютона, начальное приближение: ({x1_0}, {x2_0})")
    
    while True:
        iterations += 1
        f1_val = f1(x1, x2)
        f2_val = f2(x1, x2)
        
        # Якобиан
        df1_dx1 = 2*x1
        df1_dx2 = 2*x2
        df2_dx1 = 6*x1 - 2*x2 - 5
        df2_dx2 = -2*x1 - 6*x2
        
        det = df1_dx1*df2_dx2 - df1_dx2*df2_dx1
        if abs(det) < 1e-10:
            print("Определитель близок к нулю")
            return None, None, iterations
        
        dx1 = (f1_val*df2_dx2 - f2_val*df1_dx2) / det
        dx2 = (f2_val*df1_dx1 - f1_val*df2_dx1) / det
        
        x1_new = x1 - dx1
        x2_new = x2 - dx2
        
        if abs(x1_new - x1) < eps and abs(x2_new - x2) < eps:
            x1, x2 = x1_new, x2_new
            break
        
        x1, x2 = x1_new, x2_new
        
        if iterations > 1000:
            print("Превышено максимальное число итераций")
            return None, None, iterations
    
    print(f"Решение: x1 = {x1:.4f}, x2 = {x2:.4f}")
    print(f"Проверка: f1 = {f1(x1, x2):.6f}, f2 = {f2(x1, x2):.6f}")
    print(f"Итераций: {iterations}")
    return x1, x2, iterations

# Метод простой итерации
def simple_iteration(x1_0, x2_0, eps):
    print(f"\nМетод простой итерации, начальное приближение: ({x1_0}, {x2_0})")
    # Используем релаксацию: x_new = x_old + tau * (phi(x) - x_old)
    # x1 = x1 + tau * (sqrt(8 - x2^2) - x1)
    # x2 = x2 + tau * ((3x1^2 - 5x1 + 3) / (2x1 + 3x2) - x2)
    
    print("Проверка сходимости:")
    check_convergence("simple", x1_0, x2_0)

    x1, x2 = x1_0, x2_0
    iterations = 0
    tau = 0.3
    
    while True:
        iterations += 1
        
        val = 8 - x2**2
        if val < 0:
            val = 0
        phi1 = val**0.5 if x1 > 0 else -val**0.5
        
        denom = 2*x1 + 3*x2
        if abs(denom) < 1e-10:
            denom = 1e-10
        phi2 = (3*x1**2 - 5*x1 + 3) / denom
        
        x1_new = x1 + tau * (phi1 - x1)
        x2_new = x2 + tau * (phi2 - x2)
        
        if abs(x1_new - x1) < eps and abs(x2_new - x2) < eps:
            x1, x2 = x1_new, x2_new
            break
        
        x1, x2 = x1_new, x2_new
        
        if iterations > 10000:
            print("Превышено максимальное число итераций")
            return None, None, iterations
    
    print(f"Решение: x1 = {x1:.4f}, x2 = {x2:.4f}")
    print(f"Проверка: f1 = {f1(x1, x2):.6f}, f2 = {f2(x1, x2):.6f}")
    print(f"Итераций: {iterations}")
    return x1, x2, iterations

# Метод Зейделя
def seidel_method(x1_0, x2_0, eps):
    print(f"\nМетод Зейделя, начальное приближение: ({x1_0}, {x2_0})")
    # Используем релаксацию с обновлением x1, затем сразу используем новое x1
    # для x2
    
    print("Проверка сходимости:")
    check_convergence("seidel", x1_0, x2_0)
    
    x1, x2 = x1_0, x2_0
    iterations = 0
    tau = 0.3
    
    while True:
        iterations += 1
        
        val = 8 - x2**2
        if val < 0:
            val = 0
        phi1 = val**0.5 if x1 > 0 else -val**0.5
        x1_new = x1 + tau * (phi1 - x1)
        
        denom = 2*x1_new + 3*x2
        if abs(denom) < 1e-10:
            denom = 1e-10
        phi2 = (3*x1_new**2 - 5*x1_new + 3) / denom
        x2_new = x2 + tau * (phi2 - x2)
        
        if abs(x1_new - x1) < eps and abs(x2_new - x2) < eps:
            x1, x2 = x1_new, x2_new
            break
        
        x1, x2 = x1_new, x2_new
        
        if iterations > 10000:
            print("Превышено максимальное число итераций")
            return None, None, iterations
    
    print(f"Решение: x1 = {x1:.4f}, x2 = {x2:.4f}")
    print(f"Проверка: f1 = {f1(x1, x2):.6f}, f2 = {f2(x1, x2):.6f}")
    print(f"Итераций: {iterations}")
    return x1, x2, iterations

# Решение для всех начальных приближений
print("="*60)
print("РЕШЕНИЕ СИСТЕМЫ НЕЛИНЕЙНЫХ УРАВНЕНИЙ")
print("="*60)

for guess in initial_guesses:
    print("\n" + "="*60)
    print(f"Начальное приближение: ({guess[0]}, {guess[1]})")
    print("="*60)
    newton_method(guess[0], guess[1], eps)
    simple_iteration(guess[0], guess[1], eps)
    seidel_method(guess[0], guess[1], eps)
