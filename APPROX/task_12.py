import matplotlib.pyplot as plt

xi = [-4.71, -3.63, -2.55, -1.47, -0.39, 0.69, 1.77, 2.85, 3.93, 5.01, 6.09]
yi = [2.9853, 1.7984, 1.0358, 0.6987, 0.7116, 1.3214, 1.5742, 1.6371, 1.3936, 1.0309, 0.9754]
x_star = -0.926

def power_sum(x_arr, power):
    return sum([x**power for x in x_arr])

def power_sum_y(x_arr, y_arr, power_x):
    return sum([y * (x**power_x) for x, y in zip(x_arr, y_arr)])

def solve_gaussian(A, B):
    n = len(B)
    for i in range(n):
        for k in range(i + 1, n):
            c = A[k][i] / A[i][i]
            for j in range(i, n):
                A[k][j] -= c * A[i][j]
            B[k] -= c * B[i]

    x = [0 for _ in range(n)]
    for i in range(n - 1, -1, -1):
        x[i] = B[i] / A[i][i]
        for k in range(i - 1, -1, -1):
            B[k] -= A[k][i] * x[i]
    return x

def evaluate_poly(coeffs, x):
    res = 0
    for i, a in enumerate(coeffs):
        res += a * (x ** i)
    return res

def least_squares(x_nodes, y_nodes, degree):
    m = degree
    A = [[0.0] * (m + 1) for _ in range(m + 1)]
    B = [0.0] * (m + 1)
    
    for row in range(m + 1):
        for col in range(m + 1):
            A[row][col] = power_sum(x_nodes, row + col)

        B[row] = power_sum_y(x_nodes, y_nodes, row)
        
    coeffs = solve_gaussian(A, B)
    return coeffs


degrees = [1, 2, 3]
results = {}

print("Метод наименьших квадратов")

for deg in degrees:
    coeffs = least_squares(xi, yi, deg)
    
    phi = sum([(evaluate_poly(coeffs, x) - y)**2 for x, y in zip(xi, yi)])
    
    val_star = evaluate_poly(coeffs, x_star)
    
    results[deg] = {
        'coeffs': coeffs,
        'phi': phi,
        'val_star': val_star
    }
    
    print(f"\nМногочлен {deg}-й степени:")
    print(" ")
    poly_str = f"F_{deg}(x) = "
    for i, a in enumerate(coeffs):
        if i == 0:
            poly_str += f"{a:.4f}"
        else:
            sign = "+" if a >= 0 else "-"
            poly_str += f" {sign} {abs(a):.4f}x^{i}"
    print(poly_str)
    
    print("Коэффициенты (a0, a1...):")
    print(f"{[round(c, 4) for c in coeffs]}")
    print(f"Сумма квадратов ошибок Φ_{deg} = {phi:.4f}")
    print(f"Значение F_{deg}(x*) = {val_star:.4f}")


x_plot = [min(xi) + i*(max(xi)-min(xi))/100 for i in range(101)]

plt.figure(figsize=(10, 7))

plt.scatter(xi, yi, color='black', label='Исходные данные (yi)', zorder=5)

colors = ['green', 'blue', 'red']
styles = ['--', '-.', '-']

for deg in degrees:
    y_plot = [evaluate_poly(results[deg]['coeffs'], x) for x in x_plot]
    plt.plot(x_plot, y_plot, color=colors[deg-1], linestyle=styles[deg-1], 
             label=f'F_{deg}(x) (Степень {deg})')
    
    y_star_plot = results[deg]['val_star']
    plt.scatter([x_star], [y_star_plot], color=colors[deg-1], s=50, marker='x')

plt.axvline(x=x_star, color='gray', linestyle=':', alpha=0.5, label=f'x*={x_star}')
plt.title("Аппроксимация функций методом МНК")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.show()
