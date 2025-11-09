def check_thomas_stability(a, b, c):
    n = len(b)
    
    diagonal_ok = all(abs(b[i]) >= abs(a[i]) + abs(c[i]) for i in range(n))
    
    boundary_ok = (abs(c[0]/b[0]) < 1) and (abs(a[n-1]/b[n-1]) < 1)
    
    print(f"Диагональное преобладание: {'✓' if diagonal_ok else '✗'}")
    print(f"Граничные условия: {'✓' if boundary_ok else '✗'}")
    print(f"Метод устойчив: {'✓' if diagonal_ok and boundary_ok else '✗'}")
    
    return diagonal_ok and boundary_ok

def thomas_with_logs(matrix, d, eps=1e-12, round_digits=10):
    n = len(matrix)
    a = [0.0] * n
    b = [0.0] * n
    c = [0.0] * n
    
    for i in range(n):
        b[i] = matrix[i][i]
        
        if i > 0:
            a[i] = matrix[i][i-1]
        
        if i < n-1:
            c[i] = matrix[i][i+1]

    if not check_thomas_stability(a, b, c):
        print("Завершение...")
        exit()

    cp = [0.0]*n
    dp = [0.0]*n
    den = [0.0]*n

    den[0] = b[0]

    cp[0] = c[0] / den[0]
    dp[0] = d[0] / den[0]
    print(f"i=0: den={den[0]:f}, cp={cp[0]:f}, dp={dp[0]:f}")

    for i in range(1, n):
        den[i] = b[i] - a[i] * cp[i-1]
        
        cp[i] = (c[i] / den[i]) if i != n-1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i-1]) / den[i]
        print(f"i={i}: den={den[i]:f}, cp={cp[i]:f}, dp={dp[i]:f}")

    det = 1.0
    for v in den:
        det *= v

    x = [0.0]*n
    x[-1] = dp[-1]
    print(f"x[{n-1}] = {x[-1]:f}")
    for i in range(n-2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i+1]
        print(f"x[{i}] = dp[{i}] - cp[{i}]*x[{i+1}] = {dp[i]:.{round_digits}f} - {cp[i]:.{round_digits}f}*{x[i+1]:.{round_digits}f} = {x[i]:.{round_digits}f}")

    Ax = [( (a[i]*x[i-1] if i>0 else 0) + b[i]*x[i] + (c[i]*x[i+1] if i<n-1 else 0) ) for i in range(n)]
    residual = [Ax[i] - d[i] for i in range(n)]
    max_err = max(abs(v) for v in residual)

    x_clean = [0.0 if abs(val) < 1e-10 else val for val in x]
    
    return {
        "den": den,
        "cp": cp,
        "dp": dp,
        "x": x_clean,  # используем очищенный вектор
        "det": det,
        "Ax": Ax,
        "residual": residual,
        "max_error": max_err
    }

A = [
    [8, 7, 0, 0, 0, 0, 0, 0],
    [7, 13, -4, 0, 0, 0, 0, 0],
    [0, 2, 11, 8, 0, 0, 0, 0],
    [0, 0, -4, 11, -7, 0, 0, 0],
    [0, 0, 0, 6, -13, 5, 0, 0],
    [0, 0, 0, 0, 5, 11, 4, 0],
    [0, 0, 0, 0, 0, 3, 9, 5],
    [0, 0, 0, 0, 0, 0, 4, 5],
    ]
d = [16.0, 26.0, 31.0, 51.0, -53.0, 33.0, 54.0, 35.0]

res = thomas_with_logs(A,d)
print("\nРешение x:", [v for v in res["x"]])
print("Определитель ≈", res["det"])
