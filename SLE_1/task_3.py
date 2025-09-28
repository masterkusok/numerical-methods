def thomas_full(a, b, c, d, eps=1e-12, round_digits=3):
    n = len(b)
    cp = [0.0]*n
    dp = [0.0]*n
    den = [0.0]*n

    den[0] = b[0]
    cp[0] = c[0] / den[0]
    dp[0] = d[0] / den[0]

    for i in range(1, n):
        den[i] = b[i] - a[i] * cp[i-1]
        if abs(den[i]) < eps:
            raise ZeroDivisionError(f"Zero (or tiny) pivot at index {i} -> need pivoting")
        cp[i] = c[i] / den[i] if i != n-1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i-1]) / den[i]

    x = [0.0]*n
    x[-1] = dp[-1]
    for i in range(n-2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i+1]

    det = 1.0
    for val in den:
        det *= val

    res = [0.0]*n
    for i in range(n):
        s = b[i]*x[i]
        if i>0: s += a[i]*x[i-1]
        if i<n-1: s += c[i]*x[i+1]
        res[i] = s - d[i]

    return {
        "x": x,
        "det": det,
        "residual": res,
        "x_rounded": [round(v, round_digits) for v in x],
        "det_rounded": round(det, round_digits)
    }

a = [0.0, 7.0, 2.0, -4.0, 6.0, 5.0, 3.0, 4.0]
b = [8.0, 13.0, 11.0, 11.0, -13.0, 11.0, 9.0, 5.0]
c = [7.0, -4.0, 8.0, -7.0, 5.0, 4.0, 5.0, 0.0]
d = [16.0, 26.0, 31.0, 51.0, -53.0, 33.0, 54.0, 35.0]

res = thomas_full(a,b,c,d, round_digits=3)
print("x (rounded):", res["x_rounded"])
print("determinant (rounded):", res["det_rounded"])
print("residual (max abs):", max(abs(v) for v in res["residual"]))
