# Cumulative class counts N(T) for n=3, from the enumeration in n3_psl.py,
# against the Eskin-Mirzakhani form e^{hR}/(hR) with R = log(lambda) and
# h = 2n-4 = 2, and against li(lambda^2), the prime geodesic theorem for
# the modular group.  Primitive classes are the words that are not proper
# powers.
import math, sys
import mpmath
sys.argv = ["psl.py", "150", "/dev/null"]
src = open(__file__.replace("n3_growth.py","n3_psl.py")).read()
exec(src.split("T = int(sys.argv[1])")[0])
cls = enumerate_classes(150)
def primitive(w):
    n = len(w)
    return all(w != w[:d]*(n//d) for d in range(1,n) if n % d == 0)
print("   T   lambda    N(T)  prim  lam^2/(2 log lam)  prim/that  li(lam^2)  prim/li")
for T in (10, 20, 40, 80, 150):
    lam = 0.5*(T + math.sqrt(T*T-4))
    N = sum(1 for t in cls.values() if t <= T)
    P = sum(1 for w,t in cls.items() if t <= T and primitive(w))
    em = lam**2/(2*math.log(lam))
    li = float(mpmath.li(lam**2))
    print("%4d %8.2f %7d %6d %17.1f %10.3f %10.1f %8.3f" % (T, lam, N, P, em, P/em, li, P/li))
