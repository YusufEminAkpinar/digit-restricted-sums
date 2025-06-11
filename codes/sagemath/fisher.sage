from  sys import argv
from functools import lru_cache
from time import time

prec = 150

R = RealField(prec * 4)

def h(k):    
    return R(1)+sum(R(1)/s^k for s in range(2, 9))

alphas = [10, 360/91]
# Computing alpha(n) in definition (6) by taking the sum from 1 to n+1
def alpha(n):
    if n < 2:
        return alphas[n]
    print(f"Computing alpha({n})")
    numerator = R(10) - sum(binomial(n+1, k)*(10**(n-k+2) - 9**k)*alphas[n+1-k] for k in range(2, n+2))
    denom = (n+1)*(10^(n+1) - 9)
    alphas.append(numerator/denom)
    return numerator/denom

Poly.<x> = PolynomialRing(QQ)

def T(n):
    f = 1/2 * ( (x+sqrt(x^2-1))^n + (x-sqrt(x^2 - 1))^n)
    return f.expand().polynomial(Poly)

def Q(n):
    numerator = (T(n+1)(x=3) - T(n+1)(x=1-2*x))
    return numerator.quo_rem(1+x)[0]

def q(n):
    func, _ = Q(n).quo_rem(T(n+1)(x=3))
    return [ R(c[0]) for c in func.coefficients() ]

def test(n):
    smallCoef = R(1) - R(1)/(T(n+1).subs(x=3))
    bigCoef = R(1) + R(1)/(T(n+1).subs(x=3))
    middleValue = sum( alpha(k) * q(n)[k] * h(k+1) for k in range(n+1))
    lowerBound = middleValue/smallCoef
    upperBound = middleValue/bigCoef
    print(f"Lower Bound is:\n{lowerBound.n(digits=prec)}")
    print(f"Upper Bound is:\n{upperBound.n(digits=prec)}")
    # print(smallCoef.n())
    # print(middleValue.n())
    # print(bigCoef.n())

test(prec)

