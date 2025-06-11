import sys
from functools import lru_cache
from time import time

prec = 1000
R = RealField(prec * 4)

@lru_cache(maxsize=None)
# Computing beta(n) in definition (7) by taking the sum from 1 to n+1
def beta(n):
    if n == 0:
        return 10
    if n == 1:
        return 550/91
    # print(f"Calculating beta({n})")
    numerator = R(10) * (11^(n+1) - 10^(n+1)) - sum(binomial(n+1, k)*(10^(n-k+2) - 10^k + 1)*beta(n+1-k) for k in range(2, n+2))
    denom = (n+1) * (10^(n+1) - 9)
    return (numerator/denom).n(digits=prec)


def sumUpTo(upTo: int):
    s = 0
    for k in range(2, upTo):
        if not (k%100):
            print("Trying", k)
        s += 10^(-k)*beta(k-1)*zeta(k).n(digits=prec)
    retval = beta(0) * ln(10) - s
    print(f"Until {upTo}:\n{retval.n(digits=prec)}")
    retval += 10^(-upTo)*beta(upTo-1)*zeta(upTo)
    return retval
    # return beta(0) * ln(10) - sum([10^(-k)*beta(k-1)*zeta(k) for k in range(2, upTo)])

print(sumUpTo(prec).n(digits=prec))


