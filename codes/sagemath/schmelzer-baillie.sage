from functools import lru_cache
from time import time
from bisect import bisect_left, bisect_right
import sys
sys.setrecursionlimit(2000)


K = 200
limit = K
epsilon = 10**(-(max(limit, K)))
prec = K
R = RealField(prec*4)

restricted = 314
# restricted = 9
# restricted = 42
# restricted = 31415
# restricted = 7

strnum = str(restricted)
digits = [int(d) for d in strnum]
digit_num = len(digits)
T = Matrix([ [1 for _ in range(10)] for _ in range(digit_num) ])


# S has size len(num)
def construct_S(num=restricted):
    S = []
    # We are dealing with i=0 separately
    tmp = []
    for i in range(1, 1001):
        conditions = []
        for j in range(digit_num):
            conditions.append( i%(10**(j+1)) != int(strnum[0:j+1]) )
        conditions.append(strnum not in str(i))
        if all(conditions):
            tmp.append(i)

    S.append(tmp)

    for j in range(1, digit_num):
        tmp = []
        for i in range(1, 1001):
            if i%(10**(j)) == int(strnum[0:j]) and strnum not in str(i):
                tmp.append(i)
        S.append(tmp)

    return S

def construct_matrix(num=restricted):
    T = []
    S = construct_S(num)
    for i in range(digit_num):
        L = [1 for _ in range(10)]
        for k in range(10):
            cur_digit = k
            new_num = S[i][0]*10 + cur_digit

            value = 1
            for j in range(1, digit_num):
                if new_num%(10**(j)) == int(strnum[0:j]):
                    value = j+1
                    break
                if strnum in str(new_num):
                    value = 0
                    break

            L[cur_digit] = value
        T.append(L)
    return Matrix(T)

S = construct_S()
T = construct_matrix()

# Implement with epsilon logic. Now it's fast...
# Try to use JIT
print(T)

# 75 chiffre  => 10 sec
# 150 chiffre => 30 sec
# 300 chiffre => 93 sec

# 2 -> 3
# 4 -> 3^2
# 8 -> 3^3


# 600 chiffre approx 270 sec
# 600 chiffre = 364.42 sec

# 100 chiffre => 20-30 sec
# 20-30 fois 
# 1000 chiffre 10-15 minute
# 1_000_000 chiffre => 30 * 30 * 30 * 30 minute

def c(k, w):
    return (1 - 2*(w % 2)) * binomial(k+w-1, w)


def a(k, w, m):
    if m+w == 0:
        return 10^(-k-w) * c(k, w)
    return 10^(-k-w) * c(k, w) * m^w


def f(j, l, m):
    return int(j == T[l-1, m])


@lru_cache(maxsize=None)
def psi(i, j, k):
    if (i <= 3):
        lower = 10**(i-1)-1
        left = bisect_right(S[j-1], lower)
        upper = 10**(i)
        right = bisect_left(S[j-1], upper)
        return sum( (R(1/s).n(digits=prec))**k for s in S[j-1][left:right])
    value = R(0)
    for m in range(10):
        for l in range(1, digit_num+1):
            skip = False
            ft = f(j, l, m)
            for w in range(limit+1):
                p = psi(i-1, l, k+w)
                value += ft * a(k, w, m) * p
                if p < epsilon:
                    skip = True
                    break
            if skip:
                continue
    print(f"psi({i}, {j}, {k}) = {value.n()}")
    return value


def B():
    An = Matrix(QQ, digit_num, digit_num)
    for m in range(10):
        for row in range(digit_num):
            for col in range(digit_num):
                An[row, col] += f(row+1, col+1, m)
    An /= 10
    eigens = An.eigenvalues()
    if (max(abs(e) for e in eigens) >= 1):
        print("Eigenvalues does not fall inside of the unit disc.", file=sys.stderr)
        sys.exit(-1)
    print("An is:")
    print(An)
    idm = identity_matrix(digit_num)
    return (idm - An).inverse() - idm

def finalResult():
    vect = vector([psi(K, j, 1) for j in range(1, digit_num+1)])
    print("Vector of psi's is computed. It is:")
    print(vect.n())
    approxSum = (B() * vect).norm(1)
    value = sum( 
              sum( R(psi(i, j, 1))
                  for i in range(1, K+1)) 
              for j in range(1, digit_num+1)) + approxSum
    print(f"Approximated part is: {approxSum.n()}")
    return value.n(digits=prec)

print(finalResult())
