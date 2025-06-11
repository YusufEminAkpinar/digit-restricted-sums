import sys
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from time import time

"""
for i <= 4
estimation at 20: 22.92067 66192 64150 34821 34152 00056 80054 51733 09591 222 

for i <= 3
estimation at 30: 22.92067 66192 64150 34816 35179 52007 83243 68295 76567 010

for i <= 4
estimation at 30: 22.92067 66192 64150 34816 36570 94377 65471 98756 35155 534
22.920676619264150348163657094377654719875635155534
22.920676619264150348163657094377666872811876482985

22.920676619264150348163657094375931914944822931241


22.9206766192641503481636570943759319149447624369984815685419983565721563381899111294456260374482018989909
"""

R = RealField(700)

# calculates (-1)^n * (j+n-1) choose n
def c_jn(j, n): 
    return ((-1)**n) * binomial(j+n-1, n)

# 1^n + 2^n + ... + 8^n
def b_n(n):
    if not n:
        return 9
    return sum(k**n for k in range(1, 9))

def a_jn(j, n):
    return c_jn(j, n) * b_n(n) * 10**(-j-n)

def is_admissible(n):
    return '9' not in str(n)

# s_ij = sum of a_jn * s(i-1, j+n) for n = 0 to inf
@lru_cache(maxsize=None)
def s_ij(i, j):
    # print(f"Calculating s({i}, {j})")
    retval = R(0.0)

    if i <= 4:
        return sum( R(1/k**j) for k in range(10**(i-1), 10**i) if is_admissible(k))

    limit = 20
    for n in range(0, limit):
        retval += a_jn(j, n) * s_ij(i-1, j+n)
    return retval


def estimateSum(n):
    retval = R(0.0)
    total_time = 0
    for i in range(1, n):
        s = time()
        print("Trying to sum for i =", i)
        sij = s_ij(i, 1)
        retval += sij
        iteration_time = time()-s
        print(f"Found, s_{i},{1} = {numerical_approx(sij)}, took {iteration_time} second.")
        total_time += iteration_time
    s = time()
    print("Trying to sum for i =", n)
    sij = s_ij(n, 1)
    retval += 10*sij
    total_time += (time()-s)
    print(f"Totally, it took {total_time} seconds.")
    return retval


argc = len(sys.argv)
if argc < 2:
    n = 3
else:
    n = int(sys.argv[1])
if argc < 3:
    k = 3
else:
    k = int(sys.argv[2])

value = estimateSum(n)
print(f"estimation at {n}: {value.n(digits=50)}")

