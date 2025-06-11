from functools import lru_cache 
from time import time
from multiprocessing import Pool, cpu_count
import sys
from pprint import pprint
from math import ceil, floor
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyval
import psutil
import os


print("Usage:\n\
        sage burnol.sage <precision | int> <verbose | 0 or 1>")

start_t = time()

process = psutil.Process(os.getpid())

verbose = False
if (sys.argv[2] == '1'):
    verbose = True

l = 3

b = 10
N = 9
A = [0, 1, 2, 3, 4, 5, 6, 7, 8]
restricted_digits = [9]

prec = int(sys.argv[1])
t_first = time()
limit = int(prec/(l-1))+1

beta_limit = 10**(-prec*1.1)

bit_guard = 10
bits = ceil((prec+1)*log(10,2)) + bit_guard

RR = RealField(bits)
Rs = [RR]

result = R(0)

if verbose:
    print(f"l is {l} and precision is {prec}.")

n = int(prec/2000)
est_mem = 250 + 50*n*(n+1)

total_mem = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  
mem_mib = total_mem/(1024.**2)
if verbose:
    # print(f"Estimated time is: {est_time}")
    print(f"Estimated memory usage is: {est_mem}MB, and total memory is {mem_mib}MB")



binomMax = 1000
binomPerLine = 1000
if verbose:
    print("Pré-calcul des coefficients du binôme jusque n =", binomPerLine)

combi = [[1],[1, 1]]
for j in range(2, binomPerLine):
    # print(f"Calculating binomial for j={j}.")
    L = [1]
    for i in range(1, ceil(j/2)+1):
        cur_comb =  combi[j-1][ min(i, j-i) ] + combi[j-1][ min(i-1, j-i-1) ] 
        L.append(cur_comb)

    combi.append(L)

if verbose:
    print("Pré-calcul des coefficients du binôme jusque n =", binomPerLine)
def get_comb(x, y):
    global binomMax
    if (x < binomMax):
        return combi[x%1000][min(y, x-y)]

    if verbose:
        print(f"At least {binomPerLine} more line of binomial coefficients are needed.")
    for j in range(binomMax, binomMax+binomPerLine):
        L = [1]
        for i in range(1, ceil(j/2)+1):
            cur_comb = combi[(j-1)%1000][ min(i, j-i) ] + combi[(j-1)%1000][ min(i-1, j-i-1) ] 
            L.append(cur_comb)
        combi[j%1000] = L
    binomMax += binomPerLine
    return combi[x%1000][min(y, x-y)]


if verbose:
    print("Calculation finished.")

gamma_list = [(sum(i ** j for i in A)) for j in range(limit + 1)]


if verbose:
    print("u_m's are calculating.")
R = RR
u = [10]
append_to_u = u.append
t = time()
for m in range(1, limit+1):
    append_to_u( sum( get_comb(m, j)*gamma_list[j]*u[m-j] for j in range(1, m+1) ) / (b**(m+1)-R(N)) )
    if m%100 == 0:
        bits -= 100
        R = RealField(bits)
    if verbose and (m%100 == 0):
        print(f"u_{m} is calculated in {time() - t} seconds. With precision {int(bits)}.")
        t = time()
print(f"Total execution took {time() - start_t} seconds. Calculated {limit} amount of u.")
# sys.exit(0)


def lam(m):
    return m * R(9/8)**(m-1) * u[m-1]/u[0]

for i in range(limit+1):
    print(i, ": ", lam(i))

sys.exit(0)


if verbose:
    print("Calculation finished.")

def is_admissible(x):
    return '9' not in str(x)

def beta(l, m):
    s = 0
    if verbose:
        if m%100 == 0:
            print(f"beta({l}, {m}) computed.")
    for n in range(b^(l-1), b^l):
        if is_admissible(n):
            term = (1/(n**m)).n(digits=prec)
            if term < beta_limit:
                s += term
                break
            s += term
    return s

def compute_beta(m):
    return (m, beta(l, m+1))

def sommeFinie():
    s1 = RR(0)
    for n in range(1, b**(l-1)):
        if is_admissible(n):
            s1 += (1/n).n(digits=prec)
    if verbose:
        print("s1 finished.")
    s2 = RR(0)
    for n in range(b**(l-1), b**l):
        if is_admissible(n):
            s2 += (10/n).n(digits=prec)
    if verbose:
        print("SommeFinie finished.")
    return (s1+s2)#.n(digits=prec)

# Check if first $limit term is same for alternate serie here.
def sommeInfinie(limit=limit):
    num_workers = cpu_count()
    t = time()
    with Pool(num_workers) as pool:
        beta_values = dict(pool.imap_unordered(compute_beta, range(1, limit+1)))
    print(f"Beta's took {time()-t} seconds total.")
    #               (-1)^m
    return sum((1 - 2 * (m % 2)) * RR(u[m]) * beta_values[m] for m in range(1, limit+1))

def K(limit=limit):
    return sommeFinie() + sommeInfinie(limit)


if __name__ == '__main__':
    # for i in range(0, limit, 100):
    #     beta_i = beta(l, i)
    #     print(f"beta_{i}\n{beta_i}\nu[{i}]\n{u[i]}")
    value = K()
    # print("\nOur final value is:\n\n", value)
    # with open(f"k_prec_{prec}.txt", "w") as file:
    #     file.write(str(value)+"\n")

    # value = sommeFinie()
    print(value)
    print(f"Memory usage: {process.memory_info().rss / int(1024**2):.2f} MB")
    print(f"It took {time()-t_first} seconds.")

