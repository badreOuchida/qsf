import math
from typing import List, Tuple, Dict
import numpy as np
from itertools import chain



def is_quadratic_residue(a, p) :
    """
    Vérifie si a est un résidu quadratique modulo p utilisant le symbole de Legendre
    """
    if p == 2:
        return True

    # Calcul du symbole de Legendre
    return pow(a, (p - 1) // 2, p) == 1

def generate_factor_base(n ,factor_base_size) :
    """
    Génère la base de facteurs premiers p tels que n est un résidu quadratique modulo p
    """
    factor_base = [2]  # 2 est toujours dans la base
    p = 3

    while len(factor_base) < factor_base_size:
        if is_quadratic_residue(n, p):
            factor_base.append(p)
        p = nextprime(p)

    return factor_base

def calculate_matrix_length(n):
    """
    Calculate appropriate matrix length based on the size of n.
    This uses empirical formulas based on number field sieve complexity.

    bits = n.bit_length()

    # For smaller numbers (less than 100 bits), use a smaller matrix
    if bits < 100:
        return max(100, bits * 3)

    # For medium size numbers
    if bits < 200:
        return int(math.exp(math.sqrt(bits * math.log(bits))) / 2)

    # For larger numbers, use the L(1/2) formula from number field sieve
    # L(1/2, c) = exp((c + o(1))(log n)^(1/2)(log log n)^(1/2))
    # where c is typically around 1"""
    c = 1.0
    log_n = math.log(n)
    log_log_n = math.log(log_n)
    return int(math.exp(c * math.sqrt(log_n * log_log_n)))

def sieve(n, factor_base) :
    """
    Applique le crible sur un intervalle pour trouver des relations
    """
    length_matrix = int(calculate_matrix_length(n))

    # Initialisation du tableau de criblage
    sieve_array = np.zeros((length_matrix, len(factor_base)), dtype=int)
    values = np.zeros(length_matrix, dtype=object)

    # Calcul des valeurs de la fonction polynomiale
    M = int(math.sqrt(n)) + 1
    for i in range(length_matrix):
        x = M + i
        values[i] = x * x - n



    # Criblage avec la base de facteurs
    for i, p in enumerate(factor_base):
        # Trouve le premier x tel que x^2 ≡ n (mod p)
        x = tonelli_shanks(n, p)
        if x is None:
            continue

        # Crible avec p
        # xj = M + j and xj = root mod p
        # and so root mod p = M + j and so j = root - M mod p
        for start in [(x - M) % p, (-x - M) % p]:
            for j in range(start, length_matrix, p):
                tmp = values[j]
                while tmp % p == 0:
                    tmp //= p
                    sieve_array[j][i] += 1
                values[j] = tmp


    # Collecte des relations trouvées
    relations = []
    for i in range(length_matrix):
        if values[i] == 1:  # Nombre complètement factorisé
            relations.append((M + i, list(sieve_array[i])))

    return relations


def tonelli_shanks(n, p) :

    """
    Implémentation de l'algorithme de Tonelli-Shanks pour trouver la racine carrée modulo p
    """

    if p == 2:
        return n % 2

    if not is_quadratic_residue(n, p):
        return None

    # Factorise p-1 = q * 2^s
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1

    # Trouve un non-résidu quadratique
    z = 2
    while is_quadratic_residue(z, p):
        z += 1

    # Initialisation
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)

    # Boucle principale
    while t != 1:
        # Trouve le plus petit i tel que t^(2^i) ≡ 1 (mod p)
        i = 0
        temp = t
        while temp != 1 :
            temp = pow(temp , 2 , p)
            i += 1

            if i >= m:
                return None

        # Calcule b = c^(2^(m-i-1))
        b = pow(c, 1 << (m - i - 1), p)

        # Mise à jour des variables
        r = (r * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        m = i

    return r


def binary_gaussian_elimination(A):

    A = A % 2
    rows, cols = A.shape
    pivots = [-1] * rows

    for col in range(cols):
        for row in range(rows):
            if A[row, col] == 1:
                pivots[row] = col
                for r in range(row + 1, rows):
                    if A[r, col] == 1:
                        A[r] = (A[r] + A[row]) % 2
                break

    pivot_cols = set(p for p in pivots if p != -1)
    free_cols = [c for c in range(cols) if c not in pivot_cols]

    solutions = []
    if free_cols:
        for free_col in free_cols:
            solution = np.zeros(cols, dtype=int)
            solution[free_col] = 1
            for row, pivot_col in enumerate(pivots):
                if pivot_col != -1:
                    solution[pivot_col] = (np.dot(A[row], solution) % 2)
            solutions.append(solution)
    else:
        solutions.append(np.zeros(cols, dtype=int))

    return solutions







if __name__ == '__main__':

    n = 13 * 17

    factor_base_size = int(math.exp(math.sqrt(math.log(n) * math.log(math.log(n)))) ** 0.5)
    print(factor_base_size)

    factor_base = generate_factor_base(n ,factor_base_size)
    print(factor_base)

    #print(sieve(n , factor_base))
    relations = sieve(n , factor_base)

    #print(relations)


    matrix = np.array([rel[1] for rel in relations]) % 2
    #print(matrix)
    matrix = np.transpose(matrix)
    #print(matrix)
    kernel = binary_gaussian_elimination(matrix)
    #print(kernel)
    #print(relations)




    for solution in kernel :
      x = 1
      y = 1
      for index , xi in enumerate(solution) :
        if xi != 0 :
          x *= relations[index][0]
        for j, exp in enumerate(relations[xi][1]):
          if xi != 0 :
            y *= int(math.pow(factor_base[j] , exp))


      y = int(math.sqrt(y))

      p = math.gcd(x + y, n)

      q = math.gcd(x - y, n)

      print(f"x : {x} and y : {y} and pcgd(x + y , n) : {p} and pgcd(x - y ,n) : {q}")








