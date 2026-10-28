import numpy as np

A = np.array([
    [2.0, 1.0, 1.0],
    [4.0, 3.0, 3.0],
    [8.0, 7.0, 9.0]
], dtype=float)

B = np.array([5.0, 13.0, 33.0], dtype=float)




def decomposicao_lu(A):
    n = A.shape[0]
    L = np.eye(n, dtype=float)
    U = A.copy().astype(float)

    for k in range(n - 1):
        for i in range(k + 1, n):
            fator = U[i, k] / U[k, k]  # Multiplicador do elemento pivô
            L[i, k] = fator           # Armazena o fator na matriz L
            
            for j in range(k + 1, n):
                U[i, j] -= fator * U[k, j]  # Atualiza a linha da matriz U

            U[i, k] = 0.0             # Zera explicitamente abaixo do pivô

    return L, U


# SUBSTITUIÇÃO PROGRESSIVA (Forward Substitution)
# Resolve o sistema triangular inferior: L * y = b
# ===================================================
def substituicao_progressiva(L, B):
    n = len(B)
    y = np.zeros(n, dtype=float)

    for i in range(n):
        soma = B[i]
        for j in range(i):
            soma -= L[i, j] * y[j]
        y[i] = soma / L[i, i]

    return y


# =====================================================
# SUBSTITUIÇÃO REGRESSIVA (Back Substitution)
# Resolve o sistema triangular superior: U * x = y
# ======================================================
def substituicao_regressiva(U, y):
    n = len(y)
    x = np.zeros(n, dtype=float)

    x[n - 1] = y[n - 1] / U[n - 1, n - 1]
    for i in range(n - 2, -1, -1):
        soma = y[i]
        for j in range(i + 1, n):
            soma -= U[i, j] * x[j]
        x[i] = soma / U[i, i]

    return x

L, U = decomposicao_lu(A)
y = substituicao_progressiva(L, B)
x = substituicao_regressiva(U, y)


print("Matriz L: \n", L)
print("Matriz U: \n", U)
print("\n Vetor Intermediario y:\n", y)
print("\n Vetor Solução x: \n", x)









a = 0.5
b = 3.5

def custo(q):
    return x[0] + x[1] * q + x[2] * (q**2)

def receita(q):
    return 12.0 * q - 1.5 * q**2

def lucro(q):
    return receita(q) - custo(q)


def razao_aurea(a, b, tol=1e-6, max_iter=200):
    tau = (np.sqrt(5) - 1) / 2.0  # Razão Áurea (~0.618034)

    # 1. Escolha inicial dos dois pontos internos
    x1 = b - tau * (b - a)
    x2 = a + tau * (b - a)

    # 2. Avaliação inicial da função objetivo nos dois pontos
    f1 = lucro(x1)
    f2 = lucro(x2)

    eval_count = 2
    iter_count = 0

    # 3. Laço de redução do intervalo de busca
    while abs(b - a) > tol and iter_count < max_iter:
        if f1 < f2:
            # O mínimo está no subintervalo [x1, b]
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + tau * (b - a)
            f2 = lucro(x2)
        else:
            # O mínimo está no subintervalo [a, x2]
            b = x2
            x2 = x1
            f2 = f1
            x1 = b - tau * (b - a)
            f1 = lucro(x1)

        eval_count += 1
        iter_count += 1

    # Estimativa final da posição do ponto mínimo (ponto médio do intervalo final)
    Q_otimo = (a + b) / 2.0
    eval_count += 1

    return Q_otimo, lucro(Q_otimo), iter_count, eval_count



# 2. Chamando a função de otimização
Q_otimo, lucro_max, iters, evals = razao_aurea(a, b)

# 3. Printando os resultados
print("=" * 40)
print(" RESULTADOS DA OTIMIZAÇÃO (RAZÃO ÁUREA)")
print("=" * 40)
print(f"Vazão ótima (q*):      {Q_otimo:.5f} m³/h")
print(f"Lucro máximo L(q*):    $ {lucro_max:.2f}")
print(f"Custo no ponto C(q*):  $ {custo(Q_otimo):.2f}")
print(f"Receita no ponto R(q*):$ {receita(Q_otimo):.2f}")
print("=" * 40)