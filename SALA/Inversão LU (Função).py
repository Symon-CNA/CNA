import numpy as np

# ==============================================================================
# DECOMPOSIÇÃO LU (Eliminação de Gauss sem pivoteamento)
# Recebe: Matriz A (n x n)
# Retorna: Matriz L (Triangular Inferior) e Matriz U (Triangular Superior)
# ==============================================================================
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


# ==============================================================================
# SUBSTITUIÇÃO PROGRESSIVA (Forward Substitution)
# Resolve o sistema triangular inferior: L * y = b
# ==============================================================================
def substituicao_progressiva(L, b):
    n = len(b)
    y = np.zeros(n, dtype=float)

    for i in range(n):
        soma = b[i]
        for j in range(i):
            soma -= L[i, j] * y[j]
        y[i] = soma / L[i, i]

    return y


# ==============================================================================
# SUBSTITUIÇÃO REGRESSIVA (Back Substitution)
# Resolve o sistema triangular superior: U * x = y
# ==============================================================================
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


# ==============================================================================
# INVERSÃO DE MATRIZ VIA DECOMPOSIÇÃO LU
# Calcula a matriz inversa (C = A^-1) coluna por coluna resolvendo A * c_j = e_j
# ==============================================================================
def inverter_matriz_lu(A):
    n = A.shape[0]
    
    # 1. Fatora a matriz A em L e U (uma única vez)
    L, U = decomposicao_lu(A)
    
    C = np.zeros((n, n), dtype=float)

    # 2. Resolve o sistema para cada vetor coluna unitário e_j
    for coluna in range(n):
        e = np.zeros(n, dtype=float)
        e[coluna] = 1.0  # Vetor da base canônica (coluna j da matriz identidade)

        # Passo 1: L * y = e
        y = substituicao_progressiva(L, e)
        
        # Passo 2: U * x = y
        x = substituicao_regressiva(U, y)

        # Armazena a solução x na coluna correspondente da matriz inversa C
        C[:, coluna] = x

    return C










import numpy as np

# ==============================================================================
# MONTAGEM GENÉRICA DE MATRIZ DE RIGIDEZ TRIDIAGONAL (N Graus de Liberdade)
# Recebe: Lista ou array com os valores de rigidez das molas [k0, k1, k2, ..., kn]
# Retorna: Matriz de rigidez K de tamanho (n x n)
# ==============================================================================
def montar_matriz_rigidez(k_valores):
    # Quantidade de molas = n + 1 (onde n é o número de graus de liberdade)
    n = len(k_valores) - 1
    K = np.zeros((n, n), dtype=float)

    for i in range(n):
        # Termos da diagonal principal: soma das duas molas conectadas ao nó i
        K[i, i] = k_valores[i] + k_valores[i + 1]
        
        # Termos da acoplamento (diagonais secundárias superior e inferior)
        if i > 0:
            K[i, i - 1] = -k_valores[i]
        if i < n - 1:
            K[i, i + 1] = -k_valores[i + 1]

    return K


# ==============================================================================
# SOLVER LU GENÉRICO
# Resolve A * x = b (ou inverte A se b não for fornecido) usando Decomposição LU
# Parâmetros:
# A : Matriz quadrada de coeficientes (n x n)
# b : Vetor de termos independentes/forças (Opcional - se omitido, calcula a Inversa)
# ==============================================================================
def resolver_sistema_lu(A, b=None):
    n = A.shape[0]

    # Decomposição LU interna
    L = np.eye(n, dtype=float)
    U = A.copy().astype(float)

    for k in range(n - 1):
        for i in range(k + 1, n):
            fator = U[i, k] / U[k, k]
            L[i, k] = fator
            for j in range(k + 1, n):
                U[i, j] -= fator * U[k, j]
            U[i, k] = 0.0

    # Subfunções internas para resolver os sistemas triangulares
    def fwd(L, vec):
        y = np.zeros(n, dtype=float)
        for i in range(n):
            soma = vec[i]
            for j in range(i):
                soma -= L[i, j] * y[j]
            y[i] = soma / L[i, i]
        return y

    def bwd(U, vec):
        x = np.zeros(n, dtype=float)
        x[n - 1] = vec[n - 1] / U[n - 1, n - 1]
        for i in range(n - 2, -1, -1):
            soma = vec[i]
            for j in range(i + 1, n):
                soma -= U[i, j] * x[j]
            x[i] = soma / U[i, i]
        return x

    # MODALIDADE 1: Se o vetor de carga b for fornecido, resolve diretamente u = A^-1 * b
    if b is not None:
        y = fwd(L, b)
        x = bwd(U, y)
        return x

    # MODALIDADE 2: Se b for omitido, calcula e retorna a Matriz Inversa (A^-1)
    else:
        A_inv = np.zeros((n, n), dtype=float)
        for col in range(n):
            e = np.zeros(n, dtype=float)
            e[col] = 1.0
            y = fwd(L, e)
            x = bwd(U, y)
            A_inv[:, col] = x
        return A_inv