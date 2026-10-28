import numpy as np

# ==============================================================================
# MÉTODO DE NEWTON-RAPHSON 2x2 (FORMATO FIXO)
# Resolve o sistema não linear de vazões em tubulações paralelas.
# ==============================================================================
def newton_raphson_2x2(Q1_init, Q2_init, R1, R2, Q, tol=1e-8, max_iter=50):
    Q1 = Q1_init
    Q2 = Q2_init

    # Cabeçalho da tabela de acompanhamento
    print(
        f"{'k':<3} | {'Q1 (m³/s)':<11} | {'Q2 (m³/s)':<11} | {'F1 (m³/s)':<12} | {'F2 (m²/s²)':<12} | {'dQ1 (m³/s)':<11} | {'dQ2 (m³/s)':<11}"
    )
    print("-" * 88)

    convergiu = False

    for k in range(max_iter):
        # 1. AVALIAÇÃO DOS RESÍDUOS (SISTEMA F(Q) = 0)
        # F1: Continuidade da massa (Q1 + Q2 - Q = 0)
        # F2: Igualdade da perda de carga nos dois ramos (R1*Q1² - R2*Q2² = 0)
        F1 = Q1 + Q2 - Q
        F2 = R1 * (Q1**2) - R2 * (Q2**2)

        # 2. MONTAGEM DA MATRIZ JACOBIANA J = [[dF1/dQ1, dF1/dQ2], [dF2/dQ1, dF2/dQ2]]
        a = 1.0            # dF1/dQ1
        b = 1.0            # dF1/dQ2
        c = 2.0 * R1 * Q1  # dF2/dQ1
        d = -2.0 * R2 * Q2 # dF2/dQ2

        # Resíduos com sinal trocado para a equação J * dQ = -F
        r1 = -F1
        r2 = -F2

        # 3. SOLUÇÃO DO SISTEMA LINEAR 2x2 VIA REGRA DE CRAMER
        D = a * d - b * c                       # Determinante da Jacobiana

        if abs(D) < 1e-12:
            print("Erro: Matriz Jacobiana singular ou quase singular.")
            return None, None, k

        deltaQ1 = (r1 * d - r2 * b) / D         # Correção para Q1
        deltaQ2 = (r2 * a - r1 * c) / D         # Correção para Q2

        # Imprime os resultados da iteração k
        print(
            f"{k:<3} | {Q1:<11.8f} | {Q2:<11.8f} | {F1:<12.4e} | {F2:<12.4e} | {deltaQ1:<11.4e} | {deltaQ2:<11.4e}"
        )

        # 4. CRITÉRIO DE PARADA: Maior correção individual menor que a tolerância
        if max(abs(deltaQ1), abs(deltaQ2)) < tol:
            convergiu = True
            print("-" * 88)
            print(f"Convergência atingida na iteração {k}!")
            return Q1 + deltaQ1, Q2 + deltaQ2, k

        # 5. ATUALIZAÇÃO DAS VARIÁVEIS (CORREÇÃO APLICADA)
        Q1 += deltaQ1
        Q2 += deltaQ2

    print("-" * 88)
    print("O método não convergiu dentro do número máximo de iterações.")
    return Q1, Q2, max_iter














# ==============================================================================
# FUNÇÃO RESÍDUO GENÉRICA DO PROBLEMA DE TUBULAÇÕES
# Recebe o vetor de incógnitas X = [Q1, Q2] e os parâmetros (R1, R2, Q)
# Retorna o vetor de resíduos F(X) que deve zerar
# ==============================================================================
def residuos_tubulacao(X, R1, R2, Q):
    Q1, Q2 = X[0], X[1]
    
    f1 = Q1 + Q2 - Q               # Continuidade
    f2 = R1 * (Q1**2) - R2 * (Q2**2) # Igualdade de perda de carga
    
    return np.array([f1, f2], dtype=float)


# ==============================================================================
# MÉTODO DE NEWTON-RAPHSON GENÉRICO (MULTIVARIÁVEL N x N)
# Parâmetros:
# func     : Função que calcula o vetor de resíduos F(X, *params)
# X_init   : Vetor com o chute inicial [x0, x1, ..., xn]
# params   : Tupla de parâmetros passados para func
# tol      : Tolerância de convergência
# max_iter : Número máximo de iterações
# ==============================================================================
def newton_raphson_generico(func, X_init, params=(), tol=1e-8, max_iter=50):
    X = np.array(X_init, dtype=float)
    n = len(X)
    h = 1e-6  # Perturbação para derivada numérica (Jacobiana)

    print(f"{'k':<3} | {'Vetor X':<30} | {'Erro (max|dX|)':<15}")
    print("-" * 55)

    for k in range(max_iter):
        F = func(X, *params)

        # 1. CÁLCULO DA MATRIZ JACOBIANA POR DIFERENÇAS FINITAS
        J = np.zeros((n, n), dtype=float)
        for j in range(n):
            X_pert = X.copy()
            X_pert[j] += h
            F_pert = func(X_pert, *params)
            J[:, j] = (F_pert - F) / h  # Derivada parcial dFi/dXj

        # 2. RESOLUÇÃO DO SISTEMA LINEAR: J * dX = -F
        try:
            deltaX = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            print("Erro: Matriz Jacobiana é singular!")
            return X, k

        # Exibe o progresso da iteração
        erro = np.max(np.abs(deltaX))
        print(f"{k:<3} | {str(np.round(X, 6)):<30} | {erro:<15.4e}")

        # 3. ATUALIZAÇÃO DAS VARIÁVEIS
        X += deltaX

        # 4. CRITÉRIO DE PARADA
        if erro < tol:
            print("-" * 55)
            print(f"Convergência atingida na iteração {k}!")
            return X, k

    print("-" * 55)
    print("O método não convergiu no número máximo de iterações.")
    return X, max_iter