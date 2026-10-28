import math

# ==============================================================================
# MÉTODO DA BISSEÇÃO
# Divide o intervalo ao meio sucessivamente até encontrar a raiz.
# ==============================================================================
def bissecao(Hl, Hu, max_iter, tol):
    # 1. TESTE DE BOLZANO: Verifica se existe uma raiz no intervalo [Hl, Hu].
    # Se f(Hl) e f(Hu) tiverem o mesmo sinal, o produto será positivo (>= 0),
    # o que significa que o método não garante encontrar uma raiz.
    if f(Hl) * f(Hu) >= 0:
        print("Erro: f(Hl) * f(Hu) >= 0 (Não há troca de sinal no intervalo)")
        return None, None

    # Cabeçalho para formatar a tabela de resultados do processo iterativo
    print("\nMÉTODO DA BISSEÇÃO")
    print(f"{'k':<5}{'Hl':<12}{'Hu':<12}{'Hr':<12}{'f(Hr)':<15}")
    print("-" * 60)

    # 2. LAÇO ITERATIVO: Executa até o número máximo de iterações (max_iter)
    for k in range(1, max_iter + 1):
        # Ponto médio do intervalo atual
        Hr = (Hl + Hu) / 2.0
        
        # Avalia a função no ponto médio
        f_Hr = f(Hr)

        # Imprime os valores da iteração atual formatados
        print(f"{k:<5}{Hl:<12.6f}{Hu:<12.6f}{Hr:<12.6f}{f_Hr:<15.6e}")

        # CRITÉRIO DE PARADA (TOLERÂNCIA):
        # Se o valor absoluto de f(Hr) for menor que a tolerância, encontramos a raiz!
        if abs(f_Hr) < tol:
            print(f"Convergência atingida na iteração {k}!")
            return Hr, k  # Retorna a raiz encontrada e quantas iterações demorou

        # ATUALIZAÇÃO DO INTERVALO:
        # Se f(Hl) e f(Hr) têm sinais opostos, a raiz está na metade esquerda [Hl, Hr]
        if f(Hl) * f_Hr < 0:
            Hu = Hr  # O limite superior vira o ponto médio
        else:
            Hl = Hr  # Caso contrário, a raiz está na metade direita [Hr, Hu]

    # Caso chegue ao fim do laço sem atingir a tolerância exigida
    print("Número máximo de iterações atingido sem convergência.")
    return Hr, max_iter


# ==============================================================================
# MÉTODO DA FALSA POSIÇÃO
# Usa uma reta secante entre f(Hl) e f(Hu) para estimar a raiz (interseção com Y=0).
# ==============================================================================
def falsa_pos(Hl, Hu, max_iter, tol):
    # 1. TESTE DE BOLZANO: Garante troca de sinal inicial
    if f(Hl) * f(Hu) >= 0:
        print("Erro: f(Hl) * f(Hu) >= 0 (Não há troca de sinal no intervalo)")
        return None, None

    # Cabeçalho da tabela de acompanhamento
    print("\nMÉTODO DA FALSA POSIÇÃO")
    print(f"{'k':<5}{'Hl':<12}{'Hu':<12}{'Hr':<12}{'f(Hr)':<15}")
    print("-" * 60)

    # 2. LAÇO ITERATIVO
    for k in range(1, max_iter + 1):
        # Fórmula da Falsa Posição (interseção da reta secante com o eixo X)
        num = Hl * f(Hu) - Hu * f(Hl)
        deno = f(Hu) - f(Hl)
        Hr = num / deno
        
        # Avalia o erro/resíduo da função na nova estimativa
        f_Hr = f(Hr)

        # Exibe os dados da iteração k
        print(f"{k:<5}{Hl:<12.6f}{Hu:<12.6f}{Hr:<12.6f}{f_Hr:<15.6e}")

        # CRITÉRIO DE PARADA: Verifica se o erro atende à tolerância desejada
        if abs(f_Hr) < tol:
            print(f"Convergência atingida na iteração {k}!")
            return Hr, k

        # ATUALIZAÇÃO DO INTERVALO: Mantém a raiz cercada
        if f(Hl) * f_Hr < 0:
            Hu = Hr  # A raiz está entre Hl e Hr
        else:
            Hl = Hr  # A raiz está entre Hr e Hu

    print("Número máximo de iterações atingido sem convergência.")
    return Hr, max_iter









# ==============================================================================
# FUNÇÃO OBJETIVO (Equação de Manning adaptada)
# f(H) = Vazão Calculada - Vazão Desejada (Q) = 0
# Parâmetros:
# H = Profundidade do canal (m) [Variável a ser encontrada]
# S = Declividade do canal (m/m)
# B = Largura da base (m)
# n = Coeficiente de rugosidade de Manning
# Q = Vazão de projeto (m³/s)
# ==============================================================================
def f(H, S, B, n, Q):
    coef = math.sqrt(S) * (B * H)**(5/3)    # Numerador da equação de Manning
    den = n * (B + 2 * H)**(2/3)            # Denominador (Rugosidade x Perímetro Molhado^(2/3))
    vazao_calc = coef / den
    return vazao_calc - Q                   # Queremos zerar essa diferença


# ==============================================================================
# MÉTODO DA BISSEÇÃO (PARAMETRIZADO)
# ==============================================================================
def bissecao(Hl, Hu, S, B, n, Q, max_iter=100, tol=0.0001):
    # Passamos os parâmetros (S, B, n, Q) para cada avaliação de f()
    if f(Hl, S, B, n, Q) * f(Hu, S, B, n, Q) >= 0:
        print("Erro: f(Hl) * f(Hu) >= 0 (Não há troca de sinal no intervalo)")
        return None, None

    print("\nMÉTODO DA BISSEÇÃO")
    print(f"{'k':<5}{'Hl':<12}{'Hu':<12}{'Hr':<12}{'f(Hr)':<15}")
    print("-" * 60)

    for k in range(1, max_iter + 1):
        Hr = (Hl + Hu) / 2.0                 # Ponto médio do intervalo
        f_Hr = f(Hr, S, B, n, Q)             # Calcula f(Hr) com os parâmetros do canal

        print(f"{k:<5}{Hl:<12.6f}{Hu:<12.6f}{Hr:<12.6f}{f_Hr:<15.6e}")

        # Se f(Hr) for suficientemente próximo de zero, encerra
        if abs(f_Hr) < tol:
            print(f"Convergência atingida na iteração {k}!")
            return Hr, k

        # Redução do intervalo mantendo a troca de sinal
        if f(Hl, S, B, n, Q) * f_Hr < 0:
            Hu = Hr
        else:
            Hl = Hr

    print("Número máximo de iterações atingido sem convergência.")
    return Hr, max_iter


# ==============================================================================
# MÉTODO DA FALSA POSIÇÃO (PARAMETRIZADO E OTIMIZADO)
# ==============================================================================
def falsa_pos(Hl, Hu, S, B, n, Q, max_iter=100, tol=0.0001):
    # Armazena os valores nos limites para evitar recomputar f() repetidamente
    f_Hl = f(Hl, S, B, n, Q)
    f_Hu = f(Hu, S, B, n, Q)

    # Teste inicial de troca de sinal
    if f_Hl * f_Hu >= 0:
        print("Erro: f(Hl) * f(Hu) >= 0 (Não há troca de sinal no intervalo)")
        return None, None

    print("\nMÉTODO DA FALSA POSIÇÃO")
    print(f"{'k':<5}{'Hl':<12}{'Hu':<12}{'Hr':<12}{'f(Hr)':<15}")
    print("-" * 60)

    for k in range(1, max_iter + 1):
        # Avalia a função nos extremos do intervalo atual
        f_Hl = f(Hl, S, B, n, Q)
        f_Hu = f(Hu, S, B, n, Q)

        # Fórmula da Falsa Posição
        num = Hl * f_Hu - Hu * f_Hl
        deno = f_Hu - f_Hl
        Hr = num / deno
        
        f_Hr = f(Hr, S, B, n, Q)

        print(f"{k:<5}{Hl:<12.6f}{Hu:<12.6f}{Hr:<12.6f}{f_Hr:<15.6e}")

        # Teste de convergência
        if abs(f_Hr) < tol:
            print(f"Convergência atingida na iteração {k}!")
            return Hr, k

        # Atualização das fronteiras
        if f_Hl * f_Hr < 0:
            Hu = Hr
        else:
            Hl = Hr

    print("Número máximo de iterações atingido sem convergência.")
    return Hr, max_iter