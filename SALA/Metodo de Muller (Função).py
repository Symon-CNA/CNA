# ==============================================================================
# MÉTODO DE MÜLLER (FORMATO FIXO)
# Encontra raízes aproximando a função por um polinômio do 2º grau através de
# 3 pontos iniciais (H0, H1, H2).
# ==============================================================================
def muller(Q, H0, H1, H2, tol=1e-8, max_iter=100, verbose=False):
    # Imprime o cabeçalho da tabela caso o rastreamento (verbose) esteja ativo
    if verbose:
        print(
            f"{'k':<5}{'H0':<10}{'H1':<10}{'H2':<10}{'H3 (Nova)':<12}{'f(H3)':<15}{'Erro':<12}"
        )

    # Laço principal de iterações
    for k in range(1, max_iter + 1):
        # 1. Avaliação da função nos três pontos atuais
        f0 = f(H0, Q)
        f1 = f(H1, Q)
        f2 = f(H2, Q)

        # 2. Diferenças entre os pontos no eixo X
        h0 = H1 - H0
        h1 = H2 - H1

        # 3. Primeiras diferenças divididas (declividades das secantes)
        delta0 = (f1 - f0) / h0
        delta1 = (f2 - f1) / h1

        # 4. Coeficientes da parábola a*(H - H2)^2 + bM*(H - H2) + c = 0
        # CORREÇÃO: (h1 + h0) no denominador
        a = (delta1 - delta0) / (h1 + h0)
        bM = (a * h1) + delta1
        c = f2  # O termo constante é o valor da função no último ponto (H2)

        # 5. Cálculo do discriminante (com tratamento de números complexos)
        # CORREÇÃO: bM**2 - 4*a*c
        D = (bM**2 - 4 * a * c + 0j)**0.5

        # 6. Escolha do sinal do denominador para MAXIMIZAR o valor absoluto
        # Isso evita o cancelamento numérico e escolhe a raiz mais próxima de H2
        if abs(bM + D) > abs(bM - D):
            E = bM + D
        else:
            E = bM - D

        # 7. Incremento delta H e nova estimativa da raiz (H3)
        deltah = (-2 * c) / E
        H3 = H2 + deltah

        # 8. Cálculo do erro relativo/absoluto entre a nova estimativa e H2
        erro = abs(H3 - H2)

        # Imprime o progresso da iteração atual
        if verbose:
            f_H3 = f(H3, Q)
            print(
                f"{k:<5}{H0:<10.5f}{H1:<10.5f}{H2:<10.5f}{H3.real:<12.6f}{f_H3.real:<15.4e}{erro.real:<12.4e}"
            )

        # CRITÉRIO DE PARADA: Verifica se o erro atinge a tolerância desejada
        if erro < tol:
            return H3.real, k

        # Atualização dos pontos para a próxima iteração (desloca a janela de 3 pontos)
        H0 = H1
        H1 = H2
        H2 = H3

    raise ValueError("O método não convergiu no número máximo de iterações.")












# ==============================================================================
# FUNÇÃO OBJETIVO (Equação da Energia / Profundidade Crítica)
# f(H) = 0
# Parâmetros:
# H = Profundidade do canal (m) [Variável a encontrar]
# Q = Vazão (m³/s)
# b = Largura do fundo do canal trapezoidal (m)
# z = Talude lateral (H:V)
# g = Aceleração da gravidade (m/s²)
# ==============================================================================
def f_manning(H, Q, b=4.0, z=1.5, g=9.81):
    termo1 = g * (z**3) * (H**6)
    termo2 = 3 * b * g * (z**2) * (H**5)
    termo3 = 3 * g * (b**2) * z * (H**4)
    termo4 = g * (b**3) * (H**3)
    termo5 = 2 * (Q**2) * z * H
    termo6 = (Q**2) * b
    
    return termo1 + termo2 + termo3 + termo4 - termo5 - termo6


# ==============================================================================
# MÉTODO DE MÜLLER PARAMETRIZADO E GENÉRICO
# Parâmetros:
# func     : Função Python que calcula f(x, *params)
# x0, x1, x2 : Três estimativas iniciais para cercar/aproximar a raiz
# params   : Tupla com parâmetros adicionais passados para func (ex: (Q, b, z, g))
# tol      : Tolerância para convergência
# max_iter : Número máximo de iterações permitidas
# verbose  : Exibe a tabela de iterações se True
# ==============================================================================
def muller_generico(func, x0, x1, x2, params=(), tol=1e-8, max_iter=100, verbose=False):
    if verbose:
        print(
            f"{'k':<5}{'x0':<10}{'x1':<10}{'x2':<10}{'x3 (Nova)':<12}{'f(x3)':<15}{'Erro':<12}"
        )

    for k in range(1, max_iter + 1):
        # 1. Avalia a função genérica aplicando os parâmetros (*params)
        f0 = func(x0, *params)
        f1 = func(x1, *params)
        f2 = func(x2, *params)

        # 2. Variações de X e primeiras diferenças divididas
        h0 = x1 - x0
        h1 = x2 - x1
        delta0 = (f1 - f0) / h0
        delta1 = (f2 - f1) / h1

        # 3. Coeficientes do polinômio quadrático interpolador
        a = (delta1 - delta0) / (h1 + h0)
        bM = (a * h1) + delta1
        c = f2

        # 4. Cálculo do discriminante (suporta raízes complexas se necessário)
        D = (bM**2 - 4 * a * c + 0j)**0.5

        # 5. Escolha do denominador de maior módulo para estabilidade numérica
        if abs(bM + D) > abs(bM - D):
            E = bM + D
        else:
            E = bM - D

        # 6. Atualização do ponto x3
        deltax = (-2 * c) / E
        x3 = x2 + deltax
        erro = abs(x3 - x2)

        # 7. Exibição das iterações
        if verbose:
            f_x3 = func(x3, *params)
            print(
                f"{k:<5}{x0:<10.5f}{x1:<10.5f}{x2:<10.5f}{x3.real:<12.6f}{f_x3.real:<15.4e}{erro.real:<12.4e}"
            )

        # 8. Critério de convergência
        if erro < tol:
            return x3.real, k

        # Deslocamento dos pontos para o próximo ciclo
        x0, x1, x2 = x1, x2, x3

    raise ValueError("O método não convergiu no número máximo de iterações.")