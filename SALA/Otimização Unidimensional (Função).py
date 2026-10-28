import numpy as np

# ==============================================================================
# MÉTODO DA RAZÃO ÁUREA (FORMATO FIXO)
# Encontra o mínimo de uma função unimodal no intervalo [a, b]
# ==============================================================================
def razao_aurea(a, b, tol=1e-6, max_iter=200):
    tau = (np.sqrt(5) - 1) / 2.0  # Razão Áurea (~0.618034)

    # 1. Escolha inicial dos dois pontos internos
    x1 = b - tau * (b - a)
    x2 = a + tau * (b - a)

    # 2. Avaliação inicial da função objetivo nos dois pontos
    f1 = apenas_custo(x1)
    f2 = apenas_custo(x2)

    eval_count = 2
    iter_count = 0

    # 3. Laço de redução do intervalo de busca
    while abs(b - a) > tol and iter_count < max_iter:
        if f1 > f2: # Esta Minimizando a função, caso < estaria maximizando
            
            # O mínimo está no subintervalo [x1, b]
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + tau * (b - a)
            f2 = apenas_custo(x2)
        else:
            # O mínimo está no subintervalo [a, x2]
            b = x2
            x2 = x1
            f2 = f1
            x1 = b - tau * (b - a)
            f1 = apenas_custo(x1)

        eval_count += 1
        iter_count += 1

    # Estimativa final da posição do ponto mínimo (ponto médio do intervalo final)
    emin = (a + b) / 2.0
    cmin = apenas_custo(emin)
    eval_count += 1

    return emin, cmin, iter_count, eval_count


# ==============================================================================
# MÉTODO DA INTERPOLAÇÃO QUADRÁTICA (FORMATO FIXO)
# Ajusta uma parábola por 3 pontos (x1 < x2 < x3) para estimar o vértice (mínimo)
# ==============================================================================
def interpolacao_quadratica(x1=0.04, x2=0.07, x3=0.11, tol=1e-6, max_iter=200):
    f1 = apenas_custo(x1)
    f2 = apenas_custo(x2)
    f3 = apenas_custo(x3)

    eval_count = 3
    iter_count = 0

    while (x3 - x1) > tol and iter_count < max_iter:
        # 1. Coeficientes do polinômio do 2º grau: P(x) = A*x^2 + B*x + C
        d1 = (f2 - f1) / (x2 - x1)
        A = ((f3 - f2) / (x3 - x2) - d1) / (x3 - x1)
        B = d1 - A * (x1 + x2)

        # Proteção contra concavidade nula ou negativa (A <= 0 não tem mínimo)
        if A <= 1e-12:
            break

        # 2. Ponto de mínimo da parábola ajustada (vértice)
        x_novo = -B / (2.0 * A)

        # Proteção contra ponto fora das fronteiras atuais
        if x_novo <= x1 or x_novo >= x3:
            break

        f_novo = apenas_custo(x_novo)
        eval_count += 1

        # 3. Regras de atualização mantendo o mínimo cercado
        if x_novo < x2:
            if f_novo < f2:
                x3, f3 = x2, f2
                x2, f2 = x_novo, f_novo
            else:
                x1, f1 = x_novo, f_novo
        else:
            if f_novo < f2:
                x1, f1 = x2, f2
                x2, f2 = x_novo, f_novo
            else:
                x3, f3 = x2, f2

        iter_count += 1

    emin = x2
    cmin = f2
    return emin, cmin, iter_count, eval_count














# ==============================================================================
# FUNÇÃO DE CUSTO TOTAL GENÉRICA (Isolamento de Tubulações)
# ==============================================================================
def custo_isolamento(e, r1, L, Ti, Tinf, k, h, t, pe, ci):
    # Resistências térmicas de condução e convecção
    Rcond = np.log((r1 + e) / r1) / (2.0 * np.pi * k)
    Rconv = 1.0 / (2.0 * np.pi * h * (r1 + e))
    qlinha = (Ti - Tinf) / (Rcond + Rconv)

    # Custo da energia perdida + Custo de instalação do isolante
    Ce = pe * ((L * qlinha * t) / 1000.0)
    Ci = ci * np.pi * L * ((r1 + e)**2 - r1**2)
    
    return Ce + Ci


# ==============================================================================
# MÉTODO DA RAZÃO ÁUREA GENÉRICO
# Parâmetros:
# func   : Função objetivo a ser MINIMIZADA: func(x, *params)
# a, b   : Limites do intervalo de busca [a, b]
# params : Tupla com parâmetros adicionais passados para func
# ==============================================================================
def razao_aurea_generico(func, a, b, params=(), tol=1e-6, max_iter=200):
    tau = (np.sqrt(5) - 1) / 2.0

    x1 = b - tau * (b - a)
    x2 = a + tau * (b - a)

    f1 = func(x1, *params)
    f2 = func(x2, *params)

    eval_count = 2
    iter_count = 0

    while abs(b - a) > tol and iter_count < max_iter:
        if f1 > f2:
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + tau * (b - a)
            f2 = func(x2, *params)
        else:
            b = x2
            x2 = x1
            f2 = f1
            x1 = b - tau * (b - a)
            f1 = func(x1, *params)

        eval_count += 1
        iter_count += 1

    xmin = (a + b) / 2.0
    fmin = func(xmin, *params)
    eval_count += 1

    return xmin, fmin, iter_count, eval_count


# ==============================================================================
# MÉTODO DA INTERPOLAÇÃO QUADRÁTICA GENÉRICO
# Parâmetros:
# func       : Função objetivo a ser MINIMIZADA: func(x, *params)
# x1, x2, x3 : Três pontos iniciais (deve cumprir x1 < x2 < x3)
# params     : Tupla com parâmetros adicionais passados para func
# ==============================================================================
def interpolacao_quadratica_generico(func, x1, x2, x3, params=(), tol=1e-6, max_iter=200):
    f1 = func(x1, *params)
    f2 = func(x2, *params)
    f3 = func(x3, *params)

    eval_count = 3
    iter_count = 0

    while (x3 - x1) > tol and iter_count < max_iter:
        d1 = (f2 - f1) / (x2 - x1)
        A = ((f3 - f2) / (x3 - x2) - d1) / (x3 - x1)
        B = d1 - A * (x1 + x2)

        # Se A <= 0, a curvatura não possui mínimo relativo válido
        if A <= 1e-12:
            break

        x_novo = -B / (2.0 * A)

        if x_novo <= x1 or x_novo >= x3:
            break

        f_novo = func(x_novo, *params)
        eval_count += 1

        if x_novo < x2:
            if f_novo < f2:
                x3, f3 = x2, f2
                x2, f2 = x_novo, f_novo
            else:
                x1, f1 = x_novo, f_novo
        else:
            if f_novo < f2:
                x1, f1 = x2, f2
                x2, f2 = x_novo, f_novo
            else:
                x3, f3 = x2, f2

        iter_count += 1

    xmin = x2
    fmin = f2
    return xmin, fmin, iter_count, eval_count