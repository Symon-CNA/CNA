import numpy as np

# ==============================================================================
# MÉTODO DE NEWTON-RAPHSON (1D - FORMATO FIXO)
# Requer a função f(T) e sua derivada analítica df(T)
# ==============================================================================
def newton_raphson(T_init, tol=1e-6, max_iter=100):
    T_k = T_init

    print("=" * 60)
    print("MÉTODO DE NEWTON-RAPHSON")
    print("=" * 60)
    print(f"{'k':<5} | {'T_k (K)':<12} | {'f(T_k)':<15} | {'Erro':<15}")
    print("-" * 55)

    for k in range(max_iter):
        f_val = f(T_k)
        df_val = df(T_k)

        # Evita divisão por zero caso a derivada seja nula
        if abs(df_val) < 1e-12:
            print("Erro: Derivada próxima de zero!")
            return None, k

        # Atualização pela tangente
        T_next = T_k - f_val / df_val
        erro = abs(T_next - T_k)

        print(f"{k:<5} | {T_k:<12.6f} | {f_val:<15.6f} | {erro:<15.6e}")

        # Critério de convergência
        if erro < tol:
            print(f"\nConvergência alcançada na iteração {k + 1}!")
            print(f"Temperatura de equilíbrio T = {T_next:.6f} K\n")
            return T_next, k + 1

        T_k = T_next

    print("\nO método não convergiu no número máximo de iterações.")
    return T_k, max_iter


# ==============================================================================
# MÉTODO DA SECANTE (FORMATO FIXO)
# Aproxime a derivada usando dois pontos iniciais (T0, T1)
# ==============================================================================
def secante(T0_init, T1_init, tol=1e-6, max_iter=100):
    T0 = T0_init
    T1 = T1_init

    print("=" * 60)
    print("MÉTODO DA SECANTE")
    print("=" * 60)
    print(f"{'k':<5} | {'T_k (K)':<12} | {'f(T_k)':<15} | {'Erro (K)':<15}")
    print("-" * 55)

    # Impressão do estado inicial para k = 0
    print(f"{0:<5} | {T0:<12.6f} | {f(T0):<15.6f} | {'-':<15}")

    for k in range(1, max_iter):
        f0 = f(T0)
        f1 = f(T1)

        # Evita divisão por zero se f(T1) == f(T0)
        if abs(f1 - f0) < 1e-12:
            print("Erro: f(T1) - f(T0) muito próximo de zero!")
            return None, k

        # Fórmula da reta secante
        T2 = T1 - f1 * (T1 - T0) / (f1 - f0)
        erro = abs(T2 - T1)

        print(f"{k:<5} | {T1:<12.6f} | {f1:<15.6f} | {erro:<15.6e}")

        # Critério de parada
        if erro < tol:
            print(f"\nConvergência alcançada na iteração {k}!")
            print(f"Temperatura de equilíbrio T = {T2:.6f} K\n")
            return T2, k

        # Deslocamento dos pontos para a próxima iteração
        T0, T1 = T1, T2

    print("\nO método não convergiu no número máximo de iterações.")
    return T1, max_iter













# ==============================================================================
# FUNÇÃO RESÍDUO DO BALANÇO TÉRMICO (Convecção + Radiação = Potência)
# f(T) = h*A*(T - T_inf) + epsi*Sig*A*(T^4 - T_inf^4) - P = 0
# ==============================================================================
def f_termico(T, A, h, epsi, T_inf, P, Sig=5.67e-8):
    q_conv = h * A * (T - T_inf)
    q_rad = epsi * Sig * A * (T**4 - T_inf**4)
    return q_conv + q_rad - P


# ==============================================================================
# MÉTODO DE NEWTON-RAPHSON GENÉRICO (com derivada por diferença finita opcional)
# Parâmetros:
# func     : Função objetivo f(x, *params)
# x_init   : Chute inicial
# params   : Tupla contendo os parâmetros da função (A, h, epsi, T_inf, P, Sig)
# df_func  : Função de derivada analítica f'(x, *params) [Opcional]
# ==============================================================================
def newton_raphson_generico(func, x_init, params=(), df_func=None, tol=1e-6, max_iter=100):
    x_k = float(x_init)

    print("=" * 60)
    print("MÉTODO DE NEWTON-RAPHSON (GENÉRICO)")
    print("=" * 60)
    print(f"{'k':<5} | {'x_k':<12} | {'f(x_k)':<15} | {'Erro':<15}")
    print("-" * 55)

    for k in range(max_iter):
        f_val = func(x_k, *params)

        # Se não houver função de derivada analítica, calcula numericamente
        if df_func is not None:
            df_val = df_func(x_k, *params)
        else:
            dx = 1e-6
            df_val = (func(x_k + dx, *params) - func(x_k - dx, *params)) / (2 * dx)

        if abs(df_val) < 1e-12:
            print("Erro: Derivada nula ou insignificante.")
            return None, k

        x_next = x_k - f_val / df_val
        erro = abs(x_next - x_k)

        print(f"{k:<5} | {x_k:<12.6f} | {f_val:<15.6f} | {erro:<15.6e}")

        if erro < tol:
            print(f"\nConvergência alcançada na iteração {k + 1}!")
            return x_next, k + 1

        x_k = x_next

    print("\nO método não convergiu no número máximo de iterações.")
    return x_k, max_iter


# ==============================================================================
# MÉTODO DA SECANTE GENÉRICO
# Parâmetros:
# func     : Função objetivo f(x, *params)
# x0_init  : Primeiro ponto inicial
# x1_init  : Segundo ponto inicial
# params   : Tupla contendo os parâmetros da função
# ==============================================================================
def secante_generico(func, x0_init, x1_init, params=(), tol=1e-6, max_iter=100):
    x0 = float(x0_init)
    x1 = float(x1_init)

    print("=" * 60)
    print("MÉTODO DA SECANTE (GENÉRICO)")
    print("=" * 60)
    print(f"{'k':<5} | {'x_k':<12} | {'f(x_k)':<15} | {'Erro':<15}")
    print("-" * 55)

    print(f"{0:<5} | {x0:<12.6f} | {func(x0, *params):<15.6f} | {'-':<15}")

    for k in range(1, max_iter):
        f0 = func(x0, *params)
        f1 = func(x1, *params)

        if abs(f1 - f0) < 1e-12:
            print("Erro: Divisão por zero em f(x1) - f(x0).")
            return None, k

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        erro = abs(x2 - x1)

        print(f"{k:<5} | {x1:<12.6f} | {f1:<15.6f} | {erro:<15.6e}")

        if erro < tol:
            print(f"\nConvergência alcançada na iteração {k}!")
            return x2, k

        x0, x1 = x1, x2

    print("\nO método não convergiu no número máximo de iterações.")
    return x1, max_iter