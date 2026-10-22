import math as math 
import numpy as np

r1 = 0.04
L = 20.0
Ti = 180.0
Tinf = 40.0
k = 0.045
h = 12.0
t = 30000
pe = 0.70
ci = 9000.0

e = 0.0
max_e = 0.12

#DEFINIR PERDA DE CALOR
def taxa_perda_calor(e):
    Rcond = (np.log((r1+e)/r1)) / (2*np.pi*k)
    Rconv = 1 / (2*np.pi*h*(r1+e))
    qlinha = (Ti - Tinf) / (Rcond + Rconv)
    return qlinha

# DEFINIR CUSTO TOTAL
def custo_total(e):
    q_e = taxa_perda_calor(e)

    Ce = pe*((L*q_e*t)/1000)
    Ci = ci*np.pi*L*((r1+e)**2 - r1**2)
    Ct = Ce + Ci
    return Ce, Ci, Ct


# PARTE 1: EXPLORAÇÃO NUMÉRICA
print(f"{'e (m)':<8} | {'Ce (R$)':<12} | {'Ci (R$)':<12} | {'Custo Total (R$)':<15}")
print("-" * 55)

# Explorar numericamente 
while e <= max_e + 0.00001:
    Ce, Ci, Ct = custo_total(e)
    print(f"{e:<8.2f} | {Ce:<12.2f} | {Ci:<12.2f} | {Ct:<15.2f}")
    e = round(e + .02, 2)


# PARTE 2: MÉTODO DA RAZÃO ÁUREA
a = 0.0
b = 0.12
tol = 0.000001
max_iter = 200
tau = (np.sqrt(5)-1) / 2

def apenas_custo(e):
    Ce, Ci, Ct = custo_total(e)
    return Ct

def razao_aurea(a, b, tol, max_iter):
    x1 = b - tau*(b-a)
    x2 = a + tau*(b-a)

    f1 = apenas_custo(x1)
    f2 = apenas_custo(x2)

    eval_count = 2 
    iter = 0

    while abs(b-a) > tol and iter < max_iter:
        if (f1 > f2):
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + tau*(b-a)
            f2 = apenas_custo(x2)
            eval_count += 1

        else:
            b = x2
            x2 = x1
            f2 = f1
            x1 = b - tau*(b-a)
            f1 = apenas_custo(x1)
            eval_count += 1

        iter = iter + 1

    emin = (a + b) / 2.0
    cmin = apenas_custo(emin)
    eval_count += 1

    return emin, apenas_custo(emin), iter, eval_count

emin, cmin, total_iter, total_eval = razao_aurea(a, b, tol, max_iter)

print(f"Espessura ótima: {emin:.6f} m")
print(f"Custo mínimo: R$ {cmin:.2f}")
print(f"Iterações: {total_iter}")
print(f"Avaliações de C(e): {total_eval}")


# PARTE 3: INTERPOLAÇÃO QUADRÁTICA
def interpolacao_quadratica(x1=0.04, x2=0.07, x3=0.11, tol=1e-6, max_iter=200):
    f1 = apenas_custo(x1)
    f2 = apenas_custo(x2)
    f3 = apenas_custo(x3)
    eval_count = 3
    iter_count = 0

    while (x3 - x1) > tol and iter_count < max_iter:
        d1 = (f2 - f1) / (x2 - x1)
        denom_A = (x3 - x1)
        A = ((f3 - f2) / (x3 - x2) - d1) / denom_A
        B = d1 - A * (x1 + x2)

        # Proteção contra curvatura inválida (A <= 0 ou valor muito próximo de zero)
        if A <= 1e-12:
            break

        x_novo = -B / (2 * A)

        # Proteção se o vértice estiver fora do intervalo
        if x_novo <= x1 or x_novo >= x3:
            break

        f_novo = apenas_custo(x_novo)
        eval_count += 1

        # Regras de atualização do intervalo
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

# EXECUÇÃO E TABELA COMPARATIVA
e_ra, c_ra, iter_ra, eval_ra = razao_aurea(a, b, tol, max_iter)
e_iq, c_iq, iter_iq, eval_iq = interpolacao_quadratica(0.04, 0.07, 0.11, tol, max_iter)

print("\n" + "=" * 55)
print("=== PARTE 3: COMPARAÇÃO DOS MÉTODOS ===")
print("=" * 55)
print(f"{'Resultado':<25} | {'Razão Áurea':<15} | {'Interp. Quadrática':<15}")
print("-" * 61)
print(f"{'Espessura ótima (m)':<25} | {e_ra:<15.6f} | {e_iq:<15.6f}")
print(f"{'Custo mínimo (R$)':<25} | {c_ra:<15.2f} | {c_iq:<15.2f}")
print(f"{'Número de iterações':<25} | {iter_ra:<15d} | {iter_iq:<15d}")
print(f"{'Avaliações de C(e)':<25} | {eval_ra:<15d} | {eval_iq:<15d}")
