import numpy as np
import math as math

# Parametros 
S = 0.0002
B = 20
n = 0.03
Q = 5

Hl = 0.5
Hu = 1.0
max_iter = 100
tol = 0.0001

def f(H):
    coef = math.sqrt(S) * (B*H)**(5/3)
    den = n * (B+2*H)**(2/3)
    vazao_calc = coef / den
    return vazao_calc - Q


#PARTE 1:
#Metodo da Bisseção
def bissecao(Hl, Hu, max_iter, tol):
    #Verificação Inicial
    if f(Hl) * f(Hu) >= 0:
        print("Erro: F(Hl) * F(Hu) >= 0")
        return None, None

    print("METODO DA BISSEÇÃO")
    print(f"{'k':<5}{'Hl':<12}{'Hu':<12}{'Hr':<12}{'f(Hr)':<15}")
    print("-" * 60)

    for k in range (1, max_iter + 1):
        Hr = (Hl + Hu)/2.0
        f_Hr = f(Hr)

    print(f"{k:<5}{Hl:<12.6f}{Hu:<12.6f}{Hr:<12.6f}{f_Hr:<15.6e}")

    if abs(f_Hr) < tol:
        print(f"Convergencia atingida na interação {k}!")
        return Hr, k

    if f(Hl) * f_Hr < 0:
        Hu = Hr
    else:
        Hl = Hr

    print("Número máximo de iterações atingido sem convergência. ")
    return Hr, max_iter

#PARTE 2
#Metodo da Falsa Posição
def falsa_pos(Hl, Hu, max_iter, tol):
    # Verificação do sinal inicial
    if f(Hl) * f(Hu) >= 0:
        print("Erro: F(Hl) * F(Hu) >= 0")
        return None, None

    print("\n" + "=" * 60)
    print("MÉTODO DA FALSA POSIÇÃO")
    print(f"{'k':<5}{'Hl':<12}{'Hu':<12}{'Hr':<12}{'f(Hr)':<15}")

    
    for k in range (1, max_iter + 1):
        num = Hl*f(Hu) - Hu*f(Hl)
        deno = f(Hu) - f(Hl)
        Hr = num / deno
        f_Hr = f(Hr)

    print(f"{k:<5}{Hl:<12.6f}{Hu:<12.6f}{Hr:<12.6f}{f_Hr:<15.6e}")

    if abs (f_Hr) < tol:
        print(f"Convergencia máxima atingida na iteração {k}!")
        return Hr, k
    
    if f(Hl) * f(Hr) < 0.0:
        Hu = Hr
    else:
        Hl = Hr

    print("Número máximo de iterações atingido sem convergencia. ")
    return Hr, max_iter 


#EXECUÇÃO E RESULTADOS 
H_bissecao, iter_bissecao = bissecao(Hl, Hu, max_iter, tol)
H_falsa, iter_falsa = falsa_pos(Hl, Hu, max_iter, tol)

#Determinação da velocidade média U
if H_bissecao is not None and H_falsa is not None:
    U_bissecao = Q / (B * H_bissecao)
    U_falsa = Q / (B * H_falsa)

    print("\n" + "=" * 60)
    print("RESUMO DOS RESULTADOS")
    print("=" * 60)
    print(f"Bisseção:     H = {H_bissecao:.6f} m | Iterações = {iter_bissecao} | U = {U_bissecao:.6f} m/s")
    print(f"Falsa Posição:     H = {H_falsa:.6f} m | Iterações = {iter_falsa} | U = {U_falsa:.6f} m/s")
