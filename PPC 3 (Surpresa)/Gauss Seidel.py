import math as math
import numpy as np

n = 5
iter = 0 
max_iter = 100000
erro = 1000
tol = 0.000001

# fazer a temperatura de cada nó interno (média das temperaturas)
# T = (ti0 + Ti1 + tj0 + Tj1) / 4

T = np.full((n + 2, n + 2), 20.0)
# Ajuste a borda inferior para 100 (da coluna 1 até n):
T[n+1, 1:n+1] = 100.0 
print(T)


while erro > tol and iter < max_iter:
    erro = 0.0
    for j in range (1, n + 1):
        for i in range (1, n + 1):

            antigo = T[j, i]
            T[j, i] = (T[j, i-1] + T[j, i+1] + T[j-1, i] + T[j+1, i]) / 4

            diferenca = abs(T[j,i] - antigo)
            erro = max(erro, diferenca)

    iter = iter + 1

residuo = 0
for j in range (1, n + 1):
    for i in range (1, n + 1):
        r_ij = 4 * T[j, i] - (T[j-1, i] + T[j+1, i] + T[j, i-1] + T[j, i+1])
        residuo = max(residuo, abs(r_ij))


print("\n Tabela Atualizada")
print(T)
print('\n VERIFICAÇÃO:')
print('número de iterações: 'f"{iter}")
print('erro: ' f"{erro}")
print('temperatura do nó central: 'f"{T[3,3]}")
print("Residuo máximo: " f'{residuo}')

# a) Por que as temperaturas das bordas não são atualizadas? São condições de contorno prescritas (Dirichlet), de valor fixo conhecido, servindo como dados de entrada do sistema.

# b) O que diferencia a atualização de Gauss-Seidel da de Jacobi? Gauss-Seidel usa os valores recém-calculados imediatamente na mesma iteração, enquanto o método de Jacobi utiliza apenas os valores da iteração anterior k.

# c) Como o número de iterações muda quando a malha é refinada? O número de iterações aumenta significativamente, pois o sistema de equações fica maior e o passo espacial diminui, tornando a propagação das condições de contorno pelo domínio mais lenta.
