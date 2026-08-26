
# Cálculo Numérico Aplicado - [PPC1]

# Descrição
Este projeto analisa numericamente o processo de sedimentação de uma partícula esférica imersa em um fluido, movimentando-se sob a ação da gravidade e de forças hidrodinâmicas.

A modelagem física abrange desde o regime assintótico puro de Stokes (onde as forças viscosas predominam e os efeitos inerciais são desprezíveis) até cenários em que pequenos efeitos inerciais do fluido passam a atuar sobre a partícula, aplicando a correção não-linear proposta por Oseen. O comportamento do sistema é determinado pela relação entre o Número de Stokes (massa e relaxação da partícula) e o Número de Reynolds de partícula (influência das forças inerciais do fluido).

# Tecnologias Utilizadas
- **Linguagem:** Python 3
- **Ambiente:** Google Colab / Jupyter Notebook
- **Bibliotecas:** `numpy`, `matplotlib`

# Como Executar o Código

## Google Colab
1. Clique no notebook [`PPC_1.ipynb`](./PPC_1.ipynb) presente neste repositório.
2. Clique no botão **"Open in Colab"** no topo da página.
3. No Colab, vá no menu superior e clique em **Ambiente de execução > Executar tudo** (ou `Ctrl + F9`).
