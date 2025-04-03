"""
Criado na Sexta, dia 17 de Julho de 2020, às 13:15:47.

Dados de entrada para a geração da treliça:
    - coordenadas dos nós: coordsNos
    - incidêcia dos elementos: incElem
    - materiais das barras: materiais
    - seções transversais das barras: secoes
    - cargas dos nós da estrutura: cargas
    - apoios dos nós da estrutura: apoios

* as variáveis coordsNos, incElem, materiais, secoes, cargas e apoios deem existir obrigatoriamente neste arquivo pois são elas que levam os dados de entrada para serem calculados;
* a simples modificação destas variáveis troca completamente a estrutura sem a necessidade de alterar o arquio de cálculo geral de cálculo;
* se as variáeis não forem corretamente definidas e correlatas, o arquivo de cálculo geral não vai resolver a estrutura

Unidades adotadas: kN e cm


** a numeração dos nós e dos elementos pode ser qualquer!!

@autor: argenta
"""
import numpy as np
import calculoTrelica as cal
import funcoesAuxiliares as fa

# Coordenadas dos nós: dicionário contendo número do nó como chave e a tupla (coordenada X, coordenada Y) como valor
coordNos = {1:(100., 0.),
             20:(0., 100.),
             3:(0., 0.)}

# Incidência dos elementos: dicionário com o número do elemento como chave e a tupla (número do nó inicial, número do nó final) como valor
incElems = {1:(1, 20),
           2:(3, 20),
           30:(3, 1)}

# Material das barras: dicionário com o número do elemento da barra como chave e o valor do módulo de elasticidade como valor
materiais = {1:20000.,
             2:6900.,
             30:6900.}

# Seções transversais das barras: dicionário com o número do elemento da barra como chave e o valor da área como valor
area = np.pi*11.3**2/4 #cm2, todas são iguais
secoes = {1:area,
          2:area,
          30:area}

# Cargas aplicadas aos nós na estrutura: dicionário contendo o número do nó onde existe carga aplicada e uma tupla (valor carga X, valor carga Y) como valor
# * somente precisa conter os nós com carga
cargas = {1:(0., -10.)}

# Apoios aplicados aos nós na estrutura: dicionário contendo o número do nó onde existe o apoio aplicada e uma tupla (condição em X, condição em Y) como valor
# Se a condição em uma direção for igual a 1, significa restringido, se for igual a zero, significa lire
# * somente precisa conter os nós com apoio
apoios = {20:(1, 0),
          3:(1, 1)}

# # Gerando a visualização da malha com cargas e apoios
# fa.visual_TP(coordNos, incElems, cargas, apoios)


# Resolvendo a estrutura
resultados = cal.calculoTrelicaPlana(coordNos, incElems, materiais, 
                                     secoes, cargas, apoios)

# Abrindo os resultados
deslocamentos, reacoes, deformacoes, tensoes, normais = resultados


# # Gerando a visualização dos resultados na malha
# fa.visual_TP(coordNos, incElems, cargas, apoios, deslocamentos=deslocamentos)
# fa.visual_TP(coordNos, incElems, cargas, apoios, deformacoes=deformacoes)
# fa.visual_TP(coordNos, incElems, cargas, apoios, tensoes=tensoes)
# fa.visual_TP(coordNos, incElems, cargas, apoios, normais=normais)