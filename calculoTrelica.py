#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Criado na Segunda, dia 23 de Março de 2020, às 07:12:41

Módulo para a implementação dos procedimentos de cálculo de elementos finitos de treliça plana lineares.

@autor: argenta
"""
import numpy as np
import funcoesAuxiliares as fa

import sympy as sp #<-- para mostrar os valores para o texto!!!



# Definição da função de solução
def calculoTrelicaPlana(coordNos, incElems, materiais, secoes, cargas, apoios):
    '''
    Função para a solução de quaisquer treliças planas lineares pelo método dos
    elementos finitos conforme os argumentos que são os dados de entrada.
    
    Entrada
    -------
        * coordenadas dos nós: coordNos
        * incidência dos elementos: incElems
        * materiais das barras: materiais
        * seções transversais das barras: secoes
        * cargas dos nós da estrutura: cargas
        * apoios dos nós da estrutura: apoios
    
    Saída
    -----
        * deslocamentos da estrutura: deslocamentos
        * reações de apoio da estrutura: RXY
        * Deformações dos elementos da estrutura: deformacoes
        * tensões nos elementos da estrutura: tensoes
        * esforços normais nos elementos da estrutura: normais
    
    '''
    ### MALHA DE ELEMENTOS FINITOS ------------------------------------------------
    # Determinação dos comprimentos, cossenos, senos e matrizes de rigidez dos 
    # elementos e armazenamendo em dicionários
    
    # Criação dos dicionários vazios
    comps = {} #comprimentos
    sens = {} #senos
    coss = {} #cossenos
    kegs = {} #matrizes de rigidez
    
    # Correndo um laço nos elementos
    for elem in incElems:
        #buscando as coordenadas dos nós inicial e final
        xi, yi = coordNos[incElems[elem][0]]
        xj, yj = coordNos[incElems[elem][1]]
        
        #calculando o comprimento, seno e cosseno
        comp, sen, cos = fa.compSenCos(xi, xj, yi, yj)
        
        #armazenando os valores calculados para o respectivo elemento
        comps[elem] = comp
        sens[elem] = sen
        coss[elem] = cos
        
        #levantando os dados de material e seção do elemento
        mat = materiais[elem] #módulo de elasticidade
        sec = secoes[elem] #área da seção transversal
        
        #calculado a matriz de rigidez do elemento
        keg = fa.matRig_TP(mat, sec, comp, sen, cos)
        
        #armazenando a matriz de rigidez do respectivo elemento
        kegs[elem] = keg
    
    
    ### ELEMENTOS FINITOS DA ESTRUTURA --------------------------------------------
    # Criação automática dos indexadores para cada elemento, montagem da matriz
    # de rigidez da estrutura, do vetor de forças nodais da estrutura e a
    # criação dos graus de liberdade livres e restringidos
    
    # Determinando a quantidade de graus de liberdade totais da estrutura
    nGLs = len(coordNos)*2 #2 graus de liberdade por nó da estrutura
    
    # Criando um dicionário vazio para o armazenamento dos indexadores
    indis = {}
    
    # Inicializando o a matriz de rigidez e o vetor de forças nodais com zeros
    Kest = np.zeros((nGLs, nGLs))
    Fest = np.zeros(nGLs)
    
    # Montagem da matriz da estrutura com as matrizes dos elementos
    # Criação da listagem dos nós para a obtenção das posições
    nosList = list(coordNos.keys())
    for elem in incElems:
        #busca dos nós inicial e final
        noi, noj = incElems[elem]
        
        #definição das posições dos nós: deve iniciar em 1
        posi = nosList.index(noi) + 1
        posj = nosList.index(noj) + 1
        
        #criação dos indexadores
        indi = np.array([2*posi - 1,
                         2*posi,
                         2*posj - 1,
                         2*posj]) - 1 #já pythonizado: inicia em 0
        #armazenando o indexados para o respectivo elemento
        indis[elem] = indi
        
        #montagem da matriz de rigidez da estrutura
        for lin in range(4): #4 linhas na matriz de rigidez de cada elemento
            for col in range(4): #4 colunas na matriz de rigidez de cada elemento
                #gerando a matriz de rigidez da estrutura
                Kest[ indi[lin], indi[col] ] += kegs[elem][lin, col]
    
    # #valores para o texto:
    # print(sp.latex(sp.Matrix(np.around(Kest, 3))).replace('.', ','))
    
    
    # Montagem do etor de forças nodais da estrutura
    for no in cargas:
        #posicionamento correto das cargas já pythonizado no índice: inicia em 0
        Fest[ 2*(nosList.index(no) + 1) - 2 ] += cargas[no][0]
        Fest[ 2*(nosList.index(no) + 1) - 1 ] += cargas[no][1]
    
    # Definição dos graus de liberdade livres e restringidos já pythonizados: inicia em zero
    # Criação de listas vazias para o armazenamento de valores
    GLslivr = []
    GLsrest = []
    
    #correndo em todos os nós
    for no in coordNos:
        #buscando se existe algum apoio no nó e armazenando se o grau de liberdade
        #restringido ou livre
        if no in apoios.keys():
            restX = apoios[no][0]
            restY = apoios[no][1]
            if restX == 1:
                GLsrest.append(2*(nosList.index(no) + 1) - 2)
            elif restX == 0:
                GLslivr.append(2*(nosList.index(no) + 1) - 2)
            else:
                raise ValueError('Somente graus de liberdade livres ou ' +\
                                 'restringidos em X! Verifique as entradas em apoio.')
            if restY == 1:
                GLsrest.append(2*(nosList.index(no) + 1) - 1)
            elif restY == 0:
                GLslivr.append(2*(nosList.index(no) + 1) - 1)
            else:
                raise ValueError('Somente graus de liberdade livres ou ' +\
                                 'restringidos em Y! Verifique as entradas em apoio.')
        #caso não exista armazena os graus de liberdade livres do nó
        else:
            GLslivr.append(2*(nosList.index(no) + 1) - 2)
            GLslivr.append(2*(nosList.index(no) + 1) - 1)
    
    ### RESOLUÇÃO E SEPARAÇÃO -----------------------------------------------------
    # Cálculo dos deslocamentos, reações de apoio e separação da solução nos elementos
    
    # Separação da matriz de rigidez para o cálculo dos deslocamentos das reações 
    # de apoio
    Ku = Kest[:, GLslivr]
    
    # #valores para o texto:
    # print(sp.latex(sp.Matrix(np.around(Ku, 3))).replace('.', ','))
    
    Ku = Ku[GLslivr, :]
    
    # #valores para o texto:
    # print(sp.latex(sp.Matrix(np.around(Ku, 3))).replace('.', ','))
    
    Kr = Kest[:, GLslivr]
    Kr = Kr[GLsrest, :]
    
    # #valores para o texto:
    # print(sp.latex(sp.Matrix(np.around(Kr, 3))).replace('.', ','))
    
    # Separação do vetor de forças nodais para cálculo dos deslocamentos das reações 
    # de apoio
    Fu = Fest[GLslivr]
    Fr = Fest[GLsrest]
    
    # Determinação dos deslocamentos
    Us = np.linalg.solve(Ku, Fu)
    
    # Determinação das reações de apoio
    Re = np.matmul(Kr, Us) - Fr
    
    # Montagem do vetor de deslocamentos completo: com os deslocamentos iguais a zero
    Ug = np.zeros(nGLs)
    Ug[GLslivr] = Us
    
    # Montagem do vetor de reações completo: com as reações iguais a zero
    Rg = np.zeros(nGLs)
    Rg[GLsrest] = Re
    
    # Convertendo os deslocamentos para o formato (deslocamento X, deslocamento U) e
    # associando ao respectivo nó, as reações para (reação X, reação Y) e associando
    # ao respectivo nó e criando os vetores de deslocamentos dos elementos usando 
    # os indexadores e a matriz de decomposição iniciando com um dicionários vazios
    # para armazenamento
    UXY = {}
    RXY = {}
    ues = {}
    for elem in incElems:
        #separando os deslocamentos no elemento no sistema global
        ug = Ug[indis[elem]] #deslocamentos no elemento
        
        #separando as reações no elemento no sistema global
        rg = Rg[indis[elem]] #reações no elemento
        
        #armazenando os deslocamentos do elemento nos respetivos graus de liberdade
        # X, Y dos respectivos nós
        UXY[incElems[elem][0]] = (ug[0], ug[1]) #nó inicial
        UXY[incElems[elem][1]] = (ug[2], ug[3]) #nó final
        
        #armazenando as reações do elemento nos respetivos graus de liberdade
        # X, Y dos respectivos nós
        RXY[incElems[elem][0]] = (rg[0], rg[1]) #nó inicial
        RXY[incElems[elem][1]] = (rg[2], rg[3]) #nó final
        
        #conversão dos deslocamentos para o sistema local com a matriz de 
        #decomposição transposta
        DecT = np.array([[coss[elem], sens[elem], 0, 0],
                         [0, 0, coss[elem], sens[elem]]])
        
        # print(DecT) #--> para conferência!
        
        #armazenando os deslocamentos no elementodeslocamentos
        ues[elem] = np.matmul(DecT, ug)
    
    
    ### RESULTADOS FINAIS NOS ELEMENTOS -------------------------------------------
    # Determinação das deformações, tensões e esforços normais nos elementos, iniciando
    # com dicionários vazios para armazenamento
    defos = {}
    tenss = {}
    norms = {}
    
    #correndo em todos os elementos
    for elem in incElems:
        #criando a matriz das derivadas das funções de interpolação
        B = np.array([-1./comps[elem], 1./comps[elem]])
        
        #calculando e armazenando a deformação no elemento
        defos[elem] = np.matmul(B, ues[elem])
        
        #calculando e armazenando a tensão no elemento usando o módulo de elasticidade
        tenss[elem] = materiais[elem]*defos[elem]
        
        #calculando e armazenando os esforços normais no elemento usando a área
        #da seção transversal
        norms[elem] = secoes[elem]*tenss[elem]
    
    
    # Finalizando a função e retornado
    return UXY, RXY, defos, tenss, norms

# testando o módulo de cálculo da treliça
if __name__ == '__main__':
    ### DADOS DE ENTRADA
    # Coordenadas dos nós: dicionário contendo número do nó como chave e a tupla (coordenada X, coordenada Y) como valor
    coordNos = {1:(100., 0.),
                 2:(0., 100.),
                 3:(0., 0.)}

    # Incidência dos elementos: dicionário com o número do elemento como chave e a tupla (número do nó inicial, número do nó final) como valor
    incElems = {1:(1, 2),
               2:(3, 2),
               3:(3, 1)}

    # Material das barras: dicionário com o número do elemento da barra como chave e o valor do módulo de elasticidade como valor
    materiais = {1:20000.,
                 2:6900.,
                 3:6900.}

    # Seções transversais das barras: dicionário com o número do elemento da barra como chave e o valor da área como valor
    area = np.pi*11.3**2/4 #cm2, todas são iguais
    secoes = {1:area,
              2:area,
              3:area}

    # Cargas aplicadas aos nós na estrutura: dicionário contendo o número do nó onde existe carga aplicada e uma tupla (valor carga X, valor carga Y) como valor
    # * somente precisa conter os nós com carga
    cargas = {1:(0., -10.)}

    # Apoios aplicados aos nós na estrutura: dicionário contendo o número do nó onde existe o apoio aplicada e uma tupla (condição em X, condição em Y) como valor
    # Se a condição em uma direção for igual a 1, significa restringido, se for igual a zero, significa lire
    # * somente precisa conter os nós com apoio
    apoios = {2:(1, 0),
              3:(1, 1)}
    
    deslocamentos, RXY, deformacoes, tensoes, normais = calculoTrelicaPlana(coordNos, incElems, materiais, secoes, cargas, apoios)



