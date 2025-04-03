#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Criado na Quinta, dia 20 de Março de 2020, às 11:17:21

Módulo para as funções auxiliares:
    - compSenCos(x1, x2, y1, y2): cálculo do comprimento, seno e cosseno de elementos com 2 nós;
    - matRig_TP(E, A, L, s, c): cálculo da matriz de rigidez do elemento de treliça plano;

@autor: argenta
"""
import numpy as np
import matplotlib.pyplot as plt

import sympy as sp #<-- para mostrar os valores para o texto!!!

### FUNÇÕES AUXILIARES DE TRELIÇAS PLANAS -------------------------------------

def compSenCos(xi, xj, yi, yj):
    '''
    Função para calcular o comprimento de elementos, o seno e o cosseno do ângulo
    do elemento na estrutura, entre os eixos X e r no sentido anti-horário.
    
    Entrada
    -------
        * xi: coordenada x do nó inicial do elemento na estrutura;
        * xj: coordenada x do nó final do elemento na estrutura;
        * yi: coordenada y do nó inicial do elemento na estrutura;
        * yj: coordenada y do nó final do elemento na estrutura.
    
    Saída
    -----
        * comp: comprimento total do elemento;
        * sen: seno do ângulo do eixo local longitudinal do elemento r com o eixo x;
        * cos: cosseno do ângulo do eixo local longitudinal do elemento r com o eixo x;
    
    '''
    comp = np.sqrt((xj - xi)**2 + (yj - yi)**2)
    sen = (yj - yi)/comp
    cos = (xj - xi)/comp
    return comp, sen, cos

def matRig_TP(E, A, L, s, c):
    '''
    Função para calcular a matriz de rigidez global de um elemento de treliça plana.
    
    Entradas
    --------
        * E: módulo de elasticidade do material do elemento;
        * A: área da seção transversal da barra do elemento;
        * L: comprimento do elemento (consequentemente da barra);
        * s: seno do ângulo do eixo local longitudinal do elemento r com o eixo x;
        * c: cosseno do ângulo do eixo local longitudinal do elemento r com o eixo x;
    
    Saída
    -----
        * keg: a matriz de digidez do elemento de treliça representativo da barra no sistema global da estrutura.
    
    '''
    keg = E*A/L*np.array([[c**2, c*s, -c**2, -c*s],
                          [c*s, s**2, -c*s, -s**2],
                          [-c**2, -c*s, c**2, c*s],
                          [-c*s, -s**2, c*s, s**2]])
    
    # #valores para o texto:
    # print('Matriz de rigidez do elemento só com os sen e cos:')
    # parcial = np.around(np.array([[c**2, c*s, -c**2, -c*s],
    #                       [c*s, s**2, -c*s, -s**2],
    #                       [-c**2, -c*s, c**2, c*s],
    #                       [-c*s, -s**2, c*s, s**2]]), 3)
    # print(sp.latex(sp.Matrix(parcial)).replace('.', ','))
    
    # print('Matriz de rigidez do elemento completa:')
    # print(sp.latex(sp.Matrix(np.around(keg, 3))).replace('.', ','))
    
    return keg

def visual_TP(coordNos, incElems, cargas, apoios, deslocamentos=None, deformacoes=None, tensoes=None, normais=None, 
              escala_carga=2, escala_apoio=5, escala_desloc=1000, escala_deform=1e6, escala_tensao=100, escala_normal=1):
    '''
    Função para gerar a visualização da malha da treliça, cargas, apoios e
    resultados de deslocamentos, deformações, tensões e esforços normais.
    
    * Somente deve ser chamado um resultado por vez e os outros devem permanecer None. Se chamar mais de um, somente o primeiro será mostrado.

    Entrada
    -------
        coordNos: dict, Dicionário com o número do nó como chave e tupla (x, y) como valor (coordenadas).
        incElems: dict, Dicionário com o número do elemento como chave e tupla (nó inicial, nó final) como valor.
        cargas: dict, Dicionário com o número do nó como chave e tupla (fx, fy) como valor (cargas aplicadas).
        apoios: dict, Dicionário com o número do nó como chave e tupla (rest_x, rest_y) como valor (0 = livre, 1 = restrito).
        deslocamentos: dict, opcional, Dicionário com o número do nó como chave e tupla (dx, dy) como valor (deslocamentos). Padrão é None.
        deformacoes: dict, opcional, Dicionário com o número do elemento como chave e valor da deformação. Padrão é None.
        tensoes: dict, opcional, Dicionário com o número do elemento como chave e valor da tensão. Padrão é None.
        normais: dict, opcional, Dicionário com o número do elemento como chave e valor do esforço normal. Padrão é None.
        escala_carga: float, opcional, Fator de escala para as setas de carga. Padrão é 2.
        escala_apoio: float, opcional, Fator de escala para os triângulos de apoio. Padrão é 10.
        escala_desloc: float, opcional, Fator de escala para os deslocamentos. Padrão é 1000.
        escala_deform: float, opcional, Fator de escala para os retângulos de deformações. Padrão é 1e6.
        escala_tensao: float, opcional, Fator de escala para os retângulos de tensões. Padrão é 10.
        escala_normal: float, opcional, Fator de escala para os retângulos de normais. Padrão é 1.

    Saída
    -----
        * None: Gera e exibe a visualização da treliça plana.
    '''
    # Configurando a figura e definido a transparência para visualização dos deslocamentos
    plt.figure(figsize=(19.2, 10.8), dpi=100)
    alpha = 0.3 if deslocamentos else 1.0

    # Estrutura indeformada: elementos
    for elem, (ni, nf) in incElems.items():
        p = np.array([coordNos[ni], coordNos[nf]])
        plt.plot(p[:, 0], p[:, 1], 'b-', lw=2, alpha=alpha, label='Elementos' if elem == 1 else "")
        pm = p.mean(axis=0) + np.array([5, 5])
        plt.text(pm[0], pm[1], f'E{elem}', color='b', fontsize=10, alpha=alpha)
    
    # Estrutura indeformada: nós
    for no, p in coordNos.items():
        plt.plot(p[0], p[1], 'ko', ms=8, alpha=alpha, label='Nós' if no == 1 else "")
        if not deslocamentos:
            plt.text(p[0] + 5, p[1] + 5, f'N{no}', fontsize=10, alpha=alpha)
    
    # Cargas nos nós
    for no, (fx, fy) in cargas.items():
        p = coordNos[no]
        primeiro = True #fazendo apenas uma entrada na legenda
        if fx:
            plt.arrow(p[0], p[1], fx*escala_carga, 0, head_width=3, head_length=6, fc='r', ec='r', 
                      lw=abs(fx)*0.1, alpha=alpha, label='Cargas' if primeiro else "")
            plt.text(p[0] + fx*escala_carga + 2, p[1] + 2, f'{fx}', color='r', fontsize=8, alpha=alpha)
            primeiro = False
        if fy:
            plt.arrow(p[0], p[1], 0, fy*escala_carga, head_width=3, head_length=6, fc='r', ec='r', 
                      lw=abs(fy)*0.1, alpha=alpha, label='Cargas' if primeiro else "")
            plt.text(p[0] + 2, p[1] + fy*escala_carga + 2, f'{fy}', color='r', fontsize=8, alpha=alpha)
    
    # Apoios
    for no, (rx, ry) in apoios.items():
        x, y = coordNos[no]
        b, h = escala_apoio * 0.5, escala_apoio * 0.75
        if rx:
            plt.plot([x-h, x-h], [y-b, y+b], 'c-', lw=2, alpha=alpha)
            plt.plot([x-h, x], [y-b, y], 'c-', lw=2, alpha=alpha)
            plt.plot([x-h, x], [y+b, y], 'c-', lw=2, alpha=alpha, 
                     label='Apoios' if no == list(apoios.keys())[0] and not rx else "")
        if ry:
            plt.plot([x-b, x+b], [y-h, y-h], 'c-', lw=2, alpha=alpha)
            plt.plot([x-b, x], [y-h, y], 'c-', lw=2, alpha=alpha)
            plt.plot([x+b, x], [y-h, y], 'c-', lw=2, alpha=alpha, label='Apoios' if no == list(apoios.keys())[0] else "")

    # Estrutura deformada
    if deslocamentos:
        coords_desloc = {no: np.array(coordNos[no]) + np.array(d) * escala_desloc for no, d in deslocamentos.items()}
        for elem, (ni, nf) in incElems.items():
            p = np.array([coords_desloc[ni], coords_desloc[nf]])
            plt.plot(p[:, 0], p[:, 1], 'g-', lw=2, label='Elementos Deslocados' if elem == 1 else "")
        for no, p in coords_desloc.items():
            dx, dy = deslocamentos[no]
            plt.plot(p[0], p[1], 'go', ms=8, label='Nós Deslocados' if no == 1 else "")
            plt.text(p[0] + 5, p[1] + 5, f'N{no}\n({dx:.3e}, {dy:.3e})', fontsize=8, color='g')

    # Resultados: diagramas de deformações, tensões e esforços normais
    resultado = None
    if deformacoes:
        resultado = (deformacoes, 'red', 'Deformações', escala_deform)
    elif tensoes:
        resultado = (tensoes, 'orange', 'Tensões', escala_tensao)
    elif normais:
        resultado = (normais, 'blue', 'Esforços Normais', escala_normal)

    if resultado:
        valores, cor, label, escala = resultado
        for elem, (ni, nf) in incElems.items():
            pi, pf = np.array(coordNos[ni]), np.array(coordNos[nf])
            v_barra = pf - pi
            comp = np.linalg.norm(v_barra)
            u_barra = v_barra / comp
            v_normal = np.array([-u_barra[1], u_barra[0]])
            valor = valores[elem]
            largura = abs(valor) * escala
            lado = v_normal if valor >= 0 else -v_normal
            # Definindo os vértices do retângulo
            p0 = pi if valor >= 0 else pi - lado * largura
            p1 = p0 + v_barra
            p2 = p1 + lado * largura
            p3 = p0 + lado * largura
            rect_coords = np.array([p0, p1, p2, p3, p0])
            plt.plot(rect_coords[:, 0], rect_coords[:, 1], color=cor, alpha=0.5, label=label if elem == 1 else "")
            plt.fill(rect_coords[:, 0], rect_coords[:, 1], color=cor, alpha=0.5)
            # Centro do retângulo para o texto
            pm = (p0 + p2) / 2  # Média entre cantos opostos
            plt.text(pm[0], pm[1], f'{valor:.3e}', color=cor, fontsize=8, ha='center', va='center')
    
    # Configurações finais da visualização
    plt.grid(True)
    plt.legend()
    plt.title("Visualização da Treliça Plana")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.show()
    
    return None



#------------------------------------------------------------------------------






# # Testando as funções
# if __name__ == '__main__':
#     ### DADOS DE ENTRADA
#     # Coordenadas dos nós: dicionário contendo número do nó como chave e a tupla (coordenada X, coordenada Y) como valor
#     coordNos = {1:(100., 0.),
#                  2:(0., 100.),
#                  3:(0., 0.)}

#     # Incidência dos elementos: dicionário com o número do elemento como chave e a tupla (número do nó inicial, número do nó final) como valor
#     incElems = {1:(1, 2),
#                2:(3, 2),
#                3:(3, 1)}

#     # Material das barras: dicionário com o número do elemento da barra como chave e o valor do módulo de elasticidade como valor
#     materiais = {1:20000.,
#                  2:6900.,
#                  3:6900.}

#     # Seções transversais das barras: dicionário com o número do elemento da barra como chave e o valor da área como valor
#     area = np.pi*11.3**2/4 #cm2, todas são iguais
#     secoes = {1:area,
#               2:area,
#               3:area}

#     # Cargas aplicadas aos nós na estrutura: dicionário contendo o número do nó onde existe carga aplicada e uma tupla (valor carga X, valor carga Y) como valor
#     # * somente precisa conter os nós com carga
#     cargas = {1:(0., -10.)}

#     # Apoios aplicados aos nós na estrutura: dicionário contendo o número do nó onde existe o apoio aplicada e uma tupla (condição em X, condição em Y) como valor
#     # Se a condição em uma direção for igual a 1, significa restringido, se for igual a zero, significa lire
#     # * somente precisa conter os nós com apoio
#     apoios = {2:(1, 0),
#               3:(1, 1)}
    
#     # Gerando a visualização da malha com cargas e apoios
#     visual_TP(coordNos, incElems, cargas, apoios)
