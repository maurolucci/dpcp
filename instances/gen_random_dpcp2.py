# Crear instancias aleatorias de DPCP
# Se define un grafo aleatorio a partir del número de vértices y la densidad de aristas.
# Se especifican las cardinalidades de la primera y segunda partición.
# Se asignan aleatoriamente pares (pi,qj) a un vértice del grafo, asegurando que cada elemento de ambas particiones (es decir, cada pi y cada qj) aparezca al menos una vez. 
# Se asegura además que no se repitan las etiquetas (pi,qj) en los vértices del grafo.

# Parámetros de la línea de comandos:
# n : número de vértices del grafo
# p : densidad de aristas del grafo
# nP : cardinalidad de la primera partición
# nQ : cardinalidad de la segunda partición
# path : ruta hacia la carpeta de salida
# Ejemplo de uso:
# python3 aleatorio_dpcp2.py 300 0.7 30 20 ./aleatorias
# se sugiere que n < nP * nQ

import argparse
import random   

parser = argparse.ArgumentParser(description="Genera instancias aleatorias de DPCP.")
parser.add_argument("n", help="Número de vértices del grafo", type=int)
parser.add_argument("p", help="Densidad de aristas", type=float)
parser.add_argument("nP", help="Cardinalidad de la primera partición", type=int)
parser.add_argument("nQ", help="Cardinalidad de la segunda partición", type=int)
parser.add_argument("path", help="Ruta hacia la carpeta de salida")
args = parser.parse_args()

n = args.n
p = args.p
nP = args.nP
nQ = args.nQ

path = args.path
nombre_base = f'{path}/n{n}_p{p}_nA{nP}_nB{nQ}_unique'

# Creación de la instancia de DPCP

# Nodos
# Representamos los nodos como enteros desde 0 en la lista V2
V2 = list(range(n))
n2 = len(V2)    

# Aristas
# Se generan aleatoriamente las aristas del grafo con densidad p
G2 = {i:[] for i in V2}
m2 = 0
for i1 in V2:
    for i2 in V2:
        if i2 <= i1:
            continue
        if random.random() < p:
            G2[i1].append(i2)
            m2 += 1

# Etiquetas (pi,qj)
# Se crea una lista con todos los pares (pi,qj) posibles
Ltodos = [(pi, qj) for a in range(nP) for b in range(nQ)]
# Se crea una lista La con todos los elementos de la primera partición
La = list(range(nP))
# Se procede de forma similar para la segunda partición
Lb = list(range(nQ))
# Si nP < nQ, se añaden elementos aleatorios a La hasta completar nQ
if nP < nQ:
    La += [random.randint(0, nP - 1) for _ in range(nQ - nP)]
    random.shuffle(La) 
# Si nQ < nP, se añaden elementos aleatorios a Lb hasta completar nP
elif nQ < nP:
    Lb += [random.randint(0, nQ - 1) for _ in range(nP - nQ)]
    random.shuffle(Lb)
# Ahora ambos tienen la misma longitud, se combinan
Lab = list(zip(La, Lb))
# Se eliminan de Ltodos los pares ya usados en Lab
for par in Lab:
    Ltodos.remove(par)
# Se extraen aleatoriamente n2 - len(Lab) pares de Ltodos y se añaden a Lab utilizando sample 
if len(Lab) < n2:
    Lab += random.sample(Ltodos, n2 - len(Lab))

# Se asignan las etiquetas a los vértices
V2dict = {i : Lab[i] for i in V2}  # (pi,qj)

# Particiones
P = [[i for i in V2 if V2dict[i][0]==pi] for pi in range(nP)]
Q = [[i for i in V2 if V2dict[i][1]==qj] for qj in range(nQ)]
# Verificar que son particiones
for i in V2:
    assert(len([L for L in P if i in L])==1) 
    assert(len([L for L in Q if i in L])==1) 
# Verificar que cada a en la primera partición aparece al menos una vez
for pi in range(nP):
    assert(len(P[pi]) > 0)
# Verificar que cada b en la segunda partición aparece al menos una vez
for qj in range(nQ):
    assert(len(Q[qj]) > 0)
# Verificar si las etiquetas (pi,qj) se repiten
etiquetas = list(V2dict.values())
etiquetas_set = set(etiquetas)   
assert len(etiquetas) == len(etiquetas_set), "Error: Hay etiquetas (pi,qj) que se repiten en varios vértices."
    
# Exportar archivos
# Archivo del grafo
nombre_base += '.dpcp'
archivo_grafo = nombre_base + '.graph'
f = open(archivo_grafo, "w")
f.write(f'{n2}:{m2}\n')
for i in V2:
    for j in G2[i]:
        f.write(f'{i} {j}\n')

# Archivo de claves de nodos
archivo_diccionario = nombre_base + '.dict'
f = open(archivo_diccionario, "w")
f.write(f'{n2}:{nP}:{nQ}\n')
for i in V2:
    f.write(f'{i} {V2dict[i][0]} {V2dict[i][1]}\n')

# Archivo de primera partición
archivo_partP = nombre_base + '.partP'
f = open(archivo_partP, "w")
f.write(f'{n2}:{nP}\n')
for pi in range(len(P)):
    L = [str(i) for i in P[pi]]
    s = f'{pi} ' + f'{len(P[pi])} ' + ' '.join(L) + '\n'
    f.write(s)

# Archivo de segunda partición
archivo_partQ = nombre_base + '.partQ'
f = open(archivo_partQ, "w")
f.write(f'{n2}:{nQ}\n')
for qj in range(len(Q)):
    L = [str(i) for i in Q[qj]]
    s = f'{qj} ' + f'{len(Q[qj])} ' + ' '.join(L) + '\n'
    f.write(s)

print(f'Archivos {archivo_grafo}, {archivo_diccionario}, {archivo_partP} y {archivo_partQ} escritos con éxito.')