# Crear instancias aleatorias de DPCP
# Se define un grafo aleatorio a partir del número de vértices y la densidad de aristas.
# Se especifican las cardinalidades de la primera y segunda partición.
# Se asignan aleatoriamente pares (a,b) a un vértice del grafo, asegurando que cada elemento de ambas particiones (es decir, cada a y cada b) aparezca al menos una vez. 
# Se asegura además que no se repitan las etiquetas (a,b) en los vértices del grafo.

# Parámetros de la línea de comandos:
# n : número de vértices del grafo
# p : densidad de aristas del grafo
# nA : cardinalidad de la primera partición
# nB : cardinalidad de la segunda partición
# path : ruta hacia la carpeta de salida
# Ejemplo de uso:
# python3 aleatorio_dpcp2.py 300 0.7 30 20 ./aleatorias
# se sugiere que n < nA * nB

import argparse
import random   

parser = argparse.ArgumentParser(description="Genera instancias aleatorias de DPCP.")
parser.add_argument("n", help="Número de vértices del grafo", type=int)
parser.add_argument("p", help="Densidad de aristas", type=float)
parser.add_argument("nA", help="Cardinalidad de la primera partición", type=int)
parser.add_argument("nB", help="Cardinalidad de la segunda partición", type=int)
parser.add_argument("path", help="Ruta hacia la carpeta de salida")
args = parser.parse_args()

n = args.n
p = args.p
nA = args.nA
nB = args.nB

path = args.path
nombre_base = f'{path}/n{n}_p{p}_nA{nA}_nB{nB}_unique'

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

# Etiquetas (a,b)
# Se crea una lista con todos los pares (a,b) posibles
Ltodos = [(a, b) for a in range(nA) for b in range(nB)]
# Se crea una lista La con todos los elementos de la primera partición
La = list(range(nA))
# Se procede de forma similar para la segunda partición
Lb = list(range(nB))
# Si nA < nB, se añaden elementos aleatorios a La hasta completar nB
if nA < nB:
    La += [random.randint(0, nA - 1) for _ in range(nB - nA)]
    random.shuffle(La) 
# Si nB < nA, se añaden elementos aleatorios a Lb hasta completar nA
elif nB < nA:
    Lb += [random.randint(0, nB - 1) for _ in range(nA - nB)]
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
V2dict = {i : Lab[i] for i in V2}  # (a,b)

# Particiones
PA = [[i for i in V2 if V2dict[i][0]==a] for a in range(nA)]
PB = [[i for i in V2 if V2dict[i][1]==b] for b in range(nB)]
# Verificar que son particiones
for i in V2:
    assert(len([L for L in PA if i in L])==1) 
    assert(len([L for L in PB if i in L])==1) 
# Verificar que cada a en la primera partición aparece al menos una vez
for a in range(nA):
    assert(len(PA[a]) > 0)
# Verificar que cada b en la segunda partición aparece al menos una vez
for b in range(nB):
    assert(len(PB[b]) > 0)
# Verificar si las etiquetas (a,b) se repiten
etiquetas = list(V2dict.values())
etiquetas_set = set(etiquetas)   
assert len(etiquetas) == len(etiquetas_set), "Error: Hay etiquetas (a,b) que se repiten en varios vértices."
    
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
f.write(f'{n2}:{nA}:{nB}\n')
for i in V2:
    f.write(f'{i} {V2dict[i][0]} {V2dict[i][1]}\n')

# Archivo de primera partición
archivo_partA = nombre_base + '.partP'
f = open(archivo_partA, "w")
f.write(f'{n2}:{nA}\n')
for a in range(len(PA)):
    L = [str(i) for i in PA[a]]
    s = f'{a} ' + f'{len(PA[a])} ' + ' '.join(L) + '\n'
    f.write(s)

# Archivo de segunda partición
archivo_partB = nombre_base + '.partQ'
f = open(archivo_partB, "w")
f.write(f'{n2}:{nB}\n')
for b in range(len(PB)):
    L = [str(i) for i in PB[b]]
    s = f'{b} ' + f'{len(PB[b])} ' + ' '.join(L) + '\n'
    f.write(s)

print(f'Archivos {archivo_grafo}, {archivo_diccionario}, {archivo_partA} y {archivo_partB} escritos con éxito.')