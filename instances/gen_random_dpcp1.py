# Generador de instancias aleatorias de DPCP
#
# Parámetros de la línea de comandos:
# n : número de vértices del grafo
# p : densidad de aristas del grafo
# nA : cardinalidad de la primera partición
# nB : cardinalidad de la segunda partición
# ninst : número de instancias a generar
# path : ruta hacia la carpeta de salida
# Ejemplo de uso:
# python3 gen_random_dpcp1.py 300 0.7 30 20 3 ./aleatorias
#
# Metodologia:
# 1. Se genera un grafo aleatorio con parámetros n y p.
# 2. A los primeros nA-1 vértices se le asigna una partición distinta de P, a los restantes una aleatoria.
# 3. A los primeros nB-1 vértices se le asigna una partición distinta de Q, a los restantes una aleatoria.
#
# Warnings:
# 1. Se sugiere que n < nA * nB
# 2. La generación permite que dos vertices pertenezcan a los mismos conjuntos de P y Q.

import argparse
import random

parser = argparse.ArgumentParser(description="Genera instancias aleatorias de DPCP.")
parser.add_argument("n", help="Número de vértices del grafo", type=int)
parser.add_argument("p", help="Densidad de aristas", type=float)
parser.add_argument("nA", help="Cardinalidad de la primera partición", type=int)
parser.add_argument("nB", help="Cardinalidad de la segunda partición", type=int)
parser.add_argument("ninst", help="Número de instancias a generar", type=int)
parser.add_argument("path", help="Ruta hacia la carpeta de salida")
args = parser.parse_args()

n = args.n
p = args.p
nA = args.nA
nB = args.nB
ninst = args.ninst
path = args.path

# Set a fixed seed
random.seed(0)

for id in range(ninst):

    nombre_base = f"{path}/r_n{n}_p{p}_nA{nA}_nB{nB}_i{id}"

    # Creación de la instancia de DPCP

    # Nodos
    # Representamos los nodos como enteros desde 0 en la lista V2
    V2 = list(range(n))
    n2 = len(V2)

    # Aristas
    # Se generan aleatoriamente las aristas del grafo con densidad p
    G2 = {i: [] for i in V2}
    m2 = 0
    for i1 in V2:
        for i2 in V2:
            if i2 <= i1:
                continue
            if random.random() < p:
                G2[i1].append(i2)
                m2 += 1

    # Etiquetas (a,b)
    # Se crea una lista La de longitud n2 (número de vértices del grafo), con elementos de la primera partición extraidos aleatoriamente, asegurando que cada elemento de la partición aparece al menos una vez
    La = list(range(nA))
    La += [random.randint(0, nA - 1) for _ in range(n2 - nA)]
    random.shuffle(La)
    # Se procede de forma similar para la segunda partición
    Lb = list(range(nB))
    Lb += [random.randint(0, nB - 1) for _ in range(n2 - nB)]
    random.shuffle(Lb)
    # Se asignan las etiquetas a los vértices
    V2dict = {i: (La[i], Lb[i]) for i in V2}  # (a,b)

    # Particiones
    PA = [[i for i in V2 if V2dict[i][0] == a] for a in range(nA)]
    PB = [[i for i in V2 if V2dict[i][1] == b] for b in range(nB)]
    # Verificar que son particiones
    for i in V2:
        assert len([L for L in PA if i in L]) == 1
        assert len([L for L in PB if i in L]) == 1
    # Verificar que cada a en la primera partición aparece al menos una vez
    for a in range(nA):
        assert len(PA[a]) > 0
    # Verificar que cada b en la segunda partición aparece al menos una vez
    for b in range(nB):
        assert len(PB[b]) > 0
    # Verificar si las etiquetas (a,b) se repiten
    etiquetas = set()
    for i in V2:
        etiqueta = V2dict[i]
        if etiqueta in etiquetas:
            print(f"Advertencia: La etiqueta {etiqueta} se repite en varios vértices.")
        else:
            etiquetas.add(etiqueta)

    # Exportar archivos
    # Archivo del grafo
    nombre_base += ".dpcp"
    archivo_grafo = nombre_base + ".graph"
    f = open(archivo_grafo, "w")
    f.write(f"{n2}:{m2}\n")
    for i in V2:
        for j in G2[i]:
            f.write(f"{i} {j}\n")

    # Archivo de claves de nodos
    archivo_diccionario = nombre_base + ".dict"
    f = open(archivo_diccionario, "w")
    f.write(f"{n2}:{nA}:{nB}\n")
    for i in V2:
        f.write(f"{i} {V2dict[i][0]} {V2dict[i][1]}\n")

    # Archivo de primera partición
    archivo_partA = nombre_base + ".partP"
    f = open(archivo_partA, "w")
    f.write(f"{n2}:{nA}\n")
    for a in range(len(PA)):
        L = [str(i) for i in PA[a]]
        s = f"{a} " + f"{len(PA[a])} " + " ".join(L) + "\n"
        f.write(s)

    # Archivo de segunda partición
    archivo_partB = nombre_base + ".partQ"
    f = open(archivo_partB, "w")
    f.write(f"{n2}:{nB}\n")
    for b in range(len(PB)):
        L = [str(i) for i in PB[b]]
        s = f"{b} " + f"{len(PB[b])} " + " ".join(L) + "\n"
        f.write(s)

    print(
        f"Archivos {archivo_grafo}, {archivo_diccionario}, {archivo_partA} y {archivo_partB} escritos con éxito."
    )
