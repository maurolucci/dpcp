# Generador de instancias aleatorias de DPCP
#
# Parámetros de la línea de comandos:
# n : número de vértices del grafo
# p : densidad de aristas del grafo
# nP : cardinalidad de la primera partición
# nQ : cardinalidad de la segunda partición
# ninst : número de instancias a generar
# path : ruta hacia la carpeta de salida
# Ejemplo de uso:
# python3 gen_random_dpcp1.py 300 0.7 30 20 3 ./aleatorias
#
# Metodologia:
# 1. Se genera un grafo aleatorio con parámetros n y p.
# 2. A los primeros nP-1 vértices se le asigna una partición distinta de P, a los restantes una aleatoria.
# 3. A los primeros nQ-1 vértices se le asigna una partición distinta de Q, a los restantes una aleatoria.
#
# Warnings:
# 1. Se sugiere que n < nP * nQ
# 2. La generación permite que dos vertices pertenezcan a los mismos conjuntos de P y Q.

import argparse
import random

parser = argparse.ArgumentParser(description="Genera instancias aleatorias de DPCP.")
parser.add_argument("n", help="Número de vértices del grafo", type=int)
parser.add_argument("p", help="Densidad de aristas", type=float)
parser.add_argument("nP", help="Cardinalidad de la primera partición", type=int)
parser.add_argument("nQ", help="Cardinalidad de la segunda partición", type=int)
parser.add_argument("ninst", help="Número de instancias a generar", type=int)
parser.add_argument("path", help="Ruta hacia la carpeta de salida")
args = parser.parse_args()

n = args.n
p = args.p
nP = args.nP
nQ = args.nQ
ninst = args.ninst
path = args.path

# Set a fixed seed
random.seed(0)

for id in range(ninst):

    nombre_base = f"{path}/r_n{n}_p{p}_nA{nP}_nB{nQ}_i{id}"

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

    # Etiquetas (pi,qj)
    # Se crea una lista La de longitud n2 (número de vértices del grafo), con elementos de la primera partición extraidos aleatoriamente, asegurando que cada elemento de la partición aparece al menos una vez
    La = list(range(nP))
    La += [random.randint(0, nP - 1) for _ in range(n2 - nP)]
    random.shuffle(La)
    # Se procede de forma similar para la segunda partición
    Lb = list(range(nQ))
    Lb += [random.randint(0, nQ - 1) for _ in range(n2 - nQ)]
    random.shuffle(Lb)
    # Se asignan las etiquetas a los vértices
    V2dict = {i: (La[i], Lb[i]) for i in V2}  # (pi,qj)

    # Particiones
    P = [[i for i in V2 if V2dict[i][0] == pi] for pi in range(nP)]
    Q = [[i for i in V2 if V2dict[i][1] == qj] for qj in range(nQ)]
    # Verificar que son particiones
    for i in V2:
        assert len([L for L in P if i in L]) == 1
        assert len([L for L in Q if i in L]) == 1
    # Verificar que cada a en la primera partición aparece al menos una vez
    for pi in range(nP):
        assert len(P[pi]) > 0
    # Verificar que cada b en la segunda partición aparece al menos una vez
    for qj in range(nQ):
        assert len(Q[qj]) > 0
    # Verificar si las etiquetas (pi,qj) se repiten
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
    f.write(f"{n2}:{nP}:{nQ}\n")
    for i in V2:
        f.write(f"{i} {V2dict[i][0]} {V2dict[i][1]}\n")

    # Archivo de primera partición
    archivo_partP = nombre_base + ".partP"
    f = open(archivo_partP, "w")
    f.write(f"{n2}:{nP}\n")
    for pi in range(len(P)):
        L = [str(i) for i in P[pi]]
        s = f"{pi} " + f"{len(P[pi])} " + " ".join(L) + "\n"
        f.write(s)

    # Archivo de segunda partición
    archivo_partQ = nombre_base + ".partQ"
    f = open(archivo_partQ, "w")
    f.write(f"{n2}:{nQ}\n")
    for qj in range(len(Q)):
        L = [str(i) for i in Q[qj]]
        s = f"{qj} " + f"{len(Q[qj])} " + " ".join(L) + "\n"
        f.write(s)

    print(
        f"Archivos {archivo_grafo}, {archivo_diccionario}, {archivo_partP} y {archivo_partQ} escritos con éxito."
    )
