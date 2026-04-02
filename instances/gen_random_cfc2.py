# Crear instancias aleatorias de CFC y convertirlas a DPCP
# Se definen instancias de CFC a partir de hipergrafos aleatorios
# Para cada hipergrafo, se especifica el número de vértices n, el número de hiperaristas m,
# el tamaño promedio de las hiperaristas mu, y la desviación estándar sigma.
# Se establecen los tamaños de las hiperaristas distribuidos normalmente con media mu y desviación sigma.
# Para cada hiperarista, se seleccionan aleatoriamente los vértices que la componen (se ordenan aleatoriamente y se toman los primeros)
# A partir de la instancia CFC, se genera la instancia DPCP correspondiente.

# Parámetros de la línea de comandos:
# n : número de vértices del hipergrafo
# m : número de hiperaristas
# mu : tamaño promedio de las hiperaristas
# sigma : desviación estándar del tamaño de las hiperaristas
# path : ruta hacia la carpeta de salida
# Ejemplo de uso:
# python3 aleatorio_cfc1.py 50 100 10 2 ./instancias


import argparse
import random   

# Función para generar una instancia aleatoria de CFC
def random_cfc_2(n, m, mu, sigma):
    V = list(range(n))
    EE = []
    for k in range(m):
        t = int(round(random.normalvariate(mu, sigma)))
        t = max(1, min(n, t))  # Asegurar que 1 <= t <= n
        e = random.sample(V, t)
        EE.append(e)
    return V, EE

parser = argparse.ArgumentParser(description="Procesa archivos PCP.")
parser.add_argument("n", help="Número de vértices del hipergrafo", type=int)
parser.add_argument("m", help="Número de aristas del hipergrafo", type=int)
parser.add_argument("mu", help="Tamaño promedio de una hiperarista", type=float)
parser.add_argument("sigma", help="Desviación estándar del tamaño de una hiperarista", type=float)
parser.add_argument("path", help="Ruta hacia la carpeta de salida")
args = parser.parse_args()

n = args.n
m = args.m
mu = args.mu
sigma = args.sigma
path = args.path
nombre_base = f'{path}/n{n}_m{m}_mu{mu}_sig{sigma}'
archivo = nombre_base + ".cfc"

# Generación de la instancia CFC
V, EE = random_cfc_2(n, m, mu, sigma)

# Escritura de la instancia CFC en archivo
f = open(archivo, "w")
f.write(f'{n} {m}\n')
for e in EE:
    L = [str(i) for i in e]
    s = f'{len(e)} ' + ' '.join(L) + '\n'
    f.write(s)
print(f'Archivo {archivo} escrito con éxito.')  

# Creación de la instancia de DPCP

# Nodos
# Los nodos son pares (E,v) donde E es una hiperarista y v es un vértice en E
# Representamos los nodos como enteros desde 0 en la lista V2
# Almacenamos el par (E,v) en un diccionario V2dict 
V2dict = {}
for e in EE:
    for v in e:
        V2dict[len(V2dict)] = (EE.index(e), v)  # (E,v)
V2 = list(V2dict.keys())
n2 = len(V2)    

# Aristas
# Hay una arista entre i1=(E1,v1) y i2=(E2,v2) si i1<i2, v1!=v2, y ((v1 in E2) or (v2 in E1))
G2 = {i:[] for i in V2}
m2 = 0
for i1 in V2:
    E1, v1 = V2dict[i1]
    for i2 in V2:
        if i2 <= i1:
            continue
        E2, v2 = V2dict[i2]
        if v1 != v2 and (v1 in EE[E2] or v2 in EE[E1]):
            G2[i1].append(i2)
            m2 += 1


# Particiones
nA = m
PA = [[i for i in V2 if V2dict[i][0]==a] for a in range(nA)]
nB = n
PB = [[i for i in V2 if V2dict[i][1]==b] for b in range(nB)]
# Verificar que son particiones
for i in V2:
    assert(len([L for L in PA if i in L])==1) 
    assert(len([L for L in PB if i in L])==1) 
    
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