# Convertir instancias de VCP a DPCP

# Lee instancias de VCP en formato .col de DIMAC y las convierte a instancias
# de CFC sobre el hipergrafo de vecindades abiertas o cerradas.
# Escribe las instancias en formato .cfc y las exporta luego a DPCP,
# generando los archivos .graph, .dict, .partP y .partQ

# Ejemplo de uso:
# python dimacs2dpcp.py nombre_base ruta_salida [--cerrada]
# donde nombre_base es la ruta del archivo .col
# y ruta_salida es la ruta donde se guardarán los archivos de salida.   
# La opción --cerrada indica que se use la vecindad cerrada en lugar de la abierta.
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="Procesa archivos en formato COL.")
parser.add_argument("ruta_completa_base", help="Ruta base del archivo de entrada")
parser.add_argument("ruta_salida", help="Ruta para los archivos de salida")
parser.add_argument("--cerrada", help="Usar vecindad cerrada en lugar de abierta", action="store_true")

args = parser.parse_args()

ruta_completa_base = Path(args.ruta_completa_base)
nombre_base = ruta_completa_base.stem
ruta_base = ruta_completa_base.parent
ruta_salida = Path(args.ruta_salida)
vecindad_cerrada = args.cerrada

# Crear las carpetas de salida si no existen
ruta_salida.mkdir(parents=True, exist_ok=True)

# Lectura de la instancia en formato .col
archivo = ruta_base / (nombre_base + ".col")
try:
    with open(archivo, "r", encoding="utf-8") as f:
        while True:
            s = f.readline()
            if not s:  # condición de parada
                break
            L = s.split(' ')
            if L[0] == 'p' and L[1] == 'edge':
                # Inicializar grafo
                n = int(L[2])
                m = int(L[3])
                V = list(range(n))
                G = [[] for i in V]
            elif L[0] == 'e':
                i = int(L[1]) - 1  # ajustar índice a partir de 0
                j = int(L[2]) - 1
                G[i].append(j)
except:
    print(f'Error al abrir el archivo {archivo}')
    exit(1)
else:
    print(f'Archivo {archivo} leído con éxito.')
    # print(f'{n}:{m}')
    # for i in range(n):
    #    for j in G[i]:
    #        print(f'{i},{j}')

# Creación del hipergrafo de vecindades
EE = [[] for i in V]
for i in V:
    EE[i] = G[i].copy()
    if vecindad_cerrada:
        EE[i].append(i)  # usar vecindades cerradas
    EE[i].sort()  # ordenar las hiperaristas
    
# Escritura de la instancia CFC en archivo
archivo = ruta_salida / (nombre_base + ".cfc")
f = open(archivo, "w")
f.write(f'{n} {n}\n')
for e in EE:
    L = [str(i) for i in e]
    s = f'{len(e)} ' + ' '.join(L) + '\n'
    f.write(s)
print(f'Archivo {archivo} escrito con éxito.')  

# Creación de la instancia de DPCP

# Nodos
# Los nodos son pares (e,v) donde e es una hiperarista y v es un vértice en E
# Representamos los nodos como enteros desde 0 en la lista V2
# y las hiperaristas como los índices de la lista EE
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
nA = n
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
archivo_grafo = ruta_salida / (nombre_base + '.graph')
f = open(archivo_grafo, "w")
f.write(f'{n2}:{m2}\n')
for i in V2:
    for j in G2[i]:
        f.write(f'{i} {j}\n')

# Archivo de etiquetas de nodos
archivo_diccionario = ruta_salida / (nombre_base + '.dict')
f = open(archivo_diccionario, "w")
f.write(f'{n2}:{nA}:{nB}\n')
for i in V2:
    f.write(f'{i} {V2dict[i][0]} {V2dict[i][1]}\n')

# Archivo de primera partición
archivo_partP = ruta_salida / (nombre_base + '.partP')
f = open(archivo_partP, "w")
f.write(f'{n2}:{nA}\n')
for a in range(len(PA)):
    L = [str(i) for i in PA[a]]
    s = f'{a} ' + f'{len(PA[a])} ' + ' '.join(L) + '\n'
    f.write(s)

# Archivo de segunda partición
archivo_partQ = ruta_salida / (nombre_base + '.partQ')
f = open(archivo_partQ, "w")
f.write(f'{n2}:{nB}\n')
for b in range(len(PB)):
    L = [str(i) for i in PB[b]]
    s = f'{b} ' + f'{len(PB[b])} ' + ' '.join(L) + '\n'
    f.write(s)

print(f'Archivos {archivo_grafo}, {archivo_diccionario}, {archivo_partP} y {archivo_partQ} escritos con éxito.')
