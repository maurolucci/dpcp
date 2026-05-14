# Convertir instancias de CFC a DPCP

# Tomar primero el nombre del archivo de la línea de comandos
import argparse

parser = argparse.ArgumentParser(description="Procesa archivos CFC.")
parser.add_argument("nombre_base", help="Ruta base del archivo sin extensión .cfc")
args = parser.parse_args()

nombre_base = args.nombre_base
archivo = nombre_base + ".cfc"

# Lectura de la instancia CFC
try:
    f = open(archivo)
    s = f.readline()
    L = s.split(' ')
    n = int(L[0])
    m = int(L[1])
    V = list(range(n))
    EE = []
    for i in range(m):
        s= f.readline()
        if s=='':
            break
        L = s.split(' ')
        EE.append([int(L[i]) for i in range(1,len(L))]) # Ignorar el primer elemento que es el tamaño de la hiperarista
except:
    print(f'Error en el formato del archivo {archivo}')
else:
    print(f'Archivo {archivo} leído con éxito.')
    # print(f'{n}:{m}')
    # for i in range(n):
    #    for j in G[i]:
    #        print(f'{i},{j}')


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
nP = m
P = [[i for i in V2 if V2dict[i][0]==pi] for pi in range(nP)]
nQ = n
Q = [[i for i in V2 if V2dict[i][1]==qj] for qj in range(nQ)]
# Verificar que son particiones
for i in V2:
    assert(len([L for L in P if i in L])==1) 
    assert(len([L for L in Q if i in L])==1) 
    
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
