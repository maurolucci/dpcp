# -*- coding: utf-8 -*-
# Crear instancias infactibles de DPCP
# A partir de instancias de LCP sobre grafos multipartitos completos Kn,n,...,n

# Parámetros del grafo Kn,n,...,n
n = 4  # número de nodos en cada partición

V = list(range(n*(n-1)))
G = [[] for i in V]
Lv = [[] for i in V] # listas de colores permitidas

# Construcción del grafo Kn,n,...,n
# Cada partición tiene n nodos, y hay (n-1) particiones
# Los nodos se numeran de 0 a n*(n-1)-1
# Los nodos de la partición k son k*n, k*n+1, ..., k*n+n-1
for k in range(n-1):
    for i in range(n):
        u = k*n + i # índice del nodo i en la partición k
        # Conectar con todos los nodos que no estén en la partición k
        for k2 in range(k+1, n-1):
            for i2 in range(n):
                v = k2*n + i2
                assert v>u
                G[u].append(v)
        # La lista de colores permitidas para el nodo u es {0,1,...,n-1}\{i}
        Lv[u] = [j for j in range(n) if j!=i]
m = sum(len(L) for L in G)

# Creación de la instancia de DPCP
# Nodos
r = 0
V2dict = dict()
for i in range(len(Lv)):
    for k in Lv[i]:
        V2dict[r] = (i,k)
        r += 1

# print(V2dict)
V2 = list(V2dict.keys())
n2 = len(V2)
# print(V2)

# Aristas
G2 = [[] for i in V2]
for i in V2:
    for j in V2:
        if j<=i:
            continue
        (u1, k1) = V2dict[i]
        (u2, k2) = V2dict[j]
        # Agregar arista en G2 si k1!=k2 o u1,u2 in E 
        if (not k1==k2) or (u2 in G[u1]) or (u1 in G[u2]):
            G2[i].append(j)
m2 = sum(len(L) for L in G2)

# Particiones
nA = len(V) # número de nodos originales
PA = [[i for i in V2 if V2dict[i][0]==a] for a in V]
nB = n # número de colores
PB = [[i for i in V2 if V2dict[i][1]==b] for b in range(n)]
# Verificar que son particiones
for i in V2:
    assert(len([L for L in PA if i in L])==1) 
    assert(len([L for L in PB if i in L])==1) 

# Exportar archivos
nombre_base = f'infeasibles/K{n}'

# Archivo del grafo multipartito completo Kn,n,...,n
archivo = nombre_base + '.graph'
f = open(archivo, "w")
f.write(f'{len(V)}:{m}\n')
for i in V:
    for j in G[i]:
        f.write(f'{i},{j}\n')

# Archivo con listas de colores para los nodos
archivo = nombre_base + '.list'
f = open(archivo, "w")
f.write(f'{len(V)}:{n}\n')
for i in V:
    L = [str(c) for c in Lv[i]]
    s = f'{len(L)} ' + ' '.join(L) + '\n'
    f.write(s)

# Archivo de la instancia de DPCP
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




         

