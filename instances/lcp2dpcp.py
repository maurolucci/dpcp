# Convertir instancias de LCP a DPCP

# Lectura de la instancia LCP
nombre_base = 'infeasibles/K333'

# Lectura del grafo
archivo = nombre_base + '.graph'

try:
    f = open(archivo)
    s = f.readline()
    L = s.split(':')
    n = int(L[0])
    m = int(L[1])
    V = list(range(n))
    G = [[] for i in V]
    while True:
        s= f.readline()
        if s=='':
            break
        L = s.split(',')
        i = int(L[0])
        j = int(L[1])
        G[i].append(j)
except:
    print(f'Error en el formato del archivo {archivo}')
else:
    print(f'Archivo {archivo} leído con éxito.')
    # print(f'{n}:{m}')
    # for i in range(n):
    #    for j in G[i]:
    #        print(f'{i},{j}')

# Lectura del archivo de listas de colores
archivo = nombre_base + '.list'
try:
    f = open(archivo)
    s = f.readline()
    L = s.split(':')
    n2 = int(L[0])
    c = int(L[1])
    assert(n==n2)
    Lv = [[] for i in range(n)]
    i = 0 # índice del nodo actual
    while True:
        s= f.readline()
        if s=='':
            break
        L = s[:-1].split(' ')
        # print(L)
        assert(int(L[0])==(len(L)-1))
        for j in range(len(L)):
            if j>0:
                Lv[i].append(int(L[j]))
        i += 1
except:
    print(f'Error en el formato del archivo {archivo}')
else:
    print(f'Archivo {archivo} leído con éxito.')
    # print(f'{n}:{c}')
    # for i in range(n):
    #     L = [len(Lv[i])] + [k for k in Lv[i]]
    #     print(L)

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
nA = n
PA = [[i for i in V2 if V2dict[i][0]==a] for a in V]
nB = c
PB = [[i for i in V2 if V2dict[i][1]==b] for b in range(c)]
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
archivo_partA = nombre_base + '.partA'
f = open(archivo_partA, "w")
f.write(f'{n2}:{nA}\n')
for a in range(len(PA)):
    L = [str(i) for i in PA[a]]
    s = f'{a} ' + f'{len(PA[a])} ' + ' '.join(L) + '\n'
    f.write(s)

# Archivo de segunda partición
archivo_partB = nombre_base + '.partB'
f = open(archivo_partB, "w")
f.write(f'{n2}:{nB}\n')
for b in range(len(PB)):
    L = [str(i) for i in PB[b]]
    s = f'{b} ' + f'{len(PB[b])} ' + ' '.join(L) + '\n'
    f.write(s)




         

