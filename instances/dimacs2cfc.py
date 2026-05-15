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
                if j not in G[i]:
                    G[i].append(j)
                if i not in G[j]:
                    G[j].append(i)
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
    EE[i].sort()  # ordenar los vertices de las hiperaristas

# Remover hiperaristas vacías o repetidas
EE = [e for e in EE if len(e) > 0]  # eliminar hiperaristas vacías
EE = list(set(tuple(e) for e in EE))  # eliminar hiperaristas repetidas

# Ajustar el valor de m al número de hiperaristas
m = len(EE)

# Escritura de la instancia CFC en archivo
archivo = ruta_salida / (nombre_base + ".cfc")
f = open(archivo, "w")
f.write(f'{n} {m}\n')
for e in EE:
    L = [str(i) for i in e]
    s = f'{len(e)} ' + ' '.join(L) + '\n'
    f.write(s)
f.close()
print(f'Archivo {archivo} escrito con éxito.')  
