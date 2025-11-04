# cfcol

Implementación del algoritmo Branch & Price para resolver coloraciones parciales sin conflictos (conflict-free partial coloring).

## Requisitos
- Compilador C++ con soporte C++20 (g++ recomendado)
- IBM ILOG CPLEX (paths usados en la línea de compilación)
- Librerías y objetos en el directorio `exactcolors/`
- Sistema operativo Linux

## Compilación

Primero, compilar las dependencias de `exactcolors/`:

```bash
cd exactcolors/
make
cd ..
```

Luego, ejecutar desde la raíz del repositorio:

```bash
make
```

**Nota:** Ajustar rutas de include y librería de CPLEX en el Makefile según la instalación local.

## Uso

Tras compilar se generará el ejecutable `dpcp`. Ejecutar con:

```bash
./dpcp [opciones] <archivos_entrada>
```

Revisar `main.cpp` para detalles sobre opciones y formato de entrada.

## Estructura del proyecto

```
cfcol/
├── main.cpp              # Punto de entrada del programa
├── src/                  # Código fuente principal
│   ├── graph.cpp
│   ├── col.cpp
│   └── lp.cpp
|   └── ...
├── include/              # Cabeceras del proyecto
│   └── bp.hpp
│   └── ...
├── exactcolors/          # Dependencias externas (objetos compilados)
├── hglib/                # Librería auxiliar
├── input/                # Archivos de entrada
├── instances/            # Instancias de prueba
├── Makefile              # Script de compilación
└── README.md             # Este archivo
```

## Contribuciones

Reportar issues o enviar pull requests. Incluir descripción del cambio y pruebas mínimas.

## Licencia

Ver archivo LICENSE para detalles.
