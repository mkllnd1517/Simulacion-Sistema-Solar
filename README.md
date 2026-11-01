
# Simulación de Sistemas Gravitatorios

Proyecto en Python para estudiar distintos problemas de dinámica gravitatoria mediante métodos de integración numérica y visualización tridimensional interactiva con Ursina.

**Estado: en desarrollo.** El proyecto comienza con la simulación del Sistema Solar y está diseñado para ampliarse mediante módulos independientes que estudien otros fenómenos y configuraciones gravitatorias, reutilizando y adaptando los métodos numéricos desarrollados en la primera parte.

## 1. Simulación del Sistema Solar

La primera parte del proyecto implementa un modelo gravitatorio newtoniano para estudiar la evolución de un sistema de cuerpos sometidos a su interacción gravitatoria mutua.

La simulación permite visualizar el movimiento de los cuerpos en tres dimensiones, mientras que los módulos de análisis permiten estudiar la precisión de las soluciones numéricas y comparar distintos métodos de integración.

### Métodos numéricos

Se implementan cuatro métodos para resolver las ecuaciones del movimiento:

- **Euler:** método explícito de primer orden, utilizado como referencia inicial.
- **Euler-Cromer:** variante de Euler que actualiza la velocidad antes de calcular la nueva posición.
- **Velocity Verlet:** método de segundo orden con propiedades simplécticas, especialmente útil para estudiar sistemas mecánicos a largo plazo.
- **Runge-Kutta de cuarto orden (RK4):** método de cuarto orden que permite obtener una buena precisión por paso a cambio de un mayor número de evaluaciones de las ecuaciones.

La comparación de estos métodos permite estudiar cómo influyen el integrador y el paso temporal en los resultados obtenidos.

## 2. Estructura del proyecto

```text
.
├── data/
│   └── ci_sistema_solar.csv
├── src/
│   ├── simple_2D.py
│   ├── simulacion_3D.py
│   ├── analisis.py
│   └── metodos_numericos/
│       ├── euler.py
│       ├── euler_cromer.py
│       ├── verlet.py
│       └── rk4.py
│
├── resultados/
│   ├── analisis.md
│   ├── energia_mecanica.png
│   ├── error_energia.png
│   ├── error_momento_angular.png
│   └── momento_angular.png
│
├── pyproject.toml
├── uv.lock
└── README.md
```

### Archivos y directorios

- **`src/simple_2D.py`**: primera versión de la simulación del Sistema Solar, desarrollada en dos dimensiones.

- **`src/simulacion_3D.py`**: implementa la simulación tridimensional mediante Ursina. Gestiona la representación de los cuerpos, la cámara, las trayectorias orbitales y la evolución visual del sistema. Permite seleccionar el método numérico utilizado para integrar las ecuaciones del movimiento.

- **`src/metodos_numericos/euler.py`**: implementa el método de Euler.

- **`src/metodos_numericos/euler_cromer.py`**: implementa el método de Euler-Cromer.

- **`src/metodos_numericos/verlet.py`**: implementa el método Velocity Verlet.

- **`src/metodos_numericos/rk4.py`**: implementa el método Runge-Kutta de cuarto orden.

- **`src/analisis.py`**: permite analizar y comparar los resultados de los integradores, incluyendo la estimación de periodos orbitales y el estudio de las variaciones de la energía mecánica y el momento angular.

- **`data/ci_sistema_solar.csv`**: contiene las condiciones iniciales de los cuerpos, incluyendo sus masas, posiciones, velocidades y propiedades utilizadas para la representación visual.

- **`resultados/`**: almacena las gráficas generadas y documentos Markdown que describen los experimentos numéricos, interpretan los resultados y recogen sus conclusiones.

- **`pyproject.toml`**: define la configuración del proyecto y sus dependencias.

- **`uv.lock`**: registra las versiones resueltas de las dependencias para facilitar la reproducción del entorno.

La estructura separa la integración numérica, el análisis y la visualización, facilitando la incorporación de nuevos métodos y módulos sin concentrar toda la implementación en un único archivo.

## 3. Instalación y ejecución

El proyecto utiliza [uv](https://docs.astral.sh/uv/) para gestionar el entorno de Python y sus dependencias.

### Instalación

Clona el repositorio y accede a su directorio raíz:

```bash
git clone https://github.com/mkllnd1517/Simulaci-n-Sistema-Solar.git
cd Simulaci-n-Sistema-Solar
```

Instala las dependencias definidas en el proyecto:

```bash
uv sync
```

### Simulación 2D

Para ejecutar la primera versión de la simulación:

```bash
uv run python src/simple_2D.py
```

### Simulación 3D

La simulación tridimensional permite seleccionar el método de integración numérica desde el momento de su ejecución.

Desde el directorio raíz del proyecto, ejecuta:

```bash
uv run python -c "import sys; sys.path.insert(0, 'src'); from simulacion_3D import main; from metodos_numericos.verlet import paso_verlet; main(paso_verlet)"
```

En este ejemplo se utiliza Velocity Verlet.

Para utilizar otro integrador, basta con cambiar el método importado y el argumento de `main()`:

| Método | Función |
|---|---|
| Euler | `paso_euler` |
| Euler-Cromer | `paso_euler_cromer` |
| Velocity Verlet | `paso_verlet` |
| RK4 | `paso_rk4` |

Por ejemplo, para ejecutar la simulación con RK4:

```bash
uv run python -c "import sys; sys.path.insert(0, 'src'); from simulacion_3D import main; from metodos_numericos.rk4 import paso_rk4; main(paso_rk4)"
```

La simulación utiliza un paso temporal físico de una hora y una velocidad de visualización equivalente a un día simulado por segundo de tiempo real. La cámara permite explorar la escena desde diferentes perspectivas.

### Análisis numérico

Para ejecutar el módulo de análisis:

```bash
uv run python src/analisis.py
```

Este módulo permite realizar las comparaciones numéricas y estudiar la evolución de las magnitudes físicas sin necesidad de ejecutar la visualización tridimensional.



