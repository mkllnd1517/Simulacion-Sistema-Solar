"""
Simulación 3D del sistema solar interno. Se busca mejorar la versión 2D
añadiendo una simulación en 3D, con opción de importar los datos
directamente y poder visualizar el sistema solar desde cualquier ángulo.
Se observa que, importando otros datos, se puede simular cualquier sistema planetario
y gran variedad de sistemas gravitatorios, siempre que no haya colisiones entre cuerpos,
la aproximación de cuerpos puntuales sea válida y no haga falta tener en cuenta la relatividad.
"""

import math
import csv
from pathlib import Path

from ursina import *
from metodos_numericos.euler import paso_euler

ANCHO = 800
ALTO = 800  # para el tamaño de la ventana de la simulación
VELOCIDAD_SIMULACION = (
    24 * 3600
)  # segundos simulados por segundo real, un día cada segundo

# Variables que se usan en update() y en main().
planetas = []
lineas_orbita = {}
luz_solar = None
acumulador = 0.0
acumulador_orbitas = 0.0
pasos_desde_traza = 0


class Planeta:
    UA = 149.6e6 * 1000  # Unidades astronómicas en metros
    G = 6.67428e-11  # Constante gravitacional
    ESCALA = 200 / UA
    DELTA_T = 60 * 60  # un paso físico de una hora

    def __init__(self, posicion, velocidad, radio, rgb: tuple[int, int, int], masa):
        self.x, self.y, self.z = map(float, posicion)
        self.vx, self.vy, self.vz = map(float, velocidad)

        self.radio = float(radio)
        self.color = rgb
        self.masa = float(masa)

        self.orbita = [(self.x, self.y, self.z)]
        self.sol = False
        # para distinguir el sol de los planetas para pintar órbita o no

        self.distancia_al_sol = 0
        self.entidad = None

    def dibujar(self):
        self.entidad = Entity(
            model="sphere",
            color=color.rgb32(*self.color),
            scale=(4 * (1 + math.log1p(self.radio / 6.371e6)) / (1 + math.log(2))),
            position=(
                self.x * self.ESCALA,
                self.y * self.ESCALA,
                self.z * self.ESCALA,
            ),
        )  # el radio se toma con una compresión logarítmica para que se vea mejor la simulación

    def gravedad(self, otro):
        distancia_x = otro.x - self.x
        distancia_y = otro.y - self.y
        distancia_z = otro.z - self.z
        distancia_cuadrado = distancia_x**2 + distancia_y**2 + distancia_z**2

        # Evita dividir entre cero si dos cuerpos tienen la misma posición.
        if distancia_cuadrado == 0:
            return 0, 0, 0

        factor = (
            self.G
            * self.masa
            * otro.masa
            / (distancia_cuadrado * math.sqrt(distancia_cuadrado))
        )

        if otro.sol:
            self.distancia_al_sol = math.sqrt(distancia_cuadrado)

        f_x = factor * distancia_x
        f_y = factor * distancia_y
        f_z = factor * distancia_z
        return f_x, f_y, f_z

    def actualizar_visual(self):
        self.entidad.position = (
            self.x * self.ESCALA,
            self.y * self.ESCALA,
            self.z * self.ESCALA,
        )


def pinta_orbitas():
    """
    Actualiza las líneas que representan las órbitas de los cuerpos.
    """

    for planeta in planetas:
        if len(planeta.orbita) < 2:
            continue

        vertices = [
            (x * Planeta.ESCALA, y * Planeta.ESCALA, z * Planeta.ESCALA)
            for x, y, z in planeta.orbita
        ]

        linea_anterior = lineas_orbita.get(planeta)

        if linea_anterior is not None:
            destroy(linea_anterior)

        lineas_orbita[planeta] = Entity(
            model=Mesh(
                vertices=vertices,
                mode="line",
                thickness=2,
            ),
            color=color.rgb(*planeta.color),
        )


def update():
    """
    Actualiza la simulación en cada fotograma de Ursina.
    El cálculo físico se delega al método numérico seleccionado.
    """
    global acumulador, acumulador_orbitas, pasos_desde_traza

    # Convertir el tiempo real transcurrido en tiempo simulado.
    acumulador += time.dt * VELOCIDAD_SIMULACION

    # Ejecutar tantos pasos físicos como sean necesarios.
    while acumulador >= Planeta.DELTA_T:
        paso_euler(planetas, Planeta.DELTA_T)
        acumulador -= Planeta.DELTA_T

        pasos_desde_traza += 1

        # Guardar un punto de la órbita cada seis horas simuladas.
        if pasos_desde_traza >= 6:
            for planeta in planetas:
                planeta.orbita.append((planeta.x, planeta.y, planeta.z))

            pasos_desde_traza = 0

    # Actualizar la representación visual de los cuerpos.
    for planeta in planetas:
        planeta.actualizar_visual()

        if planeta.sol and luz_solar is not None:
            luz_solar.position = planeta.entidad.position

    # Actualizar las líneas orbitales cada medio segundo real.
    acumulador_orbitas += time.dt

    if acumulador_orbitas >= 0.5:
        pinta_orbitas()
        acumulador_orbitas = 0.0


def main():
    global planetas, lineas_orbita, luz_solar

    app = Ursina(
        title="Simulación 3D del Sistema Solar",
        size=(ANCHO, ALTO),
        fullscreen=False,
        vsync=True,
    )
    EditorCamera()
    luz_solar = PointLight(parent=scene, color=color.white, position=(0, 0, 0))

    # Busca el CSV junto al script (src/data) y después en la carpeta del proyecto (data).
    carpeta_script = Path(__file__).resolve().parent
    ruta_csv = carpeta_script / "data" / "ci_sistema_solar.csv"
    if not ruta_csv.exists():
        ruta_csv = carpeta_script.parent / "data" / "ci_sistema_solar.csv"

    if not ruta_csv.exists():
        raise FileNotFoundError(
            f"No se encontró el archivo CSV. Se buscó en: {ruta_csv}"
        )

    with open(ruta_csv, "r", newline="", encoding="utf-8") as archivo:
        lector = csv.DictReader(archivo)

        for fila in lector:
            planeta = Planeta(
                posicion=(
                    float(fila["x_UA"]) * Planeta.UA,
                    float(fila["y_UA"]) * Planeta.UA,
                    float(fila["z_UA"]) * Planeta.UA,
                ),
                velocidad=(
                    float(fila["vx_m_s"]),
                    float(fila["vy_m_s"]),
                    float(fila["vz_m_s"]),
                ),
                radio=float(fila["radio_m"]),
                rgb=tuple(int(c) for c in fila["color"].split(";")),
                masa=float(fila["masa_kg"]),
            )

            planeta.sol = fila["nombre"] == "Sol"
            planetas.append(planeta)
            planeta.dibujar()

    lineas_orbita = {planeta: None for planeta in planetas}
    app.run()


if __name__ == "__main__":
    main()
