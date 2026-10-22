"""
Análisis de la precisión de distintos métodos numéricos para integrar
las ecuaciones de movimiento del sistema solar.
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from metodos_numericos.euler import paso_euler as euler
from metodos_numericos.euler_cromer import paso_euler_cromer as euler_cromer
from metodos_numericos.verlet import paso_verlet as verlet
from metodos_numericos.rk4 import paso_rk4 as rk4
from simulacion_3D import Planeta


G = 6.67428e-11
SEGUNDOS_DIA = 24 * 60 * 60
DELTA_T = 3600

PERIODOS_REFERENCIA = {
    "Mercurio": 87.969 * SEGUNDOS_DIA,
    "Venus": 224.701 * SEGUNDOS_DIA,
    "Tierra": 365.256 * SEGUNDOS_DIA,
    "Marte": 686.980 * SEGUNDOS_DIA,
}


def cargar_planetas():
    """
    Carga el CSV y devuelve una lista nueva de planetas.
    """

    planetas = []

    carpeta_script = Path(__file__).resolve().parent
    ruta_csv = carpeta_script / "data" / "ci_sistema_solar.csv"

    if not ruta_csv.exists():
        ruta_csv = carpeta_script.parent / "data" / "ci_sistema_solar.csv"

    if not ruta_csv.exists():
        raise FileNotFoundError(f"No se ha encontrado el archivo de datos: {ruta_csv}")

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

            planeta.nombre = fila["nombre"]
            planeta.sol = planeta.nombre == "Sol"
            planeta.periodo = 0

            planetas.append(planeta)

    if not any(p.sol for p in planetas):
        raise ValueError(
            "No se ha encontrado el Sol en el CSV. "
            "Comprueba el valor de la columna 'nombre'."
        )

    return planetas


def periodo_orbital(metodo, planetas):
    """
    Calcula los periodos orbitales y los guarda en p.periodo.

    El periodo se mide mediante el cambio angular acumulado respecto
    al Sol. Si no se completa una órbita dentro del tiempo máximo,
    p.periodo permanece en cero.
    """

    sol = next(p for p in planetas if p.sol)

    angulos_previos = {}
    angulos_acumulados = {}

    for p in planetas:
        p.periodo = 0

        if p.sol:
            continue

        angulos_previos[p] = np.arctan2(
            p.y - sol.y,
            p.x - sol.x,
        )
        angulos_acumulados[p] = 0.0

    # Evita que una simulación inestable se ejecute indefinidamente.
    tiempo_maximo = 2 * max(PERIODOS_REFERENCIA.values())
    tiempo = 0

    while tiempo < tiempo_maximo:
        if all(p.periodo != 0 for p in planetas if not p.sol):
            break

        metodo(planetas, DELTA_T)
        tiempo += DELTA_T

        for p in planetas:
            if p.sol or p.periodo != 0:
                continue

            phi = np.arctan2(
                p.y - sol.y,
                p.x - sol.x,
            )

            delta = phi - angulos_previos[p]

            # Corrige el salto angular entre -pi y pi.
            delta = np.arctan2(
                np.sin(delta),
                np.cos(delta),
            )

            angulos_acumulados[p] += delta
            angulos_previos[p] = phi

            if abs(angulos_acumulados[p]) >= 2 * np.pi:
                p.periodo = tiempo

    # Aviso de los planetas que no completaron una órbita.
    for p in planetas:
        if not p.sol and p.periodo == 0:
            print(
                f"Advertencia: no se ha podido calcular el periodo "
                f"de {p.nombre} en {tiempo / SEGUNDOS_DIA:.1f} días."
            )


def energia(planetas):
    """
    Calcula la energía mecánica total del sistema, en julios.
    """

    energia_cinetica = 0.0
    energia_potencial = 0.0

    for p in planetas:
        energia_cinetica += 0.5 * p.masa * (p.vx**2 + p.vy**2 + p.vz**2)

    # Cada pareja de cuerpos se cuenta una sola vez.
    for i in range(len(planetas)):
        p1 = planetas[i]

        for j in range(i + 1, len(planetas)):
            p2 = planetas[j]

            dx = p1.x - p2.x
            dy = p1.y - p2.y
            dz = p1.z - p2.z

            distancia = np.sqrt(dx**2 + dy**2 + dz**2)

            if distancia == 0:
                raise ValueError(f"Distancia nula entre {p1.nombre} y {p2.nombre}.")

            energia_potencial -= G * p1.masa * p2.masa / distancia

    return energia_cinetica + energia_potencial


def momento_angular(planetas):
    """
    Calcula el vector de momento angular total respecto al origen.
    """

    Lx = 0.0
    Ly = 0.0
    Lz = 0.0

    for p in planetas:
        Lx += p.masa * (p.y * p.vz - p.z * p.vy)
        Ly += p.masa * (p.z * p.vx - p.x * p.vz)
        Lz += p.masa * (p.x * p.vy - p.y * p.vx)

    return np.array([Lx, Ly, Lz])


def main():
    metodos = (
        ("Euler", euler),
        ("Euler-Cromer", euler_cromer),
        ("Velocity Verlet", verlet),
        ("RK4", rk4),
    )

    for nombre_metodo, metodo in metodos:
        # Cada método comienza con las mismas condiciones iniciales.
        planetas = cargar_planetas()

        periodo_orbital(metodo, planetas)

        print(f"\nMétodo: {nombre_metodo}")
        print("-" * 65)

        print(
            f"{'Planeta':<12}"
            f"{'Simulado (días)':>18}"
            f"{'Referencia (días)':>21}"
            f"{'Error (%)':>14}"
        )

        for p in planetas:
            if p.sol:
                continue

            referencia = PERIODOS_REFERENCIA.get(p.nombre)

            if referencia is None:
                print(f"{p.nombre:<12} Sin valor de referencia")
                continue

            if p.periodo == 0:
                print(f"{p.nombre:<12} Periodo no calculado")
                continue

            periodo_dias = p.periodo / SEGUNDOS_DIA
            referencia_dias = referencia / SEGUNDOS_DIA

            error = abs(p.periodo - referencia) / referencia * 100

            print(
                f"{p.nombre:<12}"
                f"{periodo_dias:>18.3f}"
                f"{referencia_dias:>21.3f}"
                f"{error:>13.5f}"
            )


main()
