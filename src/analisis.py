"""
En este módulo se va a analizar la precisión de la simulación al usar distintos
métodos numéricos para la integración de las ecuaciones de movimiento.
De esta forma, se reutilizará el código del método más eficiente para las partes posteriores
del proyecto.
"""

import numpy as np
import matplotlib.pyplot as plt
from metodos_numericos.euler import paso_euler as euler
from metodos_numericos.euler_cromer import paso_euler_cromer as euler_cromer
from metodos_numericos.verlet import paso_verlet as verlet
from metodos_numericos.rk4 import paso_rk4 as rk4


DELTA_T = 3600

planetas = []

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


def periodo_orbital():
    """
    Compruebas cuánto tarda en completar un giro de 360º en el
    plano de las órbitas.
    """        
