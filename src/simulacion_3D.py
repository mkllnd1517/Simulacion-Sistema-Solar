"""
Simulación 3D del sistema solar interno. Se busca mejorar la versión 2D
añadiendo una simulación en 3D, con opción de importar los datos
directamente y poder visualizar el sistema solar desde cualquier ángulo.
También se estudiará la deformación del espacio-tiempo y el efecto de
las lentes gravitacionales en haces de luz.
"""

import math

from ursina import *

import csv

ANCHO = 800
ALTO = 800  # para el tamaño de la ventana de la simulación


class Planeta:
    UA = 149.6e6 * 1000  # Unidades astronómicas en metros
    G = 6.67428e-11  # Constante gravitacional
    ESCALA = 200 / UA
    DELTA_T = 3600 * 24  # Actualización = 1 día

    def __init__(self, posicion, velocidad, radio, color, masa):
        self.x, self.y, self.z = posicion
        self.vx, self.vy, self.vz = velocidad

        self.radio = radio
        self.color = color
        self.masa = masa

        self.orbita = []
        self.sol = (
            False  # para distinguir el sol de los planetas para pintar órbita o no
        )
        self.distancia_al_sol = 0

    def dibujar(self, ventana):
        """
        Dibuja el planeta en la ventana
        """
        Entity(
            model="sphere",
            color=self.color,
            scale=4 * (1 + math.log1p(self.radio / 6.371e6)) / (1 + math.log(2)),
            position=(self.x * self.ESCALA, self.y * self.ESCALA, self.z * self.ESCALA),
        )  # el radio se toma con una compresión logarítmica para que se vea mejor la simulación

    def gravedad(self, otro):
        """
        Calcula la fuerza gravitatoria entre este planeta y otro objeto.
        """
        otro_x, otro_y, otro_z = otro.x, otro.y, otro.z
        distancia_x = otro_x - self.x
        distancia_y = otro_y - self.y
        distancia_z = otro_z - self.z
        distancia_cuadrado = distancia_x**2 + distancia_y**2 + distancia_z**2

        factor = (
            self.G
            * self.masa
            * otro.masa
            / distancia_cuadrado
            * math.sqrt(distancia_cuadrado)
        )

        if otro.sol:
            self.distancia_al_sol = math.sqrt(distancia_cuadrado)

        fuerza = self.G * self.masa * otro.masa / distancia_cuadrado

        f_x = factor * distancia_x
        f_y = factor * distancia_y
        f_z = factor * distancia_z
        return f_x, f_y, f_z

    def actualizar_posicion(self, planetas):
        """
        Actualiza la posición del planeta en función de las fuerzas gravitatorias de otros planetas.
        """
        fx_total = fy_total = f_z_total = 0
        for planeta in planetas:
            if self == planeta:
                continue
            fx, fy, fz = self.gravedad(planeta)
            fx_total += fx
            fy_total += fy
            f_z_total += fz

        self.vx += fx_total / self.masa * self.DELTA_T  # 2ª ley de Newton + v=at
        self.vy += fy_total / self.masa * self.DELTA_T
        self.vz += f_z_total / self.masa * self.DELTA_T

        self.x += self.vx * self.DELTA_T
        self.y += self.vy * self.DELTA_T
        self.z += self.vz * self.DELTA_T
        self.orbita.append((self.x, self.y, self.z))


def main():
    app = Ursina(title="Simulación 3D del Sistema Solar", size=(ANCHO, ALTO), fullscreen=False, vsync=True)
    EditorCamera()
    luz_solar = PointLight(
    parent=scene,
    color=color.white,
    position=(0, 0, 0),
)

    planetas = []

    with open("data/ci_sistema_solar.csv", "r") as archivo:
        lector = csv.DictReader(archivo)

        for fila in lector:
            p = Planeta(
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
                color=tuple(
                    int(c) for c in fila["color"].split(";")
                ),
                masa=float(fila["masa_kg"]),
            )

            p.sol = fila["nombre"] == "Sol"
            planetas.append(p)
            p.dibujar()

    acumulador = 0.0  # Acumulador para integrar con pasos de tiempo fijos.


    def actualizar():

        # me queda pendiente poner la actualización de la posición de los planetas
        # actualizar la posición de la bombilla con la posición del sol en el bucle
        # añadir algo que pinte la órbita de los planetas