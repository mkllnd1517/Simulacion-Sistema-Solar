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
ALTO = 800


BLANCO = (255, 255, 255)
AMARILLO = (255, 255, 0)
AZUL = (100, 149, 237)
ROJO = (188, 39, 50)
GRIS_OSCURO = (80, 78, 81)


class Planeta:
    UA = 149.6e6 * 1000  # Unidades astronómicas en metros
    G = 6.67428e-11  # Constante gravitacional
    ESCALA = 200 / UA
    DELTA_T = 3600 * 24  # Actualización = 1 día

    def __init__(self, x, y, radio, color, masa):
        self.x = x
        self.y = y
        self.radio = radio
        self.color = color
        self.masa = masa

        self.orbita = []
        self.sol = (
            False  # para distinguir el sol de los planetas para pintar órbita o no
        )
        self.distancia_al_sol = 0

        self.x_vel = 0
        self.y_vel = 0

    def dibujar(self, ventana):
        """
        Dibuja el planeta en la ventana
        """

    def gravedad(self, otro):
        """
        Calcula la fuerza gravitatoria entre este planeta y otro objeto.
        """
        otro_x, otro_y = otro.x, otro.y
        distancia_x = otro_x - self.x
        distancia_y = otro_y - self.y
        distancia = math.sqrt(distancia_x**2 + distancia_y**2)

        if otro.sol:
            self.distancia_al_sol = distancia
        fuerza = self.G * self.masa * otro.masa / distancia**2
        theta = math.atan2(
            distancia_y, distancia_x
        )  # sacamos el ángulo para descomponer F
        fuerza_x = math.cos(theta) * fuerza
        fuerza_y = math.sin(theta) * fuerza
        return fuerza_x, fuerza_y

    def actualizar_posicion(self, planetas):
        """
        Actualiza la posición del planeta en función de las fuerzas gravitatorias de otros planetas.
        """
        fx_total = fy_total = 0
        for planeta in planetas:
            if self == planeta:
                continue
            fx, fy = self.gravedad(planeta)
            fx_total += fx
            fy_total += fy

        self.x_vel += fx_total / self.masa * self.DELTA_T  # 2ª ley de Newton + v=at
        self.y_vel += fy_total / self.masa * self.DELTA_T

        self.x += self.x_vel * self.DELTA_T
        self.y += self.y_vel * self.DELTA_T
        self.orbita.append((self.x, self.y))


def main(data):
    with open("data/ci_sistema_solar.csv", "r") as archivo:
        data = archivo.read().splitlines()
