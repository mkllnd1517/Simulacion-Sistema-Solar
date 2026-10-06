"""
Simulación simple 2D del sistema solar utilizando Pygame.
Es la primera versión del proyecto, basado en la interacción
gravitatoria de Newton.
"""

import math

import pygame

pygame.init()

ANCHO = 800
ALTO = 800
VENTANA = pygame.display.set_mode((ANCHO, ALTO))
pygame.display.set_caption("Sistema Solar 2D")

BLANCO = (255, 255, 255)
AMARILLO = (255, 255, 0)
AZUL = (100, 149, 237)
ROJO = (188, 39, 50)
GRIS_OSCURO = (80, 78, 81)


FUENTE = pygame.font.SysFont("comicsans", 16)


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
        Dibuja el planeta en la ventana de Pygame.
        """
        x = self.x * self.ESCALA + ANCHO / 2
        y = self.y * self.ESCALA + ALTO / 2
        # centrando y escalando el planteta en la ventana de pygame

        if len(self.orbita) > 2:
            puntos_actualizados = []
            for punto in self.orbita:
                x, y = punto
                x = x * self.ESCALA + ANCHO / 2
                y = y * self.ESCALA + ALTO / 2
                puntos_actualizados.append((x, y))

            pygame.draw.lines(ventana, self.color, False, puntos_actualizados, 2)

        if not self.sol:
            distancia_texto = FUENTE.render(
                f"{self.distancia_al_sol / 1000:.1f} km", 1, BLANCO
            )
            VENTANA.blit(
                distancia_texto,
                (
                    x - distancia_texto.get_width() / 2,
                    y - distancia_texto.get_height() / 2,
                ),
            )

        pygame.draw.circle(ventana, self.color, (x, y), self.radio)

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


def main():
    run = True
    reloj = pygame.time.Clock()  # para regular los FPS

    sol = Planeta(0, 0, 30, AMARILLO, 1.98892 * 10**30)
    sol.sol = True

    tierra = Planeta(-1 * Planeta.UA, 0, 16, AZUL, 5.9742 * 10**24)
    tierra.y_vel = 29.783 * 1000  # velocidad inicial de la tierra

    marte = Planeta(-1.524 * Planeta.UA, 0, 12, ROJO, 6.39 * 10**23)
    marte.y_vel = 24.077 * 1000  # velocidad inicial de marte

    mercurio = Planeta(0.387 * Planeta.UA, 0, 8, GRIS_OSCURO, 3.30 * 10**23)
    mercurio.y_vel = -47.4 * 1000  # velocidad inicial de mercurio

    venus = Planeta(0.723 * Planeta.UA, 0, 14, BLANCO, 4.8685 * 10**24)
    venus.y_vel = -35.02 * 1000  # velocidad inicial de venus

    planetas = [sol, tierra, marte, mercurio, venus]

    while run:
        reloj.tick(60)  # limita a 60 FPS máx
        VENTANA.fill((0, 0, 0))  # limpia la ventana
        for evento in pygame.event.get():
            if evento.type == pygame.QUIT:
                run = False

        for planeta in planetas:
            planeta.actualizar_posicion(planetas)
            planeta.dibujar(VENTANA)

        pygame.display.update()

    pygame.quit()


main()
