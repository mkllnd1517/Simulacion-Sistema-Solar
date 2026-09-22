"""
Simulación simple 2D del sistema solar utilizando Pygame.
Es la primera versión del proyecto, basado en la interacción
gravitatoria de Newton.
"""

import pygame
import math


pygame.init()

ANCHO = 800
ALTO = 800
VENTANA = pygame.display.set_mode((ANCHO, ALTO))
pygame.display.set_caption("Sistema Solar 2D")


def main():
    run = True

    while run:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False

    pygame.quit()
