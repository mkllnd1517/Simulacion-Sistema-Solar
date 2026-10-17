"""
Método de Runge-Kutta de cuarto orden para la simulación del sistema solar.

Este módulo se encarga exclusivamente de avanzar el estado físico
de los cuerpos. No depende de Ursina ni de la visualización 3D.
"""


def paso_rk4(planetas, dt):
    """
    Avanza todos los cuerpos un paso temporal mediante Runge-Kutta
    de cuarto orden.

    Parámetros
    ----------
    planetas : list[Planeta]
        Lista de cuerpos que interactúan gravitatoriamente.
    dt : float
        Paso temporal en segundos.
    """

    estado_inicial = [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in planetas]

    def aceleraciones(estado):
        """
        Calcula las aceleraciones correspondientes a un estado dado.
        """

        posiciones_originales = [(p.x, p.y, p.z) for p in planetas]

        try:
            # Situa los cuerpos en el estado que queremos calcular
            for planeta, valores in zip(planetas, estado):
                planeta.x, planeta.y, planeta.z = valores[:3]

            aceleraciones_estado = []

            # Calcular las aceleraciones de todos los cuerpos.
            for i, planeta in enumerate(planetas):
                fx_total = fy_total = fz_total = 0.0

                for j, otro in enumerate(planetas):
                    if i == j:
                        continue

                    fx, fy, fz = planeta.gravedad(otro)

                    fx_total += fx
                    fy_total += fy
                    fz_total += fz

                aceleraciones_estado.append(
                    (
                        fx_total / planeta.masa,
                        fy_total / planeta.masa,
                        fz_total / planeta.masa,
                    )
                )

            return aceleraciones_estado

        finally:
            # Restaurar las posiciones originales.
            for planeta, posicion in zip(planetas, posiciones_originales):
                planeta.x, planeta.y, planeta.z = posicion

    def derivadas(estado):
        """
        Devuelve las derivadas de las posiciones y velocidades.
        """

        aceleraciones_estado = aceleraciones(estado)

        return [
            (vx, vy, vz, ax, ay, az)
            for (_, _, _, vx, vy, vz), (ax, ay, az) in zip(estado, aceleraciones_estado)
        ]

    def estado_intermedio(estado, k, factor):
        """
        Calcula un estado intermedio a partir de una pendiente.
        """

        return [
            tuple(
                valor + factor * derivada
                for valor, derivada in zip(valores, derivadas_cuerpo)
            )
            for valores, derivadas_cuerpo in zip(estado, k)
        ]

    # Calcular las cuatro pendientes de RK4.
    k1 = derivadas(estado_inicial)

    k2 = derivadas(estado_intermedio(estado_inicial, k1, dt / 2))

    k3 = derivadas(estado_intermedio(estado_inicial, k2, dt / 2))

    k4 = derivadas(estado_intermedio(estado_inicial, k3, dt))

    # Combinar las pendientes para obtener el estado final.
    nuevos_estados = []

    for estado, d1, d2, d3, d4 in zip(estado_inicial, k1, k2, k3, k4):
        nuevo_estado = tuple(
            valor
            + (dt / 6) * (derivada_1 + 2 * derivada_2 + 2 * derivada_3 + derivada_4)
            for valor, derivada_1, derivada_2, derivada_3, derivada_4 in zip(
                estado, d1, d2, d3, d4
            )
        )

        nuevos_estados.append(nuevo_estado)

    # Aplicar simultáneamente los estados finales.
    for planeta, estado in zip(planetas, nuevos_estados):
        (
            planeta.x,
            planeta.y,
            planeta.z,
            planeta.vx,
            planeta.vy,
            planeta.vz,
        ) = estado
