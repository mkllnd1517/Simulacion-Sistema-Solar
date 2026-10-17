"""
Método de Leapfrog/Velocity Verlet para la simulación del sistema solar.

Este módulo se encarga exclusivamente de avanzar el estado físico
de los cuerpos. No depende de Ursina ni de la visualización 3D.
"""


def paso_verlet(planetas, dt):
    """
    Avanza todos los cuerpos un paso temporal mediante Velocity Verlet.

    Parámetros
    ----------
    planetas : list[Planeta]
        Lista de cuerpos que interactúan gravitatoriamente.
    dt : float
        Paso temporal en segundos.
    """

    estado_inicial = [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in planetas]

    aceleraciones_iniciales = []

    for i, planeta in enumerate(planetas):
        fx_total = fy_total = fz_total = 0.0

        for j, otro in enumerate(planetas):
            if i == j:
                continue

            fx, fy, fz = planeta.gravedad(otro)

            fx_total += fx
            fy_total += fy
            fz_total += fz

        aceleraciones_iniciales.append(
            (
                fx_total / planeta.masa,
                fy_total / planeta.masa,
                fz_total / planeta.masa,
            )
        )

    # Calcular las nuevas posiciones.
    nuevas_posiciones = []

    for estado, aceleracion in zip(estado_inicial, aceleraciones_iniciales):
        x, y, z, vx, vy, vz = estado
        ax, ay, az = aceleracion

        x_nueva = x + vx * dt + 0.5 * ax * dt**2
        y_nueva = y + vy * dt + 0.5 * ay * dt**2
        z_nueva = z + vz * dt + 0.5 * az * dt**2

        nuevas_posiciones.append((x_nueva, y_nueva, z_nueva))

    # Aplicar simultáneamente las nuevas posiciones.
    for planeta, posicion in zip(planetas, nuevas_posiciones):
        planeta.x, planeta.y, planeta.z = posicion

    # Calcular las nuevas aceleraciones con las nuevas posiciones.
    aceleraciones_nuevas = []

    for i, planeta in enumerate(planetas):
        fx_total = fy_total = fz_total = 0.0

        for j, otro in enumerate(planetas):
            if i == j:
                continue

            fx, fy, fz = planeta.gravedad(otro)

            fx_total += fx
            fy_total += fy
            fz_total += fz

        aceleraciones_nuevas.append(
            (
                fx_total / planeta.masa,
                fy_total / planeta.masa,
                fz_total / planeta.masa,
            )
        )

    # Actualizar las velocidades usando ambas aceleraciones.
    for planeta, estado, a_inicial, a_nueva in zip(
        planetas,
        estado_inicial,
        aceleraciones_iniciales,
        aceleraciones_nuevas,
    ):
        _, _, _, vx, vy, vz = estado

        ax_inicial, ay_inicial, az_inicial = a_inicial
        ax_nueva, ay_nueva, az_nueva = a_nueva

        planeta.vx = vx + 0.5 * (ax_inicial + ax_nueva) * dt
        planeta.vy = vy + 0.5 * (ay_inicial + ay_nueva) * dt
        planeta.vz = vz + 0.5 * (az_inicial + az_nueva) * dt
