"""
Método de Euler-Cromer para la simulación del sistema solar.

Este módulo se encarga exclusivamente de avanzar el estado físico
de los cuerpos. No depende de Ursina ni de la visualización 3D.
"""


def paso_euler_cromer(planetas, dt):
    """
    Avanza todos los cuerpos un paso temporal mediante Euler-Cromer.

    Parámetros
    ----------
    planetas : list[Planeta]
        Lista de cuerpos que interactúan gravitatoriamente.
    dt : float
        Paso temporal en segundos.
    """
    estado_inicial = [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in planetas]

    # 2. Calcular todas las fuerzas con las posiciones iniciales.
    fuerzas = []

    for i, planeta in enumerate(planetas):
        fx_total = fy_total = fz_total = 0.0

        for j, otro in enumerate(planetas):
            if i == j:
                continue

            fx, fy, fz = planeta.gravedad(otro)

            fx_total += fx
            fy_total += fy
            fz_total += fz

        fuerzas.append((fx_total, fy_total, fz_total))

    # 3. Calcular las nuevas posiciones y velocidades.
    # Euler-Cromer utiliza las velocidades y aceleraciones iniciales.
    nuevos_estados = []

    for planeta, estado, fuerza in zip(planetas, estado_inicial, fuerzas):
        x, y, z, vx, vy, vz = estado
        fx, fy, fz = fuerza

        ax = fx / planeta.masa
        ay = fy / planeta.masa
        az = fz / planeta.masa

        vx_nueva = vx + ax * dt
        vy_nueva = vy + ay * dt
        vz_nueva = vz + az * dt

        x_nueva = x + vx_nueva * dt
        y_nueva = y + vy_nueva * dt
        z_nueva = z + vz_nueva * dt

        nuevos_estados.append((x_nueva, y_nueva, z_nueva, vx_nueva, vy_nueva, vz_nueva))

    # 4. Aplicar los nuevos estados una vez calculados todos.
    for planeta, estado in zip(planetas, nuevos_estados):
        (planeta.x, planeta.y, planeta.z, planeta.vx, planeta.vy, planeta.vz) = estado
