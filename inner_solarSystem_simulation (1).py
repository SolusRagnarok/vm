import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#hello
# Orbital parameters for Inner Solar System
# [January 1, 2000, 11:58:55.816 UTC -> January 1, 2000, 12:00 UTC](https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf)
# Data format: Name: [a (AU), e(Eccentricity), i (deg relative to Earth equatorial plane), Ω (Longitude of Ascending Node), ω (Argument of Perihelion/periapsis ), Period (Earth Years), Color]
# e= c/a = ((a**2 - b**2)**(1/2)) / a
planets = {
    'Mercury': [0.387, 0.2056, 7.00, 48.33, 29.12, 0.2408, 'gray'],
    'Venus':   [0.723, 0.0067, 3.39, 76.68, 54.88, 0.6152, 'orange'],
    'Earth':   [1.000, 0.0167, 0.00, -11.26, 114.20, 1.0000, 'dodgerblue'],
    'Mars':    [1.524, 0.0934, 1.85, 49.55, 286.50, 1.8808, 'red'],
}

planetsA = {
    'Mercury A': [0.375, 0,0,0,0, 0.241, 'gray'],
    'Venus A':  [0.725,0,0,0,0, 0.615, 'orange'],
    'Earth A':  [1.000,0,0,0,0, 1.0000, 'dodgerblue'],
    'Mars A': [1.538,0,0,0,0, 1.881, 'red'],
}

def solve_kepler(M, e, tolerance=1e-6):
    E = M
    for _ in range(10):
        delta_E = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= delta_E
        if np.all(np.abs(delta_E) < tolerance): break
    return E

def get_orbital_position(a, e, i_deg, Omega_deg, omega_deg, E):
    i, Omega, omega = np.radians([i_deg, Omega_deg, omega_deg])
    x_p = a * (np.cos(E) - e)
    y_p = a * np.sqrt(1 - e**2) * np.sin(E)
    # Rotate to 3D ecliptic coordinates [source: https://en.wikipedia.org/wiki/Orbital_elements#Position_and_velocity_from_orbital_elements]
    # 1. Rotate by argument of periapsis (omega)
    # 2. Rotate by inclination (i)
    # 3. Rotate by longitude of ascending node (Omega)
    x = (np.cos(Omega)*np.cos(omega) - np.sin(Omega)*np.sin(omega)*np.cos(i))*x_p + \
        (-np.cos(Omega)*np.sin(omega) - np.sin(Omega)*np.cos(omega)*np.cos(i))*y_p
    y = (np.sin(Omega)*np.cos(omega) + np.cos(Omega)*np.sin(omega)*np.cos(i))*x_p + \
        (-np.sin(Omega)*np.sin(omega) + np.cos(Omega)*np.cos(omega)*np.cos(i))*y_p
    z = (np.sin(omega)*np.sin(i))*x_p + (np.cos(omega)*np.sin(i))*y_p
    return x, y, z

# --- Setup Figure ---
fig = plt.figure(figsize=(16, 10), facecolor='black')

ax1 = fig.add_subplot(121, projection='3d', facecolor='black')
ax1.set_title("Inner Solar System", color='white')
limit = 1.8
ax1.set_xlim([-limit, limit])
ax1.set_ylim([-limit, limit])
ax1.set_zlim([-limit, limit])
ax1.xaxis.pane.fill = False
ax1.yaxis.pane.fill = False
ax1.zaxis.pane.fill = False
ax1.xaxis.pane.set_edgecolor('black')
ax1.yaxis.pane.set_edgecolor('black')
ax1.zaxis.pane.set_edgecolor('black')
ax1.tick_params(colors='white')

ax2 = fig.add_subplot(122, projection='3d', facecolor='black')
ax2.set_title("Aryabhata's Solar System", color='white')
ax2.set_xlim([-limit, limit])
ax2.set_ylim([-limit, limit])
ax2.set_zlim([-limit, limit])
ax2.xaxis.pane.fill = False
ax2.yaxis.pane.fill = False
ax2.zaxis.pane.fill = False
ax2.xaxis.pane.set_edgecolor('black')
ax2.yaxis.pane.set_edgecolor('black')
ax2.zaxis.pane.set_edgecolor('black')
ax2.tick_params(colors='white')

points = {}
pointsA = {}

for ax, data_dict, points_dict in zip([ax1, ax2], [planets, planetsA], [points, pointsA]):
    ax.scatter([0], [0], [0], color='yellow', s=150) # The Sun
    
    E_orbit = np.linspace(0, 2*np.pi, 2000)
    for name, d in data_dict.items():
        xo, yo, zo = get_orbital_position(d[0], d[1], d[2], d[3], d[4], E_orbit)
        ax.plot(xo, yo, zo, color=d[6], alpha=0.3)
        
        p, = ax.plot([], [], [], marker='o', color=d[6], markersize=8, label=name)
        points_dict[name] = p

def update(frame):
    t = frame * 0.00274 # Speed of simulation (Earth years per frame)
    
    for name, d in planets.items():
        M = (2 * np.pi * t) / d[5]
        E = solve_kepler(M, d[1])
        x, y, z = get_orbital_position(d[0], d[1], d[2], d[3], d[4], E)
        points[name].set_data([x], [y])
        points[name].set_3d_properties([z])

    for name, d in planetsA.items():
        M = (2 * np.pi * t) / d[5]
        E = solve_kepler(M, d[1])
        x, y, z = get_orbital_position(d[0], d[1], d[2], d[3], d[4], E)
        pointsA[name].set_data([x], [y])
        pointsA[name].set_3d_properties([z])
        
    return list(points.values()) + list(pointsA.values())

ani = animation.FuncAnimation(fig, update, frames=1000, interval=30, blit=False)
plt.tight_layout()
plt.show()