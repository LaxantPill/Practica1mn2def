
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

T_c = 1348.57/2
t_c = 2855.38
t_f = 0.025
N_x = 101
x = np.linspace(0, 2, N_x)
g = [0.25,0.49]


def sim(c, N_x, T_c, punt):
    dx = 1 / (N_x - 1)
    dt = c * dx**2
    N_t = int(t_f / dt) + 1

    Tn = np.zeros((N_x, N_t))
    Tn[:, 0] = 36.5 / T_c
    Tn[0, :] = 36.5 / T_c
    Tn[-1, :] = 36.5 / T_c

    pas = None

    for j in range(N_t - 1):
        for i in range(1, N_x - 1):
            Tn[i, j+1] = Tn[i, j] + (dt / dx**2) * (Tn[i+1, j] - 2*Tn[i, j] + Tn[i-1, j]) + dt

        if Tn[punt, j+1] * T_c >= 50:
            pas = j + 1
            break

    if pas is not None:
        Tn = Tn[:, :pas + 1]

    Tn *= T_c
    return Tn, pas, dt

for c in g:
 T, pas, dt = sim(c, N_x, T_c, punt=38)

 fig, ax = plt.subplots(figsize=(8, 5))
 ax.plot(x, T[:, 0], label=r'Condició inicial ($t=0$)', linestyle='--', color='gray')
 ax.axhline(y=50, color='blue', linestyle=':', linewidth=1.5, label=r'$T=50\,^\circ\mathrm{C}$')
 ax.axvspan(0.75, 1.25, color='blue', alpha=0.25, label=r'$0.75 \leq x \leq 1.25\,\mathrm{cm}$')

 linia, = ax.plot(x, T[:, 0], color='crimson', linewidth=2, label='Perfil de temperatura')

 ax.set_xlim(0, 2)
 ax.set_ylim(35, np.max(T) + 2)
 ax.set_xlabel('x (cm)', fontsize=14)
 ax.set_ylabel('Temperatura (°C)', fontsize=14)
 ax.legend(loc='upper right')
 ax.grid(alpha=0.2)


 n_passos = T.shape[1]
 salt = max(1, n_passos // 100)
 frames = np.unique(np.append(np.arange(0, n_passos, salt), n_passos - 1))

 def actualitzar(frame):
    linia.set_ydata(T[:, frame])
    temps_s = frame * dt * t_c
    ax.set_title(f'Evolució tèrmica | pas = {frame} | t = {temps_s:.2f} s')
    return linia,

ani = FuncAnimation(fig, actualitzar, frames=frames, interval=40, repeat=True)

ani.save(f"animacio_c_{c}.gif", writer="pillow", fps=25)

plt.show()