import matplotlib.pyplot as plt

# Constants del problema en unitats SI

dx = 1/100 # m
c_v = 3698 # J/(kg·K)
p = 1081 # kg/m^3
k = 0.56 # W/(m·K)
s = 0.472 # S / m
T_0 = 36.5 # C
t_a = 0.025 # s
t = 0 # s
l = 0.02 # m

# Constants de normalitzacio:

T_c = s*0.5*(40**2)/k
t_c = p*c_v*(l**2)/k
x_c = l


# Funció per calcular T_i^j en funció de temperatures en el temps i-1

def sim(T_0, T_c, t_c, C):

  t = 0
  t_a = 0.025
  T_i = T_0/T_c
  T_0 = T_0/T_c
  T = []

  for i in range(0,101):
        T.append(T_i)

  dt = C*(0.01**2)
  steps = 0

  while t < t_a:
    steps += 1
    T_next = []
    i = 0
    while i<101:
      if i == 0 or i==100:
        T_next.append(T_0)
      else:
        T_next.append(T[i]+((T[i-1]-2*T[i]+T[i+1])/(0.01**2)+1)*dt)
      i += 1
    T = T_next
    t += dt

  for e in range(101):
    T[e] *= T_c

  print(steps)

  return T, t

# Generem l'eix horitzontal

X = []
for i in range(0,101):
  X.append(2*(i/100))

# Fem plot per cada valor de C

C = [0.51, 0.49, 0.25]
for e in C:
    Y = sim(T_0, T_c, t_c, e)
    fig, ax = plt.subplots()
    ax.plot(X, Y[0])
    a = 0.75
    b = 1.25
    ax.axvspan(a, b, color='blue', alpha=0.2)
    ax.axhline(y=36.5, color='grey', linestyle=':')
    ax.axhline(y=50, color='green', linestyle=':')
    if max(Y[0]) > 80:
        ax.axhline(y=80, color='red', linestyle=':')
    plt.title(f"Euler explícit, C = {e}")
    plt.xlabel("Posició (cm)")       # Eje X
    plt.ylabel("Temperatura (ºC)")
    fig.subplots_adjust(bottom=0.20)
    fig.text(0.5, 0.02, f"Temps final: {Y[1] * t_c:.4f} s, Temp. màx.: {max(Y[0]):.4f} ºC",
             ha='center', fontsize=10)
    if e != C[-1]:
        plt.show(block=False)
    else:
        plt.show()