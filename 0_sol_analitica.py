import numpy as np
import matplotlib.pyplot as plt

# Constants del problema en unitats SI

dx = 1/100 # m
c_v = 3698 # J/(kg·K)
p = 1081 # kg/m^3
k = 0.56 # W/(m·K)
s = 0.472 # S/m
T_0 = 36.5 # C
t_a = 0.025 #temps adimensional
t = 0 # s
l = 0.02 # m

# Constants de normalitzacio:

T_c = s*0.5*(40**2)/k
t_c = p*c_v*(l**2)/k
x_c = l


def T_analitica(x,t,M):
  suma = 0
  for m in range(1, M+1):
      n = 2*m - 1
      coef = 4*np.exp(-(n**2)*(np.pi**2)*t)/((n**3)*(np.pi**3))
      suma += coef*np.sin(n*np.pi*x)
  Tadim = T_0/T_c + x*(1-x)/2 - suma
  return Tadim*T_c

xn = []          # posició normalitzada (de 0 a 1)
X = []           # posició en cm (0 a 2, teixit d'estudi)
Y = []           # T_analitica per cada xn
N = 101

for i in range(N):
    xn.append(i/100)
    X.append(2*i/100)

for x in xn:
    Y.append(T_analitica(x, t_a, M=3))   #amb M=3 l'error és de l'ordre de 10^{-6}, menor a l'error dels mètodes numèrics (10^{-3})

T_max = T_analitica(0.5, t_a, M=3)          #sabem que per a tot t el màxim de temperatura està al centre, x=1/2
# i que la temperatura en x=1/2 creix amb el temps fins a un límit $T_{\max}=T_0+\frac{T_c}{8}$ quan t tendeix a infinit
# és a dir que l'instant on s'assoleix el màxim de temperatura en l'interval [0,t_a] ha de ser t_a

print(f"T màxima: {T_max:.4f}ºC a t = {t_a}")

fig, ax = plt.subplots()
ax.plot(X, Y)
a = 0.75
b = 1.25
ax.axvspan(a, b, color='blue', alpha=0.2)
ax.axhline(y=36.5, color='grey', linestyle=':')
ax.axhline(y=50, color='green', linestyle=':')
plt.title(f"Solució analítica")
plt.xlabel("Posició (cm)")
plt.ylabel("Temperatura (ºC)")
plt.show()

# Ara busquem l'instant en què x=0.75 cm i x=1.25 cm superen els 50 ºC per poder comparar-ho amb Crank-Nicolson
# x=0.75 cm i x=1.25 cm són xn=0.375 i xn=0.625 (són els punts on comença la regió de teixit sa)

dt = 10**(-6)
t = 0

while T_analitica(0.375, t, M=20) < 50 or T_analitica(0.625, t, M=20) < 50:         #aquí fem servir una M més gran, M=20, perquè la solució convergeix més ràpid per a M més grans, tot i que amb M=3 l'error de truncament ja és prou petit en comparació a l'error dels mètodes numèrics que fem servir
    t += dt

t_max = t*t_c
error_t = dt*t_c          # error màxim pel pas dt (el valor real de t_max es troba entre t_max-error_t i t_max)

print(f"El temps màxim que pot durar el procediment sense perjudicar el teixit sa és de t={t_max:.3f} amb un error màxim de {error_t:.3f}s, ja que si superem aquest temps la regió sana sobrepassa els 50ºC.")