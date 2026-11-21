import math

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

dx = 1/100
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


# Generem el vector T de l'instant inicial:
C = [0.5,1]

for coef in C:
    dt = coef*(0.01**2)
    n = 99
    T = []
    T_0 = 36.5

    for i in range(n):
      T.append(T_0/T_c)

    S = []
    for e in range(n):
      if e == 0 or e == n-1:
        S.append(dt*(1 + (T_0/T_c)/(dx**2)))
      else:
        S.append(dt)


    # Programem una funció que trobi la inversa d'una matriu donada pel mètode de Gauss:

    def Minversa(A):
      n = len(A)

      M = []
      for fila in A:
        M.append(fila.copy())

      # Per a cada columna i busquem una fila j = i amb l'element i diferent de zero.

      # Generem la matriu identitat que despres de les combinacions lineals de files i columnes sera la matriu inversa

      I = []
      for a in range(n):
        f = []
        for b in range(n):
          if a == b:
            f.append(1)
          else:
            f.append(0)
        I.append(f)

      i = 0 # Comencem per la columna i = 0
      for i in range(n): # Columna a columna

        pivot = 0        # Definimos el tipo de variable
        pivot_valor = 0
        for p in range(i,n):  # Averiguamos el valor de la componente i mas grande y su pivot
          if abs(M[p][i])>pivot_valor:
            pivot_valor = abs(M[p][i])
            pivot = p

        if pivot != i:  # si el pivot no esta en la diagonal lo movemos
          M[i], M[pivot] = M[pivot], M[i]  # Posem el pivot a la diagonal
          I[i], I[pivot] = I[pivot], I[i]  # Apliquem el canvi a la matriu I

        l = pivot_valor

        for e in range(n):
                M[i][e] /= l
                I[i][e] /= l

        for c in range(n):
                if c != i:

                    k = M[c][i]

                    for e in range(n):
                        M[c][e] -= M[i][e]*k
                        I[c][e] -= I[i][e]*k
      # for row in I:
        #print(row)
      return I

    n = 99
    A = []


    for i in range(n):
      fila = []
      for j in range(n):
        if i == j + 1 and i >= 1:
          fila.append(-(dt/(dx**2)))
        elif i == j:
          fila.append(1+(dt/(dx**2)*2))
        elif i == j-1:
          fila.append(-(dt/(dx**2)))
        else:
          fila.append(0)
      A.append(fila)

    A_inversa = Minversa(A)

    t = 0

    while t < t_a:
      B = []

      for e in range(n):
        B.append(T[e]+S[e])


      T_next = []

      for f in range(n):
        sum = 0
        for c in range(n):
          sum += A_inversa[f][c]*B[c]
        T_next.append(sum)

      T = T_next

      t += dt

    X = []
    for e in range(n+2):
      X.append(2 * e/n)

    Tf = [T_0]

    for e in range(n+1):
      if  e == n:
        Tf.append(T_0)
      else:
        Tf.append(T[e] * T_c)


    print(max(Tf))
    print(t*t_c)
    print(Tf[0])

    fig, ax = plt.subplots()
    ax.plot(X, Tf)
    a = 0.75
    b = 1.25
    ax.axvspan(a, b, color='blue', alpha=0.2)
    ax.axhline(y=36.5, color='grey', linestyle=':')
    ax.axhline(y=50, color='green', linestyle=':')
    if max(Tf) > 80:
        ax.axhline(y=80, color='red', linestyle=':')
    plt.title(f"Euler implícit, C = {coef}")
    plt.xlabel("Posició (cm)")       # Eje X
    plt.ylabel("Temperatura (ºC)")
    fig.subplots_adjust(bottom=0.20)
    fig.text(0.5, 0.02, f"Temps final: {t * t_c:.4f} s, Temp. màx.: {max(Tf):.4f} ºC",
             ha='center', fontsize=10)


plt.show()


