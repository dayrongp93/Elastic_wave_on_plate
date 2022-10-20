#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 11:41:48 2022

En este archivo vamos a llamar a la funcion que ejecuta la simulacion

@author: dayron
"""

import wave_plate
import matplotlib.pyplot as plt

velocidades = []
frecuencias = []

# Primer caso
# ----------------------------------------------------------------------------
frec = 100000.0
alfa = 2.5e-2
amp = 1e-3
T0 = 2.3e-5
T = 2e-6

ny = 5

dt = 1e-7

num_steps = 500

vel_aux = wave_plate.simulacion(frec, alfa, amp, T0, T, ny, dt, num_steps)

frecuencias.append(frec)
velocidades.append(vel_aux)
# ----------------------------------------------------------------------------


# Segundo caso
# ----------------------------------------------------------------------------
frec = 200000.0
alfa = 1.5e-1
amp = 1e-3
T0 = 1.2e-5
T = 2e-6

ny = 5

dt = 1e-7

num_steps = 300

vel_aux = wave_plate.simulacion(frec, alfa, amp, T0, T, ny, dt, num_steps)

frecuencias.append(frec)
velocidades.append(vel_aux)
# ----------------------------------------------------------------------------

# Tercer caso
# ----------------------------------------------------------------------------
frec = 300000.0
alfa = 8e-1
amp = 1e-3
T0 = 8e-6
T = 2e-6

ny = 6

dt = 1e-8

num_steps = 2200

vel_aux = wave_plate.simulacion(frec, alfa, amp, T0, T, ny, dt, num_steps)

frecuencias.append(frec)
velocidades.append(vel_aux)
# ----------------------------------------------------------------------------


# Cuarto caso
# ----------------------------------------------------------------------------
frec = 600000.0
alfa = 1.1
amp = 1e-3
T0 = 4e-6
T = 2e-6

ny = 6

dt = 1e-8

num_steps = 1500

vel_aux = wave_plate.simulacion(frec, alfa, amp, T0, T, ny, dt, num_steps)

frecuencias.append(frec)
velocidades.append(vel_aux)
# ----------------------------------------------------------------------------


# Quinto caso
# ----------------------------------------------------------------------------
frec = 700000.0
alfa = 1.5
amp = 1e-3
T0 = 3e-6
T = 2e-6

ny = 6

dt = 1e-8

num_steps = 1500

vel_aux = wave_plate.simulacion(frec, alfa, amp, T0, T, ny, dt, num_steps)

frecuencias.append(frec)
velocidades.append(vel_aux)
# ----------------------------------------------------------------------------

# Sexto caso
# ----------------------------------------------------------------------------
frec = 900000.0
alfa = 2
amp = 1e-3
T0 = 2.6e-6
T = 2e-6

ny = 6

dt = 1e-8

num_steps = 1500

vel_aux = wave_plate.simulacion(frec, alfa, amp, T0, T, ny, dt, num_steps)

frecuencias.append(frec)
velocidades.append(vel_aux)
# ----------------------------------------------------------------------------

# Septimo caso
# ----------------------------------------------------------------------------
frec = 1100000.0
alfa = 3.2
amp = 1e-3
T0 = 2.3e-6
T = 2e-6

ny = 6

dt = 1e-8

num_steps = 1500

vel_aux = wave_plate.simulacion(frec, alfa, amp, T0, T, ny, dt, num_steps)

frecuencias.append(frec)
velocidades.append(vel_aux)
# ----------------------------------------------------------------------------

# Octavo caso
# ----------------------------------------------------------------------------
frec = 2000000.0
alfa = 10
amp = 1e-3
T0 = 1e-6
T = 2e-6

ny = 6

dt = 1e-8

num_steps = 1500

vel_aux = wave_plate.simulacion(frec, alfa, amp, T0, T, ny, dt, num_steps)

frecuencias.append(frec)
velocidades.append(vel_aux)
# ----------------------------------------------------------------------------

plt.plot(frecuencias, velocidades, '-o')

