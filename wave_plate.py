#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 22:20:58 2022

@author: dayron
"""

from dolfin import *
from ufl import nabla_div
import numpy as np

def simulacion(frec, alfa, amp, T0, T, ny, dt, num_steps):
    
    # Caracteristicas de la placa
    Lx = 5e-2
    Ly = 1e-3
    
    pho = 2700
    
    mu = 2.624e10
    lambda_ = 5.279e10
    
    # Definimos el dominio
    nx = int(ny*(Lx/Ly))
    
    mesh = RectangleMesh(Point(0,0), Point(Lx, Ly), nx, ny)
    V = VectorFunctionSpace(mesh, 'P', 1)
    
    # Definimos el pulso en omega4
    g_t = Expression(('0.0','amp*sin(2*pi*frec*t)*exp(-alfa*(pow(t-T0,2)/pow(T,2)))'), 
                     degree=1, amp=amp, frec=frec, alfa=alfa, T0=T0, T=T, t=0)
    
    # Definimos las condiciones de frontera
    tol = 1e-14
    
    def boundary(x, on_boundary):
        if on_boundary:
            if near(x[0], 0, tol):
                return True
            else:
                return False
        else:
            return False
        
    bc = DirichletBC(V, g_t, boundary)
    
    
    # Definimos las condiciones iniciales
    u_i1 = Function(V)
    u_i2 = Function(V)
    u_i3 = Function(V)
    
    c = Expression(('0.0','0.0'), degree=1)
    
    u_i1 = interpolate(c, V)
    u_i2 = interpolate(c, V)
    u_i3 = interpolate(c, V)
    
    fi = Function(V)
    fi = 5*u_i1 - 4*u_i2 + u_i3
    
    # Definicion del problema variacional
    u = TrialFunction(V)
    d = u.geometric_dimension() 
    v = TestFunction(V)
    
    # Define strain and stress
    def strain(u):
        return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    
    
    a = 2*pho*dot(u,v)*dx + 2*mu*dt**2*inner(strain(u), strain(v))*dx + lambda_*dt**2*nabla_div(u)*nabla_div(v)*dx
    L = pho*dot(fi, v)*dx
    
    vtkfile = File('plate/solution.pvd')
    
    # Solucion aproximada
    u = Function(V)
    
    # Desplazamientos en los puntos pi
    p1 = []
    p2 = []
    p3 = []
    p4 = []
    
    set_log_active(False)
    
    # Paso del tiempo
    t = 0
    times = []
    for n in range(num_steps):
        # Actualizamos el tiempo
        t += dt
        
        # Actualizamos la frontera
        g_t.t = t
        
        # Resolvemos el problema
        problem = LinearVariationalProblem(a, L, u, bc)
        solver = LinearVariationalSolver(problem)
        solver.parameters["linear_solver"] = "gmres"
        solver.parameters["preconditioner"] = "ilu"
        solver.solve()
        
        #solve(a == L, u, bc)
        
        # Actualizamos las soluciones previas
        u_i3.assign(u_i2)
        u_i2.assign(u_i1)
        u_i1.assign(u)
        
        # Save solution
        
        if (n % 10 == 0):
            vtkfile << (u, t)
        
            
        # Vamos a guardar los desplazamientos en los puntos pi
        p1.append(u(1e-2, Ly)[1])
        p2.append(u(1.3e-2, Ly)[1])
        p3.append(u(1.6e-2, Ly)[1])
        p4.append(u(1.9e-2, Ly)[1])
        times.append(t)
    
    '''
    # Ploteo de los desplazamientos
    fig = plt.figure(figsize=(14,8))
    
    ax_5 = fig.add_subplot(122) 
    ax_1 = fig.add_subplot(421)
    ax_2 = fig.add_subplot(423)
    ax_3 = fig.add_subplot(425)
    ax_4 = fig.add_subplot(427)

    ax_1.plot(times, p1)
    ax_2.plot(times, p2)
    ax_3.plot(times, p3)
    ax_4.plot(times, p4)
    '''
    
    # Vamos a hallar el maximo de los desplazamientos
    max1 = np.max(p1)
    index1 = np.argmax(p1)
    ind01 = np.argmin(p1[index1:]) + index1
    
    max2 = np.max(p2[ind01:])
    index2 = np.argmax(p2[ind01:]) + ind01
    ind02 = np.argmin(p2[index2:]) + index2
    
    max3 = np.max(p3[ind02:])
    index3 = np.argmax(p3[ind02:]) + ind02
    ind03 = np.argmin(p3[index3:]) + index3
    
    max4 = np.max(p4[ind03:])
    index4 = np.argmax(p4[ind03:]) + ind03
    
    '''
    ax_1.plot([times[index1]],[max1],'or')
    ax_2.plot([times[index2]],[max2],'or')
    ax_3.plot([times[index3]],[max3],'or')
    ax_4.plot([times[index4]],[max4],'or')
    '''
    
    # Ya podemos calcular la velocidad de fase C(f0)
    tiempos = [times[index1], times[index2], times[index3], times[index4]]
    desplazamientos = [1e-2, 1.3e-2, 1.6e-2, 1.9e-2]
    suma = 0
    
    '''
    ax_5.plot(tiempos, desplazamientos, '-ro')
    '''
    
    for i in range(4):
        suma += desplazamientos[i]/tiempos[i] 
        
    C_f0 = suma/4
    
    print("La velocidad de fase para %d Hz es %.4f m/s" % (frec, C_f0))
    
    return C_f0
    
    
    

    


    
    





