# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 00:07:16 2021

@author: TELMA
"""

import numpy as np

"""Funciones auxiliares"""

def CabeceraInit(T):
    variables = [f'x{i+1}' for i in range(T.shape[1]-1)]
    return variables

def BaseInit(base):
    variables = [f'x{i}' for i in base]
    return variables

def PrintTabla(T, variables, base):
    T = T.astype(float)
    #variables = [f'x{i+1}' for i in range(T.shape[1]-1)]
    row_format ="{:>10}|" * (len(variables + ['Sol']) + 1)
    header = row_format.format("", *variables + ['Sol'])
    separador = len(header)*'-'
    func_objetivo = row_format.format("Z", *T[-1])
    print(header)
    print(separador)
    
    for base, row in zip(base, T[0:-1]):
        print(row_format.format(base, *row))
        print(separador)
    
    print(func_objetivo)
    print(separador)
    print()
    
def TabInit(z, A, b):
    T = np.concatenate((A, z), axis=0)
    T = np.concatenate((T, np.array([np.append(b[0], [0])]).T), axis=1)
    #print(T)
    return T

def InicTb(T, base):
    base = [x-1 for x in base]
    Tb = T[:, base]
    return Tb

"""Funciones para usar el método Dantzig"""

def Dantzig(T):
    T = T.astype(float)
    j = np.where(T[-1] == min(T[-1, 0:-1]))
    m = T[:, -1]
    n = T[:, j[0][0]]
    
    arr = m/n
    arr = arr.astype(float)
    arr = np.delete(arr, -1, 0)
    filter_arr = arr > 0
    arrfilt = arr[filter_arr]
    
    i = np.where(arr == min(arrfilt))
    print(i[0][0], j[0][0])
    return i[0][0], j[0][0]

def FinDantzig(T):
    arr = T[-1]
    arr = np.delete(arr, -1, 0)
    filter_arr = arr < 0
    print("prueba fin dantzig")
    print(filter_arr)
    if True in filter_arr:
        return True
    return False

def FinDantzig2(T):
    arr = T[-1]
    arr = np.delete(arr, -1, 0)
    print("valor arr")
    print(arr)
    filter_arr = arr < 0
    print("prueba fin dantzig2")
    print(filter_arr)
    if True in filter_arr:
        return False
    return True


"""Funcion auxiliar, pivotar en cada paso"""

def Pivot(T, i, j):
    T = T.astype(float)
    if T[i, j]!=1:
        T[i] = (1/T[i, j])*T[i]
        #print(T[i])
    
    indice = np.delete(np.arange(T.shape[0]), i)
    #print(indice)
    for k in indice:
        #print(T[k])
        T[k] = T[k] - (T[k, j]/T[i, j])*T[i]
        """print("en la iteracion pasa a")
        print(T[k])
        print("--------------")"""
    return np.around(T, decimals= 4)

"""Metodo Simplex"""

def AuxSimplex(T, base):
    while FinDantzig(T):
        i, j = Dantzig(T)
        T = Pivot(T, i, j)
        base[i] = i+1
        PrintTabla(T, CabeceraInit(T), BaseInit(base))
        
    print(base)
    #gauss
    for it in range(len(base)):
        i = np.where(np.array([T[:, base[it]-1]])[0] == 1)
        print(np.array([T[:, base[it]-1]])[0])
        print(i[0][0])
        T = Pivot(T, i[0][0], base[it]-1)
        print(base[it]-1, it)
        PrintTabla(T, CabeceraInit(T), BaseInit(base))
    
    return T, base
    
def Simplex(z, A, b, base):
    T = TabInit(z, A, b)
    T = T.astype(float)
    PrintTabla(T, CabeceraInit(T), BaseInit(base))
    
    T, base = AuxSimplex(T, base)
    
    print(FinDantzig2(T))    
    while FinDantzig2(T)==False:
        T, base = AuxSimplex(T, base)
        print(FinDantzig2(T))    

"""Dos Fases"""

def PrintTablaDosFases(T, variables, base):
    T = T.astype(float)
    artif = [f'R{i+1}' for i in range(T.shape[0] - 1)]
    variables = variables + artif
    #print(variables)
    row_format ="{:>8}|" * (len(variables + ['Sol']) + 1)
    header = row_format.format("", *variables + ['Sol'])
    separador = len(header)*'-'
    func_objetivo = row_format.format("Z", *T[-1])
    print(header)
    print(separador)
    
    for base, row in zip(base, T[0:-1]):
        print(row_format.format(base, *row))
        print(separador)
    
    print(func_objetivo)
    print(separador)
    print()

def DosFases(z, A, b):
    T = np.concatenate((A, np.identity(A.shape[0])), axis=1)
    z = np.concatenate((np.array([np.zeros(A.shape[1])]), np.array([np.ones(T.shape[0])])), axis=1)
    T = TabInit(z, T, b)
    variables = [f'x{i+1}' for i in range(A.shape[1])]
    T = T.astype(float)
    
    base = [f'R{i}' for i in np.arange(1, A.shape[0]+1)]    
    iterbase = [None]*(A.shape[0])

    PrintTablaDosFases(T, variables, base)
    
    for i in range(A.shape[0]):
        T = Pivot(T, i, A.shape[1] + i)
            
    PrintTablaDosFases(T, variables, base)
    
    while FinDantzig(T):
        i, j = Dantzig(T)
        T = Pivot(T, i, j)
        base[i] = "x{num}".format(num = j+1)
        iterbase[i] = j+1
        #print(base[i])
        PrintTablaDosFases(T, variables, base)
    
    for i in range(len(iterbase)):
        T = np.delete(T, A.shape[1], 1)
    
    return T, iterbase

def TtoAb(t):
    A = np.delete(T, -1, 1)
    A = np.delete(A, -1, 0)
    b = np.array([T[:, -1]])
    b = np.delete(b, -1, 1)
    
    return A, b

"""---Parametros---
    Para un problema en su forma estandar
    z: objetivo
    A: leq
    b: lrig
    base: indices de la base
    """
"""z = np.array([[-5, -4, 0, 0, 0, 0]])

A = np.array([[6, 4, 1, 0, 0, 0],
              [1, 2, 0, 1 ,0 ,0],
              [-1, 1, 0, 0, 1 ,0],
              [0, 1, 0, 0, 0, 1]])

b = np.array([[24, 6, 1, 2]])

base = [3, 4, 5, 6]

T = TabInit(z, A, b)

T.astype(float)

#print("-----------prueba func simplex")
Simplex(z, A, b, base)"""

#print("-----------prueba dos fases")
#k = DosFases(z, A, b)
#Simplex(z, A, b, k)

"""-----Ejercicio 1-----"""
print("---------------------")
print("-----Ejercicio 1-----")
print("---------------------")

z = np.array([[4, 1, 1]])

A = np.array([[2, 1, 2],
              [3, 3, 1]])

b = np.array([[4, 3]])

base = [2, 3]

T = TabInit(z, A, b)

T, k = DosFases(z, A, b)
A, b = TtoAb(T)
T = TabInit(z, A, b)
Simplex(z, A, b, k)

"""-----Ejercicio 2-----"""
print("---------------------")
print("-----Ejercicio 2-----")
print("---------------------")

z = np.array([[-3, -2, 0, 0, 0]])

A = np.array([[-1, 3, 1, 0, 0],
              [1, 1, 0, 1, 0],
              [2, -1, 0, 0, 1]])

b = np.array([[12, 8, 10]])

T = TabInit(z, A, b)

T, k = DosFases(z, A, b)
A, b = TtoAb(T)
print(A)
print(b)
T = TabInit(z, A, b)
Simplex(z, A, b, k)

"""-----Ejercicio 3-----"""
print("---------------------")
print("-----Ejercicio 3-----")
print("---------------------")

z = np.array([[-6, -8, -3, -9, 0, 0]])

A = np.array([[2, 1, 1, 3, 1, 0],
              [1, 3, 1, 2, 0, 1]])

b = np.array([[5, 3]])

#base = [2, 3]

T = TabInit(z, A, b)

T, k = DosFases(z, A, b)
A, b = TtoAb(T)
T = TabInit(z, A, b)
Simplex(z, A, b, k)

print("---------------------")
print("-----Ejercicio 3 - ejercicio franklin-----")
print("---------------------")

z = np.array([[-6, -8, -3, -9, 0, 0]])

A = np.array([[2, 1, 1, 3, 1, 0],
              [1, 3, 1, 2, 0, 1]])

b = np.array([[5, 3]])

base = [1, 2]

T = TabInit(z, A, b)

Simplex(z, A, b, base)