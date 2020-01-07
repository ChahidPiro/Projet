#!/usr/bin/env python
# coding: utf-8

# # Travail personnel MATHF314 par Chahid Pironnet

# On importe les modules nécessaires.

# In[12]:


import numpy as np
import scipy
from scipy import linalg
from matplotlib import pyplot as plt
import random
from math import *
from numpy.linalg import inv


# On aura besoin de ça dans la suite.

# In[13]:


N=101
L=1.0
y=numpy.linspace(-L,L,N)
dy=2*L/(N-1)


# ### Question 1

# En y=-1, grâce à https://en.wikipedia.org/wiki/Finite_difference_coefficient#Forward_finite_difference et aux conditions aux bords, on a que v'(0)=0=-3/2v(0)+2v(1)-1/2v(2) et comme v(0)=0 que v(1)=1/4v(2).
# 
#     En y=1, grâce à https://en.wikipedia.org/wiki/Finite_difference_coefficient#Backward_finite_difference et aux conditions aux bords, on a que v'(N-1)=0=1/2v(N-3)-2v(N-2)+3/2v(N-1) et comme v(N-1) = 0 que v(N-2)=1/4v(N-3)

# In[14]:


v=numpy.empty(N-4)


# ### Question 2

#     On a complété la routine afin qu'elle retourne une discrétisation de l'opérateur $D^2$ : en effet, on obtient une matrice nulle partout sauf sur les 3 diagonales centrales, composée respectivement des éléments [1, -2, 1] obtenus à https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference mais en prenant en compte les conditions aux bords à la première entrée et dernière de la matrice où on trouve grâce à l'exercice précédent le coefficient de -7/4.

# In[15]:


def D2_v(N,dy) :
    D2=np.zeros((N-4,N-4))
    for i in range (N-4):
        if i == 0 :                      #première ligne
            D2[0][0]=(-7/4)/(dy**2)
            D2[0][1]=(1)/(dy**2)
        elif i == N-5 :                  #dernière ligne
            D2[N-5][N-5]=(-7/4)/(dy**2)
            D2[N-5][N-6]=1/(dy**2)
        else :                           #les 3 diagonales centrales
            D2[i][i]=(-2)/(dy**2)
            D2[i][i-1]=(1)/(dy**2)
            D2[i][i+1]=(1)/(dy**2)
    return D2


# ### Question 3

# On a complété la routine afin qu'elle retourne une discrétisation de l'opérateur $D^4$ : en effet, on obtient une matrice nulle partout sauf sur les 5 diagonales centrales, composée respectivement des éléments [1, -4, 6, -4, 1] obtenus à https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference mais en prenant en compte les conditions aux bords pour deux éléments des 2 premières lignes et des deux dernières lignes de la matrice où on trouve grâce à l'exercice 1 les coefficient de 5, -15/4, -15/4, 5.

# In[18]:


def D4_v(N,dy) :
    D4=np.zeros((N-4,N-4))
    for i in range (N-4):
        if i == 0 : #première ligne
            D4[0][0]=(5)/(dy**4)
            D4[0][1]=(-4)/(dy**4)
            D4[0][2]=(1)/(dy**4)
        elif i == 1 : #deuxième ligne
            D4[1][0]=(-15/4)/(dy**4)
            D4[1][1]=(6)/(dy**4)
            D4[1][2]=(-4)/(dy**4)
            D4[1][3]=(1)/(dy**4)
        elif i == N-5 : #dernière ligne
            D4[N-5][N-5]=(5)/(dy**4)
            D4[N-5][N-6]=(-4)/(dy**4)
            D4[N-5][N-7]=(1)/(dy**4)
        elif i == N-6 : # avant dernière ligne
            D4[N-6][N-5]=(-15/4)/(dy**4)
            D4[N-6][N-6]=(6)/(dy**4)
            D4[N-6][N-7]=(-4)/(dy**4)
            D4[N-6][N-8]=(1)/(dy**4)
        else :
            D4[i][i]=(6)/(dy**4)
            D4[i][i-1]=(-4)/(dy**4)
            D4[i][i+1]=(-4)/(dy**4)
            D4[i][i-2]=(1)/(dy**4)
            D4[i][i+2]=(1)/(dy**4)
    return D4


# ### Question 4

# Voici l'opérateur L définit comme suit :

# In[17]:


def L_v(N,y,dy,R,alpha):
   U = np.matrix(np.zeros((N-4,N-4))) 
   for i in range(N-4): 
    U[i,i]=(1-((y[i+2])**2)) 
   Us=np.matrix((-2)*np.identity(N-4))
   D2=np.matrix(D2_v(N,dy))
   D4=np.matrix(D4_v(N,dy))
   Alpha = np.matrix(alpha * np.identity(N-4))
   L = inv(np.matrix(D2-1*Alpha**2))*(-1j*alpha*(U)*(D2-Alpha**2)+1j*alpha*(Us)+1/R*(D4-2*(Alpha**2)*D2+Alpha**4))
   return L
R=500
alpha=0.3
print(L_v(N,y,dy,R,alpha))


# Un problème récurrant dans la suite du travail me contraint d'arrêter ici.

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




