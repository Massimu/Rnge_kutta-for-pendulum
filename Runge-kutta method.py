#A program to study the driven pendulum under damping via the fourth-order Runge-Kutta algorithm.
#modules :
import math as m
import matplotlib.pyplot as plt

#Method to provide the generalized velocity vector- define g(y,t) for main algorithm
def g(y,t):
#    import math as m // when we want to use it separately
    l=len(y)
    q=0.5
    b=1.15
    omeg=2.0/3.0
    v=[{}]*l
    v[0]=y[1]
    v[1]=- m.sin(y[0])+b*m.cos(omeg*t)-q*y[1]
    return v

#Method to complete one Runge-Kutta step.
def RK(y,t,dt):
#    import math as m // when we want to use it separately
    l=len(y)
    c1=[{}]*l
    c2=[{}]*l
    c3=[{}]*l
    c4=[{}]*l
    
    c1 = g(y,t)
    for i in range(l):
         c2[i] = y[i]+dt*c1[i]/2
    c2 = g(c2,t+dt/2)
    for i in range(l):
         c3[i] = y[i]+dt*c2[i]/2
    c3 = g(c3,t+dt/2)
    for i in range(l):
         c4[i] = y[i] + dt*c3[i]
    c4 = g(c4, t+dt)
    for i in range(l):
         c1[i] = y[i] + dt*(c1[i]+2*(c2[i]+c3[i])+c4[i])/6
    return c1  
     
#starting point 
#variables
n=1000
nt=10
s=5
y1=[{}]*(n+1)
y2=[{}]*(n+1)
y=[{}]*2
   
#Set up time step and initial values
dt=(3*m.pi)/nt
y1[0]=y[0]=0
y2[0]=y[1]=2

#Perform the 4th-order Runge-Kutta integration
for i in range(n):
    t=dt*i
    y=RK(y,t,dt)
    y1[i+1] = y[0]
    y2[i+1] = y[1]
    np=int(y1[i+1]/(2*m.pi)+0.5)
    y1[i+1] -= 2*m.pi*np
plt.scatter(y2,y1)
plt.xlabel("$\\theta$")
plt.ylabel("$\omega$")
plt.grid()
plt.axhline(0)
plt.axvline(0)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    