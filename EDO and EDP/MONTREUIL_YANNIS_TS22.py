

import numpy as np
import matplotlib.pyplot as plt
import math

#Exercice Non Linear

""" ODE function """
def f(t,y):
    return 1+y*y

""" Euler Algorithm"""
def Euler(y0,a,b,nb):
    h=(b-a)/nb
    X=[a]
    Y=[y0]
    for i in range(1,nb):

        y0=Y[i-1]
        t0=X[i-1]
        t1=t0+h
        y1=y0+h*f(t0,y0)
        Y.append(y1)
        X.append(t1)
    return X,Y

def Runge2(y0,a,b,nb):
    h=(b-a)/nb
    t=a
    X=[a]
    Y=[y0]
    for i in range(1,nb):
        
        t=X[i-1]+h
        k1=Y[i-1]+(h/2)*f(X[i-1],Y[i-1])
        k2=f(X[i-1]+(h/2),k1)
        y2=Y[i-1]+h*k2
        
        Y.append(y2)
        X.append(t)
    return X,Y


def Runge4(y0,a,b,nb):
    h=(b-a)/nb
    t=a
    X=[a]
    Y=[y0]
    for i in range(1,nb):
        
        t=X[i-1]+h
        k1=f(X[i-1],Y[i-1])
        k2=f(X[i-1]+(h/2),k1*(h/2)+Y[i-1])
        k3=f(X[i-1]+(h/2),k2*(h/2)+Y[i-1])
        k4=f(X[i-1]+h,Y[i-1]+h*k3)
        y2=Y[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        
        Y.append(y2)
        X.append(t)
    return X,Y

def Adam(y0,a,b,nb):
    h=(b-a)/nb
    t1=a
    X=[a]
    Y=[y0]
    y1=y0
    #On doit calculer avec une méthode explicite la première itération
    y2=y1+h*f(X[0],Y[0])
    Y.append(y2)
    X.append(a+h)
    for i in range(2,nb):
        
        t=X[i-1]+h
        
        y3=Y[i-1]+(3*h/2)*f(X[i-1],Y[i-1])-(1/2)*h*f(X[i-2],Y[i-2])
        Y.append(y3)
        X.append(t)
    return X,Y

#Initialisation des paramètres
y0= 0.0
a = 0.0
b = 1.0
niter = 150

#Euler
x,y = Euler(y0,a,b,niter)
plt.plot(x,y,'r', label=r"Euler method", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)

#RK2
x,y = Runge2(y0,a,b,niter)
plt.plot(x,y,'b', label=r"Runge Kutta2", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)

#RK4
x,y = Runge4(y0,a,b,niter)
plt.plot(x,y,'g', label=r"Runge Kutta4", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)

#Adams
x_a,y_a = Adam(y0,a,b,niter)
plt.plot(x_a,y_a,'k', label=r"Adams", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)
plt.title("Comparaison des méthodes Equation Linéaire")
plt.grid()
plt.show()

    

#-----------Exercice Numerical Instability--------------

def Runge2_f(y0,a,b,nb,f): #Fonction avec une fonction en paramètre
    h=(b-a)/nb
    t=a
    X=[a]
    Y=[y0]
    y1=y0
    for i in range(1,nb):
        
     
        t=X[i-1]
        y1=Y[i-1]
        k1=y1+(h/2)*f(t,y1)
        k2=f(t+(h/2),k1)
        y2=y1+h*k2
        
        Y.append(y2)
        X.append(t+h)
        
    return X,Y

""" ODE function """
def f_1(t,y):
    return 3*y-4*math.exp(-t)

y0= 1.0
a = 0.0
b = 10
niter = 100

x_num,y_num = Runge2_f(y0,a,b,niter,f_1)
plt.plot(x_num,y_num,'r', label=r"Runge", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)

y_th=[y0]
t_th=np.linspace(0,10,niter)
for j in range(1,len(t_th)):
    y_th.append(math.exp(-t_th[j]))

plt.plot(t_th,y_th,'g', label=r"Théo", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)
plt.title("Comparaison théorique/numérique exercice non linear 2")
plt.grid()
plt.show()


"""On observe que la méthode ne converge pas, en effet il nous faut augmenter fortement le nombre d'itération pour un même intervalle
pour que cette méthode converge (meme chose pour RK4)
L'erreur peut aussi venir de la précision de nos flottants qui peuvent être trop faible"""

error=[]
for m_p in range(0,len(y_th)):
    error.append(math.exp(-x_num[m_p])-y_num[m_p])
    
plt.semilogy(x_num,error,'g', label=r"Théo", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)
plt.title("Comparaison théorique/numérique exercice non linear 2")
plt.grid()
plt.show()
#-------------- Exercice Non Linear 2------------------------


""" Euler Vector Algorithm"""

def F(t,V,m,l,g,A):
    L=[V[1], (A*math.sin(t)/(m*l*l))-(g/l)*math.sin(V[0])]
    return np.array(L)

def Euler_1(V0,a,b,nb,m,l,g,A):
    h=(b-a)/nb
    t=a
    X=[a]
    Y=[V0[0]]
    V1=V0
    for i in range(1,nb):
        
        V2=V1+h*F(t,V1,m,l,g,A)
        t=t+h
        Y.append(V2[0])
        X.append(t)
        V1=V2
    return X,Y

def Runge2_1(V0,a,b,nb,m,l,g,A):
    h=(b-a)/nb
    t=a
    X=[a]
    Y=[V0[0]]
    V1=V0
    for i in range(1,nb):
        
        k1=V1+(h/2)*F(t,V1,m,l,g,A)
        k2=F(t+h/2,V1,m,l,g,A)
        V2=V1+h*F(t+h/2,V1+(h/2)*F(t,V1,m,l,g,A),m,l,g,A)
        t=t+h
        Y.append(V2[0])
        X.append(t)
        V1=V2
    return X,Y

def Runge4_1(V0,a,b,nb,m,l,g,A):
    h=(b-a)/nb
    t=a
    X=[a]
    Y=[V0[0]]
    V1=V0
    for i in range(1,nb):
        
        k1=F(t,V1,m,l,g,A)
        k2=F(t+h/2,V1+(h/2)*k1,m,l,g,A)
        k3=F(t+h/2,V1+(h/2)*k2,m,l,g,A)
        k4=F(t+h,V1+h*k2,m,l,g,A)
        
        V2=V1+(h/6)*(k1+2*k2+2*k3+k4)
        t=t+h
        Y.append(V2[0])
        X.append(t)
        V1=V2
    return X,Y

#Initialisation des paramètres
teta0=0
teta_p0=1
T=20
m=0.1
l=0.3
A=0.1
g=9.81

#Initialisation du Vecteur
V0=[teta0,teta_p0]

a = 0.0
b = T
niter = 1000


#Tracer des données
x,y = Euler_1(V0,a,b,niter,m,l,g,A)
plt.plot(x,y,'k', label=r"Pendule_Euler", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)

x,y = Runge2_1(V0,a,b,niter,m,l,g,A)
plt.plot(x,y,'r', label=r"Pendule_Runge2", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)

x,y = Runge4_1(V0,a,b,niter,m,l,g,A)
plt.plot(x,y,'b', label=r"Pendule_Runge4", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)

plt.ylim(-1,1)
plt.xlim(0,20)
plt.title("Comparaison des méthodes pendule")
plt.grid()
plt.show()

""" En augmentant le nombre d'itération, on peut permettre à la méthode d'Euler de converger"""


#----------------Exercice Mass Spring -------------------

""" Euler Vector Algorithm"""

def F2(t,V,m,l,g,mu,k,V_p):
    if V_p<0:
        L=[V[1], mu*g-(k/m)*V[0]]
    else:
        L=[V[1], -mu*g-(k/m)*V[0]]   
    return np.array(L)

def Runge2_2(V0,a,b,nb,m,l,g,mu,k):
    h=(b-a)/nb
    t=a
    X=[a]
    Y=[V0[0]]
    V_p=V0[0]
    V1=V0
    for i in range(1,nb):
        
        k1=V1+(h/2)*F2(t,V1,m,l,g,mu,k,V_p)
        k2=F2(t+h/2,V1,m,l,g,mu,k,V_p)
        
        V2=V1+h*F2(t+h/2,V1+(h/2)*F2(t,V1,m,l,g,mu,k,V_p),m,l,g,mu,k,V_p)
        t=t+h
        Y.append(V2[0])
        X.append(t)
        V1=V2
        V_p=V2[1]
    return X,Y

#Initialisation des paramètres
Y_0=0.1
Y_p0=0
k=3000
mu=0.5
m=6.0
V0=[Y_0,Y_p0]
a = 0.0
b = 5
niter = 10000

xms,yms = Runge2_2(V0,a,b,niter,m,l,g,mu,k)
plt.plot(xms,yms,'b', label=r"Mass Spring", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)
plt.title("Déplacement en fonction du temps - Mass Spring")
plt.grid()
plt.show()
#On trouve que le prochain pic se trouve à y(562)=0,06076 or nous avions bien y0-4mu*m*g/k=0,06076



#------- Iron Block-------------------------

""" ODE function """

def F3(t,V,c,m,k):
    L=[V[1], -(k/m)*V[0]+c/(m*V[0]*V[0])]
    return np.array(L)

def Runge2_3(V0,a,b,nb,c,m,k):
    h=(b-a)/nb
    t=a
    T=[a]
    X=[V0[0]]
    V1=V0
    for i in range(1,nb):
        
        k1=V1+(h/2)*F3(t,V1,c,m,k)
        k2=F3(t+h/2,V1,c,m,k)
        
        V2=V1+h*F3(t+h/2,V1+(h/2)*F3(t,V1,c,m,k),c,m,k)
        t=t+h
        X.append(V2[0])
        T.append(t)
        V1=V2
    return T,X


#Initialisation 
k=120
m=1.0
L=0.2
c=5.0
x0=L
v_p0=0
V0=[x0,v_p0]


a = 0.0
b = 1
niter = 10000

#Traçage des données
t,x = Runge2_3(V0,a,b,niter,c,m,k)
plt.plot(t,x,'b', label=r"Iron block", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)
plt.title("Déplacement en fonction du temps - Iron Block")
plt.grid()
plt.show()


#Amplitude=0.35
Amp=max(x)-min(x)

#Période : On trouve T=0.33s



#---------- Magnus--------------------


""" Fonction pour les vitesses """

def F4x(t,vx,vy,vz,k,P,W):
    return -k*math.sqrt(vx*vx+vy*vy+vz*vz)*vx-P*W[2]*vy

def F4y(t,vx,vy,vz,k,P,W,m,g):
    return -k*math.sqrt(vx*vx+vy*vy+vz*vz)*vy+P*W[2]*vz-m*g

def F4z(t,vx,vy,vz,k,P,W,m,g):
    return -k*math.sqrt(vx*vx+vy*vy+vz*vz)*vz

""" Fonction pour les déplacements """

def f4(t,vx):
    return vx

#RK4
    
def Runge4_mag(V0,a,b,nb,k,P,W,m,g,x0):
    
    #Initialisation des paramètres
    h=(b-a)/nb
    t=a
    T=[a]
    Vx=[V0[0]]
    Vy=[V0[1]]
    Vz=[V0[2]]
    x=[x0[0]]
    y=[x0[0]]
    z=[x0[0]]
    
    #On va résoudre axe par axe, pour ensuite les réutiliser en incrémentant
    
    for i in range(1,nb):
        
        t=T[i-1]+h
        
        k1=F4x(T[i-1],Vx[i-1],Vy[i-1],Vz[i-1],k,P,W)
        k2=F4x(T[i-1]+(h/2),k1*(h/2)+Vx[i-1],Vy[i-1],Vz[i-1],k,P,W)
        k3=F4x(T[i-1]+(h/2),k2*(h/2)+Vx[i-1],Vy[i-1],Vz[i-1],k,P,W)
        k4=F4x(T[i-1]+h,Vx[i-1]+h*k3,Vy[i-1],Vz[i-1],k,P,W)
        Vx2=Vx[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        
        
        k1=F4y(T[i-1],Vx[i-1],Vy[i-1],Vz[i-1],k,P,W,m,g)
        k2=F4y(T[i-1]+(h/2),Vx[i-1],k1*(h/2)+Vy[i-1],Vz[i-1],k,P,W,m,g)
        k3=F4y(T[i-1]+(h/2),Vx[i-1],k2*(h/2)+Vy[i-1],Vz[i-1],k,P,W,m,g)
        k4=F4y(T[i-1]+h,Vx[i-1],Vy[i-1]+h*k3,Vz[i-1],k,P,W,m,g)
        Vy2=Vy[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        
        
        k1=F4z(T[i-1],Vx[i-1],Vy[i-1],Vz[i-1],k,P,W,m,g)
        k2=F4z(T[i-1]+(h/2),Vx[i-1],Vy[i-1],k1*(h/2)+Vz[i-1],k,P,W,m,g)
        k3=F4z(T[i-1]+(h/2),Vx[i-1],Vy[i-1],k2*(h/2)+Vz[i-1],k,P,W,m,g)
        k4=F4z(T[i-1]+h,Vx[i-1],Vy[i-1],Vz[i-1]+h*k3,k,P,W,m,g)
        Vz2=Vz[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        
        Vz.append(Vz2)
        Vx.append(Vx2)
        Vy.append(Vy2)
        T.append(t)
        
        #on détermine la position en réalisant à nouveau RK4
        
        k1=f4(T[i-1],Vx[i-1])
        k2=f4(T[i-1]+(h/2),k1*(h/2)+Vx[i-1])
        k3=f4(T[i-1]+(h/2),k2*(h/2)+Vx[i-1])
        k4=f4(T[i-1]+h,Vx[i-1]+h*k3)
        x2=x[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        x.append(x2)
        
        k1=f4(T[i-1],Vy[i-1])
        k2=f4(T[i-1]+(h/2),k1*(h/2)+Vy[i-1])
        k3=f4(T[i-1]+(h/2),k2*(h/2)+Vy[i-1])
        k4=f4(T[i-1]+h,Vy[i-1]+h*k3)
        y2=y[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        y.append(y2)
        
        k1=f4(T[i-1],Vz[i-1])
        k2=f4(T[i-1]+(h/2),k1*(h/2)+Vz[i-1])
        k3=f4(T[i-1]+(h/2),k2*(h/2)+Vz[i-1])
        k4=f4(T[i-1]+h,Vz[i-1]+h*k3)
        z2=z[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        z.append(z2)
        
    return T,Vx,Vy,Vz,x,y,z

#Adams
    

def Adams4_mag(V0,a,b,nb,k,P,W,m,g,x0):
    
    #Initialisation des paramètres
    h=(b-a)/nb
    t=a
    T=[a,a+h]
    Vx=[V0[0]]
    Vy=[V0[1]]
    Vz=[V0[2]]
    x=[x0[0]]
    y=[x0[0]]
    z=[x0[0]]
    
    #On calcul le deuxième termes à l'aide de la méthode RK4 pour obtenir la meilleure précision et ainsi les comparer
    k1=F4x(T[0],Vx[0],Vy[0],Vz[0],k,P,W)
    k2=F4x(T[0]+(h/2),k1*(h/2)+Vx[0],Vy[0],Vz[0],k,P,W)
    k3=F4x(T[0]+(h/2),k2*(h/2)+Vx[0],Vy[0],Vz[0],k,P,W)
    k4=F4x(T[0]+h,Vx[0]+h*k3,Vy[0],Vz[0],k,P,W)
    Vx2=Vx[0]+(h/6)*(k1+2*k2+2*k3+k4)    
    Vx.append(Vx2)
    
    k1=F4y(T[0],Vx[0],Vy[0],Vz[0],k,P,W,m,g)
    k2=F4y(T[0]+(h/2),Vx[0],k1*(h/2)+Vy[0],Vz[0],k,P,W,m,g)
    k3=F4y(T[0]+(h/2),Vx[0],k2*(h/2)+Vy[0],Vz[0],k,P,W,m,g)
    k4=F4y(T[0]+h,Vx[0],Vy[0]+h*k3,Vz[0],k,P,W,m,g)
    Vy2=Vy[0]+(h/6)*(k1+2*k2+2*k3+k4)
    Vy.append(Vy2)
    
    k1=F4z(T[0],Vx[0],Vy[0],Vz[0],k,P,W,m,g)
    k2=F4z(T[0]+(h/2),Vx[0],Vy[0],k1*(h/2)+Vz[0],k,P,W,m,g)
    k3=F4z(T[0]+(h/2),Vx[0],Vy[0],k2*(h/2)+Vz[0],k,P,W,m,g)
    k4=F4z(T[0]+h,Vx[0],Vy[0],Vz[0]+h*k3,k,P,W,m,g)
    Vz2=Vz[0]+(h/6)*(k1+2*k2+2*k3+k4)
    Vz.append(Vz2)
    
    #Même chose pour les positions

    
    #On va résoudre axe par axe, pour ensuite les réutiliser en incrémentant
    
    for i in range(2,nb):
        
        t=T[i-1]+h
        Vx3=Vx[i-1]+(3*h/2)*F4x(T[i-1],Vx[i-1],Vy[i-1],Vz[i-1],k,P,W)-(1/2)*h*F4x(T[i-2],Vx[i-2],Vy[i-2],Vz[i-2],k,P,W)
        Vy3=Vy[i-1]+(3*h/2)*F4y(T[i-1],Vx[i-1],Vy[i-1],Vz[i-1],k,P,W,m,g)-(1/2)*h*F4y(T[i-2],Vx[i-2],Vy[i-2],Vz[i-2],k,P,W,m,g)
        Vz3=Vz[i-1]+(3*h/2)*F4z(T[i-1],Vx[i-1],Vy[i-1],Vz[i-1],k,P,W,m,g)-(1/2)*h*F4z(T[i-2],Vx[i-2],Vy[i-2],Vz[i-2],k,P,W,m,g)
    
        Vz.append(Vz3)
        Vx.append(Vx3)
        Vy.append(Vy3)
        T.append(t)
        
    return T,Vx,Vy,Vz


#Initialiation des paramètres
m=0.430
rho=1.2
R=0.11
alpha=(math.pi)/2
Cf=0.45
P=1/(2*m)*math.pi*rho*R**3*math.sin(alpha)
k=1/(2*m)*Cf*math.pi*R**2*rho
Vx0=0
Vy0=0
Vz0=0
V0=[Vx0,Vy0,Vy0]
Wx=0
Wy=0
Wz=10
W=[Wx,Wy,Wz]
X0=[0,0,0]

a = 0.0
b = 10
niter = 10000

#Comparaison des deux méthodes

t,Vx,Vy,Vz,x,y,z = Runge4_mag(V0,a,b,niter,k,P,W,m,g,X0)
plt.plot(t,Vx,'r', label=r"Vx RK4", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)

t,Vx,Vy,Vz = Adams4_mag(V0,a,b,niter,k,P,W,m,g,X0)
plt.plot(t,Vx,'b', label=r"Vx Adams", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)
plt.title("Comparaison des méthodes - MAgnus")
plt.grid()
plt.show()

#On observe que les deux méthodes semblent converger de la même manière

#On regarde la NOrme

OM=[]
for j in range(0,len(x)):
    OM.append(math.sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]))
    
plt.plot(t,OM,'m', label=r"Norme des vecteurs déplacements RK4", linestyle='--')
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 0.95),fancybox=True, shadow=True)
plt.title("Déplacement en fonction du temps - MAgnus")
plt.grid()

plt.show()
