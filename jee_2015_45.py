import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import sympy as sy
plt.rcParams.update({'font.size': 16})
from matplotlib.ticker import MultipleLocator

v1=np.array([0,0])
v2=np.array([0,1])

V=np.vstack((v1,v2)).T

u=np.array([-2,0])

a=-1*(2*u[0])/(4.0) #focal length and x co-ordinate of focus

#The intersection points of the parabola and latus rectum are of the form (a,y)

#Get end points of Latus rectum
y=sy.symbols("y",real=True)
m=sy.Matrix([a,y])
u1=sy.Matrix(u)

sol=sy.solve(m.T*V*m+2*u1.T*m,y)

#r is the end of LR
r=np.array([a,float((sol[1])[0])])
xr=a*np.ones(100)
yr=np.linspace(-r[1],r[1],100)


fig, ax = plt.subplots()
plt.plot(xr,yr,color='green',label='Latus rectum')
#suppose equation of tangent is (n^T)x=c
n=np.matmul(r,V)+u


#Equation of normal is (m^T)x=z
inv=np.vstack(([0,1],[-1,0]))
m=np.matmul(inv,n)

m=m/(abs(m[0])) #normalize
#The normal also passes thorugh r=(1,2). So, z=(m^T)r
z=np.matmul(m.T,r)

#Now, the normal is a tangent to the circle. The tangent to circle at point p is
#((p^T)B+w^T)x+(p^T)w+F=0 with B=I

#using the equation of the circle, x.Tx-2(3,-2)=r^2-(3,-2)((3,-2).T)
#Here, F=-r^2+(3,-2)((3,-2).T)
w=np.array([-3,2])

#Comparing the above equation with the equation of normal,
#(p.T)I+w.T=m.T
#p.Tw+F=-z
p=m-w

rsq=np.matmul(p,w.T)+np.matmul(w,w.T)+z
r=np.sqrt(rsq)

print(rsq)
print(r)

x1=[]
y1=[]



xl=np.linspace(0,10,100)

y=sy.symbols("y",real=True)



y1l=[]

for i in range(0,len(xl)):
    q=sy.Matrix([xl[i],y])
    m1=sy.Matrix(m)
    
    sol=sy.solve((m1.T*q)[0]-z,y)
    
    y1l.append(float(sol[0]))

plt.plot(xl,y1l,color='black',label="Normal to parabola")
x=np.linspace(0,5,100)
for i in range(0,len(x)):
    y=sy.symbols("y",real=True)
    m=sy.Matrix([x[i],y])
    u1=sy.Matrix(u)

    sol=sy.solve(m.T*V*m+2*u1.T*m,y)
    
    x1.append(x[i])
    #x1.append(x[i])
    y1.append(float((sol[0])[0]))
    #y1.append((sol[1])[0])
    


xc=np.linspace(3.0-np.sqrt(2)-0.1,3.0+np.sqrt(2)+0.1,300)
x1c=[]
y1c=[]
x2c=[]
y2c=[]

y=sy.symbols("y",real=True)

for i in range(0,len(xc)):
    
    q=sy.Matrix([xc[i],y])
    u1=sy.Matrix(np.array([3,-2]))

    sol=sy.solve((q.T*q)[0]-(2*u1.T*q)[0]-rsq+np.matmul(w,w.T),y)
    
    
    if(len(sol)==0):
        continue
    x1c.append(xc[i])
    
    y1c.append(float(sol[0]))
    if(len(sol)>1):
        x2c.append(xc[i])
        y2c.append(float(sol[1]))




plt.plot(x1,y1,color='blue',label="parabola")
plt.plot(x1c,y1c,color="red",label="circle")
plt.plot(x2c,y2c,color="red")

plt.plot(x1,-np.array(y1),color='blue')
plt.legend(loc='upper left')
plt.axis("equal")
plt.xlabel('x')
plt.ylabel('y')
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(MultipleLocator(1))
plt.grid()
plt.show()
