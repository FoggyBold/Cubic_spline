import matplotlib.pyplot as plt
import numpy as np
import random as rn
import math
 
def function (x):
  return x**3+1               

class SplineTuple:
    def __init__(self, a, b, c, d, x):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.x = x

def BuildSpline(x, y, n,a,b):
    splines = [SplineTuple(0, 0, 0, 0, 0) for _ in range(0, len(x))]
    for i in range(0, len(x)):
        splines[i].x = x[i]
        splines[i].a = y[i]
        
    splines[0].b = a
    splines[n - 1].b = b

    alpha = [0.0 for _ in range(0, n - 1)]
    beta  = [0.0 for _ in range(0, n - 1)]
 
    for i in range(1, n - 1):
        hi  = -x[i] + x[i-1] 
        hi1 = -x[i+1] + x[i] 
        A = 1/hi #при bi-1
        C = 4.0/hi1 + 2/hi #при bi
        B = 2/hi1 #при bi+1
        F = -6.0 * (y[i+1] -y[i]) / (hi1*hi1) -3.0 * (y[i]-y[i-1])/(hi*hi) #значение после равно
        z = (A * alpha[i - 1] + C) #знаменатель
        alpha[i] = -B / z
        beta[i] = (F - A * beta[i - 1]) / z
    
    for i in range(n - 2, 1, -1):
        splines[i].b = alpha[i] * splines[i + 1].b + beta[i]
    
    # По известным коэффициентам b[i] находим значения d[i] и c[i]
    for i in range(1, n):
        hi = - x[i] + x[i - 1]
        splines[i].d = 2*(y[i]-y[i-1])/(hi**3)+(splines[i-1].b+splines[i].b)/(hi*hi)
    for i in range(1, n):
        hi = - x[i] + x[i - 1]
        splines[i].c = (splines[i-1].b-splines[i].b)/(2*hi)-3*splines[i].d*hi/2
    return splines
    
def printSpline(spline,x,i):
    return spline[i].a + spline[i].b*(x-spline[i].x)+spline[i].c*(x-spline[i].x)**2+spline[i].d*(x-spline[i].x)**3

def xx_in(x,xx,n):
  i=1
  while ((i<n) and (xx>x[i])):
    i+=1
  return i

#равномерное разбиение 
def uniform_partition(N,a,b): 
  x=np.zeros(N)
  for i in range(N):
    x[i] = a+(b-a)/(len(x)-1)*i
  return x


#разбиение Чебышева
def chebyshev_partition(N,a,b): 
  x=np.zeros(N)
  for i in range(N):
    x[i] = (a+b)/2 - (b-a)/2*math.cos((2*i+1)/(2*N+2)*math.pi)
  return x

def aproved(spline):
  for i in range(1,len(spline)):
    if(spline[i].d == 0):
      return False
  return True

def printSpl(new_x, y_new, x, y, spline_n):
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), dpi=380)
    axes.plot(new_x, y_new, 'red', label='f(x)',linestyle='dashed')
    axes.plot(x, y, '*', label='Узлы интерполяции') 
    axes.plot(new_x, spline_n, label='Сплайн')
    axes.legend(loc='upper left');
    axes.grid()
    fig.tight_layout()
    plt.show()

def countError(y_new, spline_n):
    error = 0
    for i in range(len(new_x)):
        if (abs(y_new[i]-spline_n[i]) > error):
            error = abs(y_new[i]-spline_n[i])
    return error

A=1
B=4
a=0
b=1
n=20
new_n=5 #число для гладкого построения ф-ии
n+=1
x = chebyshev_partition(n,a,b+(b-a)/n)
#print(x)
y = np.arange(a,b,(b-a)/n)
for i in range(len(x)):
  y[i]=function(x[i])
new_x = chebyshev_partition(n*new_n,a,b+(b-a)/(n*new_n))
#print(new_x)
spline = BuildSpline(x, y, len(x),A,B)
spline_n = np.arange(a,b,(b-a)/(new_n*n))

k=0
for i in range(len(new_x)):
  k=xx_in(x,new_x[i],n)
  spline_n[i]=printSpline(spline,new_x[i],k)

y_new=np.arange(a,b,(b-a)/(new_n*n))
for i in range (len(y_new)):
  y_new[i]=function(new_x[i])

printSpl(new_x, y_new, x, y, spline_n)
#print(aproved(spline))
print(countError(y_new, spline_n))    