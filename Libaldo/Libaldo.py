


from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image, display
init_printing()

from  libaldo_math2 import *
from libaldo_algorith import *


from  lib_MyEq import *
from  lib_MyEqEq import *

n,m,m1,m2,m3,m4,M,g,x,x1,x2,x3,x4,x5,y,y1,y2,X,Y,a,a1,a2,a3,v,v1,v2,M1,M2,M3, V ,V1 ,V2= symbols('n m m1 m2 m3 m4 M g x x1 x2 x3 x4 x5 y y1 y2 X Y a a1 a2 a3 v v1 v2 M1 M2 M3 V V1 V2 ')
w,w1,w2,aw,aw1,aw2,F,F1,F2,Rx,Ry,r,r1,r2,R,ax,ax1,ax2,ay,ay1,ay2= symbols('w w1 w2 aw aw1 aw2 F F1 F2 Rx Ry r r1 r2 R ax ax1 ax2 ay ay1 ay2')
mu,mu1,mu2,fr,fr1,fr2,f1, f2, f3,N1,N2,N3,Nm, L,L1,L2,h,h1,h2,b,H= symbols('mu  mu1 mu2  fr fr1 fr2 f1 f2 f3 N1 N2 N3 Nm  L L1 L2 h h1 h2 b H')
R1,R2,R3,rx,ry,fx1,fx2,fx3,fy1,fy2,fy3=symbols('R1 R2 R3 rx ry fx1 fx2 fx3 fy1 fy2 fy3')
alpha,tetha=symbols('alpha tetha')

at,an=symbols('a_t a_n')
vx,vy,In ,Fc,ac,aw =symbols('v_x v_y I_n F_c a_c a_w')
t,vxy,Po,Ti,I_n, a_c,a_t,  Io=symbols('t vxy Po Ti In  ac    at Io')
T,T1,T2,t1,t2,t3 =symbols('T T1 T2 t1 t2 t3')
K,k,X1,X2,X0,d,W,P=symbols('K k X1 X2 X0 d W P')
rho,z,z1,z2,A,p=symbols('rho z  z1 z2 A p') 
xo,yo,zo,beta,xi,xf=symbols('xo yo zo beta xi xf')
yp,xp,pp=symbols("y' x' p'") 

 
A,A1,A2,B,B1,B2,C,C1,C2,D,D1,D2,Q,Q1,Q2,S,S1,S2,Z,Z1,Z2=symbols('A A1 A2 B B1 B2 C C1 C2 D D1 D2 Q Q1 Q2 S S1 S2 Z Z1 Z2')
# diff variables

dm,ds,dx,dy,dz,dt,dr,dh,dL,da,dA,dv,dV,dM,=symbols('dm ds dx dy dz dt dr dh dL da dA dv dV dM')   

Vo,Lo,Xo,Yo,Zo=symbols('Vo Lo Xo Yo Zo')
Fx1,Fy1,To1=symbols('Fx1 Fy1 To1')

# constantes
A53=atan(frs(4,3)) # ang 53 = atan(4/3)
A37=atan(frs(3,4)) # ang 37 = atan(3/4)

 
