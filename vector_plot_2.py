from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from scipy.spatial import distance
import math
import seaborn as sns
import pandas as pd

a0=2.54

#Lower Hexagon
x_lh=[]
y_lh=[]
z_lh=[]
atom_lh=[]
cnt_atom_lh=0
with open("hexagon_data_AA_down_temp300_t0.xyz") as f:
	lines=f.readlines()

for line in lines:
	cnt_atom_lh+=1	
	atom_lh.append(line.split()[0])
	x_lh.append(float(line.split()[1]))
	y_lh.append(float(line.split()[2]))
	z_lh.append(float(line.split()[3]))


num_hexagon_l1=int(cnt_atom_lh/6)

print(num_hexagon_l1)

u_lh=[]
v_lh=[]

cnt=0
for i in range(num_hexagon_l1):
	u_lh.append(float(x_lh[cnt+1]-x_lh[cnt]))
	u_lh.append(float(x_lh[cnt+2]-x_lh[cnt+1]))
	u_lh.append(float(x_lh[cnt+3]-x_lh[cnt+2]))
	u_lh.append(float(x_lh[cnt+4]-x_lh[cnt+3]))
	u_lh.append(float(x_lh[cnt+5]-x_lh[cnt+4]))
	u_lh.append(float(x_lh[cnt+1]-x_lh[cnt+5]))
	cnt=cnt+6
	#print(i,cnt,len(u_lh),u_lh)

cnt=0
for i in range(num_hexagon_l1):
	v_lh.append(float(y_lh[cnt+1]-y_lh[cnt]))
	v_lh.append(float(y_lh[cnt+2]-y_lh[cnt+1]))
	v_lh.append(float(y_lh[cnt+3]-y_lh[cnt+2]))
	v_lh.append(float(y_lh[cnt+4]-y_lh[cnt+3]))
	v_lh.append(float(y_lh[cnt+5]-y_lh[cnt+4]))
	v_lh.append(float(y_lh[cnt+1]-y_lh[cnt+5]))
	cnt=cnt+6
	#print(i,cnt,len(v_lh),v_lh)


X1,Y1 = x_lh,y_lh
U1,V1 = u_lh,v_lh


fig1,ax1=plt.subplots()
ax1.set_title('Order_Parameter')
M1=np.hypot(U1,V1)
Q1=ax1.quiver(X1,Y1,U1,V1,M1)
qk=ax1.quiverkey (Q1, 1, 1, 1, r'order_param')

plt.savefig("h-BN_Moire_Layer1_vector_plot.png")
plt.show()

#print(cnt_atom_lh)
#print(atom_lh)
#print(x_lh)
#print(y_lh)
#print(z_lh)

#xy_lh=list(zip(x_lh,y_lh))
#print(xy_lh)

#Upper Hexagon
x_uh=[]
y_uh=[]
z_uh=[]
atom_uh=[]
cnt_atom_uh=0
with open("hexagon_data_AA_up_temp300_t0.xyz") as f:
	lines=f.readlines()

for line in lines:
	cnt_atom_uh+=1	
	atom_uh.append(line.split()[0])
	x_uh.append(float(line.split()[1]))
	y_uh.append(float(line.split()[2]))
	z_uh.append(float(line.split()[3]))


num_hexagon_l2=int(cnt_atom_uh/6)

print(num_hexagon_l2)

u_uh=[]
v_uh=[]

cnt=0
for i in range(num_hexagon_l2):
	u_uh.append(float(x_uh[cnt+1]-x_uh[cnt]))
	u_uh.append(float(x_uh[cnt+2]-x_uh[cnt+1]))
	u_uh.append(float(x_uh[cnt+3]-x_uh[cnt+2]))
	u_uh.append(float(x_uh[cnt+4]-x_uh[cnt+3]))
	u_uh.append(float(x_uh[cnt+5]-x_uh[cnt+4]))
	u_uh.append(float(x_uh[cnt+1]-x_uh[cnt+5]))
	cnt=cnt+6
	print(i,cnt,len(u_uh))

cnt=0
for i in range(num_hexagon_l2):
	v_uh.append(float(y_uh[cnt+1]-y_uh[cnt]))
	v_uh.append(float(y_uh[cnt+2]-y_uh[cnt+1]))
	v_uh.append(float(y_uh[cnt+3]-y_uh[cnt+2]))
	v_uh.append(float(y_uh[cnt+4]-y_uh[cnt+3]))
	v_uh.append(float(y_uh[cnt+5]-y_uh[cnt+4]))
	v_uh.append(float(y_uh[cnt+1]-y_uh[cnt+5]))
	cnt=cnt+6
	print(i,cnt,len(v_uh))

#print(cnt_atom_uh)
#print(atom_uh)
#print(x_uh)
#print(y_uh)
#print(z_uh)

print(len(x_uh),len(y_uh),len(u_uh),len(v_uh))
del x_uh[-5:-1] 
del y_uh[-5:-1]

X2,Y2 = x_uh,y_uh
U2,V2 = x_uh,x_uh


fig2,ax2=plt.subplots()
ax2.set_title('Order_Parameter')
M2=np.hypot(U2,V2)
Q2=ax2.quiver(X2,Y2,U2,V2,M2)
qk=ax2.quiverkey (Q2, 0.1, 0.1, 0.1, r'order_param')

plt.savefig("h-BN_Moire_Layer2_vector_plot.png")

plt.show()

x_ul=[]
y_ul=[]
u_ul=[]
v_ul=[]
for i in range(len(x_lh)):
	x_ul.append(x_lh[i])
	y_ul.append(y_lh[i])
	u_ul.append(u_lh[i])
	v_ul.append(v_lh[i])

for i in range(len(x_uh)):
	x_ul.append(x_uh[i])
	y_ul.append(y_uh[i])
	u_ul.append(x_uh[i])
	v_ul.append(x_uh[i])

print(len(x_lh),len(x_uh),len(x_ul))
print(len(y_lh),len(y_uh),len(y_ul))
print(len(u_lh),len(u_uh),len(u_ul))
print(len(v_lh),len(v_uh),len(v_ul))

X3,Y3 = x_ul,y_ul
U3,V3 = u_ul,v_ul


fig3,ax3=plt.subplots()
ax3.set_title('Order_Parameter')
M3=np.hypot(U3,V3)
Q3=ax3.quiver(X3,Y3,U3,V3,M3)
qk=ax3.quiverkey (Q3, 0.1, 0.1, 0.1, r'order_param')

plt.savefig("h-BN_Moire_Layer_both_vector_plot.png")

plt.show()


		

