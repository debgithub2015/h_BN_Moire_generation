from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math as m

a0=2.54


cnt_hexagon=0
hexagon_num=[]
rcomx_lower=[]
rcomy_lower=[]
rcomz_lower=[]
rcomx_upper=[]
rcomy_upper=[]
rcomz_upper=[]
color_index=[]
#with open("hexagon_data_AA_rcom_temp300_t0.dat") as f:
with open("initial_hexagon_com_test.xyz") as f:
	lines=f.readlines()

for line in lines:
	cnt_hexagon+=1
	hexagon_num.append(int(line.split()[0]))	
	rcomx_lower.append(float(line.split()[1]))
	rcomy_lower.append(float(line.split()[2]))
	rcomz_lower.append(float(line.split()[3]))
	rcomx_upper.append(float(line.split()[4]))
	rcomy_upper.append(float(line.split()[5]))
	rcomz_upper.append(float(line.split()[6]))


com_vec=[]

#Define the vector between the COM of upper and lower hexagon

for i in range(cnt_hexagon):
	com_vec_x=rcomx_upper[i]-rcomx_lower[i]
	com_vec_y=rcomy_upper[i]-rcomy_lower[i]
	com_vec_z=rcomz_upper[i]-rcomz_lower[i]
	com_vec.append([com_vec_x,com_vec_y,com_vec_z])

#print(com_vec)

com_dist=[]
for i in range(cnt_hexagon):
	com_dist.append(m.sqrt(com_vec[i][0]**2+com_vec[i][1]**2))
	if (com_dist[i]==0.):
		color_index[i]=0
#print(com_dist)

com_cosphi=[]
com_sinphi=[]
for i in range(cnt_hexagon):
	if (com_dist[i] != 0.) :
		com_cosphi.append(com_vec[i][0]/com_dist[i])
		com_sinphi.append(com_vec[i][1]/com_dist[i])

#print(i,com_cosphi,com_sinphi)

com_phi=[]

for i in range(cnt_hexagon):
	if (com_cosphi[i]>=0):
		if(com_sinphi[i]>= 0):
#Ist quadrant
			com_phi.append(m.acos(com_vec[i][0]/com_dist[i])*180./m.pi) 
		elif(com_sinphi[i]<0):
#IVth quadrant
			com_phi.append(-(m.acos(com_vec[i][0]/com_dist[i])*180./m.pi)) 
	elif(com_cosphi[i]<0):
		if(com_sinphi[i]>=0):
#IInd quadrant 
			com_phi.append(m.acos(com_vec[i][0]/com_dist[i])*180./m.pi) 
		elif(com_sinphi[i]<0):
#IIIrd quadrant
			com_phi.append(-(m.acos(com_vec[i][0]/com_dist[i])*180./m.pi))
	print(i,com_cosphi[i],com_sinphi[i],com_phi[i]) 
	
a0_limit=a0/m.sqrt(3.)

a0_norm=m.sqrt(3.)/2.

print(len(com_dist),len(com_phi))
for i in range(cnt_hexagon):
	print (i,com_dist[i],com_phi[i])

for i in range(cnt_hexagon):
#Case I:
	if(-15.<com_phi[i]<=0. or 0.<com_phi[i]<=15. or 105.<com_phi[i]<=135. or -75.<com_phi[i]<=-35.):
		if(0.<=com_dist[i]<=a0_limit):
			print(i,com_dist[i],com_phi[i])
			color_index.append(com_dist[i]/a0_limit)
			print(color_index[i])
#Case II:
	if(45.<com_phi[i]<=75. or 165.<com_phi[i]<=180. or -180.< com_phi[i]<=-165 or -135.<com_phi[i]<=-105.):
		if(0.<=com_dist[i]<=a0_limit):
			print(i,com_dist[i],com_phi[i])
			color_index.append(-com_dist[i]/a0_limit)
			print(color_index[i])
#Case III:
	if(15.<com_phi[i]<=45. or 135.<com_phi[i]<=165. or -105.<com_phi[i]<=-75.):
		if(0.<=(com_dist[i]*a0_norm) <=(a0/2.)):
			print(i,com_dist[i],com_phi[i])
			color_index.append((com_dist[i]*a0_norm)/a0_limit)
			print(color_index[i])
#Case IV:
	if(75.<com_phi[i]<=105. or -45.<com_phi[i]<=-15. or -165.<com_phi[i]<=-135.):
		if(0.<=(com_dist[i]*a0_norm)<=(a0/2.)):
			print(i,com_dist[i],com_phi[i])
			color_index.append(-(com_dist[i]*a0_norm)/a0_limit)
			print(color_index[i])


print(len(com_dist),len(com_phi),len(color_index))
for i in range(cnt_hexagon):
	print(i,com_dist[i],com_phi[i],color_index[i])

f1 = open("color_index_initial.txt", "w")
for i in range(cnt_hexagon):
	f1.write("%i %8.4f %8.4f %8.4f %8.4f %8.4f\n" %(hexagon_num[i],rcomx_lower[i],rcomy_lower[i],rcomx_upper[i],rcomy_upper[i],color_index[i]))
f1.close()

