import matplotlib.pyplot as plt
import numpy as np
import math as m
import statistics as stat
import glob


list_atom=[]
natom=674
maxatom=67
for i in range(maxatom):
	list_atom.append(m.floor(natom/maxatom)*i)

#print(list_atom)

def numeric_key(string):
    splitted = string.split('.')[0].split('_')
    if splitted[3].isdigit():
        return int(splitted[3])
    return -1


filename=[]
for file in glob.glob("config_nve_init_xyz*"):
	filename.append(file)

filename=sorted(filename,key=numeric_key)

#print(filename)

x_fix=[]
y_fix=[]
z_fix=[]
x_com=[]
y_com=[]
z_com=[]

for file in filename:
	myfile=open(file,"r")
	lines=myfile.readlines()
	x_sum=[]
	y_sum=[]
	z_sum=[]
	for i in range(natom):
		x_sum.append(float(lines[4336+i].split()[1]))
		y_sum.append(float(lines[4336+i].split()[2]))
		z_sum.append(float(lines[4336+i].split()[3]))

	for j in list_atom:
		x_fix.append(float(lines[4336+j].split()[1]))
		y_fix.append(float(lines[4336+j].split()[2]))
		z_fix.append(float(lines[4336+j].split()[3]))

	x_com.append(sum(x_sum)/natom)
	y_com.append(sum(y_sum)/natom)
	z_com.append(sum(z_sum)/natom)

x_fix=x_fix[:-maxatom]                
y_fix=y_fix[:-maxatom]                
z_fix=z_fix[:-maxatom]                

x_com=x_com[:-1]
y_com=y_com[:-1]
z_com=z_com[:-1]

#print(x_fix,len(x_fix))
#print(y_fix,len(y_fix))
#print(z_fix,len(z_fix))

#print(x_com,len(x_com))
#print(y_com,len(y_com))
#print(z_com,len(z_com))


filename=[]
for file in glob.glob("config_npt_xyz*"):
	filename.append(file)

filename=sorted(filename,key=numeric_key)
	
#print(filename)

for file in filename:
	myfile=open(file,"r")
	lines=myfile.readlines()
	x_sum=[]
	y_sum=[]
	z_sum=[]
	for i in range(natom):
		x_sum.append(float(lines[4336+i].split()[1]))
		y_sum.append(float(lines[4336+i].split()[2]))
		z_sum.append(float(lines[4336+i].split()[3]))

	for j in list_atom:
		x_fix.append(float(lines[4336+j].split()[1]))
		y_fix.append(float(lines[4336+j].split()[2]))
		z_fix.append(float(lines[4336+j].split()[3]))
	
	x_com.append(sum(x_sum)/natom)
	y_com.append(sum(y_sum)/natom)
	z_com.append(sum(z_sum)/natom)

x_fix=x_fix[:-maxatom]                
y_fix=y_fix[:-maxatom]                
z_fix=z_fix[:-maxatom]                

x_com=x_com[:-1]
y_com=y_com[:-1]
z_com=z_com[:-1]

#print(x_fix,len(x_fix))
#print(y_fix,len(y_fix))
#print(z_fix,len(z_fix))
	
#print(x_com,len(x_com))
#print(y_com,len(y_com))
#print(z_com,len(z_com))


filename=[]
for file in glob.glob("config_nve_fin_xyz*"):
        filename.append(file)

filename=sorted(filename,key=numeric_key)

#print(filename)

for file in filename:
	myfile=open(file,"r")
	lines=myfile.readlines()
	x_sum=[]
	y_sum=[]
	z_sum=[]
	for i in range(natom):
		x_sum.append(float(lines[4336+i].split()[1]))
		y_sum.append(float(lines[4336+i].split()[2]))
		z_sum.append(float(lines[4336+i].split()[3]))
	
	for j in list_atom:
		x_fix.append(float(lines[4336+j].split()[1]))
		y_fix.append(float(lines[4336+j].split()[2]))
		z_fix.append(float(lines[4336+j].split()[3]))
	
	x_com.append(sum(x_sum)/natom)
	y_com.append(sum(y_sum)/natom)
	z_com.append(sum(z_sum)/natom)

#print(x_fix,len(x_fix))
#print(y_fix,len(y_fix))
#print(z_fix,len(z_fix))
	
#print(x_com,len(x_com))
#print(y_com,len(y_com))
#print(z_com,len(z_com))

steps=int(len(x_fix)/maxatom)
print(steps)

vec_lay3=[]
for s in range(steps):
	for i in range(maxatom):
		#print(s,i,s*60+i,s*60)
		vec_lay3.append([(x_fix[(s*maxatom)+i]-x_com[s]),(y_fix[(s*maxatom)+i]-y_com[s]),(z_fix[(s*maxatom)+i]-z_com[s])])

vec_lay4=[]
for s in range(steps):
	for i in range(maxatom):
		#print(s,i,s*60+i,s*60)
		vec_lay4.append([(x_fix[(s*maxatom)+i]-x_com[0]),(y_fix[(s*maxatom)+i]-y_com[0]),(z_fix[(s*maxatom)+i]-z_com[0])])

		

#for s in range(steps):
#	for i in range(1,60):
#		#print(s,i,s*60+i,s*60)
#		vec_lay3.append([(x_fix[s*60+i]-x_fix[s*60]),(y_fix[s*60+i]-y_fix[s*60]),(z_fix[s*60+i]-z_fix[s*60])])

#print(vec_lay3,len(vec_lay3))
	
vec_lay3_cosphi=[]
vec_lay3_sinphi=[]
vec_lay3_dist=[]
for i in range(len(vec_lay3)):
	vec_lay3_dist.append(m.sqrt(vec_lay3[i][0]**2+vec_lay3[i][1]**2))

#print(vec_lay3_dist)

for i in range(len(vec_lay3)):
	vec_lay3_cosphi.append(vec_lay3[i][0]/vec_lay3_dist[i])
	vec_lay3_sinphi.append(vec_lay3[i][1]/vec_lay3_dist[i])

#print(vec_lay3_cosphi)
#print(vec_lay3_sinphi)


# Find angles for each vectors for each snapshot

vec_lay3_phi=[]
for i in range(len(vec_lay3)):
	if(vec_lay3_cosphi[i]>=0):
		if(vec_lay3_sinphi[i]>=0):
#Ist quadrant
			vec_lay3_phi.append(m.acos(vec_lay3[i][0]/vec_lay3_dist[i])*180./m.pi) 
		elif(vec_lay3_sinphi[i]<0):
#IVth quadrant
			vec_lay3_phi.append(-(m.acos(vec_lay3[i][0]/vec_lay3_dist[i])*180./m.pi)) 
	elif(vec_lay3_cosphi[i]<0):
		if(vec_lay3_sinphi[i]>=0):
#IInd quadrant 
			vec_lay3_phi.append(m.acos(vec_lay3[i][0]/vec_lay3_dist[i])*180./m.pi) 
		elif(vec_lay3_sinphi[i]<0):
#IIIrd quadrant
			vec_lay3_phi.append(-(m.acos(vec_lay3[i][0]/vec_lay3_dist[i])*180./m.pi))

#print(vec_lay3_phi,len(vec_lay3_phi))


mean_lay3_phi=[]
var_lay3_phi=[]
for s in range(steps):
	mean_lay3_phi.append(stat.mean(vec_lay3_phi[s*(maxatom):(s+1)*(maxatom)]))	
	var_lay3_phi.append(stat.variance(vec_lay3_phi[s*(maxatom):(s+1)*(maxatom)]))	
	
#print(mean_lay3_phi,len(mean_lay3_phi))
#print(var_lay3_phi,len(var_lay3_phi))

mean_diff_phi=[]
var_diff_phi=[]
time=[]
for i in range(len(mean_lay3_phi)-1):
	mean_diff_phi.append(mean_lay3_phi[i+1]-mean_lay3_phi[0])
	var_diff_phi.append(var_lay3_phi[i+1]-var_lay3_phi[0])
	time.append(i+1)

print(mean_diff_phi,len(mean_diff_phi))
print(var_diff_phi,len(var_diff_phi))
print(time)

fig=plt.figure()
plt.subplot(2,1,1)
plt.plot(time,mean_diff_phi,'r')
plt.xlabel("time(in ps)")
plt.ylabel("<Phi>-<Phi0> degrees")
plt.title("Average change of angle from COM (3rd layer) vs time at 0.01K")

plt.subplot(2,1,2)
plt.plot(time,var_diff_phi,'r')
plt.xlabel("time(in ps)")
plt.ylabel("Var(Phi)-Var(Phi0)")
plt.title("Variance of change of angle from COM (3rd layer) vs time at 0.01K")

plt.subplots_adjust(hspace=0.5,wspace=0.5)

plt.savefig('Rotation_COM_0.01K_TS.png')

plt.show()

vec_lay4_cosphi=[]
vec_lay4_sinphi=[]
vec_lay4_dist=[]
for i in range(len(vec_lay4)):
	vec_lay4_dist.append(m.sqrt(vec_lay4[i][0]**2+vec_lay4[i][1]**2))

#print(vec_lay4_dist)

for i in range(len(vec_lay4)):
	vec_lay4_cosphi.append(vec_lay4[i][0]/vec_lay4_dist[i])
	vec_lay4_sinphi.append(vec_lay4[i][1]/vec_lay4_dist[i])

#print(vec_lay4_cosphi)
#print(vec_lay4_sinphi)


# Find angles for each vectors for each snapshot

vec_lay4_phi=[]
for i in range(len(vec_lay4)):
	if(vec_lay4_cosphi[i]>=0):
		if(vec_lay4_sinphi[i]>=0):
#Ist quadrant
			vec_lay4_phi.append(m.acos(vec_lay4[i][0]/vec_lay4_dist[i])*180./m.pi) 
		elif(vec_lay4_sinphi[i]<0):
#IVth quadrant
			vec_lay4_phi.append(-(m.acos(vec_lay4[i][0]/vec_lay4_dist[i])*180./m.pi)) 
	elif(vec_lay4_cosphi[i]<0):
		if(vec_lay4_sinphi[i]>=0):
#IInd quadrant 
			vec_lay4_phi.append(m.acos(vec_lay4[i][0]/vec_lay4_dist[i])*180./m.pi) 
		elif(vec_lay4_sinphi[i]<0):
#IIIrd quadrant
			vec_lay4_phi.append(-(m.acos(vec_lay4[i][0]/vec_lay4_dist[i])*180./m.pi))

#print(vec_lay4_phi,len(vec_lay4_phi))


mean_lay4_phi=[]
var_lay4_phi=[]
for s in range(steps):
	mean_lay4_phi.append(stat.mean(vec_lay4_phi[s*(maxatom):(s+1)*(maxatom)]))	
	var_lay4_phi.append(stat.variance(vec_lay4_phi[s*(maxatom):(s+1)*(maxatom)]))	
	
#print(mean_lay4_phi,len(mean_lay4_phi))
#print(var_lay4_phi,len(var_lay4_phi))

mean_diff_phi_1=[]
var_diff_phi_1=[]
time=[]
for i in range(len(mean_lay4_phi)-1):
	mean_diff_phi_1.append(mean_lay4_phi[i+1]-mean_lay4_phi[0])
	var_diff_phi_1.append(var_lay4_phi[i+1]-var_lay4_phi[0])
	time.append(i+1)

print(mean_diff_phi_1,len(mean_diff_phi_1))
print(var_diff_phi_1,len(var_diff_phi_1))
print(time)

fig=plt.figure()
plt.subplot(2,1,1)
plt.plot(time,mean_diff_phi_1,'r')
plt.xlabel("time(in ps)")
plt.ylabel("<Phi>-<Phi0> degrees")
plt.title("Average change of angle from initial COM (3rd Layer) vs time at 0.01K")

plt.subplot(2,1,2)
plt.plot(time,var_diff_phi_1,'r')
plt.xlabel("time(in ps)")
plt.ylabel("Var(Phi)-Var(Phi0)")
plt.title("Variance of change of angle from initial COM (3rd Layer) vs time at 0.01K")

plt.subplots_adjust(hspace=0.5,wspace=0.5)

plt.savefig('Rotation_COM0_0.01K_TS.png')

plt.show()

