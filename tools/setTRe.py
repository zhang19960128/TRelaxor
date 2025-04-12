#!/usr/bin/env python
import random
import sys
import math
import numpy as np
def changeindex(index,cell):
	'''change the index to 3D array[index_x,index_y,index_z]'''
	re=np.zeros(3);
	re[2]=math.floor(index/cell/cell);
	re[1]=math.floor((index-re[2]*cell*cell)/cell);
	re[0]=index-re[2]*cell*cell-re[1]*cell;
	return re;
def backindex(indexlist,cell):
	'''change the three dimensional index into a one dimentional one'''
	return indexlist[0]%cell+(indexlist[1]%cell)*cell+(indexlist[2]%cell)*cell*cell;
def typein(index,Asite):
	a=[];
	for i in index:
		a.append(Asite[int(i)]);
	return a;
## Temperature Stable Relaxor x*Ba0.75Ca0.25TiO3-(1-x)*Bi(Mg0.5Ti0.5)O3
cell=2;'''here are the period on x y z cell length'''
ba_valen=1.34730;
ca_valen=1.34730;
bi_valen=1.482504;
ti_valen=1.289050;
mg_valen=1.018642;
ox_valen=-0.878783;
cocentx=float(sys.argv[1]);
bacalist=random.sample(range(cell*cell*cell),int(cell*cell*cell*cocentx));
balist=list(random.sample(bacalist,int(0.75*len(bacalist))));
calist=list(set(bacalist)-set(balist));
bilist=list(set(range(cell*cell*cell))-set(bacalist));
mglist=list(random.sample(range(cell*cell*cell),int(len(range(cell*cell*cell))*(1-cocentx)*0.5)));
tilist=list(set(range(cell*cell*cell))-set(mglist));
ca=open("cadata.txt","w");
for i in calist:
	ca.write(str(i)+"\n");
ca.close();
data=open("mixdata.BTO","w");
data.write("#LAMMPS 100 water\n");
data.write("\n");
data.write(str(5*cell*cell*cell)+" atoms\n");
data.write("6 atom types\n");
data.write("\n");
cell_length=4.04;
data.write("0.0 "+str(cell_length*cell)+" "+"xlo xhi\n");
data.write("0.0 "+str(cell_length*cell)+" "+"ylo yhi\n");
data.write("0.0 "+str(cell_length*cell)+" "+"zlo zhi\n");
data.write("\n");
data.write("\n");
data.write("Masses\n");
data.write("\n");
data.write("1 137.327  #Ba\n");
data.write("2 40.078  #Ca\n");
data.write("3 208.98 #Bi\n");
data.write("4 91.224 #Ti\n");
data.write("5 24.305 #Mg\n");
data.write("6 15.9993  #o\n");
data.write("\n");
data.write("Atoms\n");
data.write("\n");
Asite=np.zeros(cell*cell*cell);
for i in range(cell*cell*cell):
	if i in balist:
		index=changeindex(i,cell);
		data.write(str(i+1)+" "+"1"+" "+"1"+" "+str(ba_valen)+" "+str(cell_length*index[0])+" "+str(index[1]*cell_length)+" "+str(index[2]*cell_length)+"\n");
	elif i in calist:
		index=changeindex(i,cell);
		data.write(str(i+1)+" "+"1"+" "+"2"+" "+str(ca_valen)+" "+str(cell_length*index[0])+" "+str(index[1]*cell_length)+" "+str(index[2]*cell_length)+"\n");
	elif i in bilist:
		index=changeindex(i,cell);
		data.write(str(i+1)+" "+"1"+" "+"3"+" "+str(bi_valen)+" "+str(cell_length*index[0])+" "+str(index[1]*cell_length)+" "+str(index[2]*cell_length)+"\n");
for i in range(cell*cell*cell):
  if i in tilist:
    index=changeindex(i,cell);
    data.write(str(i+cell*cell*cell*1+1)+" "+"1"+" "+"4"+" "+str(ti_valen)+" "+str(cell_length*index[0]+cell_length/2)+"  "+str(cell_length*index[1]+cell_length/2)+" "+str(cell_length*index[2]+cell_length/2)+"\n");
  elif i in mglist:
    index=changeindex(i,cell);
    data.write(str(i+cell*cell*cell*1+1)+" "+"1"+" "+"5"+" "+str(mg_valen)+" "+str(cell_length*index[0]+cell_length/2)+"  "+str(cell_length*index[1]+cell_length/2)+" "+str(cell_length*index[2]+cell_length/2)+"\n");
for i in range(cell*cell*cell):
	'''bottom to up'''
	index=changeindex(i,cell);
	data.write(str(i+2*cell*cell*cell+1)+" "+"1"+" "+"6"+" "+str(ox_valen)+" "+str(cell_length*index[0]+cell_length/2)+" "+str(cell_length*index[1]+cell_length/2)+" "+str(cell_length*index[2])+"\n");
for i in range(cell*cell*cell):
	'''front to back'''
	index=changeindex(i,cell);
	data.write(str(i+3*cell*cell*cell+1)+" "+"1"+" "+"6"+" "+str(ox_valen)+" "+str(cell_length*index[0]+cell_length/2)+" "+str(cell_length*index[1])+" "+str(cell_length*index[2]+cell_length/2)+"\n");
for i in range(cell*cell*cell):
	'''left to right'''
	index=changeindex(i,cell);
	data.write(str(i+4*cell*cell*cell+1)+" "+"1"+" "+"6"+" "+str(ox_valen)+" "+str(cell_length*index[0])+" "+str(cell_length*index[1]+cell_length/2)+" "+str(cell_length*index[2]+cell_length/2)+"\n");
