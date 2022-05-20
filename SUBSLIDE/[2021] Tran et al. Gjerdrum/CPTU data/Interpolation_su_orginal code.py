# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 10:03:12 2022

@author: ivandep
"""

import numpy as np
from scipy.interpolate import griddata, RBFInterpolator, interp1d
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import time
import pandas as pd

start=time.time()

# Folder with files
fileLoc="C:\\Users\\agnet\Documents\\NTNU\Master\Gjerdrum boreholes\\"


# Open the file with borFile
borFile=open(fileLoc+'Boreholes_moved_origin.txt','r') # File with coordinates of boreholes

# Counter
numBor=-1

# Borehole data
boreholeData=[]

for line in borFile:#Read all data from borFile to boreholeData
    if(numBor==-1):
        # Skip the first line
        # Update the counter
        numBor+=1

    else:
        # Split the string
        boreholeData.append(line.split())
        # Update the counter
        numBor+=1

# Close the file
borFile.close()

# Define a list with observations
observations=[]

# Interpolation of borehole data vertically
nVer=200 #200 points for su profile for all boreholes

# Read the remaining files to build a database
for i in range(numBor):
    # Open the borehole file
    dataIn=open(fileLoc+boreholeData[i][0]+'.txt','r')
    # Skip the first line
    dataIn.readline()
    
    # Vertical values of su (list)
    su=[]
    
    for data in dataIn:
        
        su.append([data.split()[0],data.split()[1]])
        
    dataIn.close()
    
    # Interpolate su values vertically
    su=np.array(su,dtype=np.float64) #Array
    
    suInterp=interp1d(su[:,0],su[:,1])#makes a function su=f(z)
    
    zVal=np.linspace(np.min(su[:,0]),np.max(su[:,0]),num=nVer)#makes 200 points between top and bottom
    
    suVal=suInterp(zVal)#makes 200 su-values based on function above
        
    for j in range(nVer):
        observations.append([boreholeData[i][1],boreholeData[i][2],zVal[j],suVal[j]])# x y z su for all 200 points for each borehole

# Convert list to an array
observations=np.array(observations,dtype=np.float64)

# Define the bounding box
lBound=np.min(observations[:,0:3],0)#lower x y z

uBound=np.max(observations[:,0:3],0)#upper x y z

#Bolken under her må endres når jeg får en liste over materialpunkter fra QA
#--------------------------------------
# Discretize the domain
# nx=50; ny=50; nz=50 #Number of points in each direction 

# # Define 1D discretization
# x=np.linspace(lBound[0]*0.9,uBound[0]*1.1,num=nx)
# y=np.linspace(lBound[1]*0.9,uBound[1]*1.1,num=ny)
# z=np.linspace(lBound[2]*0.5,uBound[2],num=nz)

# # Generate mesh of points
# xx,yy,zz=np.meshgrid(x,y,z)

#----------------------------
result_file = open(r'C:\Users\agnet\Documents\NTNU\Master\4.1 - Grunnundersøkelser\Su\layer3_su.txt','w')

# def writeResults(points, su):
#     points = np.insert(points,3, su, axis=1)
#     np.savetxt(result_file, points, delimiter ='\t')
    


filename='C:\\Users\\agnet\\Documents\\NTNU\Master\\4.1 - Grunnundersøkelser\Su\layer3.pts'
points=[]
maxrow=10**6 #Read 1 million points each iteration 
fmt = '%1.2f','%1.2f','%1.2f', '%.0f' #Define format
for i in range(99):
    points=np.loadtxt(filename, skiprows=i*maxrow, max_rows=maxrow)
    predictionRBF=RBFInterpolator(observations[:,0:3],observations[:,3],kernel='linear',epsilon=10,degree=1)(points)
    points = np.insert(points,3, predictionRBF, axis=1)
    np.savetxt(result_file, points, fmt=fmt, delimiter ='\t')
    #writeResults(points, predictionRBF)
    

result_file.close()


print('Elapsed time:%f s' % (time.time()-start))

# idz=(points[:,1]>600)*(points[:,1]<620)==1

# idx=(points[:,0]>180)*(points[:,0]<200)==1

# Plot predictions
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# p = ax.scatter(points[:,0], points[:,1], points[:,2], c=predictionRBF[:], cmap='viridis', linewidth=0.5, vmin=20,vmax=220)
# plt.xlabel('X')
# plt.ylabel('Y')
# ax.set_zlabel('Z')
# fig.colorbar(p)

# #Plot observations plus some prediction
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# p = ax.scatter(points[idx,0], points[idx,1], points[idx,2], c=predictionRBF[idx], cmap='viridis', linewidth=0.5, vmin=20,vmax=220)
# p = ax.scatter(observations[:,0], observations[:,1], observations[:,2], c=observations[:,3], cmap='viridis', linewidth=0.5, vmin=20,vmax=220)
# plt.xlabel('X')
# plt.ylabel('Y')
# ax.set_zlabel('Z')
# fig.colorbar(p)


# #Plot observations plus some prediction
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# p = ax.scatter(points[idz,0], points[idz,1], points[idz,2], c=predictionRBF[idz], cmap='viridis', linewidth=0.5, vmin=20,vmax=220)
# p = ax.scatter(observations[:,0], observations[:,1], observations[:,2], c=observations[:,3], cmap='viridis', linewidth=0.5, vmin=20,vmax=220)
# plt.xlabel('X')
# plt.ylabel('Y')
# ax.set_zlabel('Z')
# fig.colorbar(p)

# # Compare predictions with observations
# plt.scatter(observations[:,3],RBFInterpolator(observations[:,0:3],observations[:,3],kernel='linear',epsilon=0.1,degree=1)(observations[:,0:3]))
