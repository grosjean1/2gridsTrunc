
# -*- coding: utf-8 -*-
## NIRB parabolic test with OFFLINE/ONLINE DECOMPOSITION

## Elise Grosjean
## 01/2022



import numpy as np
import sys
import pickle
import os
import os.path as osp

import Readers as MR
import SolutionVTKWriter as SVTKW

from BasicTools.FE import FETools as FT
#import pickle

import Greedy as GD
#import SVD

from scipy import linalg
from scipy import interpolate
#from scipy.interpolate import griddata 
from scipy.interpolate import interp1d #time interpolation

from scipy.spatial import cKDTree #space interpolation
from scipy.sparse import coo_matrix
        
############################################################
"""          Initialization                              """
############################################################
onlineParameter1=str(sys.argv[5])
print( "online Parameter", onlineParameter1)
## Directories
currentFolder=os.getcwd()
dataFolder=currentFolder
FinedataFolder=osp.join(dataFolder,'FineSnapshots')#+sys.argv[2]) #for fine snapshots
CoarsedataFolder=osp.join(dataFolder,'CoarseSnapshotsTrunc3')
print("fine folder: ", FinedataFolder)
print("coarse folder: ", CoarsedataFolder)


ns=0 #number of snapshots
count1=0
for _, folder, _ in os.walk(FinedataFolder): #number of folders in FineData
    count1 += len(folder)

ns=6 
print("Number of snapshots: ",ns)

nev=int(sys.argv[1])   #number of modes
time=0.0 #init
dimension=3 #3D
           
TF=1
print("Number of fine time steps: ",TF)
TG=1
print("Number of coarse time steps: ",TG)


"""
-------------------
###  Read fine mesh
------------------- 
"""

meshFileName = FinedataFolder + "/Snapshot_1.vtu";
mesh=MR.Readmesh(meshFileName)
#mesh.nodes= mesh.nodes[:,:2] # in case 2D
print("Fine mesh defined in " + meshFileName + " has been read")

nbeOfComponentsPrimal = 3 # 1 field, 3 components 
numberOfNodes = mesh.GetNumberOfNodes()
print("DoF fine mesh ", numberOfNodes)
"""
-------------------
###  Read coarse mesh
------------------- 
"""
meshFileName2 = CoarsedataFolder + "/Snapshot_1.vtu";#//FineMesh/mesh1.msh"
mesh2=MR.Readmesh(meshFileName2)
#mesh2.nodes = mesh2.nodes[:,:2] #CAS 2D

print("Coarse mesh defined in " + meshFileName2 + " has been read")
numberOfNodes2 = mesh2.GetNumberOfNodes()
print("DoF coarse mesh ", numberOfNodes2)


"""
-------------------
###  mesh space interpolation
------------------- 
"""
inputnodes=mesh2.nodes
outputnodes=mesh.nodes
kdt = cKDTree(inputnodes)
nbtp = outputnodes.shape[0]
_, ids = kdt.query(outputnodes)
cols=ids
row = np.arange(nbtp)
data = np.ones(nbtp)
Operator=coo_matrix((data, (row, cols)), shape=(nbtp , inputnodes.shape[0]))
        
"""
-------------------
###  read all snapshots ...
------------------- 
"""

parameters=['0','1','2','3','4','5']
print("parameters :",parameters)

snapshots=[]

for e,i in enumerate(parameters):
    snapshot =MR.VTKReadToNp("Velocity",FinedataFolder+"/Snapshot_",i)
    snapshot=snapshot.flatten()
    snapshots.append(snapshot)


snapshotsH=[]

for e,i in enumerate(parameters):
    snapshotH =MR.VTKReadToNp("Velocity",CoarsedataFolder+"/Snapshot_",i).flatten()
    snapshotsH.append(snapshotH)
    

############################################################
"""          Greedy                                      """
############################################################

#print("ComputeL2ScalarProducMatrix ...")

l2ScalarProducMatrix = FT.ComputeL2ScalarProducMatrix( mesh, nbeOfComponentsPrimal)
l2ScalarProducMatrixCoarse = FT.ComputeL2ScalarProducMatrix( mesh2, nbeOfComponentsPrimal)
#h1ScalarProducMatrix = FT.ComputeH10ScalarProductMatrix(mesh, nbeOfComponentsPrimal)

##### ALGO GREEDY
reducedOrderBasisPhiCoarse,GlobalIndices=GD.Greedy(snapshotsH,l2ScalarProducMatrixCoarse,NumberOfModes=nev) # on coarse grid...
reducedOrderBasisPhi=GD.GreedyNew(snapshots,l2ScalarProducMatrix,GlobalIndices,NumberOfModes=nev) #on fine grid with same parameters..

nev2=nev

print("Number of modes after greedyPhi",nev2)
############################################################
"""          Rectification  sur Phi                      """
############################################################

RI=np.zeros((TF,nev2,nev2)) #global rectification over the time steps
for time in range(TF):

    alpha=np.zeros((nev2,ns)) #fine coefficients
    beta=np.zeros((ns,nev2)) #coarse coefficients
    betainv=np.zeros((nev2,ns)) 
    R=np.zeros((nev2,nev2)) #rectification matrix
    
    for j,elt in enumerate(parameters):

        u1PT = snapshots[j]
        u1T = snapshotsH[j]
        
        for i in range(nev2):
            alpha[i,j]=u1PT@(l2ScalarProducMatrix@reducedOrderBasisPhi[i,:])
        for i in range(nev2):
            beta[j,i]=u1T@(l2ScalarProducMatrixCoarse@reducedOrderBasisPhiCoarse[i,:])
            betainv[i,j]=beta[j,i]
        
    lambd=1e-10 #Thikonov regularization parameter
    Rns=np.zeros((nev2,ns))
    
    for j in range(ns):
        Rns[:,j]=np.linalg.inv(beta.transpose()@beta+lambd*np.eye(nev2))@betainv[:,j]
        
    for i in range(nev2):
        R[i,:]=Rns@alpha[i,:]

    RI[time,:,:]=R


############################################################
"""          Online part                                 """
############################################################

for iparam in [onlineParameter1]:
    
    snapshotH =MR.VTKReadToNp("Velocity",currentFolder+"/snapshotH3ex_","H0").flatten() # read coarse solution 
   
    for time in range(1):    
        R=RI[time,:,:] #nev1
       
        u1PT=snapshotH
        coef=np.zeros(nev2)
        CompressedSolutionUj=np.zeros(nev2)
        for j in range(nev2):
            CompressedSolutionUj[j]=u1PT@(l2ScalarProducMatrixCoarse@reducedOrderBasisPhiCoarse[j,:])

        for i in range(nev2):
            coef[i]=0
            for j in range(nev2):
                coef[i]+=R[i,j]*CompressedSolutionUj[j]

        reconstructedCompressedSolution = np.dot(coef, reducedOrderBasisPhi) #rectified nirb
        ##################################################
        #######   saving solution in VTK ############

        savedata=reconstructedCompressedSolution
        savedata=savedata.reshape((numberOfNodes, nbeOfComponentsPrimal))
        VTKBase = MR.VTKReadmesh(meshFileName)
        SVTKW.numpyToVTKWrite(VTKBase,savedata,"NIRB_approximation_"+str(time)+"_"+str(nev)+".vtu")

        ##################################################

        
