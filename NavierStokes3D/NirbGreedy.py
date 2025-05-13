
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
onlineParameter2="1"#str(sys.argv[6])
print( "online Parameter", onlineParameter1,' ', onlineParameter2)
## Directories
currentFolder=os.getcwd()
dataFolder=currentFolder
FinedataFolderU=osp.join(dataFolder,'FineSnapshots')#+sys.argv[2]) #for fine snapshots
#CoarsedataFolderU=osp.join(dataFolder,'CoarseSnapshots')#+sys.argv[3]+'/'+sys.argv[4]) #for coarse snapshots
CoarsedataFolderU=osp.join(dataFolder,'CoarseSnapshotsTrunc3')#+sys.argv[3]+'/'+sys.argv[4]) #for coarse snapshots
FinedataFolder=FinedataFolderU #osp.join(dataFolder,'FineSnapshotsPhi/'+sys.argv[2]) #for fine snapshots
CoarsedataFolder=CoarsedataFolderU#osp.join(dataFolder,'CoarseSnapshotsPhi/'+sys.argv[3]+'/'+sys.argv[4]) #for coarse snapshots
print("fine folder: ", FinedataFolder)
print("coarse folder: ", CoarsedataFolder)


ns=0 #number of snapshots
count1=0
for _, folder, _ in os.walk(FinedataFolder): #number of folders in FineData
    count1 += len(folder)

ns=6 # 19-1#count1-1 #-1 because of the online parameter not included in offline snapshots
print("Number of snapshots: ",ns)
#ns=18


nev=int(sys.argv[1])   #nombre de modes

time=0.0 #init
dimension=3 #2D
           
TF=1#len([name for name in os.listdir(FinedataFolder+"/1-1/")])
print("Number of fine time steps: ",TF)

TG=1#len([name for name in os.listdir(CoarsedataFolder+"/1-1/")])
print("Number of coarse time steps: ",TG)

    

"""
-------------------
###  Read fine mesh
------------------- 
"""

meshFileName = FinedataFolder + "/Snapshot_1.vtu";
mesh=MR.Readmesh(meshFileName)
#mesh.nodes= mesh.nodes[:,:2] #2D

print("Fine mesh defined in " + meshFileName + " has been read")
nbeOfComponentsPrimal = 3 # 1 field 
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

parameters=[]
#parameters=['0.5','1','1.5','2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7','7.5','8','8.5','9','9.5']
parameters=['0','1','2','3','4','5']


#for i in range(1,ns+2):
#    if(float(0.5*i)%1>1e-3):
#        parameters.append(str(0.5*i))
#    else:
#        parameters.append(str(int(0.5*i)))

#print("param:",parameters)
#parameters.remove(onlineParameter) #online parameter mu=1
print("parameters :",parameters)


snapshots=[]

for e,i in enumerate(parameters):
    #snapshotsTime=[]
    #for time in range(0,TF):   
    snapshot =MR.VTKReadToNp("Velocity",FinedataFolder+"/Snapshot_",i)

    snapshot=snapshot.flatten()
    #    snapshotsTime.append(snapshot)
    snapshots.append(snapshot)


snapshotsH=[]

for e,i in enumerate(parameters):
    #snapshotsHTime=[]
    #for time in range(0,TG):


    snapshotH =MR.VTKReadToNp("Velocity",CoarsedataFolder+"/Snapshot_",i).flatten()
    snapshotsH.append(snapshotH)
    

############################################################
"""          Greedy                                      """
############################################################

#print("ComputeL2ScalarProducMatrix ...")

l2ScalarProducMatrix = FT.ComputeL2ScalarProducMatrix( mesh, nbeOfComponentsPrimal)
l2ScalarProducMatrixCoarse = FT.ComputeL2ScalarProducMatrix( mesh2, nbeOfComponentsPrimal)
#h1ScalarProducMatrix = FT.ComputeH10ScalarProductMatrix(mesh, nbeOfComponentsPrimal)

##### ALGO (full) GREEDY
reducedOrderBasisPhi=GD.Greedy(snapshots,l2ScalarProducMatrix,NumberOfModes=nev) #GD.greedy_algorithm(snapshots,1,l2ScalarProducMatrix,nev)
reducedOrderBasisPhiCoarse=GD.Greedy(snapshotsH,l2ScalarProducMatrixCoarse,NumberOfModes=nev) #GD.greedy_algorithm(snapshots,1,l2ScalarProducMatrix,nev)#greedy_algorithm(snapshotsH,1,l2ScalarProducMatrixCoarse,nev)

nev2=nev
nev3=nev
#assert(nev2==nev3)

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
"""          Save data for online part                   """
############################################################
### save reduced basis

#outputNamePhi = "reducedOrderBasisPhi.pkl"
#outputPhi = open(outputNamePhi, "wb")
#pickle.dump(reducedOrderBasisPhi, outputPhi)
#outputPhi.close()

#OperatorOutput = "Operator.pkl"
#outputOp=open(OperatorOutput,"wb")
#pickle.dump(Operator, outputOp)
#outputOp.close()

#REGOutput = "Rectification.pkl"
#outputReg=open(REGOutput,"wb")
#pickle.dump(RI, outputReg)
#outputReg.close()
############################################################
"""          Online part                                 """
############################################################

jparam=onlineParameter2
for iparam in [onlineParameter1]:

    snapshotsHTime=[]
    #for time in range(TG):

    snapshotH =MR.VTKReadToNp("Velocity",currentFolder+"/snapshotH3ex_","H0").flatten()
   
    for time in range(1):    
        R=RI[time,:,:] #nev1 nev2
       
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
        #reconstructedCompressedSolution =MR.VTKReadToNp("Velocity",currentFolder+"/snapshot_exh_","0")
        #reconstructedCompressedSolution=MR.VTKReadToNp("Velocity",currentFolder+"/snapshotH3ex_","H0")
        ##################################################
        #######   saving solution in VTK ############

        savedata=reconstructedCompressedSolution#[:,0]#
        savedata=savedata.reshape((numberOfNodes,3))# nbeOfComponentsPrimal))
        
        VTKBase = MR.VTKReadmesh(meshFileName)
        #SVTKW.numpyToVTKWrite(VTKBase,reconstructedCompressedSolution.reshape((numberOfNodes,3)),"NIRB_approximation_"+str(time)+"_"+str(nev)+".vtu")
        SVTKW.numpyToVTKWrite(VTKBase,savedata,"NIRB_approximation_"+str(time)+"_"+str(nev)+".vtu")

        #SVTKW.numpyToVTKWrite(VTKBase,savedata,SolutionName="NIRB_approximation_"+str(time)+"_"+str(nev)+".vtu",FieldName="U")
        ##################################################

        
