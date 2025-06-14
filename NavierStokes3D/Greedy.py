# -*- coding: utf-8 -*-
## NIRB script python Greedy Algorithm

## Elise Grosjean
## 01/2025


from BasicTools.FE import FETools as FT
import numpy as np
from scipy import linalg

def orthogonality_check(Matrix,CorrelationMatrix):
    """
    This fucntion check for the pairwise orthogonality of the new basis
    """
    list_ = list(Matrix)
    dot_matrix = np.array([[np.dot(CorrelationMatrix.dot(item1), item2) for item1 in list_] for item2 in list_])
    if (dot_matrix - np.eye(dot_matrix.shape[0]) < 1e-10).all():
        return True
    else:
        error = dot_matrix - np.eye(dot_matrix.shape[0])
        print("max error with identity: ",np.max(error), "min error : ",np.min(error))
        return False

def CheckLinearIndependence(Vectors): #check obvious linear dependence on the vectors
    import sympy
    _, inds = sympy.Matrix(Vectors).T.rref()
    
    if len(inds)<np.shape(Vectors)[0]:
        return inds
    else:
        return True

def Norm(Matrix,v):

    #print(v.shape)
    #print(Matrix.shape)
    
    return np.sqrt(np.dot(Matrix.dot(v),v))



##### ALGO GREEDY ######
# Number Of Modes can be a priori given or retrieve thanks to the tolerance threshold if no h1ScalarProducMatrix yields to L2 orthonormal basis

def Greedy(snapshots,snapshotCorrelationOperator,h1ScalarProducMatrix=None,NumberOfModes=0,Tol=1e-6):
    """
    Greedy algorithm for the construction of the reduced basis
    orthogonal basis in H1 et L2 et orthonormalized in L2
    #Algo as in https://hal.archives-ouvertes.fr/hal-01897395
    """
    if snapshotCorrelationOperator is None:
        snapshotCorrelationOperator = sparse.eye(collectionProblemData.GetSolutionsNumberOfDofs(solutionName))

    SnapshotsNorm=[]
    for s in snapshots:
        SnapshotsNorm.append(Norm(snapshotCorrelationOperator,s))
    snapshots = np.array(snapshots)
    
    DegreesOfFreedom=np.shape(snapshotCorrelationOperator)[0]
    NumberOfSnapshots=np.shape(snapshots)[0]

    if CheckLinearIndependence(snapshots) != True: ##If dependence, remove the vectors
        print("snapshots linearly dependent, removing the corresponding vectors")
        #Inds=list(CheckLinearIndependence(snapshots)) #check obvious linear dependence on the vectors
        #SnapshotsNorm=SnapshotsNorm[Inds]
        #NumberOfSnapshots=len(Inds)
        
    if NumberOfModes==0:
        
        reducedOrderBasisU=np.zeros((NumberOfSnapshots,DegreesOfFreedom)) #ns, nbd
    else:
        reducedOrderBasisU=np.zeros((NumberOfModes,DegreesOfFreedom)) #nev, nbd
    Index=SnapshotsNorm.index(max(SnapshotsNorm)) #first mode
    print("Mode 0: ", Index)

    reducedOrderBasisU[0,:]=snapshots[Index]/SnapshotsNorm[Index] #first mode
    ListeIndex=[Index] #first snapshot 
    BasisNorm=[SnapshotsNorm[Index]]
    
    Basis=[snapshots[Index]]
    MatrixBasisProduct=[snapshotCorrelationOperator.dot(Basis[0])]
    if NumberOfModes>0:
        for n in range(1,NumberOfModes):
            print("Mode ",n)
        
            TestVector=dict() # dictionnary: vector in the reduced basis if maxTest if maximum
            for j in range(NumberOfSnapshots):
                #print(j, " list: ", ListeIndex)
                if not (j in ListeIndex):
                    w=snapshots[j]-np.sum((b*np.dot(MatrixBasisProduct[k],snapshots[j])/BasisNorm[k]**2 for k,b in enumerate(Basis)),axis=0)#potential vector to add in the reduced basis
                    if (w > 1e-8).any() and SnapshotsNorm[j]>1e-10: #if index not yet in the basis:
                        NormW=Norm(snapshotCorrelationOperator,w)#np.sqrt(np.dot((l2ScalarProducMatrix.dot(w)),w))
                        GreedyMaximumTest=NormW/SnapshotsNorm[j] #we seek the max
                        TestVector[j]=[GreedyMaximumTest,w,NormW]
                        #
            Index=max(TestVector, key = lambda k: TestVector[k][0]) #index of the snapshot used
            #print(TestVector[Index])
            print("index",Index)
            ListeIndex.append(Index) #adding in the list
        
            Basis.append(TestVector[Index][1])
            BasisNorm.append(TestVector[Index][2])
            MatrixBasisProduct.append(snapshotCorrelationOperator.dot(Basis[n]))
                                
            reducedOrderBasisU[n,:]=(TestVector[Index][1]/TestVector[Index][2]) #orthonormalization in L2
            
    Check=orthogonality_check(reducedOrderBasisU,snapshotCorrelationOperator)
    print("orthogonality ", Check) #if no orthogonality, it may be due to linear dependence of vectors/ non-stability of GS
    if Check==False: #redo Gram-Schmidt procedure 
        NewReducedOrderBasisU=np.zeros((NumberOfModes,DegreesOfFreedom))
        BasisNorm=[Norm(snapshotCorrelationOperator,reducedOrderBasisU[0])]
        NewReducedOrderBasisU[0]=reducedOrderBasisU[0]/BasisNorm[0]
        
        for i in range(1,NumberOfModes):
            NewReducedOrderBasisU[i]=reducedOrderBasisU[i]-np.sum((NewReducedOrderBasisU[k]*np.dot(snapshotCorrelationOperator.dot(NewReducedOrderBasisU[k]),reducedOrderBasisU[i])/(BasisNorm[k]**2) for k in range(i)),axis=0)#potential vector to add in the reduced basis
            BasisNorm.append(Norm(snapshotCorrelationOperator,NewReducedOrderBasisU[i]))
        for i in range(NumberOfModes):
            reducedOrderBasisU[i]=NewReducedOrderBasisU[i]/BasisNorm[i]#Norm(snapshotCorrelationOperator,NewReducedOrderBasisU[i])
    Check=orthogonality_check(reducedOrderBasisU,snapshotCorrelationOperator)
    print("orthogonality ", Check) #if no orthogonality, it may be due to linear dependence of vectors/ non-stability of GS

    """
    for i in range(NumberOfModes):
        for j in range(i+1):
            t=snapshotCorrelationOperator.dot(reducedOrderBasisU[i,:])
            norm=t.dot(reducedOrderBasisU[j,:])
            print(i,j," ",norm)
    """
    ### H1 Orthogonalization
      
    if h1ScalarProducMatrix!=None:
        normRed=[]
        K=np.zeros((NumberOfModes,NumberOfModes)) #rigidity matrix
        M=np.zeros((NumberOfModes,NumberOfModes)) #mass matrix
        for i in range(NumberOfModes):
            matVecH1=h1ScalarProducMatrix.dot(reducedOrderBasisU[i,:])
            matVecL2=snapshotCorrelationOperator.dot(reducedOrderBasisU[i,:])
            for j in range(NumberOfModes):
                if i>=j:
                   
                    K[i,j]=np.dot(matVecH1,reducedOrderBasisU[j,:])
                    M[i,j]=np.dot(matVecL2,reducedOrderBasisU[j,:])
                    K[j,i]=K[i,j]
                    M[j,i]=M[i,j]
    
    
        # on resoud Kv=lambd Mv
        #mpiReducedCorrelationMatrixM = np.zeros((nev, nev))
        #MPI.COMM_WORLD.Allreduce([M,  MPI.DOUBLE], [mpiReducedCorrelationMatrixM,  MPI.DOUBLE])
        eigenValues,vr=linalg.eig(K, b=M) #eigenvalues + right eigenvectors
        idx = eigenValues.argsort()[::-1]
        eigenValues = eigenValues[idx]
        eigenVectors = vr[:, idx]
        reducedOrderBasisU=np.dot(eigenVectors.transpose(),reducedOrderBasisU)

        for i in range(NumberOfModes):
            reducedOrderBasisNorm=np.sqrt(reducedOrderBasisU[i,:]@(snapshotCorrelationOperator@reducedOrderBasisU[i,:]))
            reducedOrderBasisU[i,:]/=reducedOrderBasisNorm#np.sqrt(M[i,i]) #L2 orthonormalization
    
    
    return reducedOrderBasisU,ListeIndex










def GreedyNew(snapshots,snapshotCorrelationOperator,ListeIndex,h1ScalarProducMatrix=None,NumberOfModes=0,Tol=1e-30): ##create a RB from a priori given parameters

    
    """
    Greedy algorithm for the construction of the reduced basis
    orthogonal basis in H1 et L2 et orthonormalized in L2
    #Algo as in https://hal.archives-ouvertes.fr/hal-01897395
    """
    if snapshotCorrelationOperator is None:
        snapshotCorrelationOperator = sparse.eye(collectionProblemData.GetSolutionsNumberOfDofs(solutionName))

    SnapshotsNorm=[]
    for s in snapshots:
        SnapshotsNorm.append(Norm(snapshotCorrelationOperator,s))
    snapshots = np.array(snapshots)
    
    DegreesOfFreedom=np.shape(snapshotCorrelationOperator)[0]
    NumberOfSnapshots=np.shape(snapshots)[0]

    if CheckLinearIndependence(snapshots) != True: ##If dependence, remove the vectors
        print("snapshots linearly dependent, removing the corresponding vectors")
        
    if NumberOfModes==0:        
        reducedOrderBasisU=np.zeros((NumberOfSnapshots,DegreesOfFreedom)) #ns, nbd
    else:
        reducedOrderBasisU=np.zeros((NumberOfModes,DegreesOfFreedom)) #nev, nbd
    Index=ListeIndex[0] #SnapshotsNorm.index(max(SnapshotsNorm)) #first mode
    #print("Mode 0: ", Index)

    reducedOrderBasisU[0,:]=snapshots[Index]/SnapshotsNorm[Index] #first mode
    
    BasisNorm=[SnapshotsNorm[Index]]
    
    Basis=[snapshots[Index]]
    MatrixBasisProduct=[snapshotCorrelationOperator.dot(Basis[0])]
    if NumberOfModes>0:
        for n in range(1,NumberOfModes):
            #print("Mode: ",n)     
            TestVector=dict() # dictionnary: vector in the reduced basis if maxTest if maximum
            Index=ListeIndex[n] #max(TestVector, key = lambda k: TestVector[k][0]) #index of the snapshot used
            j=Index
            w=snapshots[j]-np.sum((b*np.dot(MatrixBasisProduct[k],snapshots[j])/BasisNorm[k]**2 for k,b in enumerate(Basis)),axis=0)#potential vector to add in the reduced basis
            if (w > 1e-8).any() and SnapshotsNorm[j]>1e-10: #if index not yet in the basis:

                NormW=Norm(snapshotCorrelationOperator,w)#np.sqrt(np.dot((l2ScalarProducMatrix.dot(w)),w))
                GreedyMaximumTest=NormW/SnapshotsNorm[j] #we seek the max
                TestVector[j]=[GreedyMaximumTest,w,NormW]
                
            #print(TestVector[Index])
            #print("index",Index)
            
            Basis.append(TestVector[Index][1])
            BasisNorm.append(TestVector[Index][2])
            MatrixBasisProduct.append(snapshotCorrelationOperator.dot(Basis[n]))
                                
            reducedOrderBasisU[n,:]=(TestVector[Index][1]/TestVector[Index][2]) #orthonormalization in L2
            
    Check=orthogonality_check(reducedOrderBasisU,snapshotCorrelationOperator)
    #print("orthogonality ", Check) #if no orthogonality, it may be due to linear dependence of vectors/ non-stability of GS
    if Check==False: #redo Gram-Schmidt procedure 
        NewReducedOrderBasisU=np.zeros((NumberOfModes,DegreesOfFreedom))
        BasisNorm=[Norm(snapshotCorrelationOperator,reducedOrderBasisU[0])]
        NewReducedOrderBasisU[0]=reducedOrderBasisU[0]/BasisNorm[0]
        
        for i in range(1,NumberOfModes):
            NewReducedOrderBasisU[i]=reducedOrderBasisU[i]-np.sum((NewReducedOrderBasisU[k]*np.dot(snapshotCorrelationOperator.dot(NewReducedOrderBasisU[k]),reducedOrderBasisU[i])/(BasisNorm[k]**2) for k in range(i)),axis=0)#potential vector to add in the reduced basis
            BasisNorm.append(Norm(snapshotCorrelationOperator,NewReducedOrderBasisU[i]))
        for i in range(NumberOfModes):
            reducedOrderBasisU[i]=NewReducedOrderBasisU[i]/BasisNorm[i]#Norm(snapshotCorrelationOperator,NewReducedOrderBasisU[i])
    Check=orthogonality_check(reducedOrderBasisU,snapshotCorrelationOperator)
    print("orthogonality ", Check) #if no orthogonality, it may be due to linear dependence of vectors/ non-stability of GS
    ### H1 Orthogonalization
      
    if h1ScalarProducMatrix!=None:
        normRed=[]
        K=np.zeros((NumberOfModes,NumberOfModes)) #rigidity matrix
        M=np.zeros((NumberOfModes,NumberOfModes)) #mass matrix
        for i in range(NumberOfModes):
            matVecH1=h1ScalarProducMatrix.dot(reducedOrderBasisU[i,:])
            matVecL2=snapshotCorrelationOperator.dot(reducedOrderBasisU[i,:])
            for j in range(NumberOfModes):
                if i>=j:
                   
                    K[i,j]=np.dot(matVecH1,reducedOrderBasisU[j,:])
                    M[i,j]=np.dot(matVecL2,reducedOrderBasisU[j,:])
                    K[j,i]=K[i,j]
                    M[j,i]=M[i,j]
    
    
        # Kv=lambd Mv
        # mpiReducedCorrelationMatrixM = np.zeros((nev, nev))
        # MPI.COMM_WORLD.Allreduce([M,  MPI.DOUBLE], [mpiReducedCorrelationMatrixM,  MPI.DOUBLE])
        eigenValues,vr=linalg.eig(K, b=M) #eigenvalues + right eigenvectors
        idx = eigenValues.argsort()[::-1]
        eigenValues = eigenValues[idx]
        eigenVectors = vr[:, idx]
        reducedOrderBasisU=np.dot(eigenVectors.transpose(),reducedOrderBasisU)

        for i in range(NumberOfModes):
            reducedOrderBasisNorm=np.sqrt(reducedOrderBasisU[i,:]@(snapshotCorrelationOperator@reducedOrderBasisU[i,:]))
            reducedOrderBasisU[i,:]/=reducedOrderBasisNorm#np.sqrt(M[i,i]) #L2 orthonormalization
    
    return reducedOrderBasisU



