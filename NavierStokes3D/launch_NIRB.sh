#! /bin/bash

#Elise Grosjean
#lance NIRB + POD Greedy Parabolique avec mu=1, a comparer avec erreur FEM fine, ex: nev=10 (H=h*2 et h=H²) finetime=0.01 0.02 0.05 0.1
# uses :

#           - meshio-convert (vtk) to convert vtu files to vtk...

finetime="0.01 0.02 0.05 0.1"
finemesh="140 70 30 15"

#grossier
coarsetimediv="0.02 0.04 0.1 0.2" #hdiv2
coarsetimesqrt="0.1 0.1414 0.22 0.32" #sqrth

coarsemeshdiv="70 36 15 7"
coarsemeshsqrt="15 10 7 5"
# 1-9 5-1 5-2 5-5
for param in 1 #2 5 8 9
do
	
    for nev in 2 #number of modes
    do
    
	for values in 4 #in 1 2 3 4
	do
	    ## for div2 ##
	    
	    ## for sqrth ##
	    tau=$(echo $finetime | cut -d ' ' -f $values)
	    echo nev: $nev
	    taucoarse=$(echo $coarsetimesqrt | cut -d ' ' -f $values)
	    
	    python3.11 NirbGreedy.py $nev $tau sqrth $taucoarse $param1
	    
	    echo !!!!!!!!!! sqrt ... time $tau !!!!!!!!!!!!!
	    
	
	done
    done
done
