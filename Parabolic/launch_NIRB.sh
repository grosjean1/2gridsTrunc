#! /bin/bash

#Elise Grosjean
#lance NIRB + POD Greedy Parabolique avec mu=1, a comparer avec erreur FEM fine, ex: nev=10 (H=h*2 et h=H²) finetime=0.01 0.02 0.05 0.1
# uses :
#           - meshio-convert (vtk) to convert vtu files to vtk...

# FINE
finetime="0.01 0.02 0.05 0.1"
finemesh="140 70 30 15"

# COARSE
coarsetimediv="0.02 0.04 0.1 0.2" #hdiv2
coarsetimesqrt="0.1 0.1414 0.22 0.32" #sqrth

coarsemeshdiv="70 36 15 7"
coarsemeshsqrt="15 10 7 5"

param1=1
for param2 in 8  #2 #9 8
do
    for nev in 5  
    do
	if [[ $param -ne 3 && $param2 -ne 4 ]]; then    
	    for values in 4 #1 2 3 4  
	    do
       
		## for sqrth ##
		tau=$(echo $finetime | cut -d ' ' -f $values)
		echo nev: $nev
		taucoarse=$(echo $coarsetimesqrt | cut -d ' ' -f $values)
		python3.11 NirbGreedy.py $nev $tau sqrth $taucoarse $param1 $param2

		FileList=$(ls NIRB_approximation*) 

		#nev=$(echo $FileList|cut -d ' ' -f 1|cut -d '_' -f 4|cut -d '.' -f 1)
		echo nev: $nev

		nbFile=$(ls NIRB_approximation*|wc -l)
        	echo !!!!!!!!!! sqrt ... time $tau !!!!!!!!!!!!!
		
	    done
	fi
    done
done
