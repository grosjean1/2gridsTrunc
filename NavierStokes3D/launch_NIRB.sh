#! /bin/bash

#Elise Grosjean
#lance NIRB + POD Greedy Parabolique avec mu=1, a comparer avec erreur FEM fine, ex: nev=10 (H=h*2 et h=H²) finetime=0.01 0.02 0.05 0.1
# uses :
#           - script  NirbGreedyTestNIRBOKParabolique.py with mu=1
#           - meshio-convert (vtk) to convert vtu files to vtk...

finetime="0.01 0.02 0.05 0.1"
finemesh="140 70 30 15"

#grossier
coarsetimediv="0.02 0.04 0.1 0.2" #hdiv2
coarsetimesqrt="0.1 0.1414 0.22 0.32" #sqrth

coarsemeshdiv="70 36 15 7"
coarsemeshsqrt="15 10 7 5"
# 1-9 5-1 5-2 5-5
for param2 in 1  
do
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
		#python3.11 NirbTradiGreedy.py $nev $tau sqrth $taucoarse $param1 $param2
		python3.11 NirbGreedy.py $nev $tau sqrth $taucoarse $param1 $param2

		FileList=$(ls NIRB_approximation*) 

		#nev=$(echo $FileList|cut -d ' ' -f 1|cut -d '_' -f 4|cut -d '.' -f 1)
		echo nev: $nev

		nbFile=$(ls NIRB_approximation*|wc -l)
		
		for (( c=0; c<$nbFile; c++ ))
		do  
		    echo meshio convert, file number: $c
		    #/home/grosjean/.local/bin/meshio-convert
		    meshio-convert NIRB_approximation_${c}_${nev}.vtu NIRB_approximation_${c}_${nev}.vtk
		    meshio-ascii NIRB_approximation_${c}_${nev}.vtk
		done
	
		## error computation
		
		nnref=$(echo $finemesh |cut -d ' ' -f $values)
		nnrefc=$(echo $coarsemeshsqrt |cut -d ' ' -f $values)
		echo sizemesh $nnref
		echo " \n sqrt h=H^2: \n" >> error.txt

		#FreeFem++-nw Crank_EulerFinNIRB.edp -tau $tau -nnref $nnref -Param1 $param1 -Param2 $param2 -nev $nev
		#FreeFem++-nw Crank_EulerFindiff.edp -tau $tau -nnref $nnref -Param1 $param1 -Param2 $param2 -theta $theta
		#FreeFem++-nw Crank_EulerFindiff.edp -tau $taucoarse -nnref $nnrefc -Param1 $param1 -Param2 $param2 -theta $theta
       		FreeFem++-nw Crank_EulerFinNIRB.edp -tau $tau -nnref $nnref -Param 1 -nev $nev
		#rm NIRB_app*
		echo !!!!!!!!!! sqrt ... time $tau !!!!!!!!!!!!!
		
	    done
	done
    done
done
