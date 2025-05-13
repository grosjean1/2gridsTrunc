#! /bin/bash
# Elise Grosjean
# Initialization

#fine setting

finetime="0.01 0.02 0.05 0.1"
finemesh="140 70 30 15"

#coarse setting
coarsetime="0.02 0.04 0.1 0.2" #hdiv2
coarsemesh="70 36 15 7"

################################
#fine setting

finetime="0.0025"
finemesh="560"

theta=1. #euler

for values in 1 
do

    for param1 in 8
    do

	for param2 in 8 9
	do
	    
	tau=$(echo $finetime | cut -d ' ' -f $values)
	echo time $tau
	nnref=$(echo $finemesh |cut -d ' ' -f $values)
	echo sizemesh $nnref
	mkdir -p  /Users/elisegrosjean/Domaine_tronque_parabolic/newcase/FineSnapshotsRef/$tau/$param1-$param2/
	FreeFem++-nw Crank_Eulerinit.edp -tau $tau -nnref $nnref -Param1 $param1 -Param2 $param2 -theta $theta -FineRef 1
	mv Snapshot* /Users/elisegrosjean/Domaine_tronque_parabolic/newcase/FineSnapshotsRef/$tau/$param1-$param2/
	done
    done
done

