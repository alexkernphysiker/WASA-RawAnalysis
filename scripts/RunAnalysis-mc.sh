for X in `seq 1 1 100`; do
    if [ -e $WMC_DATA/$1-$X.wmc.data ]; then
	echo "number $X"
	if [ -e ${WASA_OUTPUT_DATA}/MC$1-$X.root ];then
	    echo "was already done"
	else
            scriptname="run_mc-$X.sh"
            if [ -e ${scriptname} ]; then
        	echo "already in process"
    	    else
		echo "starting"
                echo "#!/bin/bash" >> ${scriptname}
                echo "./rawanalysis MC_$1 -local -mode mc -fin file:$WMC_DATA/$1-$X.wmc.data -n ${WASA_OUTPUT_DATA}/MC$1-$X -abort" >> ${scriptname}
                echo >> ${scriptname}
                echo "rm -f $PWD/${scriptname}" >> ${scriptname}
                chmod u+x ${scriptname}
                ./${scriptname}
            fi
        fi
    fi
done

