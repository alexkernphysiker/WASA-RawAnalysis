for X in `seq 1 1 100`; do
    input="${WMC_DATA}/$1-$X.wmc.data"
    output="${WASA_OUTPUT_DATA}/MC$1$X"
    if [ -e ${input} ]; then
	echo "MC number $X..."
	if [ -e ${output}.root ];then
	    echo "...was already done."
	else
	    wmc_script="run_wmc-$1-$X.sh"
	    if [ -e ${wmc_script} ]; then
		echo "... WMC is stil processing."
	    else
	            scriptname="mc$1$X.sh"
        	    if [ -e ${scriptname} ]; then
	        	echo "...already in process."
	    	    else
			echo "...PROCESSING MC $1..."
	                echo "#!/bin/bash" >> ${scriptname}
	                echo "./rawanalysis MC$1 _ -local -mode mc -fin file:${input} -n ${output} -abort" >> ${scriptname}
	                echo >> ${scriptname}
	                echo "rm -f $PWD/${scriptname}" >> ${scriptname}
	                chmod u+x ${scriptname}
	                ./${scriptname} &> ${output}.log
			echo "...done."
		    fi
            fi
        fi
    fi
done

