for X in `seq 1 1 1`; do
    input="$WMC_DATA/$1-$X.wmc.data"
    output="${WASA_OUTPUT_DATA}/RE$1$X"
    if [ -e ${input} ]; then
	echo "RE number $X..."
	if [ -e ${output}.root ];then
	    echo "...was already done"
	else
	    echo "...PROCESSING RE $1..."
	    scriptname="re$1$X.sh"
	    rm -f ${scriptname}
	    echo "#!/bin/bash" >> ${scriptname}
	    echo "./rawanalysis RE$1 _ -local -mode mc -fin file:${input} -n ${output} -abort" >> ${scriptname}
	    echo >> ${scriptname}
	    echo "rm -f $PWD/${scriptname}" >> ${scriptname}
	    chmod u+x ${scriptname}
	    ./${scriptname} &> ${output}.log
	    echo "...done."
	fi
    fi
done

