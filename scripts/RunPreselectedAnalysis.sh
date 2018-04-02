for X in `seq 45934 1 46884`; do
	echo "run #${X}..."
	output="${WASA_OUTPUT_DATA}/DataAll${X}"
	if [ -e ${output}.root ]; then
		echo "... is already done."
	else
		zezeze=${WASA_PRESELECTION}/run_${X}.bz2
		echo ${zezeze}
		if [ -e ${zezeze} ]; then
			scriptname="dataall_${X}.sh"
			lock="preselection_${X}.sh"
			if [ -e ${lock} ]; then
				echo "...still locked"
			else
				if [ -e ${scriptname} ]; then
					echo "...already in process"
				else
					echo "...PROCESSING PRESELECTED DATA ..."
					echo "#!/bin/bash" >> ${scriptname}
					echo "bzcat ${zezeze}|./rawanalysis DataAll -local -fin cluster: -r ${X} -n ${output} -abort" >> ${scriptname}
					echo "rm -f $PWD/${scriptname}" >> ${scriptname}
					chmod u+x ${scriptname}
					./${scriptname} &> ${output}.log
					echo "...done."
				fi
			fi
		fi
	fi
done
