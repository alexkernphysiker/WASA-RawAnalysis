for X in `seq 45934 1 46884`; do
	echo "run #${X}..."
	output="${WASA_OUTPUT_DATA}/DataAll${X}"
	if [ -e ${output}.root ]; then
		echo "... is already done."
	else
		RUNCAT=""
		#command find is in cycle because file list may change after each analysis
		zezeze=`find ${WASA_PRESELECTION}|grep ${X}.root`
		if [ "${zezeze}" != "" ]; then
			if [ -e ${zezeze} ]; then
				RUNCAT="cat ${zezeze}"
			fi
		fi
		if [ "${RUNCAT}" != "" ]; then
			echo "Analysis of run ${X}"
			scriptname="dataall_${X}.sh"
			if [ -e ${scriptname} ]; then
				echo "...already in process"
			else
				echo "...PROCESSING PRESELECTED DATA ..."
				echo "#!/bin/bash" >> ${scriptname}
				echo "${RUNCAT}|./rawanalysis DataAll -local -fin cluster: -r ${X} -n ${output} -abort" >> ${scriptname}
				echo "rm -f $PWD/${scriptname}" >> ${scriptname}
				chmod u+x ${scriptname}
				./${scriptname} &> ${output}.log
				echo "...done."
			fi
		fi
	fi
done
