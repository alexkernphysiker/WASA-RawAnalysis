for X in `seq 45934 1 46884`; do
	echo "run #${X}..."
	output="${WASA_PRESELECTION}/preselection_${X}"
	if [ -e ${output}.root ]; then
		echo "... is already done."
	else
		RUNCAT=""
		#command find is in cycle because file list may change after each analysis
		zezeze=`find ${RUNS_DATA}|grep ${X}.xz`
		if [ "${zezeze}" != "" ]; then
			if [ -e ${zezeze} ]; then
				RUNCAT="xzcat ${zezeze}"
			fi
		else
			zezeze=`find ${RUNS_DATA}|grep ${X}.bz2`
			if [ "${zezeze}" != "" ]; then
				if [ -e ${zezeze} ]; then
					RUNCAT="bzcat ${zezeze}"
				fi
			fi
		fi
		if [ "${RUNCAT}" != "" ]; then
			echo "Analysis of run ${X}"
			scriptname="preselection_${X}.sh"
			if [ -e ${scriptname} ]; then
				echo "...already in process"
			else
				echo "... DATA PRESELECTION ..."
				echo "#!/bin/bash" >> ${scriptname}
				echo "${RUNCAT}|./rawanalysis PRESELECTION _ -local -fin cluster: -r ${X} -n ${output} -abort" >> ${scriptname}
				echo "rm -f $PWD/${scriptname}" >> ${scriptname}
				chmod u+x ${scriptname}
				./${scriptname} &> ${output}.log
				echo "...done."
			fi
		fi
	fi
done
