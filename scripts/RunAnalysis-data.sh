for X in `seq 45934 1 46884`; do
	echo "run #${X}..."
	output="${WASA_OUTPUT_DATA}/Data$1${X}"
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
					RUNCAT="bzcat ${RUNS_DATA}/run_${X}.bz2"
				fi
			fi
		fi
		if [ "${RUNCAT}" != "" ]; then
			echo "Analysis of run ${X}"
			scriptname="data$1_${X}.sh"
			if [ -e ${scriptname} ]; then
				echo "...already in process"
			else
				echo "...PROCESSING DATA $1..."
				echo "#!/bin/bash" >> ${scriptname}
				echo "${RUNCAT}|./rawanalysis Data$1 -local -fin cluster: -r ${X} -n ${output} -abort" >> ${scriptname}
				echo "rm -f $PWD/${scriptname}" >> ${scriptname}
				chmod u+x ${scriptname}
				./${scriptname} &> ${output}.log
				echo "...done."
			fi
		fi
	fi
done
