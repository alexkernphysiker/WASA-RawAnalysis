for X in `seq 45934 1 46884`; do
	echo "run #${X}..."
	if [ -e ${WASA_OUTPUT_DATA}/Data_run_${X}.root ]; then
		echo "...has been analyzed"
	else
		RUNCAT=""
		if [ -e ${RUNS_DATA}/run_${X} ]; then
			RUNCAT="cat ${RUNS_DATA}/run_${X}"
		fi
		if [ -e ${RUNS_DATA}/run_${X}.xz ]; then
			RUNCAT="xzcat ${RUNS_DATA}/run_${X}.xz"
		fi
		if [ -e ${RUNS_DATA}/run_${X}.bz2 ]; then
			RUNCAT="bzcat ${RUNS_DATA}/run_${X}.bz2"
		fi
		echo "${RUNCAT}"
		if [ "${RUNCAT}" == "" ]; then
			echo "none."
		else
			echo "exists."
			scriptname="run_${X}$1.sh"
			if [ -e ${scriptname} ]; then
				echo "...already in process"
			else
				echo "...starting..."
				echo "#!/bin/bash" >> ${scriptname}
				echo "${RUNCAT}|./rawanalysis Data -local -fin cluster: -r ${X} -n ${WASA_OUTPUT_DATA}/Data_run_${X}$1 -lf run_${X}$1.log -abort" >> ${scriptname}
				echo "rm -f $PWD/${scriptname}" >> ${scriptname}
				chmod u+x ${scriptname}
				./${scriptname}
			fi
		fi
	fi
done
echo "FINISHED"

