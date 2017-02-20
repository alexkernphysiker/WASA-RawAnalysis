for X in `seq 45873 1 46884`; do
	echo "run #${X}..."
	if [ -e ${RUNS_DATA}/run_${X}.xz ]; then
		echo "...present"
		if [ -e ${WASA_OUTPUT_DATA}/Data_run_${X}.root ]; then
			echo "...has been analyzed"
		else
			scriptname="run_${X}.sh"
			if [ -e ${scriptname} ]; then
				echo "...already in process"
			else
				echo "...starting..."
				echo "#!/bin/bash" >> ${scriptname}
				echo "rm -f ${RUNS_TMP}/run_${X}" >> ${scriptname}
				echo "xzcat ${RUNS_DATA}/run_${X}.xz > ${RUNS_TMP}/run_${X}" >> ${scriptname}
				echo "./rawanalysis Data -local -fin cluster:${RUNS_TMP}/run_${X} -r ${X} -n ${WASA_OUTPUT_DATA}/Data_run_${X} -lf run_${X}.log -abort" >> ${scriptname}
				echo "rm -f ${RUNS_TMP}/run_${X}" >> ${scriptname}
				echo >> ${scriptname}
				echo "rm -f $PWD/${scriptname}" >> ${scriptname}
				chmod u+x ${scriptname}
				./${scriptname}
			fi
		fi
	else
		echo "...not present"
	fi
done
echo "FINISHED"

