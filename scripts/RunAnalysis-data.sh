for X in `seq 45934 1 46884`; do
	echo "run #${X}..."
	if [ -e ${RUNS_DATA}/run_${X}.xz ]; then
		echo "...present"
		if [ -e ${WASA_OUTPUT_DATA}/Data_run_${X}.root ]; then
			echo "...has been analyzed"
		else
			scriptname="run_${X}$1.sh"
			if [ -e ${scriptname} ]; then
				echo "...already in process"
			else
				echo "...starting..."
				echo "#!/bin/bash" >> ${scriptname}
				echo "xzcat ${RUNS_DATA}/run_${X}.xz|./rawanalysis Data -local -fin cluster: -r ${X} -n ${WASA_OUTPUT_DATA}/Data_run_${X}$1 -lf run_${X}$1.log -abort" >> ${scriptname}
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

