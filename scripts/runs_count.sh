rm -f runlist-*.txt
for X in `seq 45883 1 46884`; do
	RUN="NO"
	zezeze=`find ${RUNS_DATA}|grep ${X}.xz`
	if [ "${zezeze}" != "" ]; then
		if [ -e ${zezeze} ]; then
			RUN="YES"
			echo "${X}" >> runlist-present-xz.txt
		fi
	else
		zezeze=`find ${RUNS_DATA}|grep run_${X}.bz2`
		if [ "${zezeze}" != "" ]; then
			if [ -e ${zezeze} ]; then
				RUN="YES"
				echo "${X}" >> runlist-present-bz2.txt
			fi
		fi
	fi
	echo "${X} - ${RUN}" >> runlist-all.txt
	if [ "${RUN}" == "NO" ]; then
		echo "${X}" >> runlist-absent.txt
	else
		echo "${X}" >> runlist-present-all.txt
	fi
done
echo "Present runs"
cat runlist-present-all.txt|wc -l
echo "Absent runs"
cat runlist-absent.txt|wc -l
