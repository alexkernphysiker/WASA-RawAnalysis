rm -f runlist-*.txt
RUN_FILES=`find ${RUNS_DATA}`
for X in `seq 45934 1 46884`; do
	RUN="NO"
	zezeze=`echo ${RUN_FILES}|grep ${X}.xz`
	if [ "${zezeze}" != "" ]; then
		RUN="YES"
	else
		zezeze=`echo ${RUN_FILES}|grep run_${X}.bz2`
		if [ "${zezeze}" != "" ]; then
			RUN="YES"
		fi
	fi
	if [ "${RUN}" == "NO" ]; then
		echo "${X}" >> runlist-absent.txt
	else
		echo "${X}" >> runlist-present.txt
	fi
done
echo "Present runs"
cat runlist-present.txt|wc -l
echo "Absent runs"
cat runlist-absent.txt|wc -l
