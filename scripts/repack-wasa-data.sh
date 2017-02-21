for X in `seq 45873 1 46884`; do
	echo "run #${X}..."
	if [ -e run_${X}.bz2 ]; then
		if [ -e run_${X}.xz ]; then
			echo "is being packed"
		else
			echo "packing..."
			bzcat run_${X}.bz2 | xz > run_${X}.xz
			rm run_${X}.bz2
		fi
	else
		echo "not present"
	fi
done
echo "FINISHED"

