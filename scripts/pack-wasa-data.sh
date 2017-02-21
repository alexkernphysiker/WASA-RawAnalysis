for X in `seq 45873 1 46884`; do
	echo "run #${X}..."
	if [ -e run_${X} ]; then
		if [ -e run_${X}.xz ]; then
			echo "is being packed"
		else
			echo "packing..."
			xz -T0 -9 run_${X}
		fi
	else
		echo "not present"
	fi
done
echo "FINISHED"

