for X in `seq 1 1 10`; do
    if [ -e $PLUTO_OUTPUT/$1-$X.root ]; then
	echo "WMC number $X..."
	if [ -e $WMC_DATA/$1-$X.wmc.data ]; then
	    echo "...was already done."
	else
	    scriptname="run_wmc-$1-$X.sh"
	    if [ -e ${scriptname} ]; then
		echo "...already in process."
	    else
		echo "...PROCESSING WMC $1..."
		echo "#!/bin/bash" >> ${scriptname}
		echo "mkdir $PWD/wmc_tmp_$1-$X" >> ${scriptname}
		echo "cp $PWD/WMC/* $PWD/wmc_tmp_$1-$X" >> ${scriptname}
		echo "cd $PWD/wmc_tmp_$1-$X" >> ${scriptname}
		echo "./wmc.sh $PLUTO_OUTPUT/$1-$X.root $WMC_DATA/$1-$X.wmc.data" >> ${scriptname}
		echo "cd $PWD" >> ${scriptname}
		echo >> ${scriptname}
		echo "rm -r $PWD/wmc_tmp_$1-$X" >> ${scriptname}
		echo "rm -f $PWD/${scriptname}" >> ${scriptname}
		chmod u+x ${scriptname}
		./${scriptname}
		echo "...done."
	    fi
	fi
    fi
done
