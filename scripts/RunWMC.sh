mkdir $PWD/wmc_tmp_$1
cp $PWD/WMC/* $PWD/wmc_tmp_$1
for X in `seq 1 1 100`; do
    if [ -e $PLUTO_OUTPUT/$1-$X.root ]; then
	echo "number $X"
	if [ -e $WMC_DATA/$1-$X.wmc.data ]; then
	    echo "was already done"
	else
	    scriptname="run_wmc-$1-$X.sh"
	    if [ -e ${scriptname} ]; then
		echo "already in process"
	    else
		echo "#!/bin/bash" >> ${scriptname}
		echo "cd $PWD/wmc_tmp_$1" >> ${scriptname}
		echo "./wmc.sh $PLUTO_OUTPUT/$1-$X.root $WMC_DATA/$1-$X.wmc.data" >> ${scriptname}
		echo >> ${scriptname}
		echo "rm -f $PWD/${scriptname}" >> ${scriptname}
		chmod u+x ${scriptname}
		./${scriptname}
	    fi
	fi
    fi
done
rm -rf $PWD/wmc_tmp_$1