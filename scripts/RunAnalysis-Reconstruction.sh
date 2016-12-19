for X in `seq 1 1 100`; do
    if [ -e $WMC_DATA/$1-$X.wmc.data ]; then
	echo "number $X"
	if [ -e $PWD/../Sorting/RE$1-$X.root ];then
	    echo "was already done"
	else
	    scriptname="run_reconstr-$X.sh"
	    rm -f ${scriptname}
	    echo "#!/bin/bash" >> ${scriptname}
	    echo "cd $PWD/../Sorting" >> ${scriptname}
	    echo "./main RE_$1 -mode mc -fin file:$WMC_DATA/$1-$X.wmc.data -n RE$1-$X -abort" >> ${scriptname}
	    echo >> ${scriptname}
	    echo "rm -f $PWD/${scriptname}" >> ${scriptname}
	    chmod u+x ${scriptname}
	    ./${scriptname}
	fi
    fi
done

