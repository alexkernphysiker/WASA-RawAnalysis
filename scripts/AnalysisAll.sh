for R in bound1-2g bound2-2g bound3-2g bound1-6g bound2-6g bound3-6g He3eta-gg He3eta-6g ppn_qf_ He3pi0 He3pi0pi0 He3pi0pi0pi0; do
for X in `seq 1 1 10`; do
    input="${WMC_DATA}/$R-$X.wmc.data"
    output="${WASA_OUTPUT_DATA}/MC$R$X$1"
    if [ -e ${input} ]; then
        echo "MC number $X..."
        if [ -e ${output}.root ];then
            echo "...was already done."
        else
            wmc_script="run_wmc-$R-$X-$1.sh"
            if [ -e ${wmc_script} ]; then
                echo "... WMC is stil processing."
            else
                    scriptname="mc$R$X$1.sh"
                    if [ -e ${scriptname} ]; then
                        echo "...already in process."
                    else
                        echo "...PROCESSING MC $R $X $1..."
                        echo "#!/bin/bash" >> ${scriptname}
                        echo "./rawanalysis MC$R $1 -local -mode mc -fin file:${input} -n ${output} -abort" >> ${scriptname}
                        echo >> ${scriptname}
                        echo "rm -f $PWD/${scriptname}" >> ${scriptname}
                        chmod u+x ${scriptname}
                        ./${scriptname} &> ${output}.log
                        echo "...done."
                    fi
            fi
        fi
    fi
done
done
for X in `seq 45934 1 46884`; do
	echo "run #${X}..."
	output="${WASA_OUTPUT_DATA}/DataAll${X}$1"
	if [ -e ${output}.root ]; then
		echo "... is already done."
	else
		zezeze=${WASA_PRESELECTION}/run_${X}.presel
		echo ${zezeze}
		if [ -e ${zezeze} ]; then
			scriptname="data_${X}$1.sh"
			lock="preselection_${X}.sh"
			if [ -e ${lock} ]; then
				echo "...still locked"
			else
				if [ -e ${scriptname} ]; then
					echo "...already in process"
				else
					echo "...PROCESSING PRESELECTED DATA ..."
					echo "#!/bin/bash" >> ${scriptname}
	echo "./rawanalysis DataAll $1 -local -mode raw -fin file:${zezeze} -n ${output} -abort" >> ${scriptname}
					echo "rm -f $PWD/${scriptname}" >> ${scriptname}
					chmod u+x ${scriptname}
					./${scriptname} &> ${output}.log
					echo "...done."
				fi
			fi
		fi
	fi
done
