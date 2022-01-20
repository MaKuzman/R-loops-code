#!/bin/bash
#PBS -N strands

#PBS -m a
#PBS -M maja@kuzman.org

#PBS -q MASTER
#PBS -l select=1:ncpus=2:mem=50G

cd $PBS_O_WORKDIR

SPAN=10
S=$x
R=${x%%.sorted.nodup.forward.q30.bedgraph.sortd}.forward.wig
rm -f ${R}
head -45 ${S} | egrep "^browser|^track" > ${R}
grep "^chr" ${S} | cut -f1 | sort -u > ${x%%.sorted.nodup.forward.q30.bedgraph.sortd}fchr.list
cat ${x%%.sorted.nodup.forward.q30.bedgraph.sortd}fchr.list | while read C
do
    echo "variableStep chrom=${C} span=${SPAN}" >> ${R}
    awk '{if (match($1,"^'"${C}"'$")) { print } }' ${S} | sort -k2n | awk '
{
    printf "%d\t%g\n", $2+1, $4
}
' >> ${R}
done
rm ${x%%.sorted.nodup.forward.q30.bedgraph.sortd}fchr.list 
