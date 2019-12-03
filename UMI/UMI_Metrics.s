PICARD=/iwong/software/picard/build/libs/picard.jar
TMP=/iwong/tmp/


FILES=*sam
for f in $FILES
do
       echo $f
       java -Djava.io.tmpdir=$TMP -Xmx64g -jar $PICARD CollectInsertSizeMetrics I=$f O=${f}_insert_size_metrics.txt H=${f}_insert_size_histogram.pdf
done

FILES=*bam
for f in $FILES
do
       echo $f
        java -Djava.io.tmpdir=$TMP -Xmx64g -jar $PICARD CollectInsertSizeMetrics I=$f O=${f}_insert_size_metrics.txt H=${f}_insert_size_histogram.pdf
done

FILES=*groupedUMI_step07.bam
for f in $FILES
do
       echo $f
       samtools sort $f -o ${f}_sorted.bam
       java -Djava.io.tmpdir=$TMP -Xmx64g -jar $PICARD MarkDuplicates ASO=coordinate I=${f}_sorted.bam O=${f}_marked_duplicates.bam M=${f}_marked_dup_metrics.txt
done


FILES=*_unaligned_step01.bam
for f in $FILES
do
        echo $f
        samtools sort -n $f -o ${f}_readnamesorted.bam
        java -Djava.io.tmpdir=$TMP -Xmx64g -jar $PICARD MarkDuplicates ASO=queryname I=${f}_readnamesorted.bam O=${f}_marked_duplicates.bam M=${f}_marked_dup_metrics.txt
done


FILES=*_unaligned_umi_markedAdpt_step03.bam
for f in $FILES
do
        echo $f
        samtools sort $f -o ${f}_readnamesorted.sam
        java -Djava.io.tmpdir=$TMP -Xmx64g -jar $PICARD MarkDuplicates ASO=queryname I=${f}_readnamesorted.sam O=${f}_marked_duplicates.sam M=${f}_marked_dup_metrics.txt
done


FILES=*_aligned_umi_markedAdpt_bwa_mem_step05.sam
for f in $FILES
do
        echo $f
        samtools sort $f -o ${f}_sorted.sam
        java -Djava.io.tmpdir=$TMP -Xmx64g -jar $PICARD MarkDuplicates ASO=coordinate I=${f}_sorted.sam O=${f}_marked_duplicates.sam M=${f}_marked_dup_metrics.txt
done

FILES=$(ls -lhrt | grep -P ".sam$" | awk '{print $9}')
for f in $FILES
do
        echo $f
        samtools view -S -b $f > ${f}.bam
        rm -rf $f
done
