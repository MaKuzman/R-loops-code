### This folder contains code used in ROC analysis  

Files were produced using sambamba 0.6.1 for a pallette of bases varying form 1k - 100k ( 1 000, 2 000, 5 000, 10 000, 25 000, 50 000, 100 000).

```
sambamba depth window -F 'mapping_quality > 30 and not duplicate and not failed_quality_control' -o DRIPcA_allsamples_${bases}basescoverage.txt -t 8 --combined -w ${bases} DRIPcAD1770h.sorted.nodup.bam DRIPcAD1780h.sorted.nodup.bam DRIPcAD1790h.sorted.nodup.bam
sambamba depth window -F 'mapping_quality > 30 and not duplicate and not failed_quality_control' -o DRIPcA_D177_${bases}basescoverage.txt -t 8 -w ${bases} DRIPcAD1770h.sorted.nodup.bam
sambamba depth window -F 'mapping_quality > 30 and not duplicate and not failed_quality_control' -o DRIPcA_D178_${bases}basescoverage.txt -t 8 -w ${bases} DRIPcAD1780h.sorted.nodup.bam
sambamba depth window -F 'mapping_quality > 30 and not duplicate and not failed_quality_control' -o DRIPcA_D179_${bases}basescoverage.txt -t 8 -w ${bases} DRIPcAD1790h.sorted.nodup.bam
```
