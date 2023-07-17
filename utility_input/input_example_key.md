# Example Input Legend

- **dmel-all-r6.17.gtf:** A Single GTF reference file for _Drosophila melanogaster_. Taken from the
[Precomputed Files](https://github.com/McIntyre-Lab/TranD/wiki/Precomputed-Files) section of the wiki. Often used for utilities with **1 GTF input.**

- **dmel_pairwise_transcript_distance.csv:** The pairwise transcript distance file from TranD output for _Drosophila melanogaster_. Taken from the [Precomputed Files](https://github.com/McIntyre-Lab/TranD/wiki/Precomputed-Files) section of the wiki. A **1 GTF** TranD output file. 

- **flair_vs_isoseq.csv:** The pairwise transcript distance file from TranD output comparing FLAIR and Isoseq3 Cluster long read methods. A **2 GTF** minimum pairwise distance TranD output file. It was generated from _D. sim_ data, so the full suffix is a `sim_FLAIR` and `sim_isoseq3cluster`  It has been subset to allow for storage in the github (<100 MB).

- **subset_gtf_include_list.txt:** The include list used in the example for [subset_gtf.py](https://github.com/McIntyre-Lab/TranD/wiki/Utility-Descriptions-(with-Examples)#subset_gtfpy).

- **subset_pd_include_list1.txt** and **subset_pd_include_list2.txt:** The include list used in the example for [subset_trand_pairwise_transcript_distance.csv](https://github.com/McIntyre-Lab/TranD/wiki/Utility-Descriptions-(with-Examples)#subset_trand_pairwise_transcript_distancepy).
