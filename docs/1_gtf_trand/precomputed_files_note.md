## NOTE ON SINGLE GTF PRECOMPUTED FILES:

The single GTF TranD pre-computed output was created via the following commands. More information on how to run TranD can be found in the [User Guide](https://github.com/McIntyre-Lab/TranD/wiki/User-Guide).

The "input.gtf" is located on the pre-computed files page next to the output.



```
##GENE:
trand \
input.gtf \
-o outputdirectory \
-e gene \
-k \
-f \
-n num_cpus

```

```
##PAIRWISE:
trand \
input.gtf \
-o outputdirectory \
-f \
-n num_cpus

```
