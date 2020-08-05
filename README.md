# snap_atac_snakemake

snakemake -p -j 99 --cluster-config cluster.json \
--cluster "sbatch  -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &
