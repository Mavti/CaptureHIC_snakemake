bsub -J "CHICprocessing" -oo "runme_0-CHICprocessing_pipeline.o" "snakemake -j 999 --cluster-config script/cluster_config.json  -s script/0-CHIC_Processing.snakefile --use-conda --cluster 'bsub -o {cluster.output} -e {cluster.error} -n {threads} -M {cluster.memory} -R {cluster.resources}'"