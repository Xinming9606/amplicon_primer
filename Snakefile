OUTDIR = "ncbi_genome"
DE_DIR = "ncbi_genome_de"

rule all:
    input:
        directory(OUTDIR),
        directory(DE_DIR),
        "gene_seq.fasta",
        "centroid_seq.fasta",
        "aligned.fasta",
        "aligned.aln",
        "diversity.png"
  
checkpoint ncbi_download:
    output:
        directory(OUTDIR)
    params:
        domain=config["domain"],
        genera=config["genera"]
    shell:
        """
        ncbi-genome-download {params.domain} -F 'cds-fasta' -l 'complete' \
        --species-taxids {params.genera} --flat-output -o {output}
        """

rule decompress:
    input:
        directory(OUTDIR)
    output:
        directory(DE_DIR)
    shell:
        """
        mkdir -p {output}
        for gz_file in {input}/*.gz; do
            gzip -d -c "$gz_file" > "{output}/$(basename "$gz_file" .gz)"
        done
        rm -rf {input}
        """   

rule sequence_extract:
    input:
        directory(DE_DIR)
    output:
        "gene_seq.fasta"
    params:
        gene=str(config["gene"])
    shell:
        """
        set +e
        for gcf_file in {input}/*.fna; do
            for gene in {params.gene}; do
                res=$(cat "$gcf_file" | rg -U ">lcl.+gene={params.gene}.+[\w\n]+")
                if [ "$res" != "" ]; then
                echo "$res" >> {output}
                else
                echo "no {params.gene}" >&2
            fi
            done
        done
        """
        
rule vsearch:
  input:
       "gene_seq.fasta"
  output:
       "centroid_seq.fasta"
  shell:
      """
      vsearch -cluster_fast {input} -strand both --id 1 --centroids {output}
      """

rule muscle:
  input:
       "centroid_seq.fasta"
  output:
       "aligned.fasta"
  shell:
      """
      muscle -align {input} -output {output}
      """

rule mv:
   input:
       "aligned.fasta"
   output:
       "aligned.aln"
   shell:
       """
       mv {input} {output}
       """

rule primers:
    input:
        "aligned.aln"
    output:
        "diversity.png"
    params:
        len=config["primer_len"],
        max=config["amplicon_max"],
        min=config["amplicon_min"],
        div=config["div_cut"],
        gc=config["GC_tol"]
    shell:
        "Rscript MSA-primers.R {input} {output} {params.len} {params.max} {params.min} {params.div} {params.gc}"
