import os
from os import path


SNAKEDIR = path.dirname(path.abspath(workflow.snakefile))
configfile: path.join(SNAKEDIR,"config/config.yaml")
if workdir := workflow.overwrite_workdir:
    WORKDIR = path.abspath(workdir)
else:
    WORKDIR = path.abspath(os.curdir)

##### check input fastq start...
if "input_fastq" not in config:
    raise Exception("Must specified \"input_fastq\" param in your config file")

INPUT_FASTQ = config.get("input_fastq")
if not path.exists(path.abspath(INPUT_FASTQ)) or len(INPUT_FASTQ.strip()) <= 1:
    raise Exception(f"no such file or directory: {INPUT_FASTQ}")

if path.isdir(INPUT_FASTQ):
    files = os.listdir(path.abspath(INPUT_FASTQ))
    fq_files = list(filter(lambda x: x.endswith("fastq") or x.endswith("fq"),files))
    fqgz_files = list(filter(lambda x: x.endswith("fastq.gz") or x.endswith("fq.gz"),files))
    if len(fq_files) == len(files):
        COMMAND = f"cat {path.abspath(INPUT_FASTQ)}/*"
    elif len(fqgz_files) == len(files):
        COMMAND = f"zcat {path.abspath(INPUT_FASTQ)}/*"
    else:
        raise Exception(f"input_fastq directory must has only fastq or fastq.gz files, check your input_fastq: {INPUT_FASTQ}")

elif path.isfile(INPUT_FASTQ):
    if INPUT_FASTQ.endswith("fastq") or INPUT_FASTQ.endswith("fq"):
        COMMAND = f"cat {path.abspath(INPUT_FASTQ)}"
    elif INPUT_FASTQ.endswith("fastq.gz") or INPUT_FASTQ.endswith("fq.gz"):
        COMMAND = f"zcat {path.abspath(INPUT_FASTQ)}"
    else:
        raise Exception("Couldn't understand your input_fastq file, is it fastq or fq or fastq.gz or fq.gz ?")
else:
    raise Exception()
##### check input fastq end...


##### check primer bed file and trim bases params start...
primer_bed = config.get("primer_bed")
abs_primer_bed = path.abspath(primer_bed)

# if primer was specified, no trim; else trim both end 100bp default
PRIMER = path.join(SNAKEDIR,".noprimer.txt")
TRIM_BASES = 0
if path.exists(abs_primer_bed) and path.isfile(abs_primer_bed):
    PRIMER = abs_primer_bed
else:
    TRIM_BASES = config.get("trim_bases",100)
##### check primer bed file and trim bases params end...


##### init the other params start...
# now = datetime.now().strftime("%Y%m%d%H%M%S")
MIN_LEN = config.get("min_read_len",200)
MAX_LEN = config.get("max_read_len",1000)
MIN_READ_Q = config.get("min_read_qual",9)
MIN_MAPQ = config.get("min_mapq",60)
MAP_THREAD = config.get("map_thread",8)
MIN_DP = config.get("min_dp",30)
MODEL = config.get("model","r941_min_hac_g507")
SAMPLE = config.get("sample_name","test001")
GENOME = path.join(SNAKEDIR,"genome",config.get("what_sample","sars-cov-2"),"sequences.fa")
SNPEFF_DATA = path.join(SNAKEDIR,"genome")
SNPEFF_CONF = path.join(SNAKEDIR,"config/snpEff.config")

##### init the other params end...


# ============================
#       DEFINE RULES
# ============================
rule consensus:
    input:
        a=f"{SAMPLE}/consensus/{SAMPLE}.consensus.fasta",
        b=f"{SAMPLE}/figures/{SAMPLE}.coverage.pdf",
        c=f"{SAMPLE}/stat/{SAMPLE}.stat"


rule annotate:
    input:
        c=f"{SAMPLE}/variants/{SAMPLE}.pass.variants.annotate.txt",
        d=f"{SAMPLE}/variants/{SAMPLE}.fail.variants.annotate.txt",
        e=f"{SAMPLE}/pangolin/{SAMPLE}.lineage_report.csv"

rule report:
    input:
        f=f"{SAMPLE}/report/{SAMPLE}.report.pdf"


rule fastp:
    input:
        raw_data=INPUT_FASTQ
    output:
        clean_data=f"{SAMPLE}/clean_data/{SAMPLE}.clean.fastq.gz",
        js=temporary(f"{SAMPLE}/clean_data/{SAMPLE}.json"),
        html=temporary(f"{SAMPLE}/clean_data/{SAMPLE}.html")
    params:
        length_required=MIN_LEN,
        length_limit=MAX_LEN,
        average_qual=MIN_READ_Q,
        awk_query=f"{{if(NR%4==0 || NR%4==2){{print substr($0, {TRIM_BASES + 1}, length($0) - {2 * TRIM_BASES})}} else {{print $0}}}}"
    shell:
        f"{COMMAND} | awk '{{params.awk_query}}' | fastp --stdin --out1 {{output.clean_data}}  --length_required {{params.length_required}} "
        f"--length_limit {{params.length_limit}} -q 0 -w 8 --average_qual {{params.average_qual}} --html {{output.html}} --json {{output.js}}"

rule map:
    input:
        clean_data=rules.fastp.output.clean_data,genome=GENOME
    output:
        raw_bam=f"{SAMPLE}/aligns/{SAMPLE}.raw.bam",
        sorted_bam=f"{SAMPLE}/aligns/{SAMPLE}.sorted.bam",
        tempbam=temporary(f"{SAMPLE}/aligns/{SAMPLE}.temp.sorted.bam")
    params:
        min_mapq=MIN_MAPQ,
        threads=MAP_THREAD
    shell:
        f"minimap2 -ax map-ont -R \"@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}\" -t {{params.threads}} --MD {{input.genome}} {{input.clean_data}} | "
        f"samtools view -@ {{params.threads}} -bS > {{output.raw_bam}};\n"
        f"samtools view  -@ {{params.threads}} -q {{params.min_mapq}} -F 4 -bS  {{output.raw_bam}} > {{output.tempbam}} && "
        f"samtools sort -@ {{params.threads}} -O bam {{output.tempbam}} > {{output.sorted_bam}} && "
        f"samtools index {{output.sorted_bam}}"

rule trim_primers:
    input:
        sorted_bam=rules.map.output.sorted_bam,
        primer_bed=PRIMER
    output:
        trimmed_bam=f"{SAMPLE}/aligns/{SAMPLE}.primer_trimmed.bam",
        temp_bam=temporary(f"{SAMPLE}/aligns/{SAMPLE}.primer_trimmed.temp.bam")
    shell:
        f"python {SNAKEDIR}/scripts/trim_primer.py {{input.sorted_bam}} {{output.temp_bam}} {{input.primer_bed}};\n"
        f"samtools sort {{output.temp_bam}} -o {{output.trimmed_bam}} && "
        f"samtools index {{output.trimmed_bam}}"

rule medaka_consensus:
    input:
        bamfile=rules.trim_primers.output.trimmed_bam
    output:
        hdf=f"{SAMPLE}/variants/{SAMPLE}.medaka.hdf"
    params:
        model=MODEL,threads=8
    shell:
        f"medaka consensus --model {{params.model}} --threads {{params.threads}} {{input.bamfile}} {{output.hdf}}"

rule medaka_variant:
    input:
        hdf=rules.medaka_consensus.output.hdf,genome=GENOME
    output:
        vcf=f"{SAMPLE}/variants/{SAMPLE}.medaka.vcf"
    shell:
        f"medaka variant {{input.genome}} {{input.hdf}} {{output.vcf}}"


rule mosdepth:
    input:
        rules.trim_primers.output.trimmed_bam
    output:
        per_base_dp=f"{SAMPLE}/mosdepth/{SAMPLE}.per-base.bed.gz"
    shell:
        f"mosdepth {SAMPLE}/mosdepth/{SAMPLE} {{input}}"

rule longshot:
    input:
        vcf=rules.medaka_variant.output.vcf,bam=rules.trim_primers.output.trimmed_bam
    output:
        vcf=f"{SAMPLE}/variants/{SAMPLE}.longshot.vcf"
    shell:
        f"longshot -F -P 0 --bam {{input.bam}} --ref {GENOME} --out {{output.vcf}} -v {{input.vcf}}"

rule filt_vcf:
    input:
        vcf=rules.longshot.output.vcf
    output:
        pass_vcf=f"{SAMPLE}/variants/{SAMPLE}.pass.vcf",
        fail_vcf=f"{SAMPLE}/variants/{SAMPLE}.fail.vcf"
    shell:
        f"python {SNAKEDIR}/scripts/filt_vcf.py {{input.vcf}} {MIN_DP} {{output.pass_vcf}} {{output.fail_vcf}}; \n"

rule index:
    input:
        pass_vcf=rules.filt_vcf.output.pass_vcf,
        fail_vcf=rules.filt_vcf.output.fail_vcf
    output:
        pass_vcf=f"{SAMPLE}/variants/{SAMPLE}.pass.vcf.gz",
        fail_vcf=f"{SAMPLE}/variants/{SAMPLE}.fail.vcf.gz"
    shell:
        f"bgzip {{input.pass_vcf}}  && tabix -p vcf {{output.pass_vcf}} && "
        f"bgzip {{input.fail_vcf}}  && tabix -p vcf {{output.fail_vcf}}"

rule get_low_cover_region:
    input:
        per_base_dp=f"{SAMPLE}/mosdepth/{SAMPLE}.per-base.bed.gz"
    output:
        low_cover_bed=f"{SAMPLE}/mosdepth/{SAMPLE}.low_cover.bed",
        tempbed1=temporary(f"{SAMPLE}/modepth/{SAMPLE}.low_cover.temp1.bed"),
    params:
        min_dp=MIN_DP
    shell:
        f"zcat {{input.per_base_dp}} | awk '$4 < {{params.min_dp}}' | sort -k2 -n > {{output.tempbed1}} &&  "
        f"bedtools merge -i {{output.tempbed1}} > {{output.low_cover_bed}}"

rule pre_mask:
    input:
        low_cover_bed=rules.get_low_cover_region.output.low_cover_bed,
        genome=GENOME,
        fail_vcf=rules.index.output.fail_vcf
    output:
        pre_consensus=f"{SAMPLE}/consensus/{SAMPLE}.pre_consensus.fasta"
    shell:
        f"python {SNAKEDIR}/scripts/mask_reference.py {{input.genome}} {{input.low_cover_bed}} {{input.fail_vcf}} {{output.pre_consensus}}"

rule get_consensus:
    input:
        pre_consensus=rules.pre_mask.output.pre_consensus,
        pass_vcf=rules.index.output.pass_vcf,
        mask_region=rules.get_low_cover_region.output.low_cover_bed
    output:
        consensus=f"{SAMPLE}/consensus/{SAMPLE}.consensus.fasta"
    shell:
        f"bcftools consensus -f {{input.pre_consensus}} -m {{input.mask_region}} {{input.pass_vcf}} > {{output.consensus}}"

rule stat:
    input:
        input_fastq=INPUT_FASTQ,sorted_bam=rules.trim_primers.output.trimmed_bam
    output:
        stat_file=f"{SAMPLE}/stat/{SAMPLE}.stat"
    params:
        awk_query="awk '{if(NR%4==2){print $0}}'"
    shell:
        f"{COMMAND} | {{params.awk_query}} | wc -lm > {{output.stat_file}};\n"
        f"samtools fastq {{input.sorted_bam}} | {{params.awk_query}} | wc -lm >> {{output.stat_file}}"

rule plot:
    input:
        per_base_dp=rules.mosdepth.output.per_base_dp
    output:
        fig=f"{SAMPLE}/figures/{SAMPLE}.coverage.pdf"
    shell:
        f"python {SNAKEDIR}/scripts/plot.py {{input.per_base_dp}} {{output.fig}}"

rule snpEff:
    input:
        pass_vcf=rules.index.output.pass_vcf,
        fail_vcf=rules.index.output.fail_vcf
    output:
        pass_ann_vcf=f"{SAMPLE}/variants/{SAMPLE}.pass.snpEff.annotate.vcf",
        fail_ann_vcf=f"{SAMPLE}/variants/{SAMPLE}.fail.snpEff.annotate.vcf"
    shell:
        f"snpEff -c {SNPEFF_CONF} -dataDir {SNPEFF_DATA} -no-downstream -noStats -no-upstream {config.get('what_sample','sars-cov-2')} {{input.pass_vcf}} > {{output.pass_ann_vcf}};\n"
        f"snpEff -c {SNPEFF_CONF} -dataDir {SNPEFF_DATA} -no-downstream -noStats -no-upstream {config.get('what_sample','sars-cov-2')} {{input.fail_vcf}} > {{output.fail_ann_vcf}};"

rule handle_snpEff_result:
    input:
        pass_ann_vcf=rules.snpEff.output.pass_ann_vcf,
        fail_ann_vcf=rules.snpEff.output.fail_ann_vcf
    output:
        pass_variants_ann=f"{SAMPLE}/variants/{SAMPLE}.pass.variants.annotate.txt",
        fail_variants_ann=f"{SAMPLE}/variants/{SAMPLE}.fail.variants.annotate.txt"
    shell:
        f"python {SNAKEDIR}/scripts/handle_snpeff_results.py {{input.pass_ann_vcf}} {{output.pass_variants_ann}};\n"
        f"python {SNAKEDIR}/scripts/handle_snpeff_results.py {{input.fail_ann_vcf}} {{output.fail_variants_ann}};\n"


rule pangolin:
    input:
        fasta=rules.get_consensus.output.consensus
    output:
        lineage_report=f"{SAMPLE}/pangolin/{SAMPLE}.lineage_report.csv"
    shell:
        f"pangolin {{input.fasta}} --outfile {{output.lineage_report}}"

rule make_report_tex:
    input:
        consensus=rules.get_consensus.output.consensus,
        pass_vcf=rules.index.output.pass_vcf
    output:
        f"{SAMPLE}/report/{SAMPLE}.report.tex"
    params:
        sample_name=SAMPLE,
        what_sample=config.get("what_sample", "sars-cov-2")
    shell:
        f"python {SNAKEDIR}/scripts/report.py  {{params.sample_name}} {{params.what_sample}}  {{output}}"

rule make_report_pdf:
    input:
        report_tex=rules.make_report_tex.output
    output:
        report_pdf=f"{SAMPLE}/report/{SAMPLE}.report.pdf"
    shell:
        f"touch {{output.report_pdf}};\n"
        f"{SNAKEDIR}/texlive/2021/bin/x86_64-linux/xelatex -output-directory={SAMPLE}/report {{input.report_tex}};\n"
        f"rm {SAMPLE}/report/{SAMPLE}.report.log;\n"
        f"rm {SAMPLE}/report/{SAMPLE}.report.aux;\n"






