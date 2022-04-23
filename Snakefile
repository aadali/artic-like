import os
import re
from os import path
import yaml

SNAKEDIR = path.dirname(path.abspath(workflow.snakefile))
configfile: path.join(SNAKEDIR,"config/config.yaml")
with open(path.join(SNAKEDIR,"config/config.yaml"),'r',encoding='utf-8') as f:
    yaml_text = f.read()
    dict_text = yaml.load(yaml_text,Loader=yaml.FullLoader)
    config_paras = dict_text.keys()

if workdir := workflow.overwrite_workdir:
    WORKDIR = path.abspath(workdir)
else:
    WORKDIR = path.abspath(os.curdir)

##### init the other params start...
# now = datetime.now().strftime("%Y%m%d%H%M%S")
MIN_LEN = int(config.get("min_read_len",200))
MAX_LEN = int(config.get("max_read_len",1000))
MIN_READ_Q = int(config.get("min_read_qual",9))
MIN_MAPQ = int(config.get("min_mapq",60))
MIN_DP = int(config.get("min_dp",30))
MIN_QUAL = int(config.get("min_qual",20))
MIN_OVERLAP = int(config.get("min_overlap",300))
HET_SITE = config.get("het_site","more")
MODEL = config.get("model","r941_min_hac_g507")
ANALYSIS_NAME = config.get("analysis_name","test001")
GENOME = path.join(SNAKEDIR,"genome",config.get("what_sample","sars-cov-2"),"sequences.fa")
SNPEFF_DATA = path.join(SNAKEDIR,"genome")
SNPEFF_CONF = path.join(SNAKEDIR,"config/snpEff.config")
FILES_PER_BAR = int(config.get("files_per_bar",2))


##### init the other params end...

##### check input fastq start...


def get_command(directory):
    directory = path.abspath(directory)
    assert path.isdir(directory), print(f"{directory} is not directory")
    contents = os.listdir(directory)
    fq_files = list(filter(lambda x: x.endswith("fastq") or x.endswith("fq"),contents))
    fqgz_files = list(filter(lambda x: x.endswith("fastq.gz") or x.endswith("fq.gz"),contents))
    if len(fq_files) == len(contents):
        return "cat"
    elif len(fqgz_files) == len(contents):
        return f"zcat"
    else:
        raise Exception()


def detect_input(input_fastq, analysis_name=ANALYSIS_NAME):
    if path.isdir(input_fastq):
        abs_path = path.abspath(input_fastq)
        contents = os.listdir(abs_path)
        dirs = list(filter(lambda x: path.isdir(path.join(abs_path,x)),
            contents))

        files = list(filter(lambda x: path.isfile(path.join(abs_path,x)),contents))
        if len(dirs) == len(contents):
            samples = []
            for sample in dirs:
                if (len(os.listdir(path.join(abs_path,sample))) > FILES_PER_BAR
                        and
                        not re.search("unclassified",sample,re.IGNORECASE)):
                    samples.append(sample)
            abs_first_sample = path.join(abs_path,samples[0])
            command = get_command(abs_first_sample)

            # for fastq_pass that contain barcodes directory containing fastq or fastq.gz files
            # abs_path will match the abs path of fastq_pass
            # SAMPLE will match the barcode01 barcode02 or something else directory name
            return (f"{command}", f"{abs_path}/{{SAMPLE}}", "/*", "> ", samples)

        elif len(files) == len(contents):
            samples = [analysis_name]
            command = get_command(abs_path)
            # for one directory which contains fastq or fastq.gz files
            return (f"{command}", f"{abs_path}", "/*", " > ", samples)
        else:
            raise Exception()
    else:
        # if the input_fastq is one fastq or gz file
        samples = [analysis_name]
        return "ln -s", f"{path.abspath(input_fastq)}", " ", "", samples


if "input_fastq" not in config:
    raise Exception("Must specified \"input_fastq\" param in your config file")

INPUT_FASTQ = config.get("input_fastq")
if not path.exists(INPUT_FASTQ) or not INPUT_FASTQ.startswith("/"):
    raise Exception(f"No such file or directory: {INPUT_FASTQ} OR "
                    f"\"input_fastq\" shoud be absolute path")

COMMAND, INPUT, INNER, REDIRECTION, SAMPLES = detect_input(INPUT_FASTQ,ANALYSIS_NAME)
THREADS = int(workflow.cores / len(SAMPLES)) if workflow.cores > len(SAMPLES) else 1
##### check input fastq end...


##### check primer bed file and trim bases params start...
primer_bed = config.get("primer_bed")
abs_primer_bed = path.abspath(primer_bed)

# if primer was specified, no trim; else trim both end 60bp default
PRIMER = path.join(SNAKEDIR,".noprimer.txt")
TRIM_BASES = 0
if path.exists(abs_primer_bed) and path.isfile(abs_primer_bed):
    PRIMER = abs_primer_bed
else:
    TRIM_BASES = config.get("trim_bases",25)
##### check primer bed file and trim bases params end...


#### check config params, the params name from command line may be wrong because of typing
for para in config:
    if para not in config_paras:
        raise Exception("Couldn't understand the parameter: %s" % (para))

# ============================
#       DEFINE RULES
# ============================

rule all:
    input:
        a=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/consensus/{{SAMPLE}}.consensus.fasta",SAMPLE=SAMPLES),
        b=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/figures",SAMPLE=SAMPLES),
        c=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/stat/{{SAMPLE}}.stat",SAMPLE=SAMPLES),
        d=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}.qc.summary.txt",SAMPLE=SAMPLES),
        e=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.longshot.ann.vcf",SAMPLE=SAMPLES),
        f=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.primer_trimmed.report",SAMPLE=SAMPLES),
        g=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.report.snpEff.annotate.txt",SAMPLE=SAMPLES),
        h=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/pangolin/{{SAMPLE}}.lineage_report.csv",SAMPLE=SAMPLES),
        i=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/report/{{SAMPLE}}.report.pdf",SAMPLE=SAMPLES)

rule consensus:
    input:
        a=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/consensus/{{SAMPLE}}.consensus.fasta",SAMPLE=SAMPLES),
        b=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/figures",SAMPLE=SAMPLES),
        c=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/stat/{{SAMPLE}}.stat",SAMPLE=SAMPLES),
        d=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}.qc.summary.txt",SAMPLE=SAMPLES),
        e=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.longshot.ann.vcf",SAMPLE=SAMPLES),
        f=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.primer_trimmed.report",SAMPLE=SAMPLES)

rule annotate:
    input:
        g=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.report.snpEff.annotate.txt",SAMPLE=SAMPLES),
        h=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/pangolin/{{SAMPLE}}.lineage_report.csv",SAMPLE=SAMPLES)

rule report:
    input:
        i=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/report/{{SAMPLE}}.report.pdf",SAMPLE=SAMPLES)

rule fastp:
    # using fastp to filter the raw fastq dropping the low qual reads
    input:
        raw_data=INPUT
    output:
        clean_data=f"{ANALYSIS_NAME}/{{SAMPLE}}/clean_data/{{SAMPLE}}.clean.fastq.gz",
        js=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/clean_data/{{SAMPLE}}.json"),
        html=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/clean_data/{{SAMPLE}}.html"),
        tmp=temporary(f"/tmp/{ANALYSIS_NAME}/{{SAMPLE}}.tmp.fastq.gz") if INPUT.endswith("gz") else temporary(f"/tmp/{ANALYSIS_NAME}/{{SAMPLE}}.tmp.fastq")
    params:
        length_required=MIN_LEN,
        length_limit=MAX_LEN,
        average_qual=MIN_READ_Q
    threads:
        THREADS
    shell:
        f"{COMMAND} {{input.raw_data}}{INNER}  {REDIRECTION} {{output.tmp}} && "
        f"fastp -i {{output.tmp}} -o {{output.clean_data}} "
        # f"--trim_front1 {TRIM_BASES} --trim_tail1 {TRIM_BASES} "
        f"--length_required {{params.length_required}} --length_limit {{params.length_limit}} "
        f"-q 0 -w {{threads}} --average_qual {{params.average_qual}} --html {{output.html}} --json {{output.js}} "

rule nanoplot:
    # nanoplot was used to stat the raw fastq and the clean fastq and make plot
    input:
        raw_data=rules.fastp.output.tmp,
        clean_data=rules.fastp.output.clean_data,
    output:
        outdir=directory(f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot"),
        qc_summary=f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}.qc.summary.txt",
        raw_stats_file=f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}_raw.NanoStats.txt",
        clean_stats_file=f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}_clean.NanoStats.txt",
        finish=f"{ANALYSIS_NAME}/{{SAMPLE}}.nanoplot.finish"  # if nanoplot ruls finished, the file will be touched
    params:
        raw_prefix=f"{{SAMPLE}}_raw.",
        clean_prefix=f"{{SAMPLE}}_clean."
    threads:
        THREADS
    shell:
        f"NanoPlot -t {{threads}} --fastq {{input.raw_data}} --tsv_stats --downsample 2000 --prefix {{params.raw_prefix}}  -o {{output.outdir}} --dpi 300;\n"
        f"NanoPlot -t {{threads}} --fastq {{input.clean_data}} --tsv_stats --downsample 2000 --prefix {{params.clean_prefix}}  -o {{output.outdir}} --dpi 300;\n"
        f"paste {{output.raw_stats_file}} {{output.clean_stats_file}} | cut -f 1,2,4 | sed  '1cMetrics\\traw\\tclean' > {{output.qc_summary}};\n"
        f"touch {{output.finish}}  "


rule map:
    # minimap2 is used to mao reads, the low map_qual align will be dropped
    input:
        clean_data=rules.fastp.output.clean_data,genome=GENOME
    output:
        raw_bam=f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.raw.bam",
        sorted_bam=f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.sorted.bam",
        tempbam=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.temp.sorted.bam")
    params:
        min_mapq=MIN_MAPQ
    threads:
        THREADS
    shell:
        f"minimap2 -ax map-ont -R \"@RG\\tID:{{wildcards.SAMPLE}}\\tSM:{{wildcards.SAMPLE}}\" -t 4 --MD {{input.genome}} {{input.clean_data}} | "
        f"samtools view -@ {{threads}} -bS > {{output.raw_bam}};\n"
        f"samtools view  -@ {{threads}} -q {{params.min_mapq}} -F 4 -bS  {{output.raw_bam}} > {{output.tempbam}} && "
        f"samtools sort -@ {{threads}} -O bam {{output.tempbam}} > {{output.sorted_bam}} && "
        f"samtools index {{output.sorted_bam}}"

rule trim_primers:
    # if primer_bed is used, the alignments in bam file will be trimmed depended on the primer position
    # else just reserve the good alignment whose map_qual is bigger than the specified qual
    input:
        sorted_bam=rules.map.output.sorted_bam,
        primer_bed=PRIMER
    output:
        trimmed_bam=f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.primer_trimmed.bam",
        temp_bam=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.primer_trimmed.temp.bam"),
        trim_log=f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.primer_trimmed.report"
    shell:
        f"python {SNAKEDIR}/scripts/nano_trim_align.py "
        f"{{input.sorted_bam}} {{output.temp_bam}} {{input.primer_bed}} {{output.trim_log}} {MIN_OVERLAP} {TRIM_BASES};\n"
        f"samtools sort {{output.temp_bam}} -o {{output.trimmed_bam}} && "
        f"samtools index {{output.trimmed_bam}}"

rule igv_tools:
    input:
        bamfile=rules.trim_primers.output.trimmed_bam,
        genome=GENOME,
    output:
        per_base_counts=f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.igv.per-pos.bases.depth"
    shell:
        "set +eu && . $(conda info --base)/etc/profile.d/conda.sh  && "
        "conda activate artic-like-igvtools && set -eu && "
        f"igvtools count {{input.bamfile}} stdout {{input.genome}} -w 1 --bases > {{output.per_base_counts}};\n"

rule medaka_consensus:
    # medaka will be called to make hdf file
    input:
        bamfile=rules.trim_primers.output.trimmed_bam
    output:
        hdf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.medaka.hdf"
    params:
        model=MODEL,threads=8
    shell:
        "medaka consensus --model {params.model} --threads {params.threads} {input.bamfile} {output.hdf}"

rule medaka_variant:
    # get gvcf file
    input:
        hdf=rules.medaka_consensus.output.hdf,genome=GENOME,bamfile=rules.trim_primers.output.trimmed_bam
    output:
        vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.medaka.gvcf"
    shell:
        "medaka variant --gvcf {input.genome} {input.hdf} {output.vcf};"

rule bcftools_filter:
    # drop the record that doesn't has alt or whose QUAL is less than specified qual
    input:
        annotate_vcf=rules.medaka_variant.output.vcf
    output:
        vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.medaka.vcf"
    shell:
        f"bgzip -f {{input.annotate_vcf}};"
        f"tabix -p vcf {{input.annotate_vcf}}.gz;"
        f"bcftools filter -e \"ALT='.' || QUAL < {{MIN_QUAL}}\" {{input.annotate_vcf}}.gz  > {{output.vcf}};"

rule medaka_annotate:
    # medaka annotate vcf to get the SR in INFO
    input:
        vcf=rules.bcftools_filter.output.vcf,genome=GENOME,bam=rules.trim_primers.output.trimmed_bam
    output:
        medaka_ann_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.medaka.ann.vcf"
    shell:
        "medaka tools annotate --pad 25 {input.vcf} {input.genome} {input.bam} {output.medaka_ann_vcf}"


rule longshot_annotate:
    # longshot call variants from the potential variants from medaka variants to get AD info
    input:
        vcf=rules.bcftools_filter.output.vcf,genome=GENOME,bam=rules.trim_primers.output.trimmed_bam
    output:
        longshot_ann_vcf=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.longshot.ann.temp.vcf")
    shell:
        f"longshot -F -P 0 --max_cigar_indel 200 --bam {{input.bam}} --ref {{input.genome}} --out {{output.longshot_ann_vcf}} -v {{input.vcf}}"

rule get_per_base_depth:
    input:
        bam=rules.trim_primers.output.trimmed_bam
    output:
        per_base_depth=f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.per-base.depth"
    shell:
        f"samtools depth -a -J {{input.bam}} > {{output.per_base_depth}}"


rule filt_vcf:
    # get the pass vcf and the fail vcf
    input:
        vcf1=rules.medaka_annotate.output.medaka_ann_vcf
    output:
        pass_vcf=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.pass.tmp.vcf"),
        fail_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.fail.vcf"
    shell:
        f"python {SNAKEDIR}/scripts/nano_filt_vcf.py  --min-depth {MIN_DP} "
        f"--min-qual {MIN_QUAL}  {{input.vcf1}}  "
        f"{{output.pass_vcf}} {{output.fail_vcf}} ; \n"

rule calibrate_vcf:
    input:
        longshot_ann_vcf=rules.longshot_annotate.output.longshot_ann_vcf,
        medaka_pass_vcf=rules.filt_vcf.output.pass_vcf,
        per_base_counts=rules.igv_tools.output.per_base_counts,
        genome=GENOME
    output:
        longshot_ann_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.longshot.ann.vcf",
        medaka_pass_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.pass.vcf"
    shell:
        f"python {SNAKEDIR}/scripts/nano_calibrate_vcf.py {{input.per_base_counts}} {{input.genome}} {{input.longshot_ann_vcf}} {{output.longshot_ann_vcf}} "
        f"{{input.medaka_pass_vcf}} {{output.medaka_pass_vcf}}"


rule index:
    input:
        pass_vcf=rules.calibrate_vcf.output.medaka_pass_vcf,
        fail_vcf=rules.filt_vcf.output.fail_vcf
    output:
        pass_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.pass.vcf.gz",
        fail_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.fail.vcf.gz"
    shell:
        f"bgzip {{input.pass_vcf}}  && tabix -p vcf {{output.pass_vcf}} && "
        f"bgzip {{input.fail_vcf}}  && tabix -p vcf {{output.fail_vcf}}"

rule get_low_cover_region:
    # bedtools will be used to get the low coverage region, this region will be N in consensus fasta
    input:
        per_base_dp=f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.per-base.depth"
    output:
        low_cover_bed=f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.low_cover.bed",
        tempbed1=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.low_cover.temp1.bed")
    params:
        # min_dp=MIN_DP,
        awk_query=f"BEGIN{{OFS=\"\\t\"}}{{if($3< {MIN_DP}) {{print $1,$2-1,$2}} }}"
    shell:
        f"awk '{{params.awk_query}}' {{input.per_base_dp}} > {{output.tempbed1}} && "
        f"bedtools merge -i {{output.tempbed1}} > {{output.low_cover_bed}}"

rule pre_mask:
    # the position in low coverage region and fail vcf of reference will be masked by N firstly
    input:
        low_cover_bed=rules.get_low_cover_region.output.low_cover_bed,
        genome=GENOME,
        fail_vcf=rules.index.output.fail_vcf
    output:
        pre_consensus=f"{ANALYSIS_NAME}/{{SAMPLE}}/consensus/{{SAMPLE}}.pre_consensus.fasta"
    shell:
        f"python {SNAKEDIR}/scripts/nano_mask_reference.py {{input.genome}} {{input.low_cover_bed}} {{input.fail_vcf}} {{output.pre_consensus}}"

rule get_consensus:
    # get the final consensus depends on the pass vcf, the alt in pass vcf will be applied ignoring the genotype
    input:
        pre_consensus=rules.pre_mask.output.pre_consensus,
        pass_vcf=rules.index.output.pass_vcf,
        mask_region=rules.get_low_cover_region.output.low_cover_bed
    output:
        consensus=f"{ANALYSIS_NAME}/{{SAMPLE}}/consensus/{{SAMPLE}}.consensus.fasta",
    # finish_log=f"{ANALYSIS_NAME}/{{SAMPLE}}/{{SAMPLE}}.consensus.finished"
    shell:
        f"bcftools consensus -f {{input.pre_consensus}} -m {{input.mask_region}} {{input.pass_vcf}} > {{output.consensus}};\n"
        f"echo ===========================================================;\n"
        f"echo \"            {{wildcards.SAMPLE}} consensus finished        \";\n"
        f"echo ===========================================================;"


rule stat:
    # a simple stat
    input:
        input_fastq=rules.fastp.output.tmp,
        sorted_bam=rules.trim_primers.output.trimmed_bam
    output:
        stat_file=f"{ANALYSIS_NAME}/{{SAMPLE}}/stat/{{SAMPLE}}.stat"
    params:
        awk_query="awk '{if(NR%4==2){print $0}}'"
    shell:
        f"{{params.awk_query}} {{input.input_fastq}} | wc -lm > {{output.stat_file}};\n"
        f"samtools fastq {{input.sorted_bam}} | {{params.awk_query}} | wc -lm >> {{output.stat_file}}"

rule plot:
    # plot the coverage fig
    input:
        per_base_dp=rules.get_per_base_depth.output.per_base_depth
    output:
        fig=directory(f"{ANALYSIS_NAME}/{{SAMPLE}}/figures")
    shell:
        f"mkdir -p {{output.fig}}; python {SNAKEDIR}/scripts/nano_plot.py {{input.per_base_dp}} {{output.fig}} {{wildcards.SAMPLE}}"

rule snpEff:
    # snpEff annotate if {what_sample} is collected in the snpEff database
    input:
        vcf=rules.calibrate_vcf.output.longshot_ann_vcf
    output:
        vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.report.vcf"
    shell:
        f"snpEff -c {SNPEFF_CONF} -dataDir {SNPEFF_DATA} -no-downstream -noStats -no-upstream {config.get('what_sample','sars-cov-2')} "
        f"{{input.vcf}} > {{output.vcf}};\n"

rule handle_snpEff_result:
    # convert the report vcf to variants list
    input:
        report_ann_vcf=rules.snpEff.output.vcf
    output:
        report_variants_list=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.report.snpEff.annotate.txt"
    shell:
        f"python {SNAKEDIR}/scripts/nano_handle_snpeff_results.py {{input.report_ann_vcf}} {{output.report_variants_list}};\n"


rule pangolin:
    # pangolin
    input:
        fasta=rules.get_consensus.output.consensus
    output:
        lineage_report=f"{ANALYSIS_NAME}/{{SAMPLE}}/pangolin/{{SAMPLE}}.lineage_report.csv"
    shell:
        "set +eu && . $(conda info --base)/etc/profile.d/conda.sh && "
        "conda activate artic-like-pangolin && set -eu && "
        f"pangolin {{input.fasta}} --outfile {{output.lineage_report}}"

rule make_report_tex:
    # get the report tex file
    input:
        finish=rules.nanoplot.output.finish,
        consensus=rules.get_consensus.output.consensus,
        pass_vcf=rules.index.output.pass_vcf,
        coverage_fig=rules.plot.output.fig,
        snpeff_vcf=rules.snpEff.output.vcf,
        handle_vcf_result=rules.handle_snpEff_result.output.report_variants_list
    output:
        f"{ANALYSIS_NAME}/{{SAMPLE}}/report/{{SAMPLE}}.report.tex"
    params:
        sample_name="{SAMPLE}",
        what_sample=config.get("what_sample","sars-cov-2")
    shell:
        f"python {SNAKEDIR}/scripts/nano_report.py  {ANALYSIS_NAME} {{params.sample_name}} {{params.what_sample}}  {{output}}"

rule make_report_pdf:
    # xetex will be called to generate the final pdf
    input:
        report_tex=rules.make_report_tex.output
    output:
        report_pdf=f"{ANALYSIS_NAME}/{{SAMPLE}}/report/{{SAMPLE}}.report.pdf"
    shell:
        f"touch {{output.report_pdf}};\n"
        f"{SNAKEDIR}/texlive/2021/bin/x86_64-linux/xelatex -output-directory={ANALYSIS_NAME}/{{wildcards.SAMPLE}}/report {{input.report_tex}};\n"
        f"rm {ANALYSIS_NAME}/{{wildcards.SAMPLE}}/report/{{wildcards.SAMPLE}}.report.log;\n"
        f"rm {ANALYSIS_NAME}/{{wildcards.SAMPLE}}/report/{{wildcards.SAMPLE}}.report.aux;\n"
        f"rm {ANALYSIS_NAME}/{{wildcards.SAMPLE}}.nanoplot.finish"
