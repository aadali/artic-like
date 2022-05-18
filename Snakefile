import os
import re
import gzip
import yaml
from os import path
from subprocess import run

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
MODEL = config.get("model","r941_min_hac_g507")
ANALYSIS_NAME = config.get("analysis_name","test001")
GENOME = path.join(SNAKEDIR,"genome",config.get("what_sample","sars-cov-2"),"sequences.fa")
SNPEFF_DATA = path.join(SNAKEDIR,"genome")
SNPEFF_CONF = path.join(SNAKEDIR,"config/snpEff.config")
FILES_PER_BAR = int(config.get("files_per_bar",2))

if "input_fastq" not in config:
    raise Exception("Must specified \"input_fastq\" param in your config file")

INPUT_FASTQ = config.get("input_fastq")


##### init the other params end...

##### check input fastq start...

def compare_reads_filename(name1, name2):
    """
    confirm which two files is paired reads
    :param name1:
    :param name2:
    :return: a string: {name}@@{read1_file_name}@@{read2_file_name}
    """
    if name1 == name2 or len(name1) != len(name2):
        return ""
    name1 = name1.split("/")[-1]
    name2 = name2.split("/")[-1]
    diff_chars_number = 0
    diff_idx = 0
    for idx, (char1, char2) in enumerate(zip(name1,name2)):
        if char1 == char2:
            continue
        elif char1 == '1' and char2 == "2":
            diff_chars_number += 1
            diff_idx = idx
            read1, read2 = name1, name2
        elif char1 == "2" and char2 == "1":
            diff_chars_number += 1
            diff_idx = idx
            read1, read2 = name2, name1
        else:
            # name1 and name2 are not paired reads because of at the idx position, the chars is not equal
            return ""

    if diff_chars_number == 1:
        return f"{name1[:diff_idx]}@@{read1}@@{read2}"


def check_paired_reads(files):
    """
    for ngs pick paired reads file from a few files
    :param files: a list contains different file names
    :return: one list contains the return of compare_reads_filename
    """
    assert isinstance(files,list)
    sample_paired_reads = set()
    all_reads = []
    files.sort()
    for idx, file1 in enumerate(files):
        if file1 in all_reads:
            continue
        current_len = len(sample_paired_reads)
        find_paired = False
        for file2 in files:
            paired_reads = compare_reads_filename(file1,file2)
            if paired_reads:
                find_paired = True
                sample_paired_reads.add(paired_reads)
                read1_read2 = paired_reads.split("@@")[1:]
                all_reads.extend(read1_read2)
            if find_paired:
                continue
        if len(sample_paired_reads) == current_len:
            raise Exception(f"There is no mate reads file for {file1}")
    return sample_paired_reads


def is_ont_data(file):
    """
    if this reads file is ngs or ont data, check the first 250 reads
    :param file: a read file
    :return:
    """
    assert isinstance(file,str)
    count = 0
    open_file = gzip.open(file,'rt') if file.endswith("gz") else open(file,'r')

    for idx, line in enumerate(open_file):
        if len(line) > 500:
            count += 1
        if count > 0:
            return True
        if count == 0 and idx == 1000:
            return False


def get_command(directory):
    """
    for ont data, if the files are compressed, use zcat or cat.
    :param directory:
    :return:
    """
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
    """
    when input is a directory, there will be mayn situations
    :param input_fastq:
    1. directory
        1). a directory that contains sub directorys, and there are ont reads in the sub diretorys, e.g. the fastq_pass
            that contains some barcode directorys
        2). a directory that contains directly fastqs. These fastqs may be many paired ngs data or some ont fastqs just
            like the barcode in the fastq_pass
    2. file
        a fastq file of merged ont data
    3. two files separated by one comma, this means input is a paired ngs reads. The former is read1 and the later is read2
    :param analysis_name: this analysis task name, get from config.yaml
    :return: a dict depends on the input_fastq
    """
    ont = True
    if not input_fastq.startswith("/"):
        raise Exception(f"input_fastq must be absolute path")

    # when input is a directory, there will be many situations
    if path.isdir(input_fastq):
        abs_path = path.abspath(input_fastq)
        contents = os.listdir(abs_path)
        dirs = list(filter(lambda x: path.isdir(path.join(abs_path,x)),contents))
        files = list(filter(lambda x: path.isfile(path.join(abs_path,x)),contents))

        # for fastq_pass that contain barcodes directory containing fastq or fastq.gz files
        # abs_path will match the abs path of fastq_pass
        # SAMPLE will match the barcode01 barcode02 or something else directory name
        # return (f"{command}", f"{abs_path}/{{SAMPLE}}", "/*", "> ", samples)
        if len(dirs) == len(contents):
            samples = []
            for sample in dirs:
                if (len(os.listdir(path.join(abs_path,sample))) > FILES_PER_BAR
                        and
                        not re.search("unclassified",sample,re.IGNORECASE)):
                    samples.append(sample)
            print(samples)
            abs_first_sample = path.join(abs_path,samples[0])
            command = get_command(abs_first_sample)
            return dict(
                ont=ont,
                command=command,
                fastq_dir=f"{abs_path}/{{SAMPLE}}",
                fastqs="/*",
                direction=">",
                samples=samples
            )

        #  one directory contains fastqs, maybe ont data, maybe ngs paired reads
        elif len(files) == len(contents):
            ont = is_ont_data(path.join(input_fastq,files[0]))

            # for one directory which contains ont fastq or fastq.gz files
            if ont:
                samples = [analysis_name]
                command = get_command(abs_path)
                # return (f"{command}", f"{abs_path}", "/*", " > ", samples)
                return dict(
                    ont=ont,
                    command=command,
                    fastq_dir=abs_path,
                    fastqs="/*",
                    direction=">",
                    samples=samples
                )

            # for one or more ngs paired reads in one directory
            else:
                sample_paired_reads = check_paired_reads(files)
                run(f"mkdir -p /tmp/{analysis_name}",shell=True)
                samples = []
                for paired_reads in sample_paired_reads:
                    sample, read1, read2 = paired_reads.split("@@")
                    samples.append(sample)
                    suffix = "fastq.gz" if read1.endswith("gz") else "fastq"
                    wildcards_name1 = f"/tmp/{ANALYSIS_NAME}/{sample}.reads1.{suffix}"
                    wildcards_name2 = f"/tmp/{ANALYSIS_NAME}/{sample}.reads2.{suffix}"
                    run(f"ln -s -f {path.join(input_fastq,read1)} {wildcards_name1}",shell=True)
                    run(f"ln -s -f {path.join(input_fastq,read2)} {wildcards_name2}",shell=True)

                return dict(
                    ont=ont,
                    samples=samples,
                    reads1=f"/tmp/{ANALYSIS_NAME}/{{SAMPLE}}.reads1.{suffix}",
                    reads2=f"/tmp/{ANALYSIS_NAME}/{{SAMPLE}}.reads2.{suffix}"
                )

    # for one ngs sample input_fastq=read1,read2
    elif "," in input_fastq and not path.isfile(input_fastq):
        ont = False
        run(f"mkdir -p /tmp/{analysis_name}",shell=True)
        read1, read2 = input_fastq.split(",")
        if not (path.isfile(read1) and path.isfile(read2)):
            raise Exception(f"input_fastq: \"{input_fastq}\" is not correct, please check the paired reads")
        suffix = "fastq.gz" if read1.endswith("gz") else "fastq"
        wildcards_name1 = f"/tmp/{ANALYSIS_NAME}/{analysis_name}.reads1.{suffix}"
        wildcards_name2 = f"/tmp/{ANALYSIS_NAME}/{analysis_name}.reads2.{suffix}"
        samples = [analysis_name]
        run(f"ln -s -f {read1} {wildcards_name1}",shell=True)
        run(f"ln -s -f {read2} {wildcards_name2}",shell=True)
        return dict(
            ont=ont,
            samples=samples,
            reads1=f"/tmp/{ANALYSIS_NAME}/{{SAMPLE}}.reads1.{suffix}",
            reads2=f"/tmp/{ANALYSIS_NAME}/{{SAMPLE}}.reads2.{suffix}"
        )

    # one ont sample input_fastq
    elif path.isfile(input_fastq):
        samples = [analysis_name]
        return dict(
            ont=ont,
            command="ln -s",
            fastq_dir=path.abspath(input_fastq),
            fastqs=" ",
            direction="",
            samples=samples
        )

    else:
        raise Exception(f"No such file or directory: {input_fastq}")


input = detect_input(INPUT_FASTQ,ANALYSIS_NAME)
if input['ont']:
    COMMAND = input['command']
    INPUT = input['fastq_dir']
    INNER = input['fastqs']
    REDIRECTION = input['direction']
    SAMPLES = input['samples']
    ngs_label = ""
    READ1 = READ2 = "placeholder"
else:
    COMMAND = INPUT = INNER = REDIRECTION = "placeholder"
    READ1 = input['reads1']
    READ2 = input['reads2']
    SAMPLES = input['samples']
    ngs_label = ".ngs"

THREADS = int(workflow.cores / len(SAMPLES)) if workflow.cores > len(SAMPLES) else 1
##### check input fastq end...


##### check primer bed file and trim bases params start...
primer_bed = config.get("primer_bed")
abs_primer_bed = path.abspath(primer_bed)

# if primer was specified, no trim; else trim both end 60bp default
if path.exists(abs_primer_bed) and path.isfile(abs_primer_bed):
    PRIMER = abs_primer_bed
    TRIM_BASES = 0
else:
    PRIMER = path.join(SNAKEDIR,".noprimer.txt")
    TRIM_BASES = config.get("trim_bases",20) if input['ont'] else 0
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
        a=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/consensus/{{SAMPLE}}{ngs_label}.consensus.fasta",SAMPLE=SAMPLES),
        b=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/figures/{{SAMPLE}}{ngs_label}.coverage.pdf",SAMPLE=SAMPLES),
        c=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/stat/{{SAMPLE}}{ngs_label}.stat",SAMPLE=SAMPLES),
        d=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}{ngs_label}.qc.summary.txt",SAMPLE=SAMPLES),
        # e=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.longshot.ann.vcf",SAMPLE=SAMPLES),
        f=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}{ngs_label}.primer_trimmed.report",SAMPLE=SAMPLES),
        g=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}{ngs_label}.report.snpEff.annotate.txt",SAMPLE=SAMPLES),
        h=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/pangolin/{{SAMPLE}}{ngs_label}.lineage_report.csv",SAMPLE=SAMPLES),
        i=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/report/{{SAMPLE}}{ngs_label}.report.pdf",SAMPLE=SAMPLES)

rule consensus:
    input:
        a=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/consensus/{{SAMPLE}}{ngs_label}.consensus.fasta",SAMPLE=SAMPLES),
        b=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/figures/{{SAMPLE}}{ngs_label}.coverage.pdf",SAMPLE=SAMPLES),
        c=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/stat/{{SAMPLE}}{ngs_label}.stat",SAMPLE=SAMPLES),
        d=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}{ngs_label}.qc.summary.txt",SAMPLE=SAMPLES),
        # e=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.longshot.ann.vcf",SAMPLE=SAMPLES),
        f=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}{ngs_label}.primer_trimmed.report",SAMPLE=SAMPLES),
        g=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}{ngs_label}.report.snpEff.annotate.txt",SAMPLE=SAMPLES),
        h=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/pangolin/{{SAMPLE}}{ngs_label}.lineage_report.csv",SAMPLE=SAMPLES)

rule report:
    input:
        i=expand(f"{ANALYSIS_NAME}/{{SAMPLE}}/report/{{SAMPLE}}{ngs_label}.report.pdf",SAMPLE=SAMPLES)


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
        raw_fig=f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}_raw.LengthvsQualityScatterPlot_dot.png",
        clean_fig=f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}_clean.LengthvsQualityScatterPlot_dot.png",
        raw_stats_file=f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}_raw.NanoStats.txt",
        clean_stats_file=f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}_clean.NanoStats.txt",
        finish=f"{ANALYSIS_NAME}/{{SAMPLE}}.nanoplot.finish"  # if nanoplot ruls finished, the file will be touched
    shell:
        f"python {SNAKEDIR}/scripts/nanoplot.py {{input.raw_data}} {{output.raw_fig}} {{output.raw_stats_file}};\n"
        f"python {SNAKEDIR}/scripts/nanoplot.py {{input.clean_data}} {{output.clean_fig}} {{output.clean_stats_file}};\n"
        f"paste {{output.raw_stats_file}} {{output.clean_stats_file}} | cut -f 1,2,4 | sed  '1cMetrics\\traw\\tclean' > {{output.qc_summary}};\n"
        f"touch {{output.finish}}  "


rule map:
    # minimap2 is used to map reads, the low map_qual align will be dropped
    input:
        clean_data=rules.fastp.output.clean_data,genome=GENOME
    output:
        raw_bam=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.raw.bam"),
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
        f"python {SNAKEDIR}/scripts/trim_align.py "
        f"{{input.sorted_bam}} {{output.temp_bam}} {{input.primer_bed}} {{output.trim_log}} {MIN_OVERLAP} {TRIM_BASES};\n"
        f"samtools sort {{output.temp_bam}} -o {{output.trimmed_bam}} && "
        f"samtools index {{output.trimmed_bam}}"

rule call_snp:
    input:
        bamfile=rules.trim_primers.output.trimmed_bam,
        genome=GENOME,
    output:
        tmp_depth=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.tmp.per-bos.each.bases.depth"),
        per_base_counts=f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.per-pos.each.bases.depth"
    params:
        awk_query="BEGIN{OFS=\"\\t\";print \"#POS\\tA\\tC\\tG\\tT\"}{if(NR>5){print $1,$2,$3,$4,$5}}"
    shell:
        f"python {SNAKEDIR}/scripts/call_snp.py {{input.bamfile}} {{input.genome}} {{output.tmp_depth}};\n"
        f"awk '{{params.awk_query}}' {{output.tmp_depth}}  > {{output.per_base_counts}};\n"

# "set +eu && source $(conda info --base)/etc/profile.d/conda.sh  && "
# "conda activate artic-like-igvtools && set -eu && "
# f"igvtools count {{input.bamfile}} stdout {{input.genome}} -w 1 --bases > {{output.per_base_counts}};\n"

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
        f"python {SNAKEDIR}/scripts/filt_vcf.py  --min-depth {MIN_DP} "
        f"--min-qual {MIN_QUAL}  {{input.vcf1}}  "
        f"{{output.pass_vcf}} {{output.fail_vcf}} ; \n"

rule get_merged_vcf:
    input:
        longshot_ann_vcf=rules.longshot_annotate.output.longshot_ann_vcf,
        medaka_pass_vcf=rules.filt_vcf.output.pass_vcf,
        per_base_counts=rules.call_snp.output.per_base_counts,
        genome=GENOME
    output:
        longshot_ann_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.longshot.ann.vcf",
        medaka_pass_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.pass.vcf"
    shell:
        f"python {SNAKEDIR}/scripts/nano_merge_vcf.py {{input.per_base_counts}} {{input.genome}} {{input.longshot_ann_vcf}} {{output.longshot_ann_vcf}} "
        f"{{input.medaka_pass_vcf}} {{output.medaka_pass_vcf}}"


rule index:
    input:
        pass_vcf=rules.get_merged_vcf.output.medaka_pass_vcf,
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
        f"python {SNAKEDIR}/scripts/mask_reference.py {{input.genome}} {{input.low_cover_bed}} {{input.fail_vcf}} {{output.pre_consensus}}"

rule get_consensus:
    # get the final consensus depends on the pass vcf, the alt in pass vcf will be applied ignoring the genotype
    input:
        pre_consensus=rules.pre_mask.output.pre_consensus,
        pass_vcf=rules.index.output.pass_vcf,
        mask_region=rules.get_low_cover_region.output.low_cover_bed
    output:
        consensus=f"{ANALYSIS_NAME}/{{SAMPLE}}/consensus/{{SAMPLE}}.consensus.fasta",
    params:
        seq_name=f"{ANALYSIS_NAME}-{{SAMPLE}}"
    # finish_log=f"{ANALYSIS_NAME}/{{SAMPLE}}/{{SAMPLE}}.consensus.finished"
    shell:
        f"bcftools consensus -f {{input.pre_consensus}} -m {{input.mask_region}} {{input.pass_vcf}} > {{output.consensus}};\n"
        f"sed -i '1c>{{params.seq_name}}' {{output.consensus}};\n"
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
        awk_query="awk '{if(NR%4==2){print $0}}'",
        command="zcat" if INPUT.endswith("gz") else "cat"
    shell:
        f"{{params.command}} {{input.input_fastq}} | {{params.awk_query}} | wc -lm > {{output.stat_file}};\n"
        f"samtools fastq {{input.sorted_bam}} | {{params.awk_query}} | wc -lm >> {{output.stat_file}}"

rule plot:
    # plot the coverage fig
    input:
        per_base_dp=rules.get_per_base_depth.output.per_base_depth
    output:
        fig=f"{ANALYSIS_NAME}/{{SAMPLE}}/figures/{{SAMPLE}}.coverage.pdf"
    shell:
        f"python {SNAKEDIR}/scripts/plot.py {{input.per_base_dp}} {{output.fig}} {{wildcards.SAMPLE}}"

rule snpEff:
    # snpEff annotate if {what_sample} is collected in the snpEff database
    input:
        vcf=rules.get_merged_vcf.output.longshot_ann_vcf
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
        f"python {SNAKEDIR}/scripts/handle_snpeff_results.py {{input.report_ann_vcf}} {{output.report_variants_list}};\n"


rule pangolin:
    # pangolin
    input:
        fasta=rules.get_consensus.output.consensus
    output:
        lineage_report=f"{ANALYSIS_NAME}/{{SAMPLE}}/pangolin/{{SAMPLE}}.lineage_report.csv"
    shell:
        "set +eu && source $(conda info --base)/etc/profile.d/conda.sh && "
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
        handle_vcf_result=rules.handle_snpEff_result.output.report_variants_list,
        lineage_report=rules.pangolin.output.lineage_report
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

#####################################################################################
################################# FOR NGS DATA ######################################
#####################################################################################
rule ngs_fastp:
    input:
        read1=READ1,
        read2=READ2
    output:
        clean1=f"{ANALYSIS_NAME}/{{SAMPLE}}/clean_data/{{SAMPLE}}.clean.1.fastq.gz",
        clean2=f"{ANALYSIS_NAME}/{{SAMPLE}}/clean_data/{{SAMPLE}}.clean.2.fastq.gz",
        json=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/clean_data/{{SAMPLE}}.json"),
        html=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/clean_data/{{SAMPLE}}.html"),
        qc_summary=f"{ANALYSIS_NAME}/{{SAMPLE}}/nanoplot/{{SAMPLE}}.ngs.qc.summary.txt"
    shell:
        "fastp -i {input.read1} -I {input.read2} -o {output.clean1} -O {output.clean2} --html {output.html} "
        "--json {output.json} --average_qual 20 2> {output.qc_summary}"

rule ngs_map:
    input:
        read1=rules.ngs_fastp.output.clean1,
        read2=rules.ngs_fastp.output.clean2,
        genome=GENOME
    output:
        raw_bam=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.raw.bam"),
        tempbam=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.temp.bam"),
        sorted_bam=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.sorted.bam"),
        rmdup_bam=f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.sorted.rmdup.bam",
    shell:
        f"bwa mem -R \"@RG\\tID:{{wildcards.SAMPLE}}\\tSM:{{wildcards.SAMPLE}}\" {GENOME} {{input.read1}} {{input.read2}} | "
        f"samtools view  -q 60 -F 4 -bS > {{output.raw_bam}};\n "
        f"samtools fixmate -m {{output.raw_bam}}  {{output.tempbam}} && "
        f"samtools sort -O bam {{output.tempbam}} > {{output.sorted_bam}} && "
        f"samtools index {{output.sorted_bam}} && "
        f"samtools markdup -r {{output.sorted_bam}} {{output.rmdup_bam}} && "
        f"samtools index {{output.rmdup_bam}} "

rule ngs_trim_primers:
    input:
        rmdup_bam=rules.ngs_map.output.rmdup_bam,
        primer_bed=PRIMER
    output:
        trimmed_bam=f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.ngs.primer_trimmed.bam",
        temp_bam=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.ngs.primer_trimmed.temp.bam"),
        trim_log=f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.ngs.primer_trimmed.report"
    shell:
        f"python {SNAKEDIR}/scripts/trim_align.py {{input.rmdup_bam}}  {{output.temp_bam}}  {{input.primer_bed}} "
        f"{{output.trim_log}} 50 {TRIM_BASES};\n"
        f"samtools sort {{output.temp_bam}} -o {{output.trimmed_bam}} && "
        f"samtools index {{output.trimmed_bam}}"

rule ngs_call_snp:
    input:
        trimmed_bam=rules.ngs_trim_primers.output.trimmed_bam,
        genome=GENOME
    params:
        awk_query="BEGIN{OFS=\"\\t\";print \"#POS\\tA\\tC\\tG\\tT\"}{if(NR>5){print $1,$2,$3,$4,$5}}"
    output:
        per_base_counts=f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.ngs.per-pos.bases.depth"
    shell:
        "set +eu && source $(conda info --base)/etc/profile.d/conda.sh  && "
        "conda activate artic-like-igvtools && set -eu && "
        f"igvtools count {{input.trimmed_bam}} stdout {{input.genome}} -w 1 --bases | awk '{{params.awk_query}}'  > {{output.per_base_counts}};\n"

rule ngs_call_indel:
    input:
        trimmed_bam=rules.ngs_trim_primers.output.trimmed_bam,
        genome=GENOME
    output:
        indel_bam=f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.ngs.indel.bam",
        vcf=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/aligns/{{SAMPLE}}.vcf")
    threads:
        THREADS
    shell:
        "lofreq indelqual -f {input.genome} --dindel -o {output.indel_bam} {input.trimmed_bam};\n"
        "samtools index {output.indel_bam};\n"
        "lofreq call-parallel --pp-threads {threads} -f {input.genome} --call-indels -o {output.vcf} {output.indel_bam}"

rule ngs_get_merged_vcf:
    input:
        vcf=rules.ngs_call_indel.output.vcf,
        per_pos_bases=rules.ngs_call_snp.output.per_base_counts,
        ref_genome=GENOME,
    output:
        outvcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.all.variants.vcf"
    shell:
        f"python {SNAKEDIR}/scripts/merge_vcf.py {{input.vcf}} {{input.per_pos_bases}}  {{input.ref_genome}} {{output.outvcf}}"

rule ngs_filter_vcf:
    input:
        vcf=rules.ngs_get_merged_vcf.output.outvcf
    output:
        pass_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.pass.vcf",
        faile_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.fail.vcf"
    shell:
        f"python {SNAKEDIR}/scripts/filt_vcf.py --min-depth {MIN_DP} "
        f"--min-qual {MIN_QUAL} {{input.vcf}} {{output.pass_vcf}} {{output.faile_vcf}};"

rule ngs_index:
    input:
        pass_vcf=rules.ngs_filter_vcf.output.pass_vcf,
        fail_vcf=rules.ngs_filter_vcf.output.faile_vcf
    output:
        pass_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.pass.vcf.gz",
        fail_vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.fail.vcf.gz"
    shell:
        f"bgzip {{input.pass_vcf}}  && tabix -p vcf {{output.pass_vcf}} && "
        f"bgzip {{input.fail_vcf}}  && tabix -p vcf {{output.fail_vcf}}"

rule ngs_per_base_depth:
    input:
        trimmed_bam=rules.ngs_trim_primers.output.trimmed_bam
    output:
        per_base_depth=f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.ngs.per-base.depth"
    shell:
        "samtools depth -a -J {input.trimmed_bam} > {output.per_base_depth}"

rule ngs_get_low_cover_region:
    input:
        per_base_dp=rules.ngs_per_base_depth.output.per_base_depth
    output:
        low_cover_bed=f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.low_cover.bed",
        tembed1=temporary(f"{ANALYSIS_NAME}/{{SAMPLE}}/depth/{{SAMPLE}}.low_cover.temp1.bed")
    params:
        awk_query=f"BEGIN{{OFS=\"\\t\"}}{{if($3< {MIN_DP}) {{print $1,$2-1,$2}} }}"
    shell:
        f"awk '{{params.awk_query}}' {{input.per_base_dp}} > {{output.tembed1}} && "
        f"bedtools merge -i {{output.tembed1}} > {{output.low_cover_bed}}"

rule ngs_pre_mask:
    input:
        low_cover_bed=rules.ngs_get_low_cover_region.output.low_cover_bed,
        genome=GENOME,
        fail_vcf=rules.ngs_index.output.fail_vcf
    output:
        pre_consensus=f"{ANALYSIS_NAME}/{{SAMPLE}}/consensus/{{SAMPLE}}.pre_consensus.fasta"
    shell:
        f"python {SNAKEDIR}/scripts/mask_reference.py {{input.genome}} {{input.low_cover_bed}} {{input.fail_vcf}} {{output.pre_consensus}}"

rule ngs_get_consensus:
    input:
        pre_consensus=rules.ngs_pre_mask.output.pre_consensus,
        pass_vcf=rules.ngs_index.output.pass_vcf,
        mask_region=rules.ngs_get_low_cover_region.output.low_cover_bed
    params:
        seq_name=f"{ANALYSIS_NAME}-{{SAMPLE}}"
    output:
        consensus=f"{ANALYSIS_NAME}/{{SAMPLE}}/consensus/{{SAMPLE}}.ngs.consensus.fasta"
    shell:
        f"bcftools consensus -f {{input.pre_consensus}} -m {{input.mask_region}} {{input.pass_vcf}} > {{output.consensus}};\n"
        f"sed -i '1c>{{params.seq_name}}' {{output.consensus}};\n"

rule ngs_plot:
    input:
        per_base_dp=rules.ngs_per_base_depth.output.per_base_depth
    output:
        fig=f"{ANALYSIS_NAME}/{{SAMPLE}}/figures/{{SAMPLE}}.ngs.coverage.pdf"
    shell:
        f"python {SNAKEDIR}/scripts/plot.py {{input.per_base_dp}} {{output.fig}} {{wildcards.SAMPLE}}"

rule ngs_snpEff:
    input:
        vcf=rules.ngs_get_merged_vcf.output.outvcf
    output:
        vcf=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.ngs.report.vcf"
    shell:
        f"snpEff -c {SNPEFF_CONF} -dataDir {SNPEFF_DATA} -no-downstream -noStats -no-upstream {config.get('what_sample','sars-cov-2')} "
        f"{{input.vcf}} >{{output.vcf}}"


rule ngs_stat:
    # a simple stat for ngs
    input:
        read1=READ1,
        read2=READ2,
        bam=rules.ngs_trim_primers.output.trimmed_bam
    output:
        stat_file=f"{ANALYSIS_NAME}/{{SAMPLE}}/stat/{{SAMPLE}}.ngs.stat"
    params:
        awk_query="awk '{if(NR%4==2){print $0}}'",
        command="zcat" if READ1.endswith("gz") else "cat"
    shell:
        f"{{params.command}} {{input.read1}} {{input.read2}} | {{params.awk_query}}  | wc -lm > {{output.stat_file}};\n"
        f"samtools fastq {{input.bam}} | {{params.awk_query}} | wc  -lm >>{{output.stat_file}}"


rule ngs_pangolin:
    input:
        fasta=rules.ngs_get_consensus.output.consensus
    output:
        lineage_report=f"{ANALYSIS_NAME}/{{SAMPLE}}/pangolin/{{SAMPLE}}.ngs.lineage_report.csv"
    shell:
        "set +eu && source $(conda info --base)/etc/profile.d/conda.sh &&"
        "conda activate artic-like-pangolin &&set -eu &&"
        "pangolin {input.fasta} --outfile {output.lineage_report}"

rule ngs_handle_snpEff_result:
    input:
        report_ann_vcf=rules.ngs_snpEff.output.vcf
    output:
        report_variants_list=f"{ANALYSIS_NAME}/{{SAMPLE}}/variants/{{SAMPLE}}.ngs.report.snpEff.annotate.txt"
    shell:
        f"python {SNAKEDIR}/scripts/handle_snpeff_results.py {{input.report_ann_vcf}} {{output.report_variants_list}}"

rule ngs_placehorder:
    output:
        i=f"{ANALYSIS_NAME}/{{SAMPLE}}/report/{{SAMPLE}}.ngs.report.pdf"
    shell:
        "touch  {output.i}"