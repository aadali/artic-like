include {
        trimPrimer;
        getPerBaseDp; 
        plot;
        getLowCoverRegion;
        igvtools;
        snpEff;
        getConsensus;
        pangolin
        } from  './ont.nf'
process fastp {
    storeDir    "$params.directory/$params.analysis_name/$params.analysis_name/01.clean_data"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "${name}"

    input:
        tuple   val(name), val(params)
    output:
        tuple   val(name), path("${name}_clean.*.fastq.gz"), path("${name}.log"), path("${name}.json"), path("${name}.html")
    script:
    def reads2_para = params.reads2==null ? "" : " -I $params.reads2 -O ${name}_clean.reads2.fastq.gz " 
    """fastp -i $params.reads1 -o ${name}_clean.reads1.fastq.gz ${reads2_para} --html ${name}.html --json ${name}.json 2> ${name}.log"""
}

process map {
    storeDir    "$params.directory/$params.analysis_name/$params.analysis_name/02.aligns"
    conda       "$params.conda_path/artic-like"
    tag         "${name}"

    input:
        tuple   val(name), val(clean_reads)
    output:
        tuple   val(name), val(clean_reads), path("${name}.sorted.bam"), path("${name}.sorted.bam.bai")
    script:
    def reads2_para = clean_reads.clean_reads2 == null ? "" : clean_reads.clean_reads2
    """bwa mem -t $params.map_cpu -R "@RG\\tID:${name}\\tSM:${name}"  $projectDir/genomes/$params.virus/sequences.fa  $clean_reads.clean_reads1 ${reads2_para} | \
    samtools view -bS > ${name}.raw.bam && \
    samtools view -h -F 4 -F 2048 -F 256   ${name}.raw.bam | \
    samtools sort -o ${name}.sorted.bam - && \
    samtools index ${name}.sorted.bam"""
}


process ngsStat {
    storeDir    "$params.directory/$params.analysis_name/$params.analysis_name/03.stats"
    conda       "$params.conda_path/artic-like"
    tag         "${name}"

    input:
        tuple   val(name), path(trimmed_bam), path(trimmed_bam_bai), path(fastp_log)
    output:
        tuple   val(name), path("${name}.summary.txt"), path("${name}.bam.stats.txt")
    script:
    """samtools stats -@ 2 ${trimmed_bam} > ${name}.bam.stats.txt && \
    python $projectDir/scripts/ngs_stats.py ${fastp_log} ${name}.bam.stats.txt ${name}.summary.txt"""
}

process ngsGetVariants {
    storeDir    "$params.directory/$params.analysis_name/$params.analysis_name/05.variants"
    conda       "$params.conda_path/freebayes"
    tag         "${name}"

    input:
        tuple   val(name), path(trimmed_bam), path(trimmed_bam_bai)
    output:
        tuple   val(name), path("${name}.raw.vcf")
    script:
    """freebayes -f $projectDir/genomes/$params.virus/sequences.fa -m $params.min_mapq --min-coverage $params.min_dp -F $params.min_snv_freq ${trimmed_bam} > ${name}.raw.vcf"""
}

process ngsFiltVcf {
    storeDir    "$params.directory/$params.analysis_name/$params.analysis_name/05.variants"
    conda       "$params.conda_path/artic-like"
    tag         "${name}"

    input:
        tuple   val(name), path(raw_vcf)
    output:
        tuple   val(name), path("${name}.pass.vcf.gz"), path("${name}.pass.vcf.gz.tbi")
    script:
    """bcftools norm -a -m - -f $projectDir/genomes/$params.virus/sequences.fa -o ${name}.norm.vcf ${name}.raw.vcf && \
    python $projectDir/scripts/ngs_filt_vcf.py --infile ${name}.norm.vcf --outfile ${name}.pass.vcf --min_dp $params.min_dp --min_alt_freq $params.min_snv_freq && \
    bgzip ${name}.pass.vcf && tabix -p vcf ${name}.pass.vcf.gz"""
}

workflow ngs {
    main:
        fastp([params.analysis_name, params]).set { fastpOut } // fastpOut: [name, [clean.*.fastq.gz], log, json, html] or [name, clean.1.fastq.gz, log, json, html]
        fastpOut.map {
            it -> {
                clean = it[1]
                if (clean.size() == 2) {
                    clean_reads = [clean_reads1: clean[0], clean_reads2: clean[1]]
                } else {
                    clean_reads = [clean_reads1: clean, clean_reads2: null]
                }
            }
            return [it[0], clean_reads]
        }.set { mapIn } // mapIn: [name, clean_reads]

        map(mapIn).set{ mapOut }  // mapOut: [name, clean_reads, sorted.bam, sorted.bam.bai]

        mapOut.map { it -> it[0, 2, 3]}.set{ trimPrimerIn }

        trimPrimer(trimPrimerIn).map{it -> it[0..2]}.set { trimPrimerOut }

        trimPrimerOut
            .join(fastpOut.map{it -> it[0, 2]})
            .set { ngsStatIn }

        ngsStat(ngsStatIn).set{ ngsStatOut } // ngsStatOut: [name, summary.txt, bam.stats]

        (ngsGetVariants(trimPrimerOut)  | ngsFiltVcf).set{ ngsFiltVcfOut } // ngsFiltVcfOut: [name, pass.vcf.gz, pass.vcf.gz.tbi]

        getPerBaseDp(trimPrimerOut).set { getPerBaseDpOut } // getPerBaseDpOut: [name, pbd]

        getLowCoverRegion(getPerBaseDpOut).set { getLowCoverRegionOut } // getLowCoverRegionOut: [name, low_cover.bed]

        plot(getPerBaseDpOut).set { plotOut } // plotOut: [name, raw.coverage.png, modified.coverage.png]

        ngsFiltVcfOut
            .join(getLowCoverRegionOut)
            .set { getConsensusIn }

        getConsensus(getConsensusIn).set{ getConsensusOut } // getConsensusOut: [name, consensus.fasta]

        snpEffbin = file("$projectDir/genomes/$params.virus/snpEffectPredictor.bin")

        snpEffOut = snpEffbin.exists() ? snpEff(ngsFiltVcfOut) : ngsFiltVcfOut.map{ it -> it[0, 1]}

        if (params.pango_virus.contains(params.virus)) {
            pangolinOut = pangolin(getConsensusOut)
        } else {
            pangolinOut = getConsensusOut.map { it -> [it[0], null]}
        }

    emit:
        summary = ngsStatOut.map{it -> it[0, 1]} // [name, summary.txt]
        figures = plotOut // [name, rawCov.png, modifedCov.png]
        vcfRpt = snpEffOut // [name, report.vcf[.gz]]
        consensus = getConsensusOut // [name, consensus.fasta]
        pangoRpt = pangolinOut // [name, pagnRpt] || [name, null]
    
}