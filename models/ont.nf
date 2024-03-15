include {q2p; detectLongReads} from './utils.nf'

process ontPreprocess {
    storeDir    "$params.directory/$params.analysis_name/${input_fastq['name']}/01.clean_data"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "${input_fastq['name']}"

    input:
        val input_fastq
    output:
        tuple   val("${input_fastq['name']}"), path("${input_fastq['name']}.raw.${input_fastq['ext']}")
    script:
    """${input_fastq['cmd']} ${input_fastq['fastqs']}  ${input_fastq['direction']}  ${input_fastq['name']}.raw.${input_fastq['ext']}"""
}

process ontFiltlong {
    storeDir    "$params.directory/$params.analysis_name/$name/01.clean_data"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "$name"

    input:
        tuple   val(name), path(raw_fastq)
    output:
        tuple   val(name), path("${name}.clean.fastq.gz")
    script:
    def pred = q2p(params.min_read_qual) // the --min_mean_q of filtlong should be phred
    """ filtlong --min_length $params.min_read_len  --max_length $params.max_read_len --min_mean_q $pred $raw_fastq | gzip -c > ${name}.clean.fastq.gz """ 
}

process ontMap {
    storeDir    "$params.directory/$params.analysis_name/$name/02.aligns"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "$name"

    input:
        tuple   val(name), path(clean_fastq)
    output:
        tuple   val(name), path("${name}.sorted.bam"), path("${name}.sorted.bam.bai")
    script:
    """minimap2 -ax map-ont -R "@RG\\tID:${name}\\tSM:${name}" -t $params.map_cpu --MD $projectDir/genomes/$params.virus/sequences.fa ${clean_fastq} | \
    samtools view  -@ 4 -bS > ${name}.raw.bam && \
    samtools sort -@ 4 -o ${name}.sorted.bam ${name}.raw.bam && \
    samtools index ${name}.sorted.bam && \
    rm ${name}.raw.bam"""
}

process trimPrimer {
    storeDir    "$params.directory/$params.analysis_name/$name/02.aligns"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "$name"

    input:
        tuple   val(name), path(bam), path(bai)
    output:
        tuple   val(name), path("${name}.primer_trimmed.bam"), path("${name}.primer_trimmed.bam.bai"), path("${name}.primer_trimmed.report.tsv")
    script:
    def PRIMER = params.primer_bed == null ? "./.noprimer.txt" : params.primer_bed
    def OVERLAP = params.long_reads == null ? 20 : params.min_overlap
    def TRIMBASES = params.long_reads == null ? 0 : params.trim_bases
    """python $projectDir/scripts/trim_align.py $bam tmp.bam $PRIMER ${name}.primer_trimmed.report.tsv $OVERLAP $TRIMBASES && \
    samtools sort -o ${name}.primer_trimmed.bam tmp.bam && \
    samtools index ${name}.primer_trimmed.bam && \
    rm tmp.bam"""
}

process ontStat {
    storeDir    "$params.directory/$params.analysis_name/$name/03.stats"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "$name"

    input:
        tuple   val(name), path(bam), path(bai), path(raw_fastq), path(clean_fastq)
    output:
        tuple   val(name), 
                path("${name}.raw.LengthvsQualityScatterPlot_dot.png"), 
                path("${name}.clean.LengthvsQualityScatterPlot_dot.png"), 
                path("${name}.mapped.LengthvsQualityScatterPlot_dot.png"),
                path("${name}.summary.txt")
    script:
    """samtools fastq $bam > ${name}.mapped.fastq && \
    python $projectDir/scripts/nanoplot.py ${raw_fastq} ${name}.raw.LengthvsQualityScatterPlot_dot.png ${name}.raw.stats.txt && \
    python $projectDir/scripts/nanoplot.py ${clean_fastq} ${name}.clean.LengthvsQualityScatterPlot_dot.png ${name}.clean.stats.txt && \
    python $projectDir/scripts/nanoplot.py ${name}.mapped.fastq  ${name}.mapped.LengthvsQualityScatterPlot_dot.png ${name}.mapped.stats.txt && \
    echo -e "SampleName\t${params.analysis_name}\t$name" > a.txt && \
    paste ${name}.raw.stats.txt ${name}.clean.stats.txt ${name}.mapped.stats.txt | cut -f 1,2,4,6 | sed '1cMetrics\\traw\\tclean\\tmapped' > b.txt && \
    cat a.txt b.txt > ${name}.summary.txt && rm a.txt b.txt """
}

process getPerBaseDp {
    storeDir    "$params.directory/$params.analysis_name/$name/04.depth"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "$name"

    input:
        tuple   val(name), path(bam), path(bai)
    output:
        tuple   val(name), path("${name}.per-base.depth")
    script:
    """samtools depth -a -J --min-MQ $params.min_mapq $bam > ${name}.per-base.depth"""
}

process plot {
    storeDir    "$params.directory/$params.analysis_name/$name/04.depth"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "$name"

    input:
        tuple   val(name), path(pbd)
    output:
        tuple   val(name), path("*.raw.coverage.png"), path("*.modified.coverage.png")
    script:
    """python $projectDir/scripts/plot.py $pbd . $name"""
}

process ontCallIndel {
    storeDir    "$params.directory/$params.analysis_name/$name/05.variants"
    conda       "$params.conda_path/clair3"
    tag         "$name"

    input:
        tuple   val(name), path(bam), path(bai)
    output:
        tuple   val(name), path("${name}.indel.raw.vcf.gz")
    script:
    """if [ ! -e $projectDir/genomes/$params.virus/sequences.fa.fai ]; then samtools faidx $projectDir/genomes/$params.virus/sequences.fa; fi && \
    run_clair3.sh  \
    --bam_fn $bam \
    --ref_fn $projectDir/genomes/$params.virus/sequences.fa \
    --model_path $projectDir/clair3_model/$params.model \
    -p ont \
    -o clair3_output \
    --include_all_ctgs \
    --chunk_size=5000 \
    --no_phasing_for_fa \
    --fast_mode \
    --min_coverage=$params.min_dp \
    --threads $params.map_cpu  
    cp ./clair3_output/merge_output.vcf.gz  ${name}.indel.raw.vcf.gz
    # bcftools filter -i "TYPE='indel' && FILTER='PASS' && FORMAT/DP >= $params.min_dp && (GT='1/1' || FORMAT/AF > 0.6)" ./clair3_output/merge_output.vcf.gz > ${name}.indel.vcf"""
}

process ontFiltIndel {
    storeDir    "$params.directory/$params.analysis_name/$name/05.variants"
    conda       "$params.conda_path/artic-like"
    tag         "$name"

    input:
        tuple   val(name), val(indel_raw_vcf)
    output:
        tuple   val(name), path("${name}.indel.vcf")
    script:
    // TODO: bcftools filter is based on experience, especially (GT='1/1' || FORMAT/AF > 0.6). So some false negative variants may occur.
    """bcftools filter -i "TYPE='indel' && FILTER='PASS' && FORMAT/DP >= $params.min_dp && (GT='1/1' || FORMAT/AF > 0.6)" ${indel_raw_vcf} > ${name}.indel.vcf"""
}

process getLowCoverRegion {
    storeDir    "$params.directory/$params.analysis_name/$name/04.depth"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "$name"

    input:
        tuple   val(name), path(pbd)
    output:
        tuple   val(name), path("${name}.low_cover.bed")
    script:
    """awk -v min_dp=$params.min_dp 'BEGIN{OFS="\\t"}{if(\$3<min_dp) {print \$1,\$2-1,\$2}}' $pbd > 1.bed && \
    bedtools merge -i 1.bed > ${name}.low_cover.bed"""
}

process igvtools {
    storeDir    "$params.directory/$params.analysis_name/$name/04.depth"
    conda       "$params.conda_path/artic-like"
    debug       true
    tag         "$name"

    input:
        tuple   val(name), path(bam), path(bai)
    output:
        tuple   val(name), path("${name}.per-pos.each.bases.depth")
    script:
    // the first 4 lines of result of igvtools are description, so using sed to delete this. And for each contig of reference, there always is a header formatted by
    // "variableStep chrom={congitName} span=N", so this line will be deleted and new column of contig name be added by awk
    """igvtools count $bam stdout $projectDir/genomes/$params.virus/sequences.fa -w 1 --bases > a.txt && \
    sed '1,4d' a.txt | awk 'BEGIN{OFS="\\t"}{if (\$0 ~ /^variableStep/) {split(\$0, a, " "); sub("chrom=", "", a[2]); contig=a[2]} else {print contig,\$0} }' > \
    ${name}.per-pos.each.bases.depth"""
}

process ontGetVariants {
    storeDir    "$params.directory/$params.analysis_name/$name/05.variants"
    conda       "$params.conda_path/artic-like"
    tag         "$name"

    input:
        tuple   val(name), path(low_cover), path(ppebd), path(indel)
    output:
        tuple   val(name), path("${name}.vcf.gz"), path("${name}.vcf.gz.tbi")
    script:
    // result of igvtools will be used to call SNV, then merge it and indel result of clair3
    """python $projectDir/scripts/get_variants.py --snv_freq $params.min_snv_freq --outfile ${name}.vcf $projectDir/genomes/$params.virus/sequences.fa $low_cover $ppebd $indel && \
    bgzip ${name}.vcf && tabix -p vcf ${name}.vcf.gz"""
}

process snpEff {
    storeDir    "$params.directory/$params.analysis_name/$name/05.variants"
    conda       "$params.conda_path/artic-like"
    tag         "$name"

    input:
        tuple   val(name), path(vcf), path(vcf_tbi)
    output:
        tuple   val(name), path("${name}.report.vcf")
    script:
    """snpEff -c $projectDir/config/snpEff.config -dataDir $projectDir/genomes -no-downstream -noStats -no-upstream $params.virus $vcf > ${name}.report.vcf"""
}

process getConsensus {
    storeDir    "$params.directory/$params.analysis_name/$name/06.consensus"
    conda       "$params.conda_path/artic-like"
    tag         "$name"
    
    input:
        tuple   val(name), path(vcf), path(vcf_tbi), path(low_cover)
    output:
        tuple   val(name), path("${name}.consensus.fasta")
    script:
    // to distinguish every nucleic acid sequence, the name of sequence must be updated by sed
    """bcftools consensus -f $projectDir/genomes/$params.virus/sequences.fa -m ${low_cover} $vcf > ${name}.consensus.fasta && \
    sed -i  '/^>/s/>/>${name}-/' ${name}.consensus.fasta && \
    echo && \
    echo =============================================================================================== && \
    echo "                         ${name} consensu finished                                           " && \
    echo =============================================================================================== """
}

process pangolin {
    storeDir    "$params.directory/$params.analysis_name/$name/07.pangolin"
    conda       "$params.conda_path/artic-like-pangolin"
    tag         "$name"

    input:
        tuple   val(name), path(consensus)
    output:
        tuple   val(name), path("${name}.pangolin.report")
    script:
    """if [[ $consensus == "empty.file" ]]; then touch ${name}.pangolin.report; else pangolin $consensus --outfile  ${name}.pangolin.report; fi"""
}

workflow ont {
    main:
        input = Channel.fromList(detectLongReads())
        ontPreprocess(input).set{ontPreprocessOut} // ontPreprocessOut: [name, raw.fastq]

        ontFiltlong(ontPreprocessOut).set{ontFiltlongOut} // ontFiltlongOut: [name, clean.fastq]

        ontMap(ontFiltlongOut).set{ontMapOut} // ontMapOut: [name, sorted.bam, sorted.bam.bai]

        trimPrimer(ontMapOut).set{trimPrimerOut} // trimPrimerOut: [name, trimmed.bam, trimmed.bam.bai, trimmed.report]

        trimPrimerOut.map{it -> it[0..2]}
            .join(ontPreprocessOut)
            .join(ontFiltlongOut)
            .set{ontStatIn} // ontStatIn: [name, trimmed.bam, timmed.bam.bai, raw.fastq, clean.fastq]

        ontStat(ontStatIn).set{ontStatOut} // ontStatOut: [name, raw.png, clean.png, mapped.png, summary.txt]

        getPerBaseDp(trimPrimerOut.map{it -> it[0..2]}).set{getPerBaseDpOut} // getPerBaseOut: [name, pbd]


        (ontCallIndel(trimPrimerOut.map{it -> it[0..2]}) | ontFiltIndel).set{ontFiltIndelOut} // ontFiltIndelOut: [name, indel.vcf]

        igvtools(trimPrimerOut.map{it -> it[0..2]}).set{igvtoolsOut} // igvtoolsOut: [name, ppebd]

        plot(getPerBaseDpOut).set{plotOut} // plotOut: [name, rawCov.png, modifiedCov.png]
        // TODO: if the reference is multi segment e.g. influenza, so the second and third element in plotOut will be list

        getLowCoverRegion(getPerBaseDpOut).set{getLowCoverRegionOut} // getLowCoverRegion: [name, low_cover.bed]

        getLowCoverRegionOut
            .join(igvtoolsOut)
            .join(ontFiltIndelOut)
            .set{ontGetVariantsIn}
        
        ontGetVariants(ontGetVariantsIn).set{ontGetVariantsOut} // ontGetVariantsOut: [name, vcf.gz, vcf.gz.tbi]

        ontGetVariantsOut
            .join(getLowCoverRegionOut)
            .set{getConsensusIn}
        
        getConsensus(getConsensusIn).set{getConsensusOut} // getConsensusOut: [name, consensus.fasta]

        snpEffbin = file("$projectDir/genomes/$params.virus/snpEffectPredictor.bin")

        snpEffOut = snpEffbin.exists() ? snpEff(ontGetVariantsOut) : ontGetVariantsOut.map{it -> it[0, 1]} // snpEffOut: [name, report.vcf] or [name, vcf.gz]

        if (params.pango_virus.split(",").contains(params.virus)) {
            pangolinOut = pangolin(getConsensusOut)
        } else {
            pangolinOut = getConsensusOut.map{it ->[it[0], file("empty.file")]} | pangolin
        }

    emit:
        summary = ontStatOut // [name, raw.png. clean.png, mapped.png, summary.txt]
        figures = plotOut // [name, rawCov.png, modifiedCov.png]
        vcfRpt = snpEffOut // [name, report.vcf]
        consensus = getConsensusOut // [name, consensus.fasta]
        pangoRpt = pangolinOut // [name, pangoRpt]

}