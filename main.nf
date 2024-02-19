include { 
        printError; 
        checkArgs;
        detectLongReads
        } from './models/utils.nf'
include {ont} from './models/ont.nf'
include {ngs} from './models/ngs.nf'

allowedArgs = [
    'long_reads',
    'reads1',
    'reads2',
    'min_read_len',
    'max_read_len',
    'min_read_qual',
    'map_cpu',
    'min_mapq',
    'min_dp',
    'min_snv_freq',
    'trim_bases',
    'min_overlap',
    'primer_bed',
    'files_per_bar',
    'virus',
    'model',
    'conda_path',
    'analysis_name',
    'directory',
    'pango_virus',
    'help'
]

checkArgs(allowedArgs)
// a = detectLongReads()

process makeExcelReport {
    storeDir    "$params.directory/$params.analysis_name/$params.analysis_name/08.report"
    conda       "$params.conda_path/artic-like"
    tag         "${name}"

    input:
        tuple   val(name), path(summary)
    output:
        tuple   val(name), path("${name}.results_file.txt"), path("${name}.report.xlsx")
    script:
    """cp ${summary} ${name}.results_file.txt &&  python $projectDir/scripts/make_excel_report.py ${name}.results_file.txt ${name}.report.xlsx"""
}

workflow {
    if (params.long_reads == null) {
        results = ngs()
        results.summary
            .join(results.figures)
            .join(results.vcfRpt)
            .join(results.consensus)
            .join(results.pangoRpt)
            .map {
                it -> {
                    a = [['summary', 'rawCovPng', 'modifiedCovPng', 'vcfRpt', 'consensus', 'pangoRpt', "sample_name"], it[1..<it.size()] + it[0]].transpose()
                    return [it[0], a]
                }
            }.set {makeExcelReportTmp}
    } else {
        results = ont()
        results.summary
            .join(results.figures)
            .join(results.vcfRpt)
            .join(results.consensus)
            .join(results.pangoRpt)
            .map {
                it -> {
                    content = []
                    a = [
                        ['rawPng', 'cleanPng', 'mappedPng', 'summary', 'rawCovPng', 'modifiedCovPng', 'vcfRpt', 'consensus', 'pangoRpt', 'sample_name'],
                        it[1..<it.size()] + it[0]
                    ].transpose()
                    return [it[0], a]
                }
            }.set {makeExcelReportTmp}
    }
    makeExcelReportTmp.map {
        it -> {
            content = []
            it[1].each {
                k,v -> {
                    line = k == "sample_name" ? "params\t$k\t$v" : "results\t$k\t$v"
                    content << line
                }
            }
            params.each {
                k, v -> {
                    if (['files_per_bar', 'trim_bases', 'min_overlap', 'model', 'min_read_len', 'max_read_len', 'min_read_qual'].contains(k)){
                        line = "params\t#$k\t$v"
                    } else {
                        line = "params\t$k\t$v"
                    }
                    content << line
                }
            }
            results_summary_file = file("/tmp/$params.analysis_name-${it[0]}.txt")
            results_summary_file.text = content.join("\n") + '\n'
            return [it[0], results_summary_file]
        }
    }.set { makeExcelReportIn }
    makeExcelReport(makeExcelReportIn)
}