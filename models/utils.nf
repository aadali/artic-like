def printError(msg) {
    line = "=" * msg.length()
    log.error("\n\033[31m${line}\n${msg}\n${line}\033[0m")
    exit(2)
}

def q2p(qual) {
    return (1-10**(-qual/10))*100
}

def checkArgs(args) {
    badArgs= params.keySet() - args
    if (badArgs) {
        printError("Could not understand params: ${badArgs}")
    }

    if (params.help) {
        help()
        exit(0)
    }

    if (params.long_reads == null && params.reads1==null) {
        printError("one of input params: [long_reads, reads1] must be set")
    }

    if (params.long_reads != null && params.reads1 != null) {
        printError("--long_reads and --reads1 couldn't be set simultaneously")
    }

    for (String read in [params.long_reads, params.reads1, params.reads2]) {
        if (read != null) {
            if (!read.startsWith("/")) { printError("input reads  must be a absolute path") }
            if (!file(read).exists()) { printError("No such file for: ${read}") }
        }
    }

    if ( !(params.conda_path != null && params.conda_path.startsWith("/") && file(params.conda_path).exists()) ) {
        printError("conda_path msut be a absolute path and exists")
    }
    
    if (!params.directory.startsWith("/")) { printError("--directory must be a absolute path") }

    allowedVirus = file("$projectDir/genomes").list()
    if (!allowedVirus.contains(params.virus)) {
        printError("--virus must be one of ${allowedVirus}")
    }

    clair3_models = file("$projectDir/clair3_model").list()
    if (!clair3_models.contains(params.model)) {
        printError("--model must be one of ${clair3_models}")
    }
}

def getExtension(input_fastq) {
    lens = input_fastq.size()
    gz_lens = input_fastq.findAll{it -> it.endsWith(".fastq.gz") || it.endsWith(".fq.gz")}.size()
    fq_lens = input_fastq.findAll{it -> it.endsWith(".fastq") || it.endsWith(".fq")}.size()
    if (lens == gz_lens) return "fastq.gz"
    if (lens == fq_lens) return "fastq"
    printError("compressed and uncompressed files were found simultaneously")
}

def detectLongReads() {
    /*
    return ArrayList<LinkedHashMap<String, String>>, each element is a sample
        name: the sample_name of this sample
        fastqs: the fastq[.gz] file path
        cmd: how handle this fastqs, if fastqs is a file, make symbolic link; if it's files, cat them into new file
        direction: '' when fastqs is a file or '>' when fastqs is files path
        ext: extension of the symbolic link name and merged file depends on the file's extension of the --long_reads
    */
    input_fastq = file(params.long_reads, checkIfExists: true)
    if (input_fastq.isEmpty()) {printError("--long_reads is empty")}
    if (input_fastq.isFile()) {
        // --long_reads is a fastq[.gz]
        ext = getExtension([params.long_reads])
        return [[name:params.analysis_name, fastqs:params.long_reads, cmd: 'ln -sf ', direction: '', ext: ext]]
    } 
    if (input_fastq.isDirectory()) {
        // --long_reads is a directory, maybe a barcode directory or fastq_pass directory
        contents = input_fastq.list()
        contents_len = contents.size()
        subdirs_len = contents.findAll{it -> file("${params.long_reads}/$it").isDirectory()}.size()
        subfiles_len = contents.findAll{it -> file("${params.long_reads}/$it").isFile()}.size()

        if (contents_len != subdirs_len && contents_len != subfiles_len) {
            printError("file and directory found in ${params.long_reads} simultaneously")
        }

        if (contents_len == subfiles_len) {
            // --long_reads is a barcode directory
            ext = getExtension(contents.collect{ it -> "${params.long_reads}/$it"})
            a = [[name: params.analysis_name, fastqs: "${params.long_reads}/*", cmd: 'cat', direction: ' > ', ext: ext]]
            return a
        } 

        if (contents_len == subdirs_len) {
            // --long_reads is a fastq_pass directory
            a = []
            for (def each_name in contents) {
                fqs = file("${params.long_reads}/${each_name}").list()
                if (fqs.size() <= params.files_per_bar || each_name == "unclassified") continue // files number in subdirectory or name == "unclassified" will be ignored
                ext = getExtension(fqs)
                a << [name: each_name, fastqs: "${params.long_reads}/${each_name}/*", cmd: "cat", direction: " > ", ext: ext]
            }
            return a
        }
    }
}

def help() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    c_bold = "\033[1m"
    msg = """
    ________________________________________________________________________________________________________________________

    ${c_green}This is a pipeline for analysising virus amplicon reads of ont or ngs to generate variants and consensus.${c_reset}
    ________________________________________________________________________________________________________________________

    ${c_yellow}Usage example: ${c_reset}
    nextflow run main.nf --long_reads /the/path/to/nanopore/reads --directory /the/output/directory --analysis_name testName 
    nextflow run main.nf --reads1 /the/path/to/ngs/reads1 --reads2 /the/path/to/ngs/reads2 --analysis_name testNgs
    nextflow run main.nf --reads1 /the/path/to/ngs/single/reads --analsysi_name testSingleNgs

    ${c_yellow}Input:${c_reset}
    --long_reads    ${c_dim}[string]${c_reset}  the path to nanopore reads, there are three situations. ${c_dim}[default: null]${c_reset}
                                  1. a fastq[.gz] file
                                  2. a directory containing some fastq[.gz] files e.g. barcode directory of fastq_pass
                                  3. a directory containing some subdirectories that whose contents is fastq[.gz] e.g. fastq_pass of nanopore data.
                                     In ths case, multi samples will be analysised parallelly.
    --reads1        ${c_dim}[string]${c_reset}  the path to ngs reads1. one of {reads1, long_reads1} must be specified. ${c_dim}[default: null]${c_reset}
    --reads2        ${c_dim}[string]${c_reset}  the path to ngs reads2. ${c_dim}[default: null]${c_reset}

    ${c_yellow}Nanopore options:${c_reset}
    --min_read_len  ${c_dim}[integer]${c_reset} the min read length of nanopore reads. ${c_dim}[default: 300]${c_reset}
    --max_read_len  ${c_dim}[integer]${c_reset} the max read length of nanopore reads. ${c_dim}[default: 1900]${c_reset}
    --min_read_qual ${c_dim}[integer]${c_reset} the min read quality of nanopore reads. ${c_dim}[default: 9]${c_reset}
    --trim_bases    ${c_dim}[integer]${c_reset} the number of bases trimmed from alignment after mapping in order to reduce the impact of primer if there is no --primer_bed specified for nanopore reads
                              For ngs reads, this value is invalid. ${c_dim}[default: 300]${c_reset}
    --min_overlap   ${c_dim}[integer]${c_reset} if the alignment length is less than this value, this alignment will be dropped, ${c_dim}[default: 300]${c_reset} for nanopore reads and ${c_dim}[fixed: 20]${c_reset} for ngs reads
    --files_per_bar ${c_dim}[integer]${c_reset} when long_reads is situation 3, if the files number in some subdirectories is less than this value, these subdirectories will not be analysised. ${c_dim}[default: 2]${c_reset}
                              This pipeline will no analysis the reads in unclassified in fastq_pass 
    --model         ${c_dim}[string]${c_reset}  the model used to call variants for clair3, must be one of directory names in projectDir/clair3_model. ${c_dim}[default: r941_prom_sup_g501]${c_reset}

    ${c_yellow}Common options:${c_reset}
    --map_cpu       ${c_dim}[integer]${c_reset} how many threads will be used when mapping. ${c_dim}[default: 4]${c_reset}
    --min_mapq      ${c_dim}[integer]${c_reset} the min map quality. ${c_dim}[default: 60]${c_reset}
    --min_dp        ${c_dim}[integer]${c_reset} the min depth, if the depth of a site < this value, this site will be N in consensus fasta. ${c_dim}[default: 20]${c_reset}
    --min_snv_freq  ${c_dim}[float]${c_reset}   if snv freq is less than this value, this position will be reference base in consensus fasta. ${c_dim}[default: 0.55]${c_reset}
    --primer_bed    ${c_dim}[string]${c_reset}  a file of infomation of primer, pipeline will trim primer depends on this file. ${c_dim}[default: null]${c_reset}
    --virus         ${c_dim}[string]${c_reset}  what virus is this sample, must be one of directory names in projectDir/genomes. ${c_dim}[default: sars-cov-2]${c_reset}
    --conda_path    ${c_dim}[string]${c_reset}  where is the path of conda environment, user should set this value in command line or directly update this value in nextflow.config. ${c_dim}[default: null]${c_reset}
    --directory     ${c_dim}[string]${c_reset}  the path to store the results. ${c_dim}[default: $projectDir]${c_reset}
    --analysis_name ${c_dim}[string]${c_reset}  the results will be stored in the --directory, named --analysis_name. ${c_dim}[default: test001]${c_reset}
    --pango_virus   ${c_dim}[string]${c_reset}  if --virus in this virus list separated by comma, pangolin will be used to identify typing
    --help          ${c_dim}[bool]${c_reset}    show this help and exits
    """
    println(msg.stripIndent())
}