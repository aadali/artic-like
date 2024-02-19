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

    if ( !(params.conda_path.startsWith("/") && (file(params.conda_path).exists())) ) {
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
    input_fastq = file(params.long_reads, checkIfExists: true)
    if (input_fastq.isEmpty()) {printError("--long_reads is empty")}
    if (input_fastq.isFile()) {
        ext = getExtension([params.long_reads])
        return [[name:params.analysis_name, fastqs:params.long_reads, cmd: 'ln -sf ', direction: '', ext: ext]]
    } 
    if (input_fastq.isDirectory()) {
        contents = input_fastq.list()
        contents_len = contents.size()
        subdirs_len = contents.findAll{it -> file("${params.long_reads}/$it").isDirectory()}.size()
        subfiles_len = contents.findAll{it -> file("${params.long_reads}/$it").isFile()}.size()

        if (contents_len != subdirs_len && contents_len != subfiles_len) {
            printError("file and directory found in ${params.long_reads} simultaneously")
        }

        if (contents_len == subfiles_len) {
            ext = getExtension(contents.collect{ it -> "${params.long_reads}/$it"})
            a = [[name: params.analysis_name, fastqs: "${params.long_reads}/*", cmd: 'cat', direction: ' > ', ext: ext]]
            return a
        } 

        if (contents_len == subdirs_len) {
            a = []
            for (def each_name in contents) {
                fqs = file("${params.long_reads}/${each_name}").list()
                if (fqs.size() <= params.files_per_bar || each_name == "unclassified") continue
                ext = getExtension(fqs)
                a << [name: each_name, fastqs: "${params.long_reads}/${each_name}/*", cmd: "cat", direction: " > ", ext: ext]
            }
            return a
        }
    }
}