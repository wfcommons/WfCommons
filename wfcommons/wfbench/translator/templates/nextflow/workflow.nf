params.simulate = false
params.pwd = null
params.help = null
pwd = null

def printUsage(error_msg, exit_code) {
    def usage_string = """
Usage: nextflow run workflow.nf --pwd /path/to/directory [--simulate] [--help]

    Required parameters:
      --pwd         Working directory (where the workflow.nf file is located)

    Optional parameters:
      --help        Show this message and exit.
      --simulate    Use a "sleep 1" for all tasks instead of the WfBench benchmark.
"""
    if (error_msg) {
        def RED = '\u001B[31m'
        def RESET = '\u001B[0m'
        System.err.println "${RED}Error: ${RESET}" + error_msg
    }
    System.err.println usage_string
    exit exit_code
}

def validateParams() {
    if (params.help) {
        printUsage(msg = "", exit_code=0)
    }
    if (params.pwd == null) {
        printUsage(msg = "Missing required parameter: --pwd", exit_code=1)
    }
    pwd = file(params.pwd).toAbsolutePath().toString()
    if (!file(pwd).exists()) {
        printUsage(msg = "Directory not found: ${pwd}", exit_code=1)
    } 
}

// Call validation at the start
validateParams()

# Generated code goes here