# Travis CI errors, when build logs gets too long. Therefore, suppress build
# logs for dependencies by default:
SILENT() {
    logfile="logs/$(basename ${1%.sh}).log"
    "$@" >$logfile || (cat "$logfile"; false)
}

mkdir -p logs
