echo true

params.bucket = "s3://pwyming-tmp-us-east-1"
meats = file("${params.bucket}/meats.txt")

process read_meats {
    // directives
    // a container image is required
    container "ubuntu:latest"

    // compute resources for the Batch Job
    cpus 1
    memory '512 MB'

    script:
    """
    echo "${meats.text}"
    """
}