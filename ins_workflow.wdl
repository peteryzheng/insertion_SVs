version 1.0

workflow ins_workflow {
    input {
        String kmer
        String width = '20'
        String gap_open
        String gap_epen
        String mismatch_pen
        String match_pen
        String downsample_num
        String seed
        File input_file

        Int threads = 4
        Int mem_gb = 12
        Int preemptible_tries = 3
        Int output_disk_size = 100
        Int boot_disk_size = 20
    }

    call run_align_nearby_util {
        input:
            kmer = kmer,
            width = width,
            gap_open = gap_open,
            gap_epen = gap_epen,
            mismatch_pen = mismatch_pen,
            match_pen = match_pen,
            downsample_num = downsample_num,
            seed = seed,
            input_file = input_file,
            threads = threads,
            mem_gb = mem_gb,
            preemptible_tries = preemptible_tries,
            output_disk_size = output_disk_size,
            boot_disk_size = boot_disk_size
    }

    output {
        File kmer_zip = run_align_nearby_util.kmer_zip
    }
}

task run_align_nearby_util {
    input {
        String kmer
        String width
        String gap_open
        String gap_epen
        String mismatch_pen
        String match_pen
        String downsample_num
        String seed
        File input_file

        Int threads
        Int mem_gb
        Int preemptible_tries
        Int output_disk_size
        Int boot_disk_size
    }

    command <<<
        conda run -n insertions --live-stream \
            Rscript /insertion_SVs/align_nearby_utils.R \
            -i ~{kmer} -o . \
            -w ~{width} -d ~{input_file} \
            -g ~{gap_open} -e ~{gap_epen} \
            -m ~{mismatch_pen} -t ~{match_pen} \
            -n ~{downsample_num} -s ~{seed} \
            -c ~{threads}

        tar cvzf ~{kmer}.tar.gz ./~{kmer}
    >>>

    runtime {
        docker: 'zhengy1/ins_workflow:latest'
        memory: mem_gb + "GB"
        cpu: threads
        disks: "local-disk " + output_disk_size + " HDD"
        preemptible: preemptible_tries
        bootDiskSizeGb: boot_disk_size
    }

    output {
        File kmer_zip = "~{kmer}.tar.gz"
    }

}

