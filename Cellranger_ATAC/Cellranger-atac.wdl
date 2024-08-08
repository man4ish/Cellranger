# Task to process raw sequencing data using Cell Ranger ATAC
task count {
    input {
        # Optional parameters
        Boolean? nosecondary     # Flag to indicate whether secondary analysis is to be skipped
        Int? localcores          # Number of local cores to use
        String? localmem         # Local memory to allocate
        String? localvmem        # Local virtual memory to allocate
        
        # Required parameters
        String sid               # Sample ID for Cell Ranger ATAC
        File fastq_tar           # Tarball of gzipped FASTQ files
        File ref_tar             # Tarball of gzipped reference genome
        Int memory_count         # Memory to allocate for this task
        Int cpu_count            # Number of CPUs to allocate for this task
        Int storage_count        # Disk space required for this task
        
        # Derived parameter
        String base = basename(fastq_tar, ".tar.gz")  # Base name for result files
    }

    parameter_meta {
        sid: "id"
        fastq_tar: "fastqs gzipped file"
        ref_tar: "gzipped reference genome"
        nosecondary: "no secondary"
        localcores: "number of local cores"
        localmem: "local memory"
        localvmem: "local virtual memory"
    }
    
    command {
        # Exit on error and fail on any command in the pipeline that fails
        set -exo pipefail
        
        # Unpack FASTQ files
        mkdir fastq_bundle
        echo ~{fastq_tar}
        tar -xzf ~{fastq_tar} -C fastq_bundle --no-same-owner
        fastq_path=$(realpath fastq_bundle)
        echo $fastq_path
        mv fastq_bundle/*/* fastq_bundle
        fastq_file_path="$(dirname $fastq_path)"
        echo $fastq_file_path

        # Unpack reference genome
        mkdir ref_bundle
        echo ~{ref_tar}
        tar -xzf ~{ref_tar} -C ref_bundle --no-same-owner
        ref_path=$(realpath ref_bundle)
        echo $ref_path
        mv ref_bundle/*/* ref_bundle
        ref_file_path="$(dirname $ref_path)"
        echo $ref_file_path

        # Construct Cell Ranger ATAC command
        cmd="atac --id=${sid} --fastqs=$fastq_path --reference=$ref_path"
        
        # Add optional parameters to the command if specified
        if ${nosecondary}; then cmd+=" --nosecondary"; fi
        if [ -n "${localcores}" ]; then cmd+=" --localcores ${localcores}"; fi
        if [ -n "${localmem}" ]; then cmd+=" --localmem ${localmem}"; fi
        if [ -n "${localvmem}" ]; then cmd+=" --localvmem ${localvmem}"; fi
        echo $cmd

        # Run the Cell Ranger ATAC command
        /software/reboot-utils/cellranger/bin/cellranger $cmd
        
        # Package the results
        tar -czvf ${base}_result.tar.gz "$ref_file_path/${sid}"
    }
    
    output { 
        # Output result tarball
        File out = '${base}_result.tar.gz'
    }
    
    runtime {
        # Resource requirements
        memory: memory_count + "G"
        cpu: cpu_count
        disk: storage_count + "GB"
        docker: "docker.io/man4ish/cellranger:latest"
    }
}

# Task to aggregate results from multiple Cell Ranger ATAC outputs
task aggr {
    input {
        # Required parameters
        String aggr_id           # Aggregation ID for Cell Ranger ATAC
        File data_tar            # Tarball of archive files for ATAC-seq data
        File aggr_csv_file       # CSV file listing samples and their paths
        String? norm_type        # Optional normalization type for aggregation
        
        # Resource requirements
        Int memory_aggr          # Memory to allocate for this task
        Int cpu_aggr             # Number of CPUs to allocate for this task
        Int storage_aggr         # Disk space required for this task
        
        # Derived parameter
        String base = basename(aggr_csv_file, ".csv")  # Base name for result files
    }

    parameter_meta {
        aggr_id: "aggr id"
        data_tar: "archive files for ATAC-seq data"
        aggr_csv_file: "csv file"
        norm_type: "norm type"
    }

    command {
        # Exit on error and fail on any command in the pipeline that fails
        set -exo pipefail
        
        # Unpack data files
        mkdir data_bundle
        echo ~{data_tar}
        tar -xzf ~{data_tar} -C data_bundle --no-same-owner
        data_path=$(realpath data_bundle)
        echo $data_path
        mv data_bundle/*/* data_bundle
        data_file_path="$(dirname $data_path)"
        echo $data_file_path
        
        # Process the CSV file to update paths
        while IFS= read -r line; do
            if [[ "$line" == "sample_id"* ]]; then
                echo "$line" >> ${aggr_csv_file}.tmp
            else
                A="$(cut -d',' -f1 <<<"$line")"
                echo "$A,$data_path/$A/outs/molecule_info.h5" >> ${aggr_csv_file}.tmp
            fi
        done < ${aggr_csv_file}
        
        # Replace original CSV with updated CSV
        mv ${aggr_csv_file}.tmp ${aggr_csv_file}

        # Construct Cell Ranger ATAC aggregation command
        ref_path=$(realpath ${aggr_csv_file})
        ref_file_path="$(dirname $ref_path)"
        cmd="aggr --id=${aggr_id} --csv=${aggr_csv_file}"
        
        # Add optional normalization type to the command if specified
        if [ -n "${norm_type}" ]; then cmd+=" --norm_type ${norm_type}"; fi
        
        # Run the Cell Ranger ATAC aggregation command
        /software/reboot-utils/cellranger/bin/cellranger $cmd
        
        # Package the results
        tar -czvf ${base}_result.tar.gz "$data_file_path/${aggr_id}"
    }

    output {
        # Output result tarball
        File out = '${base}_result.tar.gz'
    }

    runtime {
        # Resource requirements
        memory: memory_aggr + "G"
        cpu: cpu_aggr
        disk: storage_aggr + "GB"
        docker: "docker.io/man4ish/cellranger:latest"
    }
}

# Workflow to execute the Cell Ranger ATAC count and aggregation tasks
workflow CellrangerATAC {
    input {
        # Inputs for both tasks
        Boolean? nosecondary
        Int? localcores
        String? localmem
        String? localvmem
        String sid
        File fastq_tar
        File ref_tar
        String aggr_id
        File data_tar
        File aggr_csv_file
        String? norm_type
        Int memory_count
        Int cpu_count
        Int storage_count
        Int memory_aggr
        Int cpu_aggr
        Int storage_aggr
    }

    call count {
        input: 
            nosecondary=nosecondary, 
            localcores=localcores, 
            localmem=localmem,  
            localvmem=localvmem, 
            sid=sid, 
            fastq_tar=fastq_tar, 
            ref_tar=ref_tar,
            memory_count=memory_count,
            cpu_count=cpu_count,
            storage_count=storage_count
    }
    
    call aggr {
        input: 
            aggr_id=aggr_id, 
            data_tar=data_tar, 
            aggr_csv_file=aggr_csv_file, 
            norm_type=norm_type,
            memory_aggr=memory_aggr,
            cpu_aggr=cpu_aggr,
            storage_aggr=storage_aggr
    }
}

