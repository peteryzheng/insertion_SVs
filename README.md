# Identifying Statistically Significant Templates for SV Junctional Insertions

## Introduction

Briefly, the insertion template significance is measured as follows:

1.  We have a set of SVs with insertion sequences (size N). Each SV breakpoint has four possible directions where the adjacent sequences could have given rise to the insertion sequences.
    
2. For each adjacent directional sequence, 
    
    1.  For each insertion sequence (length l), 
    
        1.  We will use the other SVs with insertion sequences as the background.

        2.  We will obtain alignment scores between insertion sequences and breakpoint adjacent sequences of varying lengths (l - 20 x l) to get a N x 19l matrix of alignment scores.

        3.  We than calculate column wise quantile scores which represents how well the insertion sequences align to the real vs background breakpoint adjacent sequences at each offset position away from the breakpoint.

        4.  We will then get the maximum quantile score for the insertion sequence at every row (SV) which represents how well the insertion sequence aligns to the real vs background breakpoint adjacent sequences overall (at the best position away from the breakpoint). This gives us a N x 1 vector of maximum quantile scores.

        5.  Finally, we do a quantile conversion of the maximum quantile scores to get the final offset-adjusted p-values for the insertion sequences, specifically for the real breakpoint adjacent sequences of the SV where the insertion sequence came from.

    2.  FDR correction is done on the final offset-adjusted p-values to get the final insertion template significance.


Our method for detecting structural variants with statistically significant templates is based on the following assumptions:

1.  Some insertions originate from sequencesclose to the breakpoints

2.  The insertion sequences are not common in breakpoint adjacent sequences of other SVs with insertion sequences


## Usage

The code base is mainly divided into three steps:

1.  Data preprocessing (processing various SV caller output VCFs into a common format)

2.  Template significance calculation (calculating the significance of insertion templates)

3.  Configuration calling (calling the configuration of the SVs based on the template significance)

Each step is run as separate scripts. 

### Data Preprocessing

#### Scripts and Commands

The data preprocessing step is done using the `insertion_vcf_processing.R` script. This script calls helper functions in `helper_functions.R` to process the VCFs from various SV callers into a common format. The script takes in the following arguments:

```
Options:
	-d DATASET, --dataset=DATASET
		Dataset to use. Current option: [HCMI, 1KG, hcmi_cancer_models]
	-c CALLER, --caller=CALLER
		Caller to use. Current option: [manta, svaba]
	-p DEDUP, --dedup=DEDUP
		Whether to dedup the SVs.
	-o OUTPUTDIR, --outputdir=OUTPUTDIR
		Output directory to use.
```


* To run processing on HCMI datasets using Manta, you can run this on **DIPG server** with the following command:

    ```
    conda run -n rameen \
        Rscript /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/insertion_vcf_processing.R  \
        -d hcmi_cancer_models -c manta -p FALSE
    ```

* To run larger dataset especially germline datasets, you probably need to run this on a cluster with more memory. You can run this on **UGER** with the following command:

    ```
    qsub  -l h_vmem=256G,h_rt=08:00:00 -o ./ -e ./ \
        /xchip/beroukhimlab/youyun/miniconda3/condabin/conda run -n rameen \
        Rscript /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/insertion_vcf_processing.R  \
        -d 1KG -c manta -p TRUE
    ```

#### Workflow Information

By default, this generates a folder in the `outputdir` called **ins_align_total_[timestamp]**. The path of this subfolder will be the *input to the next step*. This folder will contain the following files:
1. **total_SVs_processed_[timestamp].tsv**: This file contains all SVs in the cohort in the processed format
2. **insertions_SVs_processed_[timestamp].tsv**: This file contains all SVs with insertions in the cohort in the processed format with some added columns (e.g. adjacent sequences) that we will use in the analysis. 
3. Various visualization plots (e.g. insertion length distribution, SV burden, etc.)



### Template Significance Calculation

#### Scripts and Commands

The template significance calculation step is done using the `batch_align_kmers.R` script. This script writes a qsub script into the output directory which the user can use on UGER to submit a task array for each insertion sequence. The script takes in the following arguments:

```
Options:
	-a ALIGNPARAM, --alignparam=ALIGNPARAM
		Alignment Parameters to Use. Current option: [mhe, default_bwa]
	-o OUTPUTDIR, --outputdir=OUTPUTDIR
		Output directory to use.
```

To generate the qsub script, you can run the following command and follow the instructions:

```
conda run -n rameen \
    Rscript /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/batch_align_kmers.R \
    -a mhe -o /xchip/beroukhimlab/youyun/nti/output/ins_align_total_[timestamp]
```

Not all the tasks will finish the first time around, so there is another script designed to resubmit the failed tasks. This script is called `batch_align_kmers_rerun.r`. This script takes in the following arguments:

```
Options:
	-a ALIGNPARAM, --alignparam=ALIGNPARAM
		Alignment Parameters to Use. Current option: [mhe, default_bwa]
	-o OUTPUTDIR, --outputdir=OUTPUTDIR
		Output directory to use.
	-t TASK_NUMBER, --task_number=TASK_NUMBER
		Task array number to check.
	-k KMER_FILE, --kmer_file=KMER_FILE
		Latest kmer file to check.
```

To rerun the failed tasks, you can run the following command and follow the instructions:

```
conda run -n rameen \
    Rscript /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/batch_align_kmers_rerun.R \
    -a mhe -o /xchip/beroukhimlab/youyun/nti/output/ins_align_total_[timestamp] \
    -t [job ID] -k /xchip/beroukhimlab/youyun/nti/output/ins_align_total_[timestamp]/inputs/[kmer_file]
```

Both of the above scripts will use helper functions in `align_and_config_call_helper.R` to generate the task array shell scripts. The task array shell scripts will use `align_nearby_utils.R` to align the insertion sequences to the breakpoint adjacent sequences.

#### Workflow Information

This scripts create subfolder structure within the **ins_align_[timestamp]**. The structure of this is as following (--| denotes a folder, -- denotes a file):

```
| ins_align_[timestamp]
--| inputs
    -- kmer file — input for the initial run of task array
    -- failed kmer file - input for reruns of the failed tasks
    -- task array shell scripts — shell scripts for the task arrays
--| outputs
    --| alignment_files — alignment score outpus
    --| task_array_output — task array outputs
    --| tmp_working_directory — core dump directory
    -- [job ID].tsv - task array resource usage tables
    -- runtime_plot_[job ID].pdf -- task array runtime by kmer length plots
-- total_SVs_processed_[timestamp].tsv
-- insertions_SVs_processed_[timestamp].tsv
-- various visualization plots
```

### Configuration Calling

#### Scripts and Commands

The configuration calling step is done using the `configuration_calling.R` script. This script takes in the following arguments:

```
Options:
	-o OUTPUTDIR, --outputdir=OUTPUTDIR
		Output directory to use.
	-n DATASETNAME, --datasetname=DATASETNAME
		Dataset to use. Current option: [HCMI, 1KG, hcmi_cancer_models]
```

To run the configuration calling, you can run the following command:

```
conda run -n rameen \
    Rscript /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/configuration_calling.R \
    -o /xchip/beroukhimlab/youyun/nti/analysis_files/insertions/ins_align_total_[timestamp] \
    -n hcmi_cancer_models
```

#### Workflow Information

This script generates output files in the **ins_align_total_[timestamp]** folder. We will generate several files

1. Single breakend configuration
    
    * **[datasetname]_insertions_w_sig_alignment_[timestamp].tsv** -- this file contains information on all configurations at a single breakend level
    * **[datasetname]_alignment_upset_plot_[timestamp].pdf** -- this file contains an upset plot of the single breakend configurations

2. Double breakend configuration

    * **[datasetname]_alignment_combination_[timestamp].tsv** -- this file contains information on all configurations at a breakend pair level

## Final Output

Thus the final output structure of the workflow will be as follows (--| denotes a folder, -- denotes a file):

```
| ins_align_[timestamp]
--| outputs
    -- [datasetname]_insertions_w_sig_alignment_[timestamp].tsv
    -- [datasetname]_alignment_upset_plot_[timestamp].pdf
    -- [datasetname]_alignment_combination_[timestamp].tsv
    --| alignment_files — alignment score outputs
        --| [kmer 1]
        --| [kmer 2]
        --| ...
    --| task_array_output — task array outputs
        -- [job ID]:1.error/output
        -- [job ID]:2.error/output
        -- ...
    --| tmp_working_directory — core dump directory
    -- [job ID].tsv - task array resource usage tables
    -- runtime_plot_[job ID].pdf -- task array runtime by kmer length plots
--| inputs
    -- kmer file — input for the initial run of task array
    -- failed kmer file - input for reruns of the failed tasks
    -- task array shell scripts — shell scripts for the task arrays
-- total_SVs_processed_[timestamp].tsv
-- insertions_SVs_processed_[timestamp].tsv
-- various visualization plots
```