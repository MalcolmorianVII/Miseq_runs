#!/bin/bash

path=$1
batch_abs_path=dirname "$(dirname "$(dirname "$path")")"
batch=basename $batch_abs_path

echo "Starting to run the analyze_miseq_data.sh"
echo "Are you sure you have backed up the data for this run:$batch.If not cancel using Ctrl + C and back the data to NAS"
sleep 10

echo "Proceeding with running analyze_miseq_data.sh"
echo "Getting the batch directory from the given absolute path"

echo "Batch is $batch"
echo "Confirm this? (y/n)"
read -r reply

if [[ $reply == "y" ]]; then
    echo "Batch confirmed"
else
    echo "Incorrect Batch name. Provide the script with proper path"
    exit 1
fi

echo Renaming Fastq to mlw_processing

mv $path ${batch_abs_path}/mlw_processing 

# Verify mlw_processing exists

if [ -d "${batch_abs_path}/mlw_processing" ] && [ ! -d $path ]; then
    echo "Successfully renamed Fastq to mlw_processing"
else
    echo "Directory 'mlw_processing' does not exist or 'Fastq' exists."
    exit 1
fi

# echo "Creating a fofn of samples"
# find "${batch_abs_path}/mlw_processing" -type f -name '*.fastq.gz' ! -name 'FastqSummaryF1L1.txt' ! -name 'NTC_*' ! -name 'Undetermined_*' | sed -e 's#.*/##' -e 's/_.*//' | sort -u > "${batch_abs_path}/mlw_processing/fofn.txt"

# if [ -s "${batch_abs_path}/mlw_processing/fofn.txt" ]; then
#     echo "Created a FOFN of samples!!!"
# else
#     echo "The FOFN is empty or does not exist."
#     exit 1
# fi

# echo "Ordering the fastqs into their respective sample directories"

# while IFS= read -r sample
#     do
#         mkdir -p "$sample"
#         mv "${batch_abs_path}/mlw_processing/${sample}"*.fastq.gz "${batch_abs_path}/mlw_processing/${sample}"
#         # Check if this successful per sample
#     done < "${batch_abs_path}/mlw_processing/fofn.txt"

# Run bactopia i.e ensure bactopia conda env is properly activated!!!

# prepare sample sheet

# bactopia prepare -f '_001.fastq.gz' "${batch_abs_path}/mlw_processing -r  > "${batch_abs_path}/mlw_processing/samples.txt 

# # Validation i.e samples.txt file exits and contain data
# if [ -s "${batch_abs_path}/mlw_processing/samples.txt " ]; then
#     echo "Created a samples.txt!!!"
# else
#     echo "The samples.txt is empty or does not exist."
#     exit 1
# fi

# # Run bactopia ...preferably in a screen session .... no run this script inside a screen session
# echo "Running bactopia using datasets for species Salmonella enterica"
# # echo "If different species please provide for it and run the script again"
# sleep 5
# bactopia --samples ${batch_abs_path}/mlw_processing/samples.txt --datasets /home/bkutambe/data/bactopia/datasets --species "Salmonella enterica" --outdir  ${batch_abs_path}/mlw_processing/${batch}_bactopia --max_cpus 100

# then proceed with processing
renameFastq() {
    sourceDir=$1
    targetDir=$2
    
    mv "${sourceDir}/${BATCH}" "${targetDir}"
}

generatefofn() {
    directory=$1
    echo "Creating a fofn of samples"
    find $directory -type f -name '*.fastq.gz' ! -name 'FastqSummaryF1L1.txt' ! -name 'NTC_*' ! -name 'BLANK_*' ! -name 'Undetermined_*' | sed -e 's#.*/##' -e 's/_.*//' | sort -u > "${directory}/fofn.txt"
    
    if [ -s "$directory/fofn.txt" ]; then
        echo "Created a FOFN of samples!!!"
    else
        echo "The FOFN is empty or does not exist."
        exit 1
    fi
    
}

order_fastqs() {
    directory=$1
    echo "Ordering the fastqs into their respective sample directories"

    while IFS= read -r sample
        do
            mkdir -p "${directory}/${sample}"
            mv "${directory}/${sample}"*.fastq.gz "${directory}/${sample}"
            # Check if this successful per sample
        done < "${directory}/fofn.txt"

}

prepare_samples() {
    directory=$1
    bactopia prepare -f '_001.fastq.gz' ${directory} -r  > "${directory}/samples.txt" 

    # Validation i.e samples.txt file exits and contain data
    if [ -s "${directory}/samples.txt"  ]; then
        echo "Created a samples.txt!!!"
    else
        echo "The samples.txt is empty or does not exist."
        exit 1
    fi
}

run_bactopia() {
    directory=$1
    # Run bactopia ...preferably in a screen session .... no run this script inside a screen session
    echo "Running bactopia using datasets for species Salmonella enterica"
    # echo "If different species please provide for it and run the script again"
    sleep 5
    bactopia --samples "${directory}/samples.txt" --datasets /home/bkutambe/data/bactopia/datasets --species "Salmonella enterica" --outdir  ${directory}/${batch}_bactopia --max_cpus 100
}


renameFastq
generatefofn ${batch_abs_path}/mlw_processing
order_fastqs ${batch_abs_path}/mlw_processing
# Activate bactopia env 

source ~/miniconda3/etc/profile.d/conda.sh

conda activate bactopia-dev 

prepare_samples ${batch_abs_path}/mlw_processing
run_bactopia ${batch_abs_path}/mlw_processing
# Deal with ownership on the miseq data from root to user
# unwanted = [
#     FastqSummaryF1L1.txt,
#     NTC_S51_L001_R1_001.fastq.gz,
#     Undetermined_S0_L001_R1_001.fastq.gz
# ]

find . -type f -name '*.fastq.gz' ! -name 'FastqSummaryF1L1.txt' ! -name 'NTC_*' !  -name 'BLANK_*' ! -name 'Undetermined_*' | sed -e 's#.*/##' -e 's/_.*//' | sort -u > fofn.txt


batchDir "$source_dir" "$dest_dir"
fofn batchDir.out
order_fastqs fofn.out batchDir.out
