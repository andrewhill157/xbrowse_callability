import subprocess
import os
import glob
import gzip

# Path to sample BED files
BED_FILES_PATH = '/humgen/atgu1/fs03/DM-Lab/projects/Muscle-Disease/North-64/callable/'

# Paths to exome capture interval files on Broad cluster
AGILENT_EXOME_CAPTURE_INTERVAL_PATH = '/humgen/atgu1/fs03/ahill/xbrowse_callability/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.bed'
ICE_EXOME_CAPTURE_INTERVAL_PATH = '/humgen/atgu1/fs03/ahill/xbrowse_callability/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.bed'

# Path to summary file
OUTPUT_FILE_PATH = 'summary.tsv'


def get_low_coverage_base_count(bed_gz_file):
    """
    Get list of low coverage entries from .bed.gz files
    """

    f = gzip.open(bed_gz_file)
    text = f.read()
    f.close()
    lines = [x.split('\t') for x in text.split('\n') if 'LOW_COVERAGE' in x]

    return sum([int(x[2]) - int(x[1]) for x in lines])


def get_overlapping_base_count(file1, file2):
    """
    Get number of overlapping bases between two BED files.
    """

    pipe = subprocess.Popen(['bedtools', 'intersect',
                             '-a', file1,
                             '-b', file2,
                             '-wo'], stdout=subprocess.PIPE)

    text = pipe.communicate()[0]
    lines = [x.split('\t') for x in text.split('\n')]
    return sum([int(x[9]) for x in lines if 'LOW_COVERAGE' in x])


########################################################################################################################
# Script
########################################################################################################################
if __name__ == '__main__':
    # Get all BED files on cluster
    bed_files = glob.glob(os.path.join(BED_FILES_PATH, '*.bam.bed.gz'))

    # Calculate how many bases fall outside each exome capture region for each sample BED file
    print 'Processing BED Files...'
    sample_entries = []
    for path in bed_files:
        low_coverage_base_count = get_low_coverage_base_count(path)
        agilent_overlap = get_overlapping_base_count(path, AGILENT_EXOME_CAPTURE_INTERVAL_PATH)
        ice_overlap = get_overlapping_base_count(path, ICE_EXOME_CAPTURE_INTERVAL_PATH)

        # Assemble summary information
        sample_id = os.path.splitext(os.path.basename(path))[0]
        summary = [sample_id, str(low_coverage_base_count), str(agilent_overlap), str(ice_overlap)]
        sample_entries.append(summary)
        print summary

    # Write to summary file
    with open(OUTPUT_FILE_PATH, 'w') as f:
        f.write('\n'.join(['\t'.join(x) for x in sample_entries]))