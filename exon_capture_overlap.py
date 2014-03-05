import banyan
import os
import glob
import gzip

# Path to sample BED files
BED_FILES_PATH = '/humgen/atgu1/fs03/DM-Lab/projects/Muscle-Disease/North-64/callable/'

# Paths to exome capture interval files on Broad cluster
AGILENT_EXOME_CAPTURE_INTERVAL_PATH = '/humgen/atgu1/fs03/lek/resources/gatk/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list'
ICE_EXOME_CAPTURE_INTERVAL_PATH = '/seq/references/HybSelOligos/whole_exome_illumina_coding_v1/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list'

# Path to summary file
OUTPUT_FILE_PATH = 'summary.tsv'

def get_low_coverage_intervals(bed_gz_file):
    """
    Get list of low coverage entries from .bed.gz files
    """

    f = gzip.open(bed_gz_file)
    text = f.read()
    lines = text.split('\n')
    return [x.split('\t') for x in lines if 'LOW_COVERAGE' in x]


def get_interval_list(interval_list_file):
    """
    Get list of lists where each sub-list contains each column for a given interval entry in the interval list file.
    """

    with open(interval_list_file, 'r') as f:
        file_lines = f.read().strip().split('\n')

    return [x.split('\t') for x in file_lines if '@' not in x]


def create_interval_trees(interval_list):
    """
    Creates a dictionary with one entry per chromosome. Each entry will contain an interval tree for intervals on the
    respective chromosome.
    """

    chromosome_dict = dict()

    for interval in interval_list:
        chromosome_number = interval[0]
        start = int(interval[1])
        end = int(interval[2])

        if not chromosome_dict.has_key(chromosome_number):
            chromosome_dict[chromosome_number] = banyan.SortedSet(key_type=(int, int), updator=banyan.OverlappingIntervalsUpdator)
        chromosome_dict[chromosome_number].add([start, end])

    return chromosome_dict


def get_non_overlapping_count(interval1, interval2):
    """
    Return number of bases from interval1 that do not overlap with interval2.
    """

    start1 = interval1[0]
    start2 = interval2[0]
    end1 = interval1[1]
    end2 = interval2[1]

    if start1 < start2:
        return start2 - start1
    if end1 > end2:
        return end1 - end2

    return 0

########################################################################################################################
# Script
########################################################################################################################
if __name__ == '__main__':
    # Get all BED files on cluster
    bed_files = glob.glob(os.path.join(BED_FILES_PATH, '*.bam.bed.gz'))

    # Build interval trees for exome capture intervals
    print 'Building Interval Trees...'
    agilent_exome_capture_intervals = get_interval_list(AGILENT_EXOME_CAPTURE_INTERVAL_PATH)
    ice_exome_capture_intervals = get_interval_list(ICE_EXOME_CAPTURE_INTERVAL_PATH)

    agilent_exome_capture = create_interval_trees(agilent_exome_capture_intervals)
    ice_exome_capture = create_interval_trees(ice_exome_capture_intervals)

    # Calculate how many bases fall outside each exome capture region for each sample BED file
    print 'Processing BED Files...'
    sample_entries = []
    for path in bed_files:

        low_coverage_base_count = 0
        agilent_non_overlapping_count = 0
        ice_non_overlapping_count = 0

        low_coverage_intervals = get_low_coverage_intervals(path)

        for interval in low_coverage_intervals:
            chromosome_number = interval[0]
            start = int(interval[1])
            stop = int(interval[2])
            bounds = [start, stop]

            low_coverage_base_count += stop - start + 1

            try:
                agilent_exome_capture_overlap = agilent_exome_capture[chromosome_number].overlap(bounds)[0]
                agilent_non_overlapping_count += get_non_overlapping_count(bounds, agilent_exome_capture_overlap)
            except:
                # No overlapping exome capture interval, add total length
                agilent_non_overlapping_count += stop - start + 1

            try:
                ice_exome_capture_overlap = ice_exome_capture[chromosome_number].overlap(bounds)[0]
                ice_non_overlapping_count += get_non_overlapping_count(bounds, ice_exome_capture_overlap)
            except:
                # No overlapping exome capture interval, add total length
                ice_non_overlapping_count += stop - start + 1

        # Assemble summary information
        sample_id = os.path.splitext(os.path.basename(path))[0]
        summary = [sample_id, str(low_coverage_base_count), str(agilent_non_overlapping_count), str(ice_non_overlapping_count)]
        sample_entries.append(summary)
        print summary

    # Write to summary file
    with open(OUTPUT_FILE_PATH, 'w') as f:
        f.write('\n'.join(['\t'.join(x) for x in sample_entries]))