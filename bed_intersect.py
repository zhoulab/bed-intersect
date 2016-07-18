import os
import csv
import numpy as np
from pybedtools import BedTool

BASE_DIR = os.path.dirname(os.getcwd())
BEDS_DIR = os.path.join(BASE_DIR, 'data/Beds')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
INTERSECTS_DIR = os.path.join(RESULTS_DIR, 'intersects')
if not os.path.exists(INTERSECTS_DIR):
    os.mkdir(INTERSECTS_DIR)

DM6_CHROM_SIZES = os.path.join(BASE_DIR, 'data/dm6.chrom.sizes')


def generate_bed_intersect(beds_dir, table_file, dm6_file,
                           list_results_file, threshold=None,
                           intersects_dir=None, slop=None,
                           filename_tail='_peaks.bed'):
    files = next(os.walk(beds_dir))[2]
    bed_objs = [BedTool(os.path.join(beds_dir, f))
                for f in files if f.endswith(filename_tail)]
    arr = np.zeros((len(files),) * 2)
    for (x, y), val in np.ndenumerate(arr):
        if not slop:
            a = bed_objs[x]
            b = bed_objs[y]
        else:
            a = bed_objs[x].slop(g=dm6_file, b=slop)
            b = bed_objs[y].slop(g=dm6_file, b=slop)
        intersect = a.intersect(b, u=True)
        proportion = 1.0 * intersect.count() / bed_objs[x].count()
        arr[x][y] = proportion
        if intersects_dir and x != y:
            int_file = '{}_intersect_{}.bed'.format(files[x].replace(filename_tail, ''),
                                                    files[y].replace(filename_tail, ''))
            intersect.saveas(os.path.join(intersects_dir, int_file))

    with open(table_file, 'w') as file1:
        filewriter = csv.writer(file1, delimiter='\t')
        filewriter.writerow([None] + [f.replace(filename_tail, '')
                                      for f in files])
        for i, row in enumerate(arr):
            filewriter.writerow([files[i].replace(filename_tail, '')] +
                                [round(x * 100, 2) for x in row])

    with open(list_results_file, 'w') as file2:
        filewriter = csv.writer(file2, delimiter='\t')
        filewriter.writerow(['File1', 'File2',
                             'File1.intersect(File2)_percent'])
        for (x, y), val in np.ndenumerate(arr):
            if threshold and val * 100 > threshold and x is not y:
                filewriter.writerow([files[x].replace(filename_tail, ''),
                                     files[y].replace(filename_tail, ''),
                                     round(val * 100, 2)])


def get_intersect(bed1, bed2, dm6_file, slop=None, intersects_dir=None):
    if slop:
        bed1 = bed1.slop(g=dm6_file, b=slop)
        bed2 = bed2.slop(g=dm6_file, b=slop)
    return bed1.intersect(bed2, u=True)


if __name__ == '__main__':
    print 'Generate intersect files and summary using *dm6* genome.'
    slop = ''
    while not slop.isdigit():
        slop = raw_input('Enter slop size (+/- n base pairs on each side): ')
    slop = int(slop)
    table_fname = raw_input('Enter table results filename: ')
    table_file = os.path.join(BASE_DIR, 'results', table_fname)
    threshold = ''
    while not threshold.isdigit():
        threshold = raw_input('Enter threshold for filtered results (0-100): ')
    threshold = int(threshold)
    filtered_file = raw_input('Enter filtered results filename: ')
    filtered = os.path.join(BASE_DIR, 'results', filtered_file)

    generate_bed_intersect(BEDS_DIR, table_file, DM6_CHROM_SIZES, filtered,
                           threshold=threshold, intersects_dir=INTERSECTS_DIR,
                           slop=slop)
    print 'Table file saved at "{}".'.format(table_file)
    print 'Filtered results file saved at "{}".'.format(filtered_file)
