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


def generate_bed_intersect(beds_dir, out_file, dm6_file, intersects_dir=None,
                           flank=None, filename_tail='_peaks.bed'):
    files = next(os.walk(beds_dir))[2]
    bed_objs = [BedTool(os.path.join(beds_dir, f)) for f in files]
    arr = np.zeros((len(files),) * 2)
    for (x, y), val in np.ndenumerate(arr):
        if not flank:
            a = bed_objs[x]
            b = bed_objs[y]
        else:
            a = bed_objs[x].slop(g=dm6_file, b=flank)
            b = bed_objs[y].slop(g=dm6_file, b=flank)
        intersect = a.intersect(b, u=True)
        proportion = 1.0 * intersect.count() / bed_objs[x].count()
        arr[x][y] = proportion
        if intersects_dir and x != y:
            int_file = '{}_intersect_{}.bed'.format(files[x][:files[x].index(filename_tail)],
                                                    files[y][:files[y].index(filename_tail)])
            intersect.saveas(os.path.join(intersects_dir, int_file))
    with open(out_file, 'w') as file:
        filewriter = csv.writer(file, delimiter='\t')
        filewriter.writerow([None] + [f[:f.index(filename_tail)]
                                      for f in files])
        for i, row in enumerate(arr):
            filewriter.writerow([files[i][:files[i].index(filename_tail)]] +
                                [round(x * 100, 2) for x in row])

if __name__ == '__main__':
    print 'Generate intersect files and summary using *dm6* genome.'
    flank = ''
    while not flank.isdigit():
        flank = raw_input('Enter slop size (+/- n base pairs on each side): ')
    flank = int(flank)
    filename = raw_input('Enter results filename: ')
    out_file = os.path.join(BASE_DIR, 'results', filename)
    generate_bed_intersect(BEDS_DIR, out_file, DM6_CHROM_SIZES,
                           INTERSECTS_DIR, flank=flank)
    print 'File saved at "{}".'.format(out_file)
