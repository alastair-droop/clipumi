#!/usr/bin/env python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
from sys import exit, stderr
import os.path
import gzip
from signal import signal, SIGPIPE, SIG_DFL
import logging
from collections import OrderedDict

# Handle broken pipes:
signal(SIGPIPE, SIG_DFL) 

# Get the version:
version = {}
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'version.py')) as f: exec(f.read(), version)

def main():
    # Create the CLI:
    defaults = {'verbosity':'warning', 'umi':6, 'offset':0, 'mismatch':1}
    parser = argparse.ArgumentParser(description='Identify and remove UMI tags from sequence starts')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(version))
    parser.add_argument('-V', '--verbose', dest='verbosity_level', default=defaults['verbosity'], choices=['error', 'warning', 'info', 'debug'], help='Set logging level (default {verbosity})'.format(**defaults))
    parser.add_argument('-a', '--return-all', dest='include_invalid', action='store_true', default=False, help='return all sequences, even if not mapped')
    parser.add_argument('-n', '--umi-length', dest='umi_length', type=int, metavar='n', default=defaults['umi'], help='length of the UMI sequence (default {umi})'.format(**defaults))
    parser.add_argument('-o', '--max-offset', dest='max_offset', type=int, metavar='n', default=defaults['offset'], help='maximum offset (default {offset})'.format(**defaults))
    parser.add_argument('-m', '--max-mismatch', dest='max_mismatch', type=int, metavar='n', default=defaults['mismatch'], help='maximum permissible primer mismatches (default {mismatch})'.format(**defaults))
    parser.add_argument('-f', '--output-fastq', dest='output_fastq', default=None, help='output FASTQ file')
    parser.add_argument(dest='input_primers', metavar='<primers>', help='primer FASTA file')
    parser.add_argument(dest='input_fastq', metavar='<fastq>', help='FASTQ file to process')
    args = parser.parse_args()

    # Log an error and exit:
    def error(message, exit_code=1):
        log.error(message)
        exit(exit_code)

    # Load a read (4 lines) from file:
    def readSequence(f):
        header = f.readline().decode('ASCII')
        if header == '': return None
        header = header.strip().lstrip('@')
        sequence = f.readline().decode('ASCII').strip()
        f.readline()
        quality = f.readline().decode('ASCII').strip()
        return (header, sequence, quality)

    # Return the hamming distance of two same-length strings:
    def mismatches(x, y):
        if len(x) > len(y): raise ValueError('y too short')
        diff = 0
        for i in range(len(x)):
            if x[i] != y[i]: diff += 1
        return diff

    # Set up logging based on the verbosity level set by the command line arguments:
    log = logging.getLogger()
    log_handler = logging.StreamHandler()
    log.default_msec_format = ''
    log_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    log.setLevel(args.verbosity_level.upper())
    log.addHandler(log_handler)

    # Write the settings to the log:
    log.info('using UMI length of {}'.format(args.umi_length))
    log.info('using maximum offset of {}'.format(args.max_offset))
    log.info('allowing <={} mismatches'.format(args.max_mismatch))
    if args.include_invalid: log.info('including unselected sequences')
    else: log.info('removing unselected sequences')

    # Read in the primers from FASTA format:
    log.info('reading primers from "{}"'.format(args.input_primers))
    primers = OrderedDict()
    max_primer_len = 0
    try:
        with open(args.input_primers, 'rt') as primer_file:
            while True:
                line = primer_file.readline()
                if line == '': break
                line = line.strip()
                if line.startswith('>'):
                    header = line.lstrip('>')
                    continue
                primers[header] = line
                max_primer_len = max(max_primer_len, len(line))
    except Exception as err: error(str(err))
    log.debug('read {} primers'.format(len(primers)))
    for primer in primers.keys(): log.debug('primer {}: {}'.format(primer, primers[primer]))

    # Attempt to open the input FASTQ file:
    try:
        input_fastq = gzip.open(args.input_fastq, 'r')
        input_fastq.readline().decode('ASCII').strip()
        f_type = 'compressed'
    except OSError:
        try:
            input_fastq = open(args.input_fastq, 'rb')
            input_fastq.readline().decode('ASCII').strip()
            f_type = 'plain'
        except: error('failed to open file "{}"'.format(args.input_fastq))
    input_fastq.seek(0)
    log.info('reading sequences from {} file "{}"'.format(f_type, args.input_fastq))

    # Attempt to create an output FASTQ file, if requested:
    if args.output_fastq is not None:
        try:
            if os.path.splitext(args.output_fastq)[1] == '.gz':
                output_type = 'compressed'
                output_f = gzip.open(args.output_fastq, 'wt')
            else:
                output_type = 'plain'
                output_f = open(args.output_fastq, 'wt')
        except: error('failed to open file "{}"'.format(args.output_fastq))
        log.info('writing clipped sequences to {} file "{}"'.format(output_type, args.output_fastq))
    else: output_f = None

    # A function to test a set of primers against a specific offset:
    def checkPrimers(seq, primers, offset, max_mismatch):
        best_primers = []
        best_mismatch = max_mismatch
        for primer_id, primer_seq in primers.items():
            mismatch = mismatches(primer_seq, seq[offset:])
            if mismatch < best_mismatch:
                best_primers = [primer_id]
                best_mismatch = mismatch
            elif mismatch == best_mismatch:
                best_primers.append(primer_id)
        if len(best_primers) == 0: return (None, [])
        return (best_mismatch, best_primers)    

    # Store frequent variables to save time:
    max_mismatch = args.max_mismatch
    max_offset = args.max_offset
    umi_length = args.umi_length
    write_fastq = output_f is not None
    return_invalid = args.include_invalid

    # Iterate through all the reads processing them one at a time:
    read_n = 0
    while True:
        # if read_n > 1: break
        read_n += 1
        read = readSequence(input_fastq)
        if read is None: break
        log.debug('read {}: {}'.format(read_n, read[1]))
        # Select the optimal (primer, offset) pair:
        selected_primer = None
        selected_primer_mismatches = None
        for offset in range(max_offset + 1):
            start_pos = offset + umi_length
            primer_mismatches = checkPrimers(read[1], primers, offset=start_pos, max_mismatch=max_mismatch)
            if primer_mismatches[0] is None: continue
            elif len(primer_mismatches[1]) > 1:
                log.debug('read {}: multiple primers found with {} mismatches'.format(read_n, primer_mismatches[0]))
                break
            else:
                selected_primer = primer_mismatches[1][0]
                selected_primer_mismatches = primer_mismatches[0]
                log.debug('read {}: primer {} ({} mismatches at offset {}) selected'.format(read_n, selected_primer, selected_primer_mismatches, offset))
                break
        # If we're not told to write out the invalid data, then we can skip to the end here:
        if (selected_primer is None):
            if len(primer_mismatches[1]) == 0: log.debug('read {}: no primers selected'.format(read_n))
            if (return_invalid is False): continue
        # If necessary, write out the FASTQ:
        if write_fastq is True:
            output_f.write('@{}\n'.format(read[0]))
            output_f.write('{}\n'.format(read[1][start_pos:]))
            output_f.write('+\n')
            output_f.write('{}\n'.format(read[2][start_pos:]))
        # Write out the mapping data:
        selected_UMI = read[1][offset:(offset + umi_length)]
        print('{}\t{}\t{}\t{}\t{}'.format(read[0].split(' ')[0], selected_UMI, selected_primer, selected_primer_mismatches, offset))
    log.info('processed {} reads'.format(read_n))

# If run as main, run main():
if __name__ == '__main__': main()
