# Detection and Removal of initial UMI Tags from FASTQ Files


`clip-umi` detects and removes UMI tags from the beginning of FASTQ reads using the presence of primer sequences to locate the UMI tags.

## Overview

`clip-umi` assumes that all reads have an initial UMI followed immediately by one of several primer sequences. A small offset can be specified to allow for sequence before the start of the UMI. The read is therefore:

~~~
ooUUUUUUpppppppppSSSSSSSSS
~~~

* `o` marks possible offset bases;
* `U` marks the UMI tag to be detected and clipped;
* `p` marks one of several primers; and
* `S` marks the remaining read sequence

If a primer is detected at a given offset, the preceding UMI is detected, and the UMI (plus offset) is removed. Thus, the returned read would be:

~~~
pppppppppSSSSSSSSS
~~~

### Primer Matching

Primer matching starts at an offset of zero. For each offset, each primer is assigned a mismatch score based on how well it matches the read. Primers matching with scores above the threshold are rejected. The primer(s) with the lowest match score are examined:

* If no primers match, the next offset is tried;
* If exactly one primer matches, it is selected;
* If more than one primers match; the read is rejected (as we are unable to assign a primer unambiguously)

## Usage

~~~
usage: clipumi [-h] [-v] [-V {error,warning,info,debug}] [-a] [-n n] [-o n]
                [-m n] [-f OUTPUT_FASTQ]
                <primers> <fastq>

Identify and remove UMI tags from sequence starts

positional arguments:
  <primers>             primer FASTA file
  <fastq>               FASTQ file to process

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -V {error,warning,info,debug}, --verbose {error,warning,info,debug}
                        Set logging level (default warning)
  -a, --return-all      return all sequences, even if not mapped
  -n n, --umi-length n  length of the UMI sequence (default 6)
  -o n, --max-offset n  maximum offset (default 0)
  -m n, --max-mismatch n
                        maximum permissible primer mismatches (default 1)
  -f OUTPUT_FASTQ, --output-fastq OUTPUT_FASTQ
                        output FASTQ file
~~~

## Installation

Installation should be as simple as:

~~~bash
git clone https://github.com/alastair-droop/clipumi.git
cd clipumi
python setup.py install
~~~

If you do not have admin privileges, you can install this locally using `python setup.py install --user`.

After installation, you can verify that you have the correct version using `clipumi -v`.

## Licence

These tools are released under the [GNU General Public License version 3](http://www.gnu.org/licenses/gpl.html).
