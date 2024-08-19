# PFP-eBWT
PFP-eBWT is a tool for computing the eBWT and optionally the Generalized Conjugate Array (GCA) of string collections using the Prefix-Free Parsing preprocessing (PFP).

# Usage

### Construction of the eBWT:
```
usage: pfpebwt [-h] [-w WSIZE] [-p MOD] [-t T] [-n N] [--rle] [--samples]
               [--GCA] [--reads] [--remainders] [--period] [--invert] [--keep]
               [--parsing] [--verbose]
               input

Tool to compute the eBWT and the GCA of a string collection.

positional arguments:
  input                 input fasta file name

optional arguments:
  -h, --help            show this help message and exit
  -w WSIZE, --wsize WSIZE
                        sliding window size (def. 10)
  -p MOD, --mod MOD     hash modulus (def. 100)
  -t T                  number of helper threads (def. None)
  -n N                  number of different primes when using --reads (def. 1)
  --rle                 store eBWT in RLE format (def. False)
  --samples             compute the GCA samples (def. False)
  --GCA                 compute the complete GCA (def. False)
  --reads               process input ad a reads multiset (def. False)
  --remainders          use multiple remainders instead of multiple primes (def. False)
  --period              remove periodic sequences (def. False)
  --invert              invert the eBWT (def. False)
  --keep                keep auxiliary files (debug only)
  --parsing             stop after the parsing phase (debug only)
  --verbose             verbose (def. False)
```
The default PFP algorithm will run with one prime number and one remainder to search for the trigger strings. If you need to compute the eBWT of a collection of short sequences, you can use the `--reads` flag to enable the algorithm to search for trigger strings using multiple prime numbers. With the `-n` flag, you can set the maximum number of primes allowed. With `--remainders`, you can enable the algorithm to use multiple remainders (and one prime) instead of multiple primes (and one remainder). 
You can activate the `--period` flag if your dataset contains non-primitive words; this flag will run an algorithm to filter them out. 
You can activate the `--invert` flag to invert the eBWT and write the inverted string collection to a file. This process requires computing the WT of the eBWT in the internal memory. 
You can activate the `--samples` flag to compute the GCA samples at the end and beginning of an eBWT run, the `--GCA` flag enable the computation of the complete GCA. 
The `.info` file contains more information on the output files.

# Example
### Download and Compile

```console
git clone https://github.com/davidecenzato/PFP-eBWT.git
cd PFP-eBWT
mkdir build
cd build
cmake ..
make
```

### Run on Example Data

```console
// Build the eBWT on a toy data set
python3 pfpebwt yeast.fasta 
// Build the eBWT and the GCA samples on a toy data set
python3 pfpebwt yeast.fasta --samples
// Build the eBWT and the complete GCA on a toy data set
python3 pfpebwt yeast.fasta --GCA
```
# External resources

* [gSACA-K](https://github.com/felipelouza/gsa-is.git)
* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)

# References and citations

[1] Christina Boucher, Davide Cenzato, Zsuzsanna Lipták, Massimiliano Rossi, Marinella Sciortino:
Computing the Original eBWT Faster, Simpler, and with Less Memory. SPIRE 2021: 129-142 ([go to the paper](https://link.springer.com/chapter/10.1007/978-3-030-86692-1_11)) 

[2] Christina Boucher, Davide Cenzato, Zsuzsanna Lipták, Massimiliano Rossi, Marinella Sciortino:
Computing the original eBWT faster, simpler, and with less memory. CoRR abs/2106.11191 (2021) ([go to the paper](https://arxiv.org/abs/2106.11191)) 

[3] Christina Boucher, Davide Cenzato, Zsuzsanna Lipták, Massimiliano Rossi, Marinella Sciortino:
r-indexing the eBWT. Inf. Comput. 298: 105155 (2024) ([go to the paper](https://doi.org/10.1016/j.ic.2024.105155)) 

Please, if you use this tool in an academic setting cite the following papers:

### conference version
    @inproceedings{BoucherCL0S21a,
      author    = {Christina Boucher and
                   Davide Cenzato and
                   Zsuzsanna Lipt{\'{a}}k and
                   Massimiliano Rossi and
                   Marinella Sciortino},
      title     = {Computing the Original e{BWT} Faster, Simpler, and with Less Memory},
      booktitle = {Proceedings of 28th International Symposium in String Processing and Information Retrieval {SPIRE} 2021},
      series    = {Lecture Notes in Computer Science},
      volume    = {12944},
      pages     = {129--142},
      year      = {2021}
    }

### extended version
    @article{extBoucherCL0S21a,
      author       = {Christina Boucher and
                      Davide Cenzato and
                      Zsuzsanna Lipt{\'{a}}k and
                      Massimiliano Rossi and
                      Marinella Sciortino},
      title        = {Computing the original eBWT faster, simpler, and with less memory},
      journal      = {CoRR},
      volume       = {abs/2106.11191},
      year         = {2021},
    }

### journal paper containing the GCA computation algorithm
    @article{BoucherCLRS24,
      author       = {Christina Boucher and
                      Davide Cenzato and
                      Zsuzsanna Lipt{\'{a}}k and
                      Massimiliano Rossi and
                      Marinella Sciortino},
      title        = {r-indexing the eBWT},
      journal      = {Inf. Comput.},
      volume       = {298},
      pages        = {105155},
      year         = {2024}
    }

# Authors

* Christina Boucher
* Davide Cenzato
* Zsuzsanna Lipták
* Massimiliano Rossi
* Marinella Sciortino

### Implementation and experiments:

* [Davide Cenzato](https://github.com/davidecenzato) 
* [Massimiliano Rossi](https://github.com/maxrossi91)