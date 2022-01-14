# PFP-eBWT
PFP-eBWT is a tool that computes the eBWT of string collections using the Prefix-Free Parsing method (PFP).

# Usage

### Construction of the eBWT:
```
usage: pfpebwt [-h] [-w WSIZE] [-p MOD] [-t T] [-n N] [-r] [-v] [-d] [--reads]
               [--period] [--invert] [--keep] [--parsing]
               input

positional arguments:
  input                 input fasta file name

optional arguments:
  -h, --help            show this help message and exit
  -w WSIZE, --wsize WSIZE
                        sliding window size (def. 10)
  -p MOD, --mod MOD     hash modulus (def. 100)
  -t T                  number of helper threads (def. None)
  -n N                  number of different primes (def. 2)
  -r                    store the eBWT in RLE format (def. False)
  -v                    verbose (def. False)
  -d                    use different remainders instead of different prime numbers (def. False)
  --reads               process input as a read collection (def. False)
  --period              remove sequences which are not primitive (def. False)
  --invert              invert the eBWT (def. False)
  --keep                keep auxiliary files (debug only)
  --parsing             stop after the parsing phase (debug only)
```
The default algorithm will run with one prime number and one remainder to search for the trigger strings. If you want to build the eBWT of a collection of short sequences you have to use the `--reads` flag, and set with `-n` the maximum number of different primes allowed. With `-d` you can use different remainders instead of different primes. 
You can activate the `--period` flag in case your dataset contains non-primitive words, this flag will allow to filter them out. 
You can activate the `--invert` flag in case you want to invert the eBWT and write the inverted sequences to a file, this process requires to compute the WT of the eBWT in the internal memory. 

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
python3 pfpebwt yeast.fasta -w 10 -p 100 
```
# External resources

* [gSACA-K](https://github.com/felipelouza/gsa-is.git)
* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)

# Citation 

Please, if you use this tool in an academic setting cite the following paper:

### PFP-eBWT
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

# Authors

### Theoretical results:

* Christina Boucher
* Davide Cenzato
* Zsuzsanna Lipták
* Massimiliano Rossi
* Marinella Sciortino

### Implementation and experiments:

* [Davide Cenzato](https://github.com/davidecenzato) 
* [Massimiliano Rossi](https://github.com/maxrossi91)

# References

[1] Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini: *Prefix-Free Parsing for Building Big BWTs.* In Proc. of the 18th International Workshop on Algorithms in Bioinformatics, WABI 2018.

[2] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini and Taher Mun: *Prefix-free parsing for building big BWTs.* Algorithms Mol. Biol. 14(1): 13:1-13:15 (2019)

[3] Ge Nong, Sen Zhang and Wai Hong Chan: *Two Efficient Algorithms for Linear Time Suffix Array Construction.* IEEE Trans. Computers 60(10): 1471-1484 (2011)

[4] Hideo Bannai, Juha Kärkkäinen, Dominik Köppl and Marcin Piatkowski: *Constructing the bijective and the extended Burrows-Wheeler-Transform in linear time.* In Proc. of the 32nd Annual Symposium on Combinatorial Pattern Matching, CPM 2021.