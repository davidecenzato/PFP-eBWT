# PFP-eBWT
Tool to build the eBWT of string collections using the Prefix Free Parse (PFP).

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
  -r                    store eBWT in RLE format (def. False)
  -v                    verbose (def. False)
  -d                    use remainders instead of primes (def. False)
  --reads               process input ad a reads multiset (def. False)
  --period              remove periodic sequences (def. False)
  --invert              invert the eBWT (def. False)
  --keep                keep auxiliary files (debug only)
  --parsing             stop after the parsing phase (debug only)
```
The default algorithm will run with one prime and one remainder to check if a window is a trigger string. If you want to build the eBWT of a collection of short sequences you have to use the `--reads` flag, and set with `-n` the maximum number of primes or remainders (with `-n`) allowed.
You can activate the `--period` in case your dataset contains non-primitive word, this flag will allow to filter them out. 

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

[1] Christina Boucher, Davide Cenzato, Zsuzsanna Lipták, Massimiliano Rossi and Marinella Sciortino, *"Computing the original eBWT faster, simpler, and with less memory."*, (2021).

[2] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini, and Taher Mun. *"Prefix-free parsing for building big BWTs."*, Algorithms for Molecular Biology 14, no. 1 (2019): 13.