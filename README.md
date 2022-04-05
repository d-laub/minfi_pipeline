# Dependencies
R dependencies for the `minfi` preprocessing script:

- minfi
- minfiData
- IlluminaHumanMethylationEPICanno.ilm10b4.hg19
- IlluminaHumanMethylationEPICmanifest
- IlluminaHumanMethylation450kanno.ilmn12.hg19
- IlluminaHumanMethylation450kmanifest
- sva
- EnsDb.Hsapiens.v86

Optional dependencies:

- arrow

Also need to install [maxprobes](https://github.com/markgene/maxprobes), see link for installation instructions. If installing with conda, use install_github(..., dependencies = FALSE) as R does not detect packages installed through conda and the above dependencies covers the needed packages.

# Usage
`lib_minfi.R` needs to be in the same directory as `run_minfi.R` when you call it. It can then be called as:

`Rscript run_minfi.R <idat_type> <idat_dir> <info_file> <out_dir>`

- `idat_type` should be one of: Both, 450k, or EPIC
- `idat_dir` is the directory to the .idat files
- `info_file` is the path to a csv of phenotypes which must include the columns:
    - Basename
    - sex, with the characters M,F denoting the sex of the samples
- `out_dir` is a directory to put output files

Note that this script will not exclude arbitrary files so those have be removed either beforehand by e.g. editing the `info_file` to only include the Basenames/samples of interest or afterward from the output. In other words, the samples that failed the bench QC checks should be excluded from the `info_file` csv.