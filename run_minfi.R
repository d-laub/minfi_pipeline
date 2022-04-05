# library(arrow)
source('lib_minfi.R')

args = commandArgs(TRUE)
idat_type = args[1]
idat_dir = args[2]
pheno_file = args[3]
out_dir = args[4]

if (idat_type == "Both") {
    message("Running minfi for both 450k and EPIC .idat files.")
    gmset = process_both(idat_dir, pheno_file)
} else {
    message("Running minfi for ", idat_type, " .idat files.")
    gmset = process_data(idat_dir, pheno_file)
}

cg_ann = cg_annotations(gmset)
write_feather(cg_ann, file.path(out_dir, "cg_ann.fth"))

mvals = m_values(gmset)
write_feather(mvals, file.path(out_dir, "minfi_mvals.fth"))