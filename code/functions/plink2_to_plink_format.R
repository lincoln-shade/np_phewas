#=========================================================================
# Convert plink2 .glm.logistic data.tables to a format that is backwards
# compatible with scripts that use plink1.9 .assoc.logistic data.tables
#=========================================================================

plink2_to_plink_format <- function(dt, ci = 0.95, keep_plink2_cols = TRUE) {
  require(data.table)
  tmp <- copy(dt)
  setnames(tmp, 
           c('#CHROM', 'POS', 'ID', 'OBS_CT', 'LOG(OR)_SE', 'Z_STAT'),
           c('CHR', 'BP', 'SNP', 'NMISS', 'SE', 'STAT'))
  
  if (ci) {
    CI <- as.integer(ci * 100)
    setcolorder(tmp, c('CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 'OR', 'SE',
                       paste0('L', CI), paste0('U', CI), 'STAT', 'P'))
  } else {
    setcolorder(tmp, c('CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 'OR', 'SE',
                       'STAT', 'P'))
  }
  
  if (!keep_plink2_cols) {
    tmp[, `:=`(REF = NULL, ALT = NULL, ERRCODE = NULL)]
  }
  
  return(tmp)
}
