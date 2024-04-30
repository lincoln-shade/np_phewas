##--------------------------------------------------
## Create bisulfite conversion efficiency covariate
##--------------------------------------------------


library(data.table)
library(magrittr)
library(synapser)
library(illuminaio)

synLogin()

covars <- fread(synapser::synGet("syn5843544")$path, 
                integer64 = "double")

idat_file_folders_dir <- "data/synapse/DNA_methylation/IDAT Files/"
idat.file.folders <- list.files(idat_file_folders_dir)

##--------------------------------------
## Create data.table of control labels
## Methylation450 control probes
## taken from manifest file
##--------------------------------------
control.probes <- c("22711390", "22795447", "56682500", "46651360", 
                    "24637490", "33665449", "54705438", "49720470", 
                    "26725400", "57693375", "15700381", "33635504", 
                    "43720395", "70664314", "71718498", "30724412")

# `I` and `II` mean probe type 1 and 2
# C means converted and U means unconverted
label <- c("I-C1", "I-C2", "I-C3", "I-U1", 
           "I-U2", "I-U3", "I-C4", "I-C5", 
           "I-C6", "I-U4", "I-U5", "I-U6", 
           "II-1", "II-2", "II-3", "II-4")

control.probe.labels <- data.table(rn=control.probes, label)

##-------------------------------------------------
## function for creating control probe data.table
##-------------------------------------------------
Make.Control.DT <- function(path, controls) {
  idat <- readIDAT(path)
  idat.quants <- as.data.table(idat[["Quants"]], keep.rownames=T)
  controls <- merge(idat.quants, controls, by="rn")
  return(controls)
}
##-----------------------------------------------------------
## Function for calculating bisulfite conversion efficiency
##-----------------------------------------------------------
Calc.Bisulfite.Conversion.Efficiency <- function(green, red, offset=100) {
  Calc.Type.I.Probe.Efficiency <- function(dt, offset) {
    # 1
    T1.1 <- dt[label == "I-C1", Mean] / 
      (dt[label == "I-C1", Mean] + dt[label == "I-U1", Mean] + offset)
    
    # 2
    T1.2 <- dt[label == "I-C2", Mean] / 
      (dt[label == "I-C2", Mean] + dt[label == "I-U2", Mean] + offset)
    
    # 3
    T1.3 <- dt[label == "I-C3", Mean] / 
      (dt[label == "I-C3", Mean] + dt[label == "I-U3", Mean] + offset)
    
    # 4
    T1.4 <- dt[label == "I-C4", Mean] / 
      (dt[label == "I-C4", Mean] + dt[label == "I-U4", Mean] + offset)
    
    # 5
    T1.5 <- dt[label == "I-C5", Mean] / 
      (dt[label == "I-C5", Mean] + dt[label == "I-U5", Mean] + offset)
    
    # 6
    T1.6 <- dt[label == "I-C6", Mean] / 
      (dt[label == "I-C6", Mean] + dt[label == "I-U6", Mean] + offset)
    
    return(c(T1.1, T1.2, T1.3, T1.4, T1.5, T1.6))
  }
  Calc.Type.II.Probe.Efficiency <- function(green, red, offset) {
    #1 
    T2.1 <- red[label == "II-1", Mean] / 
      (red[label == "II-1", Mean] + green[label == "II-1", Mean] + offset)
    
    #2 
    T2.2 <- red[label == "II-2", Mean] / 
      (red[label == "II-2", Mean] + green[label == "II-2", Mean] + offset)
    
    #3 
    T2.3 <- red[label == "II-3", Mean] / 
      (red[label == "II-3", Mean] + green[label == "II-3", Mean] + offset)
    
    #4 
    T2.4 <- red[label == "II-4", Mean] / 
      (red[label == "II-4", Mean] + green[label == "II-4", Mean] + offset)
    
    return(c(T2.1, T2.2, T2.3, T2.4))
  }
  control.probes.efficiency <- c(Calc.Type.I.Probe.Efficiency(green, offset), 
                                 Calc.Type.I.Probe.Efficiency(red, offset), 
                                 Calc.Type.II.Probe.Efficiency(green, red, offset))
  if (sum(control.probes.efficiency > 1) > 0) {
    stop("At least one efficiency greater than 1")
  }
  
  return(median(control.probes.efficiency))
}

##------------------------------
## Loop for creating covariate
##------------------------------
specimenID <- character(length = nrow(covars))
conversion.efficiency <- numeric(length = nrow(covars))
# for (barcode in idat.file.folders) {
#   sentrix.position <- covars[Sentrix_ID == barcode, Sentrix_Position]
#   
#   for (array.position in sentrix.position) {
#     path.green.idat <- paste0(idat_file_folders_dir, 
#                               barcode, "/", barcode, "_", array.position, "_Grn.idat")
#     path.red.idat <- paste0(idat_file_folders_dir, 
#                             barcode, "/", barcode, "_", array.position, "_Red.idat")
#     if (file.exists(path.green.idat) & file.exists(path.red.idat)) {
#       green.controls <- Make.Control.DT(path.green.idat, control.probe.labels)
#       red.controls <- Make.Control.DT(path.red.idat, control.probe.labels)
#       
#       specimenID <- c(specimenID, covars[Sentrix_ID == barcode & Sentrix_Position == array.position, 
#                                          Sample])
#       conversion.efficiency <- c(conversion.efficiency, 
#                                  Calc.Bisulfite.Conversion.Efficiency(green.controls, 
#                                                                       red.controls))
#     }
#   }
# }

for (i in 1:nrow(covars)) {
  barcode <- covars[i, Sentrix_ID]
  array.position <- covars[i, Sentrix_Position]
  idat_file_prefix <- paste0(idat_file_folders_dir, 
                             barcode, "/", barcode, "_", array.position)
  path.green.idat <- paste0(idat_file_prefix, "_Grn.idat")
  path.red.idat <- paste0(idat_file_prefix, "_Red.idat")
  if (file.exists(path.green.idat) & file.exists(path.red.idat)) {
    green.controls <- Make.Control.DT(path.green.idat, control.probe.labels)
    red.controls <- Make.Control.DT(path.red.idat, control.probe.labels)
    conversion.efficiency[i] <- Calc.Bisulfite.Conversion.Efficiency(green.controls, red.controls)
  } else {
    conversion.efficiency[i] <- NA
  }
}

bisulfite.conversion.efficiency <- data.table(covars$Sample, conversion.efficiency)
save(bisulfite.conversion.efficiency, file = "data/bisulfite.conversion.efficiency.RData")

