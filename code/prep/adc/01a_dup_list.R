#===============================================
# Create list of duplicates to remove
#===============================================

library(data.table)
library(magrittr)
library(stringi)

related <- fread("tmp/adc_dup.con")
setnames(related, c("ID1", "ID2"), c("IID1", "IID2"))
related[, PI_HAT := (2 * N_IBS2 + N_IBS1) / (2 * (N_IBS0 + N_IBS1 + N_IBS2))]
related[, pair := 1:.N]
related <- related[, .(IID1, IID2, pair, PI_HAT)]
related_long <- melt.data.table(related, 
                                measure.vars = c("IID1", "IID2"), 
                                value.name = "IID")
related_long[, NACCID := stri_replace_first_regex(IID, ".*_", "")]

# identify ids in multiple pairs
dup_iids <- related_long$IID[duplicated(related_long$NACCID)]
dup_naccids <- related_long$NACCID[duplicated(related_long$NACCID)]
setorder(related_long, pair, variable)


related_long_uniq <- related_long


#----------------------------------
# pairs with uniq iids
# add miss and demographics info
#----------------------------------

miss <- fread("tmp/adc_dup.imiss")

uds <- fread("/data_global/nacc/investigator_nacc56.csv", 
             header = T, 
             na.strings = c(-4, "999", "9999", "888")) %>% 
  .[NACCVNUM == NACCAVST] %>% 
  # setnames(., "NACCID", "IID") %>% 
  merge(related_long_uniq, ., "NACCID") %>% 
  .[, .(NACCID, IID, pair, PI_HAT, NPFORMVER, NPSEX, 
        BIRTHMO, BIRTHYR, HISPANIC, RACE)]

mds <- fread("/data_global/nacc/fardo09062019.csv", 
             na.strings = c(-4, "999", "9999", "888")) %>%
  # setnames(., "NACCID", "IID") %>% 
  merge(related_long_uniq, ., "NACCID") %>% 
  .[, .(NACCID, IID, pair, PI_HAT, NPFORMVER, NPSEX, MDSBIRTHMO, 
        MDSBIRTHYR, MDSHISPANIC, MDSRACE)] 

colnames(mds) <- stri_replace_all_fixed(colnames(mds), "MDS", "")

# examine demographics of pairs
nacc <- rbindlist(list(uds, mds)) %>% 
  merge(., miss, "IID") %>% 
  # set order by pair and miss
  setorder(., pair, F_MISS) %>% 
  .[, FID := 0] %>% 
  setcolorder(., c("FID", "IID"))

# identify pairs for which birth year and birth month are equal
nacc[, num := rep_len(1:2, .N)]
nacc[, info_sum := BIRTHYR + 10000 * BIRTHMO]
nacc_info_wide <- dcast(nacc[, .(pair, IID, info_sum, num)], 
                        pair ~ num, 
                        value.var = 'info_sum'
                        )
pairs_same_birthdate <- nacc_info_wide[`1` == `2`, pair]

# For pairs where demographic data (birth date) is not the same,
# both iids should be removed (we don't know which is correct).
iids_diff_birthdate <- nacc[!(pair %in% pairs_same_birthdate), IID]
iids_np_remove <- nacc[(pair %in% pairs_same_birthdate) & !(is.na(NPFORMVER)),
                   ][duplicated(pair), IID]
iids_no_np_remove <- nacc[(pair %in% pairs_same_birthdate) & (is.na(NPFORMVER))
                          ][!(duplicated(pair)), IID]
#---------------------
# pairs with dup iids
# check missingness
#---------------------

# dup_miss <- merge(related_long_dup, miss, by = "IID")
# setorder(dup_miss, pair)
# dup_miss[, num := rep_len(1:2, .N)]
# dup_miss_wide <- dcast(dup_miss[, .(pair, num, N_MISS)],
#                        pair ~ num,
#                        value.var = 'N_MISS')
# non_matching_pairs <- dup_miss_wide[`1` != `2`, pair]

# iids with higher missingness
# iid_miss <- dup_miss[pair %in% non_matching_pairs] %>% 
#   setorder(., pair, N_MISS) %>% 
#   .[seq(2, .N, 2), IID]
# dup_miss_remove <- dup_miss[!(pair %in% non_matching_pairs)
#                             ][grep("_dup", IID), IID]

id_remove <- data.table(FID = 0,
                         IID = c(iids_diff_birthdate, 
                                 iids_np_remove, 
                                 iids_no_np_remove))

write.table(id_remove, 
            file = "tmp/adc_dup_remove.tmp", 
            quote = F, 
            row.names = F, 
            col.names = F)
