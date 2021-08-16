#===============================================
# Create list of duplicates to remove
#===============================================

library(pacman)
p_load(data.table, magrittr, stringi)
related <- fread("data/tmp/adc_dup.con")
setnames(related, c("ID1", "ID2"), c("IID1", "IID2"))
related[, PI_HAT := (2 * N_IBS2 + N_IBS1) / (2 * (N_IBS0 + N_IBS1 + N_IBS2))]
related[, pair := 1:.N]
related <- related[, .(IID1, IID2, pair, PI_HAT)]
related_long <- melt.data.table(related, measure.vars = c("IID1", "IID2"), value.name = "IID")

# identify ids in multiple pairs
dup_iids <- related_long$IID[duplicated(related_long$IID)]
setorder(related_long, pair, variable)

# related pairs with uniq ids (i.e. one id not the other + "_dup")
# note: all "_dup" ids under IID2 column
uniq_id_pairs <- related[!(grep("_dup", IID2)), pair]
related_long_uniq <- related_long[pair %in% uniq_id_pairs]
related_long_dup <- related_long[!(pair %in% uniq_id_pairs)]
# #-------------------------------------------------------
# # examining individuals related to more than one person
# #-------------------------------------------------------
# 
# # don't need to worry about 3+ person clusters as removing person with highest 
# # miss in each pair still leaves just one person in each cluster
# # (because as long as name is flagged at least once, PLINK will remove)
# 
# # 3 sets of 2 pairs that have one shared person in both pairs.
# # - 1 PO and one 2nd-degree pair
# # - 1 PO and one 2nd-degree pair
# # - 2 2nd-degree pairs
# # In the interest of maximizing sample size, 
# # I will preferentially remove the individuals present
# # in 2 pairs (i.e. NACC903503, NACC697927, NACC696497).
# related_one_dup <- related[(IID1 %in% dup_iids & !(IID2 %in% dup_iids)) |
#                             (!(IID1 %in% dup_iids) & IID2 %in% dup_iids)
# ]
# related_one_dup.dups <- intersect(dup_iids, c(related_one_dup$IID1, related_one_dup$IID2))

#----------------------------------
# pairs with uniq iids
# add miss and demographics info
#----------------------------------

miss <- fread("data/tmp/adc_dup.imiss")

uds <- fread("/data_global/nacc/investigator_nacc53.csv", header = T, na.strings = c(-4, "999", "9999", "888")) %>% 
  .[NACCVNUM == NACCAVST] %>% 
  setnames(., "NACCID", "IID") %>% 
  merge(related_long_uniq, ., "IID") %>% 
  .[, .(IID, pair, PI_HAT, NPFORMVER, NPSEX, BIRTHMO, BIRTHYR, HISPANIC, RACE)]

mds <- fread("/data_global/nacc/fardo09062019.csv", na.strings = c(-4, "999", "9999", "888")) %>%
  setnames(., "NACCID", "IID") %>% 
  merge(related_long_uniq, ., "IID") %>% 
  .[, .(IID, pair, PI_HAT, NPFORMVER, NPSEX, MDSBIRTHMO, MDSBIRTHYR, MDSHISPANIC, MDSRACE)] 

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
pairs_not_remove <- nacc_info_wide[`1` == `2`, pair]
# second name in pairs matched pairs both with non-missing NPFORMVER
extra_iids <- nacc[(pair %in% pairs_not_remove) & !(is.na(NPFORMVER)),][duplicated(pair), IID]
iids_uniq_remove <- nacc[!(pair %in% pairs_not_remove), IID]
#---------------------
# pairs with dup iids
# check missingness
#---------------------

dup_miss <- merge(related_long_dup, miss, by = "IID")
setorder(dup_miss, pair)
dup_miss[, num := rep_len(1:2, .N)]
dup_miss_wide <- dcast(dup_miss[, .(pair, num, N_MISS)], pair ~ num, value.var = 'N_MISS')
non_matching_pairs <- dup_miss_wide[`1` != `2`, pair]
# iids with higher missingness
iid_miss <- dup_miss[pair %in% non_matching_pairs] %>% 
  setorder(., pair, N_MISS) %>% 
  .[seq(2, .N, 2), IID]
dup_miss_remove <- dup_miss[!(pair %in% non_matching_pairs)
                            ][grep("_dup", IID), IID]

iid_remove <- data.table(FID = 0,
                         IID = c(extra_iids, iids_uniq_remove, dup_miss_remove, iid_miss)
                         )
write.table(iid_remove, file = "data/tmp/adc_dup_remove.tmp", quote = F, row.names = F, col.names = F)

rm(list = ls())
p_unload(all)
