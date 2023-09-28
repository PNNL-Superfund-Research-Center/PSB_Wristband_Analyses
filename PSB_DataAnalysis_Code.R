# packages that must be installed: #
# pcaMethods
# factoextra
# trelliscopejs
# ggrepel
# ggpubr

library(tidyverse)

# read in data #
psb_data = read_csv("PSB_data.csv")

# group by study #
nested_data = psb_data %>% group_by(Study) %>% nest()

# calculate number of wristbands per study #
nested_data2 = nested_data %>% mutate(NSamps= map_dbl(data, function(x) length(unique(x$ID))))

#### Figure 1 Analyses ####

# take data with at least one observation per chemical, log transform and replace BLOD measurements with 1/2 LOD #
nested_data2 = nested_data2 %>%
  mutate(
    filtered_data_ones = map(data, function(x) {
      temp <- x %>% group_by(ChemID) %>% summarise(Nobs = sum(!is.na(Log2_PicoMoles)))
      temp <- temp %>% filter(Nobs == 1)
      x %>% filter(!(ChemID %in% temp$ChemID)) 
    }),
    log_data = filtered_data_ones,
    log_half_LOD_transform = map(log_data, function(x) {
      x %>% mutate(Log2_PicoMoles = ifelse(is.na(Log2_PicoMoles), LOD_Log2PM * 0.5, Log2_PicoMoles))
    })
  )

# set a random seed for reproducibility #
set.seed(1137)

# PCA calculation function #
pca_return <- function(dataset) {
  
  pre_matrix <- dataset %>%
    dplyr::select(ID, ChemID, Log2_PicoMoles) %>%
    # make columns the chemical and rows the wristbands #
    pivot_wider(id_cols = ID, names_from = ChemID, values_from = Log2_PicoMoles) 
  theMatrix <- as.matrix(pre_matrix[,2:ncol(pre_matrix)])
  row.names(theMatrix) <- pre_matrix$ID
  
  # do PCA with scaling the log concentrations #
  pcaMethods::pca(theMatrix, method = "ppca", scale = "uv")
  
}

# return results for non-imputed data #
log_t = pca_return(nested_data2$log_data[[1]])

# return results for imputed data #
log_t_lod = pca_return(nested_data2$log_half_LOD_transform[[1]])


## Fig 1 ##
## make the pca plots ##
temp1 = data.frame(SampleID = row.names(log_t@scores),
                   PC1 = log_t@scores[,1], 
                   PC2 = log_t@scores[,2])

# run silhouette scores on scaled PCAs #
# 4 clusters is optimal based on silhoutte #
factoextra::fviz_nbclust(log_t@scores, kmeans, method = "silhouette")

# determine clusters for non-imputed PCA data #
temp1$Cluster = kmeans(log_t@scores, centers = 4)$cluster

# compile results for PCA based on imputed data #
temp2 = data.frame(SampleID = row.names(log_t_lod@scores),
                   PC1i = log_t_lod@scores[,1], 
                   PC2i = log_t_lod@scores[,2])

# join results from the two analyses #
pc_joined = temp1 %>% left_join(temp2)
pc_joined$ID = as.character(1:22)
pc_joined$Cluster = as.character(pc_joined$Cluster)

# pull R^2 values for each PCA #
temp1_r2s = log_t@R2
temp2_r2s = log_t_lod@R2

# PCA plots #
p1 = ggplot(data = pc_joined, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = ID)) +
  theme_bw() +
  xlab(paste("PC1 (", round(100*temp1_r2s[1],2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(100*temp1_r2s[2],2), "%)", sep = "")) +
  guides(color = "none")

p2 = ggplot(data = pc_joined, aes(x = PC1i, y = PC2i, color = Cluster)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = ID)) +
  theme_bw() +
  xlab(paste("PC1 (", round(100*temp2_r2s[1],2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(100*temp2_r2s[2],2), "%)", sep = "")) +
  guides(color = "none")

# get loadings values #
load_vals = data.frame(Imputation = rep(c("None", "1/2 LOD"), each = nrow(log_t@loadings)),
                       Chemical = c(row.names(log_t@loadings), row.names(log_t_lod@loadings)),
                       Load_PC1 = c(log_t@loadings[,1], log_t_lod@loadings[,1])
)

# make data.frame for loadings in format needed for geom_segments #
load_segs = load_vals %>% pivot_wider(values_from = Load_PC1, names_from = Imputation)

# order chemicals by PCA (no imputation) PC1 laodings #
chem_ord = load_vals %>% filter(Imputation == "None") %>% arrange(Load_PC1)

# define factors for plotting #
load_vals$Chemical = factor(load_vals$Chemical, levels = chem_ord$Chemical)
load_segs$Chemical = factor(load_segs$Chemical, levels = chem_ord$Chemical)
load_vals$Imputation = factor(load_vals$Imputation, labels = c("1/2 * LOD", "None"))

# PC1 loading plot #
p3 = ggplot(data = load_vals, aes(x = Chemical, y = Load_PC1)) +
  geom_segment(data = load_segs, aes(x = Chemical, xend = Chemical, y = None, yend = `1/2 LOD`)) +
  geom_point(size = 2, aes(color = Imputation)) +
  geom_hline(yintercept = 0, lty = 2, color = "darkgrey") +
  theme_bw() +
  xlab("") +
  ylab("PC1 Loading") +
  theme(axis.title = element_text(size = 13), axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5), legend.title = element_text(size = 13), legend.text = element_text(size = 11), legend.position = "bottom")

# arrange the plots #
ggpubr::ggarrange(ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = T, legend = "right", labels = "AUTO"), p3, ncol = 1, nrow = 2, labels = c("", "C"), vjust = -1.2)

ggsave("/Users/bram489/Desktop/Fig1.png", width = 7, height = 6, units = "in")


#### Figure 2 ####

# The mice method is described here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7550746/
# And the code is pulled from here: https://github.com/BTDPhD/Frontiers-Wristbands-NHBCS/blob/master/20200731%20-%20GitHub%20-%20Wristband%20Descriptives%20-%20Analyses%20-%2012W.R

# Mice Impute: 
# RF Impute: Samples rows, Chemicals Columns 


# First filter to data with at most 40% missingness, then add mice imputing, and random forest imputing
nested_data3 = nested_data2 %>%
  mutate(
    
    # Limit to compounds where only 40% is missing
    forty_missing = map(log_data, function(x) {
      
      compounds_to_impute <- x %>% 
        group_by(ChemID) %>% 
       summarise(Missingness = sum(is.na(Log2_PicoMoles)) / length(Log2_PicoMoles)) %>% 
        dplyr::filter(Missingness <= 0.4) %>% 
        dplyr::select(ChemID) %>% 
        unlist()
      
      return(x %>% dplyr::filter(ChemID %in% all_of(compounds_to_impute)))
      
    }), 
    
    # Conduct the mice imputing 
    mice_impute = map(forty_missing, function(x) {
      
      # Create dataset, MICE impute, add whether it is BLOD or Matrix Interference
      pre_matrix <- x %>% 
        dplyr::select(ChemID, Log2_PicoMoles, ID) %>%
        tidyr::pivot_wider(names_from = ChemID, values_from = Log2_PicoMoles, id_cols = ID) 

      # Build imputation matrix
      impute <- pre_matrix %>% 
        dplyr::select(-ID) %>% 
        unlist() %>%
        matrix(nrow = nrow(pre_matrix)) %>%
        mice(m = 25, printFlag = FALSE)
      
      # Take the median per row 
      medianImp <- do.call(rbind, 1:length(impute$imp) %>% lapply(function(num) {
        
        tibble(
          ChemID = colnames(pre_matrix[2:ncol(pre_matrix)])[num],
          Log2_PicoMoles = apply(impute$imp[[num]], 1, median)
        )
        
      })) %>%
        filter(!is.na(Log2_PicoMoles)) %>%
        dplyr::mutate(
          ID = apply(pre_matrix, 2, function(x) {pre_matrix$ID[is.na(x)]}) %>% unlist(),
          Imputed = TRUE
        )
      
      # Make data.frame of missingness
      x %>% 
        dplyr::select(ChemID, ID, Result_Qualifier) %>%
        merge(medianImp, by = c("ChemID", "ID")) %>%
        dplyr::select(ChemID, Log2_PicoMoles, Result_Qualifier, ID, Imputed) %>%
        return()
      
    }),
    
    # Random forest imputation
    rf_impute = map(forty_missing, function(x) {
      
      # Create dataset, rf impute, add whether it is BLOD or Matrix Interference
      pre_matrix <- x %>% 
        dplyr::select(ChemID, Log2_PicoMoles, ID) %>%
        tidyr::pivot_wider(names_from = ChemID, values_from = Log2_PicoMoles, id_cols = ID) 
      
     # Impute
      impute <- pre_matrix %>% 
        dplyr::select(-ID) %>% 
        unlist() %>%
        matrix(nrow = nrow(pre_matrix)) %>%
        missForest()
      
      # Pull out imputation matrix
      imp_matrix <- impute$ximp
      colnames(imp_matrix) <- colnames(pre_matrix)[2:ncol(pre_matrix)]
      row.names(imp_matrix) <- pre_matrix$ID
      
      # Make data.frame of missingness
      x %>% 
        dplyr::select(ChemID, Log2_PicoMoles, Result_Qualifier, ID) %>%
        dplyr::filter(is.na(Log2_PicoMoles)) %>%
        mutate(
          Imputed = TRUE, 
          Log2_PicoMoles = map2(ID, ChemID, function(rownm, colnm) {
            imp_matrix[rownm, colnm]
          }) %>% unlist()
        ) %>% 
        return()
      
    })
    
  )


# pull results for two imputation methods #
or_imp_rf = nested_data3[3,] %>% select(rf_impute) %>% unnest(rf_impute) %>% rename(RF_Imputed = Log2_PicoMoles)
or_imp_mice = nested_data3[3,] %>% select(mice_impute) %>% unnest(mice_impute) %>% rename(Mice_Imputed = Log2_PicoMoles)

# join results #
or_imp = or_imp_rf %>% inner_join(or_imp_mice)

# pull LOD info #
lod_inf = nested_data2 %>% filter(Study == "OR") %>% unnest(data) %>% filter(Result_Qualifier == "BLOD") %>% select(ChemID, LOD_Log2PM) %>% unique() %>% ungroup() %>% summarise(LODmin = min(LOD_Log2PM), LODmed = median(LOD_Log2PM), LODmean = mean(LOD_Log2PM), LODmax = max(LOD_Log2PM)) 

ggplot(data = or_imp, aes(x = RF_Imputed, y = Mice_Imputed)) +
  geom_point(data = subset(or_imp, Result_Qualifier == "BLOD"), aes(color = "Below LOD")) +
  geom_point(data = subset(or_imp, ResultQualifier == "MI"), aes(color = "MI")) +
  geom_vline(xintercept = lod_inf$LODmed, lty = 2) +
  geom_hline(yintercept = lod_inf$LODmed, lty = 2) +
  geom_abline(a = 0, b = 1, lty = 2) +
  xlab("RF Imputed Log Concentration") +
  ylab("MICE Imputed Log Concentration") +
  scale_color_manual(name='Missing Data Type',
                     breaks=c('Below LOD', 'MI'),
                     values=c('Below LOD'='darkgrey', 'MI'='blue')) +
  theme_bw() +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13), legend.title = element_text(size = 13), legend.text = element_text(size = 11))



#### Figure 3 ####

# subset to Oregon study #
or_subset = nested_data2 %>% select(Study, data) %>% ungroup() %>% unnest() %>% filter(Study == "OR") %>% group_by(ChemID) %>% nest()

# make histograms of every chemical for exploratory purposes #
plot_fn0 = function(x){
  
  x2 = x %>% mutate(value = ifelse(is.na(Log2_PicoMoles) & Result_Qualifier == "BLOD", 2^LOD_Log2PM, 2^Log2_PicoMoles))
  
  ggplot(data = x2, aes(x = value)) +
    geom_histogram(bins = 25) +
    theme_bw() +
    xlab("Concentration (pmol)") +
    ylab("Count")
}

# make histograms of every chemical for exploratory purposes #
plot_fn = function(x){
  
  x2 = x %>% mutate(value = ifelse(is.na(Log2_PicoMoles) & Result_Qualifier == "BLOD", LOD_Log2PM, Log2_PicoMoles))
  
  ggplot(data = x2, aes(x = value)) +
    geom_histogram(bins = 25) +
    theme_bw() +
    xlab("Log Concentration (pmol)") +
    ylab("Count")
}

# map plots #
or_subset2 = or_subset %>% mutate(panels = trelliscopejs::map_plot(data, plot_fn))

# compute missing values for different types #
or_subset2 = or_subset2 %>% mutate(N_MI = map_dbl(data, function(x) sum(x$Result_Qualifier == "MI")), N_BLOD = map_dbl(data, function(x) sum(x$Result_Qualifier == "BLOD")))


plot_ids = c(grep("chemH", or_subset2$ChemID), grep("chemYB", or_subset2$ChemID), grep("chemPB", or_subset2$ChemID), grep("CB", or_subset2$ChemID))

p1 = plot_fn0(or_subset2$data[[plot_ids[1]]])
p2 = plot_fn(or_subset2$data[[plot_ids[1]]])
p3 = plot_fn(or_subset2$data[[plot_ids[3]]])
p4 = plot_fn(or_subset2$data[[plot_ids[4]]])


ggpubr::ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, labels = "AUTO")


#### Figure 4 ####
# from oregon data #
# calculate number of missing observations per chemical #
# order from fewest to most missing values #
or_nmiss = or_subset %>% mutate(Nmiss = map_dbl(data, function(x) sum(is.na(x$Log2_PicoMoles)))) %>% arrange(Nmiss)

# get data from top 3 chemicals with least missing values #
or_top3 = or_subset %>% filter(ChemID %in% or_nmiss$ChemID[1:3]) %>% unnest(data)

# format data for calculation of correlations #
or_top3_long = or_top3 %>% dplyr::select(ChemID, ID, Log2_PicoMoles) %>% pivot_wider(names_from = ChemID, values_from = Log2_PicoMoles)

# chemUB and chemT #
rmv_nas_1 = which(is.na(or_top3_long$chemUB) | is.na(or_top3_long$chemT))

# chemH and chemT #
rmv_nas_2 = which(is.na(or_top3_long$chemH) | is.na(or_top3_long$chemT))

# subset down to non-missing values for two cases
or_cor1 = or_top3_long[-rmv_nas_1, -which(names(or_top3_long) == "chemH")]
or_cor2 = or_top3_long[-rmv_nas_2, -which(names(or_top3_long) == "chemUB")]

orig_cor1 = cor(or_cor1$chemT, or_cor1$chemUB, method = "spearman")
orig_cor2 = cor(or_cor2$chemH, or_cor2$chemT, method = "spearman")

set.seed(4819)
samp_pool = 1:nrow(or_cor1)
new_cor_c1 = NULL
new_cor_ig_c1 = NULL
for(i in 1:(nrow(or_cor1) - 5)){
  temp_dat = or_cor1
  temp_dat_ig = or_cor1
  rndm_s = sample(samp_pool, i)  
  # fill in chemT with 1/2 LOD #
  temp_dat[rndm_s, "chemT"] = -2.244858*0.5
  # fill in chemUB with 1/2 LOD #
  temp_dat[rndm_s, "chemUB"] = -0.1403715*0.5
  # turn them into NAs for ignore missing val case #
  temp_dat_ig[rndm_s, "chemT"] = NA
  temp_dat_ig[rndm_s, "chemUB"] = NA
  
  new_cor_c1[i] = cor(temp_dat$chemT, temp_dat$chemUB, method = "spearman")
  new_cor_ig_c1[i] = cor(temp_dat_ig$chemT, temp_dat_ig$chemUB, use = "pairwise.complete.obs", method = "spearman")
}

cor_res_df_c1 = data.frame(Case = "Case 1", N = seq(from = nrow(or_cor1), to = 5, by = (-1)), Cor_val = c(orig_cor1, new_cor_c1), Cor_val_ig = c(orig_cor1, new_cor_ig_c1)) %>% mutate(Detection_Perc = N/N[1] * 100)


samp_pool = 1:nrow(or_cor2)
new_cor_c2 = NULL
new_cor_ig_c2 = NULL
for(i in 1:(nrow(or_cor2) - 5)){
  temp_dat = or_cor2
  temp_dat_ig = or_cor2
  rndm_s = sample(samp_pool, i)  
  temp_dat[rndm_s, "chemH"] = -2.82985*0.5
  temp_dat[rndm_s, "chemT"] = -2.244858*0.5
  temp_dat_ig[rndm_s, "chemH"] = NA
  temp_dat_ig[rndm_s, "chemT"] = NA
  
  new_cor_c2[i] = cor(temp_dat$chemT, temp_dat$chemH, method = "spearman")
  new_cor_ig_c2[i] = cor(temp_dat_ig$chemT, temp_dat_ig$chemH, use = "pairwise.complete.obs", method = "spearman")
}

cor_res_df_c2 = data.frame(Case = "Case 2", N = seq(from = nrow(or_cor2), to = 5, by = (-1)), Cor_val = c(orig_cor2, new_cor_c2), Cor_val_ig = c(orig_cor2, new_cor_ig_c2)) %>% mutate(Detection_Perc = N/N[1] * 100)

cor_res = bind_rows(cor_res_df_c1, cor_res_df_c2)

cor_res_long = cor_res %>% pivot_longer(3:4) %>% rename(Computation = name)

cor_res_long$Computation = factor(cor_res_long$Computation, levels = c("Cor_val", "Cor_val_ig"), labels = c("1/2 * LOD", "Ignore"))


p1 = ggplot(data = subset(cor_res_long, Case == "Case 1"), aes(x = Detection_Perc, y = value, color = Computation, group = Computation)) +
  geom_vline(xintercept = 70, lty = 2, color = "darkgrey", size = 1.25) +
  geom_hline(yintercept = orig_cor1, color = "darkgrey", size = 1.25) +
  geom_line() +
  scale_y_continuous(limits = c(-0.3,1), breaks = c(-0.2,0,0.2,0.4,0.6,0.8,1)) +
  theme_bw() +
  xlab("% Detection Frequency") +
  ylab("Spearman Correlation") +
  guides(color = guide_legend(title = "Missing Value\nTreatment")) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13), legend.text = element_text(size = 11), legend.title = element_text(size = 13))

p2 = ggplot(data = subset(cor_res_long, Case == "Case 2"), aes(x = Detection_Perc, y = value, color = Computation, group = Computation)) +
  geom_vline(xintercept = 70, lty = 2, color = "darkgrey", size = 1.25) +
  geom_hline(yintercept = orig_cor2, color = "darkgrey", size = 1.25) +
  geom_line() +
  theme_bw() +
  scale_y_continuous(limits = c(0.9,1)) +
  xlab("% Detection Frequency") +
  ylab("Spearman Correlation") +
  guides(color = guide_legend(title = "Missing Value\nTreatment")) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13), legend.text = element_text(size = 11), legend.title = element_text(size = 13))

ggpubr::ggarrange(p1,p2, ncol = 2, labels = "AUTO", common.legend = T, legend = "bottom")



#### Figure 5 ####
ny_miss_summ = nested_data2 %>% filter(Study == "NY Pilot") %>% select(data) %>% unnest() %>% group_by(ChemID) %>% summarise(NBLOD = sum(Result_Qualifier == "BLOD")) %>% arrange(NBLOD)

ny_z = nested_data2 %>% filter(Study == "NY Pilot") %>% select(data) %>% unnest() %>% filter(ChemID == "chemZ")

ny_z = ny_z %>% mutate(value = ifelse(is.na(Log2_PicoMoles) & Result_Qualifier == "BLOD", LOD_Log2PM, Log2_PicoMoles))

ggplot(data = ny_z, aes(x = value)) +
  geom_histogram(bins = 15) +
  theme_bw()

set.seed(16)
samp_ids = sample(1:nrow(ny_z), nrow(ny_z), replace = T)

ny_z_samp = ny_z[samp_ids,]

ggplot(data = ny_retene_samp, aes(x = value)) +
  geom_histogram(bins = 15) +
  theme_bw()

ny_both = data.frame(Data = rep(c("Original", "Sampled"), each = nrow(ny_z)), bind_rows(ny_z, ny_z_samp))


ggplot(data = ny_both, aes(x = value, group = Data, color = Data)) +
  geom_density() +
  theme_bw()

set.seed(1)
ny_z_samp2 = list()
for(i in 1:25){
  samp_ids = sample(1:nrow(ny_z), nrow(ny_z), replace = T)
  
  ny_z_samp2[[i]] = data.frame(Sample = i, ny_z[samp_ids,])
  
}

ny_z_samp2 = do.call(rbind, ny_z_samp2) 
ggplot(data = ny_z, aes(x = value)) +
  geom_density(data = ny_z_samp2, aes(x = value, group = factor(Sample), color = Sample)) +
  geom_density(color = "red", size = 2) +
  theme_bw() +
  ylab("Density") +
  xlab("Log2 Concentration") +
  theme(axis.text = element_text(size= 11), axis.title = element_text(size = 13)) +
  guides(color = "none")

#### Figure 6 ####

# subset to NY and OR data #
nyor_data = nested_data2 %>% filter(Study %in% c("NY", "OR")) %>% select(Study, log_half_LOD_transform) %>% unnest(log_half_LOD_transform)

# subset to chemCB
chm_cb_sub = nyor_data %>% filter(ChemID == "chemCB")

# calculate tertile values #
threshs = chm_cb_sub %>% summarise(T1 = quantile(Log2_PicoMoles, 1/3, na.rm = T), T2 = quantile(Log2_PicoMoles, 2/3, na.rm = T))

ggplot(data = chm_cb_sub, aes(x = Log2_PicoMoles, color = Study, group = Study)) +
  stat_ecdf(linewidth = 1.25) +
  geom_hline(yintercept = 1/3, col = "darkgrey", lty = 2) +
  geom_hline(yintercept = 2/3, col = "darkgrey", lty = 2) +
  theme_bw() +
  xlab("Log2 Concentration") +
  ylab("Percentile") +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13), legend.text = element_text(size = 11), legend.title = element_text(size = 13))
