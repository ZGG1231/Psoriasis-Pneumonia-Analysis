==============================================================================
# Script 5: Bayesian Colocalization Analysis (05_Colocalization.R)
# Description: Performs coloc.abf and coloc.susie to identify shared causal variants between exposure and outcome traits within a +/- 500kb window of lead SNPs.
==============================================================================
library(coloc)
library(data.table)

# 1. Ensure global tables are in data.table format
setDT(gwas_pso)
setDT(pneu_mr_ready)

# Standardize column names
setnames(gwas_pso, old = c("rsid", "pval"), new = c("SNP", "P"), skip_absent = TRUE)

# Get the list of target SNPs
snp_list <- unique(res_single_snps$SNP) 

# 2. Initialize an empty data frame to collect results
all_coloc_results <- data.frame()

# 3. Run batch processing loop
for (i in 1:length(snp_list)) {
  my_target <- snp_list[i]
  
  # Find coordinates
  snp_info <- pneu_mr_ready[pneu_mr_ready$SNP == my_target, ]
  if(nrow(snp_info) == 0) next 
  
  chr_val <- snp_info$CHR[1]
  bp_val <- snp_info$BP[1]
  w_start <- bp_val - 500000
  w_end <- bp_val + 500000

  # Region slicing and deduplication
  reg_pso <- gwas_pso[CHR == chr_val & BP >= w_start & BP <= w_end]
  reg_pneu <- pneu_mr_ready[CHR == chr_val & BP >= w_start & BP <= w_end]
  
  reg_pso <- reg_pso[order(P)][!duplicated(SNP)]
  reg_pneu <- reg_pneu[order(P)][!duplicated(SNP)]
  
  # Align shared SNPs
  shared_snps <- intersect(reg_pso$SNP, reg_pneu$SNP)
  if(length(shared_snps) < 10) next
  
  reg_pso <- reg_pso[SNP %in% shared_snps]
  reg_pneu <- reg_pneu[SNP %in% shared_snps]
  setorder(reg_pso, SNP)
  setorder(reg_pneu, SNP)
  
  # Package data for coloc (using true BETA and varbeta = SE^2 directly, no MAF required)
  dataset_pso <- list(
    beta = reg_pso$BETA,
    varbeta = reg_pso$SE^2,
    N = 494544,                     
    s = 36466 / 494544,                     
    type = "cc",                     
    snp = reg_pso$SNP
  )
  
  dataset_pneu <- list(
    beta = reg_pneu$BETA,
    varbeta = reg_pneu$SE^2,
    N = 500348,                     
    s = 78791 / 500348,                     
    type = "cc",
    snp = reg_pneu$SNP
  )
  
  # Run Bayesian colocalization
  coloc_run <- tryCatch({
    coloc.abf(dataset1 = dataset_pso, dataset2 = dataset_pneu)
  }, error = function(e) {
    return(NULL) 
  })
  
  if (is.null(coloc_run)) next 
  
  # Extract and combine results
  res_row <- as.data.frame(t(coloc_run$summary))
  res_row$Target_SNP <- my_target
  res_row$CHR <- chr_val
  res_row$BP <- bp_val
  
  all_coloc_results <- rbind(all_coloc_results, res_row)
}

# 4. Format final results table and sort by PP.H4 descending
if (nrow(all_coloc_results) > 0) {
  all_coloc_results <- all_coloc_results[, c("Target_SNP", "CHR", "BP", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")]
  all_coloc_results <- all_coloc_results[order(-all_coloc_results$PP.H4.abf), ]
  
  write.csv(all_coloc_results, file = "Coloc_114IVs_Beta_Varbeta_Results.csv", row.names = FALSE)
}
###Forest Plot of Mendelian Randomization & Coloc Results
data <- read.delim(file.choose(), header = TRUE, sep = "\t")
library(ggplot2)
library(dplyr)
plot_data <- data %>%
   mutate(
      # Calculate OR and 95% confidence interval
      OR = exp(b),
      LCI = exp(b - 1.96 * se),
      UCI = exp(b + 1.96 * se),
      # Format P values and FDR for display on the plot (optional)
      p_label = sprintf("%.4f", p),
      FDR_label = sprintf("%.4f", FDR)
   )
# To make the forest plot more aesthetically pleasing, reorder SNPs by OR value
plot_data$SNP <- factor(plot_data$SNP, levels = plot_data$SNP[order(plot_data$OR)])

p <- ggplot(plot_data, aes(y = SNP, x = OR)) +
   # Draw confidence interval error bars
   geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2, color = "#2c3e50", size = 0.8) +
   # Draw effect size points (squares)
   geom_point(shape = 15, size = 3.5, color = "#e74c3c") +
   # Add null line (vertical line at OR = 1)
   geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray", size = 1) +
   # Use a clean academic theme
   theme_classic() +
   # Adjust fonts and axis labels
   theme(
      axis.text.y = element_text(size = 11, face = "bold", color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
   ) +
   # Set title and axis labels
   labs(
      title = "Forest Plot of Individual SNPs",
      x = "Odds Ratio (95% CI)",
      y = "SNP"
   )
print(p)

### MR scatter plot
library(TwoSampleMR)
mr_results <- mr(dat_clean_steiger)
p1 <- mr_scatter_plot(mr_results, dat_clean_steiger)
p1[[1]]

# Use ggsave to save as a high-resolution image
library(ggplot2)
ggsave("MR_scatter_plot.pdf", plot = p1[[1]], width = 7, height = 5, dpi = 300)


### Colocalization scatter comparison plot (LocusCompare Plot)
devtools::install_github("boxiangliu/locuscomparer")
library(locuscomparer)
colnames(gwas_pso)[colnames(gwas_pso) == "SNP"] <- "rsid"
colnames(gwas_pso)[colnames(gwas_pso) == "P"] <- "pval"
colnames(gwas_pneu)[colnames(gwas_pneu) == "SNP"] <- "rsid"
colnames(gwas_pneu)[colnames(gwas_pneu) == "P"] <- "pval"
gwas_pso_clean2 <- gwas_pso %>%
   group_by(rsid) %>%
   slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
   ungroup()
# 2. Deduplicate pneumonia data: group by rsid, select the row with the minimum pval in each group
gwas_pneu_clean2 <- gwas_pneu %>%
   group_by(rsid) %>%
   slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
   ungroup()
# Take the intersection of the two datasets by rsid and merge them
merged_data <- inner_join(gwas_pso_clean2, gwas_pneu_clean2, by = "rsid", suffix = c(".pso", ".pneu"))
# Convert P values to -log10 format, which is the standard form for plotting
merged_data <- merged_data %>%
   mutate(
      logP_pso = -log10(pval.pso),
      logP_pneu = -log10(pval.pneu)
   )
# Specify the core SNP to highlight
target_snp <- "rs35741374"            ####take rs35741374 as an example
# Plot using ggplot2
p_custom <- ggplot(merged_data, aes(x = logP_pso, y = logP_pneu)) +
   geom_point(color = "darkgray", size = 2, alpha = 0.6) +
   
   geom_point(data = filter(merged_data, rsid == target_snp), 
              color = "red", size = 4, shape = 18) +
   
   geom_text(data = filter(merged_data, rsid == target_snp),
             aes(label = rsid), vjust = -1.5, color = "darkred", fontface = "bold") +
   
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue", alpha = 0.5) +
   
   theme_bw() +
   labs(
      title = paste("Colocalization at", target_snp, "Region"),
      x = "-log10(P) for Psoriasis",
      y = "-log10(P) for Pneumonia"
   ) +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold")
   )

print(p_custom)

### locuszoom
library(locuscomparer)
library(dplyr)
library(data.table)

colnames(gwas_pso)[colnames(gwas_pso) == "SNP"] <- "rsid"
colnames(gwas_pso)[colnames(gwas_pso) == "P"] <- "pval"
colnames(gwas_pneu)[colnames(gwas_pneu) == "SNP"] <- "rsid"
colnames(gwas_pneu)[colnames(gwas_pneu) == "P"] <- "pval"

gwas_pso_clean2 <- gwas_pso %>%
  filter(!is.na(rsid) & !is.na(pval)) %>%
  group_by(rsid) %>%
  slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
  ungroup()

gwas_pneu_clean2 <- gwas_pneu %>%
  filter(!is.na(rsid) & !is.na(pval)) %>%
  group_by(rsid) %>%
  slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
  ungroup()

target_snp <- "rs35741374"        ##take rs35741374 as an example
snp_info <- gwas_pneu_clean2 %>% filter(rsid == target_snp)
target_chr <- snp_info$CHR[1]
target_bp <- snp_info$BP[1]

pso_region <- gwas_pso_clean2 %>%  
  filter(CHR == target_chr & BP >= (target_bp - 500000) & BP <= (target_bp + 500000)) %>%
  select(rsid, pval)

pneu_region <- gwas_pneu_clean2 %>%  
  filter(CHR == target_chr & BP >= (target_bp - 500000) & BP <= (target_bp + 500000)) %>%
  select(rsid, pval)

p <- locuscompare(
  in_fn1 = pso_region,  
  in_fn2 = pneu_region,  
  title1 = "Psoriasis GWAS",  
  title2 = "Pneumonia GWAS",
  snp = target_snp
)

ggsave(
  filename = "LocusCompare_rs35741374_new.pdf", 
  plot = p, 
  width = 15, 
  height = 8, 
  dpi = 300
)


#Susie.coloc
args <- commandArgs(trailingOnly=TRUE)
config_file <- if(length(args)) args[1] else "config.R"
source(config_file)
suppressPackageStartupMessages({
  library(data.table); library(readxl); library(susieR); library(coloc)
})
stopifnot("annotate_susie" %in% getNamespaceExports("coloc"))
dir.create(WORK_DIR,recursive=TRUE,showWarnings=FALSE)
WORK_DIR <- normalizePath(WORK_DIR)
dirs <- file.path(WORK_DIR,c("raw","clean","regions","results","tables"))
for(d in dirs) dir.create(d,showWarnings=FALSE)
raw_dir<-dirs[1]; clean_dir<-dirs[2]; region_dir<-dirs[3]
result_dir<-dirs[4]; table_dir<-dirs[5]
setDTthreads(THREADS)
command <- function(exe, argv) {
  rc <- system2(exe,args=vapply(as.character(argv),shQuote,""))
  if(rc!=0) stop("Command failed: ",exe)
}
for(p in c(IV_FILE,LIFTOVER,CHAIN,paste0(REFERENCE,c(".bed",".bim",".fam"))))
  if(!file.exists(p)) stop("Missing input: ",p)
if(!nzchar(Sys.which(PLINK)) && !file.exists(PLINK)) stop("PLINK not found")
file.copy(config_file,file.path(WORK_DIR,"config.used.R"),overwrite=TRUE)
writeLines(capture.output(sessionInfo()),file.path(WORK_DIR,"sessionInfo.txt"))
download <- function(url,name) {
  p<-file.path(raw_dir,name)
  if(!file.exists(p)) {
    command("wget",c("-c","-O",paste0(p,".part"),url))
    command("gzip",c("-t",paste0(p,".part")))
    if(!file.rename(paste0(p,".part"),p)) stop("Download rename failed")
  } else command("gzip",c("-t",p))
  p
}
ps_file<-download("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90472001-GCST90473000/GCST90472771/harmonised/GCST90472771.h.tsv.gz","GCST90472771.h.tsv.gz")
pn_file<-download("https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_J10_PNEUMONIA.gz","finngen_R12_J10_PNEUMONIA.gz")
iv_input<-if(grepl("\\.xlsx?$",IV_FILE,ignore.case=TRUE)) as.data.table(read_excel(IV_FILE)) else fread(IV_FILE)
if(!"SNP" %in% names(iv_input)) stop("IV file needs a SNP column")
iv_ids<-unique(trimws(as.character(iv_input$SNP)))
if(!length(iv_ids)||anyNA(iv_ids)||any(!grepl("^rs[0-9]+$",iv_ids))) stop("Invalid IV IDs")
if(!is.null(EXPECTED_IV_COUNT)&&length(iv_ids)!=EXPECTED_IV_COUNT)
  stop("IV count is ",length(iv_ids),", expected ",EXPECTED_IV_COUNT)
message("Using ",length(iv_ids)," IVs; no automatic truncation")
fwrite(data.table(SNP=iv_ids),file.path(WORK_DIR,"IVs_used.tsv"),sep="\t")
read_bim<-function(path) {
  x<-fread(path,col.names=c("CHR19","SNP","CM","BP19","A1","A2"))
  x[,CHR19:=suppressWarnings(as.integer(CHR19))]
  x[,`:=`(A1=toupper(A1),A2=toupper(A2))]
  x
}
original_ref<-read_bim(paste0(REFERENCE,".bim"))
coordinate_sources<-list(original_ref[SNP %chin% iv_ids & CHR19 %in% 1:22,.(IV_SNP=SNP,CHR19,BP19)])
rm(original_ref); gc()
panel<-file.path(clean_dir,"reference.maf005")
command(PLINK,c("--bfile",REFERENCE,"--autosome","--snps-only","just-acgt",
  "--nonfounders","--maf","0.005","--keep-allele-order","--threads",THREADS,
  "--memory",PLINK_MEMORY_MB,"--make-bed","--out",panel))
ref<-read_bim(paste0(panel,".bim"))
ref<-ref[CHR19 %in% 1:22 & A1 %chin% c("A","C","G","T") &
  A2 %chin% c("A","C","G","T") & A1!=A2 & !is.na(SNP) & !SNP %chin% c(".","NA","0")]
if(anyDuplicated(ref$SNP)) stop("Duplicate reference IDs")
ref[,`:=`(K1=pmin(A1,A2),K2=pmax(A1,A2))]
ref[,KEY_N:=.N,by=.(CHR19,BP19,K1,K2)]
ref<-ref[KEY_N==1,.(CHR19,BP19,K1,K2,PANEL_SNP=SNP,A1,A2)]

prepare_trait<-function(path,trait) {
  message("Cleaning/lifting ",trait)
  mapping<-if(trait=="Psoriasis") c(SNP="rsid",CHR38="chromosome",BP38="base_pair_location",
    EA="effect_allele",OA="other_allele",BETA="beta",SE="standard_error",P="p_value",
    EAF="effect_allele_frequency",N_EFF="cum_eff_sample_size") else
    c(SNP="rsids",CHR38="#chrom",BP38="pos",EA="alt",OA="ref",BETA="beta",SE="sebeta",P="pval",EAF="af_alt")
  x<-fread(cmd=paste("gzip -dc",shQuote(path)),select=unname(mapping),na.strings=c("NA","NaN",".",""))
  setnames(x,unname(mapping),names(mapping))
  if(trait=="Pneumonia") x[,N_EFF:=4*78791*421557/(78791+421557)]
  x[,SNP:=as.character(SNP)]
  x[,CHR38:=suppressWarnings(as.integer(sub("^chr","",CHR38)))]
  for(col in c("BP38","BETA","SE","P","EAF","N_EFF"))
    set(x,j=col,value=suppressWarnings(as.numeric(x[[col]])))
  x[,`:=`(EA=toupper(EA),OA=toupper(OA))]
  x<-x[CHR38 %in% 1:22 & is.finite(BP38) & BP38>0 & BP38==floor(BP38) &
    EA %chin% c("A","C","G","T") & OA %chin% c("A","C","G","T") & EA!=OA &
    is.finite(BETA) & is.finite(SE) & SE>0 & is.finite(P) & P>=0 & P<=1]
  x[!is.finite(EAF)|EAF<0|EAF>1,EAF:=NA_real_]
  x[!is.finite(N_EFF)|N_EFF<=0,N_EFF:=NA_real_]
  x[,ROW_ID:=.I]
  bed<-file.path(clean_dir,paste0(trait,".hg38.bed"))
  mapped<-file.path(clean_dir,paste0(trait,".hg19.mapped.bed"))
  unmapped<-file.path(clean_dir,paste0(trait,".unmapped.bed"))
  fwrite(x[,.(CHR=paste0("chr",CHR38),START=BP38-1,END=BP38,ROW_ID,SCORE=0,STRAND="+")],
    bed,sep="\t",col.names=FALSE)
  command(LIFTOVER,c(bed,CHAIN,mapped,unmapped))
  m<-fread(mapped,col.names=c("CHR19","START","END","ROW_ID","SCORE","LIFT_STRAND"))
  m[,CHR19:=suppressWarnings(as.integer(sub("^chr","",CHR19)))]
  m[,N_MAP:=.N,by=ROW_ID]
  m<-m[N_MAP==1 & CHR19 %in% 1:22 & END-START==1 & LIFT_STRAND %chin% c("+","-")]
  x<-merge(x,m[,.(ROW_ID,CHR19,BP19=START+1,LIFT_STRAND)],by="ROW_ID")
  # Track IV coordinates before excluding palindromic/low-NEFF SNPs.
  coords<-x[SNP %chin% iv_ids,.(IV_SNP=SNP,CHR19,BP19)]
  multi<-which(!is.na(x$SNP)&grepl("[;,]",x$SNP))
  if(length(multi)) {
    extra<-rbindlist(lapply(multi,function(i) {
      hit<-intersect(trimws(unlist(strsplit(x$SNP[i],"[;,]"))),iv_ids)
      if(length(hit)) data.table(IV_SNP=hit,CHR19=x$CHR19[i],BP19=x$BP19[i])
    }))
    coords<-rbindlist(list(coords,extra),use.names=TRUE,fill=TRUE)
  }
  x[LIFT_STRAND=="-",`:=`(EA=chartr("ACGT","TGCA",EA),OA=chartr("ACGT","TGCA",OA))]
  x<-x[!((EA=="A"&OA=="T")|(EA=="T"&OA=="A")|(EA=="C"&OA=="G")|(EA=="G"&OA=="C"))]
  x[,`:=`(K1=pmin(EA,OA),K2=pmax(EA,OA))]
  x<-merge(x,ref,by=c("CHR19","BP19","K1","K2"),sort=FALSE)
  x[,BETA_FLIPPED:=as.integer(EA!=A1)]
  x[BETA_FLIPPED==1,`:=`(BETA=-BETA,EAF=1-EAF)]
  x[,`:=`(ORIG_SNP=SNP,SNP=PANEL_SNP,EA=A1,OA=A2)]
  setorder(x,SNP,P,ROW_ID)
  x<-unique(x,by="SNP")
  x<-x[,.(SNP,CHR19,BP19,EA,OA,BETA,SE,P,EAF,N_EFF,ORIG_SNP,CHR38,BP38,LIFT_STRAND,BETA_FLIPPED)]
  setorder(x,CHR19,BP19,SNP)
  fwrite(x,file.path(clean_dir,paste0(trait,".hg19.matched.tsv")),sep="\t",na="NA")
  list(data=x,coords=unique(coords))
}
pa<-prepare_trait(ps_file,"Psoriasis"); coordinate_sources[[2]]<-pa$coords; A<-pa$data; rm(pa)
pb<-prepare_trait(pn_file,"Pneumonia"); coordinate_sources[[3]]<-pb$coords; B<-pb$data; rm(pb,ref); gc()
coordinates<-unique(rbindlist(coordinate_sources))
coordinates[,N_COORD:=.N,by=IV_SNP]
fwrite(coordinates,file.path(WORK_DIR,"IV_coordinates.hg19.tsv"),sep="\t")
if(any(coordinates$N_COORD!=1)) stop("Conflicting IV coordinates")
if(!setequal(coordinates$IV_SNP,iv_ids))
  stop("Unresolved IVs: ",paste(setdiff(iv_ids,coordinates$IV_SNP),collapse=";"))
windows<-coordinates[,.(CHR19,START=pmax(1,BP19-FLANK),END=BP19+FLANK,IV=IV_SNP)]
setorder(windows,CHR19,START,END)
rr<-list()
for(i in seq_len(nrow(windows))) {
  w<-windows[i]; n<-length(rr)
  if(n>0 && rr[[n]]$CHR19==w$CHR19 && w$START<=rr[[n]]$END) {
    rr[[n]]$END<-max(rr[[n]]$END,w$END)
    rr[[n]]$IVS<-paste(rr[[n]]$IVS,w$IV,sep=";")
  } else rr[[n+1]]<-data.table(CHR19=w$CHR19,START=w$START,END=w$END,IVS=w$IV)
}
regions<-rbindlist(rr)
regions[,REGION:=sprintf("R%03d",.I)]
regions[,HG19_region:=paste0("chr",CHR19,":",START,"-",END)]
fwrite(regions,file.path(WORK_DIR,"regions.hg19.tsv"),sep="\t")
A<-A[is.finite(N_EFF)&N_EFF>0.9*103614]
for(i in seq_len(nrow(regions))) {
  r<-regions[i]
  a<-A[CHR19==r$CHR19 & BP19>=r$START & BP19<=r$END]
  b<-B[CHR19==r$CHR19 & BP19>=r$START & BP19<=r$END]
  a<-a[SNP %chin% intersect(a$SNP,b$SNP)]; setorder(a,BP19,SNP)
  b<-b[match(a$SNP,b$SNP)]
  stopifnot(identical(a$SNP,b$SNP),all(a$BP19==b$BP19),all(a$EA==b$EA),all(a$OA==b$OA))
  d<-file.path(region_dir,r$REGION); dir.create(d,showWarnings=FALSE)
  fwrite(a,file.path(d,"Psoriasis.tsv"),sep="\t",na="NA")
  fwrite(b,file.path(d,"Pneumonia.tsv"),sep="\t",na="NA")
  fwrite(a[,.(SNP)],file.path(d,"snps.txt"),col.names=FALSE)
  fwrite(a[,.(SNP,EA)],file.path(d,"a1.tsv"),sep="\t",col.names=FALSE)
}
rm(A,B); gc()
all_status<-list(); all_pairs<-list()
save_tables<-function() {
  ss<-rbindlist(all_status,fill=TRUE)
  fwrite(ss,file.path(WORK_DIR,"analysis_status.tsv"),sep="\t",na="NA")
  s1<-merge(regions[,.(REGION,HG19_region,MR_index_IVs=IVS)],
    ss[,.(REGION,Psoriasis_CS_count=N_CS_A,Pneumonia_CS_count=N_CS_B)],by="REGION",all.x=TRUE,sort=FALSE)
  s1<-s1[match(regions$REGION,s1$REGION)]; s1[,REGION:=NULL]
  fwrite(s1,file.path(table_dir,"Table_S1.csv"),bom=TRUE,na="NA")
  fields<-c("HG19_region","MR_index_IVs","Psoriasis_signal_lead_SNP","Pneumonia_signal_lead_SNP",
    "PP.H0","PP.H1","PP.H2","PP.H3","PP.H4","Psoriasis_signal_index","Pneumonia_signal_index")
  s2<-if(length(all_pairs)) rbindlist(all_pairs,fill=TRUE) else
    as.data.table(setNames(replicate(length(fields),character(),simplify=FALSE),fields))
  fwrite(s2[,..fields],file.path(table_dir,"Table_S2.csv"),bom=TRUE,na="NA")
}
for(i in seq_len(nrow(regions))) {
  r<-regions[i]; message("Region ",i,"/",nrow(regions),": ",r$HG19_region)
  input<-file.path(region_dir,r$REGION); out<-file.path(result_dir,r$REGION)
  dir.create(out,showWarnings=FALSE)
  status<-list(REGION=r$REGION,STATUS="STARTED",N_CS_A=NA_integer_,N_CS_B=NA_integer_)
  tryCatch({
    a<-fread(file.path(input,"Psoriasis.tsv")); b<-fread(file.path(input,"Pneumonia.tsv"))
    if(nrow(a)<2) {
      status$STATUS<-"INSUFFICIENT_SNPS"
    } else {
      geno<-file.path(out,"genotypes")
      command(PLINK,c("--bfile",panel,"--extract",file.path(input,"snps.txt"),
        "--a1-allele",file.path(input,"a1.tsv"),"2","1","--keep-allele-order",
        "--threads",THREADS,"--memory",PLINK_MEMORY_MB,"--make-bed","--out",geno))
      command(PLINK,c("--bfile",geno,"--a1-allele",file.path(input,"a1.tsv"),"2","1",
        "--keep-allele-order","--recode","A","--out",geno))
      bim<-read_bim(paste0(geno,".bim")); raw<-fread(paste0(geno,".raw"),check.names=FALSE)
      stopifnot(!anyDuplicated(bim$SNP),setequal(a$SNP,bim$SNP),setequal(b$SNP,bim$SNP))
      a<-a[match(bim$SNP,a$SNP)]; b<-b[match(bim$SNP,b$SNP)]
      stopifnot(all(a$EA==bim$A1),all(b$EA==bim$A1),all(a$OA==bim$A2),all(b$OA==bim$A2),
        identical(names(raw)[-(1:6)],paste(bim$SNP,bim$A1,sep="_")))
      G<-as.matrix(raw[,-(1:6),with=FALSE]); storage.mode(G)<-"double"
      means<-colMeans(G,na.rm=TRUE); ok<-is.finite(means)
      G<-G[,ok,drop=FALSE]; a<-a[ok]; b<-b[ok]; means<-means[ok]
      for(j in seq_len(ncol(G))) G[is.na(G[,j]),j]<-means[j]
      sd_g<-apply(G,2,sd); ok<-is.finite(sd_g)&sd_g>0
      G<-G[,ok,drop=FALSE]; a<-a[ok]; b<-b[ok]
      if(ncol(G)<2) stop("Fewer than two polymorphic variants")
      X<-scale(G); R<-crossprod(X)/(nrow(X)-1)
      R<-(R+t(R))/2; diag(R)<-1; dimnames(R)<-list(a$SNP,a$SNP)
      stopifnot(all(is.finite(R)),all(is.finite(a$BETA/a$SE)),all(is.finite(b$BETA/b$SE)))
      n_a<-median(a$N_EFF); n_b<-4*78791*421557/(78791+421557)
      status$N_USED<-nrow(a); status$N_REFERENCE<-nrow(G)
      status$N_SUSIE_A<-n_a; status$N_SUSIE_B<-n_b
      fwrite(a,file.path(out,"Psoriasis.input.tsv"),sep="\t",na="NA")
      fwrite(b,file.path(out,"Pneumonia.input.tsv"),sep="\t",na="NA")
      fit_one<-function(d,n,trait) {
        fit<-susieR::susie_rss(z=d$BETA/d$SE,R=R,n=n,L=min(10L,nrow(d)),
          coverage=0.95,min_abs_corr=0.5,estimate_residual_variance=FALSE,max_iter=1000,verbose=FALSE)
        fit<-coloc::annotate_susie(fit,d$SNP,R)
        saveRDS(fit,file.path(out,paste0(trait,".susie.rds")))
        if(length(fit$sets$cs)) {
          members<-rbindlist(lapply(names(fit$sets$cs),function(k) {
            j<-fit$sets$cs[[k]]; data.table(CS=k,SNP=d$SNP[j],PIP=fit$pip[j])
          }))
          fwrite(members,file.path(out,paste0(trait,".credible_sets.tsv")),sep="\t")
        }
        fit
      }
      fa<-fit_one(a,n_a,"Psoriasis"); fb<-fit_one(b,n_b,"Pneumonia")
      status$CONVERGED_A<-isTRUE(fa$converged); status$CONVERGED_B<-isTRUE(fb$converged)
      status$N_CS_A<-length(fa$sets$cs); status$N_CS_B<-length(fb$sets$cs)
      if(!isTRUE(fa$converged)||!isTRUE(fb$converged)) status$STATUS<-"SUSIE_NOT_CONVERGED" else
      if(!length(fa$sets$cs)||!length(fb$sets$cs)) status$STATUS<-"NO_CREDIBLE_SET_IN_ONE_OR_BOTH" else {
        co<-coloc::coloc.susie(fa,fb,p1=1e-4,p2=1e-4,p12=1e-5)
        saveRDS(co,file.path(out,"coloc.rds")); z<-as.data.table(co$summary)
        if(!nrow(z)) status$STATUS<-"NO_COLOC_COMPARISONS" else {
          fwrite(z,file.path(out,"coloc.raw_summary.tsv"),sep="\t")
          pairs<-z[,.(HG19_region=r$HG19_region,MR_index_IVs=r$IVS,
            Psoriasis_signal_lead_SNP=hit1,Pneumonia_signal_lead_SNP=hit2,
            PP.H0=PP.H0.abf,PP.H1=PP.H1.abf,PP.H2=PP.H2.abf,PP.H3=PP.H3.abf,PP.H4=PP.H4.abf,
            Psoriasis_signal_index=idx1,Pneumonia_signal_index=idx2)]
          all_pairs[[length(all_pairs)+1]]<-pairs
          status$STATUS<-"COMPLETED"; status$N_SUPPORTED_PAIRS<-sum(z$PP.H4.abf>=0.8)
        }
      }
    }
  },error=function(e) {
    status$STATUS<<-"ERROR"; status$ERROR<<-conditionMessage(e)
    message("Region failed: ",r$REGION,"; ",conditionMessage(e))
  })
  all_status[[length(all_status)+1]]<-as.data.table(status)
  save_tables(); gc()
}
message("Finished. CSV tables: ",table_dir)
ss<-rbindlist(all_status,fill=TRUE); print(ss[,.N,by=STATUS])
if(any(ss$STATUS %chin% c("ERROR","SUSIE_NOT_CONVERGED"))) {
  message("Inspect analysis_status.tsv: failed/nonconverged regions are not negative results.")
  quit(status=2)
}
print(p_combined)
