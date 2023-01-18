setwd("/Users/eilishmcmaster/Documents/ThemTria")
load("dms.RDS") # import the dms and m2 from the main rmd 

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/ea46cc026bb56cafd339f5af383c94f46e0de2dd/read_dart_counts_csv_faster_new.r?raw=TRUE")

counts2 <- read_dart_counts_csv_faster('ThemTria/dart_raw/SEQ_SNPs_counts_0_Target.csv', # import readcount data 
                                       minAlleleCount=1, 
                                       minGenotypeCount=0)

#this count file has not labelled things with nsw number but its in the meta
counts$sample_names <- counts$meta[4,]
colnames(counts$c1) <- counts$meta[4,]
colnames(counts$c2) <- counts$meta[4,]

# keep only the individuals in the dms that are assigned a species (sp)

dms <- remove.by.list(dms2, m2[!is.na(m2$sp),] %>%.$sample)


a_function <- function(dms){
  species <- unique(dms$meta$analyses[,"sp"])
  print(species)
  
  for(i in 1:length(species)){
    dmsx <- remove.by.list(dms, m2[(m2$sp %in% paste(species[i])),] %>%.$sample) %>% remove.by.maf(., 0.05)
    test2 <- count_subsetter(dmsx, counts, 0.05)
    
    tr <-  t(test2$c1)
    o <- lapply(split(tr,rownames(tr)), as.list)
    
    tr2 <-  t(test2$c2)
    o2 <- lapply(split(tr2,rownames(tr2)), as.list)
    
    nn <- mergeLists_internal(o, o2)
    
    minor <- lapply(nn, sapply, function(x) min(x)/sum(x))
    a <- do.call(rbind, minor) #make matrix
    major <- lapply(nn, sapply, function(x) max(x)/sum(x))
    b <- do.call(rbind, major) #make matrix
    c <- cbind(a,b)
    
    # only for six samples at a time
    par(mfrow = c(4, 5), mai=c(0.5,0.2,0.2,0.2))  # Set up a 2 x 2 plotting space
    
    # Create the loop.vector (all the columns)
    loop.vector <- 1:nrow(c)
    z <- paste(unique(dms$meta$analyses[,"sp"])[i])
    print(z)
    
    hist(c, main=z, xlab="", ylab="", breaks=50)
    
    for (i in loop.vector) { # Loop over loop.vector
      
      # store data in row.i as x
      x <- c[i,]
      if(sum(x, na.rm=TRUE)>0){ # skip empties
        # Plot histogram of x
        hist(x,breaks=50,
             main = paste(rownames(c)[i]),
             xlab = "",#"MAF reads/ total reads",
             ylab="",
             xlim = c(0, 1))
      }
    }
    
  }
}


# plot histograms of the readcound data for all samples at the same time
# run by species dms
# minor allele frequency 0.05
# 
# dms_sg <-  remove.by.list(dms, m2[(m2$sp %in% "sg"),] %>%.$sample) 
# dms_aus <-  remove.by.list(dms, m2[(m2$sp %in% "aus"),] %>%.$sample) 
# dms_ng <-  remove.by.list(dms, m2[(m2$sp %in% "ng"),] %>%.$sample) 
# dms_sole <-  remove.by.list(dms, m2[(m2$sp %in% "sole"),] %>%.$sample) 
# dms_nole <-  remove.by.list(dms, m2[(m2$sp %in% "nole"),] %>%.$sample) 

# par(mfrow = c(2, 3),mai=c(0.5,0.5,0.2,0.2)) 
# doitall(dms_sg, counts2, 0.05, "sg")
# doitall(dms_aus, counts2, 0.05, "aus")
# doitall(dms_ng, counts2, 0.05, "ng")
# doitall(dms_nole, counts2, 0.05, "sole")
# doitall(dms_sole, counts2, 0.05, "nole")

pdf(file="all.pdf")
a_function(dms)
dev.off()
