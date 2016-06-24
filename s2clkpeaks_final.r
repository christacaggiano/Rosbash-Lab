###felipe data analysis 
###christa caggiano 
###15 june 
###rosbash lab 

rm(list=ls()) #removes all variables in R's memory 

setwd("/Users/Christa/Documents/Rosbash lab") #sets working directory- must be changed for each computer 

sequencing = read.csv('S2CLKseq_gsig.csv', header = T, sep = ',') #loads sequencing file
binding = read.csv('S2_clk_2REP.2.peaks..csv', header = T, sep = ',')#loads CLK peaks file
complete = read.csv('updated_drosophila_genes.csv', header = T, sep = ',')#loads dm3 genome information from UCSC genome browser

##makes a file that contains only genes where CLK activates transcription relative to ctrl 
##by 1.5 in both replicates
cutoff = subset(sequencing, sequencing$CLK.WT_R1_Signed>1.5 & sequencing$CLK.WT_R2_Signed>1.5)

#finds only the genes from the cutoff list and restricts the genome to them 
intersect = intersect(cutoff$gene, complete$gene)
common_genes = complete[complete$gene %in% intersect, ]

unique = unique(common_genes)

##converts to character for comparison 
cutoff$gene = as.character(cutoff$gene)
unique$gene = as.character(unique$gene)
unique$strand = as.character(unique$strand)
cutoff$chr = as.character(cutoff$chr)

##uses genome information to determine strandedness, calculates a start 500bp of the sequenced start coordinates
##relative to what strand the gene is on 
for(i in 1:nrow(cutoff)){
  for(j in 1:nrow(unique)){ 
    if((identical(cutoff$gene[i], unique$gene[j])) & (unique$strand[j] == '+')){
      cutoff$upstream[i]=cutoff$start[i]-500
      cutoff$strand[i]=unique$strand[j]
    }else if ((identical(cutoff$gene[i], unique$gene[j])) & (unique$strand[j] == '-')){
      cutoff$upstream[i]=cutoff$start[i]+500
      cutoff$strand[i]=unique$strand[j]
    }
  }
}

#creates empty lists to store overlapping coordinates, both for the genes CLK activates transcription on, 
#and where the peaks occur on those genes 
start_trans = list()
end_trans = list() 
start_bind = list() 
end_bind = list()
strand = list() 
chromosome = list()


#searches for peaks within the gene region where CLK activates transcription, including 500bp upstream region 
#appends these coordinates to a list 
for(i in 1:nrow(cutoff)){
  for (j in 1:nrow(binding)){
    if ((cutoff$upstream[i]<=binding$cons_start[j])&(cutoff$end[i]>=binding$cons_end[j])){
      start_trans = c(start_trans, cutoff$start[i])
      end_trans = c(end_trans, cutoff$end[i])
      start_bind = c(start_bind, binding$cons_start[j])
      end_bind = c(end_bind, binding$cons_end[j])
      chromosome = c(chromosome, cutoff$chr[i])
      strand = c(strand, cutoff$strand[i])
      
    }
  }
}

#makes lists of coordinates a dataframe/multidimensioned list 
coord = do.call(rbind, Map(data.frame, start_gene=start_trans, end_gene=end_trans, 
                           start_peak = start_bind, end_peak = end_bind, chromosome = chromosome, strand = strand ))

#makes the dataframe a CSV file and saves it to whatever directory you're working in 
write.csv(coord, file = "coordinates.csv")

