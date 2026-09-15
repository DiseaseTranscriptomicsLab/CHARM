library(data.table)
library(tidyr)
library(dplyr)
library(stringr)


vttable <- fread("~/path/to/your/project/INCLUSION_LEVELS_FULL-hg38-25-v251.tab.gz")[,1:6]

##To facilitate integration into my python script, I need to separate both of these events into different dataframes

vttable_INT <- subset(vttable, grepl("^HsaINT", EVENT))
vttable_EX <- subset(vttable, grepl("^HsaEX", EVENT))


#Lets start with Intron Retention Dataframe Cleanup. It is simpler.
#With a simple cleanup by intro length. Very small introns are probably mistakes.
vttable_INT <- subset(vttable_INT, LENGTH>30)

vttable_INT <- separate(vttable_INT, COORD, into = c("CHROM", "COORD.ST", "COORD.EN"), sep = "[:-]")
vttable_INT <- extract(vttable_INT, FullCO, into = c("EX1.ST", "EX1.EN", "EX2.ST", "EX2.EN", "STRAND"), regex = ".*:(\\d+)-(\\d+)=(\\d+)-(\\d+):(.*)")

#We need to make some changes in accordance to the sense or antisense


for (row in 1:nrow(vttable_INT)){
  if (vttable_INT[row, "STRAND"]=="-"){
    vttable_INT[row,] <- vttable_INT[row, c(1:6,9:10,7:8,11)]
    colnames(vttable_INT[row,]) <- colnames(vttable_INT[row, ])[c(1:6,9:10,7:8,11)]
  }
}

#Now we define our regions of interest
vttable_INT_regions <- vttable_INT
vttable_INT_regions <- mutate_at(vttable_INT_regions, vars(4:10), as.numeric)

for (row in 1:nrow(vttable_INT_regions)){
  if (vttable_INT_regions[row, "COORD.EN"]-vttable_INT_regions[row, "COORD.ST"]<=400){
    INT_length_upstream <- floor((vttable_INT_regions[row, "COORD.EN"]-vttable_INT_regions[row, "COORD.ST"]-1)/2)
    INT_length_downstream <- ceiling((vttable_INT_regions[row, "COORD.EN"]-vttable_INT_regions[row, "COORD.ST"]-1)/2)
    if (vttable_INT_regions[row, "STRAND"]=="+"){
      vttable_INT_regions[row,"IN.ST"] <- vttable_INT_regions[row,"COORD.ST"] + INT_length_upstream
      vttable_INT_regions[row,"IN.EN"] <- vttable_INT_regions[row,"COORD.EN"] - INT_length_downstream
    }
    else{
      vttable_INT_regions[row,"IN.ST"] <- vttable_INT_regions[row,"COORD.ST"] + INT_length_downstream
      vttable_INT_regions[row,"IN.EN"] <- vttable_INT_regions[row,"COORD.EN"] - INT_length_upstream
    }
  }
  else{
    vttable_INT_regions[row,"IN.ST"] <- vttable_INT_regions[row,"COORD.ST"] + 200 #Just in case you forget again, the region is given by
    vttable_INT_regions[row,"IN.EN"] <- vttable_INT_regions[row,"COORD.EN"] - 200 #COORD.ST  --- IN.ST
  }
  vttable_INT_regions[row,"UPS.ST"] <- vttable_INT_regions[row,"EX1.EN"] - 50
  vttable_INT_regions[row,"UPS.EN"] <- vttable_INT_regions[row,"EX1.EN"]
  vttable_INT_regions[row,"DOW.ST"] <- vttable_INT_regions[row,"EX2.ST"]
  vttable_INT_regions[row,"DOW.EN"] <- vttable_INT_regions[row,"EX2.ST"] + 50
}

vttable_IR_regions <- vttable_INT_regions


saveRDS(vttable_INT_regions, "~/path/to/your/project/vttable_IR_regions.RDS")



#Lets move on to the Exon Skipping Cleanup
#Not sure if its correct or not, but I decided to remove exons smaller than 10bp
vttable_EX <- subset(vttable_EX, LENGTH>10)

vttable_EX <- separate(vttable_EX, COORD, into = c("CHROM", "COORD.ST", "COORD.EN"), sep = "[:-]")


# Assuming your column is named 'your_column' and your dataframe is named 'your_dataframe'
vttable_EX <- vttable_EX %>%
  mutate(FullCO = gsub(".*:", "", FullCO),  # Discard everything up to ":"
  )

# Separate the modified column into new columns
vttable_EX <- separate(vttable_EX, FullCO, into = c("EX1", "AE", "EX2"), sep = ",", convert = TRUE)





#In this scenario, we define coordinates in a way that maximizes intron size. 
vttable_EX_bigintron <- vttable_EX
vttable_EX_bigintron$EX1 <- sapply(strsplit(as.character(vttable_EX_bigintron$EX1), "\\+"), function(x) head(x, 1))
vttable_EX_bigintron$EX2 <- sapply(strsplit(as.character(vttable_EX_bigintron$EX2), "\\+"), function(x) tail(x, 1))

vttable_EX_bigintron <- separate(vttable_EX_bigintron, AE, into = c("AE.ST", "AE.EN"), sep = "-", convert = TRUE)

vttable_EX_bigintron$AE.ST <- sapply(strsplit(as.character(vttable_EX_bigintron$AE.ST), "\\+"), function(x) tail(x, 1))
vttable_EX_bigintron$AE.EN <- sapply(strsplit(as.character(vttable_EX_bigintron$AE.EN), "\\+"), function(x) head(x, 1))

vttable_EX_bigintron_regions <- vttable_EX_bigintron
vttable_EX_bigintron_regions <- mutate_at(vttable_EX_bigintron_regions, vars(4:10), as.numeric)


for (row in 1:nrow(vttable_EX_bigintron_regions)){
  if (vttable_EX_bigintron_regions[row, "EX2"]<vttable_EX_bigintron_regions[row, "EX1"]){
    vttable_EX_bigintron_regions[row,] <- vttable_EX_bigintron_regions[row, c(1:6,10,8:9,7)]
    colnames(vttable_EX_bigintron_regions[row,]) <- colnames(vttable_EX_bigintron_regions[row, ])[c(1:6,10,8:9,7)]
    vttable_EX_bigintron_regions[row,"STRAND"] <- "-"
  }
  else{
    vttable_EX_bigintron_regions[row,"STRAND"] <- "+"
  }
}


vttable_EX_bigintron_regions$UPS.INT.ST.ST <- (vttable_EX_bigintron_regions$EX1 + 1)
vttable_EX_bigintron_regions$UPS.INT.EN.EN <- (vttable_EX_bigintron_regions$COORD.ST - 1)
vttable_EX_bigintron_regions$DOW.INT.ST.ST <- (vttable_EX_bigintron_regions$COORD.EN + 1)
vttable_EX_bigintron_regions$DOW.INT.EN.EN <- (vttable_EX_bigintron_regions$EX2 - 1)


for (row in 1:nrow(vttable_EX_bigintron_regions)){
  #Defining the alternative exon itself
  if (vttable_EX_bigintron_regions[row, "AE.EN"]-vttable_EX_bigintron_regions[row, "AE.ST"]<=100){
    AE_length_upstream <- floor((vttable_EX_bigintron_regions[row, "AE.EN"]-vttable_EX_bigintron_regions[row, "AE.ST"]-1)/2)
    AE_length_downstream <- ceiling((vttable_EX_bigintron_regions[row, "AE.EN"]-vttable_EX_bigintron_regions[row, "AE.ST"]-1)/2)
    if (vttable_EX_bigintron_regions[row, "STRAND"]=="+"){
      vttable_EX_bigintron_regions[row,"AER.ST"] <- vttable_EX_bigintron_regions[row,"AE.ST"] + AE_length_upstream
      vttable_EX_bigintron_regions[row,"AER.EN"] <- vttable_EX_bigintron_regions[row,"AE.EN"] - AE_length_downstream
    }
    else{
      vttable_EX_bigintron_regions[row,"AER.ST"] <- vttable_EX_bigintron_regions[row,"AE.ST"] + AE_length_downstream
      vttable_EX_bigintron_regions[row,"AER.EN"] <- vttable_EX_bigintron_regions[row,"AE.EN"] - AE_length_upstream
    }
  }
  else{
    vttable_EX_bigintron_regions[row,"AER.ST"] <- vttable_EX_bigintron_regions[row,"AE.ST"] + 50
    vttable_EX_bigintron_regions[row,"AER.EN"] <- vttable_EX_bigintron_regions[row,"AE.EN"] - 50
  }
  #Defining the upstream intron
  if (vttable_EX_bigintron_regions[row, "UPS.INT.EN.EN"]-vttable_EX_bigintron_regions[row, "UPS.INT.ST.ST"]<=400){
    UI_length_upstream <- floor((vttable_EX_bigintron_regions[row, "UPS.INT.EN.EN"]-vttable_EX_bigintron_regions[row, "UPS.INT.ST.ST"]-1)/2)
    UI_length_downstream <- ceiling((vttable_EX_bigintron_regions[row, "UPS.INT.EN.EN"]-vttable_EX_bigintron_regions[row, "UPS.INT.ST.ST"]-1)/2)
    if (vttable_EX_bigintron_regions[row, "STRAND"]=="+"){
      vttable_EX_bigintron_regions[row,"UPS.INT.ST.EN"] <- vttable_EX_bigintron_regions[row,"UPS.INT.ST.ST"] + UI_length_upstream
      vttable_EX_bigintron_regions[row,"UPS.INT.EN.ST"] <- vttable_EX_bigintron_regions[row,"UPS.INT.EN.EN"] - UI_length_downstream
    }
    else{
      vttable_EX_bigintron_regions[row,"UPS.INT.ST.EN"] <- vttable_EX_bigintron_regions[row,"UPS.INT.ST.ST"] + UI_length_downstream
      vttable_EX_bigintron_regions[row,"UPS.INT.EN.ST"] <- vttable_EX_bigintron_regions[row,"UPS.INT.EN.EN"] - UI_length_upstream
    }
  }
  else{
    vttable_EX_bigintron_regions[row,"UPS.INT.ST.EN"] <- vttable_EX_bigintron_regions[row, "UPS.INT.ST.ST"]+200
    vttable_EX_bigintron_regions[row,"UPS.INT.EN.ST"] <- vttable_EX_bigintron_regions[row, "UPS.INT.EN.EN"]-200
  }
  #Defining the downstream intron
  if (vttable_EX_bigintron_regions[row, "DOW.INT.EN.EN"]-vttable_EX_bigintron_regions[row, "DOW.INT.ST.ST"]<=400){
    DI_length_upstream <- floor((vttable_EX_bigintron_regions[row, "DOW.INT.EN.EN"]-vttable_EX_bigintron_regions[row, "DOW.INT.ST.ST"]-1)/2)
    DI_length_downstream <- ceiling((vttable_EX_bigintron_regions[row, "DOW.INT.EN.EN"]-vttable_EX_bigintron_regions[row, "DOW.INT.ST.ST"]-1)/2)
    if (vttable_EX_bigintron_regions[row, "STRAND"]=="+"){
      vttable_EX_bigintron_regions[row,"DOW.INT.ST.EN"] <- vttable_EX_bigintron_regions[row,"DOW.INT.ST.ST"] + DI_length_upstream
      vttable_EX_bigintron_regions[row,"DOW.INT.EN.ST"] <- vttable_EX_bigintron_regions[row,"DOW.INT.EN.EN"] - DI_length_downstream
    }
    else{
      vttable_EX_bigintron_regions[row,"DOW.INT.ST.EN"] <- vttable_EX_bigintron_regions[row,"DOW.INT.ST.ST"] + DI_length_upstream
      vttable_EX_bigintron_regions[row,"DOW.INT.EN.ST"] <- vttable_EX_bigintron_regions[row,"DOW.INT.EN.EN"] - DI_length_downstream
    }
  }
  else{
    vttable_EX_bigintron_regions[row,"DOW.INT.ST.EN"] <- vttable_EX_bigintron_regions[row, "DOW.INT.ST.ST"]+200
    vttable_EX_bigintron_regions[row,"DOW.INT.EN.ST"] <- vttable_EX_bigintron_regions[row, "DOW.INT.EN.EN"]-200
  }
  #Defining the upstream and downstream exons
  vttable_EX_bigintron_regions[row,"UPS.EX.ST"] <- vttable_EX_bigintron_regions[row,"EX1"] - 50
  vttable_EX_bigintron_regions[row,"DOW.EX.EN"] <- vttable_EX_bigintron_regions[row,"EX2"] + 50
}

saveRDS(vttable_EX_bigintron_regions, "~/path/to/your/project/vttable_EX_bigintron_regions.RDS")

#In this scenario, we instead define coordinates in a way that maximizes the alternative exon's own size
vttable_EX_bigexon <- vttable_EX
vttable_EX_bigexon$EX1 <- sapply(strsplit(as.character(vttable_EX_bigexon$EX1), "\\+"), function(x) tail(x, 1))
vttable_EX_bigexon$EX2 <- sapply(strsplit(as.character(vttable_EX_bigexon$EX2), "\\+"), function(x) head(x, 1))

vttable_EX_bigexon <- separate(vttable_EX_bigexon, AE, into = c("AE.ST", "AE.EN"), sep = "-", convert = TRUE)

vttable_EX_bigexon$AE.ST <- sapply(strsplit(as.character(vttable_EX_bigexon$AE.ST), "\\+"), function(x) head(x, 1))
vttable_EX_bigexon$AE.EN <- sapply(strsplit(as.character(vttable_EX_bigexon$AE.EN), "\\+"), function(x) tail(x, 1))

vttable_EX_bigexon_regions <- vttable_EX_bigexon
vttable_EX_bigexon_regions <- mutate_at(vttable_EX_bigexon_regions, vars(4:10), as.numeric)


for (row in 1:nrow(vttable_EX_bigexon_regions)){
  if (vttable_EX_bigexon_regions[row, "EX2"]<vttable_EX_bigexon_regions[row, "EX1"]){
    vttable_EX_bigexon_regions[row,] <- vttable_EX_bigexon_regions[row, c(1:6,10,8:9,7)]
    colnames(vttable_EX_bigexon_regions[row,]) <- colnames(vttable_EX_bigexon_regions[row, ])[c(1:6,10,8:9,7)]
    vttable_EX_bigexon_regions[row,"STRAND"] <- "-"
  }
  else{
    vttable_EX_bigexon_regions[row,"STRAND"] <- "+"
  }
}


vttable_EX_bigexon_regions$UPS.INT.ST.ST <- (vttable_EX_bigexon_regions$EX1 + 1)
vttable_EX_bigexon_regions$UPS.INT.EN.EN <- (vttable_EX_bigexon_regions$COORD.ST - 1)
vttable_EX_bigexon_regions$DOW.INT.ST.ST <- (vttable_EX_bigexon_regions$COORD.EN + 1)
vttable_EX_bigexon_regions$DOW.INT.EN.EN <- (vttable_EX_bigexon_regions$EX2 - 1)


for (row in 1:nrow(vttable_EX_bigexon_regions)){
  #Defining the alternative exon itself
  if (vttable_EX_bigexon_regions[row, "AE.EN"]-vttable_EX_bigexon_regions[row, "AE.ST"]<=100){
    AE_length_upstream <- floor((vttable_EX_bigexon_regions[row, "AE.EN"]-vttable_EX_bigexon_regions[row, "AE.ST"]-1)/2)
    AE_length_downstream <- ceiling((vttable_EX_bigexon_regions[row, "AE.EN"]-vttable_EX_bigexon_regions[row, "AE.ST"]-1)/2)
    if (vttable_EX_bigexon_regions[row, "STRAND"]=="+"){
      vttable_EX_bigexon_regions[row,"AER.ST"] <- vttable_EX_bigexon_regions[row,"AE.ST"] + AE_length_upstream
      vttable_EX_bigexon_regions[row,"AER.EN"] <- vttable_EX_bigexon_regions[row,"AE.EN"] - AE_length_downstream
    }
    else{
      vttable_EX_bigexon_regions[row,"AER.ST"] <- vttable_EX_bigexon_regions[row,"AE.ST"] + AE_length_downstream
      vttable_EX_bigexon_regions[row,"AER.EN"] <- vttable_EX_bigexon_regions[row,"AE.EN"] - AE_length_upstream
    }
  }
  else{
    vttable_EX_bigexon_regions[row,"AER.ST"] <- vttable_EX_bigexon_regions[row,"AE.ST"] + 50
    vttable_EX_bigexon_regions[row,"AER.EN"] <- vttable_EX_bigexon_regions[row,"AE.EN"] - 50
  }
  #Defining the upstream intron
  if (vttable_EX_bigexon_regions[row, "UPS.INT.EN.EN"]-vttable_EX_bigexon_regions[row, "UPS.INT.ST.ST"]<=400){
    UI_length_upstream <- floor((vttable_EX_bigexon_regions[row, "UPS.INT.EN.EN"]-vttable_EX_bigexon_regions[row, "UPS.INT.ST.ST"]-1)/2)
    UI_length_downstream <- ceiling((vttable_EX_bigexon_regions[row, "UPS.INT.EN.EN"]-vttable_EX_bigexon_regions[row, "UPS.INT.ST.ST"]-1)/2)
    if (vttable_EX_bigexon_regions[row, "STRAND"]=="+"){
      vttable_EX_bigexon_regions[row,"UPS.INT.ST.EN"] <- vttable_EX_bigexon_regions[row,"UPS.INT.ST.ST"] + UI_length_upstream
      vttable_EX_bigexon_regions[row,"UPS.INT.EN.ST"] <- vttable_EX_bigexon_regions[row,"UPS.INT.EN.EN"] - UI_length_downstream
    }
    else{
      vttable_EX_bigexon_regions[row,"UPS.INT.ST.EN"] <- vttable_EX_bigexon_regions[row,"UPS.INT.ST.ST"] + UI_length_downstream
      vttable_EX_bigexon_regions[row,"UPS.INT.EN.ST"] <- vttable_EX_bigexon_regions[row,"UPS.INT.EN.EN"] - UI_length_upstream
    }
  }
  else{
    vttable_EX_bigexon_regions[row,"UPS.INT.ST.EN"] <- vttable_EX_bigexon_regions[row, "UPS.INT.ST.ST"]+200
    vttable_EX_bigexon_regions[row,"UPS.INT.EN.ST"] <- vttable_EX_bigexon_regions[row, "UPS.INT.EN.EN"]-200
  }
  #Defining the downstream intron
  if (vttable_EX_bigexon_regions[row, "DOW.INT.EN.EN"]-vttable_EX_bigexon_regions[row, "DOW.INT.ST.ST"]<=400){
    DI_length_upstream <- floor((vttable_EX_bigexon_regions[row, "DOW.INT.EN.EN"]-vttable_EX_bigexon_regions[row, "DOW.INT.ST.ST"]-1)/2)
    DI_length_downstream <- ceiling((vttable_EX_bigexon_regions[row, "DOW.INT.EN.EN"]-vttable_EX_bigexon_regions[row, "DOW.INT.ST.ST"]-1)/2)
    if (vttable_EX_bigexon_regions[row, "STRAND"]=="+"){
      vttable_EX_bigexon_regions[row,"DOW.INT.ST.EN"] <- vttable_EX_bigexon_regions[row,"DOW.INT.ST.ST"] + DI_length_upstream
      vttable_EX_bigexon_regions[row,"DOW.INT.EN.ST"] <- vttable_EX_bigexon_regions[row,"DOW.INT.EN.EN"] - DI_length_downstream
    }
    else{
      vttable_EX_bigexon_regions[row,"DOW.INT.ST.EN"] <- vttable_EX_bigexon_regions[row,"DOW.INT.ST.ST"] + DI_length_upstream
      vttable_EX_bigexon_regions[row,"DOW.INT.EN.ST"] <- vttable_EX_bigexon_regions[row,"DOW.INT.EN.EN"] - DI_length_downstream
    }
  }
  else{
    vttable_EX_bigexon_regions[row,"DOW.INT.ST.EN"] <- vttable_EX_bigexon_regions[row, "DOW.INT.ST.ST"]+200
    vttable_EX_bigexon_regions[row,"DOW.INT.EN.ST"] <- vttable_EX_bigexon_regions[row, "DOW.INT.EN.EN"]-200
  }
  #Defining the upstream and downstream exons
  vttable_EX_bigexon_regions[row,"UPS.EX.ST"] <- vttable_EX_bigexon_regions[row,"EX1"] - 50
  vttable_EX_bigexon_regions[row,"DOW.EX.EN"] <- vttable_EX_bigexon_regions[row,"EX2"] + 50
}

saveRDS(vttable_EX_bigexon_regions, "~/path/to/your/project/vttable_EX_bigexon_regions.RDS")


#Some minor processing before the python pipeline
#To make table understanding better, I will order the table by genomic coordinate
#

vttable_INT_final <- vttable_IR_regions[,c(1:3,11,6,15:16,4,13:14,5,17:18)] #This one is correct.
vttable_EX_bigintron_final <- vttable_EX_bigintron_regions[,c(1:3,11,22,7,12,18,19,13,8,16,17,9,14,20,21,15,10,23)]
vttable_EX_bigexon_final <- vttable_EX_bigexon_regions[,c(1:3,11,22,7,12,18,19,13,8,16,17,9,14,20,21,15,10,23)]
# # 
# # 
# # #Saving all of this
write.csv(vttable_INT_final, "~/path/to/your/project/Final_IR_Table.txt", row.names = F)
write.csv(vttable_EX_bigintron_final, "~/path/to/your/project/Final_SmallExon_Table.txt", row.names = F)
write.csv(vttable_EX_bigexon_final, "~/path/to/your/project/Final_BIGEX_Table.txt", row.names = F)
