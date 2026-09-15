library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

rmatstable <- fread("~/path/to/your/project/final_output_rMATS/SE.MATS.JC.txt")[,1:11]
vttable <- fread("~/path/to/your/project/Final_BIGEX_Table.txt")

#I need to add a metric for Length. I need to bear in mind that rmats is in base 10

rmatstable$Length <- rmatstable$exonEnd - rmatstable$exonStart_0base

rmatstable <- subset(rmatstable, Length>10)  

#My boy rmats is in base 0. Eu não papo grupos. Vou converter tudo para base 1
rmatstable$exonStart_0base <- rmatstable$exonStart_0base+1
rmatstable$exonEnd <- rmatstable$exonEnd+1
rmatstable$upstreamES <- rmatstable$upstreamES+1
rmatstable$upstreamEE <- rmatstable$upstreamEE+1
rmatstable$downstreamES <- rmatstable$downstreamES+1
rmatstable$downstreamEE <- rmatstable$downstreamEE+1

#Lets add coordenates we can infer from the provided ones
rmatstable$UPS.INT.ST.ST <- (rmatstable$upstreamEE + 1)
rmatstable$UPS.INT.EN.EN <- (rmatstable$exonStart_0base - 1)
rmatstable$DOW.INT.ST.ST <- (rmatstable$exonEnd + 1)
rmatstable$DOW.INT.EN.EN <- (rmatstable$downstreamES - 1)

#Correct till here

for (row in 1:nrow(rmatstable)){
  #Defining the alternative exon itself
  if (rmatstable[row, "exonEnd"]-rmatstable[row, "exonStart_0base"]<=100){
    AE_length_upstream <- floor((rmatstable[row, "exonEnd"]-rmatstable[row, "exonStart_0base"]-1)/2)
    AE_length_downstream <- ceiling((rmatstable[row, "exonEnd"]-rmatstable[row, "exonStart_0base"]-1)/2)
    if (rmatstable[row, "strand"]=="+"){
      rmatstable[row,"AER.ST"] <- rmatstable[row,"exonStart_0base"] + AE_length_upstream
      rmatstable[row,"AER.EN"] <- rmatstable[row,"exonEnd"] - AE_length_downstream
    }
    else{
      rmatstable[row,"AER.ST"] <- rmatstable[row,"exonStart_0base"] + AE_length_downstream
      rmatstable[row,"AER.EN"] <- rmatstable[row,"exonEnd"] - AE_length_upstream
    }
  }
  else{
    rmatstable[row,"AER.ST"] <- rmatstable[row,"exonStart_0base"] + 50
    rmatstable[row,"AER.EN"] <- rmatstable[row,"exonEnd"] - 50
  }
  #Defining the upstream intron
  if (rmatstable[row, "UPS.INT.EN.EN"]-rmatstable[row, "UPS.INT.ST.ST"]<=400){
    UI_length_upstream <- floor((rmatstable[row, "UPS.INT.EN.EN"]-rmatstable[row, "UPS.INT.ST.ST"]-1)/2)
    UI_length_downstream <- ceiling((rmatstable[row, "UPS.INT.EN.EN"]-rmatstable[row, "UPS.INT.ST.ST"]-1)/2)
    if (rmatstable[row, "strand"]=="+"){
      rmatstable[row,"UPS.INT.ST.EN"] <- rmatstable[row,"UPS.INT.ST.ST"] + UI_length_upstream
      rmatstable[row,"UPS.INT.EN.ST"] <- rmatstable[row,"UPS.INT.EN.EN"] - UI_length_downstream
    }
    else{
      rmatstable[row,"UPS.INT.ST.EN"] <- rmatstable[row,"UPS.INT.ST.ST"] + UI_length_downstream
      rmatstable[row,"UPS.INT.EN.ST"] <- rmatstable[row,"UPS.INT.EN.EN"] - UI_length_upstream
    }
  }
  else{
    rmatstable[row,"UPS.INT.ST.EN"] <- rmatstable[row, "UPS.INT.ST.ST"]+200
    rmatstable[row,"UPS.INT.EN.ST"] <- rmatstable[row, "UPS.INT.EN.EN"]-200
  }#Defining the downstream intron
  if (rmatstable[row, "DOW.INT.EN.EN"]-rmatstable[row, "DOW.INT.ST.ST"]<=400){
    DI_length_upstream <- floor((rmatstable[row, "DOW.INT.EN.EN"]-rmatstable[row, "DOW.INT.ST.ST"]-1)/2)
    DI_length_downstream <- ceiling((rmatstable[row, "DOW.INT.EN.EN"]-rmatstable[row, "DOW.INT.ST.ST"]-1)/2)
    if (rmatstable[row, "strand"]=="+"){
      rmatstable[row,"DOW.INT.ST.EN"] <- rmatstable[row,"DOW.INT.ST.ST"] + DI_length_upstream
      rmatstable[row,"DOW.INT.EN.ST"] <- rmatstable[row,"DOW.INT.EN.EN"] - DI_length_downstream
    }
    else{
      rmatstable[row,"DOW.INT.ST.EN"] <- rmatstable[row,"DOW.INT.ST.ST"] + DI_length_downstream
      rmatstable[row,"DOW.INT.EN.ST"] <- rmatstable[row,"DOW.INT.EN.EN"] - DI_length_upstream
    }
  }
  else{
    rmatstable[row,"DOW.INT.ST.EN"] <- rmatstable[row, "DOW.INT.ST.ST"]+200
    rmatstable[row,"DOW.INT.EN.ST"] <- rmatstable[row, "DOW.INT.EN.EN"]-200
  }
  #THIS IS DIFFERENT THAN WITH VT
  #Defining the upstream exon. We just need to change if exon is bigger than 50
  if (rmatstable[row, "upstreamEE"]-rmatstable[row, "upstreamES"]>=50){
    rmatstable[row,"upstreamES"] <- rmatstable[row,"upstreamEE"] - 50
  }
  #Defining the downstream exon. We just need to change if exon is bigger than 50
  if (rmatstable[row, "downstreamEE"]-rmatstable[row, "downstreamES"]>=50){
    rmatstable[row,"downstreamEE"] <- rmatstable[row,"downstreamES"] + 50
  }
}

saveRDS(rmatstable, "~/path/to/your/project/rmatstable_as.RDS")

#Now we transform the rmatstable into a vtttols table format
#First we remove GENEID and Length
rmatstable <- readRDS("~/path/to/your/project/rmatstable_as.RDS")[,-c(2,12)]


rmatstable <- rmatstable[,c(2,1,3,4,7,8,11,17,18,12,5,15,16,6,13,19,20,14,9,10)]
newnames <- colnames(vttable)
colnames(rmatstable) <- newnames

write.csv(rmatstable, "~/path/to/your/project/Final_AS_Table_Rmats.txt", row.names = F)
