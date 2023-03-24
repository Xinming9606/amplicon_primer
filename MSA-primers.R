rm(list=ls());gc()

#####
#Mess around here
#####

#adress of alignment
args <- commandArgs(TRUE)
ALNS <- args[1]
 
#global options
primer_len <- as.numeric(args[3])     # length for finding potential primers
amplion_max_len <- args[4] # max amplicon length
amplion_min_len <- args[5] # min amplicon length
div_cut <- args[6]           # initial diversity cutoff to start out with
GC_tol <- args[7]        # max difference in GC%



#####
#stop messing around
#####

#libraries and plot settings
if (!requireNamespace("seqinr", quietly = TRUE)) {
  install.packages("seqinr")
}
library("seqinr")
if (!requireNamespace("zoo", quietly = TRUE)) {
  install.packages("zoo")
}
library("zoo")
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library("pheatmap")

par(mar = c(2.5, 2.5, 1.8,.5), family="serif", mfrow=c(1,1),mgp = c(1.3, 0.3, 0), font.lab=2)

#####
#Functions
#####
pos=478
#Formatter function for printing alignments
getAln=function(pos) {
  #work out if there are wobbles
  #looks at if the table() returns more than one element on each column
  subs=paste(apply(DNA_mat[,c(pos:(pos+primer_len))], 2,
                   function(x) length(table(x))-1 ),collapse="")
  
  #fishing out the consensus sequence of the primer
  consPrimer=CONSENSUS[pos:(pos+primer_len)]
  
  #go through DNA_mat and work out how many sequences diverge from the consensus
  subs_primer=c()
  for(j in 0:primer_len){
    temp_subs=sum( (consPrimer[j+1]!=as.character(DNA_mat[,pos+j]))==T)
    subs_primer=c(subs_primer,ifelse(temp_subs>9,"+",temp_subs))
  }
  
  
  degens=unlist(strsplit(subs,""))
  degPrimer=consPrimer
  if(any(degens!="0"))  {
    for( i in which(degens!="0")) {
      degPos=i+pos
      degPrimer[i]=bma(names(table(DNA_mat[,degPos-1])))
    }
  }
  
  cat("Consensus:\t",consPrimer,"\n#wobbles:\t",subs,"\n#substitutions:\t",subs_primer, "\nPrimer:\t\t",degPrimer,"\n",sep="")
  
}

#helper function for printing primers
printPair=function(pair){
  
  cat(paste("First:\n\n"))
  print(goodPairs[pair[1],])
  cat("\n")
  getAln(c(goodPairs[pair[1],1]))
  
  cat(paste("\nSecond:\n\n"))
  
  print(goodPairs[pair[2],])
  cat("\n")
  getAln(c(goodPairs[pair[2],1]))
  
}



#address of fasta file
fasta=read.fasta(ALNS)

#turning into a matrix
DNA_mat=do.call(rbind,fasta)

#finding dimensions of data
max_len=NCOL(DNA_mat)
n_seqs=NROW(DNA_mat)

#Finding consensus
CONSENSUS=apply(DNA_mat, MARGIN = 2,function(x) names(sort(table(x),decreasing = T))[1])

#making empty vector for diversity across string
divs  = rep(-1, max_len)
sumDF = data.frame("a"= rep(-1,max_len),"g"=rep(-1,max_len),"c"=rep(-1,max_len),"t"=rep(-1,max_len),"ins"=rep(-1,max_len))

i=10
#looping across all nucleotide positions
for( i in 1:max_len) {
  #turn each position into a table with of counts for each nucleotide, while including 0s.
  per_pos=factor(DNA_mat[,i], levels=c("a","g", "c", "t", "-"))
  probs=table(per_pos)/n_seqs
  divs[i]=-sum(probs*ifelse(log(probs)==-Inf, 0, log(probs))) #shannon index
  sumDF[i,]=probs
}

#calculation rolling mean of diversity
roll_means_30=rollmean(divs, k = 10)
print(sum(divs)/max_len)

#plotting diversity
jpeg(width=800*3,height=3*600,res=72*3,file = args[2])
plot(divs, pch=16, cex=0.1, ylim=c(-0.1,max(divs)),type="p", main="Diversity", xlab="Nucleotide position", ylab="Diversity")
lines(roll_means_30, col=3)
legend("topright", legend = "Rolling mean", col=3, lty=1)
dev.off()

######
#Making primer candidates
#make kmers of diversity vector, see which ones are conserved (low diversity)
#and see if any overlap correct amplicon length
######

#making data frame for for Kmers/primers and their diversity
kmers=data.frame(pos=rep(0,length(divs)-primer_len),kmer=rep("",length(divs)-primer_len), 
                 divs=rep(0,length(divs)-primer_len),GC=rep(0,length(divs)-primer_len))
j=1
for(j in seq(1, length(divs)-primer_len, by = 1)){
  kmers$pos [j] = j
  kmers$kmer[j] = paste(CONSENSUS[c(j:(j+primer_len))],collapse ="")
  kmers$divs[j] = sum(divs[c(j:(j+primer_len))])
  kmers$GC  [j] = length(which(CONSENSUS[c(j:(j+primer_len))]=="g" | CONSENSUS[c(j:(j+primer_len))] == "c"))/primer_len
  
}  

#try primer combinations whilst progressively increasing
#diversity cutoff untill a valid pair is found
no_primers=F
keep_going=T
while(keep_going) {
  
  #working out kmer pairs having low diversity
  goodPairs=kmers[which(kmers$divs<=div_cut),]
  
  #if there are no good pairs, increase cutoff and try again
  if(NROW(goodPairs)==0){
    keep_going=T 
    div_cut=div_cut+0.05
    next
  } else {
    keep_going=F
  }
  
  #if good pairs exist, create matrix of all pairs
  #sum of pair diversity
  divMat  = matrix(-.1,nrow = NROW(goodPairs),ncol = NROW(goodPairs))
  #difference in pair GC content
  dGC     = matrix(-.1,nrow = NROW(goodPairs),ncol = NROW(goodPairs))
  #lengt of pair amplicon
  lenMat  = matrix(-.1,nrow = NROW(goodPairs),ncol = NROW(goodPairs))
  #home-made score of primer pair
  #1/sqrt(diversity^2 + 10*deltaGC^2)
  scoreMat= matrix(-.1,nrow = NROW(goodPairs),ncol = NROW(goodPairs))
  
  i=1;j=6
  #run through all primer combinations
  for(i in 1:NROW(divMat)) {
    for(j in 1:NCOL(divMat)) {
      
      #first calculate amplicon length
      lenMat[i,j]=abs(goodPairs$pos[i] - goodPairs$pos[j])
      
      #if the length is good AND difference in GC i is less than GC cutoff
      #calculate div, dGC and score
      if(lenMat[i,j]>amplion_min_len & lenMat[i,j]<amplion_max_len & (abs(goodPairs$GC[i] - goodPairs$GC[j]) <GC_tol)) {
        divMat[i,j]   = goodPairs$divs[i] + goodPairs$divs[j]
        dGC[i,j]      = abs(goodPairs$GC[i] - goodPairs$GC[j])
        scoreMat[i,j] = 1/((sqrt(divMat[i,j]^2) + 10*dGC[i,j]^2)+.01)
      }
    }
  }
  
  #if none of the candidates still wont work, then increase div cutoff and try again
  if(all(scoreMat==-.1 ) ){
    keep_going=T 
    div_cut=div_cut+0.05
  } else {
    keep_going=F
  }
  
  if(div_cut>3){
    no_primers=T
    break
  }
}

if(no_primers) {
  cat(paste("No primers found. Increase GC cutoff.\n"))
} else {
  #pheatmap(divMat,cluster_rows = F,cluster_cols = F)
  #pheatmap(dGC,cluster_rows = F,cluster_cols = F)
  pheatmap(scoreMat,cluster_rows = F,cluster_cols = F,labels_col = goodPairs$pos, labels_row = goodPairs$pos)
  
  unqScores=sort(unique(as.numeric(scoreMat)))[-1]
  counter=1
  for(i in unqScores)  {
    cat(paste("\n***pair", counter,"\n"))
    Pair=which(scoreMat==i,arr.ind = T)[1,]
    printPair(Pair)
    cat(paste("\n--------------\n\n"))
    counter=counter+1
  }
}
