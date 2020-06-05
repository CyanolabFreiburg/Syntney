
# GLASSgo postprocessing script

# dependencies: 
require(stringi)


#CALL:
#R --slave -f  GLASSgo_postprocessing_11_sqlite.r --args filename=candidates.fasta synteny_window=3000 script_path=~/media/cyano_share/data/TOOLS/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py db_path=~/media/jens@margarita/jensSicherung/GLASSgo2/mySQLiteDB_new.db"


filename<-"candidates.fasta" # result fasta file from GLASSgo
script_path<-"~/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py"
db_path<-"~/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py"


synteny_window<-3000 # number of bases upstream and downstream of the sRNA that were searched for protein coding genes for the synteny analysis
random_extension=T # make locus_tags unique by adding random extensions

args <- commandArgs(trailingOnly = TRUE) 

for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }
 
synteny_window<-as.numeric(synteny_window)



split_glassgo<-function(x){
	tmp<-strsplit(x[1], ">")[[1]]
	if(length(tmp)>2){
		tmp2<-strsplit(tmp[length(tmp)],"p.c.VAL:")[[1]]
		tmp2<-gsub(";.*","",tmp2[length(tmp2)])
		tmp<-paste(tmp[2], "p.c.VAL:", tmp2[length(tmp2)], sep="")
	}else {
			tmp<-tmp[2]
		}
	tmp<-strsplit(tmp, ":")[[1]]
	id<-tmp[1]
	taxid<-tmp[length(tmp)]
	identity<-strsplit(tmp[length(tmp)-1],"-taxID")[[1]][1]
	tmp2<-strsplit(tmp[2]," ")[[1]]
	tmp3<-paste(tmp2[2:(length(tmp2)-1)],collapse=" ")
	name<-strsplit(tmp3,",")[[1]][1]
	tmp2<-strsplit(tmp2,"-")[[1]]
	a<-tmp2[1]
	b<-tmp2[2]
	a<-strsplit(a,"c")[[1]]
	strand<-"+"
	if(length(a)==2){
		st<-b
		en<-a[2]
		strand<-"-"
	}else{
		st<-a
		en<-b
	}
	out<-c(id, strand, st,en,name)
	out
}

export_ncRNA_coordinates<-function(x){ 
	header_row <- grep(">", x)
	headers <- as.character(x[header_row])
	headers<-headers[2:length(headers)]
	seqs<-as.character(x[header_row+1])
	seqs<-seqs[2:length(seqs)]
	tmp<-do.call(rbind,lapply(headers,split_glassgo))
	tmp<-cbind(tmp,headers,seqs)
	colnames(tmp)<-c("Accesion_number", "Strand","start","end","name","Full_header","sequence")
	tmp
}


locus_tag2org<-function(out2){
	tag<-c()
	org<-c()
	for(i in 1:length(out2)){
		temp_tag<-out2[[i]][,5]
		temp_org<-rep(names(out2)[i], length(temp_tag))
		tag<-c(tag,temp_tag)
		org<-c(org,temp_org)
	}
	out<-cbind(tag,org)
	out
}


rand_extension<-function(x, Accession){
	ra<-stri_rand_strings(x,length=4,  pattern = "[A-Za-z0-9]")
	temp<-paste(Accession,ra,sep="_")
	#temp<-gsub("\"","",temp)
	
	temp
}

get_prot_fasta3<-function(out){
  fasta<-c()
  for(i in 1:length(out)){
    for(j in 1:nrow(out[[i]])){
        if(is.na(out[[i]][j,6])==F){
        temp<-as.character(out[[i]][j,6])
		na<-as.character(out[[i]][j,5])
		na<-gsub("\\\"","",na)
         na<-paste(">",na,sep="")
		 temp<-c(na,temp)
         fasta<-c(fasta,temp)
        }
      }
    }
  write.table(fasta, file="protein_fasta.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}


# identify homologous proteins using CDhit
cdhit_run<-function(fasta="protein_fasta.txt", outname="psi", thres=0.3, psi=T){
  wd<-getwd()
  di<-paste(wd,"/", "psi_out",sep="")
  dir.create(di)
  if(psi==T){
    inp<-paste("./psi-cd-hit.pl -i ", fasta,  " -d 50 -o ", di,"/",outname, " -c ", thres, sep="")
  }
  if(psi==F){
    inp<-paste("cd-hit -i ", fasta,  " -d 50 -o ", di,"/",outname, " -c ",  thres ," -n 2", " -aL 0.6", sep="")
  }
  print(inp)
  system(inp)
  cd<-paste(di, "/", outname, ".clstr", sep="")
  cd<-read.delim(cd, header=F, sep="?")
  cd<-as.character(cd[,1])
  cd<-gsub("\t"," ", cd)
  cd
}

proc_cdhit<-function(x){ 
  clustlist<-list()
  numb<-grep(">Cluster", x)
  for(i in 1:length(numb)){
    if(i<length(numb)){
      end<-numb[i+1]-1
    }
    if(i==length(numb)){
      end<-length(x)
    }
    temp<-x[(numb[i]+1):end]
    temp<-gsub(".*aa, >","",temp)
    temp<-gsub("\\.\\.\\..*","",temp)
    clustlist[[i]]<-temp
  }
  clustlist	
}


# Execute functions



fasta<-read.delim(filename, header=F, sep="\t")
fasta<-as.character(fasta[,1])
coor<-export_ncRNA_coordinates(fasta)


s<-as.numeric(coor[,3])
e<-as.numeric(coor[,4])
m<-round(s+(e-s)/2,digits=0)

wi<-rep(synteny_window,nrow(coor))

coor3<-cbind(paste(coor[,1],"_",coor[,3],sep=""),coor[,1],m,wi)

write.table(coor3,file="coordinates.txt", sep="\t", row.names=F, col.names=F, quote=F)

command<-paste("python3 ", script_path, " -s ", db_path, " -a coordinates.txt")

#print(command)
dat<-do.call(rbind,strsplit(system(command, intern=T),"\t")) 


#print(dat)
unlink("coordinates.txt")
na<-which(dat[,3]=="no annotation")
na2<-which(dat[,3]=="missing entry in LUT")
na<-c(na,na2)

if(length(na)>0){
	dat<-dat[-na,]
}
ids<-unique(dat[,1])

id2<-paste(coor[,1],coor[,3],sep="_")

out<-vector("list", length(ids))
names(out)<-ids

for(i in 1:length(ids)){
	tmp<-which(dat[,1]==ids[i])
	srna<-match(ids[i],id2)
	stra<--1
    if(coor[srna,2]=="+"){
      stra<-1
    }
	temp_out<-cbind(dat[tmp,7],dat[tmp,5],dat[tmp,6],dat[tmp,3],dat[tmp,4],dat[tmp,8])
	colnames(temp_out)<-c("strand","start","end","gene_name","locus_tag","AA_sequence")
	ma<-max(as.numeric(temp_out[,3]))
	mi<-min(as.numeric(temp_out[,2]))
	aa<-as.numeric(temp_out[,2])-mi
	bb<-as.numeric(temp_out[,3])-mi
	s_srna<-min(as.numeric(coor[srna,3:4]))-mi
	e_srna<-max(as.numeric(coor[srna,3:4]))-mi
	nan<-which(is.na(temp_out[,"locus_tag"]))
	na<-which(temp_out[,"locus_tag"]=="na")
	nan<-unique(c(na,nan))
	if(length(nan)>0){
		temp_out[nan,"locus_tag"]<-unlist(lapply(length(nan),rand_extension,Accession=ids[i]))
	}
	if(random_extension==T){
		if(length(rest)>0){
			temp_out[rest,"locus_tag"]<-unlist(lapply(length(rest),rand_extension,Accession=temp_out[rest,"locus_tag"]))
		}
	}
	temp_out<-data.frame(temp_out,aa,bb,rep(s_srna,nrow(temp_out)),rep(e_srna,nrow(temp_out)),rep(stra,nrow(temp_out)),rep(coor[srna,5],nrow(temp_out)))

	
	out[[i]]<-temp_out
	
}

tagtable<-locus_tag2org(out)
get_prot_fasta3(out)
cd<-cdhit_run(psi=F,thres=0.4)
cd<-proc_cdhit(cd)



unlink("protein_fasta.txt")
 cluster_table<-function(cd){
	out<-c()
	for( i in 1:length(cd)){
		temp<-paste(cd[[i]], collapse=",")
		out<-c(out,temp)
	}
	names(out)<-paste("cluster_",1:length(out),sep="")
	out
}

cluster<-cluster_table(cd)
write.table(cluster, file=paste(name,"cluster_table.txt",sep="_"), sep="\t", quote=F)
 
synt_table<-function(out3, coor){
	res<-matrix(,length(out3),6)
	colnames(res)<-c("ID","organism","accesion","refseqID","neighbourhood_genes","position2sRNA")
	for(i in 1:length(out3)){
		tmp<-out3[[i]][order(out3[[i]][,"aa"]),]
		pos<-which(tmp[,"aa"]<tmp[,"rep.s_srna..nrow.temp_out.."])
		pos2<-seq(1,nrow(tmp)-length(pos))
		pos<-c(rev(pos),pos2)
		res[i,1]<-names(out3)[i]
		res[i,2]<-coor[match(gsub("_.*","",names(out3)[i]), coor[,1]),"name"]
		res[i,3]<-coor[match(gsub("_.*","",names(out3)[i]), coor[,1]),1]
		res[i,4]<-coor[match(gsub("_.*","",names(out3)[i]), coor[,1]),1]
		res[i,5]<-paste(gsub("\"","",tmp[,"locus_tag"]),collapse=",")
		res[i,6]<-paste(pos,collapse=",")
	}
	res
}
 
synteny<-synt_table(out,coor)
write.table(synteny, file=paste(name,"synteny_table.txt",sep="_"), sep="\t", quote=F, row.names=F)






