
# GLASSgo postprocessing script

# dependencies: 
require(stringi)


#CALL:
#R --slave -f  ~/Syntney/packages/Rscript/Synteny_Cluster_Script_sqlite.r --args write_files=FALSE threads=10 filename=candidates.fasta synteny_window=3000 script_path=~/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py db_path=/media/cyano_share/exchange/Jens/Syntney/mySQLiteDB_new.db
#	python3 ~/Syntney/Syntney.py -i ~/Syntney/testfiles/candidates.fasta  -o ~/Syntney/testfiles/ -n cys -r off -d /media/cyano_share/exchange/Jens/Syntney/mySQLiteDB_new.db -c ~/Syntney/packages/Rscript/Synteny_Cluster_Script_sqlite.r  -s ~/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py

filename<-"candidates.fasta" # result fasta file from GLASSgo
script_path<-"~/Syntney/packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py"
db_path<-"/media/cyano_share/exchange/Jens/Syntney/mySQLiteDB_new.db"
threads<-30
name<-"sRNA"
write_files<-F

synteny_window<-3000 # number of bases upstream and downstream of the sRNA that were searched for protein coding genes for the synteny analysis

args <- commandArgs(trailingOnly = TRUE) 

for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }
 
synteny_window<-as.numeric(synteny_window)
threads<-as.numeric(threads)
write_files<-as.logical(write_files)

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
	first_line<-grep(":",headers[1])
	seqs<-as.character(x[header_row+1])
	if(length(first_line)==0){
		headers<-headers[2:length(headers)]
		seqs<-seqs[2:length(seqs)]
	}
	
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

get_prot_fasta3<-function(out, filen){
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
  write.table(fasta, file=filen, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}


# identify homologous proteins using CDhit
cdhit_run<-function(fasta="protein_fasta.txt", outname="psi", thres=0.3, psi=T, threads=2){
	tempf<-tempfile()
	wd<-getwd()
	di<-paste(wd,"/", "psi_out",sep="")
	dir.create(di)
	if(psi==T){
		inp<-paste("./psi-cd-hit.pl -i ", fasta,  " -d 50 -o ", tempf, " -c ", thres, sep="")
	 }
	if(psi==F){
		inp<-paste("cd-hit -i ", fasta,  " -d 50 -o ",tempf, " -c ",  thres ," -n 2", " -aL 0.6", " -T ", threads, sep="")
	}
	 # print(inp)
	 system(inp)
	 cd<-paste(tempf, ".clstr", sep="")
	 cd<-readLines(cd)
	 cd<-as.character(cd)
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



unlink("coordinates.txt")
na<-which(dat[,3]=="no annotation")
na2<-which(dat[,3]=="missing entry in LUT")
na<-c(na,na2)

no_anno<-dat[na,1:2]


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
	temp_out<-data.frame(temp_out,aa,bb,rep(s_srna,nrow(temp_out)),rep(e_srna,nrow(temp_out)),rep(stra,nrow(temp_out)),rep(coor[srna,5],nrow(temp_out)))
	out[[i]]<-temp_out
}
filen<-tempfile()
tagtable<-locus_tag2org(out)
get_prot_fasta3(out, filen)
cd<-cdhit_run(fasta=filen, psi=F,thres=0.4, threads=threads)
cd<-proc_cdhit(cd)

synt_table<-function(out3, coor){
	position<-c()
	res<-matrix(,length(out3),6)
	colnames(res)<-c("ID","organism","accesion","refseqID","neighbourhood_genes","position2sRNA")
	for(i in 1:length(out3)){
		tmp<-out3[[i]][order(out3[[i]][,"aa"]),]
		tmp_pos<-matrix(,nrow(tmp),2)
		tmp_pos[,2]<-"right"
		if(tmp[1,"rep.stra..nrow.temp_out.."]=="-1"){
			tmp_pos[,2]<-"left"
		}
		tmp_pos[,1]<-tmp[,"locus_tag"]
		pos<-which(tmp[,"aa"]<tmp[,"rep.s_srna..nrow.temp_out.."])
		if(tmp[1,"rep.stra..nrow.temp_out.."]=="-1"){
			tmp_pos[pos,2]<-"right"
		} else{
			tmp_pos[pos,2]<-"left"
		}
		
		genes<-tmp[,"locus_tag"]
		
		pos2<-seq(1,nrow(tmp)-length(pos))
		pos<-c(rev(pos),pos2)
		if(tmp[1,"rep.stra..nrow.temp_out.."]=="-1"){
			pos<-rev(pos)
			genes<-rev(genes)
		}		
		genes<-paste(gsub("\"","",genes),collapse=",")
		res[i,1]<-names(out3)[i]
		res[i,2]<-coor[match(gsub("_.*","",names(out3)[i]), coor[,1]),"name"]
		res[i,3]<-coor[match(gsub("_.*","",names(out3)[i]), coor[,1]),1]
		res[i,4]<-coor[match(gsub("_.*","",names(out3)[i]), coor[,1]),1]
		res[i,5]<-genes
		res[i,6]<-paste(pos,collapse=",")
		position<-rbind(position,tmp_pos)
	}
	return(list(res,position))
}
 
synteny<-synt_table(out,coor)

 cluster_table<-function(cd,dat, synteny){
	out<-c()
	anno<-matrix(,length(cd),3)
	colnames(anno)<-c("node","cluster_name","position_relative_to_sRNA")
	for( i in 1:length(cd)){
		temp<-paste(cd[[i]], collapse=",")
		out<-c(out,temp)
		tmp_pos<-na.omit(match(cd[[i]],synteny[[2]][,1]))
		if(length(tmp_pos)>0){
			tmp_pos<-synteny[[2]][tmp_pos,2]
			tmp_pos<-table(tmp_pos)
			tmp_pos<-names(tmp_pos)[which(tmp_pos==max(tmp_pos))[1]]
		} else {
			tmp_pos<-""
		}
		tmp<-na.omit(match(cd[[i]], dat[,4]))
		tmp<-dat[tmp,3]
		na<-which(tmp=="na")
		if(length(na)>0){
			tmp<-tmp[-na]			
		}
		if(length(tmp)>0){
			tmp<-table(tmp)
			tmp<-names(tmp)[which(tmp==max(tmp))[1]]
		} else {
			tmp<-cd[[i]][1]
		}
		anno[i,1]<-paste("cluster_",i,sep="")
		anno[i,2]<-tmp
		anno[i,3]<-tmp_pos
	}
	names(out)<-paste("cluster_",1:length(out),sep="")
	out<-list(out,anno)
}

cluster<-cluster_table(cd,dat, synteny)

if(write_files==TRUE){
	write.table(cluster[[1]], file=paste(name,"cluster_table.txt",sep="_"), sep="\t", quote=F)
	
	write.table(cluster[[2]], file=paste(name,"network_annotation.txt",sep="_"), sep="\t", quote=F,row.names = FALSE )
	write.table(synteny[[1]], file=paste(name,"synteny_table.txt",sep="_"), sep="\t", quote=F, row.names=F)	
} else {
	cat("#cluster_table\n")
	write.table(cluster[[1]], file=stdout(), sep="\t", quote=F,row.names=F, col.names=F)
	cat("#synteny_table\n")
	write.table(synteny[[1]], file=stdout(), sep="\t", quote=F, row.names=F)	
	cat("#network_annotation\n")
	write.table(cluster[[2]], file=stdout(), sep="\t", quote=F,row.names = FALSE )
	cat("#missing_data\n")
	x <- capture.output(write.table(no_anno, file=stdout(), sep="\t", quote=F,row.names = FALSE , col.names=F))
	cat(paste(x, collapse = "\n"))
}


cat(paste(x, collapse = "\n"), file = "df.csv") 