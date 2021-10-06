#!usr/bin/R
suppressMessages(library("optparse"))
suppressMessages(library("metagenomeSeq"))
 
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="OTUs counts.tsv", metavar="character"),
	make_option(c("-o", "--output"), type="character", default="normalized.tsv", help="output file name [default= %default]", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("I not found the input file .tsv", call.=FALSE)
}


OTUs=read.csv(opt$input, check.names=FALSE, header=TRUE, sep="\t")
OTUs=as.data.frame(OTUs)
row.names(OTUs)=OTUs$Genus
OTUs=OTUs[,-1]
Object= newMRexperiment(OTUs) 

metaSeqObject_CSS  = cumNorm(Object , p=cumNormStatFast(Object))

normalized= data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))
write.table(normalized,file=opt$output, sep="\t")
