# This script produces a matrix of VIPER regulator activities for an input matrix
#
# example files in /home/exacloud/lustre1/share_your_data_here/precepts/aracne_viper_analyses/kallisto/
# example expr = TCGA_BRCA_bygene_prepped_ensembl.txt.gz
# example network = network.txt
# tfs used to generate example network= pancancer-tfgenes_2014_04_17_filtered_brca_ensembl.txt
#
# Author: Hannah Manning & Joey Estabrook
# Date: April 9, 2018

#example of optparse usage https://gist.github.com/ericminikel/8428297
suppressPackageStartupMessages(require(optparse))
suppressWarnings(library(Biobase))
suppressWarnings(library(viper))

# function that loads expression matrix of unzipped or gzipped file
load_expr <- function(expr_file) {
	if (endsWith(expr_file, '.gz')) {
		expr <- as.matrix(read.table(gzfile(expr_file, 'rt'), row.names=1, sep='\t',
							check.names=FALSE, header=TRUE))
	} else {
		expr <- as.matrix(read.table(expr_file, row.names=1, sep='\t',
							check.names=FALSE, header=TRUE))
	}
	
	return(expr)
}

main <- function() {
	
	print("Running run_viper.R main()")
	
	# allow for commandline args
	option_list = list(
	make_option(c("-d", "--dir"), action="store", default=NA, type='character',
				help='Full path to outdir.'),
	make_option(c("-e", "--expr"), action="store", default=NA, type='character',
				help='Full path to input expression matrix used to build ARACNE network.'),
				
	make_option(c("-n", "--network"), action="store", default=NA, type='character',
				help='Full path to ARACNE network.txt.'),
	make_option(c("-s", "--smmartexpr"), action="store", default=NA, type='character',
				help='Full path to file containing expression of samples for which viper\n
				activity scores are desired.')		
	)
	
	opt = parse_args(OptionParser(option_list=option_list))
	
	# todo: make sure this works
	# todo: report which isn't working and exit
	if(any(is.na(c(opt$d, opt$e, opt$n)))) {
		cat("One or more variables have not been set!")
	}
	
	outdir=opt$d
	expr_aracne_fl=opt$e
	network_fl=opt$n
	expr_smmart_fl = opt$s
	
	# print session info
	writeLines(capture.output(sessionInfo()), paste(outdir, "RsessionInfo.txt", sep=""))
	
	print("Loading expression data from which ARACNE network was derived in R.")
	expr_a <- load_expr(expr_aracne_fl)
	
	print("Loading expression data of samples for which VIPER scores are desired.")
	expr_s <- load_expr(expr_smmart_fl)
								 
	print("Loading ARACNE network in R and converting to regulon obj")
	regul <- aracne2regulon(network_fl, expr_a)
	
	# TODO: Although regul will be generated with the tatlow data alone, vpres will be calculated using a different expr matrix
	# that expr matrix will be one that includes a UHR and smmart samples (a single smmart sample?)
	print("Generating activity scores")
	vpres <- viper(expr_s, regul, verbose=FALSE)
	
	print("Writing out activity scores")
	write.table(vpres, file=paste(outdir, 'viper_activities_ensembl.tsv', sep=''), sep='\t', quote=FALSE)
	
}


main()


