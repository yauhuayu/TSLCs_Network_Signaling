# ===========================================================================
# Extract related GeneIDs from Entrez Gene queries of a specific topic
# ===========================================================================

#!/usr/bin/perl

open(inFile, "gene_20090122_EMT.txt") || die "Cannot open input file: $!\n";; #open the data file
open(out, ">g_entrez090122_emt.txt"); #output file

#$line=<inFile>; #throw out the first line (column heading)
#print $line;


while ($line=<inFile>)   #read in data file one line at a time
{
	if (($line ne "") && ($line ne "\n"))  #not an empty line
	{ 
		chomp $line; #remove end-of-line character (\n)
		if($line =~ /^(\d+):/ ){
			$line =~ m/(\d+): (.*)/s;
			my $sym = $2;
			print out $sym."\t";
		}
		
		if($line =~ /Official Symbol/ ){
			$line =~ m/Official Symbol\s(.*)and Name:(.*)\[Homo sapiens\]/s;
			my $gene = $1."\t".$2."\t";
			print out $gene."\t";
		}
				
		if($line =~ /Chromosome:/ ){
			$line =~ s/;//;
			$line =~ m/Chromosome: (.*) Location: (.*)/s;
			my $ck = $1;
			my $loc = $2;
			my $chrom = $ck."\t".$loc;
			print out $chrom."\t";
		}
		if($line =~ /GeneID:/ ){
			$line =~ m/GeneID: (.*)/s;
			my $id = $1;
			print out $id."\n";
		}		
 	}
}



close out;

