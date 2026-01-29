# Dated back to 2007~2009
# might not working with current version of IPA Svg
#
# Extracts all edge's interactions from html


#!/usr/local/bin/perl

my $type = "atrt";	# also the name of the directory

for ($i = 1; $i <= 4; $i++) {
my $n = $i;

my $dirOut = "out_".$type;
my $inF = "net".$n.".svg";
my $outF = ">".$dirOut."\\".$type."_svg".$n."_out.txt";
print $inF."\n";
open(net, $inF) || die "Cannot open input file: $!\n";; #open the data file

mkdir($dirOut);	# create directory
open(out, $outF); #output file

while(<net>){

	
   my @lines = split(/\n/, $_);
   my $part = @lines[0];
   #print out "check:\t$part"."\n\n";

   foreach my $line (@lines) {
        
        # --- proc genes --- example:
        #  <text x="459" xml:space="preserve"  y="-93.6914" clip-path="url(#clipPath2)" stroke="none">RARRES2</text>
        if($line =~ m/<text/ ){
        	$line =~ s/\"//g;  
        	$line =~ s/\s//g;   	 	
        	$line =~ m/<text(.*)xml:space=preserve(.*)clip-path(.*)>(.*)<\/text>/s;
            my $ncordx = $1;
			my $ncordy = $2;
        	my $garb = $3;
        	my $node = $4;
        	
        	if ($ncordx) {
        	print out $node."\t".$ncordx."\t".$ncordy."\n";
        	#print $node."\t".$ncordx."\t".$ncordy."\n";
        	}
        }
        
        # --- proc genes --- example:
        # <path fill="none" stroke-dasharray="12,4" d="M176.0736 226.2548 L219.9286 131.5308" clip-path="url(#clipPath2)"/>
        if($line =~ m/<path/ ){
        	$line =~ m/d=\"(.*)\s(.*)\s(.*)\s(.*)\" clip-path(.*)\/>/s;
            	my $path = $1."\t".$2."\t".$3."\t".$4;
            	my $ck = $4;
        	      	
        	if ($ck!="Z" ) {
        		if( !($1=~ m/L/g)) {
        			print out "path:\t".$path."\n";
        			#print  out $path."\n";
        		}
        	}
        }
        
   }  # foreach

}  # while 


close out;

}  # end of for loop

