#!/usr/bin/perl1
use strict;
system("cls");
use warnings;

my @codons=();
my $min_ORF_length;
my $file;
my $name="";
my $sequence="";

print "Enter a minimum open reading length or type in default ";
chomp($min_ORF_length=<>);
if ($min_ORF_length eq "default") {$min_ORF_length=50;}

print "Enter name of FASTA file: ";
chomp($file=<>);
my $fasta_file=shift;

my $FASTA;
open($FASTA, $file) or die "can't open $fasta_file: $!\n";

my %sequences;
while (Get_Sequence_Info($FASTA, \%sequences)) {
   my $id="$sequences{header}";
   $name=$id;
   find_ORFS($sequences{sequence});
}
################ READS FASTA####################################

sub Get_Sequence_Info {
   my $fileInfo=$_[0];
   my $sequences = $_[1]; #Sequence hash
   
   # clears out previous sequence
   
   delete $sequences->{sequence};
   
   #Checks to see if it is the first header, or in other words if next_header has a value, if so it reasigns
   #header to be the value of the stored header that was found from the previous run at the end of the sequence
  
  if (exists($sequences{"next"}))
  {
	$sequences{"header"}=delete $sequences{"next"}; 
  }

   my $NotEmpty = ""; 
   while (<$fileInfo>) {
        $NotEmpty = "false";
	chomp;
       next if ($_=~/^\s*$/);  # <---skip blank lines

      # Finds the header line via regex and stores it
      #If no header line was found concatinate the sequence
      
      if ($_!~/^>/){    
         $_=&remove_white_space($_);  # remove any white space
         $sequences{"sequence"} .= $_;
      }
      elsif ($_=~/^>/) {
         
	 #When/If another header was found store it as next header for the next run so that it may be used as the primary header and return the
	 #hash with the header key and sequence value
	 #We have reached the end of the first sequence
	 
         if (exists($sequences{"header"})) {
             my $new_header = $_;    
         $new_header =$_;
	     $new_header=~ s/^>//;  #Stores the whole line minus 
	     $sequences{"next"} = $new_header;
            return $sequences;   
         }              
         
         # This is only for the First header of the file
	 
         else {
		my $header = $_;    
		$header =$_;
		$header=~ s/^>//;  #Stores the whole line minus >
            $sequences{"header"} = $header;
         }              
      }         
               
   }    
	#If the file is Empty return the hash with the last header and sequence
	#No header will be found so it must return one last time
	#Only for last sequence
   
   if ($NotEmpty ne "") {
      return $sequences;
   }    
    
}

################ORF FINDER####################################

#Sends to the ORF subroutine the sequence if frame is between 1 and 3 and the reverse complement if over 4: Scott

sub find_ORFS{
my $seq=$_[0];
my $frame=1;
while ($frame<=6){

#print $name."\n";
	
        if($frame<=3){&ORF($seq,$frame);}
	else {&ORF(&revcom($seq),$frame);}
	$frame++;
}
}

#Finds the ORFS: Scott

sub ORF(){
	my $start_position=-1;
	my $end_position=0;
	my $found=0;
	my $len;
	my $DNA=$_[0];
	my $fragment="";
	my $frame=$_[1];
        
	#Reasigns the offset for frames 4 and greater back to 1-3
	
        my $temp_frame=$frame;
	if ($frame ==4){$temp_frame=1;}
	elsif ($frame ==5){$temp_frame=2;}
	elsif ($frame ==6){$temp_frame=3;}
	for (my $i=$temp_frame-1;$i<length($DNA);$i+=3)
	{
		#makes a codon
		
                $fragment= substr $DNA,$i,3;
		
                #checks to see if it a start codon if so it records the start position
		
                if ( $fragment eq "ATG"){
			$start_position=$i+1;
		}
	
		if ($start_position!=-1){
			push( @codons, $fragment);
		}
		
                #Check to see if end codon if so sends the length and start postion along with frame to output
		
                if ($start_position!=-1 && ($fragment eq 'TAG' || $fragment eq 'TGA' || $fragment eq 'TAA')){
			$end_position=$i+3+1;
			my $length=$end_position-$start_position;
			&output($start_position,$length,$frame);
			$start_position=-1;
			undef @codons;
			
		}

	}
	
################OUTPUT####################################

#Teng

sub output(){
my ($pos, $length, $frame) = @_;

#assigns global array @codon to scalar variable $size
#The print statement outputs a headerline containing the "> name" followed by a bar "|" then Frame, Position and Length of the ORF.
#up to 15 codons per line will be output per ORF sequence

if ($length>=$min_ORF_length,){
	if($frame >=4) {$pos = $pos*-1;}
	my $format_seq = "(> $name| FRAME = $frame| POS = $pos|LEN= $length)\n";
	print "\n".$format_seq;
	my $size=@codons;
        
	#Prints 15 per line
	
        for (my $i=0; $i<$size; $i++) {
		if($i %15==0){print "\n";}
		print $codons[$i]." ";
        }
	print "\n";
}
}
}

#Tom

sub revcom() {
$_ = shift; 
$_ =~ tr/ACGTacgt/TGCAtgca/;    #complement
return scalar reverse("$_");	#reverse             
}

#Tom

sub remove_white_space() {
my $sequence=$_;

    #remove whitespace
    
    $sequence =~ s/\s+//;

    return $sequence;
}
