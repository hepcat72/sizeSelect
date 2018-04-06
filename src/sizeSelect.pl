#!/usr/bin/perl

#USAGE: Run with no options to get usage or with --help for basic details

use CommandLineInterface;
use warnings;
use strict;

our $VERSION = '2.1';

setScriptInfo(CREATED => '4/5/2018',
              VERSION => $VERSION,
              AUTHOR  => 'Robert William Leach',
              CONTACT => 'rleach@princeton.edu',
              COMPANY => 'Princeton University',
              LICENSE => 'Copyright 2018',
              HELP    => << 'END_HELP'

Filters a sequence file (fastq or fasta) for sequence size.  Example:

    sizeSelect.pl -s 20 -i sequences.fa --verbose --noheader > sequences20.fa

END_HELP
	     );

setDefaults(HEADER        => 1,
	    ERRLIMIT      => 3,
	    COLLISIONMODE => 'error', #,merge,rename (when outfile conflict)
	    DEFRUNMODE    => 'usage');

my $seqfileid =
  addInfileOption(GETOPTKEY   => 'i|seq-file|infile=s',
		  PRIMARY     => 1,
		  FLAGLESS    => 1,
		  REQUIRED    => 1,
		  SMRY_DESC   => "Sequence file (fastq or fasta).",
		  DETAIL_DESC => << 'end_detail'

A fasta or fastq sequence file to be filtered for sequence size.

end_detail
		  ,
		  FORMAT_DESC => << 'end_format'

Standard fasta or fastq format sequence file.

end_format
);

my $filetype = 'auto';
addOption(GETOPTKEY   => 't|filetype=s',
	  GETOPTVAL   => \$filetype,
	  ACCEPTS     => ['fastq','fastq','auto'],
	  DEFAULT     => $filetype,
	  DETAIL_DESC => join('',('Input file (-i) type.  Using this instead ',
				  'of auto will make file reading faster.  ',
				  '"auto" cannot be used when redirecting a ',
				  'file in.')));

my $sample_size = 0;
addOption(GETOPTKEY   => 'n|sample-size=s',
	  GETOPTVAL   => \$sample_size,
	  DEFAULT     => $sample_size,
	  DETAIL_DESC => join('',('Check at least this many sequences.  This ',
				  'is only used when -s is not provided.  ',
				  'Ignored otherwise.  If 0 is supplied, all ',
				  'sequences will be checked (unless -m has ',
				  'a value).  If both -m and -n are set, ',
				  'both requirements must be satisfied in ',
				  'order to stop checking sequence ',
				  'lengths.')));

my $order_thresh = 3;
addOption(GETOPTKEY   => 'm|order-of-mag-thresh=s',
	  GETOPTVAL   => \$order_thresh,
	  DEFAULT     => $order_thresh,
	  DETAIL_DESC => join('',
			      ('If -s is not supplied, check as many ',
			       'sequences as it takes until the most ',
			       'abundant sequence length is larger than or ',
			       'equal to the number of the rest of the ',
			       'sequences to this order of magnitude.  E.g. ',
			       'If `-m 3` is supplied, the search for the ',
			       'most abundant sequence length will stop ',
			       'when: most_abundant >= (the_remainder * ',
			       '10^3).  If both -m and -n are set, both ',
			       'requirements must be satisfied in order to ',
			       'stop checking sequence lengths.')));

my $outfileid =
  addOutfileSuffixOption(GETOPTKEY   => 'o|outfile-siffix|extension=s',
			 FILETYPEID  => $seqfileid,
			 PRIMARY     => 1,
			 SMRY_DESC   => join('',('Outfile suffix appended to ',
						 'files supplied via -i to ',
						 'create outfile names.')),
			 DETAIL_DESC => join('',
					     ('Outfile suffix appended to ',
					      'files supplied via -i to ',
					      'create outfile names.  Note, ',
					      'the output format will match ',
					      'the input format.  This ',
					      'script is not intended to ',
					      'convert between versions.  ',
					      'Refer to convertSeq.pl for ',
					      'format conversion.')),
			 FORMAT_DESC => 'Same as input file format (-i).');

my $sizes = [];
addArrayOption(GETOPTKEY   => 's|size=s',
	       GETOPTVAL   => $sizes,
	       DEFAULT     => 'auto',
	       SMRY_DESC   => 'Sequence size to select.',
	       DETAIL_DESC => join('',
				   ('Only sequences of this size will be ',
				    'output.  All other sizes will be ',
				    'skipped/ignored.  May be supplied ',
				    'multiple times (or spece-delimited) ',
				    'with different values (e.g. `-s 20 -s ',
				    '30 -s 40` or `-s "20 30 40"`).  ',
				    'Ranges are also accepted in the form of ',
				    '"20-30", "20-", or "-30" in order to ',
				    'select sizes 20-30, 20-max, and 0-20 ',
				    'respectively.  An empty string value ',
				    'will select all sequence sizes (for ',
				    'purposes of pipelining).  The same ',
				    'size(s) will be selected in each input ',
				    'file.')));

processCommandLine();

my @errs   = ();
my $all    = 0;
my $over   = '';
my $ranges = {};
foreach my $size_str (map {s/^-/0-/;$_} @$sizes)
  {
    if($size_str =~ /[^\d\-]/ || $size_str =~ /-.*-/)
      {push(@errs,$size_str)}
    elsif($size_str eq '')
      {$all = 1}
    else
      {
	if($size_str =~/^(\d+)-(\d+)$/)
	  {
	    my($l,$g) = sort {$a <=> $b} ($1,$2);
	    if(!exists($ranges->{$l}) ||
	       (exists($ranges->{$l}) && $g > $ranges->{$l}
		&& $ranges->{$l} != -1))
	      {$ranges->{$l} = $g}
	  }
	elsif($size_str =~/^(\d+)-$/)
	  {$ranges->{$1} = -1}
	elsif($size_str =~/^(\d+)$/)
	  {$ranges->{$1} = $1 unless(exists($ranges->{$1}))}
	else
	  {push(@errs,$size_str)}
      }
  }

if(scalar(@errs))
  {
    error("Invalid sizes supplied to -s: [",join(' ',@errs),"].");
    usage(1);
    exit(1);
  }

my $auto = !scalar(keys(%$ranges));

if($filetype =~ /fasta/i)
  {$filetype = 'fasta'}
elsif($filetype =~ /fastq/i)
  {$filetype = 'fastq'}
elsif($filetype =~ /auto/)
  {$filetype = 'auto'}
else
  {
    error("Unrecognized file type: [$filetype].  Must be 'fasta', 'fastq', ",
	  "or 'auto'.");
    quit(1);
  }

while(nextFileCombo())
  {
    my $infile  = getInfile($seqfileid);
    my $outfile = getOutfile($outfileid);

    if($auto)
      {
	my $s = determineAdundantSize($infile,$sample_size,$order_thresh);
	$ranges->{$s} = $s;
      }

    openIn(*INF,$infile)    || next;
    openOut(*OUTF,$outfile) || next;

    if($all)
      {print(<INF>)}
    else
      {
	while(my $rec = getNextSeqRec(*INF,1,$infile))
	  {
	    my $s = length(formatSequence($rec->[1],
					  0,0,0,undef,undef,undef,0,0,0));

	    if(scalar(grep {my($l,$g)=($_,$ranges->{$_});
			    $s >= $l && ($s <= $g || $g == -1)}
		      keys(%$ranges)))
	      {printSeqRec($rec,$infile)}
	  }
      }

    closeOut(*OUTF);
    closeIn(*INF);
  }

















#Returns the most abundant sequence length.  Looks at all sequences by default.
#Looks at $sample_size is provided (and order_thresh not provided).  Looks at
#as many sequences as it takes for the most abundant sequence to be more
#abundant than the remaining number of sequences to the indicated order of
#magnitude (if sample_size is not provided).  If both sample_size and
#order_thresh are provided, looks until both requirements are satisfied.
sub determineAbundantSize
  {
    my $seqfile      = $_[0];
    my $sample_size  = defined($_[1]) ? $_[1] : 0;
    my $def_order    = 3;
    my $order_thresh = defined($_[2]) ? $_[2] : $def_order;

    if($order_thresh < 1)
      {
	warning("Invalid order of magnitude threshold: [$order_thresh].  ",
		"Must be greater than or equal to 1.  Setting default of ",
		"[$def_order].");
	$order_thresh = $def_order;
      }

    openIn(*SIZEF,$seqfile) || return({0 => -1});

    my $count     = 0;
    my $hash      = {};
    my $max_len   = 0;
    my $max_abund = 0;

    while(my $rec = getNextSeqRec(*SIZEF,0,$seqfile))
      {
	my $len = length($rec->[1]);
	$hash->{$len}++;
	if($hash->{$len} > $max_abund)
	  {
	    $max_abund = $hash->{$len};
	    $max_len = $len;
	  }
	$count++;
	my $remainder = $count - $max_abund;
	if($sample_size > 0 && $count >= $sample_size &&
	   $max_abund >= ($remainder * 10**$order_thresh))
	  {last}
      }

    closeIn(*SIZEF);

    return($max_len);
  }

#Globals used: $filetype, $main::lastfiletype
sub printSeqRec
  {
    my $rec = $_[0];

    #Both formats have a defline and a sequence line first
    print($rec->[0],"\n",$rec->[1],"\n");

    #We can assume it's fastq if the rec array is larger than 2
    if($#{$rec} > 1)
      {print($rec->[3],"\n",$rec->[2],"\n")}
  }

#Uses global variables: lastfiletype & filetype
sub getNextSeqRec
  {
    my $input_file = $_[2];

    debug("Determining previous type");

    if(!defined($main::lastfiletype) || $filetype ne 'auto')
      {
	if($filetype eq 'fasta')
	  {$main::getnextsub = \&getNextFastaRec}
	elsif($filetype eq 'fastq')
	  {$main::getnextsub = \&getNextFastqRec}
      }
    elsif(defined($main::lastfiletype) &&
	  exists($main::lastfiletype->{$input_file}))
      {
	if($main::lastfiletype->{$input_file} eq 'fasta')
	  {$main::getnextsub = \&getNextFastaRec}
	elsif($main::lastfiletype->{$input_file} eq 'fastq')
	  {$main::getnextsub = \&getNextFastqRec}
      }

    if($filetype eq 'auto' &&
       (!defined($main::lastfiletype) ||
	!exists($main::lastfiletype->{$input_file})))
      {
	debug("Determining type");

	if($input_file eq '-')
	  {
	    error("`-t auto` cannot be used when the input file is supplied ",
		  "on standard input.  Please supply the exact file type.");
	    quit(2);
	  }

	if(!-e $input_file)
	  {
	    error("`-t auto` cannot be used when the input file does not ",
		  "exist.  Please supply the exact file type.");
	    quit(8);
	  }

	my($num_fastq_defs);
	if(-e $input_file)
	  {
	    $num_fastq_defs =
	      `head -n 50 "$input_file" | grep -c -E '^[\@\+]'`;
	    debug("System output from: [",
		  qq(head -n 50 "$input_file" | grep -c -E '^[\@\+]'),
		  "]:\n$num_fastq_defs");
	    $num_fastq_defs =~ s/^\D+//;
	    $num_fastq_defs =~ s/\D.*//;
	  }
	else
	  {$num_fastq_defs = 0}

	if($num_fastq_defs > 0)
	  {
	    $main::getnextsub = \&getNextFastqRec;
	    $main::lastfiletype->{$input_file} = 'fastq';
	  }
	else
	  {
	    my($num_fasta_defs);
	    if(-e $input_file)
	      {
		$num_fasta_defs = `head -n 50 "$input_file" | grep -c -E '^>'`;

		debug("System output from: [",
		      qq(head -n 50 "$input_file" | grep -c -E '^>'),
		      "]:\n$num_fasta_defs");

		$num_fasta_defs =~ s/^\D+//;
		$num_fasta_defs =~ s/\D.*//;
	      }
	    else
	      {$num_fasta_defs = 0}

	    if($num_fasta_defs > 0)
	      {
		$main::getnextsub = \&getNextFastaRec;
		$main::lastfiletype->{$input_file} = 'fasta';
	      }
	    else
	      {
		if(!defined($main::lastfiletype) ||
		   !exists($main::lastfiletype->{$input_file}))
		  {
		    debug("Num fasta deflines: [$num_fasta_defs].");
		    error("Unable to determine file type.  Skipping file ",
			  "[$input_file].");
		    return(undef);
		  }
		warning("Unable to determine file type.  Defaulting to ",
			"[$main::lastfiletype->{$input_file}].");
		if($main::lastfiletype->{$input_file} eq 'fasta')
		  {$main::getnextsub = \&getNextFastaRec}
		else
		  {$main::getnextsub = \&getNextFastqRec}
	      }
	  }
      }

    debug("Returning record");

    return($main::getnextsub->(@_));
  }

#Copied from DNAstiffness.pl on 2/12/2014 so as to be independent -Rob
sub getNextFastaRec
  {
#    my $self       = shift(@_);
    my $handle    = $_[0];      #File handle or file name
    my $no_format = $_[1];

    if(exists($main::{FASTABUFFER}) && exists($main::{FASTABUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTABUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTABUFFER}->{$handle}});
		@{$main::{FASTABUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTABUFFER}->{$handle}}));
	  }
	elsif(eof($handle))
	  {return(undef)}
      }

    my $parent_id_check = {};
    my $first_loop      = 0;
    my $line_num        = 0;
    my $verbose_freq    = 1000;
    my $line            = '';
    my $defline         = '';
    my $seq_lines       = 0;
    my($seq);

    #For each line in the current input file
    while(getLine($handle))
      {
	$line_num++;

	verboseOverMe("Reading line [$line_num].")
	  unless($line_num % $verbose_freq);

	$line = $_;

	next if($line !~ /\S/ || $line =~ /^\s*#/);
	if($line =~ />/)
	  {
	    if($defline)
	      {
		my $solidseq = (defined($seq) ?
				($seq_lines == 1 || $no_format ?
				 $seq : formatSequence($seq)) : '');
		chomp($solidseq);
		chomp($defline);

		push(@{$main::{FASTABUFFER}->{$handle}},[$defline,$solidseq]);
	      }
	    $defline   = $line;
	    $seq_lines = 0;

	    my $tmp_id = $defline;
	    $tmp_id =~ s/^\s*>\s*//;
	    $tmp_id =~ s/\s.*//;
	    if($tmp_id eq '')
	      {warning("No Defline ID on line: [$line_num] of current file.  ",
		       " Universal coordinates will be used if some were ",
		       "supplied either via command line arguments of via ",
		       "coordinate file with no parent sequence ID.")}
	    elsif(exists($parent_id_check->{$tmp_id}))
	      {
		error("Two sequences found with the same ID on the ",
		      "defline: [$tmp_id] in current fasta file.  The same ",
		      "pairs of coordinates will be used for each sequence.");
	      }

	    undef($seq);
	  }
	elsif($line =~ /^([^\t]+?) *\t\s*(.*)/)
	  {
	    $defline = $1;
	    $seq     = $2;

	    my $solidseq =
	      ($seq_lines == 1 || $no_format ? $seq : formatSequence($seq));
	    chomp($solidseq);
	    chomp($defline);

	    push(@{$main::{FASTABUFFER}->{$handle}},[$defline,$solidseq]);

	    undef($seq);
	  }
	else
	  {
	    $seq .= $line;
	    $seq_lines++;
	  }
      }

    #Handle the last sequence (if there were any sequences)
    if(defined($seq))
      {
	my $solidseq = (defined($seq) ?
			($no_format ? $seq :
			 formatSequence($seq)) : '');
	chomp($solidseq);
	chomp($defline);

	push(@{$main::{FASTABUFFER}->{$handle}},[$defline,$solidseq]);
      }

    #Return the first sequence (if sequence was parsed)
    if(exists($main::{FASTABUFFER}) && exists($main::{FASTABUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTABUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTABUFFER}->{$handle}});
		@{$main::{FASTABUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTABUFFER}->{$handle}}));
	  }
	else
	  {return(undef)}
      }
    else
      {return(undef)}
  }

#Copied from DNAstiffness.pl on 2/12/2014 so as to be independent -Rob
sub formatSequence
  {
    #1. Read in the parameters.
    my $sequence          = $_[0];
    my $chars_per_line    = $_[1];
    my $coords_left_flag  = $_[2];
    my $coords_right_flag = $_[3];
    my $start_coord       = $_[4];
    my $coords_asc_flag   = $_[5];
    my $coord_upr_bound   = $_[6];
    my $uppercase_flag    = $_[7];
    my $print_flag        = $_[8];
    my $nucleotide_flag   = $_[9];

    my($formatted_sequence,
       $sub_string,
       $sub_sequence,
       $coord,
       $max_num_coord_digits,
       $line_size_left,
       $lead_spaces,
       $line);
    my $coord_separator = '  ';
    my $tmp_sequence = $sequence;
    $tmp_sequence =~ s/\s+//g;
    $tmp_sequence =~ s/<[^>]*>//g;
    my $seq_len = length($tmp_sequence);

    #2. Error check the parameters and set default values if unsupplied.
    my $default_chars_per_line    = ''; #Infinity
    my $default_coords_left_flag  = 0;
    my $default_coords_right_flag = 0;
    my $default_start_coord       = (!defined($coords_asc_flag) ||
				     $coords_asc_flag ? 1 : $seq_len);
    my $default_coords_asc_flag   = 1;
    my $default_coord_upr_bound   = undef();  #infinity (going past 1 produces
    my $default_uppercase_flag    = undef();  #          negative numbers)
    my $default_print_flag        = 0;

    if(!defined($chars_per_line) || $chars_per_line !~ /^\d+$/)
      {
        if(defined($chars_per_line) &&
	   $chars_per_line !~ /^\d+$/ && $chars_per_line =~ /./)
	  {print("WARNING:seq-lib.pl:formatSequence: Invalid ",
	         "chars_per_line: [$chars_per_line] - using default: ",
		 "[$default_chars_per_line]<BR>\n")}
        #end if(chars_per_line !~ /^\d+$/)
	$chars_per_line = $default_chars_per_line;
      }
    elsif(!$chars_per_line)
      {$chars_per_line = ''}
    #end if(!defined($chars_per_line) || $chars_per_line !~ /^\d+$/)
    if(!defined($coords_left_flag))
      {$coords_left_flag = $default_coords_left_flag}
    #end if(!defined(coords_left_flag))
    if(!defined($coords_right_flag))
      {$coords_right_flag = $default_coords_right_flag}
    #end if(!defined(coords_right_flag))
    if(!defined($start_coord) || $start_coord !~ /^\-?\d+$/)
      {
        if(defined($start_coord) &&
           ($coords_left_flag || $coords_right_flag))
          {print("WARNING:formatSequence.pl:formatSequence: Invalid ",
                 "start_coord: [$start_coord] - using default: ",
                 "[$default_start_coord]\n")}
        #end if($start_coord !~ /^\d+$/)
        $start_coord = $default_start_coord;
      }
    #end if(!defined($start_coord) || $start_coord !~ /^\d+$/)
    if(!defined($coords_asc_flag))
      {$coords_asc_flag = $default_coords_asc_flag}
    #end if(!defined(coords_right_flag))
    if(defined($coord_upr_bound) && $coord_upr_bound !~ /^\d+$/)
      {undef($coord_upr_bound)}
    if(!defined($print_flag))
      {$print_flag = $default_print_flag}
    #end if(!defined($print_flag))

    if(defined($coord_upr_bound) && $start_coord < 1)
      {$start_coord = $coord_upr_bound + $start_coord}
    elsif($start_coord < 1)
      {$start_coord--}
    elsif(defined($coord_upr_bound) && $start_coord > $coord_upr_bound)
      {$start_coord -= $coord_upr_bound}

    #3. Initialize the variables used for formatting.  (See the DATASTRUCTURES
    #   section.)
    if($coords_asc_flag)
      {
        if(defined($coord_upr_bound) &&
           ($seq_len + $start_coord) > $coord_upr_bound)
          {$max_num_coord_digits = length($coord_upr_bound)}
        else
          {$max_num_coord_digits = length($seq_len + $start_coord - 1)}

        $coord = $start_coord - 1;
      }
    else
      {
        if(defined($coord_upr_bound) && ($start_coord - $seq_len + 1) < 1)
          {$max_num_coord_digits = length($coord_upr_bound)}
        elsif(!defined($coord_upr_bound) &&
              length($start_coord - $seq_len - 1) > length($start_coord))
          {$max_num_coord_digits = length($start_coord - $seq_len - 1)}
        else
          {$max_num_coord_digits = length($start_coord)}

        $coord = $start_coord + 1;
      }
    $line_size_left = $chars_per_line;
    $lead_spaces    = $max_num_coord_digits - length($start_coord);

    #5. Add the first coordinate with spacing if coords_left_flag is true.
    $line = ' ' x $lead_spaces . $start_coord . $coord_separator
      if($coords_left_flag);

    #6. Foreach sub_string in the sequence where sub_string is either a
    #   sub_sequence or an HTML tag.
    foreach $sub_string (split(/(?=<)|(?<=>)/,$sequence))
      {
        #6.1 If the substring is an HTML tag
        if($sub_string =~ /^</)
          #6.1.1 Add it to the current line of the formatted_sequence
          {$line .= $sub_string}
        #end if(sub_string =~ /^</)
        #6.2 Else
        else
          {
            $sub_string =~ s/\s+//g;

	    if($nucleotide_flag)
	      {
		my(@errors);
		(@errors) = ($sub_string =~ /([^ATGCBDHVRYKMSWNX])/ig);
		$sub_string =~ s/([^ATGCBDHVRYKMSWNX])//ig;
		if(scalar(@errors))
		  {print STDERR ("WARNING:formatSequence.pl:formatSequence:",
				 scalar(@errors),
				 " bad nucleotide characters were ",
				 "filtered out of your sequence: [",
				 join('',@errors),
				 "].\n")}
	      }

            #6.2.1 If the sequence is to be uppercased
            if(defined($uppercase_flag) && $uppercase_flag)
              #6.2.1.1 Uppercase the sub-string
              {$sub_string = uc($sub_string)}
            #end if(defined($uppercase_flag) && $uppercase_flag)
            #6.2.2 Else if the sequence is to be lowercased
            elsif(defined($uppercase_flag) && !$uppercase_flag)
              #6.2.2.1 Lowercase the sub-string
              {$sub_string = lc($sub_string)}
            #end elsif(defined($uppercase_flag) && !$uppercase_flag)

            #6.2.3 While we can grab enough sequence to fill the rest of a line
            while($sub_string =~ /(.{1,$line_size_left})/g)
              {
                $sub_sequence = $1;
                #6.2.3.1 Add the grabbed sequence to the current line of the
                #        formatted sequence
                $line .= $sub_sequence;
                #6.2.3.2 Increment the current coord by the amount of sequence
                #        grabbed
                my $prev_coord = $coord;
                if($coords_asc_flag)
                  {
                    $coord += length($sub_sequence);
                    if(defined($coord_upr_bound)      &&
                       $prev_coord <= $coord_upr_bound &&
                       $coord > $coord_upr_bound)
                      {$coord -= $coord_upr_bound}
                  }
                else
                  {
                    $coord -= length($sub_sequence);
                    if(defined($coord_upr_bound) &&
                       $prev_coord >= 1 && $coord < 1)
                      {$coord = $coord_upr_bound + $coord - 1}
                    elsif($prev_coord >= 1 && $coord < 1)
                      {$coord--}
                  }
                #6.2.3.3 If the length of the current sequence grabbed
                #        completes a line
                if($line_size_left eq '' ||
		   length($sub_sequence) == $line_size_left)
                  {
                    $lead_spaces = $max_num_coord_digits - length($coord);
                    #6.2.3.3.1 Conditionally add coordinates based on the
                    #          coords flags
                    $line .= $coord_separator . ' ' x $lead_spaces . $coord
                      if($coords_right_flag);

                    #6.2.3.3.2 Add a hard return to the current line of the
                    #          formatted sequence
                    $line .= "\n";

                    #6.2.3.3.3 Add the current line to the formatted_sequence
                    $formatted_sequence .= $line;
                    #6.2.3.3.4 Print the current line if the print_flag is true
                    print $line if($print_flag);

                    #6.2.3.3.5 Start the next line
                    $lead_spaces = $max_num_coord_digits - length($coord+1);
                    $line = '';
                    $line = ' ' x $lead_spaces
                          . ($coords_asc_flag ? ($coord+1) : ($coord-1))
                          . $coord_separator
                      if($coords_left_flag);

                    #6.2.3.3.6 Reset the line_size_left (length of remaining
                    #          sequence per line) to chars_per_line
                    $line_size_left = $chars_per_line;
                  }
                #end if(length($sub_sequence) == $line_size_left)
                #6.2.3.4 Else
                else
                  #6.2.3.4.1 Decrement line_size_left (length of remaining
                  #          sequence per line) by the amount of sequence
                  #          grabbed
                  {$line_size_left -= length($sub_sequence)}
                #end 6.2.3.4 Else
              }
            #end while($sub_string =~ /(.{1,$line_size_left})/g)
          }
        #end 6.2 Else
      }
    #end foreach $sub_string (split(/(?=<)|(?<=>)/,$sequence))
    #7. Add the last coodinate with enough leadin white-space to be lined up
    #   with the rest coordinates if the coords_right_flag is true
    $lead_spaces = $max_num_coord_digits - length($coord);
    $line .= ' ' x $line_size_left . $coord_separator . ' ' x $lead_spaces
          . $coord
      if($coords_right_flag && $line_size_left != $chars_per_line);
    $line =~ s/^\s*\d+$coord_separator\s*$// if($coords_left_flag);

    #8. Add the ending PRE tag to the last line of the formatted sequence
    $line =~ s/\n+$/\n/s;

    #9. Add the last line to the formatted_sequence
    $formatted_sequence .= $line;
    #10. Print the last line if the print_flag is true
    print "$line\n" if($print_flag);

    if($coord < 1 && ($coords_left_flag || $coords_right_flag))
      {print("WARNING: The sequence straddles the origin.  Coordinates are ",
             "inaccurate.")}

    #11. Return the formatted_sequence
    return $formatted_sequence;
  }

#Merged the above getNextFastaRec subroutine with the code from convertSeq.pl
#on 2/12/2014 - have not yet tested
sub getNextFastqRec
  {
    my $handle    = $_[0];      #File handle or file name
    my $no_format = $_[1];

    if(exists($main::{FASTQBUFFER}) && exists($main::{FASTQBUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTQBUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTQBUFFER}->{$handle}});
		@{$main::{FASTQBUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTQBUFFER}->{$handle}}));
	  }
	elsif(eof($handle))
	  {return(undef)}
      }

    my $parent_id_check  = {};
    my $first_loop       = 1;
    my $line_num         = 0;
    my $line             = '';
    my $defline          = '';
    my $qual_defline     = '';
    my $seq              = '';
    my $qual             = '';
    my $getting_sequence = 0;
    my $comment_buffer   = '';
    my $seq_lines        = 0;
    my $qual_lines       = 0;
    my $verbose_freq     = 1000;

    #For each line in the current input file
    while(getLine($handle))
      {
	$line_num++;

	verboseOverMe("Reading line [$line_num].")
	  unless($line_num % $verbose_freq);

	$line = $_;

	next if($line !~ /\S/ || ($first_loop && $line =~ /^\s*#/));

	$first_loop = 0;

	#If this is the defline, or the quality length is the same as the seq
	if(length($qual) >= length($seq) && /^\s*\@[^\n\r]*/)
	  {
	    if($defline ne '' || $seq ne '' || $qual ne '')
	      {
		my $solidseq = ($seq ne '' ?
				($seq_lines == 1 || $no_format ? $seq :
				 formatSequence($seq)) : '');
		$qual =~ s/[\s\r\n\t]+//g if(!$no_format && $qual_lines > 1);
		chomp($solidseq);
		chomp($qual);
		chomp($defline);
		chomp($qual_defline);

		push(@{$main::{FASTQBUFFER}->{$handle}},
		     [$defline,$solidseq,$qual,$qual_defline]);
	      }
	    $defline    = $line;
	    $seq_lines  = 0;
	    $qual_lines = 0;

	    my $tmp_id = $defline;
	    $tmp_id =~ s/^\s*\@\s*//;
	    $tmp_id =~ s/\s.*//;

	    if($tmp_id eq '')
	      {warning("No Defline ID on line: [$line_num] of current file.")}
	    elsif(exists($parent_id_check->{$tmp_id}))
	      {error("Two sequences found with the same ID on the ",
		     "defline: [$tmp_id] in current fastq file.")}

	    $seq              = '';
	    $qual             = '';
	    $qual_defline     = '';
	    $getting_sequence = 1;
	  }
	elsif($getting_sequence && /^\s*\+[^\n\r]*/)
	  {
	    $getting_sequence = 0;
	    $qual_defline     = $line;
	  }
	elsif($getting_sequence)
	  {
	    s/\s+//g;
	    $seq_lines++;
	    if(/^[A-Za-z\n\.~]*$/)
	      {$seq .= $_}
	    else
	      {
		error("Expected a sequence character string, but ",
		      "got: [$_].  Appending anyway.");
		$seq .= $_;
	      }
	  }
	elsif($seq =~ /./)
	  {
	    s/\s+//g;
	    $qual_lines++;
	    if(/^[\!-\~]*$/)
	      {$qual .= $_}
	    else
	      {
		error("Expected a quality character string, but ",
		      "got: [$_].  Appending anyway.");
		$qual .= $_;
	      }

	    my $solidseq = ($seq ne '' ?
			    ($seq_lines == 1 || $no_format ? $seq :
			     formatSequence($seq)) : '');
	    my $solidqual = $qual;
	    $solidqual =~ s/[\s\r\n\t]+//g if(!$no_format && $qual_lines > 1);

	    if(length($solidseq) == length($solidqual))
	      {
		chomp($defline);
		chomp($qual_defline);
		return([$defline,$solidseq,$solidqual,$qual_defline]);
	      }
	  }
	#else must be a comment, ignore it
      }

    #Handle the last sequence (if there were any sequences)
    if($seq ne '')
      {
	my $solidseq = ($seq ne '' ?
			($no_format ? $seq :
			 formatSequence($seq)) : '');
	$qual =~ s/[\s\r\n\t]+//g unless($no_format);
	chomp($solidseq);
	chomp($defline);
	chomp($qual);
	chomp($qual_defline);

	push(@{$main::{FASTQBUFFER}->{$handle}},
	     [$defline,$solidseq,$qual,$qual_defline]);
      }

    #Return the first sequence (if sequence was parsed)
    if(exists($main::{FASTQBUFFER}) && exists($main::{FASTQBUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTQBUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTQBUFFER}->{$handle}});
		@{$main::{FASTQBUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTQBUFFER}->{$handle}}));
	  }
	else
	  {return(undef)}
      }
    else
      {return(undef)}
  }
