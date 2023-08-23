#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
my $first_position = 0;

my $metaparfile = "./cnf/meta.cnf";
# make it unbuffered
select STDIN; $| = 1;
select STDOUT; $| = 1;

my %metapar;
# read metaparameters file
open (META, "<$metaparfile") or die "Cant open $metaparfile: $!\n";
foreach my $line (<META>)
  {
    chomp($line);
    my @fields = split(/\s*=\s*/, $line);
    # remove spaces;
    $fields[0] = trim_spaces($fields[0]);
    $fields[1] = trim_spaces($fields[1]);
    $metapar{$fields[0]} = $fields[1];
  }
close(META);

my $start_codon_length=$metapar{start_length};
my $start_codon_offset=0;
my $intron_short_length = $metapar{intron_short_length};
my $intergenic_length = $metapar{intergenic_length};
my $stop_codon_length= $metapar{stop_codon_length};
my $stop_codon_offset=$metapar{stop_codon_offset};
my $acceptor_length= $metapar{acceptor_length};
my $initial_pattern_length = $metapar{acceptor_offset};
my $acceptor_initial_pattern_length = $metapar{acceptor_initial_pattern_length};
my $donor_initial_pattern_length = $metapar{donor_initial_pattern_length};
my $acceptor_offset=$metapar{acceptor_offset};
my $donor_length=$metapar{donor_length};
my $donor_offset=$metapar{donor_offset};
my $start_initial_pattern_length=$metapar{start_initial_pattern_length};
my $branch_length = $metapar{branch_length};

# Fixing offsets
my $fixed_donor_offset = $donor_offset + $donor_initial_pattern_length;
my $fixed_stop_offset = $stop_codon_offset;
my $donor_signal_length = $donor_length + $donor_initial_pattern_length;
my $acceptor_signal_length = $branch_length + $acceptor_length + $acceptor_initial_pattern_length;
my $start_signal_length = $start_codon_length + 3 + $start_initial_pattern_length;
my $stop_signal_length = $stop_codon_length;
my $exon_length_start = $start_codon_length - $start_codon_offset + $start_initial_pattern_length;
my $exon_length_stop = $stop_codon_offset;
my $exon_length_acceptor =  $acceptor_length - $acceptor_offset - 2 + $acceptor_initial_pattern_length;
my $exon_length_donor = $donor_offset + $donor_initial_pattern_length;
my $exon_delta_initial = $exon_length_start + $exon_length_donor;
my $exon_delta_internal = $exon_length_acceptor + $exon_length_donor;
my $exon_delta_final = $exon_length_acceptor + $exon_length_stop;
my $exon_delta_single = $exon_length_start + $exon_length_stop;
my $intron_delta = $branch_length + $acceptor_offset + 2 + $donor_length- $donor_offset;
my $intron_short_offset_forward = $donor_length - $donor_offset;
my $intron_short_offset_reverse = $branch_length + $acceptor_offset + 2 + $intron_short_length;
my $branch_offset = $branch_length + $acceptor_offset;

my $start_offset_begin = 0;
my $start_offset_end =  7;
my $stop_offset_begin = ($fixed_stop_offset);
my $donor_offset_begin = ($fixed_donor_offset);
my $acceptor_offset_end = ($exon_length_acceptor);

my $id = 0;
while (my $seq = <STDIN>)
{
    my $name;
    my @sequence;
print STDERR $seq."\n";

    if(!($seq =~ /(<.+>),(.+)/)){
      die "ERROR: tops_to_gtf.pl: invalid entry format\n";
    }
    my $id = $1;
    my $seq2 = $2;
    $id =~ /<(.+):(\d+),(\d+)/;
    $name = $1;
    $first_position = ($2-1);
    my @seqentry = split(/:\t/, $seq2);
    my $str = $seqentry[1];
    @sequence = split(/ /, $str);
    my $current_state = $sequence[0];
    my $begin = 0;
    push @sequence, "+_+_+_+_+_+_+";
    my $frame = 0;
    my $length = 0;
    $id = 0;


    for(my $i = 1; $i < scalar(@sequence); $i++)
    {
        if(!($sequence[$i] eq $current_state)){
            my $nameid = "MYOP".".$name".".$id";
            if($current_state eq "start") {
                $frame  = 0;
                print "$name\tmyop\tstart_codon\t".($begin + $start_offset_begin+1 + $first_position)."\t".($begin + $start_offset_begin+3+ $first_position)."\t.\t+\t".$frame."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n";
            } elsif( $current_state eq "stop") {
                print "$name\tmyop\tstop_codon\t".($begin + $stop_offset_begin+1+ $first_position)."\t".($begin + $stop_offset_begin+3+ $first_position)."\t.\t+\t".((3-$frame)%3)."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n\n";
                $id ++;
            } elsif ( $current_state =~ m/^CDS$/) {
                print "$name\tmyop\tCDS\t".($begin - $start_offset_end + 1+ $first_position)."\t".($i +  $first_position )."\t.\t+\t".$frame."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n";
            }elsif ($current_state =~ m/^CDS(\d)/) {
                $frame = $1;
                print "$name\tmyop\tCDS\t".($begin + $first_position + 1)."\t".($i +  $first_position )."\t.\t+\t".$frame."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n";
            }  elsif($current_state =~m/rstop/) {
                $i = process_reverse_strand($name, \@sequence, $begin,  \$id);
            }
            $begin = $i;
            $current_state = $sequence[$i];
        }
    }
}


sub process_reverse_strand {
    my $name = shift;
    my $seq = shift;
    my $i =  shift;
    my $idref = shift;
    my @sequence = @$seq;
    my $begin_of_gene = scalar(@sequence)-1;

    for(my $j = $i; $j < scalar(@sequence); $j++)
    {
        if(!($sequence[$j] =~ /^r/))
        {
            $begin_of_gene = $j;
            last;
        }
    }

    my $length = 0;
    my $frame = 0;
    my @output;
    my $begin = $begin_of_gene;
    my $current_state = $sequence[$begin];
    for(my $j = $begin_of_gene; $j >= $i-1; $j--)
    {
        if(!($sequence[$j] eq $current_state)) {
            my $nameid = "MYOP".".$name".".$$idref";
            if($current_state eq "rstart") {
                $frame = 0;
                push @output, "$name\tmyop\tstart_codon\t".($begin - $start_offset_begin - 1+ $first_position)."\t".($begin - $start_offset_begin +1+ $first_position )."\t.\t-\t".$frame."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n\n";
            } elsif( $current_state eq "rstop") {
                $frame = 0;
                $$idref ++;
                push @output, "$name\tmyop\tstop_codon\t".($begin - $stop_offset_begin-1+ $first_position)."\t".($begin - $stop_offset_begin+1+ $first_position)."\t.\t-\t".$frame."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n";
                $length = 0;
            } elsif ( $current_state =~ m/^rEI(\d)/) {
                $frame = 0; # first cds always have frame 0
                push @output, "$name\tmyop\tCDS\t".($j - $donor_offset_begin + 2+ $first_position)."\t".($begin + $start_offset_end + 1+ $first_position)."\t.\t-\t".$frame."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n";
                $length =  - ($j - $donor_offset_begin  + 2) + ($begin + $start_offset_end + 1 ) + 1;
            }elsif ($current_state =~ m/^rE(\d)(\d)/) {
                $frame = (3 - ($length - $frame) % 3)%3;
                push @output, "$name\tmyop\tCDS\t".($j - $donor_offset_begin +2+ $first_position )."\t".($begin + $acceptor_offset_end+1+ $first_position)."\t.\t-\t".$frame."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n";
                $length = - ($j - $donor_offset_begin+2 )+ ($begin + $acceptor_offset_end +1 )+1
            } elsif ($current_state =~ m/^rET(\d)/) {
                $frame = (3 - ($length - $frame) % 3)%3;
                push @output, "$name\tmyop\tCDS\t".($j - $stop_offset_begin + 2+ $first_position)."\t".($begin + $acceptor_offset_end + 1+ $first_position)."\t.\t-\t".$frame."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n";
        } elsif ($current_state =~ m/^rES/) {
            $frame = 0;
            push @output, "$name\tmyop\tCDS\t".($j - $stop_offset_begin + 2+ $first_position )."\t".($begin + $start_offset_end+1+ $first_position)."\t.\t-\t".$frame."\tgene_id \"$nameid\"; transcript_id \"$nameid\";\n";
        }
            $current_state = $sequence[$j];
            $begin = $j;
        }
    }
    for (my $k = scalar(@output)-1; $k >= 0; $k--)
    {
        print $output[$k];
    }
    return $begin_of_gene;
}

sub trim_spaces {
  my $v = shift;
  $v =~ s/^\s+//;     $v =~ s/\s+$//;
  return $v;
}
