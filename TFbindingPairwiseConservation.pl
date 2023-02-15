#!/usr/local/bin/perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::AlignIO;
use Bio::EnsEMBL::Utils::Exception qw(throw);

# input and output
my $quer_name;
my $quer_bed;
my $ref_name;
my $out;

# global variables
#my $reg_conf;
my $genomicalign_adaptor;
my $genomedb_adaptor;
#my $method_link_species_set_adaptor;


######################## usage and options ###############################
my $USAGE=("...fill in with description... .\nValid species names are Ensembl database species names (eg. homo_sapiens).\nTakes into account pairwise alignments using the EPO alignment.\n\nRun with following options:\n./TEflankingAlignment.pl -quer_name query_species -quer_bed query.bed -ref_name reference_species -out out_file\n\n");

GetOptions('quer_name=s' => \$quer_name, 'quer_bed=s' => \$quer_bed,
	'ref_name=s' => \$ref_name, 'out=s' => \$out);

die $USAGE unless ($quer_name and $quer_bed
	and $ref_name and $out);





########################### main program ##########################


my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous',
    -db_version => 98
);
$registry->set_reconnect_when_lost(1);

# get the MethodLinkSpeciesSet adaptor of Compara
my $method_link_species_set_adaptor = $registry->get_adaptor("Multi", "Compara", "MethodLinkSpeciesSet");

my $methodLinkSpeciesSet = $method_link_species_set_adaptor->fetch_by_method_link_type_species_set_name("EPO", "mammals");


# get species *core* Adaptor for Slices
my $slice_adaptor = $registry->get_adaptor(
       $quer_name, "Core", "Slice");

 # get the GenomicAlignBlock adaptor for Compara database
my $genomic_align_block_adaptor = $registry->get_adaptor(
    "Multi", "Compara", "GenomicAlignBlock");

 # (?) get the GenomeDB adaptor of Compara database
$genomedb_adaptor = $registry->get_adaptor("Multi", "Compara", "GenomeDB");



print "Projecting query species: ", $quer_name, "\tto reference: ",
	 $ref_name, "\n\tquery file: ", $quer_bed, "\n";

open (OUT, ">$out") or die "Cannot make outfile: $out\n";
print "\tmaking output: $out\n\n";

#print OUT "#", $ref_name,"_projected_chr\t", $ref_name, "_projected_start\t",
#$ref_name, "_projected_end\t", $quer_name, "_position\n";

open (IN, $quer_bed) or die "cannot open $quer_bed\n";


while (<IN>) {
    chomp;
    my @col = split(/\t/,$_);
    my $query_slice = $slice_adaptor->fetch_by_region(
        "toplevel", $col[0], $col[1], $col[2]);
    throw("No Slice can be created with coordinates $col[0]:$col[1]-$col[2]") if
        (!$query_slice);

    # The fetch_all_by_MethodLinkSpeciesSet_Slice() returns a ref.
    # to an array of GenomicAlingBlock objects
    my $all_genomic_align_blocks = $genomic_align_block_adaptor->
        fetch_all_by_MethodLinkSpeciesSet_Slice(
            $methodLinkSpeciesSet, $query_slice);
    #print scalar $all_genomic_align_blocks;

    if ( scalar @$all_genomic_align_blocks == 1 ) { #if there's only 1 alignment block corresponding to that slice
        foreach my $block (@$all_genomic_align_blocks) {
            next if (!defined $block);

            my $restricted = $block->restrict_between_reference_positions($query_slice->start(), $query_slice->end(), undef, 1);
            next if (!defined $restricted);

			my @genomic_aligns = @{ $restricted->get_all_GenomicAligns };
    		foreach my $align (@genomic_aligns) {
				next if (!$align);
				if ($align->genome_db()->name() eq $ref_name) {
            		print OUT $align->dnafrag->name(), "\t", # ref chr
            		$align->dnafrag_start(), "\t", # ref start
            		$align->dnafrag_end(), "\t", # ref end
            		$query_slice->seq_region_name(), "\t", # query chr
            		$query_slice->start(), "\t", # query start
            		$query_slice->end(), "\t", # query end
		            $col[3], "\t", # num of CpGs
		            $col[4], "\t", # peak summit
		            $col[5], "\t", # foldE
		            $col[6], "\t", # pval
		            $col[7], "\n"; # avg Meth
				}
			}
        }
    }
}

close(OUT);
close(IN);
