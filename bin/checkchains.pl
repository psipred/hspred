#!/usr/bin/perl -w

use FileHandle;
use strict;
use English;
use Data::Dumper;
# exit(0);
my $file = $ARGV[0];
my $chains = $ARGV[1];
my $charmm = $ARGV[2];

my $aChains = [];

my $hResidues = read_charmm19();

@$aChains = split//, $chains;
my $hChains = {};
my $hChains_master = {};
my $chain_count = 0;
foreach my $chain(@$aChains)
{
	$hChains->{$chain} = 0;
	$hChains_master->{$chain} = 0;
	$chain_count++;
}

if($chain_count == 1)
{
	print "You must specify more than 1 chain\n";
	exit(70);
}

my $fhInput = new FileHandle($file, "r");

my $tmp_pdb = $file;
$tmp_pdb =~ s/\.pdb$/.tmp/;
$tmp_pdb =~ s/\.input$/.tmp/;
print $tmp_pdb."\n";
my $fhOut = new FileHandle($tmp_pdb,"w");

my $hRogue_atoms = {};
my $chains_found = 0;
my $rogue_residue_found = 0;
my $line_count =0;
my $model_found = 0;
while(my $line = $fhInput->getline)
{
	$line_count++;
	if($line =~ /^ATOM.{8}(.{4}).(.{3})\s(.)\s*\d+\s+.+$/)
	{
		my $aa_code = $2;
		my $chain = $3;
		my $atom = $1;
		$atom =~ s/\s//g;
		#print $chain."\n";
		if(exists $hChains_master->{$chain})
		{
			print $fhOut $line;
		}
		if(exists $hChains->{$chain} && $hChains->{$chain} == 0)
		{
			delete $hChains->{$chain};
			$chains_found++;
		}
		if(! exists $hResidues->{$aa_code})
		{
			$hRogue_atoms->{$aa_code} = $line_count;
			$rogue_residue_found++;
		}
	}
	elsif($line =~ /^MODEL/)
	{
		$model_found++;
	}
	else
	{
		print $fhOut $line;
	}
}

if($chains_found == $chain_count && $rogue_residue_found == 0 && $model_found == 0)
{
	print "All Chains Found\n";
	`mv $tmp_pdb $file`;
	exit;
}
else
{
	`rm $tmp_pdb`;
	if($chains_found != $chain_count)
	{
		print "The following specified chains were not found in the provided pdb file: ";
		foreach my $chain (keys %$hChains)
		{
			print $chain." ";
		}
		print "\n";
	}
	if($rogue_residue_found != 0)
	{
		print "Input PDB file contained residue types not represented in charmm, please use the standard 3 letter abbreviations to correct the following: ";
		foreach my $res (keys %$hRogue_atoms)
		{
			print $res." ";
		}
		print "\n";
	}
	if($model_found != 0)
	{
		print "HSPred can not process NMR ensembles, please submit a pdb file containing a single structure.";
		print "\n";
	}

	exit(80);
}


sub read_charmm19
{
	my $fhIn = new FileHandle($charmm,"r");
	#my $fhIn = new FileHandle("/cs/research/bioinf/home1/green/dbuchan/charmm19_ha.aaa","r");

	my $hData = {};
	while(my $line = $fhIn->getline)
	{
		if($line =~ /(\S+)\s+\S+\s+.+\s+.+/)
		{
			$hData->{$1} = 0;
		}
	}
	return($hData);
}
