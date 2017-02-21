#!/usr/bin/perl -w
#
# Takes a pdb file and associated hspred scores and spits out 2 files for either
# side of the calculated interface, also edits the temperature factor with the
# hs-pred scores for visualisation.
#

use FileHandle;
use strict;
use English;
use Data::Dumper;

my $id = $ARGV[0];
my $chainsA = $ARGV[1];
my $chainsB = $ARGV[2];
my $tmpDir = $ARGV[3];
my $pdb = '';
if(-e $tmpDir."/".$id.".pdb")
{
	$pdb = $tmpDir."/".$id.".pdb";
}
elsif(-e $tmpDir."/".$id.".input")
{
	$pdb = $tmpDir."/".$id.".input";
}
else
{
	die "can not find pdb file";
}

my $hspred = $tmpDir."/".$id."_hs-pred.out";
my $output_1 = $tmpDir."/".$id."_initial.pdb";
my $output_2 = $tmpDir."/".$id."_second.pdb";

my $hChainsA = {};
my $hChainsB = {};
my $aEntries = [];

@$aEntries = split //, $chainsA;
foreach my $entry (@$aEntries)
{
	$hChainsA->{$entry} = 1;
}
@$aEntries = split //, $chainsB;
foreach my $entry (@$aEntries)
{
	$hChainsB->{$entry} = 1;
}

my $hHspred = {};
my $fHspred = new FileHandle($hspred, "r");
while(my $line = $fHspred->getline)
{
	chomp $line;
	@$aEntries = split /\s+/, $line;
	@$aEntries[2] =~ s/^-.+/0.000/;
	@$aEntries[2] =~ s/^(.{4})./ $1/;

	$hHspred->{@$aEntries[0]}{RESIDUE} = @$aEntries[1];
	$hHspred->{@$aEntries[0]}{SCORE} = @$aEntries[2];

}
#print Dumper $hHspred;
$fHspred->close;

my $fhPDB = new FileHandle($pdb,"r");
my $fhOutput1 = new FileHandle($output_1,"w");
my $fhOutput2 = new FileHandle($output_2,"w");
my $end_found = 0;
while(my $line = $fhPDB->getline)
{
	if($line =~ /^HEADER|^EXPDATA|^REMARK/)
	{
		print $fhOutput1 $line;
		print $fhOutput2 $line;
	}
	if($line =~ /^ATOM.{17}(.)\s*(\d+)/)
	{
		my $chain = $1;
		my $res = $2;
		if(exists $hChainsA->{$chain})
		{
			if(exists $hHspred->{$chain.$res})
			{
				if($hHspred->{$chain.$res}{SCORE} =~ /\s0\.00/)
				{
					#$line =~ s/(.{61}).{5}(.+)/$1 0.50$2/;
					#we handle regular pdb files with all 80 characters
					#and we handle models from Hex (66 characters)
					if($line =~ /.{80}/)
					{
						$line =~ s/(.{61}).{5}(.+)/$1   50$2/;
					}
					elsif($line =~ /.{66}/)
					{
						$line =~ s/(.{61}).{5}/$1   50/;
					}
				}
				else
				{
					if($line =~ /.{80}/)
					{
						$line =~ s/(.{61}).{5}(.+)/$1  100$2/;
					}
					elsif($line =~ /.{66}/)
					{
						$line =~ s/(.{61}).{5}/$1  100/;
					}
					#$line =~ s/(.{61}).{5}(.+)/$1 1.00$2/;
				}
			}
			else
			{
				if($line =~ /.{80}/)
				{
					$line =~ s/(.{61}).{5}(.{14})/$1 0.00$2/;
				}
				elsif($line =~ /.{66}/)
				{
					$line =~ s/(.{61}).{5}/$1 0.00/;
				}
			}
			if($line =~ /.{80}/)
			{
				$line =~ s/(.{60}).(.{19})/$1 $2/;
			}
			elsif($line =~ /.{66}/)
			{
				$line =~ s/(.{60}).(.{5})/$1 $2/;
			}
			print $fhOutput1 $line;
		}
		if(exists $hChainsB->{$chain})
		{
			#if(exists $hHspred->{$chain.$res})
			#{
			#	if($hHspred->{$chain.$res}{SCORE} =~ /\s0\.00/)
			#	{
			#		#$line =~ s/(.{61}).{5}(.+)/$1 0.50$2/;
			#		$line =~ s/(.{61}).{5}(.+)/$1   50$2/;
			#	}
			#	else
			#	{
			#		#$line =~ s/(.{61}).{5}(.+)/$1 1.00$2/;
			#		$line =~ s/(.{61}).{5}(.+)/$1  100$2/;
			#	}
			#}
			#else
			#{
			#	$line =~ s/(.{61}).{5}(.{14})/$1 0.00$2/;
			#}
			#
			#$line =~ s/(.{60}).(.{19})/$1 $2/;
			if(exists $hHspred->{$chain.$res})
			{
				if($hHspred->{$chain.$res}{SCORE} =~ /\s0\.00/)
				{
					#$line =~ s/(.{61}).{5}(.+)/$1 0.50$2/;
					#we handle regular pdb files with all 80 characters
					#and we handle models from Hex (66 characters)
					if($line =~ /.{80}/)
					{
						$line =~ s/(.{61}).{5}(.+)/$1   50$2/;
					}
					elsif($line =~ /.{66}/)
					{
						$line =~ s/(.{61}).{5}/$1   50/;
					}
				}
				else
				{
					if($line =~ /.{80}/)
					{
						$line =~ s/(.{61}).{5}(.+)/$1  100$2/;
					}
					elsif($line =~ /.{66}/)
					{
						$line =~ s/(.{61}).{5}/$1  100/;
					}
					#$line =~ s/(.{61}).{5}(.+)/$1 1.00$2/;
				}
			}
			else
			{
				if($line =~ /.{80}/)
				{
					$line =~ s/(.{61}).{5}(.{14})/$1 0.00$2/;
				}
				elsif($line =~ /.{66}/)
				{
					$line =~ s/(.{61}).{5}/$1 0.00/;
				}
			}
			if($line =~ /.{80}/)
			{
				$line =~ s/(.{60}).(.{19})/$1 $2/;
			}
			elsif($line =~ /.{66}/)
			{
				$line =~ s/(.{60}).(.{5})/$1 $2/;
			}
			print $fhOutput2 $line;
		}

	}
	if($line =~ /^HETATM.{15}(.)/)
	{
		my $chain = $1;
		$line =~ s/(.{61}).{5}(.{14})/$1 0.00$2/;
		if(exists $hChainsA->{$chain})
		{
			print $fhOutput1 $line;
		}
		if(exists $hChainsB->{$chain})
		{
			print $fhOutput2 $line;
		}
	}

	if($line =~ /^END/)
	{
		print $fhOutput1 $line;
		print $fhOutput2 $line;
		$end_found++;
	}

}

if($end_found == 0)
{
		print $fhOutput1 "END\n";
		print $fhOutput2 "END\n";
}

$fhPDB->close;
$fhOutput1->close;
$fhOutput2->close;
