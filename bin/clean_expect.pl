#!/usr/bin/perl -w
#
# replacement expect script for the old tcl one that hs_pred was using
#
# /usr/bin/expect
#set pdbfile [lindex $argv 0]
#puts "pdb: $pdbfile";
#
#spawn /webdata/binaries/current/hs_pred/clean.ex
#expect "Enter filename containing coordinates of structure\n"
#send "$pdbfile\r"
#interact
#

#use lib "/webdata/binaries/current/hs_pred/lib/perl5/site_perl/5.8.8/";
#use lib "/webdata/binaries/current/hs_pred/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";
#use lib "/SAN/bioinf/bioserv/binaries/current/hs_pred/lib/perl5/site_perl/5.8.8/";
#use lib "/SAN/bioinf/bioserv/binaries/current/hs_pred/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";

use strict;
use Expect;
use FindBin;
#chdir('/SAN/bioinf/bioserv/tmp/NewPredServer/') or die "$!";

my $file = $ARGV[0];

#print $file."\n";
my $command = "$FindBin::Bin/clean.ex";
my @params = [];
my $exp = Expect->spawn($command, @params) or die "Cannot spawn $command: $!\n";
#print $exp "zug\n";

if( $exp->expect(10, "Enter filename containing coordinates of structure"))
{
   $exp->send("$file\n");
}
else
{
	  print "hmm...";
}
 $exp->soft_close();
