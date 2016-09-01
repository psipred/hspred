#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my (@fields,$i,$j,$k,$x,$y);
my ($pdb_code,$chains_1,$chains_2,$res_id,$aa_tmp);
my (@resnum_1,@restype_1,@xyz_1,@atom_1,@xyz_1_cb,@atom_pdbname_1,@sidechain_1,@mainchain_1);
my (@resnum_2,@restype_2,@xyz_2,@atom_2,@xyz_2_cb,@atom_pdbname_2,@sidechain_2,@mainchain_2);
my (%resseq_1,%nter_1,%cter_1,%resseq_2,%nter_2,%cter_2,@int_1,@int_2);
my (@epsilon,@r_min,@dgref,@dgfree,@fracvol,@corrlen,@vdwrad);
my (@atom_type_1,@atom_charge_1,@atom_type_2,@atom_charge_2,%hb_energy);
my (@model,@mean,@stdev,@ddg,@predictions);
my ($format);
my (@conn_mat_1,@conn_mat_2);

#
# Input name of (pdb) file to be processed.
# The chains of the two interacting monomers need also to be specified
#
# Added changes to make it portable and thread safe for the server backend nodes
#
#
use constant PI             => 4 * atan2 1, 1;

my %aa_code = (
    "ALA" => "A", "ARG" => "R", "ASN" => "N", "ASP" => "D", "CYS" => "C",
    "GLN" => "Q", "GLU" => "E", "GLY" => "G", "HIS" => "H", "ILE" => "I",
    "LEU" => "L", "LYS" => "K", "MET" => "M", "PHE" => "F", "PRO" => "P",
    "SER" => "S", "THR" => "T", "TRP" => "W", "TYR" => "Y", "VAL" => "V");

$pdb_code = $ARGV[0];
$chains_1 = $ARGV[1];
$chains_2 = $ARGV[2];
my $data_dir = $ARGV[3];
my $tmp_dir = $ARGV[4];
#open INPUT_FILE,"input_1.dat" or die "input file not found ..." ;
#while (<INPUT_FILE>) {
#    chomp;
#    if (/^(\w+)\((\w+):(\w+)\)/){
#	$pdb_code = $1;
#	$chains_1 = $2;
#	$chains_2 = $3;
#    } else {
#	die;
#    }
#}
#close INPUT_FILE;

&models;

&pdb_parser;

&interface_residues;

&force_fields;

&h_bonds;

&local_atomic_connectivity;

@predictions=();
foreach $x (@int_1){
    $res_id = $resnum_1[$x];    
    if ($restype_1[$x] =~ /([A-Z]{3})$/){
	$aa_tmp = $aa_code{$1};
    } else {
	die;
    }
    if ($aa_tmp eq 'E'){
	@ddg=();
	push(@ddg,&hb_sc(1,$x));
	push(@ddg,&hb_frag(1,$x));
	push(@ddg,&vdw_intra_sc(1,$x));
	&norm_features($aa_tmp);
	system 
	    "$data_dir/svm_classify -v 0  $tmp_dir/$pdb_code"."_input_file ${model[1]}.dat $tmp_dir/".$pdb_code."_output_file";
	open OUT_FILE, $tmp_dir."/".$pdb_code."_output_file" or die "minchia ...";
	while (<OUT_FILE>) {
	    chomp;
	    push(@predictions,[($res_id,$aa_tmp,$_)]);
	}
	close OUT_FILE;
    } elsif ($aa_tmp eq 'R'){
	@ddg=();
	push(@ddg,&vdw_sc(1,$x));
	push(@ddg,&hb_frag(1,$x));
	push(@ddg,&vdw_intra_sc(1,$x));
	push(@ddg,&coul_intra_sc(1,$x));
	&norm_features($aa_tmp);
	system 
	    "$data_dir/svm_classify -v 0  $tmp_dir/$pdb_code"."_input_file ${model[2]}.dat  $tmp_dir/".$pdb_code."_output_file";
	open OUT_FILE, $tmp_dir."/".$pdb_code."_output_file" or die "minchia ...";
	while (<OUT_FILE>) {
	    chomp;
	    push(@predictions,[($res_id,$aa_tmp,$_)]);
	}
	close OUT_FILE;
    } elsif ('CDFHIKLMNQSTVWY' =~ /$aa_tmp/){
	@ddg=();
	push(@ddg,&vdw_sc(1,$x));
	push(@ddg,&hb_sc(1,$x));
	push(@ddg,&solv_sc(1,$x));
	push(@ddg,&vdw_frag(1,$x));
	push(@ddg,&hb_frag(1,$x));
	push(@ddg,&solv_frag(1,$x));
	push(@ddg,&vdw_intra_sc(1,$x));
	&norm_features($aa_tmp);
	system 
	    "$data_dir/svm_classify -v 0  $tmp_dir/$pdb_code"."_input_file ${model[0]}.dat $tmp_dir/".$pdb_code."_output_file";
	open OUT_FILE,  $tmp_dir."/".$pdb_code."_output_file" or die "minchia ...";
	while (<OUT_FILE>) {
	    chomp;
	    push(@predictions,[($res_id,$aa_tmp,$_)]);
	}
	close OUT_FILE;
    }
}
foreach $y (@int_2){
    $res_id = $resnum_2[$y];    
    if ($restype_2[$y] =~ /([A-Z]{3})$/){
	$aa_tmp = $aa_code{$1};
    } else {
	die;
    }
    if ($aa_tmp eq 'E'){
	@ddg=();
	push(@ddg,&hb_sc(2,$y));
	push(@ddg,&hb_frag(2,$y));
	push(@ddg,&vdw_intra_sc(2,$y));
	&norm_features($aa_tmp);
	system 
	    "$data_dir/svm_classify -v 0 $tmp_dir/$pdb_code"."_input_file ${model[1]}.dat  $tmp_dir/".$pdb_code."_output_file";
	open OUT_FILE, $tmp_dir."/".$pdb_code."_output_file" or die "minchia ...";
	while (<OUT_FILE>) {
	    chomp;
	    push(@predictions,[($res_id,$aa_tmp,$_)]);
	}
	close OUT_FILE;
    } elsif ($aa_tmp eq 'R'){
	@ddg=();
	push(@ddg,&vdw_sc(2,$y));
	push(@ddg,&hb_frag(2,$y));
	push(@ddg,&vdw_intra_sc(2,$y));
	push(@ddg,&coul_intra_sc(2,$y));
	&norm_features($aa_tmp);
	system 
	    "$data_dir/svm_classify -v 0 $tmp_dir/$pdb_code"."_input_file ${model[2]}.dat $tmp_dir/".$pdb_code."_output_file";
	open OUT_FILE, $tmp_dir."/".$pdb_code."_output_file" or die "minchia ...";
	while (<OUT_FILE>) {
	    chomp;
	    push(@predictions,[($res_id,$aa_tmp,$_)]);
	}
	close OUT_FILE;
    } elsif ('CDFHIKLMNQSTVWY' =~ /$aa_tmp/){
	@ddg=();
	push(@ddg,&vdw_sc(2,$y));
	push(@ddg,&hb_sc(2,$y));
	push(@ddg,&solv_sc(2,$y));
	push(@ddg,&vdw_frag(2,$y));
	push(@ddg,&hb_frag(2,$y));
	push(@ddg,&solv_frag(2,$y));
	push(@ddg,&vdw_intra_sc(2,$y));
	&norm_features($aa_tmp);
	system 
	    "$data_dir/svm_classify -v 0 $tmp_dir/$pdb_code"."_input_file ${model[0]}.dat $tmp_dir/".$pdb_code."_output_file";
	open OUT_FILE, $tmp_dir."/".$pdb_code."_output_file" or die "minchia ...";
	while (<OUT_FILE>) {
	    chomp;
	    push(@predictions,[($res_id,$aa_tmp,$_)]);
	}
	close OUT_FILE;
    }
}

open OUT_FILE, ">$tmp_dir/$pdb_code"."_hs-pred.out" or die "minchia ...";
foreach $i (0..$#predictions){
    printf OUT_FILE "%-5s %1s %7.3f\n", 
    $predictions[$i][0],$predictions[$i][1],$predictions[$i][2];
}
close OUT_FILE;

##########################################################################

sub norm_features{
    my $aa_tmp=$_[0];
    my ($format,$i,$val,@input_feat,$i_f,$i_aa);
    
    if ('CDFHIKLMNQSTVWY' =~ /$aa_tmp/){
	$i_aa =0;
    } elsif ($aa_tmp eq 'E'){
	$i_aa =1;
    } elsif ($aa_tmp eq 'R'){
	$i_aa =2;
    }

    $format = "%2d";
    foreach $i (0..$#ddg){
	$val = ($ddg[$i]-$mean[$i_aa][$i])/$stdev[$i_aa][$i];
	push(@input_feat,$val);
	$i_f=$i+1;
	$format .= " $i_f:%-7.3f";
    }
    open TEST_FILE, ">$tmp_dir/$pdb_code"."_input_file" or die "minchia ..." ;
    printf TEST_FILE "$format\n", 0, @input_feat;
    close TEST_FILE;
}


###########################################################################

sub pdb_parser{
    my ($tmp_1a,$tmp_1b,$tmp_2,$tmp_2a,$tmp_3,$tmp_4,$tmp_5,$tmp_6,$tmp_7,$tmp_8,$tmp_9);

    $tmp_1a=$tmp_1b=-1;
    @resnum_1=@restype_1=@xyz_1=@atom_1=@xyz_1_cb=@atom_pdbname_1=@sidechain_1=@mainchain_1=();
    @resnum_2=@restype_2=@xyz_2=@atom_2=@xyz_2_cb=@atom_pdbname_2=@sidechain_2=@mainchain_2=();
    %resseq_1=%nter_1=%cter_1=%resseq_2=%nter_2=%cter_2=();    
    open PDB_FILE, $tmp_dir."/${pdb_code}.pdb" or die "pdb file not found ..." ;
    while (<PDB_FILE>) {
	chomp;
	if (/^ATOM/){
	    $tmp_2=substr($_ ,12,4);          # 2 Atom name
	    $tmp_2a=substr($_ ,13,1);         # 2a Atom (H,C,O,...)
	    if ($tmp_2 =~ /^H/){$tmp_2a = 'H'}
	    $tmp_3=substr($_ ,16,1);          # 3 Alternate location indicator
	    $tmp_4=substr($_ ,17,3);          # 4 Residue name.
	    $tmp_5=substr($_ ,21,1);          # 5 Chain identifier
	    $tmp_6=substr($_ ,22,5);          # 6 Residue sequence number
	    $tmp_7=substr($_ ,30,8);          # 7 X coordinate
	    $tmp_8=substr($_ ,38,8);          # 8 Y coordinate
	    $tmp_9=substr($_ ,46,8);          # 9 Z coordinate
	    $tmp_2 =~ s/^\s+//; $tmp_2 =~ s/\s+$// ;
	    $tmp_2a =~ s/^\s+//; $tmp_2a =~ s/\s+$// ;
	    $tmp_3 =~ s/^\s+//; $tmp_3 =~ s/\s+$// ;
	    $tmp_4 =~ s/^\s+//; $tmp_4 =~ s/\s+$// ;
	    $tmp_5 =~ s/^\s+//; $tmp_5 =~ s/\s+$// ;
	    $tmp_6 =~ s/^\s+//; $tmp_6 =~ s/\s+$// ;
	    $tmp_7 =~ s/^\s+//; $tmp_7 =~ s/\s+$// ;
	    $tmp_8 =~ s/^\s+//; $tmp_8 =~ s/\s+$// ;
	    $tmp_9 =~ s/^\s+//; $tmp_9 =~ s/\s+$// ;
	    if (($tmp_2a) && ($tmp_2a ne 'H') && (!($tmp_3) || ($tmp_3 eq 'A'))){
		if (($tmp_4 eq 'ILE') && ($tmp_2 eq 'CD')){
		    $tmp_2 = 'CD1';		
		}
		$res_id="$tmp_5"."$tmp_6";
		if ($chains_1 =~ /$tmp_5/){
		    $tmp_1a += 1;
		    if (($#resnum_1 == -1) || ($resnum_1[-1] ne "$res_id")) {
			push (@resnum_1,$res_id);
			push (@restype_1,$tmp_4);
			$resseq_1{$res_id}=$#resnum_1;
			if ( ! exists  $nter_1{$tmp_5} ){
			    $nter_1{$tmp_5}= $#resnum_1;
			}
			$cter_1{$tmp_5}= $#resnum_1;
		    }
		    push(@atom_pdbname_1,$tmp_2);
		    push(@xyz_1,[($tmp_7,$tmp_8,$tmp_9)]);
		    push (@{$atom_1[$#resnum_1]}, $tmp_1a);
		    if (! defined($xyz_1_cb[$#resnum_1])){
			$xyz_1_cb[$#resnum_1]=[($tmp_7,$tmp_8,$tmp_9)];
		    }
		    if (("$tmp_2" eq 'CB') || 
			(("$tmp_4" eq 'GLY') && ("$tmp_2" eq 'CA'))){
			$xyz_1_cb[$#resnum_1]=[($tmp_7,$tmp_8,$tmp_9)];
		    }
		    if (($tmp_2 eq 'N') || ($tmp_2 eq 'CA') || ($tmp_2 eq 'CB') || 
			($tmp_2 eq 'C') || ($tmp_2 eq 'O')){
			push(@{$mainchain_1[$#resnum_1]}, $tmp_1a);
		    } else {
			push(@{$sidechain_1[$#resnum_1]}, $tmp_1a);
		    }
		} elsif ($chains_2 =~ /$tmp_5/){
		    $tmp_1b += 1;
		    if (($#resnum_2 == -1) || ($resnum_2[-1] ne "$res_id")) {
			push (@resnum_2,$res_id);
			push (@restype_2,$tmp_4);
			$resseq_2{$res_id}=$#resnum_2;
			if ( ! exists  $nter_2{$tmp_5} ){
			    $nter_2{$tmp_5}= $#resnum_2;
			}
			$cter_2{$tmp_5}= $#resnum_2;
		    }
		    push(@atom_pdbname_2,$tmp_2);
		    push(@xyz_2,[($tmp_7,$tmp_8,$tmp_9)]);
		    push (@{$atom_2[$#resnum_2]}, $tmp_1b);
		    if (! defined($xyz_2_cb[$#resnum_2])){
			$xyz_2_cb[$#resnum_2]=[($tmp_7,$tmp_8,$tmp_9)];
		    }
		    if (("$tmp_2" eq 'CB') || 
			(("$tmp_4" eq 'GLY') && ("$tmp_2" eq 'CA'))){
			$xyz_2_cb[$#resnum_2]=[($tmp_7,$tmp_8,$tmp_9)];
		    }
		    if (($tmp_2 eq 'N') || ($tmp_2 eq 'CA') || ($tmp_2 eq 'CB') || 
			($tmp_2 eq 'C') || ($tmp_2 eq 'O')){
			push(@{$mainchain_2[$#resnum_2]}, $tmp_1b);
		    } else {  
			push(@{$sidechain_2[$#resnum_2]}, $tmp_1b);
		    }
		}
	    }
	}
    }
    close PDB_FILE;

    foreach $x (values %nter_1){
	$restype_1[$x] = 'N' . $restype_1[$x] ;
    }
    foreach $x (values %cter_1){
	$restype_1[$x] = 'C' . $restype_1[$x] ;
    }
    foreach $x (values %nter_2){
	$restype_2[$x] = 'N' . $restype_2[$x] ;
    }
    foreach $x (values %cter_2){
	$restype_2[$x] = 'C' . $restype_2[$x] ;
    }

}

##############################################################################

sub local_atomic_connectivity{
    my ($i,$j,$aa_tmp,$dist,$atm_i,$atm_j,%conn_res,%conn_peptbond,%conn_peptbond_pro);

    %conn_res=();
    open CONN_FILE, $data_dir."/conn_mat.txt" or die "minchia ..." ;
    while (<CONN_FILE>){
	chomp;
	@fields = split;
	if (shift(@fields) =~/^([A-Z]{3})_([0-9]):/){
	    $aa_tmp = $1;
	    $dist = $2;
	    foreach (@fields){
		if (/^(\w+)-(\w+)$/){
		    $conn_res{$aa_tmp}{$1}{$2}=$dist;
		    $conn_res{$aa_tmp}{$2}{$1}=$dist;
		}
	    }
	}
    }
    close CONN_FILE;

    $conn_peptbond{'C'}{'N'}=1;
    $conn_peptbond{'CA'}{'N'}=$conn_peptbond{'C'}{'CA'}=$conn_peptbond{'O'}{'N'}=2;
    $conn_peptbond{'N'}{'N'}=$conn_peptbond{'CB'}{'N'}=$conn_peptbond{'CA'}{'CA'}=3;
    $conn_peptbond{'O'}{'CA'}=$conn_peptbond{'C'}{'CB'}=$conn_peptbond{'C'}{'C'}=3;
    $conn_peptbond_pro{'C'}{'CD'}=2; 
    $conn_peptbond_pro{'CA'}{'CD'}=$conn_peptbond_pro{'O'}{'CD'}=$conn_peptbond_pro{'C'}{'CG'}=3;
    
    @conn_mat_1=();
    foreach $x (0..$#resnum_1) {
	if ($restype_1[$x] =~ /([A-Z]{3})$/){
	    $aa_tmp = $1;
	}
	foreach $i (@{$atom_1[$x]}){
	    foreach $j (@{$atom_1[$x]}){
		$atm_i=$atom_pdbname_1[$i];
		$atm_j=$atom_pdbname_1[$j];
		if (exists $conn_res{$aa_tmp}{$atm_i}{$atm_j}) {
		    $conn_mat_1[$i][$j]= $conn_res{$aa_tmp}{$atm_i}{$atm_j};
		    $conn_mat_1[$j][$i]= $conn_res{$aa_tmp}{$atm_i}{$atm_j};
		}
	    }
	}
    }
    foreach $x (0..$#resnum_1-1) {
	$dist=10;
	foreach $i (@{$atom_1[$x]}){
	    foreach $j (@{$atom_1[$x+1]}){
		$atm_i=$atom_pdbname_1[$i];
		$atm_j=$atom_pdbname_1[$j];
		if (($atm_i eq 'C') && ($atm_j eq 'N') ){
		    $dist = 0;
		    foreach $k (0..2){
			$dist += ($xyz_1[$i][$k] - $xyz_1[$j][$k])**2 ;
		    }
		    $dist=sqrt($dist);
		}
	    }
	}
	if ($dist < 1.5){
	    foreach $i (@{$atom_1[$x]}){
		foreach $j (@{$atom_1[$x+1]}){
		    $atm_i=$atom_pdbname_1[$i];
		    $atm_j=$atom_pdbname_1[$j];
		    if (exists $conn_peptbond{$atm_i}{$atm_j}){
			$conn_mat_1[$i][$j]=$conn_peptbond{$atm_i}{$atm_j};
			$conn_mat_1[$j][$i]=$conn_peptbond{$atm_j}{$atm_i};
		    }
		}
	    }
	    if ($restype_1[$x+1] eq 'PRO'){
		foreach $i (@{$atom_1[$x]}){
		    foreach $j (@{$atom_1[$x+1]}){
			$atm_i=$atom_pdbname_1[$i];
			$atm_j=$atom_pdbname_1[$j];
			if (exists $conn_peptbond_pro{$atm_i}{$atm_j}){
			    $conn_mat_1[$i][$j]=$conn_peptbond_pro{$atm_i}{$atm_j};
			    $conn_mat_1[$j][$i]=$conn_peptbond_pro{$atm_j}{$atm_i};
			}
		    }
		}
	    }
	}
    }
    foreach $x (values %cter_1){
	foreach $i (@{$atom_1[$x]}){
	    foreach $j (@{$atom_1[$x]}){
		$atm_i=$atom_pdbname_1[$i];
		$atm_j=$atom_pdbname_1[$j];
		if ($atm_j eq 'OXT') {
		    if ($atm_i eq 'OXT') {
			$conn_mat_1[$j][$j]=0;
		    } elsif ($atm_i eq 'C') {
			$conn_mat_1[$i][$j]=1;			
			$conn_mat_1[$j][$i]=1;
		    } elsif ($atm_i eq 'CA') {
			$conn_mat_1[$i][$j]=2;			
			$conn_mat_1[$j][$i]=2;
		    } elsif (($atm_i eq 'N') || ($atm_i eq 'CB')){
			$conn_mat_1[$i][$j]=3;			
			$conn_mat_1[$j][$i]=3;
		    }
		}
	    }
	}
    }

    foreach $i (0..$#atom_pdbname_1){
	foreach $j (0..$#atom_pdbname_1){
	    if (! defined($conn_mat_1[$i][$j])){
		$conn_mat_1[$i][$j]=4;
	    } 
	}	    
    }

    @conn_mat_2=();
    foreach $y (0..$#resnum_2) {
	if ($restype_2[$y] =~ /([A-Z]{3})$/){
	    $aa_tmp = $1;
	}
	foreach $i (@{$atom_2[$y]}){
	    foreach $j (@{$atom_2[$y]}){
		$atm_i=$atom_pdbname_2[$i];
		$atm_j=$atom_pdbname_2[$j];
		if (exists $conn_res{$aa_tmp}{$atm_i}{$atm_j}) {
		    $conn_mat_2[$i][$j]= $conn_res{$aa_tmp}{$atm_i}{$atm_j};
		    $conn_mat_2[$j][$i]= $conn_res{$aa_tmp}{$atm_i}{$atm_j};
		}
	    }
	}
    }
    foreach $y (0..$#resnum_2-1) {
	$dist=10;
	foreach $i (@{$atom_2[$y]}){
	    foreach $j (@{$atom_2[$y+1]}){
		$atm_i=$atom_pdbname_2[$i];
		$atm_j=$atom_pdbname_2[$j];
		if (($atm_i eq 'C') && ($atm_j eq 'N') ){
		    $dist = 0;
		    foreach $k (0..2){
			$dist += ($xyz_2[$i][$k] - $xyz_2[$j][$k])**2 ;
		    }
		    $dist=sqrt($dist);
		}
	    }
	}
	if ($dist < 1.5){
	    foreach $i (@{$atom_2[$y]}){
		foreach $j (@{$atom_2[$y+1]}){
		    $atm_i=$atom_pdbname_2[$i];
		    $atm_j=$atom_pdbname_2[$j];
		    if (exists $conn_peptbond{$atm_i}{$atm_j}){
			$conn_mat_2[$i][$j]=$conn_peptbond{$atm_i}{$atm_j};
			$conn_mat_2[$j][$i]=$conn_peptbond{$atm_j}{$atm_i};
		    }
		}
	    }
	    if ($restype_2[$y+1] eq 'PRO'){
		foreach $i (@{$atom_2[$y]}){
		    foreach $j (@{$atom_2[$y+1]}){
			$atm_i=$atom_pdbname_2[$i];
			$atm_j=$atom_pdbname_2[$j];
			if (exists $conn_peptbond_pro{$atm_i}{$atm_j}){
			    $conn_mat_2[$i][$j]=$conn_peptbond_pro{$atm_i}{$atm_j};
			    $conn_mat_2[$j][$i]=$conn_peptbond_pro{$atm_j}{$atm_i};
			}
		    }
		}
	    }
	}

    }
    foreach $y (values %cter_2){
	foreach $i (@{$atom_2[$y]}){
	    foreach $j (@{$atom_2[$y]}){
		$atm_i=$atom_pdbname_2[$i];
		$atm_j=$atom_pdbname_2[$j];
		if ($atm_j eq 'OXT') {
		    if ($atm_i eq 'OXT') {
			$conn_mat_2[$j][$j]=0;
		    } elsif ($atm_i eq 'C') {
			$conn_mat_2[$i][$j]=1;			
			$conn_mat_2[$j][$i]=1;
		    } elsif ($atm_i eq 'CA') {
			$conn_mat_2[$i][$j]=2;			
			$conn_mat_2[$j][$i]=2;
		    } elsif (($atm_i eq 'N') || ($atm_i eq 'CB')){
			$conn_mat_2[$i][$j]=3;			
			$conn_mat_2[$j][$i]=3;
		    }
		}
	    }
	}
    }
    foreach $i (0..$#atom_pdbname_2){
	foreach $j (0..$#atom_pdbname_2){
	    if (! defined($conn_mat_2[$i][$j])){
		$conn_mat_2[$i][$j]=4;
	    } 
	}	    
    }

}


##############################################################################

sub force_fields{
    # reads force fields atomic parameters (e.g. charge,radius, etc) and assigns them
    # to individual atoms in the structure
    #
    my (%tbl_atom_par,%tbl_atom_type,%tbl_atom_charge);
    my ($aa_tmp,$atm_tmp,$x,$y,$i,$j);

    %tbl_atom_par=();
    @epsilon=@r_min=@dgref=@dgfree=@fracvol=@corrlen=@vdwrad=();
    open PAR_FILE, $data_dir."/charmm19.par" or 
	die "minchia ..." ;
    $i=-1;
    while (<PAR_FILE>) {
	chomp;
	@fields = split;
	$i++;
	$tbl_atom_par{$fields[0]} = $i ;
	$epsilon[$i] = $fields[1];
	$r_min[$i]   = $fields[2];
	$dgref[$i]   = $fields[3] ;
	$dgfree[$i]  = $fields[4] ;
	$fracvol[$i] = $fields[5] ;
	$corrlen[$i] = $fields[6] ; 	
	$vdwrad[$i]  = $fields[7] ; 
    }
    close PAR_FILE;
    
    %tbl_atom_type=%tbl_atom_charge=();
    open PAR_FILE, $data_dir."/charmm19_ha.aaa" or 
	die "minchia ..." ;
    while (<PAR_FILE>) {
		chomp;
		@fields = split;
		if ($fields[0] =~ /HIS/)
		{
	    	if (scalar(@fields) == 6)
			{
				$fields[3] =  ($fields[3]+$fields[5])/2 ;
	    	}
			else 
			{
				die;
	    	}
		}
		$tbl_atom_type{$fields[0]}{$fields[1]}   = $tbl_atom_par{$fields[2]} ;
		$tbl_atom_charge{$fields[0]}{$fields[1]} = $fields[3];
		
		#DB 21/03/12 here we add some C and N terminal atom types that sometimes crystalographers 
		#or modelling packages don't correctly label in the pdb file. For instance OXT is valid 
		#and present in charmm19 but in pdb files it is usually just accompanied with the LEU aa type
		#rather than the CLEU or NLEU. The "correct" way to handle this would be to change the 
		#pdb parser in this script to correctly re-label the C and N termini residues of chains
		#from the 3 letter code to the charmm19 C and N types (LEU -> CLEU). The problem appears to
		#be that the parser here doesn't correctly handle discontigous chain segments in pdb files.
		if($fields[0] =~ /^C|^N/)
		{
			my $tmp_res_type = $fields[0];
			$tmp_res_type =~ s/^C|^N//;	
			if(exists $tbl_atom_type{$tmp_res_type}{$fields[1]})
			{
				next;
			}
			else
			{
				$tbl_atom_type{$tmp_res_type}{$fields[1]}   = $tbl_atom_par{$fields[2]} ;
				$tbl_atom_charge{$tmp_res_type}{$fields[1]} = $fields[3];
			}
		}
		
    }
    close PAR_FILE;
	#print Dumper \%tbl_atom_type;
	#exit;
    @atom_type_1=@atom_charge_1=@atom_type_2=@atom_charge_2=();
    foreach $x (0..$#resnum_1) {
	$aa_tmp = $restype_1[$x];
	foreach $i (@{$atom_1[$x]}){
	    $atm_tmp=$atom_pdbname_1[$i];
	    if (exists $tbl_atom_type{$aa_tmp}{$atm_tmp}){
		$atom_type_1[$i]    = $tbl_atom_type{$aa_tmp}{$atm_tmp};
		$atom_charge_1[$i]  = $tbl_atom_charge{$aa_tmp}{$atm_tmp};
	    } else {
		print "$aa_tmp $atm_tmp   \n";
		die;
	    }
	}
    }		
    foreach $y (0..$#resnum_2) {
	$aa_tmp = $restype_2[$y];
	foreach $j (@{$atom_2[$y]}){
	    $atm_tmp=$atom_pdbname_2[$j];
	    if (exists $tbl_atom_type{$aa_tmp}{$atm_tmp}){
		$atom_type_2[$j]    = $tbl_atom_type{$aa_tmp}{$atm_tmp};
		$atom_charge_2[$j]  = $tbl_atom_charge{$aa_tmp}{$atm_tmp};
	    } else {
		print "$aa_tmp $atm_tmp   \n";
		die;
	    }
	}
    }		
}

#############################################################################

sub h_bonds{
    # step 1: parse hbplus file and identify h-bonds at the interface 
    # step 2: calculate energy associate to each interface h-bond
    my ($tmp_1a,$tmp_1b,$tmp_2a,$tmp_2b,$tmp_3a,$tmp_3b,$tmp_4a,$tmp_4b);
    my ($tmp_5,$tmp_6,$tmp_7,$tmp_8,$tmp_9,$tmp_10);
    my (%tbl_atom_type,%hb_list,$res_id_1,$res_id_2,$class_1,$class_2,$v_hbond_ij);
    my (@solvexp_1,@solvexp_2,$donor,$acceptor);

    %tbl_atom_type=();
    open PAR_FILE, $data_dir."/charmm19_ha.aaa" or 
	die "minchia ..." ;
    while (<PAR_FILE>) {
	chomp;
	@fields = split;
	if ($fields[0] =~ /HIS/){
	    if (scalar(@fields) == 6){
		$fields[3] =  ($fields[3]+$fields[5])/2 ;
	    } else {
		die;
	    }
	}
	$tbl_atom_type{$fields[0]}{$fields[1]}   = $fields[2] ;
    }
    close PAR_FILE;

    &solvexp(\@solvexp_1,\@solvexp_2) ;

    %hb_list=();
    open HB_FILE, $tmp_dir."/${pdb_code}.hb2" or die "file ${pdb_code}.hb2 not found ..." ;   
    foreach $j (1..8){
	$_=<HB_FILE>;
    }
    while (<HB_FILE>) {
	$tmp_1a=substr($_,0,1);                                                    # chain id D
	$tmp_1b=substr($_,14,1);                                                   # chain id A 
	if ((($chains_1 =~ /$tmp_1a/) && ($chains_2 =~ /$tmp_1b/)) 
	    || (($chains_1 =~ /$tmp_1b/) && ($chains_2 =~ /$tmp_1a/))){
	    $tmp_2a=substr($_,1,5);  $tmp_2a =~ s/^0+//; $tmp_2a =~ s/-$// ;       # res number D
	    $tmp_2b=substr($_,15,5); $tmp_2b =~ s/^0+//; $tmp_2b =~ s/-$// ;       # res number A 
	    $tmp_3a=substr($_,6,3);                                                # aa D
	    $tmp_3b=substr($_,20,3);                                               # aa A
	    $tmp_4a=substr($_,9,4);  $tmp_4a =~ s/^\s+//; $tmp_4a =~ s/\s+$// ;    # atom D
	    $tmp_4b=substr($_,23,4); $tmp_4b =~ s/^\s+//; $tmp_4b =~ s/\s+$// ;    # atom A
	    if ((exists $tbl_atom_type{$tmp_3a}{$tmp_4a}) && 
		(exists $tbl_atom_type{$tmp_3b}{$tmp_4b})){
		$donor    = "$tmp_1a"."$tmp_2a"."$tmp_4a" ;
		$acceptor = "$tmp_1b"."$tmp_2b"."$tmp_4b" ;
		$tmp_5=substr($_,27,5);  $tmp_5 =~ s/^\s+//;  $tmp_5 =~ s/\s+$// ; # DA dist
		$tmp_6=substr($_,33,2);                                            # hb_type
		$tmp_7=substr($_,46,5);  $tmp_7 =~ s/^\s+//;  $tmp_7 =~ s/\s+$// ; # DHA angle
		$tmp_8=substr($_,52,5);  $tmp_8 =~ s/^\s+//;  $tmp_8 =~ s/\s+$// ; # H-A dist 
		$tmp_9=substr($_,58,5);  $tmp_9 =~ s/^\s+//;  $tmp_9 =~ s/\s+$// ; # H-A-AA angle
		$tmp_10=substr($_,64,5); $tmp_10 =~ s/^\s+//; $tmp_10 =~ s/\s+$// ;# D-A-AA angle
		
		$hb_list{$donor}{$acceptor} = 1 ;
		$hb_list{$acceptor}{$donor} = 1 ;

		if ($chains_1 =~ /$tmp_1a/){
		    $res_id_1 = "$tmp_1a"."$tmp_2a" ;
		    $res_id_2 = "$tmp_1b"."$tmp_2b" ;
		} elsif ($chains_1 =~ /$tmp_1b/){
		    $res_id_1 = "$tmp_1b"."$tmp_2b" ;
		    $res_id_2 = "$tmp_1a"."$tmp_2a" ;
		}
		$class_1 = $solvexp_1[$resseq_1{$res_id_1}] ; 
		$class_2 = $solvexp_2[$resseq_2{$res_id_2}] ;
		$v_hbond_ij=0;
		if (($tmp_7 >0) && ($tmp_8 >0) && ($tmp_9 >0)){ 
		    $v_hbond_ij=&v_hbond($tmp_8,$tmp_7,$tmp_9,$tmp_6,$class_1,$class_2) ;
		}
		$hb_energy{$donor}{$acceptor} = $v_hbond_ij ;
		$hb_energy{$acceptor}{$donor} = $v_hbond_ij ;
	    } else {
		unless (($tmp_3a eq 'HOH') || ($tmp_3b eq 'HOH')){
#		    print "$tmp_3a $tmp_4a\n";
#		    print "$tmp_3b $tmp_4b\n";
#		    die ;
		}
	    }
	}
    }
}

#################################################################################

sub v_hbond {
    my ($d_ha,$theta,$psi,$hb_type,$class_1,$class_2)=@_ ;
    my ($v_d_ha,$v_theta,$v_psi,$x_tmp,$cat_tmp,$v_hbond_ij);
    my $r_0 = 1.9  ; 
    my $w_d = 1.0  ; 
    my $w_t = 1.03 ; 
    my $w_p = 0.2  ;
    my @env_class = ([1,1,2],
		     [1,2,3],
		     [2,3,3]);
    my @w_hb = (0.49, 0.35, 0.95, 2.04) ;    
    my $r_smooth = 0.25 ;
    my $v_hbond_ij_th = 1e6;

    $v_hbond_ij = 0;
    if ($d_ha <  ($r_0 - $r_smooth)) {
	$x_tmp = ($d_ha + $r_smooth)/$r_0 ; 
    } elsif ($d_ha >  ($r_0 + $r_smooth)){
	$x_tmp = ($d_ha - $r_smooth)/$r_0 ;	
    } else {
	$x_tmp = 1 ;
    }

    $v_d_ha = (5/($x_tmp**12) - 6/($x_tmp**10)) ;

    if ($theta < 100){ 
	$v_theta = 0 ;
    } elsif ($theta < 110){
	$v_theta = -0.35 ;
    } elsif ($theta < 120){
	$v_theta = -0.48 ;
    } elsif ($theta < 130){
	$v_theta = -0.62 ;
    } elsif ($theta < 140){
	$v_theta = -0.74 ;
    } elsif ($theta < 150){
	$v_theta = -0.83 ;
    } elsif ($theta < 160){
	$v_theta = -0.91 ;
    } elsif ($theta < 170){
	$v_theta = -0.97 ;
    } elsif ($theta < 180){
	$v_theta = -1.00 ;
    }

    if ($psi < 90){ 
	$v_psi = -0.325;
    } elsif ($psi < 100){
	$v_psi = -0.565 ;
    } elsif ($psi < 110){
	$v_psi = -0.820 ;
    } elsif ($psi < 120){
	$v_psi = -0.950 ;
    } elsif ($psi < 130){
	$v_psi = -1.000 ;
    } elsif ($psi < 140){
	$v_psi = -0.935 ;
    } elsif ($psi < 150){
	$v_psi = -0.865 ;
    } elsif ($psi < 160){
	$v_psi = -0.680 ;
    } elsif ($psi < 170){
	$v_psi = -0.550 ;
    } elsif ($psi < 180){
	$v_psi = -0.355 ;
    }

    if ( $hb_type =~ /SS/i){
	$cat_tmp = $env_class[$class_1][$class_2];
    } else {
	$cat_tmp = 0;
    }

    $v_hbond_ij = $w_hb[$cat_tmp]*($w_d*$v_d_ha + $w_t*$v_theta + $w_p*$v_psi) ;    
    if ($v_hbond_ij > $v_hbond_ij_th){$v_hbond_ij = $v_hbond_ij_th}
    return $v_hbond_ij;
}

###########################################################################

sub interface_residues{
    my (@int_1_flag,@int_2_flag);

    @int_1_flag=@int_2_flag=();
    foreach $x (0..$#resnum_1) {
	foreach $y (0..$#resnum_2) {
	    if (&aa_dist_cb($x,$y) <30){
		if (&aa_ctc($x,$y)) {
		    $int_1_flag[$x]+=1;
		    $int_2_flag[$y]+=1;
		}
	    }
	}
    }

    @int_1=@int_2=();
    foreach $x (0 .. $#int_1_flag){
	if ($int_1_flag[$x]){
	    push(@int_1,$x);
	}
    }
    foreach $y (0 .. $#int_2_flag){
	if ($int_2_flag[$y]){
	    push(@int_2,$y);
	}
    }
}

###########################################################################

sub models{

    $model[0] = $data_dir.'/svm_model';                           # CDFHIKLMNQSTVWY 
    $model[1] = $data_dir.'/svm_model_E';                         # E
    $model[2] = $data_dir.'/svm_model_R';                         # R

    @mean=@stdev=();                                        # read mean and stdev for input
    foreach $i (0 .. $#model){                              # features normalization
	open PAR_FILE, "${model[$i]}.par" or 
	    die "minchia ..." ; 
	$_ = <PAR_FILE>;
	chomp;
	@fields = split;
	push(@mean,[@fields]);
	$_ = <PAR_FILE>;
	chomp;
	@fields = split;
	push(@stdev,[@fields]);
	close PAR_FILE;
    }
}

###########################################################################

sub aa_dist_cb {
    my ($x,$y)=@_ ;
    my $dist_cb;
    $dist_cb = 0;
    foreach $k (0..2){
        $dist_cb += ($xyz_1_cb[$x][$k] - $xyz_2_cb[$y][$k])**2 ;
    }
    $dist_cb=sqrt($dist_cb);
    $dist_cb;
}

###########################################################################

sub aa_ctc {
    my ($x,$y)=@_ ;
    my ($dist,$contact,$i,$j,$k);
    my $dist_th=5;
    foreach $i (0..$#{$atom_1[$x]}){
	foreach $j (0..$#{$atom_2[$y]}){
	    $dist = 0;
	    foreach $k (0..2){
		$dist += ($xyz_1[$atom_1[$x][$i]][$k] - $xyz_2[$atom_2[$y][$j]][$k])**2 ;
	    }
	    $dist=sqrt($dist);
	    if ($dist < $dist_th) {
		return 1;
	    }
	}
    }
    return 0 ;
}

###########################################################################

sub vdw_sc {
    my ($monomer,$z_m)=@_;
    my ($ddg_vdw_sc,$v_vdw_ij,$dist,$x,$y,$i,$j,$k,$ip,$iq,$sigma_ij,$epsilon_ij,$x_tmp);
    
    my $dist_th = 8;
    my $r_smooth = 0.5 ;
    my $v_vdw_th = 1e6;

    $ddg_vdw_sc=0;
    if ($monomer == 1){
	foreach $y (0..$#resnum_2) {
	    if (&aa_dist_cb($z_m,$y) <30){
		foreach $i (@{$sidechain_1[$z_m]}){		    
		    foreach $j (@{$atom_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_2[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $dist_th) {
			    $ip = $atom_type_1[$i] ;
			    $iq = $atom_type_2[$j] ;
			    $sigma_ij = $r_min[$ip] + $r_min[$iq];
			    $epsilon_ij = sqrt($epsilon[$ip]*$epsilon[$iq]);
			    if ($dist <  ($sigma_ij - $r_smooth)) {
				$x_tmp = ($dist+$r_smooth)/$sigma_ij ; 
			    } elsif ($dist >  ($sigma_ij + $r_smooth)){
				$x_tmp = ($dist-$r_smooth)/$sigma_ij ;	
			    } else {
				$x_tmp = 1 ;
			    }
			    $v_vdw_ij = $epsilon_ij * (1/($x_tmp**12) - 2/($x_tmp**6)) ;
			    if ($v_vdw_ij > $v_vdw_th){$v_vdw_ij=$v_vdw_th;}
			    $ddg_vdw_sc += $v_vdw_ij ;
			}
		    }
		}
	    }
	}
    } elsif ($monomer == 2){ 
	foreach $x (0..$#resnum_1) {
	    if (&aa_dist_cb($x,$z_m) <30){
		foreach $j (@{$sidechain_2[$z_m]}){		    
		    foreach $i (@{$atom_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_2[$j][$k])**2 ;
			}
 			$dist=sqrt($dist);
			if ($dist < $dist_th) {
			    $ip = $atom_type_1[$i] ;
			    $iq = $atom_type_2[$j] ;
			    $sigma_ij = $r_min[$ip] + $r_min[$iq];
			    $epsilon_ij = sqrt($epsilon[$ip]*$epsilon[$iq]);
			    if ($dist <  ($sigma_ij - $r_smooth)) {
				$x_tmp = ($dist+$r_smooth)/$sigma_ij ; 
			    } elsif ($dist >  ($sigma_ij + $r_smooth)){
				$x_tmp = ($dist-$r_smooth)/$sigma_ij ;	
			    } else {
				$x_tmp = 1 ;
			    }
			    $v_vdw_ij = $epsilon_ij * (1/($x_tmp**12) - 2/($x_tmp**6)) ;
			    if ($v_vdw_ij > $v_vdw_th){$v_vdw_ij=$v_vdw_th;}
			    $ddg_vdw_sc += $v_vdw_ij ;
			}
		    }
		}
	    }
	}
    }
    return $ddg_vdw_sc;
}

###########################################################################

sub hb_sc{
    my ($monomer,$z_m)=@_;
    my ($ddg_hb_sc,$atname_tmp_1,$atname_tmp_2);

    $ddg_hb_sc=0;
    if ($monomer == 1){
	foreach $y (0..$#resnum_2) {
	    if (&aa_dist_cb($z_m,$y) <30){
		foreach $i (@{$sidechain_1[$z_m]}){		    
		    foreach $j (@{$atom_2[$y]}){
			$atname_tmp_1 =
			    "$resnum_1[$z_m]"."$atom_pdbname_1[$i]" ;
			$atname_tmp_2 =
			    "$resnum_2[$y]"."$atom_pdbname_2[$j]" ;
			if (exists ($hb_energy{$atname_tmp_1}{$atname_tmp_2})){
			    $ddg_hb_sc += $hb_energy{$atname_tmp_1}{$atname_tmp_2} ;
			}
		    }
		}
	    }
	}			    
    } elsif ($monomer == 2){ 
	foreach $x (0..$#resnum_1) {
	    if (&aa_dist_cb($x,$z_m) <30){
		foreach $j (@{$sidechain_2[$z_m]}){		    
		    foreach $i (@{$atom_1[$x]}){
			$atname_tmp_1 =
			    "$resnum_1[$x]"."$atom_pdbname_1[$i]" ;
			$atname_tmp_2 =
			    "$resnum_2[$z_m]"."$atom_pdbname_2[$j]" ;
			if (exists ($hb_energy{$atname_tmp_1}{$atname_tmp_2})){
			    $ddg_hb_sc += $hb_energy{$atname_tmp_1}{$atname_tmp_2};
			}
		    }
		}
	    }
	}
    }
    return $ddg_hb_sc;
}

###########################################################################

sub solv_sc {
    my ($monomer,$z_m)=@_;
    my ($ddg_solv_sc,$v_solv_ij,$x,$y,$i,$j,$k,$dist,$ip,$iq);
    my ($f_i,$alpha_i,$x_i,$f_j,$alpha_j,$x_j);

    my $dist_th = 8;
    
    $ddg_solv_sc=0;
    if ($monomer == 1){
	foreach $y (0..$#resnum_2) {
	    if (&aa_dist_cb($z_m,$y) <30){
		foreach $i (@{$sidechain_1[$z_m]}){		    
		    foreach $j (@{$atom_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_2[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $dist_th) {
			    $ip = $atom_type_1[$i] ;
			    $iq = $atom_type_2[$j] ;

			    $x_i = ($dist - $vdwrad[$ip])/$corrlen[$ip] ;
			    $alpha_i = $dgfree[$ip]/(2*PI*sqrt(PI)*$corrlen[$ip]);
			    $f_i =  $alpha_i * exp(-$x_i**2)/($dist **2);

			    $x_j = ($dist - $vdwrad[$iq])/$corrlen[$iq] ;
			    $alpha_j = $dgfree[$iq]/(2*PI*sqrt(PI)*$corrlen[$iq]);
			    $f_j =  $alpha_j * exp(-$x_j**2)/($dist **2);

			    $v_solv_ij = - ($f_i * $fracvol[$iq] + $f_j * $fracvol[$ip]) ;
			    $ddg_solv_sc += $v_solv_ij ;
			}
		    }
		}
	    }
	}
    } elsif ($monomer == 2){ 
	foreach $x (0..$#resnum_1) {
	    if (&aa_dist_cb($x,$z_m) <30){
		foreach $j (@{$sidechain_2[$z_m]}){		    
		    foreach $i (@{$atom_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_2[$j][$k])**2 ;
			}
 			$dist=sqrt($dist);
			if ($dist < $dist_th) {
			    $ip = $atom_type_1[$i] ;
			    $iq = $atom_type_2[$j] ;

			    $x_i = ($dist - $vdwrad[$ip])/$corrlen[$ip] ;
			    $alpha_i = $dgfree[$ip]/(2*PI*sqrt(PI)*$corrlen[$ip]);
			    $f_i =  $alpha_i * exp(-$x_i**2)/($dist **2);

			    $x_j = ($dist - $vdwrad[$iq])/$corrlen[$iq] ;
			    $alpha_j = $dgfree[$iq]/(2*PI*sqrt(PI)*$corrlen[$iq]);
			    $f_j =  $alpha_j * exp(-$x_j**2)/($dist **2);

			    $v_solv_ij = - ($f_i * $fracvol[$iq] + $f_j * $fracvol[$ip]) ;
			    $ddg_solv_sc += $v_solv_ij ;
			}
		    }
		}
	    }
	}
    }
    return $ddg_solv_sc;
}

############################################################################################

sub vdw_frag {
    my ($monomer,$z_m)=@_;
    my ($ddg_vdw_frag,$v_vdw_ij,$dist,$x,$y,$i,$j,$k,$ip,$iq,$sigma_ij,$epsilon_ij,$x_tmp);
    my (@frag_1,@frag_2,$dist_cb);
    
    my $dist_th = 8;
    my $d_frag = 10;
    my $r_smooth = 0.5 ;
    my $v_vdw_th = 1e6;

    $ddg_vdw_frag=0;
    @frag_1=@frag_2=();
    if ($monomer == 1){
	foreach $x (0..$#resnum_1) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$z_m][$k] - $xyz_1_cb[$x][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		if ($x != $z_m){
		    foreach $i (@{$atom_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_1,$i);
			}
		    }
		} else {
		    foreach $i (@{$mainchain_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_1,$i);
			}
		    }
		}
	    }
	}
	foreach $y (0..$#resnum_2) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$z_m][$k] - $xyz_2_cb[$y][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		foreach $j (@{$atom_2[$y]}){
		    $dist = 0;
		    foreach $k (0..2){
			$dist += ($xyz_2[$j][$k] - $xyz_1_cb[$z_m][$k])**2 ;
		    }
		    $dist=sqrt($dist);
		    if ($dist < $d_frag) {
			push(@frag_2,$j);
		    }
		}
	    }
	}
    } elsif ($monomer == 2){ 
	foreach $x (0..$#resnum_1) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$x][$k] - $xyz_2_cb[$z_m][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		foreach $i (@{$atom_1[$x]}){
		    $dist = 0;
		    foreach $k (0..2){
			$dist += ($xyz_1[$i][$k] - $xyz_2_cb[$z_m][$k])**2 ;
		    }
		    $dist=sqrt($dist);
		    if ($dist < $d_frag) {
			push(@frag_1,$i);
		    }
		}
	    }
	}
	foreach $y (0..$#resnum_2) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_2_cb[$y][$k] - $xyz_2_cb[$z_m][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		if ($y != $z_m){
		    foreach $j (@{$atom_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$j][$k] - $xyz_2_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_2,$j);
			}
		    }
		} else {
		    foreach $j (@{$mainchain_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$j][$k] - $xyz_2_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_2,$j);
			}
		    }
		}
	    }
	}
    }
    foreach $i (@frag_1){
	foreach $j (@frag_2){
	    $dist = 0;
	    foreach $k (0..2){
		$dist += ($xyz_1[$i][$k] - $xyz_2[$j][$k])**2 ;
	    }
	    $dist=sqrt($dist);
	    if ($dist < $dist_th) {
		$ip = $atom_type_1[$i] ;
		$iq = $atom_type_2[$j] ;
		$sigma_ij = $r_min[$ip] + $r_min[$iq];
		$epsilon_ij = sqrt($epsilon[$ip]*$epsilon[$iq]);
		if ($dist <  ($sigma_ij - $r_smooth)) {
		    $x_tmp = ($dist+$r_smooth)/$sigma_ij ; 
		} elsif ($dist >  ($sigma_ij + $r_smooth)){
		    $x_tmp = ($dist-$r_smooth)/$sigma_ij ;	
		} else {
		    $x_tmp = 1 ;
		}
		$v_vdw_ij = $epsilon_ij * (1/($x_tmp**12) - 2/($x_tmp**6)) ;
		if ($v_vdw_ij > $v_vdw_th){$v_vdw_ij=$v_vdw_th;}
		$ddg_vdw_frag += $v_vdw_ij ;
	    }
	}
    }
    return $ddg_vdw_frag;
}

###########################################################################

sub hb_frag{
    my ($monomer,$z_m)=@_;
    my ($ddg_hb_frag,$atname_tmp_1,$atname_tmp_2);
    my (@frag_1,@frag_2,$dist,$dist_cb);

    my $d_frag = 10;

    $ddg_hb_frag=0;
    @frag_1=@frag_2=();
    if ($monomer == 1){
	foreach $x (0..$#resnum_1) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$z_m][$k] - $xyz_1_cb[$x][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		if ($x != $z_m){
		    foreach $i (@{$atom_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_1,[($x,$i)]);
			}
		    }
		} else {
		    foreach $i (@{$mainchain_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_1,[($x,$i)]);
			}
		    }
		}
	    }
	}
	foreach $y (0..$#resnum_2) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$z_m][$k] - $xyz_2_cb[$y][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		foreach $j (@{$atom_2[$y]}){
		    $dist = 0;
		    foreach $k (0..2){
			$dist += ($xyz_2[$j][$k] - $xyz_1_cb[$z_m][$k])**2 ;
		    }
		    $dist=sqrt($dist);
		    if ($dist < $d_frag) {
			push(@frag_2,[($y,$j)]);
		    }
		}
	    }
	}
    } elsif ($monomer == 2){ 
	foreach $x (0..$#resnum_1) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$x][$k] - $xyz_2_cb[$z_m][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		foreach $i (@{$atom_1[$x]}){
		    $dist = 0;
		    foreach $k (0..2){
			$dist += ($xyz_1[$i][$k] - $xyz_2_cb[$z_m][$k])**2 ;
		    }
		    $dist=sqrt($dist);
		    if ($dist < $d_frag) {
			push(@frag_1,[($x,$i)]);
		    }
		}
	    }
	}
	foreach $y (0..$#resnum_2) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_2_cb[$y][$k] - $xyz_2_cb[$z_m][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		if ($y != $z_m){
		    foreach $j (@{$atom_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$j][$k] - $xyz_2_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_2,[($y,$j)]);
			}
		    }
		} else {
		    foreach $j (@{$mainchain_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$j][$k] - $xyz_2_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_2,[($y,$j)]);
			}
		    }
		}
	    }
	}
    }
    foreach $i (0..$#frag_1){
	foreach $j (0..$#frag_2){
	    $atname_tmp_1 =
		"$resnum_1[$frag_1[$i][0]]"."$atom_pdbname_1[$frag_1[$i][1]]" ;
	    $atname_tmp_2 =
		"$resnum_2[$frag_2[$j][0]]"."$atom_pdbname_2[$frag_2[$j][1]]" ;
	    if (exists ($hb_energy{$atname_tmp_1}{$atname_tmp_2})){
		$ddg_hb_frag += $hb_energy{$atname_tmp_1}{$atname_tmp_2} ;
	    }
	}
    }
    return $ddg_hb_frag;
}

###############################################################################################

sub solv_frag {
    my ($monomer,$z_m)=@_;
    my ($ddg_solv_frag,$v_solv_ij,$dist,$x,$y,$i,$j,$k,$ip,$iq);
    my ($f_i,$alpha_i,$x_i,$f_j,$alpha_j,$x_j);
    my (@frag_1,@frag_2,$dist_cb);
    
    my $dist_th = 8;
    my $d_frag = 10;

    $ddg_solv_frag=0;
    @frag_1=@frag_2=();
    if ($monomer == 1){
	foreach $x (0..$#resnum_1) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$z_m][$k] - $xyz_1_cb[$x][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		if ($x != $z_m){
		    foreach $i (@{$atom_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_1,$i);
			}
		    }
		} else {
		    foreach $i (@{$mainchain_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_1,$i);
			}
		    }
		}
	    }
	}
	foreach $y (0..$#resnum_2) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$z_m][$k] - $xyz_2_cb[$y][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		foreach $j (@{$atom_2[$y]}){
		    $dist = 0;
		    foreach $k (0..2){
			$dist += ($xyz_2[$j][$k] - $xyz_1_cb[$z_m][$k])**2 ;
		    }
		    $dist=sqrt($dist);
		    if ($dist < $d_frag) {
			push(@frag_2,$j);
		    }
		}
	    }
	}
    } elsif ($monomer == 2){ 
	foreach $x (0..$#resnum_1) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$x][$k] - $xyz_2_cb[$z_m][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		foreach $i (@{$atom_1[$x]}){
		    $dist = 0;
		    foreach $k (0..2){
			$dist += ($xyz_1[$i][$k] - $xyz_2_cb[$z_m][$k])**2 ;
		    }
		    $dist=sqrt($dist);
		    if ($dist < $d_frag) {
			push(@frag_1,$i);
		    }
		}
	    }
	}
	foreach $y (0..$#resnum_2) {
	    $dist_cb = 0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_2_cb[$y][$k] - $xyz_2_cb[$z_m][$k])**2 ;
	    }
	    $dist_cb=sqrt($dist_cb);
	    if ($dist_cb <30){
		if ($y != $z_m){
		    foreach $j (@{$atom_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$j][$k] - $xyz_2_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_2,$j);
			}
		    }
		} else {
		    foreach $j (@{$mainchain_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$j][$k] - $xyz_2_cb[$z_m][$k])**2 ;
			}
			$dist=sqrt($dist);
			if ($dist < $d_frag) {
			    push(@frag_2,$j);
			}
		    }
		}
	    }
	}
    }
    foreach $i (@frag_1){
	foreach $j (@frag_2){
	    $dist = 0;
	    foreach $k (0..2){
		$dist += ($xyz_1[$i][$k] - $xyz_2[$j][$k])**2 ;
	    }
	    $dist=sqrt($dist);
	    if ($dist < $dist_th) {
		$ip = $atom_type_1[$i] ;
		$iq = $atom_type_2[$j] ;

		$x_i = ($dist - $vdwrad[$ip])/$corrlen[$ip] ;
		$alpha_i = $dgfree[$ip]/(2*PI*sqrt(PI)*$corrlen[$ip]);
		$f_i =  $alpha_i * exp(-$x_i**2)/($dist **2);
		
		$x_j = ($dist - $vdwrad[$iq])/$corrlen[$iq] ;
		$alpha_j = $dgfree[$iq]/(2*PI*sqrt(PI)*$corrlen[$iq]);
		$f_j =  $alpha_j * exp(-$x_j**2)/($dist **2);
		
		$v_solv_ij = - ($f_i * $fracvol[$iq] + $f_j * $fracvol[$ip]) ;
		$ddg_solv_frag += $v_solv_ij ;
	    }
	}
    }
    return $ddg_solv_frag;
}

###############################################################################################

sub vdw_intra_sc {
    my ($monomer,$z_m)=@_;
    my ($ddg_vdw_intra_sc,$v_vdw_ij,$dist,,$dist_cb,$x,$y,$i,$j,$k,$ip,$iq,$sigma_ij,$epsilon_ij,$x_tmp);
    
    my $dist_th = 8;
    my $r_smooth = 0.5 ;
    my $v_vdw_th = 1e6;
    my $nn_dist=3;

    $ddg_vdw_intra_sc=0;
    if ($monomer == 1){
	foreach $x (0..$#resnum_1) {
	    if ($x == $z_m){
		foreach $i (@{$sidechain_1[$z_m]}){
		    foreach $j (@{$mainchain_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if (($dist < $dist_th) && ($conn_mat_1[$i][$j] > $nn_dist)){
			    $ip = $atom_type_1[$i] ;
			    $iq = $atom_type_1[$j] ;
			    $sigma_ij = $r_min[$ip] + $r_min[$iq];
			    $epsilon_ij = sqrt($epsilon[$ip]*$epsilon[$iq]);
			    if ($dist <  ($sigma_ij - $r_smooth)) {
				$x_tmp = ($dist+$r_smooth)/$sigma_ij ; 
			    } elsif ($dist >  ($sigma_ij + $r_smooth)){
				$x_tmp = ($dist-$r_smooth)/$sigma_ij ;	
			    } else {
				$x_tmp = 1 ;
			    }
			    $v_vdw_ij = $epsilon_ij * (1/($x_tmp**12) - 2/($x_tmp**6)) ;
			    if ($v_vdw_ij > $v_vdw_th){$v_vdw_ij=$v_vdw_th;}
			    $ddg_vdw_intra_sc += $v_vdw_ij ;
			}
		    }
		    foreach $j (@{$sidechain_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if (($dist < $dist_th) && ($conn_mat_1[$i][$j] > $nn_dist)){
			    $ip = $atom_type_1[$i] ;
			    $iq = $atom_type_1[$j] ;
			    $sigma_ij = $r_min[$ip] + $r_min[$iq];
			    $epsilon_ij = sqrt($epsilon[$ip]*$epsilon[$iq]);
			    if ($dist <  ($sigma_ij - $r_smooth)) {
				$x_tmp = ($dist+$r_smooth)/$sigma_ij ; 
			    } elsif ($dist >  ($sigma_ij + $r_smooth)){
				$x_tmp = ($dist-$r_smooth)/$sigma_ij ;	
			    } else {
				$x_tmp = 1 ;
			    }
			    $v_vdw_ij = $epsilon_ij * (1/($x_tmp**12) - 2/($x_tmp**6)) ;
			    if ($v_vdw_ij > $v_vdw_th){$v_vdw_ij=$v_vdw_th;}
			    $ddg_vdw_intra_sc += ($v_vdw_ij/2) ;
			}
		    }
		}
	    } else {
		$dist_cb = 0;
		foreach $k (0..2){
		    $dist_cb += ($xyz_1_cb[$z_m][$k] - $xyz_1_cb[$x][$k])**2 ;
		}
		$dist_cb=sqrt($dist_cb);
		if ($dist_cb <30){
		    foreach $i (@{$sidechain_1[$z_m]}){		    
			foreach $j (@{$atom_1[$x]}){
			    $dist = 0;
			    foreach $k (0..2){
				$dist += ($xyz_1[$i][$k] - $xyz_1[$j][$k])**2 ;
			    }
			    $dist=sqrt($dist);
			    if (($dist < $dist_th) && ($conn_mat_1[$i][$j] > $nn_dist)){
				$ip = $atom_type_1[$i] ;
				$iq = $atom_type_1[$j] ;
				$sigma_ij = $r_min[$ip] + $r_min[$iq];
				$epsilon_ij = sqrt($epsilon[$ip]*$epsilon[$iq]);
				if ($dist <  ($sigma_ij - $r_smooth)) {
				    $x_tmp = ($dist+$r_smooth)/$sigma_ij ; 
				} elsif ($dist >  ($sigma_ij + $r_smooth)){
				    $x_tmp = ($dist-$r_smooth)/$sigma_ij ;	
				} else {
				    $x_tmp = 1 ;
				}
				$v_vdw_ij = $epsilon_ij * (1/($x_tmp**12) - 2/($x_tmp**6)) ;
				if ($v_vdw_ij > $v_vdw_th){$v_vdw_ij=$v_vdw_th;}
				$ddg_vdw_intra_sc += $v_vdw_ij ;
			    }
			}
		    }
		}
	    }
	}
    } elsif ($monomer == 2){ 
	foreach $y (0..$#resnum_2) {
	    if ($y == $z_m){
		foreach $i (@{$sidechain_2[$z_m]}){		    
		    foreach $j (@{$mainchain_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$i][$k] - $xyz_2[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if (($dist < $dist_th) && ($conn_mat_2[$i][$j] > $nn_dist)){
			    $ip = $atom_type_2[$i] ;
			    $iq = $atom_type_2[$j] ;
			    $sigma_ij = $r_min[$ip] + $r_min[$iq];
			    $epsilon_ij = sqrt($epsilon[$ip]*$epsilon[$iq]);
			    if ($dist <  ($sigma_ij - $r_smooth)) {
				$x_tmp = ($dist+$r_smooth)/$sigma_ij ; 
			    } elsif ($dist >  ($sigma_ij + $r_smooth)){
				$x_tmp = ($dist-$r_smooth)/$sigma_ij ;	
			    } else {
				$x_tmp = 1 ;
			    }
			    $v_vdw_ij = $epsilon_ij * (1/($x_tmp**12) - 2/($x_tmp**6)) ;
			    if ($v_vdw_ij > $v_vdw_th){$v_vdw_ij=$v_vdw_th;}
			    $ddg_vdw_intra_sc += $v_vdw_ij ;
			}
		    }
		    foreach $j (@{$sidechain_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$i][$k] - $xyz_2[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if (($dist < $dist_th) && ($conn_mat_2[$i][$j] > $nn_dist)){
			    $ip = $atom_type_2[$i] ;
			    $iq = $atom_type_2[$j] ;
			    $sigma_ij = $r_min[$ip] + $r_min[$iq];
			    $epsilon_ij = sqrt($epsilon[$ip]*$epsilon[$iq]);
			    if ($dist <  ($sigma_ij - $r_smooth)) {
				$x_tmp = ($dist+$r_smooth)/$sigma_ij ; 
			    } elsif ($dist >  ($sigma_ij + $r_smooth)){
				$x_tmp = ($dist-$r_smooth)/$sigma_ij ;	
			    } else {
				$x_tmp = 1 ;
			    }
			    $v_vdw_ij = $epsilon_ij * (1/($x_tmp**12) - 2/($x_tmp**6)) ;
			    if ($v_vdw_ij > $v_vdw_th){$v_vdw_ij=$v_vdw_th;}
			    $ddg_vdw_intra_sc += ($v_vdw_ij/2) ;
			}
		    }
		}
	    } else {
		$dist_cb = 0;
		foreach $k (0..2){
		    $dist_cb += ($xyz_2_cb[$z_m][$k] - $xyz_2_cb[$y][$k])**2 ;
		}
		$dist_cb=sqrt($dist_cb);
		if ($dist_cb <30){
		    foreach $i (@{$sidechain_2[$z_m]}){		    
			foreach $j (@{$atom_2[$y]}){
			    $dist = 0;
			    foreach $k (0..2){
				$dist += ($xyz_2[$i][$k] - $xyz_2[$j][$k])**2 ;
			    }
			    $dist=sqrt($dist);
			    if (($dist < $dist_th) && ($conn_mat_2[$i][$j] > $nn_dist)){
				$ip = $atom_type_2[$i] ;
				$iq = $atom_type_2[$j] ;
				$sigma_ij = $r_min[$ip] + $r_min[$iq];
				$epsilon_ij = sqrt($epsilon[$ip]*$epsilon[$iq]);
				if ($dist <  ($sigma_ij - $r_smooth)) {
				    $x_tmp = ($dist+$r_smooth)/$sigma_ij ; 
				} elsif ($dist >  ($sigma_ij + $r_smooth)){
				    $x_tmp = ($dist-$r_smooth)/$sigma_ij ;	
				} else {
				    $x_tmp = 1 ;
				}
				$v_vdw_ij = $epsilon_ij * (1/($x_tmp**12) - 2/($x_tmp**6)) ;
				if ($v_vdw_ij > $v_vdw_th){$v_vdw_ij=$v_vdw_th;}
				$ddg_vdw_intra_sc += $v_vdw_ij ;
			    }
			}
		    }
		}
	    }
	}
    }
    return $ddg_vdw_intra_sc;
}

#################################################################################################

sub coul_intra_sc {
    my ($monomer,$z_m)=@_;
    my ($ddg_coul_intra_sc,$v_coul_ij,$dist,,$dist_cb,$x,$y,$i,$j,$k,$xp,$xq);
    
    my $dist_th = 8;
    my $nn_dist=3;

    $ddg_coul_intra_sc=0;
    if ($monomer == 1){
	foreach $x (0..$#resnum_1) {
	    if ($x == $z_m){
		foreach $i (@{$sidechain_1[$z_m]}){
		    foreach $j (@{$mainchain_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if (($dist < $dist_th) && ($conn_mat_1[$i][$j] > $nn_dist)){
			    $xp = $atom_charge_1[$i] ;
			    $xq = $atom_charge_1[$j] ;
			    $v_coul_ij = 332.0637*($xp*$xq)/($dist*$dist);
			    $ddg_coul_intra_sc += $v_coul_ij ;
			}
		    }
		    foreach $j (@{$sidechain_1[$x]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_1[$i][$k] - $xyz_1[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if (($dist < $dist_th) && ($conn_mat_1[$i][$j] > $nn_dist)){
			    $xp = $atom_charge_1[$i] ;
			    $xq = $atom_charge_1[$j] ;
			    $v_coul_ij = 332.0637*($xp*$xq)/($dist*$dist);
			    $ddg_coul_intra_sc += ($v_coul_ij/2) ;
			}
		    }
		}
	    } else {
		$dist_cb = 0;
		foreach $k (0..2){
		    $dist_cb += ($xyz_1_cb[$z_m][$k] - $xyz_1_cb[$x][$k])**2 ;
		}
		$dist_cb=sqrt($dist_cb);
		if ($dist_cb <30){
		    foreach $i (@{$sidechain_1[$z_m]}){		    
			foreach $j (@{$atom_1[$x]}){
			    $dist = 0;
			    foreach $k (0..2){
				$dist += ($xyz_1[$i][$k] - $xyz_1[$j][$k])**2 ;
			    }
			    $dist=sqrt($dist);
			    if (($dist < $dist_th) && ($conn_mat_1[$i][$j] > $nn_dist)){
				$xp = $atom_charge_1[$i] ;
				$xq = $atom_charge_1[$j] ;
				$v_coul_ij = 332.0637*($xp*$xq)/($dist*$dist);
				$ddg_coul_intra_sc += $v_coul_ij;
			    }
			}
		    }
		}
	    }
	}
    } elsif ($monomer == 2){ 
	foreach $y (0..$#resnum_2) {
	    if ($y == $z_m){
		foreach $i (@{$sidechain_2[$z_m]}){		    
		    foreach $j (@{$mainchain_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$i][$k] - $xyz_2[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if (($dist < $dist_th) && ($conn_mat_2[$i][$j] > $nn_dist)){
			    $xp = $atom_charge_2[$i] ;
			    $xq = $atom_charge_2[$j] ;
			    $v_coul_ij = 332.0637*($xp*$xq)/($dist*$dist);
			    $ddg_coul_intra_sc += $v_coul_ij;
			}
		    }
		    foreach $j (@{$sidechain_2[$y]}){
			$dist = 0;
			foreach $k (0..2){
			    $dist += ($xyz_2[$i][$k] - $xyz_2[$j][$k])**2 ;
			}
			$dist=sqrt($dist);
			if (($dist < $dist_th) && ($conn_mat_2[$i][$j] > $nn_dist)){
			    $xp = $atom_charge_2[$i] ;
			    $xq = $atom_charge_2[$j] ;
			    $v_coul_ij = 332.0637*($xp*$xq)/($dist*$dist);
			    $ddg_coul_intra_sc += ($v_coul_ij/2);
			}
		    }
		}
	    } else {
		$dist_cb = 0;
		foreach $k (0..2){
		    $dist_cb += ($xyz_2_cb[$z_m][$k] - $xyz_2_cb[$y][$k])**2 ;
		}
		$dist_cb=sqrt($dist_cb);
		if ($dist_cb <30){
		    foreach $i (@{$sidechain_2[$z_m]}){		    
			foreach $j (@{$atom_2[$y]}){
			    $dist = 0;
			    foreach $k (0..2){
				$dist += ($xyz_2[$i][$k] - $xyz_2[$j][$k])**2 ;
			    }
			    $dist=sqrt($dist);
			    if (($dist < $dist_th) && ($conn_mat_2[$i][$j] > $nn_dist)){
				$xp = $atom_charge_2[$i] ;
				$xq = $atom_charge_2[$j] ;
				$v_coul_ij = 332.0637*($xp*$xq)/($dist*$dist);
				$ddg_coul_intra_sc += $v_coul_ij;
			    }
			}
		    }
		}
	    }
	}
    }
    return $ddg_coul_intra_sc;
}

#################################################################################################



sub solvexp {
    my ($solvexp_1,$solvexp_2) = @_;
    my (@ncb_1,@ncb_2,$x,$y,$dist_cb,$k);
    my $dist_cb_th = 8;

    @ncb_1=@ncb_2=();
    foreach $x (0..$#resnum_1-1) {
	foreach $y ($x+1..$#resnum_1) {
	    $dist_cb=0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$x][$k] - $xyz_1_cb[$y][$k])**2 ;		
	    }
	    $dist_cb = sqrt($dist_cb);
	    if ($dist_cb < $dist_cb_th){
		$ncb_1[$x] += 1;
		$ncb_1[$y] += 1;
	    }
	}
    }
    foreach $x (0..$#resnum_2-1) {
	foreach $y ($x+1..$#resnum_2) {
	    $dist_cb=0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_2_cb[$x][$k] - $xyz_2_cb[$y][$k])**2 ;		
	    }
	    $dist_cb = sqrt($dist_cb);
	    if ($dist_cb < $dist_cb_th){
		$ncb_2[$x] += 1;
		$ncb_2[$y] += 1;
	    }
	}
    }
    foreach $x (0..$#resnum_1) {
	foreach $y (0..$#resnum_2) {
	    $dist_cb=0;
	    foreach $k (0..2){
		$dist_cb += ($xyz_1_cb[$x][$k] - $xyz_2_cb[$y][$k])**2 ;		
	    }
	    $dist_cb = sqrt($dist_cb);
	    if ($dist_cb < $dist_cb_th){
		$ncb_1[$x] += 1;
		$ncb_2[$y] += 1;
	    }
	}
    }
    foreach $x (0..$#resnum_1) {
	if ((! $ncb_1[$x]) || ($ncb_1[$x] <= 8)) {
	    $$solvexp_1[$x]=0 ;
	} elsif ($ncb_1[$x] <= 14) {
	    $$solvexp_1[$x]=1 ;
	} else {
	    $$solvexp_1[$x]=2 ;
	}
    }
    foreach $y (0..$#resnum_2) {
	if ((! $ncb_2[$y]) ||($ncb_2[$y] <= 8)) {
	    $$solvexp_2[$y]=0 ;
	} elsif ($ncb_2[$y] <= 14) {
	    $$solvexp_2[$y]=1 ;
	} else {
	    $$solvexp_2[$y]=2 ;
	}
    }
}
