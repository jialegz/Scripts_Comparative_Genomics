#!/usr/bin/perl 
# Copyright (c) BMK 2018/9/28
# Writer:         liufuyan <liufuyan@biomarker.com.cn>
# Program Date:   2018/9/28.
# Modifier:       liufuyan <liufuyan@biomarker.com.cn>
# Last Modified:  2018/9/28.

use Cwd;
use FindBin qw($Bin $Script);
BEGIN{
     unshift (@INC,"$Bin/Soft/");
}
use  Bio::SeqIO;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);

my $programe_dir=basename($0);
my $path=dirname($0);

my $ver    = "1.0";
my $Writer = "liufuyan <liufuyan\@biomarker.com.cn>";
my $Data   = "2019/1/28";
my $BEGIN=time();

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($pep,$gff,$od);
GetOptions(
			"h|?" =>\&help,
			"od:s"=>\$od,
			"pep:s"=>\$pep,
			"gff:s"=>\$gff,

			) || &help;
&help unless ($pep && $od);

sub help
{
	print <<"	Usage End.";
    Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: 
    Usage:
        -pep        pep file from gffread or JGI longest pep file
        -gff        gff3 file   
        -od         outdir
        -h          Help document
	Usage End.
	exit;
}
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";
################
my %hash;
`mkdir $od`;
my %hash_long;
my $out_longest=basename($pep);
open (OUT, ">$od/$out_longest.longest.pep.fa")|| die "cannot open :$!";
my $ina = Bio::SeqIO->new(-file => $pep, -format => 'fasta');
while(my $obj = $ina->next_seq()){
        my $id = $obj->id;
        my $seq = $obj->seq;
	my $des= $obj->desc;
	if(!defined $des||$des!~/\w+/){
		$des=$id;
	}
        my @new=($id,length$seq,$seq);
        push @{$hash{$des}},[@new];

}
foreach my $key(keys %hash){
  my @array=@{$hash{$key}};
  if($#array>=1){
     @array=sort{$b->[1]<=>$a->[1]} @array;
   }
    $hash_long{$array[0][0]}=1;
   $array[0][2]=~s/\.//g;
  print OUT ">$array[0][0]\n$array[0][2]\n";
   }
close OUT;
my %have;
open(GFF,$gff);
open(OUT2,">$od/$out_longest.gff_for_MCSCANX.txt");
open(OUT3,">$od/$out_longest.gff_for_WGD.txt");
open(OUT4,">$od/$out_longest.bed");
while(<GFF>){
   chomp;
   next if(/^\#/);
   my @tmp=split/\t+/;
   my ($id)=$tmp[8]=~/ID=([^;]+)/;
   if(defined $hash_long{$id}&&!defined $have{"$tmp[0]\t$id\t$tmp[3]\t$tmp[4]"}){
       print OUT2 "$tmp[0]\t$id\t$tmp[3]\t$tmp[4]\n";
       print OUT3 "$tmp[0]\t$tmp[3]\t$tmp[4]\t$id\t0\t$tmp[6]\n";
	print OUT4 "$tmp[0]\t$tmp[3]\t$tmp[4]\t$id\t0\t$tmp[6]\n";
	$have{"$tmp[0]\t$id\t$tmp[3]\t$tmp[4]"}=1;
   }
}
close GFF;




###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);


###############Subs
sub sub_format_datetime #Time calculation subroutine
{
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Runtime # &Runtime($BEGIN);
{
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total $programe_dir elapsed time : [",&sub_time($t),"]\n";
}
sub sub_time
{
	my ($T)=@_;chomp $T;
	my $s=0;my $m=0;my $h=0;
	if ($T>=3600) {
		my $h=int ($T/3600);
		my $a=$T%3600;
		if ($a>=60) {
			my $m=int($a/60);
			$s=$a%60;
			$T=$h."h\-".$m."m\-".$s."s";
		}else{
			$T=$h."h-"."0m\-".$a."s";
		}
	}else{
		if ($T>=60) {
			my $m=int($T/60);
			$s=$T%60;
			$T=$m."m\-".$s."s";
		}else{
			$T=$T."s";
		}
	}
	return ($T);
}

sub AbsolutePath
{		#获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;
	$/="\n";

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'f_dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/";
                chdir($pwd);
        }
		 elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}
