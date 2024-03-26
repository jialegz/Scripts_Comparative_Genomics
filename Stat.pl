
#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($glist,$Sdir,$result,$fkey,$outdir,$p);
GetOptions(
				"help|?" =>\&USAGE,
				"glist:s"=>\$glist,
				"result:s"=>\$result,
				"k:s"=>\$fkey,
				"Sdir:s"=>\$Sdir,
				"p:f"=>\$p,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($glist and $result and $fkey);
$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
$glist=AbsolutePath("file",$glist);
$result=AbsolutePath("file",$result);
############################键入程序###################################################
my %hash;
open(GLIST,$glist)or die $!;
#1       c25855T1c0orf1
#2       c29889T1c0orf1
while(<GLIST>)
{
	chomp;
	next if(/^$/);
	my($id,$geneId)=split(/\s+/,$_);	
	$hash{$id}=$geneId;
}
close(GLIST);
############################################
open(OUT,">$outdir/Paml_branch_site_all.stat")or die $!;
print OUT "#gene\tbranch_pvalue\tsites\n";
open(RES,$result)or die $!;
while(<RES>)
{
	chomp;
	next if(/^$/);
	my @array=split(/\s+/,$_);
	$array[0]=~ s/_mlc20//;
	$array[0]=$hash{$array[0]};
	my $line=join("\t",@array);
	next if($array[1]>=$p||!defined $array[2]);
	print OUT join("\t",@array)."\n";
}
close(OUT);

#&qsub("$outdir/iprscan/work.sh"," --reqsub --queue middle.q  --lines 1 --maxproc $maxproc ");
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zhanghailun<zhanghl\@biomarker.com.cn> 

Usage:
  Options:
  -glist  <file>   Input gene list file, forced
  -Sdir              single dir
   -p     float     p value
  -result	<file>  Input result file,forced
   -k  <str>    Key of output file, forced
  -od <dir>    Dir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
