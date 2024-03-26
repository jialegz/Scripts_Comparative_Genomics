
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
my ($zero,$model,$fkey,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"zero:s"=>\$zero,
				"model:s"=>\$model,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($zero and $model and $fkey);
$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
############################键入程序###################################################
$zero=AbsolutePath("file",$zero);
$model=AbsolutePath("file",$model);
my $dir=dirname($zero);
##############################################################################
open(ZERO,$zero)or die $!;
#lnL(ntime: 12  np: 14):   -378.757310      +0.000000
my $lnL1=0;
while(<ZERO>)
{
        chomp;
        next if(/^$/);
        if(/^lnL\(/)
        {
         my(undef,undef,undef,$info)=split(/:/,$_);
         $info=~ s/^\s+//g;
         ($lnL1)=split(/\s+/,$info);
        }
}
close(ZERO);
open(CHOSE,$model)or die $!;
#lnL(ntime: 12  np: 14):   -378.757310      +0.000000
my $lnL2=0;
my $index=0;
my %sites;
while(<CHOSE>)
{
        chomp;
        next if(/^$/);
        if(/^lnL\(/)
        {
         my(undef,undef,undef,$info)=split(/:/,$_);
         $info=~ s/^\s+//g;
         ($lnL2)=split(/\s+/,$info);
        }
		if(/^Bayes/)
		{
			$index++;
		}
		if(/^Positive/)
		{
			$index++;
		}
		if($_=~ /^The\s+grid\s+\(/)
		{
		 $index=1;
		}
		if($index eq 3)
		{
		my $siteinfo=$_;
		$siteinfo=~ s/^\s+//g;
		$siteinfo=~ s/\s+$//g;
		push @{$sites{$model}},$siteinfo;
		}
}
close(CHOSE);
#print $lnL1,"\t",$lnL2,"\n";
my $pvalue=0;
if($lnL1 != 0 && $lnL1 !=0)
{
my $ltr=2*abs($lnL1-$lnL2);
my $info =`/home/jiale/Bioinformatics/Software/PAML/paml-4.10.6/bin/chi2 2 $ltr `;
#df =  2  prob = 0.812280139 = 8.123e-01
(undef,undef,$pvalue)=split(/=/,$info);
$pvalue=~ s/\s+//g;
}
my $name=basename($model);
open(OUT,">$dir/$name.result");
	print OUT "$name\t$pvalue";
foreach my $key (keys %sites)
{	next if(! exists $sites{$key});
	my @array=@{$sites{$key}};
	foreach my $line (@array)
	{
	next if($line =~ /^Positive/);
	my ($pos,$type,$star)=split(/\s+/,$line);
	if($star=~ /\*/)
	{
		print OUT "\t",$pos,",",$type,",",$star;
	}
	}
}
print OUT "\n";
close(OUT);
#&qsub("$outdir/iprscan/work.sh"," --reqsub --queue middle.q  --lines 1 --maxproc $maxproc ");
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub qsub($$){
	my $shell = shift;
	my $qsub_para = shift;
	my $note = `hostname`;
	my $cmd;
	if($note=~ m/cluster/){
		$cmd = "qsub-sge.pl $shell $qsub_para";
		`$cmd`;
	}else{
		$cmd = "ssh cluster -Y qsub-sge.pl $shell $qsub_para";
		`$cmd`;
	}
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
  -zero  <file>   Input zero mlc file, forced
  -model <file>   Input model mlc file,forced  
  -k  <str>    Key of output file, forced
  -od <dir>    Dir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
