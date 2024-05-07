#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.2";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($list,$tree,$mafft,$codeml,$outdir);
my $p=0.05;
GetOptions(
				"help|?" =>\&USAGE,
				"list:s"=>\$list,
				"tree:s"=>\$tree,
				"p:f"=>\$p,
				"mafft:s"=>\$mafft,
				"codeml:s"=>\$codeml,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($list and $tree );
$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
`mkdir $outdir/conf/`;
$outdir=AbsolutePath("dir",$outdir);
`mkdir $outdir/mlc/`;
my $config1;
my $config2;
my @array1;
my @array2;
my @genelist ;
my @index;
print "creat the config file\n";
open(IN,"$list");
<IN>;
while(<IN>){
chomp;
my @tmp=split/\s+/;
push @index,$tmp[0];
}
close IN;

open(CTL , ">$outdir/ctl.list")or die $!;
for(my $i=0 ;$i<=$#index;$i++)
{       my $seq ="$mafft/$index[$i].fa.mafft.fa.codon.fa";
	$seq=AbsolutePath("file",$seq );
        my $tree ="$tree/Grouptree/$index[$i].tree";
	$tree =AbsolutePath("file",$tree);
        open(CFGM00,">$outdir/conf/$index[$i].ctl00")or die $!;
        open(CFGM20,">$outdir/conf/$index[$i].ctl20")or die $!;
$config1 =<< "CFG1";
      seqfile = $seq
      treefile = $tree   
      outfile = /home/jiale/Bioinformatics/Project/06.Ka_Ks/20230404/Results/mlc/$index[$i]\_mlc00
      ndata =1
      
      noisy = 3
      verbose = 1
      runmode = 0
      seqtype = 1
      CodonFreq = 2
      clock = 0
      aaDist = 0
      model = 2
      NSsites = 2
      icode = 0
      Mgene = 0
      fix_kappa = 0
      kappa = .3
      fix_omega = 1
      omega = 1 
      ncatG = 2  
      getSE = 0 
      RateAncestor = 0
      Small_Diff = .45e-6
      cleandata = 1  
      fix_blength = 0    
 
CFG1

 $config2 =<< "CFG2";
      seqfile = $seq
      treefile = $tree
      outfile = /home/jiale/Bioinformatics/Project/06.Ka_Ks/20230404/Results/mlc/$index[$i]\_mlc20
      ndata =1
  
  noisy = 3   
  verbose = 1   
  runmode = 0   
  seqtype = 1   
  CodonFreq = 2  
  clock = 0   
  aaDist = 0 
  model = 2
  NSsites = 2
  icode = 0
  Mgene = 0
  fix_kappa = 0 
  kappa = .3 
  fix_omega = 0
  omega = 1.5
  ncatG = 2
  getSE = 0
  RateAncestor = 0
  Small_Diff = .45e-6
  cleandata = 1
  fix_blength = 0
CFG2
	 print CFGM00 $config1;
	 print CFGM20 $config2;
	print CTL ">$outdir/conf/$index[$i].ctl00\n";
         print CTL ">$outdir/conf/$index[$i].ctl20\n";
         if(-f "$seq")
        {
         push @array1,"$outdir/conf/$index[$i].ctl00";
         push @array2,"$outdir/conf/$index[$i].ctl20";

        }
}
############################################qsub###########################################
        open(SH ,">$outdir/psg.sh")or die $!;
	open(SH2 ,">$outdir/psg2.sh");
        for (my $i=0;$i<@array1;$i++)
        {
         my $ii=$array1[$i];
                $ii=~ s/$outdir//g;
                $ii=~ s/conf//g;
                $ii=~ s/\.ctl00//g;
                my $out1="$outdir/mlc/$ii\_mlc00";
                my $out2="$outdir/mlc/$ii\_mlc20";
                my $out="$outdir/mlc/$ii\_mlc.result";

        print SH "cd $outdir && $codeml  $array1[$i] >$out1\n";
	print SH "$codeml $array2[$i] >$out2\n";
	print SH2 "perl $Bin/Extract_site.pl -zero $out1  -model  $out2 -k z\n";
        }
     `sh   $outdir/psg.sh`;
	`sh $outdir/psg2.sh`;
`cat $outdir/mlc/*.result > $outdir/all.result`;
`perl $Bin/Stat.pl -p $p -result  $outdir/all.result -glist  $tree/glist.txt  -k test -Sdir $mafft/  -od $outdir `;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub AbsolutePath
{		#»ñÈ¡Ö¸¶¨Ä¿Â¼»òÎÄ¼þµÄ¾ö¶¨Â·¾¶
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
  -list    <file>   Input geneID list file, forced
  -tree    <file>   Input the tree file,forced
  -mafft     <dir>   mafft dir
  -p      <float>    p value
  -codeml  <str>    codemls,forced
  -od      <dir>    Dir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
