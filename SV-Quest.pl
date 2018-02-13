#!/usr/bin/perl
use strict;
use List::Util qw(max);
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);

&title;

my $fastq1 = "";
my $fastq2 = "";
my $bam = "";
my $fasta = "";
my $threshould = 0.1;
my $size = 25;
my $call_threshould = "0";
my $coverage = 0.7;
my $coverage2 = 1;
my $outputdir = "summary";
my $output = "InDel.txt";
my $max = "100000";
my $check = "1";
my $wide = "12";
my $Lborder = "300";
my $cpu = "35";
my $seed = "25";
my $seedpenalty  = "1";
my $loose = "275";
my $lowcoverage = "0.1";
my $duplicateRemove = "1";
my $map = "yes";
my $compress = "2";
my $radius = "0.65";
my $slope_threshould = "3";

GetOptions('i=s' => \$bam, 'f=s' => \$fasta, 'm=f' => \$threshould, 'c=i' => \$size, 't=f' => \$call_threshould, 'C=f' => \$coverage, 'y=i' => \$coverage2, 'o=s' => \$output, 'O=s' => \$outputdir, 'L=i' => \$max, 'r=f' => \$check, 'w=i' => \$wide, 'U=i' => \$Lborder, 'p=i' => \$cpu, 's=i' => \$seed, 'M=i' => \$seedpenalty, 'e=i' => \$loose,'1=s' => \$fastq1,'2=s' => \$fastq2,'cov=f' => \$lowcoverage,'d=i' => \$duplicateRemove,'k=s' => \$map,'z=i' => \$compress, 'R=f' => \$radius, 'S=i' => \$slope_threshould);
die "\nMapped.bam or fastq file is required !\n\n\n" if($bam eq "" && $fastq1 eq "" && $check == 1);
die "\ninput fasta fileis required !\n\n\n" if($fasta eq "");
die "\nThe threshould must be 0-1 real number ! default 0.3\n\n\n" unless($threshould <= 1 or $threshould >= 0);
die "\nThe call_size must be plus Integer ! default 30\n\n\n" unless($size >= 1);
die "\nThe call threshould must be plus Integer ! default 30\n\n\n" unless($call_threshould >= 0);
die "\nThe coverage threshould for calling deletion must be 0-1 real number ! default 0.3\n\n\n" unless($coverage <= 1 or $coverage >= 0);
die "\nThe maximum length for searching deletion must be greater than 101 ! default 100000\n\n\n" unless($max > 101);

print " - Mode is $check\n - Threshould is $threshould\n - Wides is $wide\n - Coverage threshould $coverage\n\n";
system("mkdir temp $outputdir");
open TEMP0, ">temp/temp0" or die "cant save \n"&&exit;
print TEMP0 "$fasta";
open TEMP2, ">temp/temp2" or die "cant save \n"&&exit;
print TEMP2 "$size";
open TEMP8, ">temp/temp8" or die "cant save \n"&&exit;
print TEMP8 "$output";
system("mkdir temp/circos_config");
close TEMP0;close TEMP2;close TEMP8;


&return if($check == 1);
system("sleep 3s");
&count if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tProcessing fasta\n" if($check == 1);
&mapping if($check == 1 && $bam eq "");
&sam_convert if($fastq1 eq "" && $check == 1);
system("rm $output") if($check == 1);
system("rm $output $outputdir/log.txt 6F.txt 6R.txt F_mismatch_mountain.txt R_mismatch_mountain.txt temp/F_mismatch_mountain_all.txt temp/R_mismatch_mountain_all.txt $output") unless($check == 1);
&split if($check == 1);#split the reads by the orientation of the reads
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools view of F.sam\n" if($check == 1);
system("samtools view -@ $cpu -bS temp/F.sam | samtools sort -@ $cpu - > temp/F_sorted.bam") if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools view of R.sam\n" if($check == 1);
system("samtools view -@ $cpu -bS temp/R.sam | samtools sort -@ $cpu - > temp/R_sorted.bam") if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools index of F.bam\n" if($check == 1);
system("samtools index temp/F_sorted.bam") if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools index of R.bam\n" if($check == 1);
system("samtools index temp/R_sorted.bam") if($check == 1);
system("rm temp/F.sam temp/R.sam temp/input.sam") if($check == 1);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\t" if($check == 1);
system("samtools mpileup -f temp/reference.fa temp/F_sorted.bam -O -Q 0.1 > temp/F_mpileup") if($check == 1 or $check == 3);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\t" if($check == 1 or $check == 3);
system("samtools mpileup -f temp/reference.fa temp/R_sorted.bam -O -Q 0.1 > temp/R_mpileup") if($check == 1 or $check == 3);
my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts calculation of F&R mismatch\n";
&extract;#calclulate mismatch base from mpileup
&circos if($map eq "yes");
my @now = localtime;print "\nINFO $now[2]:$now[1]:$now[0]\tProgram end.\n";
system("rm temp/temp* *sai .non-redundant.txt");
system("zip -r $output.zip */ InDel.$output $outputdir/log.txt circos* *html $output") if($compress == 1);
print "\nRessults are saved as file name \"$output\"\n\n" unless($compress == 1);
print "\nRessults are saved as file name \"$output\.zip\"\n\n" if($compress == 1);
exit;


#######################################################################################################################################################################################
### SUBROUTINES ###

sub title {
	####################################################################################################################################
	#
	#						SV-Quest version 1.0
	#
	#						Kazuma Uesaka
	#						University of Nagoya
	#						9 February 2018
	#		
	#						A Perl scripts to call SV position from mapped.bam.
	#
	#	
	#		SV Quest: Sensitive Structure Variation detection tool.
	#		Kazuma Uesaka, Hiroshi Kouno, Kazuki Terauchi, Yuichi Fujita, Tatsuo Omata, and Kunio Ihara
	#
	#		
	#						Input:
	#							bam file and reference.fasta for the mapping
	#
	#						Outnput:
	#							Insertion and deletion position printed to STDOUT
	#						
	#						Usage:
	#						perl SV-Quest.pl -f reference.fa -1 forward.fq -2 reverse.fq
	#
	#				The mapped.bam and it's reference.fasta should be included in the same folder.
	#
	#
	#
	#	ver0.4 2017-02-23 version0.4 solid mapped.bam are supported.illumina paired-end fastq are supported.
	#	ver0.5 2017-05-15 coverage ratio from 0.3 to 0.7
	#	ver0.6 check coverage (dont call low coverage region)
	#	ver0.7 mismatch mountain height supported.
	#	ver0.8 redundant call removed.
	####################################################################################################################################


	print "\n\n############################################################################################################################################################\n";
	print "Program: SV-Quest\n";
	print "version 1.0\n\n";
	print "\nUsage:	perl SV-Quest.pl <options>\n\n";
	print "Input/output options:\n\n";
	print "\t -1	input fastq (Required if .bam file doesn't assigned)\n";
	print "\t -2	input pair fastq (for paired end reads)\n";
	print "\t -i	input BAM (Required if .fastq file doesn't assigned)\n";
	print "\t -f	input fasta (Required)\n";
	print "\t -o	output file name (default Indel.txt)\n";
	print "\t -O	output file directory (default summary)\n";

	print "\nMismatch calclulation options:\n\n";
	print "\t -y	Threshould of mismatch ratio coverage. Dont count mismatch ratio if supprting reads are smaller than value (default 1)\n";
	print "\t -C	Threshould of coverage  (default 0.7; stop searching when the candidate read coverage are bigger than 70% of genomic average coverage.)\n";
	print "\t -L	Maximum length (bp) for searching deletion  (default 100000)\n";
	print "\t -U	The border width (bp) of Pf and Pr to call insertion or deletion (default average read length)\n";
	print "\t -c	Calculation size (bp) of mismatch value (default 25)\n";
	print "\t -m	Threshould of mismatch ratio (default 0.1; stop searching when the mismatch ratio was bigger than 0.1)\n";
	print "\t -t	Threshould of call (default 0)\n\n";
	print "\t -S	Mismatch mountain slope threshould to call indel (default 3)\n";
	print "\t -r	Run or skip mapping analysis (default 1)\n";
	print "\t			1;runs\n";
	print "\t			2;skip (re-perform SV-Quest with same mapping parameter. use previous mpileup results)\n";
	print "\t			3;skip and perform mpileup (re-perform SV-Quest with same mapping parameter, but do mpileup)\n";	
	print "\t -w	The size of mismatch value (bp)  (default 12: 12 means ± 12bp mismatch at each postion are included)\n";
	print "\t -p	CPU thread (default 20)\n";
	print "\t -s	seed length (default 25)\n";
	print "\t -M	maximum differences in the seed (default 1)\n";
	print "\t -e	do not put an indel within INT bp towards the ends (default 275)\n";
	print "\t -d	remove duplicate call (default 1)\n";
	print "\t			1;runs\n";
	print "\t			2;skip\n";
	print "\t -k	draw indel map using circos (default yes)\n";
	print "\t			yes;runs\n";
	print "\t			no;skip\n";
	print "\t -R	diameter of map (default 0.65)\n";
	print "\t -z	compress results and related folder with zip (default 2)\n";
	print "\t			1;runs\n";
	print "\t			2;skip\n";
	#print "\t -cov	remove low coverage region. 0.1 means the possible Pf/Pr that have 1/10 times small than average coverage are removed (default 0.1)\n";
	print "############################################################################################################################################################\n\n";
	my @now = localtime;print "\nINFO $now[2]:$now[1]:$now[0]\t";
	print "Starts SV-Quest\n";
	system("sleep 2s");
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub return {
	open INPUT1, "<$fasta" or die "cant open sam file!\n";
	open INPUT2, "<$fasta";
	open (OUT, '>temp/temp0');

	my $line2 = <INPUT2>;
	while (my $line1 = <INPUT1>) {
		my $line2 = <INPUT2>;#1行先読みファイルの入力
		chomp($line1);#改行を除く
		print OUT "$line1\n" if($line1 =~ "\>");#先頭行を出力
		next if($line1 =~ "\>");
		print OUT "$line1";
		print OUT "\n"  if($line2 =~ "\>");#1行先読みファイルを識別に利用している
	}
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub count {
	open INPUT3, "<temp/temp0";
	open (INDEX, '>temp/fasta.index');
	open (REFERENCE, '>temp/reference.fa');
	my $a = 1;
	while (my $line1 = <INPUT3>) {
	chomp($line1);
	$line1 =~ s/\>//;
	my $line2 = <INPUT3>;
	my $contig_size = length($line2);
	print INDEX "$a\t$line1\t$contig_size\n";
	print REFERENCE "\>$a\n$line2";
	$a++;
	}
close INPUT3;close INDEX;
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub mapping {
	die "\nnFastq file is required !\n\n\n" if($fastq1 eq "");
	system("bwa index -p temp/reference.fa -a is temp/reference.fa");
	system("bwa aln -l $seed -k $seedpenalty -n $loose -i $loose -t $cpu temp/reference.fa $fastq1 | bwa samse temp/reference.fa - $fastq1 > temp/input.sam");
	
	next if($fastq2 eq "");
	system("bwa aln -l $seed -k $seedpenalty -n $loose -i $loose -t $cpu temp/reference.fa $fastq2 | bwa samse temp/reference.fa - $fastq2 > temp/inputR2.sam");
	system("sleep 3s");
	
	open INPUT, "<temp/inputR2.sam" or die "cant open temp/inputR2.sam file!\n";
	open (OUT, '>temp/inputR2');
	while (my $line1 = <INPUT>) {
		next if($line1 =~ /^\@/);
		print OUT "$line1";
	}
	system("sleep 3s");
	system("cat temp/inputR2 >> temp/input.sam");
	system("sleep 3s");
	system("rm temp/inputR2 temp/inputR2.sam");
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub sam_convert {
	die "Mapped.bam file is required !\n\n\n" if($bam eq "");
	my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tStarts samtools view\n" if($check == 1);
	system("samtools view -@ $cpu -h $bam > temp/input.sam") if($check == 1);
	my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tDivides input.sam into F.sam and R.sam\n\n" if($check == 1);
}

#-----------------------------------------------------------------------------------------------------------------------------------
sub split {
	open INPUT, "<temp/input.sam" or die "cant open sam file!\n";
	open (OUTF, '>temp/F.sam');
	open (OUTR, '>temp/R.sam');
	open (NAME , '>temp/fasta');
	open (LOG , '>$outputdir/log.txt');
	my $totalp = 0;my $totalm = 0;my $unmap = 0;

	while (my $line1 = <INPUT>) {
		my $head = substr($line1,0,1);
		my @array = split(/\t/, $line1);
		if($head =~ /^\@/){
			print OUTF "$line1";
			print OUTR "$line1";
			my @name2 = "";
			my $head2 = substr($array[1],0,2);
			my @array2 = split(/\:/, $array[1]);
			print NAME "$array2[1]\t" if($head2 =~ /SN/);
		}
		next if($head =~ /^\@/);
		
		print OUTF "$line1" if($array[1] == 0);#forward reads
		print OUTF "$line1" if($array[1] == 256);#forward reads
		$totalp++ if($array[1] == 0);#forward reads
		$totalp++ if($array[1] == 256);#forward reads

		print OUTR "$line1" if($array[1] == 16);#reverse reads
		print OUTR "$line1" if($array[1] == 272);#reverse reads
		$totalm++ if($array[1] == 16);#reverse reads
		$totalm++ if($array[1] == 272);#reverse reads
	}
print LOG "Forward reads $totalp\nReverse reads $totalm\n\n";
close OUTF;close OUTR;close LOG;close NAME;close INPUT;
system("sleep 5s");
}#&split end

#-----------------------------------------------------------------------------------------------------------------------------------
sub extract {
	open CONTIG, "<temp/fasta.index" or die;
	my $counterA = 0;
	my $border = 0;
	while(my $cycle = <CONTIG>){#while1 start
		$counterA++;
		chomp($cycle);
		my @box =  split(/\t/,$cycle);
		#system("mkdir $box[1]");
		print "Chromosome name is $box[1]\n";
		my @now = localtime;print "\nINFO $now[2]:$now[1]:$now[0]\tStarts calculation of $box[1] mismatch\n";
		
		open F, "<temp/F_mpileup" or die ("File can\'t find.") && exit;
		open OUTF, ">temp/mismatch_ratio_F" or die "cant save \n"&&exit;
		my $count = 0;my @CoverageBoxF = "";
		while(my $line = <F>) {
			my @array = split(/\t/, $line);
			$CoverageBoxF[$array[1]] = $array[3];
			next unless($array[0] == $box[0]);
			my @seq = split(//, $array[4]);
			my $mis = 0; my $all = 0;my $nomiscov = 0;
			$count++;
			unless ($count == $array[1]) {
				while ($count<$array[1]) {
					print OUTF "$count\t0\t0\n";$count++;
				}
			}
			foreach (@seq) {#count mismatched bases
				$mis++ if ((/A/i) or (/T/i) or (/G/i) or (/C/i) or (/N/i));
				$all++ if ((/\./) or (/\,/) or (/A/i) or (/T/i) or (/G/i) or (/C/i) or (/N/i));
			}
			$CoverageBoxF[$array[1]] = $all - $mis;#coverageF
			my $mismuch_ratioF = 0;
			$mismuch_ratioF = $mis / $all unless($all == "0");
			$mismuch_ratioF = 0 if($mismuch_ratioF < $threshould);
			$mismuch_ratioF = 0 if($all <= $coverage2);
			$nomiscov = $all - $mis;#ver1.7追加部位
			print OUTF "$array[1]\t$nomiscov\t$mismuch_ratioF\n";#position, coverage, mismatch ratio of F #ver1.7修正部位
		}
		close F; close OUTF;


		open R, "<temp/R_mpileup" or die ("File can\'t find.") && exit;
		open OUTM, ">temp/mismatch_ratio_R" or die "cant save \n"&&exit;
		my $count = 0;my @CoverageBoxR = "";
		while(my $line = <R>) {
			my @array = split(/\t/, $line);
			$CoverageBoxR[$array[1]] = $array[3];#coverageR
			next unless($array[0] == $box[0]);
			my @seq = split(//, $array[4]);
			my $mis = 0; my $all = 0;my $nomiscov = 0;
			$count++;
			unless ($count == $array[1]) {
				while ($count<$array[1]) {
					print OUTM "$count\t0\t0\n";$count++;
				}
			}
			foreach (@seq) {#count mismatched bases
				$mis++ if ((/A/i) or (/T/i) or (/G/i) or (/C/i) or (/N/i));
				$all++ if ((/\./) or (/\,/) or (/A/i) or (/T/i) or (/G/i) or (/C/i) or (/N/i));
			}
		$CoverageBoxR[$array[1]] =  $all - $mis;
		my $mismuch_ratioR = 0;
		$mismuch_ratioR = $mis / $all unless($all == "0");
		$mismuch_ratioR = 0 if($mismuch_ratioR < $threshould);
		$mismuch_ratioR = 0 if($all <= $coverage2);
		$nomiscov = $all - $mis;#ver1.7追加部位
		print OUTM "$array[1]\t$nomiscov\t$mismuch_ratioR\n";#position, coverage, mismatch ratio of R #ver1.7修正部位
		}
		close R; close OUTM;
		print "mismatch ratio calculation end\n";
		
		my $half = int($size / 2);
		open F1, "<temp/mismatch_ratio_F" or die;
		open F2, "<temp/mismatch_ratio_F" or die;
		open FOUT1, ">temp/mismatch_value_F" or die;
		my @cov = (); my @misratioF = (); my @pos = ();my $i = 1;
	
		while (my $f1 = <F1>){
			chomp($f1);
			my @array = split(/\t/, $f1);
			$pos[$i] = $i;
			$cov[$i] = $array[1];
			$misratioF[$i] = $array[2];
			$i++;
		}
		my $ii = 1;
		while (my $f2 = <F2>){
			my @array2 = split(/\t/,$f2);
			my $right = $ii + $wide;
			my $left = $ii - $wide;
			$left = 1 if($left < 1);
			my @group = ();my $mismuch_valueF = 0;
			for (my $iii = $left; $iii <= $right;$iii++){#calclulate mismatch value in the loop
				$mismuch_valueF = $mismuch_valueF + $misratioF[$iii];#calclulate mismatch value of F
			}
		print FOUT1 "$pos[$ii]\t$array2[1]\t$mismuch_valueF\n";
		$ii++;
		}
		close F1;close F2;close FOUT1;

		open R1, "<temp/mismatch_ratio_R" or die;
		open R2, "<temp/mismatch_ratio_R" or die;
		open ROUT1, ">temp/mismatch_value_R" or die;
		my @cov = (); my @misratioR = (); my @pos = ();my $i = 1;
		while (my $r1 = <R1>){
			chomp($r1);
			my @array = split(/\t/, $r1);
			$pos[$i] = $i;
			$cov[$i] = $array[1];
			$misratioR[$i] = $array[2];
			$i++;
		}

		my $ii = 1;
		while (my $r2 = <R2>){
			my @array2 = split(/\t/,$r2);
			my $left = $ii - $wide;
			my $right = $ii + $wide;
			$left = 1 if($left < 1);
			my @group = (); my $mismuch_valueR = 0; my $region = 0;
			for (my $iii = $right; $iii >= $left;$iii--){
				$mismuch_valueR = $mismuch_valueR + $misratioR[$iii];#calclulate mismatch value of R
			}
		print ROUT1 "$pos[$ii]\t$array2[1]\t$mismuch_valueR\n";
		$ii++;
		}
		close R1;close R2;close ROUT1;
		print "mismatch value calculation end\n";
		
		
		open (LOG , '>>$outputdir/log.txt');
		open PLUSM, "temp/mismatch_value_F" or die;
		open MINUSM, "<temp/mismatch_value_R" or die;
		my $i = 1; my @plus_cov = "";#average read count of F
		my $p = <PLUSM>;
		while (my $p = <PLUSM>){
			chomp($p); my @boxP = split(/\t/,$p);
			$plus_cov[$i] = $boxP[1] if($i <= $box[2]);
			$i++;
		}
		
		my $i = 1; my @minus_cov = ""; my @cov = "";my $total = 0.1;#average read count of R
		my $m = <MINUSM>;
		while (my $m = <MINUSM>){
			chomp($m); my @boxM = split(/\t/,$m);
			$minus_cov[$i] = $boxM[1] if($i <= $box[2]);
			$cov[$i] = $plus_cov[$i] + $minus_cov[$i]  if($i <= $box[2]);
			$total = $total + $cov[$i];#total read count
			$i++;
		}
		
		
		$total = $total / $box[2];#average read count
		print LOG "average read count in $box[1] is $total\n";
		
		open F3, "<temp/mismatch_value_F" or die;
		open F4, "<temp/mismatch_value_F" or die;
		open FOUT3, ">temp/F_mismatch_mountain.txt" or die;
		open FOUT4, ">>temp/F_mismatch_mountain_all.txt" or die;
		
		my $counter = 0;my $sum = 0;my $start = 0;
		my $f4 = <F4>;
		my $mismatch_wideness = 0;
		while (my $f3 = <F3>){
			my $f4 = <F4>;
			my @array3 = split(/\t/,$f3);
			my @array4 = split(/\t/,$f4);
			chomp($array3[2]);chomp($array4[2]);
			my $CoverageThisPositionF = $CoverageBoxF[$array3[0]] / $total / 2;
			#next if($CoverageThisPositionF < $lowcoverage);
			if ($array3[2] > 0){
				$counter++;
				$sum += $array3[2];
				$start = $array3[0] if($counter == 1);
				$mismatch_wideness++;
			}
			if ($array3[2] > 0 && $array4[2] == 0){
				my $endness = $mismatch_wideness;
				$mismatch_wideness = 0;
				my $end = $array3[0];
				my $center = ($start + $end) / 2;
				#my $sizew = $end - $start;
				$center = int($center);
				my $slope = 0;
				$slope = $sum / $endness unless($endness == 0);
				print FOUT3 "$center\t$sum\t$endness\n" if($sum >= $call_threshould && $slope >= $slope_threshould);#version1.7修正部位
				print FOUT4 "$box[0]\t$center\t$sum\t$endness\n" if($sum >= $call_threshould && $slope >= $slope_threshould);#version1.7修正部位
				$counter = 0;$sum = 0;
			}
		}
		close F3;close F4;


		open R3, "<temp/mismatch_value_R" or die;
		open R4, "<temp/mismatch_value_R" or die;
		open ROUT3, ">temp/R_mismatch_mountain.txt" or die;
		open ROUT4, ">>temp/R_mismatch_mountain_all.txt" or die;

		my $counter = 0;my $sum = 0;my $start = 0;
		my $r4 = <R4>;
		my $mismatch_wideness = 0;
		while (my $r3 = <R3>){
			my $r4 = <R4>;
			my @array3 = split(/\t/,$r3);
			my @array4 = split(/\t/,$r4);
			chomp($array3[2]);chomp($array4[2]);
			my $CoverageThisPositionR = $CoverageBoxR[$array3[0]] / $total / 2;
			#next if($CoverageThisPositionR < $lowcoverage);
			if ($array3[2] > 0){
				$counter++;
				$sum += $array3[2];
				$start = $array3[0] if($counter == 1);
				$mismatch_wideness++;
			}
		
			if ($array3[2] > 0 && $array4[2] == 0){
				my $endness = $mismatch_wideness;
				$mismatch_wideness = 0;
				my $end = $array3[0];
				my $center = ($start + $end) / 2;
				$center = int($center);
				#my $sizew = $end - $start;
				my $slope = 0;
				$slope = $sum / $endness unless($endness == 0);
				print ROUT3 "$center\t$sum\t$endness\n" if($sum >= $call_threshould && $slope >= $slope_threshould);#version1.7修正部位
				print ROUT4 "$box[0]\t$center\t$sum\t$endness\n" if($sum >= $call_threshould && $slope >= $slope_threshould);#version1.7修正部位
				$counter = 0;$sum = 0;
			}
		}
		close R3;close R4;
		system("sleep 5s");
		
		my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\tMismatch calculation of $box[1] was end. Now starts indel call of $box[1] from mismatch information\n";
		#call Insertion and small deletion candidate
		open PLUS2, "<temp/F_mismatch_mountain.txt" or die;
		open MINUS2, "<temp/R_mismatch_mountain.txt" or die;
		open OUTPUT3I, ">temp/In.txt" or die;
		open OUTPUT3IC, ">>$output" or die;

		print OUTPUT3I "\#chrom\tstart\tend\tmismatch\n";
		print OUTPUT3IC "\#chrom\ttype\tpositon\tmismatch\n" if($counterA eq 1);
		my @mismuchMP =();my @mismuchMM =();my $i = 0;
		while (my $minus2 = <MINUS2>){
			my @array2 = split(/\t/,$minus2);chomp($array2[1]);
			$mismuchMP[$i] = $array2[0];
			$mismuchMM[$i] = $array2[1];
			$i++;
		}
		my $length = @mismuchMP;
		
		while (my $plus2 = <PLUS2>){
			
			my @array1 = split(/\t/,$plus2); chomp($array1[1]);
			for (my $i = 0; $i < $length; $i++){
			
				my $right =  $mismuchMP[$i] + $Lborder;#≤ 300 bp
				my $left =  $mismuchMP[$i] -  $Lborder;#≤ 300 bp
				$array1[0] =~ s/\s+//g;
				
				if ($array1[0] >= $left && $array1[0] <= $right){#Both Pf and Pr are found (≤ 300 bp)
					my $value = ($array1[1] + $mismuchMM[$i]) / 2;
					$value  = int($value);
					print OUTPUT3I "$box[1]\t$array1[0]\t$mismuchMP[$i]\t$value\n" if($array1[0] <= $mismuchMP[$i]);
					print OUTPUT3IC "$box[1]\tType_I_SV\t$array1[0]\-$mismuchMP[$i]\t$value\n" if($array1[0] <= $mismuchMP[$i]);
					
					print OUTPUT3I "$box[1]\t$mismuchMP[$i]\t$array1[0]\t$value\n" if($array1[0] > $mismuchMP[$i]);
					print OUTPUT3IC "$box[1]\tType_I_SV\t$mismuchMP[$i]\-$array1[0]\t$value\n" if($array1[0] > $mismuchMP[$i]);
					
				}
			}
		}
		
		#call large deletion candidate
		open COVF, "<temp/mismatch_ratio_F" or die;
		open COVR, "<temp/mismatch_ratio_R" or die;
		my @coverageF = "";	
		while(my $f = <COVF>){
		my @fbox = split(/\t/,$f);
		$coverageF[$fbox[0]] = $fbox[1];
		}
		my @coverageR = "";	
		while(my $r = <COVR>){
			my @rbox = split(/\t/,$r);
			$coverageR[$rbox[0]] = $rbox[1];
		}

		open PLUS2, "<temp/F_mismatch_mountain.txt" or die;
		open OUTPUT4D, ">temp/Del.txt" or die;
		open OUTPUT4DC, ">>$output" or die;
		print OUTPUT4D "\#chrom\tstart\tend\trelative_coverage\n";
			while (my $checkp = <PLUS2>){
			chomp($checkp);
			my @inputp = split(/\t/,$checkp);
			my $right =  $inputp[0] + $max;
			my $left = $inputp[0] + $Lborder;##Pf (read length ≥ , ≤ 100 kbp)
			
			for (my $i = 0; $i < $length; $i++){
				if ($left <= $mismuchMP[$i] && $mismuchMP[$i] <= $right){

					my $sumC = 0;my $ii = 0;
					for (my $k = $inputp[0]; $k <= $mismuchMP[$i]; $k++){
						$sumC = $sumC + $coverageF[$k] + $coverageR[$k];
						$ii++;
					}
					my $AveC = $sumC / $ii;
					my $read_ratio = 0;
					$read_ratio = $AveC / $total;
					print OUTPUT4D "$box[1]\t$inputp[0]\t$mismuchMP[$i]\t$read_ratio\n" if($read_ratio < $coverage);
					print OUTPUT4DC "$box[1]\tType_II_SV\t$inputp[0]\-$mismuchMP[$i]\t$read_ratio\n" if($read_ratio < $coverage);
				}
			}
		}
		my @now = localtime;print "INFO $now[2]:$now[1]:$now[0]\t$box[1] end.\n";
		system("sleep 2s");
		if($duplicateRemove == 1){
			open LAST1, "<temp/In.txt" or die;
			open LAST2, "<temp/In.txt" or die;
			open OUTP, ">temp/In_non-redundant.txt" or die;
			open OUTPP, ">>.non-redundant.txt" or die;
			print OUTP "\#chrom\tstart\tend\tmismatch\n";
			print OUTPP "\#chrom\ttype\positon\tmismatch\n";
			my $input1 = <LAST1>;
			my $input2 = <LAST2>;my $input2 = <LAST2>;
			while($input1 = <LAST1>){
				$input2 = <LAST2>;
				next if($input1 eq "");
				chomp($input1);
				chomp($input2);
				my @lbox1 = split(/\t/,$input1);
				my @lbox2 = split(/\t/,$input2);
				if(($lbox1[1] + 300) > $lbox2[1]){
					my $firstP = $lbox1[1];
					my $endP = $lbox2[2];
					$lbox1[3] = $lbox2[3] if($lbox1[3] < $lbox2[3]);
					for(my $e=1; $e < 1000; $e++){
						$input1 = <LAST1>;
						$input2 = <LAST2>;
						chomp($input1);
						chomp($input2);
						my @lbox3 = split(/\t/,$input1);
						my @lbox4 = split(/\t/,$input2);
						$lbox1[3] = $lbox4[3] if(($lbox1[1] + 300) > $lbox4[1] && $lbox1[3] < $lbox4[3]);
						$endP = $lbox4[2] if(($lbox1[1] + 300) > $lbox4[1]);
						print OUTP "$lbox1[0]\t$firstP\t$lbox1[2]\t$lbox1[3]\n" if((($lbox1[1] + 300) < $lbox4[1]) or ($lbox4[1] eq ""));
						print OUTPP "$lbox1[0]\tType_I_SV\t$firstP\-$lbox1[2]\t$lbox1[3]\n" if((($lbox1[1] + 300) < $lbox4[1]) or ($lbox4[1] eq ""));
						last if((($lbox1[1] + 300) < $lbox4[1]) or ($lbox4[1] eq ""));
					}#for
				}elsif((($lbox1[1] + 300) < $lbox2[1]) or ($lbox2[1] eq "")){
					print OUTP "$input1\n";
					print OUTPP "$lbox1[0]\tType_I_SV\t$lbox1[1]\-$lbox1[2]\t$lbox1[3]\n";
				}
			}
		}
		
		if($duplicateRemove == 1){
			open LAST3, "<temp/Del.txt" or die;
			open LAST4, "<temp/Del.txt" or die;
			open OUTD, ">temp/Del_non-redundant.txt" or die;
			open OUTDR, ">>.non-redundant.txt" or die;
			print OUTD "\#chrom\tstart\tend\trelative_coverage\n";
			my $input1 = <LAST3>;
			my $input2 = <LAST4>;my $input2 = <LAST4>;
			while($input1 = <LAST3>){
				$input2 = <LAST4>;
				next if($input1 eq "");
				chomp($input1);
				chomp($input2);
				my @lbox1 = split(/\t/,$input1);
				my @lbox2 = split(/\t/,$input2);
				if($lbox1[1] == $lbox2[1] || $lbox1[2] == $lbox2[2]){
					
					my $firstP = $lbox1[1];
					my $endP = $lbox2[2];
					$lbox1[3] = $lbox2[3] if($lbox1[3] < $lbox2[3]);
					for(my $e=1; $e < 1000; $e++){
						$input1 = <LAST3>;
						$input2 = <LAST4>;
						chomp($input1);
						chomp($input2);
						my @lbox3 = split(/\t/,$input1);
						my @lbox4 = split(/\t/,$input2);
						$lbox1[3] = $lbox4[3] if(($lbox1[1] + 300) > $lbox4[1] && $lbox1[3] < $lbox4[3]);
						$endP = $lbox4[2] if(($lbox1[1] + 300) > $lbox4[1]);
						print OUTD "$lbox1[0]\t$firstP\t$lbox1[2]\t$lbox1[3]\n" if((($lbox1[1] != $lbox4[1]) && ($lbox1[2] != $lbox4[2]))  or ($lbox4[1] eq ""));
						print OUTDR "$lbox1[0]\tType_II_SV\t$firstP\-$lbox1[2]\t$lbox1[3]\n" if((($lbox1[1] != $lbox4[1]) && ($lbox1[2] != $lbox4[2])) or ($lbox4[1] eq ""));
						last if((($lbox1[1] != $lbox4[1]) && ($lbox1[2] != $lbox4[2])) or ($lbox4[1] eq ""));
					}#for
				}else{
					print OUTD "$input1\n";
					print OUTDR "$lbox1[0]\tType_II_SV\t$lbox1[1]\-$lbox1[2]\t$lbox1[3]\n";
				}
			}
		}
	}#while1end
}#&extract end

#-----------------------------------------------------------------------------------------------------------------------------------
sub circos {

open INPUT0, "<temp/fasta.index" or die "cant open fasta.index file!\n";

my @number = "";
my @name = "";
my @length = "";
my $k = 1;
while (my $line0 = <INPUT0>){
	chomp($line0);
	my @array0 = split(/\t/, $line0);
	$number[$k] = $array0[0];
	$name[$k] = $array0[1];
	$length[$k] = $array0[2];
	$k++;
}
$k--;
print "$k\n";
open SUMMARY, ">summary.html" or die "cant open html file!\n";
open CIRCOSORF, ">chr_all.txt" or die "cant open txt file!\n";


my $top2 = 0;
for(my $i = 1; $i <= $k; $i++){
	if($i == 1){
	my $merge0 = "<HTML><HEAD><TITLE>Summary</TITLE></HEAD><BODY bgcolor=\"#FAFAFA\"><PRE><CENTER><br><br><br>\"SV-Quest: a program to detect SV\"<br><br><br><font size=\"4\"><B>Prediction Summary</B></font><br>";
	print SUMMARY "$merge0\n";
	}
	open INPUT3, "<temp/mismatch_value_F";
	open OUT3A, ">F_average_mismatch_value.txt";
	my $max = 0;
	my $counter = 1;
	my $sum = 0;
	my $nowF = 0;
	my $unit = int($length[$i] / 3 / 1000);
	$unit = 10 if($length[$i] < 10000);
	my @boxF1 = "";my @boxF2 = "";my @boxF3 = "";my @boxF4 = "";
	print "unit is $unit\n";
	while (my $line3 = <INPUT3>){
		chomp($line3);
		my @array3 = split(/\t/, $line3);
		$max = $array3[2] if($max < $array3[2]);
		$sum = $sum + $array3[1];
		if($counter == $unit){
			my $start = $array3[0] -$unit + 1;
			my $average = $sum / $unit;
			print OUT3A "$name[$i]\t$start\t$array3[0]\t$max\n";
			$boxF1[$nowF] = $name[$i];
			$boxF2[$nowF] = $start;
			$boxF3[$nowF] = $array3[0];
			$boxF4[$nowF] = $average;
			$max = 0;
			$sum = 0;
			$counter = 0;
			$nowF++;
		}elsif($array3[0] == $length[$i]){
			my $max =  $array3[0] -$counter;
			my $start = $array3[0] - $counter;
			my $average = $sum / $counter;
			print OUT3A "$name[$i]\t$start\t$array3[0]\t$max\n";
			$boxF1[$nowF] = $name[$i];
			$boxF2[$nowF] = $start;
			$boxF3[$nowF] = $array3[0];
			$boxF4[$nowF] = $average;
		}
		$counter++;
	}
	
	
	open INPUT4, "<temp/mismatch_value_R" or die "cant open R_mismatch value file!\n";
	open OUT4A, ">R_average_mismatch_value.txt" or die "cant save file!\n";
	my $max = 0;
	my $counter = 1;
	my $sum = 0;
	my $nowR = 0;
	my @boxR1 = "";my @boxR2 = "";my @boxR3 = "";my @boxR4 = "";
	while (my $line4 = <INPUT4>){
		chomp($line4);
		my @array4 = split(/\t/, $line4);
		$max = $array4[2] if($max < $array4[2]);
		$sum = $sum + $array4[1];
		if($counter == $unit){
			my $start = $array4[0] - $unit + 1;
			my $average = $sum / $unit;
			print OUT4A "$name[$i]\t$start\t$array4[0]\t$max\n";
			$boxR1[$nowR] = $name[$i];
			$boxR2[$nowR] = $start;
			$boxR3[$nowR] = $array4[0];
			$boxR4[$nowR] = $average;
			$max = 0;
			$sum = 0;
			$counter = 0;
			$nowR++;
		}elsif($array4[0] == $length[$i]){
			my $max =  $array4[0] -$counter;
			my $average = $sum / $counter;
			my $start = $array4[0] - $counter;
			print OUT4A "$name[$i]\t$start\t$array4[0]\t$max\n";
			$boxR1[$nowR] = $name[$i];
			$boxR2[$nowR] = $start;
			$boxR3[$nowR] = $array4[0];
			$boxR4[$nowR] = $average;
		}
		$counter++;
	}
	
	open OUT, ">FR_average_coverage.txt" or die "cant open file!\n";
	open LOG, ">log_coverage.txt" or die "cant open file!\n";
	my $top = 0;
	my $L = @boxF1;
	for(my $loop = 0; $loop < $L; $loop++){
		next unless($boxF2[$loop] > 0 or $boxR2[$loop] > 0);
		my $plus = $boxF4[$loop] + $boxR4[$loop];
		$plus = 0.1 if($plus == 0 or $plus eq "");
		my $plus2 = log($plus);
		$top = $plus if($plus > $top);
		$top2 = $plus2 if($plus2 > $top2);
		print OUT "$name[$i]\t$boxF2[$loop]\t$boxF3[$loop]\t$plus\n" if($boxF2[$loop] > 0 && $boxF3[$loop] > 0);
		print LOG "$name[$i]\t$boxF2[$loop]\t$boxF3[$loop]\t$plus2\n" if($boxF2[$loop] > 0 && $boxF3[$loop] > 0);
	}
	
	my $merge1A = "<font size=\"4\"><br><b>Type_I_SV</b><br></font><table border=\"1\" width=\"350\" cellspacing=\"0\"><tr align=\"center\"><th bgcolor=\"d3d3d3\">positon</th><th bgcolor=\"d3d3d3\">mismatch value</th></tr>";
	my $merge1B = "<font size=\"4\"><b>Type_II_SV</b><br></font><table border=\"1\" width=\"350\" cellspacing=\"0\"><tr align=\"center\"><th bgcolor=\"d3d3d3\">positon</th><th bgcolor=\"d3d3d3\">length (bp)</th><th bgcolor=\"d3d3d3\">relative coverage</th></tr>";
	my $merge1C = "<table border=\"1\" width=\"350\" cellspacing=\"0\"><tr align=\"center\"><th bgcolor=\"d3d3d3\">chr</th><th bgcolor=\"d3d3d3\">length (bp)</th><th bgcolor=\"d3d3d3\">Type_I_SV</th><th bgcolor=\"d3d3d3\">Type_II_SV</th></tr>" if($i == 1);
	print SUMMARY "$merge1C\n" if($i == 1);
	system("cp .non-redundant.txt $output") if($duplicateRemove == 1);
	system("rm .non-redundant.txt") unless($duplicateRemove == 1);
	open INPUTI, "<$output" or die "cant open R_mismatch value file!\n";
	open OUTPUTI, ">temp/In.txt" or die "cant open R_mismatch value file!\n";
	open OUTPUTD, ">temp/Del.txt" or die "cant open R_mismatch value file!\n";
	
	my $indel = <INPUTI>;
	my $Intotal = 0;
	my $Deltotal = 0;
	while (my $indel = <INPUTI>){
		chomp($indel);
		my @indel = split(/\t/, $indel);
		my @position = split(/\-/, $indel[2]);
		next unless($indel[0] eq $name[$i]);
		my $width = abs($position[1] - $position[0]);
		my $value = int($indel[3] * 1000) / 1000;
		my $temp1A = "</tr><tr align=\"center\"><td>$indel[2]</td><td>$indel[3]</td></tr>" if($indel[1] eq "Type_I_SV");
		my $temp1B = "</tr><tr align=\"center\"><td>$indel[2]</td><td>$width</td><td>$value</td></tr>" if($indel[1] eq "Type_II_SV");
		$merge1A = $merge1A . $temp1A;
		$merge1B = $merge1B . $temp1B;
		$Intotal++ if($indel[1] eq "Type_I_SV");
		$Deltotal++ if($indel[1] eq "Type_II_SV");
		print OUTPUTI "$name[$i]\t$position[0]\t$position[1]\t$indel[3]\n" if($indel[1] eq "Type_I_SV");
		print OUTPUTD "$name[$i]\t$position[0]\t$position[1]\t$indel[3]\n" if($indel[1] eq "Type_II_SV");
	}
	my $temp = "</table>";
	$merge1A = $merge1A . $temp;
	$merge1B = $merge1B . $temp;
	#
	my $temp1C = "</tr><tr align=\"center\"><td><a href = \"$name[$i]/$name[$i].html\" target=\"_blank\">$name[$i]</a></td><td>$length[$i]</td><td>$Intotal</td><td>$Deltotal</td></tr>";
	
	print SUMMARY "$temp1C\n";
	print SUMMARY "$temp\n" if($i == $k);
	

open OUT5A, ">IDEOGRAM.CONF" or die "cant save file!\n";
my $textA = <<'EOS';
<ideogram>
<spacing>
break   = 0.5r
default = 0u #環状になる
</spacing>
<<include ideogram.position.conf>>
<<include ideogram.label.conf>>
<<include bands.conf>>	
radius*       = 0.825r
</ideogram>
EOS
print OUT5A "$textA\n";
	
	
open OUT5B, ">IDEOGRAM.LABEL.CONF" or die "cant save file!\n";
my $textB = <<'EOS';
show_label       = yes
label_font       = bold
label_radius = dims(ideogram,radius_outer) + 0.05r
label_size       = 30
label_parallel   = yes
label_case       = lower
label_format     = eval(sprintf("%s",var(label)))
EOS
print OUT5B "$textB\n";
	
open OUT5C, ">TICKS.CONF" or die "cant save file!\n";
my $textC = <<'EOS';
show_ticks          = yes
show_tick_labels    = yes
<ticks>
radius           = dims(ideogram,radius_outer)
orientation      = out
label_multiplier = 1e-4
color            = black
size             = 15p
thickness        = 1p
label_offset     = 5p
<tick>
spacing        = .1u
show_label     = no
</tick>
<tick>
spacing        = .5u
show_label     = yes 
label_size     = 16p
format         = %.1f
</tick>
<tick>
spacing        = 5u
show_label     = yes
label_size     = 30p
format         = %d
</tick>
</ticks>
EOS
print OUT5C "$textC\n";
	
open OUT5D, ">BANDS.CONF" or die "cant save file!\n";
my $textD = <<'EOS';
	show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 0
EOS
print OUT5D "$textD\n";

open OUT5E, ">IDEOGRAM.POSITION.CONF" or die "cant save file!\n";
my $textE = <<'EOS';
radius           = 0.775r
thickness        = 10p #これが環の太さパラメータ
fill             = yes
fill_color       = black
stroke_thickness = 10
stroke_color     = black
EOS
print OUT5E "$textE\n";

open OUT5F, ">genome.txt" or die "cant save file!\n";
my $textF = <<"EOS";
chr - $name[$i] $name[$i] 0 $length[$i] $name[$i]
EOS
print OUT5F "$textF\n";
print CIRCOSORF "$textF";

my $pwd = "";
system("pwd > $pwd");
my $moji = 1;
my $tick1 = length($length[$i]);
for(my $s = 1; $s < $tick1;$s++){
	$moji = $moji . 0;
}
$moji = $moji / 10;
print "$moji\n";

open OUT5G, ">MAIN.CONF" or die "cant save file!\n";
my $textG = <<"EOS";
# 1.2 IDEOGRAM LABELS, TICKS, AND MODULARIZING CONFIGURATION
karyotype = genome.txt
chromosomes_units = $moji
<<include ideogram.conf>>
<<include ticks.conf>>

<plots>

<plot>
type      = line
thickness = 0.1
max_gap = 1u
file    = FR_average_coverage.txt
color   = while
min     = 0
max     = $top
r0      = 0.90r
r1      = 0.98r
fill_color = 76,76,76
<backgrounds>
<background>
color     = white
y0        = 0.006
</background>
<background>
color     = white
y1        = 0.002
</background>
</backgrounds>
<axes>
<axis>
color     = white
thickness = 1
spacing   = 10r
</axis>
</axes>
</plots>

<plot>
type      = line
thickness = 0.1
max_gap = 1u
file    = F_average_mismatch_value.txt
color   = while
min     = 0
max     = 20
r0      = 0.83r
r1      = 0.90r
fill_color = blue
<backgrounds>
<background>
color     = white
y0        = 0.006
</background>
<background>
color     = white
y1        = 0.002
</background>
</backgrounds>
<axes>
<axis>
color     = white
thickness = 1
spacing   = 10r
</axis>
</axes>
</plots>

<plot>
type      = line
thickness = 0.1
max_gap = 1u
file    = R_average_mismatch_value.txt
color   = while
min     = 0
max     = 20
r0      = 0.76r
r1      = 0.83r
fill_color = red
<backgrounds>
<background>
color     = white
y0        = 0.006
</background>
<background>
color     = white
y1        = 0.002
</background>
</backgrounds>
<axes>
<axis>
color     = white
thickness = 1
spacing   = 10r
</axis>
</axes>
</plots>


<plot>
type         = scatter
file    = temp/In.txt
glyph        = triangle
glyph_size   = 20
orientation  = out
fill_color   = black
r0      = 0.72r
r1      = 0.72r
</plots>


<plot>
type      = heatmap
file    = temp/Del.txt
scale_log_base = 2
min          = 0.00001
color        = red
r0      = 0.69r
r1      = 0.71r
</plots>

</plots>

<image>
<<include etc/image.conf>>
</image>
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
EOS
print OUT5G "$textG\n";

system("sleep 3s");
system("circos -conf MAIN.CONF");
open OUTH, ">$name[$i].html" or die "cant open html file!\n";
my $merge0 = "<HTML><HEAD><TITLE>Summary</TITLE></HEAD><BODY bgcolor=\"#FAFAFA\"><PRE><CENTER><br><font size=\"4\"><B>$name[$i]</B></font><br>$length[$i] bp<br>";
my $merge2 = "<img src=\"circos.png\" width=\"1000\" height=\"1000\"></CENTER></BODY></HTML>";
print OUTH "$merge0\n$merge1A\n\n$merge1B\n\n$merge2\n";
print SUMMARY "$merge2\n" if($i == $k);
system("sleep 1s");
system("cp TICKS.CONF TICKS2.CONF") if($i == 1);
system("cat log_coverage.txt >> FR_all_coverage.txt");
system("cat temp/In.txt >> In_all.txt");
system("cat temp/Del.txt >> Del_all.txt");
system("cat F_average_mismatch_value.txt >> F_all_mismatch.txt");
system("cat R_average_mismatch_value.txt >> R_all_mismatch.txt");
system("sleep 1s");
#system("mv *_average_* $name[$i].html circos.png circos.svg MAIN.CONF IDEOGRAM.CONF IDEOGRAM.LABEL.CONF IDEOGRAM.POSITION.CONF BANDS.CONF TICKS.CONF genome.txt In.txt Del.txt *mismatch_value.txt FR_average_coverage* *_average_mismatch_value.txt $name[$i]/");
system("sleep 1s");
}# for end


open INPUT, "<chr_all.txt" or die "cant save file!\n";#circos; the maximum is currently set at [200].
open SMALL, ">chr_all_small.txt" or die "cant save file!\n";
my $n = 0;
while (my $decrease = <INPUT>){
	print SMALL "$decrease\n" if($n <= 199);
	$n++;
}
system("mv chr_all_small.txt chr_all.txt");


open INPUT, "<FR_all_coverage.txt" or die "cant save file!\n";#circos; the maximum plot is currently set at [25000].
open SMALL, ">FR_all_coverage_small.txt" or die "cant save file!\n";
my @array1 = "";
my @array2 = "";
my @array3 = "";
my @array4 = "";
my $n = 0;
while (my $decrease = <INPUT>){
	chomp($decrease);
	my @decreasebox = split(/\t/, $decrease);
	$array1[$n] = $decreasebox[0];
	$array2[$n] = $decreasebox[1];
	$array3[$n] = $decreasebox[2];
	$array4[$n] = $decreasebox[3];
	$n++;
}
if($n >= 25000){
	my $ratio = $n / 25000;
	$ratio = 2 if($ratio < 2);
	int($ratio);
	my $ratio2 = $ratio * 2;
	for(my $q = 0; $q <= $n; $q++){
		print SMALL "$array1[$q]\t$array2[$q]\t$array3[$q]\t$array4[$q]\n" unless($array1[$q] eq $array1[$q+$ratio2]);
		if($array1[$q] eq $array1[$q+$ratio2]){
		print SMALL "$array1[$q]\t$array2[$q]\t$array3[$q+$ratio2]\t$array4[$q]\n";
		$q = $q + $ratio2 - 1;
		}
	}
}
system("mv FR_all_coverage_small.txt FR_all_coverage.txt") if($n >= 25000);


open INPUT, "<F_all_mismatch.txt" or die "cant open file!\n";#circos; the maximum plot is currently set at [25000].
open INPUT2, "<F_all_mismatch.txt" or die "cant open file!\n";
open OUTFALL, ">F_all_mismatch_small.txt" or die "cant save file!\n";
my $decrease2 = <INPUT2>;
my $point = 0;
while (my $decrease = <INPUT>){
	$decrease2 = <INPUT2>;
	chomp($decrease);
	chomp($decrease2);
	my @decreasebox = split(/\t/, $decrease);
	my @decreasebox2 = split(/\t/, $decrease2);
	
	$point++ if($decreasebox[3] == 0);
	unless($decreasebox[3] == 0){
		$point = 0;
		my $half = ($decreasebox[2] - $decreasebox[1]) / 2;
		my $halfA = $decreasebox[1] + $half;
		$halfA = int($halfA);
		my $halfB = $halfA + 1;
		print OUTFALL "$decreasebox[0]\t$decreasebox[1]\t$halfA\t0\n";
		print OUTFALL "$decreasebox[0]\t$halfB\t$decreasebox[2]\t$decreasebox[3]\n";
		print OUTFALL "$decreasebox2[0]\t$decreasebox2[1]\t$decreasebox2[2]\t$decreasebox2[3]\n" if($decreasebox2[3] == 0);
	}
	print OUTFALL "$decrease\n" if($point >= 10 && $point <= 15);
	print OUTFALL "$decrease\n" if($point >= 100 && $point <= 130);
	print OUTFALL "$decrease\n" if($point >= 400 && $point <= 430);
	print OUTFALL "$decrease\n" if($point >= 1500 && $point <= 1530);
	print OUTFALL "$decrease\n" if($point >= 2000 && $point <= 2030);
	print OUTFALL "$decrease\n" if($point >= 2500 && $point <= 1530);
	print OUTFALL "$decrease\n" if($point >= 3000 && $point <= 3030);
	print OUTFALL "$decrease\n" if($point >= 4000 && $point <= 4030);
	print OUTFALL "$decrease\n" if($point >= 5000 && $point <= 5030);
	print OUTFALL "$decrease\n" if($point >= 10000 && $point <= 10100);
	print OUTFALL "$decrease\n" if($point >= 20000);
}
system("sleep 2s");
system("mv F_all_mismatch_small.txt F_all_mismatch.txt");


open INPUT, "<R_all_mismatch.txt" or die "cant open file!\n";#circos; the maximum plot is currently set at [25000].
open INPUT2, "<R_all_mismatch.txt" or die "cant open file!\n";
open OUTRALL, ">R_all_mismatch_small.txt" or die "cant save file!\n";
my $decrease2 = <INPUT2>;
my $point = 0;
while (my $decrease = <INPUT>){
	$decrease2 = <INPUT2>;
	chomp($decrease);
	chomp($decrease2);
	my @decreasebox = split(/\t/, $decrease);
	my @decreasebox2 = split(/\t/, $decrease2);
	$point++ if($decreasebox[3] == 0);
	unless($decreasebox[3] == 0){
		$point = 0;
		my $half = ($decreasebox[2] - $decreasebox[1]) / 2;
		my $halfA = $decreasebox[1] + $half;
		$halfA = int($halfA);
		my $halfB = $halfA + 1;
		print OUTRALL "$decreasebox[0]\t$decreasebox[1]\t$halfA\t0\n";
		print OUTRALL "$decreasebox[0]\t$halfB\t$decreasebox[2]\t$decreasebox[3]\n";
		print OUTRALL "$decreasebox2[0]\t$decreasebox2[1]\t$decreasebox2[2]\t$decreasebox2[3]\n" if($decreasebox2[3] == 0);
	}
	print OUTRALL "$decrease\n" if($point >= 10 && $point <= 15);
	print OUTRALL "$decrease\n" if($point >= 100 && $point <= 130);
	print OUTRALL "$decrease\n" if($point >= 400 && $point <= 430);
	print OUTRALL "$decrease\n" if($point >= 1500 && $point <= 1530);
	print OUTRALL "$decrease\n" if($point >= 2000 && $point <= 2030);
	print OUTRALL "$decrease\n" if($point >= 2500 && $point <= 1530);
	print OUTRALL "$decrease\n" if($point >= 3000 && $point <= 3030);
	print OUTRALL "$decrease\n" if($point >= 4000 && $point <= 4030);
	print OUTRALL "$decrease\n" if($point >= 5000 && $point <= 5030);
	print OUTRALL "$decrease\n" if($point >= 10000 && $point <= 10100);
	print OUTRALL "$decrease\n" if($point >= 20000);
}
system("sleep 2s");
system("mv R_all_mismatch_small.txt R_all_mismatch.txt");

my $radius2 = "$radius" . "r";
my $radius3 = $radius - 0.5;
my $radius3 = "$radius3" . "r";

my $labelsize = 15;
$labelsize = 25 if($k <= 10);
$labelsize = 20 if($k <= 50 && $k >>10);
open OUT5G, ">MAIN2.CONF" or die "cant save file!\n";
my $textG = <<"EOS";
# 1.2 IDEOGRAM LABELS, TICKS, AND MODULARIZING CONFIGURATION
karyotype = chr_all.txt
chromosomes_units = 100000;
<<include ideogram2.conf>>
<<include ticks2.conf>>

<plots>

<plot>
type      = line
thickness = 0.1
max_gap = 1u
file    = FR_all_coverage.txt
color   = while
min     = 0
max     = $top2
r0      = 0.90r
r1      = 0.98r
fill_color = 76,76,76
<backgrounds>
<background>
color     = white
y0        = 0.006
</background>
<background>
color     = white
y1        = 0.002
</background>
</backgrounds>
<axes>
<axis>
color     = white
thickness = 1
spacing   = 10r
</axis>
</axes>
</plots>

<plot>
type      = line
thickness = 0.1
max_gap = 1u
file    = F_all_mismatch.txt
color   = while
min     = 0
max     = 20
r0      = 0.83r
r1      = 0.90r
fill_color = blue
<backgrounds>
<background>
color     = white
y0        = 0.006
</background>
<background>
color     = white
y1        = 0.002
</background>
</backgrounds>
<axes>
<axis>
color     = white
thickness = 1
spacing   = 10r
</axis>
</axes>
</plots>

<plot>
type      = line
thickness = 0.1
max_gap = 1u
file    = R_all_mismatch.txt
color   = while
min     = 0
max     = 20
r0      = 0.76r
r1      = 0.83r
fill_color = red
<backgrounds>
<background>
color     = white
y0        = 0.006
</background>
<background>
color     = white
y1        = 0.002
</background>
</backgrounds>
<axes>
<axis>
color     = white
thickness = 1
spacing   = 10r
</axis>
</axes>
</plots>


<plot>
type         = scatter
file    = In_all.txt
glyph        = triangle
glyph_size   = 20
orientation  = out
fill_color   = black
r0      = 0.72r
r1      = 0.72r
</plots>


<plot>
type      = heatmap
file    = Del_all.txt
scale_log_base = 2
min          = 0.00001
color        = red
r0      = 0.69r
r1      = 0.71r
</plots>

</plots>

<image>
<<include etc/image.conf>>
</image>
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
EOS
print OUT5G "$textG\n";

open OUT5A, ">IDEOGRAM2.CONF" or die "cant save file!\n";
my $textA = <<"EOS";
<ideogram>
<spacing>
break   = 0.1r
default = 0.1u
</spacing>
<<include ideogram.position.conf>>
<<include ideogram.label.conf>>
<<include bands.conf>>
radius*       = $radius2
</ideogram>
EOS
print OUT5A "$textA\n";

open OUT5B, ">IDEOGRAM.LABEL.CONF" or die "cant save file!\n";
my $textB = <<"EOS";
show_label       = yes
label_font       = bold
label_radius = dims(ideogram,radius_outer) + 0.05r
label_size       = $labelsize
label_parallel   = no
label_case       = upper
label_format     = eval(sprintf("%s",var(label)))
EOS
print OUT5B "$textB\n";


open OUT5D, ">BANDS.CONF" or die "cant save file!\n";
my $textD = <<"EOS";
	show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 0
EOS
print OUT5D "$textD\n";

open OUT5E, ">IDEOGRAM.POSITION.CONF" or die "cant save file!\n";
my $textE = <<"EOS";
radius           = $radius3
thickness        = 10p #環の太さパラメータ
fill             = yes
fill_color       = black
stroke_thickness = 10
stroke_color     = black
EOS
print OUT5E "$textE\n";

print "last circos process\n";
system("sleep 3s");
system("circos -conf MAIN2.CONF");
system("mv temp/*mismatch_mountain_all.txt circos.svg InDel.txt $outputdir");
system("mv *CONF *average_mismatch_value.txt FR_average_coverage.txt genome.txt In_all.txt Del_all.txt FR_all_coverage_small.txt log_coverage.txt chr_all.txt FR_all_coverage.txt *all_mismatch.txt temp/circos_config/");

}
####################################################################################################################################