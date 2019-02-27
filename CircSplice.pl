#!/usr/bin/perl
# usage:
# perl circ.pl chimeric.out.sam genome bed_refFlat

use strict;
use warnings;

my $chimeric_file = shift;
#my $aligned_file = shift;
my $genomefile = shift;
my $gencodefile = shift;

my $prefile = $chimeric_file . '.tmp';
preprocessing($chimeric_file, $prefile);

my $circtmpfile = $prefile . '.circ';
my $astmpfile = $prefile . '.as';
tmpsplit($prefile, $circtmpfile, $astmpfile);

my $precircfile = $chimeric_file . '.circ.pre';
my $preasfile = $chimeric_file . '.as.pre';
presplit($circtmpfile, $astmpfile, $precircfile, $preasfile);

my $precircposfile = $precircfile . '.pos';
position($precircfile, $precircposfile);
my $preasposfile = $preasfile . '.pos';
position($preasfile, $preasposfile);

my $asfile = $chimeric_file . '.as';
my $circfile = $chimeric_file . '.circ';
as($precircposfile, $preasposfile, $asfile, $circfile);

my $circbedfile = $circfile . '.bed';
tran2bed($circfile, $circbedfile);

my $circsitefile = $circfile.'.site';
getsplice($circbedfile, $genomefile, $circsitefile);

my $comcircsitefile = $circfile.'.com';
circ_gtag_ctac($circsitefile, $comcircsitefile);

my $forasfile = $asfile . '.find';
findcircas($comcircsitefile, $asfile, $forasfile);
my $forcircfile = $circfile . '.find';
findcircas($comcircsitefile, $circfile, $forcircfile);

### for CIRC ###

my $circreadsfile = $forcircfile . '.bed';
getcircreads($forcircfile, $circreadsfile);

my $circintersectfile = $circreadsfile . '.intersect';
intersect($circreadsfile, $gencodefile, $circintersectfile);

my $resultcircfile = $chimeric_file . '.result.circ';
getcircresult($circintersectfile, $chimeric_file, $resultcircfile);

### for AS ###

my $asreadsbedfile = $forasfile . '.bed';
astran2bed($forasfile, $asreadsbedfile);

my $asintersectfile = $asreadsbedfile . '.intersect';
intersect($asreadsbedfile, $gencodefile, $asintersectfile);

my $readastypefile = $asfile . '.type';
astype($asintersectfile, $readastypefile);

my $astypefastafile = $readastypefile . '.fasta';
tran2fasta($forasfile, $astypefastafile);

my $asseqfile = $asfile . '.seq';
getasseq($astypefastafile, $genomefile, $asseqfile);

my $comassitefile = $asfile.'.com';
as_gtag_ctac($asseqfile, $readastypefile, $comassitefile);

my $asreadsinfo = $chimeric_file . '.asreads';
getasreads($comassitefile, $asreadsinfo);

my $resultasfile = $chimeric_file . '.result.as';
getasresult($asreadsinfo, $chimeric_file, $resultasfile);

deltmpfiles();

sub deltmpfiles {
		system("rm " . $prefile);
		system("rm " . $circtmpfile);
		system("rm " . $astmpfile);
		system("rm " . $precircfile);
		system("rm " . $preasfile);
		system("rm " . $precircposfile);
		system("rm " . $preasposfile);
		system("rm " . $asfile);
		system("rm " . $circfile);
		system("rm " . $circbedfile);
		system("rm " . $circsitefile);
		system("rm " . $comcircsitefile);
		system("rm " . $forasfile);
		system("rm " . $forcircfile);
		system("rm " . $circreadsfile);
		system("rm " . $circintersectfile);
		system("rm " . $asreadsbedfile);
		system("rm " . $asintersectfile);
		system("rm " . $readastypefile);
		system("rm " . $astypefastafile);
		system("rm " . $asseqfile);
		system("rm " . $comassitefile);
		system("rm " . $asreadsinfo);
}


sub preprocessing {
	my ($chimeric_file, $prefile) = @_;
	
	my %id_reads;
	open(FILE, $chimeric_file) or die;
	while(<FILE>){
		s/[\r\n]//g;
		if(/^\@/){
		}else{
			my @row = split /\t/;	
			my $strand = '';
			if($row[1] & 16) {
				$strand =  "-";
			} else {
				$strand =  "+";
			}	
			push @{$id_reads{$row[0]}}, $strand . "\t" . $_;
		}
	}
	close(FILE);
	
	open(RESULT, '>'.$prefile) or die;
	foreach (sort {$a cmp $b} keys %id_reads) {
		my @num = @{$id_reads{$_}};
		if(@num <= 2){
			#for(my $i=0;$i<@num;$i++){
				#print RESULT2 "$num[$i]\n";
			#}
		}elsif(@num > 2){
			for(my $i=0; $i<@num; $i++){
				print RESULT "$num[$i]\n";
			}
		}
	}
	close(RESULT);
}


sub tmpsplit {
	my ($prefile, $circfile, $asfile) = @_;
	
	open(RESULT, '>'.$circfile) or die;
	open(RESULT2, '>'.$asfile) or die;
	
	my %id_reads;
	open(FILE, $prefile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row=split /\t/;
		push @{$id_reads{$row[1]}},$_;
	}
	close(FILE);
	
	foreach my $id (sort {$a cmp $b} keys %id_reads) {
		my @num=@{$id_reads{$id}};
		my %seq_reads;
		for(my $i=0; $i<@num; $i++){
			my @row=split /\t/, $num[$i];
	
			###反转互补序列
			my @lxm = split //, $row[10];
			@lxm = reverse @lxm;
			my $seq = join '', @lxm;
			$seq =~ tr/AGCT/TCGA/;
	
			push @{$seq_reads{$row[10]}}, $num[$i];
			push @{$seq_reads{$seq}}, $num[$i];
		}
		
		my %circ_num;
		my %splice_num;
		foreach my $seq (keys %seq_reads) {   ###序列为key
			my @tmp = @{$seq_reads{$seq}};
			#my %kao;
			#@tmp = grep { ++$kao{$_} < 2 } @tmp;  ###数组去重复
			if(@tmp == 1){
				$splice_num{$tmp[0]} = 1;
			}else{
				my $value = join "\n", @tmp;
				$circ_num{$value} = 1;
			}
		}
	
		foreach my $circread (keys %circ_num) {
			print RESULT "$circread\n";
		}
		foreach my $spliceread (keys %splice_num) {
			print RESULT2 "$spliceread\n";
		}
	}	
	close(RESULT);
	close(RESULT2);
}


sub presplit {
	my ($circtmpfile, $astmpfile, $precircfile, $preasfile) = @_;
	
	my %strand;
	my %chr;
	my %id_reads;
	open(FILE, $circtmpfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		push @{$strand{$row[1]}}, $row[0];
		push @{$chr{$row[1]}}, $row[3];
		push @{$id_reads{$row[1]}}, $_;
	}
	close(FILE);
	
	
	open(RESULT,'>'.$precircfile) or die;
	my %id_sign;
	foreach my $id (keys %id_reads) {
		my @strand = @{$strand{$id}};
		my %kao;
		@strand = grep { ++$kao{$_} < 2 } @strand;  ###数组去重复
		my @chr = @{$chr{$id}};
		my %kao1;
		@chr = grep { ++$kao1{$_} < 2 } @chr;  ###数组去重复
		if(@strand == 1 and @chr == 1){ ###考虑链和染色体
		#if(@chr==1){  ###只考虑染色体
			my @row = @{$id_reads{$id}};
			print RESULT join "\n",@row;
			print RESULT "\n";
			$id_sign{$id}=1;
		}
	}
	close(RESULT);
	
	open(RESULT,'>'.$preasfile) or die;
	open(FILE, $astmpfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		if(exists $id_sign{$row[1]}){
			print RESULT "$_\n";
		}
	}
	close(FILE);
	close(RESULT);
}


sub position {
	my ($file, $posfile) = @_;
	
	open(RESULT, '>'.$posfile) or die;
	my %pos_sign;
	open(FILE, $file) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		my @a = $row[6] =~ /(\d+)[A-Z]/g;
		my @b = $row[6] =~ /\d+([A-Z])/g;
		my $tmp = join '', @b;
		#if($tmp=~/DN/g or $tmp=~/ND/g){
			$pos_sign{$tmp}=1;
		#	print RESULT2 "$tmp\n";
		#}
	
		my @loc;
		my $co = $row[4];
		$co -= 1;  ### $co��1-based,�ĳ�0-based
		push @loc, $co;
		for(my $i=0; $i<@a; $i++){
			if($b[$i] eq 'M' or $b[$i] eq 'N' or $b[$i] eq 'D'){
				$co += $a[$i];
				push @loc, $co;
			}elsif($b[$i] eq 'I'){
				push @loc, $co;
			}
		}
	
		my %loc_sign;
		for(my $k=0; $k<@loc; $k++){
			push @{$loc_sign{$loc[$k]}},1;
		}
		my @new_loc;
		foreach my $loc (keys %loc_sign) {
			my @dd = @{$loc_sign{$loc}};
			if(@dd == 1){
				push @new_loc, $loc;
			}
		}
	
		@new_loc = sort {$a <=> $b} @new_loc;
		print RESULT join ",", @new_loc;
		print RESULT "\t$_\n";
	}
	close(FILE);
	close(RESULT);
}


sub as {
	my ($circposfile, $spliceposfile, $asfile, $ascircfile) = @_;
	
	my %id_poses;
	my %id_reads;
	open(FILE, $circposfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		push @{$id_poses{$row[2]}}, $row[0];
		push @{$id_reads{$row[2]}}, $_;
	}
	close(FILE);
	
	my %as_sign;
	my %circas_sign;
	open(FILE, $spliceposfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		if(exists $id_poses{$row[2]}){
			my @poses = @{$id_poses{$row[2]}};
			my $loc = join ',',@poses;
			my @circ = split /,/, $loc;
			@circ = sort {$a<=>$b} @circ;
			my @lin = split /,/, $row[0];
			my $start = $circ[0] - 2;
			my $end = $circ[-1] + 2;
			if($start <= $lin[0] and $lin[-1] <= $end){  ###保证另一端在环里面
				my @reads = @{$id_reads{$row[2]}};
				for(my $i=0; $i<@reads; $i++){
					$as_sign{$reads[$i]} = 1;			###储存环端
					$circas_sign{$reads[$i]} = 1;	###储存环端
				}
				$as_sign{$_} = 1;								###储存线性端
			}
		}
	}
	close(FILE);
	
	open(RESULT,'>'.$asfile) or die;
	foreach  (sort {$a cmp $b} keys %as_sign) {
		print RESULT "$_\n";
	}
	close(RESULT);
	
	open(RESULT,'>'.$ascircfile) or die;
	foreach  (sort {$a cmp $b} keys %circas_sign) {
		print RESULT "$_\n";
	}
	close(RESULT);
}


sub tran2bed {
	my ($ascircfile, $ascircbedfile) = @_;
	
	my %id_reads;
	open(FILE, $ascircfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		push @{$id_reads{$row[2]}}, $_;
	}
	close(FILE);
	
	open(RESULT,'>'.$ascircbedfile) or die;
	foreach  (sort {$a cmp $b} keys %id_reads) {
		my @reads = @{$id_reads{$_}};
		my @loc;
		my $chr;
		my $strand;
		for(my $i=0; $i<@reads; $i++){
			my @row = split /\t/, $reads[$i];
			push @loc, $row[0];
			$chr = $row[4];
			$strand = $row[1];
		}
		if($chr =~ /chr/){
			my $loc = join ',', @loc;
			my @loc2 = split /,/, $loc;
			@loc2 = sort {$a <=> $b} @loc2;
			my $start = $loc2[0] - 2;
			my $end = $loc2[-1] + 2;
			print RESULT "$chr\t$start\t$loc2[0]\t$_\t$loc2[0]-$loc2[-1]\t$strand\n";
			print RESULT "$chr\t$loc2[-1]\t$end\t$_\t$loc2[0]-$loc2[-1]\t$strand\n";
		}
	}
	close(RESULT);
}


sub getsplice {
	my ($ascircbedfile, $genomefile, $circsitefile) = @_;
	#my $bedversion = `bedtools --version`;
	#if($bedversion =~ /2\.25/) {
	#	my $cmd = "bedtools getfasta -fi " . $genomefile . " -tab -name+ -bed " . $ascircbedfile . " > " . $circsitefile;
	#	system($cmd);
	#} elsif($bedversion =~ /2\.26/) {
	#	my $cmd = "bedtools getfasta -fi " . $genomefile . " -tab -name+ -bed " . $ascircbedfile . " -fo " . $circsitefile;
	#	system($cmd);
	#} else {
	#	print("bedtools v2.25 or v2.26 is needed for this program.\n");
	#	exit(1);
	#}
	my $cmd = "bedtools getfasta -fi " . $genomefile . " -tab -name+ -bed " . $ascircbedfile . " -fo " . $circsitefile;
	system($cmd);
}


sub circ_gtag_ctac {
	my ($circsplicfile, $combinefile) = @_;
	
	my $ref={};
	open(FILE, $circsplicfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		my @read = split /::/, $row[0];
		$row[1] = uc $row[1];
		$read[1] =~ /(\S+):(\d+)-(\d+)/;
		$ref->{$read[0]}->{$2}->{$row[1]} = $.; ###整行装在了ref的key里面
	}
	close(FILE);
	
	open(RESULT,'>'.$combinefile) or die;
	foreach my $cor (sort {$a cmp $b} keys %$ref) {
		print RESULT "$cor\t";
		foreach my $loc (sort {$b <=> $a} keys %{$ref->{$cor}}) {
			#print RESULT "\t$loc";
			foreach my $site (keys %{$ref->{$cor}->{$loc}}) {
				print RESULT "$site";
			}
		}
		print RESULT "\n";
	}
	close(RESULT);
}


sub findcircas {
	my ($combinefile, $targetfile, $resultfile) = @_;
	
	my %hash;
	open(FILE, $combinefile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		$hash{$row[0]} = $row[1];
	}
	close(FILE);
	
	open(RESULT,'>'.$resultfile) or die;
	open(FILE, $targetfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		if(exists $hash{$row[2]}){
			if($hash{$row[2]} eq 'GTAG' or $hash{$row[2]} eq 'CTAC'){
				print RESULT "$_\n";
			}
		}
	}
	close(FILE);
	close(RESULT);
}


sub getcircreads {
	my ($findcircfile, $circreadsfile) = @_;
	
	my %id_reads;
	open(FILE, $findcircfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		push @{$id_reads{$row[2]}},$_;
	}
	close(FILE);

	my %read_id;
	my %read_strand;
	my %read_chr;
	foreach  (sort {$a cmp $b} keys %id_reads) {
		my @reads = @{$id_reads{$_}};
		my @loc;
		my $chr;
		my $strand;
		for(my $i=0; $i<@reads; $i++){
			my @row = split /\t/, $reads[$i];
			push @loc, $row[0];
			$chr = $row[4];
			$strand = $row[1];
		}
		my $loc = join ',',@loc;
		my @loc2 = split /,/,$loc;
		@loc2 = sort {$a <=> $b} @loc2;
		# my $info = "$loc2[0],$loc2[-1]\t$strand\t$_\t$chr";
		my $locstr = $loc2[0].",".$loc2[-1];
		push @{$read_id{$locstr}}, $_;  					# read id
		push @{$read_strand{$locstr}}, $strand;		# read strand
		push @{$read_chr{$locstr}}, $chr;   			# read chr
	}
	open(RESULT,'>'.$circreadsfile) or die;
	foreach my $key (sort {$a cmp $b} keys %read_id) {
			my @num = @{$read_id{$key}};
			my @strand = @{$read_strand{$key}};
			my @chr = @{$read_chr{$key}};
			my %kao;
			@num = grep { ++$kao{$_} < 2 } @num;  ###数组去重复
			my %kao2;
			@strand = grep { ++$kao2{$_} < 2 } @strand;  ###数组去重复
			my %kao4;
			@chr = grep { ++$kao4{$_} < 2 } @chr;  ###数组去重复
			my $n = @num;
			my @tmp = split /,/, $key;
			print RESULT join ';',@chr;
			print RESULT "\t";
			print RESULT "$tmp[0]\t$tmp[1]\t";
			print RESULT "$key\t$n\t";
			print RESULT join ';', @strand;
			print RESULT "\t";
			print RESULT join ';', @chr;
			print RESULT "\t";
			print RESULT join ';', @num;
			print RESULT "\n";
		}
	close(RESULT);
}


sub astran2bed {
	my ($findasfile, $asreadsfile) = @_;
	
	open(RESULT, '>'.$asreadsfile) or die;
	open(FILE, $findasfile) or die;	
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		my @key = split /,/, $row[0];
		print RESULT "$row[4]\t$key[0]\t$key[-1]\t$row[2]\t$row[3]\t$row[1]\t$_\n";
	}
	
	close(FILE);
	close(RESULT); 
}


sub intersect {
	my ($readsfile, $gencodefile, $intersectfile) = @_;
	# bedtools intersect -wo -f 1 -a 11-circ.txt -b bed-gencode_all_hg38.txt > 12-circ.txt
	system("bedtools intersect -wo -f 1 -a " . $readsfile . " -b " . $gencodefile . " > " . $intersectfile);
}


sub getcircresult {
	my ($intersectfile, $chimeric_file, $resultcircfile) = @_;

	my $run = "samtools flagstat ". $chimeric_file ." | grep 'properly paired' | awk '{print \$1}' > awk.tmp";
	system($run);
	my $total_count = 0;
	open(FD, "awk.tmp") or die "$!";
	while(<FD>){
		$total_count = $_;
	}
	close(FD);
	$run = "rm awk.tmp";
	system($run);

	my %hash;
	open(FILE, $intersectfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row = split /\t/;
		push @{$hash{$row[0]}}, $_;  # 还是用chr做index
	}
	close(FILE);
	
	my %hash2;
	my %hash3;
	open(FILE,$intersectfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		if(exists $hash3{$_}){  ###每行只用一次
		}else{
			my @row = split /\t/;
			my @loc = split /,/,$row[3];
			if(exists $hash{$row[0]}){
				my @num2;
				my @num = @{$hash{$row[0]}};
				for(my $i=0; $i<@num; $i++){
					if(exists $hash3{$num[$i]}){  ###每行只用一次
					}else{
						my @line = split /\t/, $num[$i];
						my @loc2 = split /,/, $line[3];
						if(abs($loc[0]-$loc2[0]) <= 2 and abs($loc[1]-$loc2[1]) <= 2){  #设定容错范围
							push @num2, $num[$i];
							$hash3{$num[$i]} = 1;   ###每行只用一次
						}
					}
				}
				@num2 = sort {$a cmp $b} @num2;
				my $tmp = join "\t\t\t", @num2;
				$hash2{$tmp} = 1;
			}
		}
	}
	close(FILE);
	
	open(RESULT,'>'.$resultcircfile) or die;
	foreach my $key (sort {$a cmp $b} keys %hash2) {
		my @row=split /\t\t\t/, $key;
		my @as;
		my @id;
		my @gene;
		my @strand;
		my @stranda;
		my @chr;
		#my @read;
		#my @exon;
		for(my $i=0; $i<@row; $i++){
			my @line = split /\t/, $row[$i];
			push @as, $line[3];
			push @id, $line[7];
			my $shi = $line[14].','.$line[15].','.$line[25];
			push @gene, $shi;
			push @strand, $line[5];
			push @stranda, $line[13];
			push @chr, $line[0];
			#push @read,$line[25];
			#push @exon,$line[8];
		}
		my %up;
		@as = grep {++$up{$_}<2} @as;
		my $as = join '|', @as;	
		################合并id
		my $id = join ';', @id;
		my @id2 = split /\;/, $id;
		my %kao;
		@id2 = grep { ++$kao{$_} < 2 } @id2;  ###数组去重复
		my $id2 = join ';', @id2;
		my $n = @id2;  ###统计read个数
		################合并gene
		my $gene = join ';', @gene;
		my @gene2 = split /\;/, $gene;
		my %kao1;
		@gene2 = grep { ++$kao1{$_} < 2 } @gene2;  ###数组去重复
		my $gene2 = join ';', @gene2;
		################合并strand
		my $strand = join ';', @strand;
		my @strand2 = split /\;/, $strand;
		my %kao2;
		@strand2 = grep { ++$kao2{$_} < 2 } @strand2;  ###数组去重复
		my $strand2 = join ';', @strand2;
		################合并stranda
		my $stranda = join ';', @stranda;
		my @stranda2 = split /\;/, $stranda;
		my %kao6;
		@stranda2 = grep { ++$kao6{$_} < 2 } @stranda2;  ###数组去重复
		my $stranda2 = join ';', @stranda2;
		################合并chr
		my $chr = join ';', @chr;
		my @chr2 = split /\;/, $chr;
		my %kao3;
		@chr2 = grep { ++$kao3{$_} < 2 } @chr2;  ###数组去重复
		my $chr2 = join ';', @chr2;
		################合并read
		#my $read=join ';',@read;
		#my @read2=split /\;/,$read;
		#my %kao4;
		#@read2=grep { ++$kao4{$_} < 2 } @read2;  ###数组去重复
		#my $read2=join ';',@read2;
		################合并exon
		#my $exon=join ';',@exon;
		#my @exon2=split /\;/,$exon;
		#my %kao5;
		#@exon2=grep { ++$kao5{$_} < 2 } @exon2;  ###数组去重复
		#my $exon2=join ';',@exon2;
		################	
		my $norm=$n/$total_count*1000000;	
		print RESULT "$as\t$n\t$norm\t$id2\t$gene2\t$stranda2\t$chr2\n";
		#print RESULT "$as\t$n\t$id2\t$gene2\t$strand2\t$stranda2\t$chr2\n";
	}
	close(RESULT);
}


sub astype {
	my ($asintersectfile, $readastypefile) = @_;
	
	open(RESULT,'>'.$readastypefile) or die;
	open(FILE, $asintersectfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row=split /\t/;
		my $info = join "\t", @row[6..22];
		$info .= "\t";
		$info .= join "\t",@row[29..40];
		
		@row = split /\t/, $info;
	
		###构建read的坐标
		my @read = split /,/,$row[0];
		@read = sort {$a<=>$b} @read;
		
		###构建注释基因的外显子坐标
		my $loc = $row[26].''.$row[27];
		$loc =~ s/,$//;
		my @loc = split /,/, $loc;
		@loc = sort {$a<=>$b} @loc;
	
		###开始判断read中是否有坐标跟loc中至少有一个一致，容错2bp
		my %hash;
		for(my $i=0; $i<@read; $i++){
			my $value1 = $read[$i] - 2;
			my $value2 = $read[$i] - 1;
			my $value3 = $read[$i];
			my $value4 = $read[$i] + 1;
			my $value5 = $read[$i] + 2;
			$hash{$value1} = 1;
			$hash{$value2} = 1;
			$hash{$value3} = 1;
			$hash{$value4} = 1;
			$hash{$value5} = 1;
		}
		my $n = 0;
		for(my $i=0; $i<@loc; $i++){
			if(exists $hash{$loc[$i]}){
				$n++;
			}
		}
		if($n>0){
			###开始析出剪切类型SE A3SS A5SS
			my @type;
			for(my $i=1; $i<@read-1; $i=$i+2){
				my $start = $read[$i];
				my $end = $read[$i+1];
				for(my $k=0; $k<@loc; $k=$k+2){
					if(($start <= $loc[$k] and $loc[$k+1] < $end) or ($start < $loc[$k] and $loc[$k+1] <= $end)){
						if(abs($loc[$k]-$loc[$k+1])>2){  ###剪切事件坐标不能太小
							push @type,"SE,$loc[$k],$loc[$k+1]\|$start,$end";
						}
					}elsif($loc[$k] <= $start and $end <= $loc[$k+1]){
						if(abs($start-$end) > 2){
							push @type,"SE,$start,$end\|$start,$end";
						}
					}elsif($loc[$k] < $start and $start < $loc[$k+1] and $loc[$k+1] < $end){
						if(abs($start-$loc[$k+1])>2){
							if($row[20] eq '+'){
								push @type,"A5SS,$start,$loc[$k+1]\|$start,$end";
							}elsif($row[20] eq '-'){
								push @type,"A3SS,$start,$loc[$k+1]\|$start,$end";
							}elsif($row[20] eq '-;+' or $row[20] eq '+;-'){
								push @type,"A5SS,$start,$loc[$k+1]\|$start,$end";
							}													
						}
					}elsif($start < $loc[$k] and $loc[$k] < $end and $end < $loc[$k+1]){
						if(abs($loc[$k]-$end)>2){
							if($row[20] eq '+'){
								push @type,"A3SS,$loc[$k],$end\|$start,$end";
							}elsif($row[20] eq '-'){
								push @type,"A5SS,$loc[$k],$end\|$start,$end";
							}elsif($row[20] eq '-;+' or $row[20] eq '+;-'){
								push @type,"A3SS,$start,$loc[$k+1]\|$start,$end";
							}							
						}
					}
				}
			}
			if(@type > 0){
				#if($row[0] eq $row[21]){   ###链可以不考虑
					print RESULT join ';',@type;
					#####去重复
					#my %saw;
					#@saw{ @site } = ( );
					#my @saw = sort keys %saw;
					#print RESULT "\t";
					#print RESULT join ';',@saw;
					#####
					#my $o=@saw;
					#print RESULT "\t$o";
					print RESULT "\t$info\n";
				#}
			}
			###开始析出剪切类型RI
			my @type1;
			for(my $i=0;$i<@read;$i=$i+2){
				my $start=$read[$i];
				my $end=$read[$i+1];
				for(my $k=1;$k<@loc-1;$k=$k+2){
					if($start<=$loc[$k] and $loc[$k+1]<=$end){
					#if((($end-$loc[$k])>=4 and $start<=$loc[$k]) or (($loc[$k+1]-$start)>=4 and $end>=$loc[$k+1]) or ($loc[$k]<=$start and $end<=$loc[$k+1]) or ($start<=$loc[$k] and $loc[$k+1]<=$end)){
						if(abs($loc[$k]-$loc[$k+1])>2){  ###剪切事件坐标不能太小
							push @type1,"RI,$loc[$k],$loc[$k+1]";
						}
					}
				}
			}
			if(@type1>0){
				#if($row[0] eq $row[21]){   ###链可以不考虑
					print RESULT join ';',@type1;
					#print RESULT "\tNA\tNA";
					print RESULT "\t$info\n";
				#}
			}
		}
	}
	close(FILE);
}


sub tran2fasta {
	my ($readastypefile, $astypefastafile) = @_;
	
	my %hash;
	open(FILE, $readastypefile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row=split /\t/;
		my @read=split /,/, $row[0];
		if($row[4] =~ /chr/){
			for(my $i=1; $i<@read-1; $i=$i+2){
				my $start=$read[$i]+2;
				my $end=$read[$i+1]-2;
				my $key="$row[4]\t$read[$i]\t$start\t$read[$i]-$read[$i+1]";
				$hash{$key}=1;
				my $key1="$row[4]\t$end\t$read[$i+1]\t$read[$i]-$read[$i+1]";
				$hash{$key1}=1;
			}
		}
	}
	close(FILE);
	
	open(RESULT, '>'.$astypefastafile) or die;
	foreach  (sort {$a cmp $b} keys %hash) {
		print RESULT "$_\n";
	}
	close(RESULT);
}


sub getasseq {
	my ($astypefastafile, $genomefile, $asseqfile) = @_;
	
	my $cmd = "bedtools getfasta -fi " . $genomefile . " -tab -name+ -bed " . $astypefastafile . " -fo " . $asseqfile;
	system($cmd);
}


sub as_gtag_ctac {
	my ($asseqfile, $readastypefile, $comassitefile) = @_;
	
	my $ref={};
	open(FILE, $asseqfile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row=split /\t/;
		my @read=split /::/,$row[0];
		$row[1]=uc $row[1];
		$read[1]=~/(\S+):(\d+)-(\d+)/;
		my $id=$1.':'.$read[0];
	  $ref->{$id}->{$2}->{$row[1]} = $.; ###整行装在了ref的key里面
	}
	close(FILE);
	
	my %hash;
	foreach my $cor (sort {$a cmp $b} keys %$ref) {
		my $info = "$cor\t";
		foreach my $loc (sort {$a <=> $b} keys %{$ref->{$cor}}) {
			#print RESULT "\t$loc";
			foreach my $site (keys %{$ref->{$cor}->{$loc}}) {
				$info .= "$site";
			}
		}
		my @row=split /\t/, $info;
		$hash{$row[0]}=$row[1];
	}
	
	open(RESULT, '>'.$comassitefile) or die;
	open(FILE, $readastypefile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row=split /\t/;
		if($row[0]=~/\|/){
			my @num2;
			my @loc=split /\;/,$row[0];
			for(my $i=0;$i<@loc;$i++){
				my @tmp=split /\|/,$loc[$i];
				$tmp[1]=~s/,/-/;
				my $id=$row[5].':'.$tmp[1];
				if(exists $hash{$id}){
					if($hash{$id} eq 'GTAG' or $hash{$id} eq 'CTAC'){
						push @num2,$tmp[0];
					}
				}
			}
			if(@num2){  ###如果有满足条件的剪切类型
				print RESULT join ';',@num2;
				my $tmp=shift @row;
				print RESULT "\t";
				print RESULT join "\t",@row;
				print RESULT "\n";
			}
		}else{
			print RESULT "$_\n";
		}
	}
	close(FILE);
	close(RESULT);
}


sub getasreads {
	my ($comassitefile, $asreadsinfo) = @_;
	
	my %hash;
	my %hash2;
	my %strand1;
	my %strand2;
	my %chr;
	my %loc;
	my %exon;
	#my %junc;
	open(RESULT, '>'.$asreadsinfo) or die;
	open(FILE ,$comassitefile) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row=split /\t/;
		my @tmp=split /\;/,$row[0];
		my $gene=$row[18].','.$row[19].','.$row[29];
		###把外显子start和end放一起
		$row[27]=~s/,$//;
		$row[28]=~s/,$//;
		my $exon=$row[27].','.$row[28];
		my @exon=split /,/,$exon;
		@exon=sort {$a<=>$b} @exon;
		$exon=join ',',@exon;
		###
		#if($row[2] ne $row[21]){     ###要求链不一致，可以先不考虑
			for(my $i=0;$i<@tmp;$i++){
				push @{$hash{$tmp[$i]}},$row[3];  ###read id
				push @{$hash2{$tmp[$i]}},$gene;   ###gene symbol and type
				push @{$strand1{$tmp[$i]}},$row[2];   ### read strand
				push @{$strand2{$tmp[$i]}},$row[21];  ### gene strand
				push @{$chr{$tmp[$i]}},$row[5];   ### chr
				push @{$loc{$tmp[$i]}},$row[1];   ###read 坐标
				push @{$exon{$tmp[$i]}},$exon;   ### 外显子坐标
				#push @{$junc{$tmp[$i]}},$row[1];  ###junction type
			}
		#}
	}
	close(FILE);
	
	foreach my $key (sort {$a cmp $b} keys %hash) {
		my @num=@{$hash{$key}};
		my @gene=@{$hash2{$key}};
		my @strand1=@{$strand1{$key}};
		my @strand2=@{$strand2{$key}};
		my @chr=@{$chr{$key}};
		my @loc=@{$loc{$key}};
		my @exon=@{$exon{$key}};
		#my @junc=@{$junc{$key}};
		my %kao;
		@num=grep { ++$kao{$_} < 2 } @num;  ###数组去重复
		my %kao1;
		@gene=grep { ++$kao1{$_} < 2 } @gene;  ###数组去重复
		my %kao2;
		@strand1=grep { ++$kao2{$_} < 2 } @strand1;  ###数组去重复
		my %kao3;
		@strand2=grep { ++$kao3{$_} < 2 } @strand2;  ###数组去重复
		my %kao4;
		@chr=grep { ++$kao4{$_} < 2 } @chr;  ###数组去重复
		my %kao5;
		@loc=grep { ++$kao5{$_} < 2 } @loc;  ###数组去重复
		my %kao6;
		@exon=grep { ++$kao6{$_} < 2 } @exon;  ###数组去重复
		#my %kao7;
		#@junc=grep { ++$kao7{$_} < 2 } @junc;  ###数组去重复
		my $n=@num;
		print RESULT "$key\t$n\t";
		print RESULT join ';',@num;
		#print RESULT "\t";
		#print RESULT join ';',@junc;
		print RESULT "\t";
		print RESULT join ';',@gene;
		print RESULT "\t";
		print RESULT join ';',@strand1;
		print RESULT "\t";
		print RESULT join ';',@strand2;
		print RESULT "\t";
		print RESULT join ';',@chr;
		print RESULT "\t";
		print RESULT join ';',@loc;
		print RESULT "\t";
		print RESULT join ';',@exon;
		print RESULT "\n";
	}
	close(RESULT);
}


sub getasresult {
	my ($asreadsinfo, $chimeric_file, $resultasfile) = @_;
	
	my $run = "samtools flagstat ". $chimeric_file ." | grep 'properly paired' | awk '{print \$1}' > awk.tmp";
	system($run);
	my $total_count = 0;
	open(FD, "awk.tmp") or die "$!";
	while(<FD>){
		$total_count = $_;
	}
	close(FD);
	$run = "rm awk.tmp";
	system($run);
	
	my %hash;
	open(FILE, $asreadsinfo) or die;
	while(<FILE>){
		s/[\r\n]//g;
		my @row=split /\t/;
		push @{$hash{$row[6]}},$_;  ###还是用chr做index
	}
	close(FILE);
	
	my %hash2;
	my %hash3;
	open(FILE, $asreadsinfo) or die;
	while(<FILE>){
		s/[\r\n]//g;
		if(exists $hash3{$_}){  ###每行只用一次
		}else{
			my @row=split /\t/;
			my @loc=split /,/,$row[0];
			if(exists $hash{$row[6]}){
				my @num2;
				my @num=@{$hash{$row[6]}};
				for(my $i=0;$i<@num;$i++){
					if(exists $hash3{$num[$i]}){  ###每行只用一次
					}else{
						my @line=split /\t/,$num[$i];
						my @loc2=split /,/,$line[0];
						if(abs($loc[1]-$loc2[1])<=2 and abs($loc[2]-$loc2[2])<=2){  ###设定容错范围
							push @num2,$num[$i];
							$hash3{$num[$i]}=1;   ###每行只用一次
						}
					}
				}
				@num2=sort {$a cmp $b} @num2;
				my $tmp=join "\t\t\t",@num2;
				$hash2{$tmp}=1;
			}
		}
	}
	close(FILE);
	
	open(RESULT, '>'.$resultasfile) or die;
	foreach my $key (sort {$a cmp $b} keys %hash2) {
		my @row=split /\t\t\t/,$key;
		my @as;
		my @id;
		my @gene;
		my @strand;
		my @stranda;
		my @chr;
		my @read;
		my @exon;
		for(my $i=0;$i<@row;$i++){
			my @line=split /\t/,$row[$i];
			push @as,$line[0];
			push @id,$line[2];
			push @gene,$line[3];
			push @strand,$line[4];
			push @stranda,$line[5];
			push @chr,$line[6];
			push @read,$line[7];
			push @exon,$line[8];
		}
		my $as=join '|',@as;
		################合并id
		my $id=join ';',@id;
		my @id2=split /\;/,$id;
		my %kao;
		@id2=grep { ++$kao{$_} < 2 } @id2;  ###数组去重复
		my $id2=join ';',@id2;
		my $n=@id2;  ###统计read个数
		################合并gene
		my $gene=join ';',@gene;
		my @gene2=split /\;/,$gene;
		my %kao1;
		@gene2=grep { ++$kao1{$_} < 2 } @gene2;  ###数组去重复
		my $gene2=join ';',@gene2;
		################合并strand
		my $strand=join ';',@strand;
		my @strand2=split /\;/,$strand;
		my %kao2;
		@strand2=grep { ++$kao2{$_} < 2 } @strand2;  ###数组去重复
		my $strand2=join ';',@strand2;
		################合并stranda
		my $stranda=join ';',@stranda;
		my @stranda2=split /\;/,$stranda;
		my %kao6;
		@stranda2=grep { ++$kao6{$_} < 2 } @stranda2;  ###数组去重复
		my $stranda2=join ';',@stranda2;
		################合并chr
		my $chr=join ';',@chr;
		my @chr2=split /\;/,$chr;
		my %kao3;
		@chr2=grep { ++$kao3{$_} < 2 } @chr2;  ###数组去重复
		my $chr2=join ';',@chr2;
		################合并read
		my $read=join ';',@read;
		my @read2=split /\;/,$read;
		my %kao4;
		@read2=grep { ++$kao4{$_} < 2 } @read2;  ###数组去重复
		my $read2=join ';',@read2;
		################合并exon
		my $exon=join ';',@exon;
		my @exon2=split /\;/,$exon;
		my %kao5;
		@exon2=grep { ++$kao5{$_} < 2 } @exon2;  ###数组去重复
		my $exon2=join ';',@exon2;
		################	
		my $norm=$n/$total_count*1000000;
		print RESULT "$as\t$n\t$norm\t$id2\t$gene2\t$stranda2\t$chr2\t$read2\t$exon2\n";
		#print RESULT "$as\t$n\t$norm\t$id2\t$gene2\t$strand2\t$stranda2\t$chr2\t$read2\t$exon2\n";
	}
	close(RESULT);
}
