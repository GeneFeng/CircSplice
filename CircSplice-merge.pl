#!/usr/bin/perl
# Usage: merge.pl dir-as dir-circ

use strict;
use warnings;

my $as_dir=shift;
my $circ_dir=shift;

###数组去重复
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

# merge files
sub file_merge
{
	my $dir = shift;
	my $mergefile = shift;
	open(RESULT,'>'.$mergefile) or die "$!\n";
	opendir(DIR,$dir) or die "$!\n";
	while(my $file=readdir(DIR)){
		if($file !~ /^\./) {
			open(FILE,$dir.'/'.$file) or die "$!\n";
			#my @name=split /\-/,$file;
			while(<FILE>){
				s/\r//g;
				s/\n//g;
				my @row=split /\t/;
				my $title=shift @row;
				my @tmp=split /\|/,$title;
				for(my $i=0;$i<@tmp;$i++){
					print RESULT "$tmp[$i]\t";
					print RESULT join "\t",@row;
					#print RESULT "\t$name[0]\n";
					print RESULT "\t$file\n";
				}
			}
			close(FILE);
		} 
	}
	close(DIR);
	close(RESULT);
}

# Deal with AS
sub as_deal
{
	my $file = shift;
	my $result = shift;
	my %hash;
	open(FILE,$file) or die "$!\n";
	while(<FILE>){
		s/[\r\n]//g;
		my @row=split /\t/;
		push @{$hash{$row[6]}},$_;  ###还是用chr做index
	}
	close(FILE);
	
	my %hash2;
	my %hash3;
	open(FILE,$file) or die "$!\n";
	while(<FILE>){
		s/\r//g;
		s/\n//g;
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
	
	open(RESULT,'>'.$result) or die "$!\n";
	foreach my $key (sort {$a cmp $b} keys %hash2) {
		my @row=split /\t\t\t/,$key;
		my @as;
		my @id;
		my @gene;
		#my @strand;
		my @stranda;
		my @chr;
		my @read;
		my @exon;
		my @sample;
		my @norm;
		for(my $i=0;$i<@row;$i++){
			my @line=split /\t/,$row[$i];
			push @as,$line[0];
			push @id,$line[3];
			push @gene,$line[4];
			#push @strand,$line[4];
			push @stranda,$line[5];
			push @chr,$line[6];
			push @read,$line[7];
			push @exon,$line[8];
			my @tmp=split '-',$line[9];
			my $count=$line[9].','.$line[1];
			my $count2=$tmp[0].','.$line[2];
			push @sample,$count;
			push @norm,$count2;
		}
		################合并坐标
		my $as=join '|',@as;
		my @as2=split /\|/,$as;
		@as2=uniq @as2;
		my $as2=join '|',@as2;
		################合并id
		my $id=join ';',@id;
		my @id2=split /\;/,$id;
		@id2=uniq @id2;  ###数组去重复
		my $id2=join ';',@id2;
		#my $n=@id2;  ###统计read个数
		################合并gene
		my $gene=join ';',@gene;
		my @gene2=split /\;/,$gene;
		@gene2=uniq @gene2;  ###数组去重复
		my $gene2=join ';',@gene2;
		################合并strand
		#my $strand=join ';',@strand;
		#my @strand2=split /\;/,$strand;
		#@strand2=uniq @strand2;  ###数组去重复
		#my $strand2=join ';',@strand2;
		################合并stranda
		my $stranda=join ';',@stranda;
		my @stranda2=split /\;/,$stranda;
		@stranda2=uniq @stranda2;  ###数组去重复
		my $stranda2=join ';',@stranda2;
		################合并chr
		my $chr=join ';',@chr;
		my @chr2=split /\;/,$chr;
		@chr2=uniq @chr2;  ###数组去重复
		my $chr2=join ';',@chr2;
		################合并read
		my $read=join ';',@read;
		my @read2=split /\;/,$read;
		@read2=uniq @read2;  ###数组去重复
		my $read2=join ';',@read2;
		################合并exon
		my $exon=join ';',@exon;
		my @exon2=split /\;/,$exon;
		@exon2=uniq @exon2;  ###数组去重复
		my $exon2=join ';',@exon2;
		################合并sample
		@sample=uniq @sample;  ###数组去重复
		#my $sample=join ';',@sample;
		#my @sample2=split /\;/,$sample;
		#@sample2=uniq @sample2;  ###数组去重复
		#my $sample2=join ';',@sample2;
		################合并norm
		my $norm=join ';',@norm;
		my @norm2=split /\;/,$norm;
		@norm2=uniq @norm2;  ###数组去重复
		my $norm2=join ';',@norm2;
		##########################################判断样本类型
		my @sample2;
		my @type;
		#if($sample2=~/T/ and $sample2=~/N/){
		#	$type='Cancer/Normal';
		#}elsif($sample2=~/N/){
		#	$type='Normal';
		#}else{
		#	$type='Cancer';
		#}
		foreach my $samp (@sample) {
			my @tmp=split '-', $samp;
			my @tmp2=split ',', $tmp[2];
			push @sample2, $tmp[0].','.$tmp2[1];
			push @type, $tmp[1];
		}
		@sample2=uniq @sample2;
		@type=uniq @type;
		my $sample2=join ';', @sample2;
		my $type=join '/', @type;
		################	
		#print RESULT "$as2\t$sample2\t$norm2\t$id2\t$gene2\t$strand2\t$stranda2\t$chr2\t$read2\t$exon2\n";
		print RESULT "$as2\t$type\t$sample2\t$norm2\t$gene2\t$stranda2\t$chr2\t$read2\t$exon2\t$id2\n";
	}
	close(RESULT);
}

# Deal with Circ
sub circ_deal
{
	my $cfile = shift;
	my $cresult = shift;
	my %chash;
	open(CFILE,$cfile) or die;
	while(<CFILE>){
		s/\r//g;
		s/\n//g;
		my @row=split /\t/;
		push @{$chash{$row[6]}},$_;  ###还是用chr做index
	}
	close(CFILE);
	
	my %chash2;
	my %chash3;
	open(CFILE,$cfile) or die;
	while(<CFILE>){
		s/\r//g;
		s/\n//g;
		if(exists $chash3{$_}){  ###每行只用一次
		}else{
			my @row=split /\t/;
			my @loc=split /,/,$row[0];
			if(exists $chash{$row[6]}){
				my @num2;
				my @num=@{$chash{$row[6]}};
				for(my $i=0;$i<@num;$i++){
					if(exists $chash3{$num[$i]}){  ###每行只用一次
					}else{
						my @line=split /\t/,$num[$i];
						my @loc2=split /,/,$line[0];
						if(abs($loc[0]-$loc2[0])<=2 and abs($loc[1]-$loc2[1])<=2){  ###设定容错范围,circ与as这一句不一样
							push @num2,$num[$i];
							$chash3{$num[$i]}=1;   ###每行只用一次
						}
					}
				}
				@num2=sort {$a cmp $b} @num2;
				my $tmp=join "\t\t\t",@num2;
				$chash2{$tmp}=1;
			}
		}
	}
	close(CFILE);
	
	open(RESULT,'>'.$cresult) or die;
	foreach my $key (sort {$a cmp $b} keys %chash2) {
		my @row=split /\t\t\t/,$key;
		my @as;
		my @id;
		my @gene;
		#my @strand;
		my @stranda;
		my @chr;
		my @sample;
		my @norm;
		for(my $i=0;$i<@row;$i++){
			my @line=split /\t/,$row[$i];
			push @as,$line[0];
			push @id,$line[3];
			push @gene,$line[4];
			#push @strand,$line[4];
			push @stranda,$line[5];
			push @chr,$line[6];
			my @tmp=split '-',$line[7];
			my $count=$line[7].','.$line[1];
			my $count2=$tmp[0].','.$line[2];
			push @sample,$count;
			push @norm,$count2;
		}
		################合并坐标
		my $as=join '|',@as;
		my @as2=split /\|/,$as;
		@as2=uniq @as2;
		my $as2=join '|',@as2;
		################合并id
		my $id=join ';',@id;
		my @id2=split /\;/,$id;
		@id2=uniq @id2;  ###数组去重复
		my $id2=join ';',@id2;
		#my $n=@id2;  ###统计read个数
		################合并gene
		my $gene=join ';',@gene;
		my @gene2=split /\;/,$gene;
		@gene2=uniq @gene2;  ###数组去重复
		my $gene2=join ';',@gene2;
		################合并strand
		#my $strand=join ';',@strand;
		#my @strand2=split /\;/,$strand;
		#@strand2=uniq @strand2;  ###数组去重复
		#my $strand2=join ';',@strand2;
		################合并stranda
		my $stranda=join ';',@stranda;
		my @stranda2=split /\;/,$stranda;
		@stranda2=uniq @stranda2;  ###数组去重复
		my $stranda2=join ';',@stranda2;
		################合并chr
		my $chr=join ';',@chr;
		my @chr2=split /\;/,$chr;
		@chr2=uniq @chr2;  ###数组去重复
		my $chr2=join ';',@chr2;
		################合并sample
		@sample=uniq @sample;  ###数组去重复
		#my $sample=join ';',@sample;
		#my @sample2=split /\;/,$sample;
		#@sample2=uniq @sample2;  ###数组去重复
		#my $sample2=join ';',@sample2;
		################合并norm
		my $norm=join ';',@norm;
		my @norm2=split /\;/,$norm;
		@norm2=uniq @norm2;  ###数组去重复
		my $norm2=join ';',@norm2;
		##########################################判断样本类型
		my @sample2;
		my @type;
		#if($sample2=~/T/ and $sample2=~/N/){
		#	$type='Cancer/Normal';
		#}elsif($sample2=~/N/){
		#	$type='Normal';
		#}else{
		#	$type='Cancer';
		#}
		foreach my $samp (@sample) {
			my @tmp=split '-', $samp;
			my @tmp2=split ',', $tmp[2];
			push @sample2, $tmp[0].','.$tmp2[1];
			push @type, $tmp[1];
		}
		@sample2=uniq @sample2;
		@type=uniq @type;
		my $sample2=join ';', @sample2;
		my $type=join '/', @type;
		################	
		#print RESULT "$as2\t$sample2\t$norm2\t$id2\t$gene2\t$strand2\t$stranda2\t$chr2\n";
		print RESULT "$as2\t$type\t$sample2\t$norm2\t$gene2\t$stranda2\t$chr2\t$id2\n";
	}
	close(RESULT);
}

# Save AS result
sub save_result
{
	my $asfile = shift;
	my $circfile = shift;
	my $asresult = shift;
	my $circresult = shift;
	
	my %hash;
	my %hash2;
	open(FILE,$circfile) or die "$!\n";  ###环文件
	while(<FILE>){
		s/\r//g;
		s/\n//g;
		my @row=split /\t/;
		my @id=split /\;/,$row[7];
		for(my $i=0;$i<@id;$i++){
			push @{$hash{$id[$i]}},$row[0];  
		}
		my $tmp=pop @row;
		my $tmp2=join "\t",@row;
		$hash2{$tmp2}=1;  ###保存环文件
	}
	close(FILE);
	
	
	open(ASRESULT,'>'.$asresult) or die "$!\n";
	open(ASFILE,$asfile) or die "$!\n";  ###AS
	while(<ASFILE>){
		s/\r//g;
		s/\n//g;
		my @row=split /\t/;
		my @id=split /\;/,$row[9];
		my @circ;
		for(my $i=0;$i<@id;$i++){
			if(exists $hash{$id[$i]}){
				my @tmp=@{$hash{$id[$i]}};
				push @circ,@tmp;
			}
		}
		@circ=uniq @circ;
		my $circ2=join ';',@circ;
		print ASRESULT "$row[0]\t$circ2\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\n";
	}
	close(ASFILE);
	close(ASRESULT);
	
	open(CIRCRESULT,'>'.$circresult) or die;
	foreach  (keys %hash2) {
		print CIRCRESULT "$_\n";
	}
	close(CIRCRESULT);
}

# Remove tmp files
sub tmp_clear
{
	foreach my $tmpfile (@_) {
		system("rm $tmpfile");	
	}
}

# main
# AS files
my $as_mergefile = "as_merge";
file_merge($as_dir, $as_mergefile);
my $as_dealfile = $as_mergefile . '_deal';
as_deal($as_mergefile, $as_dealfile);

# Circ files
my $circ_mergefile = "circ_merge";
file_merge($circ_dir, $circ_mergefile);
my $circ_dealfile = $circ_mergefile . '_deal';
circ_deal($circ_mergefile, $circ_dealfile);

# Save result
my $as_resultfile = "as_result";
my $circ_resultfile = "circ_result";
save_result($as_dealfile, $circ_dealfile, $as_resultfile, $circ_resultfile);
tmp_clear($as_mergefile, $as_dealfile, $circ_mergefile, $circ_dealfile);

