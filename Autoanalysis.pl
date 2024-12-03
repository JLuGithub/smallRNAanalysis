#!/usr/local/bin/perl -w
# Bookmarks: 0,0 0,0 0,0 0,2333; CollapsedSubs: genome_size  truncate  datapurify  normalization_total  normalization_totalandgenomesize  normalization_microRNA  normalization_microRNAandgenomesize  normalization_number  bowtie_program  totalstore_generation  length_distribution_percentage  count_reads_number_figure  complementary  export_miRNA_list_excel  lengthdistribution_show_5basebias  hit_distribution  overlappingcheck  basebiascurve_oneside  data_process. ### Version 2.0, 2024.10 Yunpeng Dai.
sub genome_size {
    open (input_data, "<$_[0]") or die "Couldn't open: $!";
    @size=();
    @truesize=();
    @segmentID=();
    $segmentnumber=0;
    my $temp="";
    my @temp=();
    $totalsize=0;
    my $i=0;
    @jointsize=();
    $size[0]=0;
    while (1) {
           chomp ($temp=<input_data>);
           if ($temp eq "") {
                next;
           } elsif ($temp=~/^>/) {
                $totalsize=$totalsize+$size[$segmentnumber];
                $truesize[$segmentnumber]=$size[$segmentnumber];
                $size[$segmentnumber]=(int($size[$segmentnumber]/100)+1)*100;
                $segmentnumber++;
                $segmentID[$segmentnumber]=$temp;
                $size[$segmentnumber]=0;
                next;
           } else {
                $size[$segmentnumber]=$size[$segmentnumber] + length ($temp);
           }
           if (eof) {
              last;
           }
    }

    $totalsize=$totalsize+$size[$segmentnumber];
    $truesize[$segmentnumber]=$size[$segmentnumber];

    $jointsize[0]=0;
    for ($i=1;$i<=$segmentnumber;$i++) {
          $jointsize[$i]=$jointsize[$i-1]+$size[$i];
    }

    close (input_data);
    for ($i=1;$i<=$segmentnumber;$i++) {
          print ("Segment $i : $segmentID[$i] $truesize[$i] nt\n");
    }
    return $jointsize[$segmentnumber];
}

sub truncate {
    $input_file_trunc=$_[0];
    chomp ($input_file_trunc);
    open (input_data, "<$input_file_trunc") or die "Couldn't open: $!";
    ($output_file_virus_trunc=$input_file_trunc)=~s/.fastq/_adapter_removal.fastq/;
    open (output_result_virus, ">./$output_file_virus_trunc") or die "Couldn't open: $!";

    my $id1="";
    my $seq="";
    my $id2="";
    my $qual="";
    my $head5=$_[1];
    my $end3=$_[2];

    while (1){
          chomp ($id1=<input_data>);
          chomp ($seq=<input_data>);
          chomp ($id2=<input_data>);
          chomp ($qual=<input_data>);
        if ($head5 ne "none") {                                        
               if ($seq=~/^\Q$head5\E/) {
                   $seq=$';                                                           
                   $qual=substr ($qual, length($head5), length($qual)-length($head5));
                   for ($i=length($seq)-length($end3); $i>=18; $i--){ 
                       if ((substr($seq,$i,length($end3)) eq $end3)) { 
                           $seq=substr ($seq, 0, $i); 
                           $qual=substr ($qual, 0, length($seq)); 
                       if (!($seq=~/N/)) 
                       {
                                   print output_result_virus "$id1\n";
                                   print output_result_virus "$seq\n";
                                   print output_result_virus "$id2\n";
                                   print output_result_virus "$qual\n";
                       } 
                               last;
                          }
                    }
               }
          }else {                                                
              for ($i=length($seq)-length($end3); $i>=18; $i--){
                  if ((substr($seq,$i,length($end3)) eq $end3)) {
                          $seq=substr ($seq, 0, $i);
                          $qual=substr ($qual, 0, length($seq));
                          if (!($seq=~/N/)) {
                               print output_result_virus "$id1\n";
                               print output_result_virus "$seq\n";
                               print output_result_virus "$id2\n";
                               print output_result_virus "$qual\n";
                          }
                          last;
                    }
               }
          }
          if (eof) {
              last;
          }
    }

    close input_data;
    close output_result_virus;
}

sub datapurify {
    my $output_file_datapurify="";
    open (input_data, "<./$_[0]") or die "Couldn't open: $!";
    ($output_file_datapurify=$_[0])=~s/.fastq/_datapurify_\Q$description\E.fastq/;
    open (output_result, ">./$output_file_datapurify") or die "Couldn't open: $!";

    my %sequence_existed=();
    my $id1="";
    my $seq="";
    my $id2="";
    my $qual="";
    my @tempstore=();
    my $i=0;
    while (1){
           chomp ($id1=<input_data>);
           chomp ($seq=<input_data>);
           chomp ($id2=<input_data>);
           chomp ($qual=<input_data>);
        if (!defined($sequence_existed{$seq})) {
               $sequence_existed{$seq}=1;
               $tempstore[$i][0]=$id1;
               $tempstore[$i][1]=$seq;
               $tempstore[$i][2]=$id2;
               $tempstore[$i][3]=$qual;
               $i++;
           } elsif (defined($sequence_existed{$seq})) {
               $sequence_existed{$seq}++;
           }
           if (eof) {
               last;
           }
    }

    for ($i=0;$i<=$#tempstore;$i++) {
           $id1=$tempstore[$i][0];
           $seq=$tempstore[$i][1];
           $id2=$tempstore[$i][2];
           $qual=$tempstore[$i][3];
           if ($description eq "withrepeats") {
                $id1=join "", ($id1,"_",$sequence_existed{$seq});
           } elsif ($description eq "unique") {
                $id1=join "", ($id1,"_1");
           }
           $id1=~s/\s+/./g;
           print output_result "$id1\n";
           print output_result "$seq\n";
           print output_result "$id2\n";
           print output_result "$qual\n";
    }

    close (input_data);
    close (output_result);
    return ($output_file_datapurify);
}

sub normalization_total {
    my $output_file_normalization="";
    open (input_data, "<./$_[0]") or die "Couldn't open: $!";
    ($output_file_normalization=$_[0])=~s/.fastq/_normalized_by_total.fastq/;
    open (output_result, ">./$output_file_normalization") or die "Couldn't open: $!";

    my $lowlimit=$_[1];
    my $highlimit=$_[2];
    my @id1=();
    my @seq=();
    my @id2=();
    my @qual=();
    my $i=0;
    my $totalreads=0;
    while (1){
           chomp ($id1[$i]=<input_data>);
           chomp ($seq[$i]=<input_data>);
           chomp ($id2[$i]=<input_data>);
           chomp ($qual[$i]=<input_data>);
           if ((length($seq[$i])>=$lowlimit) and (length($seq[$i])<=$highlimit)) {
                $id1[$i]=~/_/;
                $totalreads=$totalreads+$';
           }
           $i++;
           if (eof) {
               last;
           }
    }

    my $temp=0;
    for ($i=0;$i<=$#qual;$i++) {
           $id1[$i]=~/_/;
           $temp=1000000*$'/$totalreads;
           $id1[$i]=join "_", ($`,$temp);
           print output_result "$id1[$i]\n";
           print output_result "$seq[$i]\n";
           print output_result "$id2[$i]\n";
           print output_result "$qual[$i]\n";
    }

    close (input_data);
    close (output_result);
    return ($output_file_normalization);
}

sub normalization_totalandgenomesize {
    my $output_file_normalization="";
    open (input_data, "<./$_[0]") or die "Couldn't open: $!";
    $reference=~/.fastq/;
    my $tempstring=$`;
    ($output_file_normalization=$_[0])=~s/.fastq/_normalized_by_totaland\Q$tempstring\E.fastq/;
    open (output_result, ">./$output_file_normalization") or die "Couldn't open: $!";

    my $lowlimit=$_[1];
    my $highlimit=$_[2];
    my @id1=();
    my @seq=();
    my @id2=();
    my @qual=();
    my $i=0;
    my $totalreads=0;
    while (1){
           chomp ($id1[$i]=<input_data>);
           chomp ($seq[$i]=<input_data>);
           chomp ($id2[$i]=<input_data>);
           chomp ($qual[$i]=<input_data>);
           if ((length($seq[$i])>=$lowlimit) and (length($seq[$i])<=$highlimit)) {
                $id1[$i]=~/_/;
                $totalreads=$totalreads+$';
           }
           $i++;
           if (eof) {
               last;
           }
    }

    my $temp=0;
    for ($i=0;$i<=$#qual;$i++) {
           $id1[$i]=~/_/;
           $temp=1000000*$'/$totalreads;
           $temp=1000*$temp/$totalsize;
           $id1[$i]=join "_", ($`,$temp);
           print output_result "$id1[$i]\n";
           print output_result "$seq[$i]\n";
           print output_result "$id2[$i]\n";
           print output_result "$qual[$i]\n";
    }

    close (input_data);
    close (output_result);
    return ($output_file_normalization);
}

sub normalization_microRNA {

    my $output_file_normalization="";
    open (input_data, "<./$_[0]") or die "Couldn't open: $!";
    ($output_file_normalization=$_[0])=~s/.fastq/_normalized_by_microRNA.fastq/;
    open (output_result, ">./$output_file_normalization") or die "Couldn't open: $!";

    my $lowlimit=$_[1];
    my $highlimit=$_[2];
    open (microRNA, "<$microRNA_file") or die "Couldn't open: $!";
    my $id1="";
    my $seq="";
    my %microRNA_ids=();
    while (1){
           chomp ($id1=<microRNA>);
           chomp ($seq=<microRNA>);
           $seq=~s/U/T/g;
           $microRNA_ids{$seq}=$id1;
           if (eof) {
               last;
           }
    }

    my @id1=();
    my @seq=();
    my @id2=();
    my @qual=();
    my $i=0;
    my $totalreads=0;
    my $totalmicroRNA=0;
    while (1){
           chomp ($id1[$i]=<input_data>);
           chomp ($seq[$i]=<input_data>);
           chomp ($id2[$i]=<input_data>);
           chomp ($qual[$i]=<input_data>);
           if ((length($seq[$i])>=$lowlimit) and (length($seq[$i])<=$highlimit)) {
                $id1[$i]=~/_/;
                $totalreads=$totalreads+$';
                if (defined($microRNA_ids{$seq[$i]})) {
                    $totalmicroRNA=$totalmicroRNA+$';
                }
           }
           $i++;
           if (eof) {
               last;
           }
    }

    my $temp=0;
    for ($i=0;$i<=$#qual;$i++) {
           $id1[$i]=~/_/;
           $temp=1000000*$'/$totalmicroRNA;
           $id1[$i]=join "_", ($`,$temp);
           print output_result "$id1[$i]\n";
           print output_result "$seq[$i]\n";
           print output_result "$id2[$i]\n";
           print output_result "$qual[$i]\n";
    }

    close (input_data);
    close (output_result);
    close (microRNA);
    return ($output_file_normalization);
}

sub normalization_microRNAandgenomesize {

    my $output_file_normalization="";
    open (input_data, "<./$_[0]") or die "Couldn't open: $!";
    $reference=~/.fastq/;
    my $tempstring=$`;
    ($output_file_normalization=$_[0])=~s/.fastq/_normalized_by_microRNAand\Q$tempstring\E.fastq/;
    open (output_result, ">./$output_file_normalization") or die "Couldn't open: $!";
	
    my $lowlimit=$_[1];
    my $highlimit=$_[2];
    open (microRNA, "<$microRNA_file") or die "Couldn't open: $!";
    my $id1="";
    my $seq="";
    my %microRNA_ids=();
    while (1){
           chomp ($id1=<microRNA>);
           chomp ($seq=<microRNA>);
           $seq=~s/U/T/g;
           $microRNA_ids{$seq}=$id1;
           if (eof) {
               last;
           }
    }

    my @id1=();
    my @seq=();
    my @id2=();
    my @qual=();
    my $i=0;
    my $totalreads=0;
    my $totalmicroRNA=0;
    while (1){
           chomp ($id1[$i]=<input_data>);
           chomp ($seq[$i]=<input_data>);
           chomp ($id2[$i]=<input_data>);
           chomp ($qual[$i]=<input_data>);
           if ((length($seq[$i])>=$lowlimit) and (length($seq[$i])<=$highlimit)) {
                $id1[$i]=~/_/;
                $totalreads=$totalreads+$';
                if (defined($microRNA_ids{$seq[$i]})) {
                    $totalmicroRNA=$totalmicroRNA+$';
                }
           }
           $i++;
           if (eof) {
               last;
           }
    }

    my $temp=0;
    for ($i=0;$i<=$#qual;$i++) {
           $id1[$i]=~/_/;
           $temp=1000000*$'/$totalmicroRNA;
           $temp=1000*$temp/$genomesize;
           $id1[$i]=join "_", ($`,$temp);
           print output_result "$id1[$i]\n";
           print output_result "$seq[$i]\n";
           print output_result "$id2[$i]\n";
           print output_result "$qual[$i]\n";
    }

    close (input_data);
    close (output_result);
    close (microRNA);
    return ($output_file_normalization);
}

sub normalization_number {

    my $output_file_normalization="";
    open (input_data, "<./$_[0]") or die "Couldn't open: $!";
    ($output_file_normalization=$_[0])=~s/.fastq/_normalized_by_number.fastq/;
    open (output_result, ">./$output_file_normalization") or die "Couldn't open: $!";

    my @id1=();
    my @seq=();
    my @id2=();
    my @qual=();
    my $i=0;
    while (1){
           chomp ($id1[$i]=<input_data>);
           chomp ($seq[$i]=<input_data>);
           chomp ($id2[$i]=<input_data>);
           chomp ($qual[$i]=<input_data>);
           $i++;
           if (eof) {
               last;
           }
    }

    my $temp=0;
    for ($i=0;$i<=$#qual;$i++) {
           $id1[$i]=~/_/;
           $temp=1000000*$'/$_[1];
           $id1[$i]=join "_", ($`,$temp);
           print output_result "$id1[$i]\n";
           print output_result "$seq[$i]\n";
           print output_result "$id2[$i]\n";
           print output_result "$qual[$i]\n";
    }

    close (input_data);
    close (output_result);
    return ($output_file_normalization);
}

sub bowtie_program {
    $input_file_genome=$reference;
    chomp ($input_file_genome);
    $output_file_virus_map=join "_", ("reads_to",$reference,$description);

    $bowtiebuildcommand=join " ", ("bowtie-build -q", $input_file_genome, "./genome.binary");
    system($bowtiebuildcommand);

    $bowtiecommand_virus=join "", ("bowtie -q -a ./genome.binary ./", $output_file_virus_trunc, " -v 0 -p 4 --quiet --suppress 6,7,8 ./", $output_file_virus_map);
    system ($bowtiecommand_virus);

}

sub totalstore_generation {
    my $id1="";
    my $seq="";
    my $id2="";
    my $qual="";

    open (total_truncated, "<./$output_file_virus_trunc") or die "Couldn't open: $!";

    @totalstore=();
    $temp=0;
    while (1){
           chomp ($id1=<total_truncated>);
           chomp ($seq=<total_truncated>);
           chomp ($id2=<total_truncated>);
           chomp ($qual=<total_truncated>);

           $id1=~/_/;
           $temp=$#totalstore+1;
           $totalstore[$temp][0]=$id1;
           $totalstore[$temp][1]="+";
           $totalstore[$temp][2]="total_reads";
           $totalstore[$temp][3]=1;
           $totalstore[$temp][4]=$seq;
           $totalstore[$temp][5]=1/$';

           if (eof) {
               last;
           }
    }

    close (total_truncated);
}

sub length_distribution_percentage {

    use GD;

    my @plus=();
    my @minus=();
    local *store=$_[1];
    my $strand=$_[2];
    my $i=0;
    my $hundredpercent=300;
    my $plus_total=0;
    my $minus_total=0;

    my $output_file=join "", ("Size_distribution_of_whole_library_showpercentage_range", $_[3], "to", $_[4], "nt.png");
    open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

    for ($i=$_[3];$i<=$_[4];$i++) {
          $plus[$i]=0;
          $minus[$i]=0;
    }

    for ($i=0;$i<=$#store;$i++) {
          if (($store[$i][1] eq "+") and ($strand ne "minus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $plus[length($store[$i][4])]=$plus[length($store[$i][4])]+1/$store[$i][5];
               $plus_total=$plus_total+1/$store[$i][5];
          }
          if (($store[$i][1] eq "-") and ($strand ne "plus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $minus[length($store[$i][4])]=$minus[length($store[$i][4])]+1/$store[$i][5];
               $minus_total=$minus_total+1/$store[$i][5];
          }
    }

    my $max_plus=0;
    my $max_minus=0;
    for ($i=$_[3];$i<=$_[4];$i++) {
          if (($strand ne "minus") and ($plus_total!=0)) {
               $plus[$i]=$plus[$i]/$plus_total*$hundredpercent;
               if ($plus[$i]>$max_plus) {
                    $max_plus=$plus[$i];
               }
          }
          if (($strand ne "plus") and ($minus_total!=0)) {
               $minus[$i]=$minus[$i]/$minus_total*$hundredpercent;
               if ($minus[$i]>$max_minus) {
                    $max_minus=$minus[$i];
               }
          }
    }

    my $width_per_group=25;
    my $img=GD::Image->new(($_[4]-$_[3]+1)*$width_per_group+60,$max_plus+$max_minus+200);

    local $white=$img->colorAllocate(255,255,255);
    local $red=$img->colorAllocate(255,0,0);
    local $black=$img->colorAllocate(0,0,0);
    local $blue=$img->colorAllocate(0,0,255);
    $img->fill(1,1,$white);

    $img->line(40,$max_plus+100,($_[4]-$_[3]+1)*$width_per_group+50,$max_plus+100,$black);
    $img->line(50,50,50,$max_plus+$max_minus+150,$black);

    my $temp=0;
    for ($i=0;$i<=0.9;$i=$i+0.1) {
          if ($max_plus>$i*$hundredpercent) {
               $img->line(45,100+$max_plus-($i+0.1)*$hundredpercent,50,100+$max_plus-($i+0.1)*$hundredpercent,$black);
               $temp=($i+0.1)*100;
               $img->string(gdSmallFont,25,95+$max_plus-($i+0.1)*$hundredpercent,"$temp%",$black);
          }
          if ($max_minus>$i*$hundredpercent) {
               $img->line(45,100+$max_plus+($i+0.1)*$hundredpercent,50,100+$max_plus+($i+0.1)*$hundredpercent,$black);
               $temp=($i+0.1)*100;
               $img->string(gdSmallFont,25,95+$max_plus+($i+0.1)*$hundredpercent,"$temp%",$black);
          }
    }

    $temp=int($plus_total);
    $img->string(gdSmallFont,5,40,"total=$plus_total",$black);

    for ($i=$_[3];$i<=$_[4];$i++) {
          $img->filledRectangle(50+($i-$_[3])*$width_per_group,$max_plus+100-$plus[$i],50+($i-$_[3]+1)*$width_per_group-5,$max_plus+100,$red);
          $img->filledRectangle(50+($i-$_[3])*$width_per_group,$max_plus+100+$minus[$i],50+($i-$_[3]+1)*$width_per_group-5,$max_plus+100,$blue);
          $img->string(gdSmallFont,50+($i-$_[3]+1/4)*$width_per_group,30,"$i",$black);
          $img->string(gdSmallFont,50+($i-$_[3]+1/4)*$width_per_group,160+$max_plus+$max_minus,"$i",$black);
    }

    binmode output_result;
    print output_result GD::Image::png($img);
    close (output_result);

}

sub count_reads_number_figure {
    local *store=$_[1];
    my $strand=$_[2];

    open (total_truncated, "<./$output_file_virus_trunc") or die "Couldn't open: $!";
    open (microRNA, "<$microRNA_file") or die "Couldn't open: $!";

    my %microRNA_existed=();
    my %microRNA_ids=();
    my $id1="";
    my $seq="";
    while (1){
           chomp ($id1=<microRNA>);
           chomp ($seq=<microRNA>);
           $seq=~s/U/T/g;
           $microRNA_existed{$seq}=0;
           $microRNA_ids{$seq}=$id1;
           if (eof) {
               last;
           }
    }

    $id1="";
    $seq="";
    my $id2="";
    my $qual="";
    my $total_reads=0;
    my $total_microRNA_reads=0;
    my $total_mapped_reads_plus=0;
    my $total_mapped_reads_minus=0;
    my $total_mapped_reads=0;
    my $i=0;
    my $temp=0;

    while (1){
           chomp ($id1=<total_truncated>);
           chomp ($seq=<total_truncated>);
           chomp ($id2=<total_truncated>);
           chomp ($qual=<total_truncated>);
           if ((length($seq)>=$_[3]) and (length($seq)<=$_[4])) {
                $id1=~/_/;
                $total_reads=$total_reads+$';
                if (defined($microRNA_existed{$seq})) {
                    $total_microRNA_reads=$total_microRNA_reads+$';
                    $microRNA_existed{$seq}=$microRNA_existed{$seq}+$';
                }
           }
           if (eof) {
               last;
           }
    }

    @mirstore=();
    $seq="";
    foreach $seq (keys (%microRNA_existed)) {
             if ($microRNA_existed{$seq}!=0) {
                  $temp=$#mirstore+1;
                  $mirstore[$temp][0]=$microRNA_ids{$seq};
                  if ($mirstore[$temp][0]=~/\*$/) {
                       $mirstore[$temp][1]="-";
                       $mirstore[$temp][4]=&complementary ($seq);
                  } else {
                       $mirstore[$temp][1]="+";
                       $mirstore[$temp][4]=$seq;
                  }
                  $mirstore[$temp][2]=$microRNA_file;
                  $mirstore[$temp][3]=1;
                  $mirstore[$temp][5]=1/$microRNA_existed{$seq};
             }
    }

    for ($i=0;$i<=$#store;$i++) {
          if (($store[$i][1] eq "+") and ($strand ne "minus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $total_mapped_reads=$total_mapped_reads+1/$store[$i][5];
               $total_mapped_reads_plus=$total_mapped_reads_plus+1/$store[$i][5];
          }
          if (($store[$i][1] eq "-") and ($strand ne "plus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $total_mapped_reads=$total_mapped_reads+1/$store[$i][5];
               $total_mapped_reads_minus=$total_mapped_reads_minus+1/$store[$i][5];
          }
    }


    $total_mapped_reads=int($total_mapped_reads);
    $total_mapped_reads_plus=int($total_mapped_reads_plus);
    $total_mapped_reads_minus=int($total_mapped_reads_minus);


    use GD;

    my $plusscale=0;
    my $minusscale=0;
    my $mirscale=0;
    my $everyscale=50;
    my $maxscale=0;
    my $minscale=0;

    my $output_file=join "", ("General_information_range", $_[3], "to", $_[4], "nt.png");
    open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

    if ($total_mapped_reads_plus!=0) {
         $plusscale=(log($total_mapped_reads_plus/$total_microRNA_reads))/(log(10));
         if ($plusscale>$maxscale) {
              $maxscale=$plusscale;
         }
         if ($plusscale<$minscale) {
              $minscale=$plusscale;
         }
    }
    if ($total_mapped_reads_minus!=0) {
         $minusscale=(log($total_mapped_reads_minus/$total_microRNA_reads))/(log(10));
         if ($minusscale>$maxscale) {
              $maxscale=$minusscale;
         }
         if ($minusscale<$minscale) {
              $minscale=$minusscale;
         }
    }
    $maxscale=int($maxscale)+1;
    $minscale=int($minscale)-1;

    my $width_per_group=100;
    my $img=GD::Image->new(4*$width_per_group+110,($maxscale-$minscale+1)*$everyscale+200);

    local $white=$img->colorAllocate(255,255,255);
    local $red=$img->colorAllocate(255,0,0);
    local $black=$img->colorAllocate(0,0,0);
    local $blue=$img->colorAllocate(0,0,255);
    local $green=$img->colorAllocate(0,255,0);
    $img->fill(1,1,$white);

    $img->line(90,($maxscale-$minscale+1)*$everyscale+100,4*$width_per_group+100,($maxscale-$minscale+1)*$everyscale+100,$black);
    $img->line(100,50,100,($maxscale-$minscale+1)*$everyscale+150,$black);

    my $x=50;
    my $y=15;
    $img->filledRectangle($x, $y, $x+10, $y+10, $green);
    $img->filledRectangle($x, $y+10, $x+10, $y+20, $red);
    $img->filledRectangle($x, $y+20, $x+10, $y+30, $blue);
    $img->string(gdSmallFont,$x+12,$y-3,"miRNA",$black);
    $img->string(gdSmallFont,$x+12,$y+7,"plus",$black);
    $img->string(gdSmallFont,$x+12,$y+17,"minus",$black);

    for ($i=$minscale;$i<=$maxscale;$i++) {
          $img->line(95,100+($maxscale-$i)*$everyscale,100,100+($maxscale-$i)*$everyscale,$black);
          $temp=10**$i;
          $img->string(gdSmallFont,50,92+($maxscale-$i)*$everyscale,"$temp",$black);
    }

    $img->string(gdSmallFont,50,155+($maxscale-$minscale+1)*$everyscale,"total reads number=$total_reads",$black);
    $img->string(gdSmallFont,50,165+($maxscale-$minscale+1)*$everyscale,"total virus reads=$total_mapped_reads",$black);

    $img->filledRectangle(120,100+($maxscale-$mirscale)*$everyscale,100+$width_per_group,($maxscale-$minscale+1)*$everyscale+100,$green);
    $img->string(gdSmallFont,120+15,85+($maxscale-$mirscale)*$everyscale,"$total_microRNA_reads",$black);
    if ($total_mapped_reads_plus!=0) {
         $img->filledRectangle(100+1.5*$width_per_group,100+($maxscale-$plusscale)*$everyscale,80+2.5*$width_per_group,($maxscale-$minscale+1)*$everyscale+100,$red);
         $img->string(gdSmallFont,100+15+1.5*$width_per_group,85+($maxscale-$plusscale)*$everyscale,"$total_mapped_reads_plus",$black);
    } else {
         $img->string(gdSmallFont,100+15+1.5*$width_per_group,85+($maxscale-$minscale+1)*$everyscale,"$total_mapped_reads_plus",$black);
    }
    if ($total_mapped_reads_minus!=0) {
         $img->filledRectangle(100+2.5*$width_per_group,100+($maxscale-$minusscale)*$everyscale,80+3.5*$width_per_group,($maxscale-$minscale+1)*$everyscale+100,$blue);
         $img->string(gdSmallFont,100+15+2.5*$width_per_group,85+($maxscale-$minusscale)*$everyscale,"$total_mapped_reads_minus",$black);
    } else {
         $img->string(gdSmallFont,100+15+2.5*$width_per_group,85+($maxscale-$minscale+1)*$everyscale,"$total_mapped_reads_minus",$black);
    }
    binmode output_result;
    print output_result GD::Image::png($img);
    close (output_result);
    close (total_truncated);
    close (microRNA);

}

sub complementary {
    my $i=0;
    my $temp="";
    my $letter="";
    for ($i=length($_[0])-1;$i>=0;$i--) {
          $letter=substr ($_[0], $i, 1);
          if ($letter eq "A") {
               $temp=join "", ($temp, "T");
          }
          if ($letter eq "T") {
               $temp=join "", ($temp, "A");
          }
          if ($letter eq "G") {
               $temp=join "", ($temp, "C");
          }
          if ($letter eq "C") {
               $temp=join "", ($temp, "G");
          }
    }
    return $temp;
}

sub export_miRNA_list_excel {
    local *store=$_[1];
    my $strand=$_[2];
    my $i=0;
    my $j=0;
    my $temp="";
    my $TtoU="";
    my $filenum=1;
    my $morethan=$_[5];
    my $k=0;

    use Spreadsheet::WriteExcel;
    my $excel=Spreadsheet::WriteExcel->new("./result/Mature_microRNAs.xls");
    my $sheet=$excel->add_worksheet();

    $sheet->write(0,0,"ID");
    $sheet->write(0,1,"sequence");
    $sheet->write(0,2,"reads");
    $sheet->write(0,3,"size");
    $sheet->write(1,5,"17 nt ");
    $sheet->write(2,5,"18 nt");
    $sheet->write(3,5,"19 nt");
    $sheet->write(4,5,"20 nt");
    $sheet->write(5,5,"21 nt");
    $sheet->write(6,5,"22 nt");
    $sheet->write(7,5,"23 nt");
    $sheet->write(8,5,"24 nt");
    $sheet->write(9,5,"25 nt");
    $sheet->write(1,6,"=SUMIF(D:D,F2,C:C)");
    $sheet->write(2,6,"=SUMIF(D:D,F3,C:C)");
    $sheet->write(3,6,"=SUMIF(D:D,F4,C:C)");
    $sheet->write(4,6,"=SUMIF(D:D,F5,C:C)");
    $sheet->write(5,6,"=SUMIF(D:D,F6,C:C)");
    $sheet->write(6,6,"=SUMIF(D:D,F7,C:C)");
    $sheet->write(7,6,"=SUMIF(D:D,F8,C:C)");
    $sheet->write(8,6,"=SUMIF(D:D,F9,C:C)");
    $sheet->write(9,6,"=SUMIF(D:D,F10,C:C)");


    my $blank=$excel->add_format();
    $blank->set_color('white');
    $j=0;
    for ($i=0;$i<=$#store;$i++) {
          if ((1/$store[$i][5]>=$_[5]) and ($store[$i][1] eq "+") and ($strand ne "minus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $j++;
               $sheet->write($j,0,"$store[$i][0]");
               $TtoU=$store[$i][4];
               $TtoU=~s/T/U/g;
               $sheet->write($j,1,"$TtoU");
               $temp=1/$store[$i][5];
               $sheet->write($j,2,"$temp");
               $temp=length($store[$i][4]);
               $sheet->write($j,3,"$temp nt");
          }
          if ($j==60000) {
               $j=0;
               $filenum++;
               $excel=Spreadsheet::WriteExcel->new("./result/Mature_microRNAs_file$filenum.xls");
               $sheet=$excel->add_worksheet();

               $sheet->write(0,0,"ID");
               $sheet->write(0,1,"sequence");
               $sheet->write(0,2,"reads");
               $sheet->write(0,3,"size");

               $blank=$excel->add_format();
               $blank->set_color('white');
          }
    }
}

sub lengthdistribution_show_5basebias {
    use GD;

    my @plus=();
    my @minus=();
    local *store=$_[1];
    my $strand=$_[2];
    my $i=0;
    my $base="";
    my $comple="";

    my $output_file=join "", ("Size_distribution_range", $_[3], "to", $_[4], "nt_scale_",$_[5], ".png");
    open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

    for ($i=$_[3];$i<=$_[4];$i++) {
          foreach ("G","C","U","A") {
               $plus[$i]{$_}=0;
               $minus[$i]{$_}=0;
          }
    }

    my $lowcount=21;
    my $highcount=23;
    my $pluscount=0;
    my $minuscount=0;
    my $pluscountu=0;
    my $minuscountu=0;

    for ($i=0;$i<=$#store;$i++) {
          if (($store[$i][1] eq "+") and ($strand ne "minus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $base=substr($store[$i][4],0,1);
               if ($base eq "T") {
                    $base="U";
               }
               $plus[length($store[$i][4])]{$base}=$plus[length($store[$i][4])]{$base}+1/$store[$i][5];
          }
          if (($store[$i][1] eq "-") and ($strand ne "plus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $comple=&complementary ($store[$i][4]);
               $base=substr($comple,0,1);
               if ($base eq "T") {
                    $base="U";
               }
               $minus[length($store[$i][4])]{$base}=$minus[length($store[$i][4])]{$base}+1/$store[$i][5];
          }

          if (($store[$i][1] eq "+") and (length($store[$i][4])>=$lowcount) and (length($store[$i][4])<=$highcount)) {
               $pluscount=$pluscount+1/$store[$i][5];
               $base=substr($store[$i][4],0,1);
               if (($base eq "T") or ($base eq "U")) {
                    $pluscountu=$pluscountu+1/$store[$i][5];
               }
          }
          if (($store[$i][1] eq "-") and (length($store[$i][4])>=$lowcount) and (length($store[$i][4])<=$highcount)) {
               $minuscount=$minuscount+1/$store[$i][5];
               $comple=&complementary ($store[$i][4]);
               $base=substr($comple,0,1);
               if (($base eq "T") or ($base eq "U")) {
                    $minuscountu=$minuscountu+1/$store[$i][5];
               }
          }
    }
    if ($pluscount>0) {
         $pluscountu=$pluscountu/$pluscount*100;
    }
    if ($minuscount>0) {
         $minuscountu=$minuscountu/$minuscount*100;
    }

    my $max_plus=0;
    my $max_minus=0;
    my $temptotal=0;
    for ($i=$_[3];$i<=$_[4];$i++) {
          $temptotal=0;
          foreach ("G","C","U","A") {
               $temptotal=$temptotal+$plus[$i]{$_};
          }
          if ($temptotal>$max_plus) {
               $max_plus=$temptotal;
          }
          $temptotal=0;
          foreach ("G","C","U","A") {
               $temptotal=$temptotal+$minus[$i]{$_};
          }
          if ($temptotal>$max_minus) {
               $max_minus=$temptotal;
          }
    }

    my $yesorno= $_[5];
    my $down= $_[3];
    my $up= $_[4];

    if ($yesorno eq "y") {
         if (($max_plus>0) and ($max_plus<$max_minus)) {
               $max_plus=$max_minus;
         } elsif (($max_minus>0) and ($max_minus<$max_plus)) {
               $max_minus=$max_plus;
         }
    }

    my $scale_plus=0;
    my $scale_minus=0;

    if ($max_plus>=1) {
         $scale_plus=10**(int(log($max_plus)/log(10)));
         if ($max_plus<=(2*$scale_plus)) {
              $scale_plus=$scale_plus/5;
         } elsif ($max_plus>=(8*$scale_plus)) {
              $scale_plus=2*$scale_plus;
         }
    } elsif ($max_plus>0) {
         $scale_plus=10**(int(log($max_plus)/log(10))-1);
         if ($max_plus<=(2*$scale_plus)) {
              $scale_plus=$scale_plus/5;
         } elsif ($max_plus>=(8*$scale_plus)) {
              $scale_plus=2*$scale_plus;
         }
    }

    if ($max_minus>=1) {
         $scale_minus=10**(int(log($max_minus)/log(10)));
         if ($max_minus<=(2*$scale_minus)) {
              $scale_minus=$scale_minus/5;
         } elsif ($max_minus>=(8*$scale_minus)) {
              $scale_minus=2*$scale_minus;
         }
    } elsif ($max_minus>0) {
         $scale_minus=10**(int(log($max_minus)/log(10))-1);
         if ($max_minus<=(2*$scale_minus)) {
              $scale_minus=$scale_minus/5;
         } elsif ($max_minus>=(8*$scale_minus)) {
              $scale_minus=2*$scale_minus;
         }
    }

    if ($scale_plus>0) {
         $max_plus=(int($max_plus/$scale_plus)+1)*$scale_plus;
    }
    if ($scale_minus>0) {
         $max_minus=(int($max_minus/$scale_minus)+1)*$scale_minus;
    }

    my $weight_per_nt_plus=0;
    my $weight_per_nt_minus=0;
    if ($max_plus>0) {
         $weight_per_nt_plus=500/$max_plus;
    }
    if ($max_minus>0) {
         $weight_per_nt_minus=500/$max_minus;
    }

    my $width_per_group=25;
    my $img=GD::Image->new(($_[4]-$_[3]+1)*$width_per_group+100,$max_plus*$weight_per_nt_plus+$max_minus*$weight_per_nt_minus+200);

    local $white=$img->colorAllocate(255,255,255);
    local $red=$img->colorAllocate(255,0,0);
    local $black=$img->colorAllocate(0,0,0);
    local $blue=$img->colorAllocate(0,128,255);
    local $green=$img->colorAllocate(0,255,0);
    local $yellow=$img->colorAllocate(255,200,0);
    local $orange=$img->colorAllocate(255,128,0);
    local $purple=$img->colorAllocate(255,0,128);
    $img->fill(1,1,$white);

    my $x=50;
    my $y=15;
    $img->filledRectangle($x, $y, $x+10, $y+10, $blue);
    $img->filledRectangle($x+10, $y, $x+20, $y+10, $green);
    $img->filledRectangle($x+20, $y, $x+30, $y+10, $purple);
    $img->filledRectangle($x+30, $y, $x+40, $y+10, $orange);
    $img->string(gdSmallFont,$x+2,$y-15,"G",$black);
    $img->string(gdSmallFont,$x+12,$y-15,"C",$black);
    $img->string(gdSmallFont,$x+22,$y-15,"A",$black);
    $img->string(gdSmallFont,$x+32,$y-15,"U",$black);

    $img->line(90,$max_plus*$weight_per_nt_plus+100,($_[4]-$_[3]+1)*$width_per_group+100,$max_plus*$weight_per_nt_plus+100,$black);
    $img->line(100,50,100,$max_plus*$weight_per_nt_plus+$max_minus*$weight_per_nt_minus+150,$black);

    my $tempstring="";
    $tempstring=join "", ("plus ",$lowcount,"~",$highcount,"nt: ",$pluscount,"(1U%=",$pluscountu,"%)");
    $img->string(gdSmallFont,20,27,"$tempstring",$black);
    $tempstring=join "", ("minus ",$lowcount,"~",$highcount,"nt: ",$minuscount,"(1U%=",$minuscountu,"%)");
    $img->string(gdSmallFont,20,155+$max_plus*$weight_per_nt_plus+$max_minus*$weight_per_nt_minus,"$tempstring",$black);

    if ($scale_plus>0) {
         for ($i=$scale_plus;$i<=$max_plus;$i=$i+$scale_plus) {
               $img->line(95,100+$max_plus*$weight_per_nt_plus-$i*$weight_per_nt_plus,100,100+$max_plus*$weight_per_nt_plus-$i*$weight_per_nt_plus,$black);
               $img->string(gdSmallFont,20,92+$max_plus*$weight_per_nt_plus-$i*$weight_per_nt_plus,"$i",$black);
         }
    }
    if ($scale_minus>0) {
         for ($i=$scale_minus;$i<=$max_minus;$i=$i+$scale_minus) {
               $img->line(95,100+$max_plus*$weight_per_nt_plus+$i*$weight_per_nt_minus,100,100+$max_plus*$weight_per_nt_plus+$i*$weight_per_nt_minus,$black);
               $img->string(gdSmallFont,20,92+$max_plus*$weight_per_nt_plus+$i*$weight_per_nt_minus,"$i",$black);
         }
    }

    my $edge=2;
    for ($i=$_[3];$i<=$_[4];$i++) {
          $img->string(gdSmallFont,100+($i-$_[3]+1/4)*$width_per_group,80,"$i",$black);
          $img->string(gdSmallFont,100+($i-$_[3]+1/4)*$width_per_group,110+$max_plus*$weight_per_nt_plus+$max_minus*$weight_per_nt_minus,"$i",$black);

          $temptotal=0;
          foreach ("G","C","U","A") {
               $temptotal=$temptotal+$plus[$i]{$_};
          }
          $img->filledRectangle(100+($i-$_[3])*$width_per_group,$max_plus*$weight_per_nt_plus+100-$temptotal*$weight_per_nt_plus,100+($i-$_[3]+1)*$width_per_group-5,$max_plus*$weight_per_nt_plus+100,$black);

          $temptotal=0;
          foreach ("G","C","U","A") {
               $temptotal=$temptotal+$minus[$i]{$_};
          }
          $img->filledRectangle(100+($i-$_[3])*$width_per_group,$max_plus*$weight_per_nt_plus+100+$temptotal*$weight_per_nt_minus,100+($i-$_[3]+1)*$width_per_group-5,$max_plus*$weight_per_nt_plus+100,$black);

          $temptotal=0;
          $img->filledRectangle(100+($i-$_[3])*$width_per_group+$edge,$max_plus*$weight_per_nt_plus+100-$temptotal*$weight_per_nt_plus-$plus[$i]{"U"}*$weight_per_nt_plus+int($edge/2),100+($i-$_[3]+1)*$width_per_group-5-$edge,$max_plus*$weight_per_nt_plus-$temptotal*$weight_per_nt_plus-int($edge/2)+100,$orange);
          $temptotal=$temptotal+$plus[$i]{"U"};
          $img->filledRectangle(100+($i-$_[3])*$width_per_group+$edge,$max_plus*$weight_per_nt_plus+100-$temptotal*$weight_per_nt_plus-$plus[$i]{"A"}*$weight_per_nt_plus+int($edge/2),100+($i-$_[3]+1)*$width_per_group-5-$edge,$max_plus*$weight_per_nt_plus-$temptotal*$weight_per_nt_plus-int($edge/2)+100,$purple);
          $temptotal=$temptotal+$plus[$i]{"A"};
          $img->filledRectangle(100+($i-$_[3])*$width_per_group+$edge,$max_plus*$weight_per_nt_plus+100-$temptotal*$weight_per_nt_plus-$plus[$i]{"C"}*$weight_per_nt_plus+int($edge/2),100+($i-$_[3]+1)*$width_per_group-5-$edge,$max_plus*$weight_per_nt_plus-$temptotal*$weight_per_nt_plus-int($edge/2)+100,$green);
          $temptotal=$temptotal+$plus[$i]{"C"};
          $img->filledRectangle(100+($i-$_[3])*$width_per_group+$edge,$max_plus*$weight_per_nt_plus+100-$temptotal*$weight_per_nt_plus-$plus[$i]{"G"}*$weight_per_nt_plus+int($edge/2),100+($i-$_[3]+1)*$width_per_group-5-$edge,$max_plus*$weight_per_nt_plus-$temptotal*$weight_per_nt_plus-int($edge/2)+100,$blue);
          $temptotal=$temptotal+$plus[$i]{"G"};

          $temptotal=0;
          $img->filledRectangle(100+($i-$_[3])*$width_per_group+$edge,$max_plus*$weight_per_nt_plus+100+$temptotal*$weight_per_nt_minus+$minus[$i]{"U"}*$weight_per_nt_minus-int($edge/2),100+($i-$_[3]+1)*$width_per_group-5-$edge,$max_plus*$weight_per_nt_plus+$temptotal*$weight_per_nt_minus+int($edge/2)+100,$orange);
          $temptotal=$temptotal+$minus[$i]{"U"};
          $img->filledRectangle(100+($i-$_[3])*$width_per_group+$edge,$max_plus*$weight_per_nt_plus+100+$temptotal*$weight_per_nt_minus+$minus[$i]{"A"}*$weight_per_nt_minus-int($edge/2),100+($i-$_[3]+1)*$width_per_group-5-$edge,$max_plus*$weight_per_nt_plus+$temptotal*$weight_per_nt_minus+int($edge/2)+100,$purple);
          $temptotal=$temptotal+$minus[$i]{"A"};
          $img->filledRectangle(100+($i-$_[3])*$width_per_group+$edge,$max_plus*$weight_per_nt_plus+100+$temptotal*$weight_per_nt_minus+$minus[$i]{"C"}*$weight_per_nt_minus-int($edge/2),100+($i-$_[3]+1)*$width_per_group-5-$edge,$max_plus*$weight_per_nt_plus+$temptotal*$weight_per_nt_minus+int($edge/2)+100,$green);
          $temptotal=$temptotal+$minus[$i]{"C"};
          $img->filledRectangle(100+($i-$_[3])*$width_per_group+$edge,$max_plus*$weight_per_nt_plus+100+$temptotal*$weight_per_nt_minus+$minus[$i]{"G"}*$weight_per_nt_minus-int($edge/2),100+($i-$_[3]+1)*$width_per_group-5-$edge,$max_plus*$weight_per_nt_plus+$temptotal*$weight_per_nt_minus+int($edge/2)+100,$blue);
          $temptotal=$temptotal+$minus[$i]{"G"};
    }

    binmode output_result;
    print output_result GD::Image::png($img);
    close (output_result);

}

sub hit_distribution {
    use GD;
    my @plus=();
    my @minus=();
    local *store=$_[1];
    my $strand=$_[2];
    my $i=0;

    my $output_file=join "", ("Genomic_coverage_range", $_[3], "to", $_[4], "nt.png");
    open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

    my $nt_per_group=1;

    for ($i=1;$i<=int($genomesize/$nt_per_group)+1;$i++) {
          $plus[$i]=0;
          $minus[$i]=0;
    }

    for ($i=0;$i<=$#store;$i++) {
          if (($store[$i][1] eq "+") and ($strand ne "minus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               for ($j=$store[$i][3];$j<=$store[$i][3]+length($store[$i][4])-1;$j++) {
                     $plus[int(($j-1)/$nt_per_group)+1]=$plus[int(($j-1)/$nt_per_group)+1]+1/$store[$i][5];
               }
          }
          if (($store[$i][1] eq "-") and ($strand ne "plus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               for ($j=$store[$i][3];$j<=$store[$i][3]+length($store[$i][4])-1;$j++) {
                     $minus[int(($j-1)/$nt_per_group)+1]=$minus[int(($j-1)/$nt_per_group)+1]+1/$store[$i][5];
               }
          }
    }

    my $max_plus=0;
    my $max_minus=0;
    for ($i=1;$i<=$#plus;$i++) {
          if ($plus[$i]>$max_plus) {
               $max_plus=$plus[$i];
          }
          if ($minus[$i]>$max_minus) {
               $max_minus=$minus[$i];
          }
    }

    my $yesorno= "N";
    if ($yesorno eq "y") {
         if (($max_plus>0) and ($max_plus<$max_minus)) {
               $max_plus=$max_minus;
         } elsif (($max_minus>0) and ($max_minus<$max_plus)) {
               $max_minus=$max_plus;
         }
    }

    my $scale_plus=0;
    my $scale_minus=0;

    if ($max_plus>=1) {
         $scale_plus=10**(int(log($max_plus)/log(10)));
         if ($max_plus<=(2*$scale_plus)) {
              $scale_plus=$scale_plus/5;
         } elsif ($max_plus>=(8*$scale_plus)) {
              $scale_plus=2*$scale_plus;
         }
    } elsif ($max_plus>0) {
         $scale_plus=10**(int(log($max_plus)/log(10))-1);
         if ($max_plus<=(2*$scale_plus)) {
              $scale_plus=$scale_plus/5;
         } elsif ($max_plus>=(8*$scale_plus)) {
              $scale_plus=2*$scale_plus;
         }
    }

    if ($max_minus>=1) {
         $scale_minus=10**(int(log($max_minus)/log(10)));
         if ($max_minus<=(2*$scale_minus)) {
              $scale_minus=$scale_minus/5;
         } elsif ($max_minus>=(8*$scale_minus)) {
              $scale_minus=2*$scale_minus;
         }
    } elsif ($max_minus>0) {
         $scale_minus=10**(int(log($max_minus)/log(10))-1);
         if ($max_minus<=(2*$scale_minus)) {
              $scale_minus=$scale_minus/5;
         } elsif ($max_minus>=(8*$scale_minus)) {
              $scale_minus=2*$scale_minus;
         }
    }

    if ($scale_plus>0) {
         $max_plus=(int($max_plus/$scale_plus)+1)*$scale_plus;
    }
    if ($scale_minus>0) {
         $max_minus=(int($max_minus/$scale_minus)+1)*$scale_minus;
    }

    my $weight_per_nt_plus=0;
    my $weight_per_nt_minus=0;
    if ($max_plus>0) {
         $weight_per_nt_plus=500/$max_plus;
    }
    if ($max_minus>0) {
         $weight_per_nt_minus=500/$max_minus;
    }
    my $width_per_group=1;
    my $img=GD::Image->new($#plus*$width_per_group+int($width_per_group/2)+120,$max_plus*$weight_per_nt_plus+$max_minus*$weight_per_nt_minus+250);

    local $white=$img->colorAllocate(255,255,255);
    local $red=$img->colorAllocate(255,0,0);
    local $black=$img->colorAllocate(0,0,0);
    local $blue=$img->colorAllocate(0,0,255);
    local $green=$img->colorAllocate(0,255,0);
    $img->fill(1,1,$white);

    $img->line(90,$max_plus*$weight_per_nt_plus+100,$#plus*$width_per_group+int($width_per_group/2)+100,$max_plus*$weight_per_nt_plus+100,$black);
    $img->line(100,90,100,$max_plus*$weight_per_nt_plus+$max_minus*$weight_per_nt_minus+110,$black);
    if ($scale_plus>0) {
         for ($i=$scale_plus;$i<=$max_plus;$i=$i+$scale_plus) {
               $img->line(95,100+$max_plus*$weight_per_nt_plus-$i*$weight_per_nt_plus,100,100+$max_plus*$weight_per_nt_plus-$i*$weight_per_nt_plus,$black);
               $img->string(gdSmallFont,20,92+$max_plus*$weight_per_nt_plus-$i*$weight_per_nt_plus,"$i",$black);
         }
    }
    if ($scale_minus>0) {
         for ($i=$scale_minus;$i<=$max_minus;$i=$i+$scale_minus) {
               $img->line(95,100+$max_plus*$weight_per_nt_plus+$i*$weight_per_nt_minus,100,100+$max_plus*$weight_per_nt_plus+$i*$weight_per_nt_minus,$black);
               $img->string(gdSmallFont,20,92+$max_plus*$weight_per_nt_plus+$i*$weight_per_nt_minus,"$i",$black);
         }
    }

    for ($i=1;$i<=$#plus;$i++) {
          $img->filledRectangle(100+$i*$width_per_group-int($width_per_group/2),$max_plus*$weight_per_nt_plus+100-$plus[$i]*$weight_per_nt_plus,100+$i*$width_per_group+int($width_per_group/2),$max_plus*$weight_per_nt_plus+100,$red);
          $img->filledRectangle(100+$i*$width_per_group-int($width_per_group/2),$max_plus*$weight_per_nt_plus+100,100+$i*$width_per_group+int($width_per_group/2),$max_plus*$weight_per_nt_plus+100+$minus[$i]*$weight_per_nt_minus,$blue);
    }

    #Other
    for ($i=1;$i<=$#size;$i++) {
          $img->filledRectangle(100+(int(($jointsize[$i-1]+1-1)/$nt_per_group)+1)*$width_per_group-int($width_per_group/2),$max_plus*$weight_per_nt_plus+$max_minus*$weight_per_nt_minus+130,100+(int(($jointsize[$i-1]+$truesize[$i]-1)/$nt_per_group)+1)*$width_per_group+int($width_per_group/2),$max_plus*$weight_per_nt_plus+$max_minus*$weight_per_nt_minus+140,$black);
    }

    binmode output_result;
    print output_result GD::Image::png($img);
    close (output_result);

}

sub overlappingcheck {
    my %m5top3=();
    my %m5top5=();
    my %m3top3=();
    my %m3top5=();
    my %m5top3plus=();
    my %m5top5plus=();
    my %m3top3plus=();
    my %m3top5plus=();
    my %m5top3minus=();
    my %m5top5minus=();
    my %m3top3minus=();
    my %m3top5minus=();
    my %m5top3total=();
    my %m5top5total=();
    my %m3top3total=();
    my %m3top5total=();

    local *store=$_[1];
    my $type=$_[2];
    my $i=0;
    my $j=0;
    my $ii=0;
    my $temp="";
    my $TtoU="";

    for ($i=(1-$_[4]);$i<=($_[6]-1);$i++) {
          $m5top3{$i}=0;
          $m5top3plus{$i}=0;
          $m5top3minus{$i}=0;
          $m5top3total{$i}=0;
    }
    for ($i=0;$i<=($_[4]+$_[6]-2);$i++) {
          $m5top5{$i}=0;
          $m5top5plus{$i}=0;
          $m5top5minus{$i}=0;
          $m5top5total{$i}=0;
    }
    for ($i=(2-$_[6]-$_[4]);$i<=0;$i++) {
          $m3top3{$i}=0;
          $m3top3plus{$i}=0;
          $m3top3minus{$i}=0;
          $m3top3total{$i}=0;
    }
    for ($i=1-$_[6];$i<=$_[4]-1;$i++) {
          $m3top5{$i}=0;
          $m3top5plus{$i}=0;
          $m3top5minus{$i}=0;
          $m3top5total{$i}=0;
    }

    use Spreadsheet::WriteExcel;
    my $excel=Spreadsheet::WriteExcel->new("./result/Duplexes_details p $_[3] $_[4] m $_[5] $_[6].xls");
    my $sheet=$excel->add_worksheet();

    $sheet->write(0,0,"strand");
    $sheet->write(0,1,"position on +");
    $sheet->write(0,2,"5'head");
    $sheet->write(0,3,"sequence");
    $sheet->write(0,4,"reads");
    $sheet->write(0,5,"size");
    $sheet->write(0,7,"strand");
    $sheet->write(0,8,"position on +");
    $sheet->write(0,9,"5'head");
    $sheet->write(0,10,"sequence");
    $sheet->write(0,11,"reads");
    $sheet->write(0,12,"size");
    $sheet->write(0,14,"m5top5");
    $sheet->write(0,15,"m5top3");
    $sheet->write(0,16,"m3top5");
    $sheet->write(0,17,"m3top3");
    $sheet->write(0,18,"-2 pairs");
    $sheet->write(0,19,"-2 dulplex pairs");
    my $blank=$excel->add_format();
    $blank->set_color('white');
    $j=0;
  
    my @plusheadlist=();
    for ($i=1;$i<=$genomesize;$i++) {
          $plusheadlist[$i]="+";
    }
    for ($i=0;$i<=$#store;$i++) {
          if (($store[$i][1] eq "+") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $plusheadlist[$store[$i][3]]=join "", ($plusheadlist[$store[$i][3]], " ", $i);
          }
    }

    my $k=0;
    my $m=0;
    my @tempdata=();
    for ($i=0;$i<=$#store;$i++) {
          if (($store[$i][1] eq "-") and (length($store[$i][4])>=$_[5]) and (length($store[$i][4])<=$_[6])) {
                for ($k=($store[$i][3]-$_[4]+1);$k<=($store[$i][3]+length($store[$i][4])-1);$k++) {
                      if ($k<1) {
                           next;
                      }
                      @tempdata=();
                      @tempdata=split /\s+/, $plusheadlist[$k];
                      for ($m=1;$m<=$#tempdata;$m++) {
                            $ii=$tempdata[$m];
                            if (
                                ($store[$ii][1] eq "+")
                                and
                                (length($store[$ii][4])>=$_[3])
                                and
                                (length($store[$ii][4])<=$_[4])
                                and
                                (
                                   (($store[$ii][3]>=$store[$i][3]) and ($store[$ii][3]<=($store[$i][3]+length($store[$i][4])-1)))
                                   or
                                   ((($store[$ii][3]+length($store[$ii][4])-1)>=$store[$i][3]) and (($store[$ii][3]+length($store[$ii][4])-1)<=($store[$i][3]+length($store[$i][4])-1)))
                                   or
                                   (($store[$i][3]>=$store[$ii][3]) and ($store[$i][3]<=($store[$ii][3]+length($store[$ii][4])-1)))
                                   or
                                   ((($store[$i][3]+length($store[$i][4])-1)>=$store[$ii][3]) and (($store[$i][3]+length($store[$i][4])-1)<=($store[$ii][3]+length($store[$ii][4])-1)))
                                )
                            ) {
                                  $j++;
                                  $sheet->write($j,0,"-");
                                  $temp=join "", ($store[$i][3], "~", $store[$i][3]+length($store[$i][4])-1);
                                  $sheet->write($j,1,"$temp");
                                  $temp=$store[$i][3]+length($store[$i][4])-1;
                                  $sheet->write($j,2,"$temp");
                                  $TtoU=$store[$i][4];
                                  $TtoU=&complementary($TtoU);
                                  $TtoU=~s/T/U/g;
                                  $sheet->write($j,3,"$TtoU");
                                  $temp=1/$store[$i][5];
                                  $sheet->write($j,4,"$temp");
                                  $temp=length($store[$i][4]);
                                  $sheet->write($j,5,"$temp nt");
                                  $sheet->write($j,7,"+");
                                  $temp=join "", ($store[$ii][3], "~", $store[$ii][3]+length($store[$ii][4])-1);
                                  $sheet->write($j,8,"$temp");
                                  $sheet->write($j,9,"$store[$ii][3]");
                                  $TtoU=$store[$ii][4];
                                  $TtoU=~s/T/U/g;
                                  $sheet->write($j,10,"$TtoU");
                                  $temp=1/$store[$ii][5];
                                  $sheet->write($j,11,"$temp");
                                  $temp=length($store[$ii][4]);
                                  $sheet->write($j,12,"$temp nt");

                                  $temp=$store[$i][3]+length($store[$i][4])-1-$store[$ii][3];
                                  $sheet->write($j,14,"$temp");
                                  $m5top5{$temp}=$m5top5{$temp}+(1/$store[$i][5])*(1/$store[$ii][5]);
                                  $m5top5plus{$temp}=$m5top5plus{$temp}+1/$store[$ii][5];
                                  $m5top5minus{$temp}=$m5top5minus{$temp}+1/$store[$i][5];
                                  $m5top5total{$temp}=$m5top5total{$temp}+1/$store[$i][5]+1/$store[$ii][5];

                                  $temp=$store[$i][3]+length($store[$i][4])-$store[$ii][3]-length($store[$ii][4]);
                                  $sheet->write($j,15,"$temp");
                                  $m5top3{$temp}=$m5top3{$temp}+(1/$store[$i][5])*(1/$store[$ii][5]);
                                  $m5top3plus{$temp}=$m5top3plus{$temp}+1/$store[$ii][5];
                                  $m5top3minus{$temp}=$m5top3minus{$temp}+1/$store[$i][5];
                                  $m5top3total{$temp}=$m5top3total{$temp}+1/$store[$i][5]+1/$store[$ii][5];

                                  $temp=$store[$i][3]-$store[$ii][3];
                                  $sheet->write($j,16,"$temp");
                                  $m3top5{$temp}=$m3top5{$temp}+(1/$store[$i][5])*(1/$store[$ii][5]);
                                  $m3top5plus{$temp}=$m3top5plus{$temp}+1/$store[$ii][5];
                                  $m3top5minus{$temp}=$m3top5minus{$temp}+1/$store[$i][5];
                                  $m3top5total{$temp}=$m3top5total{$temp}+1/$store[$i][5]+1/$store[$ii][5];

                                  $temp=$store[$i][3]-$store[$ii][3]-length($store[$ii][4])+1;
                                  $sheet->write($j,17,"$temp");
                                  $m3top3{$temp}=$m3top3{$temp}+(1/$store[$i][5])*(1/$store[$ii][5]);
                                  $m3top3plus{$temp}=$m3top3plus{$temp}+1/$store[$ii][5];
                                  $m3top3minus{$temp}=$m3top3minus{$temp}+1/$store[$i][5];
                                  $m3top3total{$temp}=$m3top3total{$temp}+1/$store[$i][5]+1/$store[$ii][5];

                            }
                      }
                }
          }
    }

    use GD;
    my $width_per_group=20;
    my $max=0;
    my $weight_per_nt=0;
    my $img;
    my $output_file="";
    my $pensize=5;
    my $scale=0;

    if ($type=~/m5top5/) {
         $output_file=join "", ($title, "_overlappingcheck m5top5 pairs +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt.png");
         open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

         $max=0;
         for ($i=0;$i<=($_[4]+$_[6]-2);$i++) {
               if ($m5top5{$i}>$max) {
                   $max=$m5top5{$i};
               }
         }

         $scale=0;
         if ($max>=1) {
              $scale=10**(int(log($max)/log(10)));
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         } elsif ($max>0) {
              $scale=10**(int(log($max)/log(10))-1);
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         }
         if ($scale>0) {
              $max=(int($max/$scale)+1)*$scale;
         }

         $weight_per_nt=0;
         if ($max!=0) {
              $weight_per_nt=500/$max;
         }

         $img=GD::Image->new(($_[4]+$_[6]-1)*$width_per_group+110,$max*$weight_per_nt+200);
         local $white=$img->colorAllocate(255,255,255);
         local $red=$img->colorAllocate(255,0,0);
         local $black=$img->colorAllocate(0,0,0);
         local $blue=$img->colorAllocate(0,0,255);
         local $green=$img->colorAllocate(0,255,0);
         local $yellow=$img->colorAllocate(255,255,0);
         local $purple=$img->colorAllocate(255,0,255);
         $img->fill(1,1,$white);

         $img->line(90,$max*$weight_per_nt+100,($_[4]+$_[6]-1)*$width_per_group+100,$max*$weight_per_nt+100,$black);
         $img->line(100,50,100,$max*$weight_per_nt+150,$black);
         if ($scale>0) {
              for ($i=$scale;$i<=$max;$i=$i+$scale) {
                    $img->line(95,100+$max*$weight_per_nt-$i*$weight_per_nt,100,100+$max*$weight_per_nt-$i*$weight_per_nt,$black);
                    $img->string(gdSmallFont,20,92+$max*$weight_per_nt-$i*$weight_per_nt,"$i",$black);
              }
         }
         $img->string(gdMediumBoldFont,30,50,"Pairs",$black);
         $temp=join "", ("overlappingcheck m5top5 pairs +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt");
         $img->string(gdMediumBoldFont,30,160+$max*$weight_per_nt,"$temp",$black);

         for ($i=0;$i<=($_[4]+$_[6]-2);$i++) {
               $img->string(gdSmallFont,97+($i+1/2)*$width_per_group,70,$i,$black);
               $img->string(gdSmallFont,97+($i+1/2)*$width_per_group,130+$max*$weight_per_nt,$i,$black);
         }
         for ($i=1;$i<=($_[4]+$_[6]-2);$i++) {
               for ($j=-$pensize;$j<=$pensize;$j++) {
                     $img->line(100+$width_per_group*($i-1/2), $max*$weight_per_nt+100-$m5top5{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+1/2), $max*$weight_per_nt+100-$m5top5{$i}*$weight_per_nt+$j, $black);
               }
         }

         binmode output_result;
         print output_result GD::Image::png($img);
         close (output_result);

         $output_file=join "", ($title, "_overlappingcheck m5top5 reads +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt.png");
         open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

         $max=0;
         for ($i=0;$i<=($_[4]+$_[6]-2);$i++) {
               if ($m5top5total{$i}>$max) {
                   $max=$m5top5total{$i};
               }
         }

         $scale=0;
         if ($max>=1) {
              $scale=10**(int(log($max)/log(10)));
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         } elsif ($max>0) {
              $scale=10**(int(log($max)/log(10))-1);
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         }
         if ($scale>0) {
              $max=(int($max/$scale)+1)*$scale;
         }

         $weight_per_nt=0;
         if ($max!=0) {
              $weight_per_nt=500/$max;
         }

         $img=GD::Image->new(($_[4]+$_[6]-1)*$width_per_group+110,$max*$weight_per_nt+200);
         local $white=$img->colorAllocate(255,255,255);
         local $red=$img->colorAllocate(255,0,0);
         local $black=$img->colorAllocate(0,0,0);
         local $blue=$img->colorAllocate(0,0,255);
         local $green=$img->colorAllocate(0,255,0);
         local $yellow=$img->colorAllocate(255,255,0);
         local $purple=$img->colorAllocate(255,0,255);
         $img->fill(1,1,$white);

         $img->line(90,$max*$weight_per_nt+100,($_[4]+$_[6]-1)*$width_per_group+100,$max*$weight_per_nt+100,$black);
         $img->line(100,50,100,$max*$weight_per_nt+150,$black);
         if ($scale>0) {
              for ($i=$scale;$i<=$max;$i=$i+$scale) {
                    $img->line(95,100+$max*$weight_per_nt-$i*$weight_per_nt,100,100+$max*$weight_per_nt-$i*$weight_per_nt,$black);
                    $img->string(gdSmallFont,20,92+$max*$weight_per_nt-$i*$weight_per_nt,"$i",$black);
              }
         }
         $img->string(gdMediumBoldFont,30,50,"Reads",$black);
         $temp=join "", ("overlappingcheck m5top5 reads +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt");
         $img->string(gdMediumBoldFont,30,160+$max*$weight_per_nt,"$temp",$black);

         for ($i=0;$i<=($_[4]+$_[6]-2);$i++) {
               $img->string(gdSmallFont,97+($i+1/2)*$width_per_group,70,$i,$black);
               $img->string(gdSmallFont,97+($i+1/2)*$width_per_group,130+$max*$weight_per_nt,$i,$black);
         }
         for ($i=1;$i<=($_[4]+$_[6]-2);$i++) {
               for ($j=-$pensize;$j<=$pensize;$j++) {
                     $img->line(100+$width_per_group*($i-1/2), $max*$weight_per_nt+100-$m5top5total{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+1/2), $max*$weight_per_nt+100-$m5top5total{$i}*$weight_per_nt+$j, $purple);
                     $img->line(100+$width_per_group*($i-1/2), $max*$weight_per_nt+100-$m5top5plus{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+1/2), $max*$weight_per_nt+100-$m5top5plus{$i}*$weight_per_nt+$j, $red);
                     $img->line(100+$width_per_group*($i-1/2), $max*$weight_per_nt+100-$m5top5minus{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+1/2), $max*$weight_per_nt+100-$m5top5minus{$i}*$weight_per_nt+$j, $blue);
               }
         }

         binmode output_result;
         print output_result GD::Image::png($img);
         close (output_result);
    }

    if ($type=~/m5top3/) {
         $output_file=join "", ($title, "_overlappingcheck m5top3 pairs +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt.png");
         open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

         $max=0;
         for ($i=(1-$_[4]);$i<=($_[6]-1);$i++) {
               if ($m5top3{$i}>$max) {
               $max=$m5top3{$i};
               }
         }

         $scale=0;
         if ($max>=1) {
              $scale=10**(int(log($max)/log(10)));
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         } elsif ($max>0) {
              $scale=10**(int(log($max)/log(10))-1);
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         }
         if ($scale>0) {
              $max=(int($max/$scale)+1)*$scale;
         }

         $weight_per_nt=0;
         if ($max!=0) {
              $weight_per_nt=500/$max;
         }

         $img=GD::Image->new(($_[4]+$_[6]-1)*$width_per_group+110,$max*$weight_per_nt+200);
         local $white=$img->colorAllocate(255,255,255);
         local $red=$img->colorAllocate(255,0,0);
         local $black=$img->colorAllocate(0,0,0);
         local $blue=$img->colorAllocate(0,0,255);
         local $green=$img->colorAllocate(0,255,0);
         local $yellow=$img->colorAllocate(255,255,0);
         local $purple=$img->colorAllocate(255,0,255);
         $img->fill(1,1,$white);

         $img->line(90,$max*$weight_per_nt+100,($_[4]+$_[6]-1)*$width_per_group+100,$max*$weight_per_nt+100,$black);
         $img->line(100,50,100,$max*$weight_per_nt+150,$black);
         if ($scale>0) {
              for ($i=$scale;$i<=$max;$i=$i+$scale) {
                    $img->line(95,100+$max*$weight_per_nt-$i*$weight_per_nt,100,100+$max*$weight_per_nt-$i*$weight_per_nt,$black);
                    $img->string(gdSmallFont,20,92+$max*$weight_per_nt-$i*$weight_per_nt,"$i",$black);
              }
         }
         $img->string(gdMediumBoldFont,30,50,"Pairs",$black);
         $temp=join "", ("overlappingcheck m5top3 pairs +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt");
         $img->string(gdMediumBoldFont,30,160+$max*$weight_per_nt,"$temp",$black);

         for ($i=(1-$_[4]);$i<=($_[6]-1);$i++) {
               $img->string(gdSmallFont,97+($i+$_[4]-1/2)*$width_per_group,70,$i,$black);
               $img->string(gdSmallFont,97+($i+$_[4]-1/2)*$width_per_group,130+$max*$weight_per_nt,$i,$black);
         }
         for ($i=(2-$_[4]);$i<=($_[6]-1);$i++) {
               for ($j=-$pensize;$j<=$pensize;$j++) {
                     $img->line(100+$width_per_group*($i+$_[4]-3/2), $max*$weight_per_nt+100-$m5top3{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[4]-1/2), $max*$weight_per_nt+100-$m5top3{$i}*$weight_per_nt+$j, $black);
               }
         }

         binmode output_result;
         print output_result GD::Image::png($img);
         close (output_result);

         $output_file=join "", ($title, "_overlappingcheck m5top3 reads +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt.png");
         open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

         $max=0;
         for ($i=(1-$_[4]);$i<=($_[6]-1);$i++) {
               if ($m5top3total{$i}>$max) {
               $max=$m5top3total{$i};
               }
         }

         $scale=0;
         if ($max>=1) {
              $scale=10**(int(log($max)/log(10)));
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         } elsif ($max>0) {
              $scale=10**(int(log($max)/log(10))-1);
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         }
         if ($scale>0) {
              $max=(int($max/$scale)+1)*$scale;
         }

         $weight_per_nt=0;
         if ($max!=0) {
              $weight_per_nt=500/$max;
         }

         $img=GD::Image->new(($_[4]+$_[6]-1)*$width_per_group+110,$max*$weight_per_nt+200);
         local $white=$img->colorAllocate(255,255,255);
         local $red=$img->colorAllocate(255,0,0);
         local $black=$img->colorAllocate(0,0,0);
         local $blue=$img->colorAllocate(0,0,255);
         local $green=$img->colorAllocate(0,255,0);
         local $yellow=$img->colorAllocate(255,255,0);
         local $purple=$img->colorAllocate(255,0,255);
         $img->fill(1,1,$white);

         $img->line(90,$max*$weight_per_nt+100,($_[4]+$_[6]-1)*$width_per_group+100,$max*$weight_per_nt+100,$black);
         $img->line(100,50,100,$max*$weight_per_nt+150,$black);
         if ($scale>0) {
              for ($i=$scale;$i<=$max;$i=$i+$scale) {
                    $img->line(95,100+$max*$weight_per_nt-$i*$weight_per_nt,100,100+$max*$weight_per_nt-$i*$weight_per_nt,$black);
                    $img->string(gdSmallFont,20,92+$max*$weight_per_nt-$i*$weight_per_nt,"$i",$black);
              }
         }
         $img->string(gdMediumBoldFont,30,50,"Reads",$black);
         $temp=join "", ("overlappingcheck m5top3 reads +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt");
         $img->string(gdMediumBoldFont,30,160+$max*$weight_per_nt,"$temp",$black);

         for ($i=(1-$_[4]);$i<=($_[6]-1);$i++) {
               $img->string(gdSmallFont,97+($i+$_[4]-1/2)*$width_per_group,70,$i,$black);
               $img->string(gdSmallFont,97+($i+$_[4]-1/2)*$width_per_group,130+$max*$weight_per_nt,$i,$black);
         }
         for ($i=(2-$_[4]);$i<=($_[6]-1);$i++) {
               for ($j=-$pensize;$j<=$pensize;$j++) {
                     $img->line(100+$width_per_group*($i+$_[4]-3/2), $max*$weight_per_nt+100-$m5top3total{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[4]-1/2), $max*$weight_per_nt+100-$m5top3total{$i}*$weight_per_nt+$j, $purple);
                     $img->line(100+$width_per_group*($i+$_[4]-3/2), $max*$weight_per_nt+100-$m5top3plus{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[4]-1/2), $max*$weight_per_nt+100-$m5top3plus{$i}*$weight_per_nt+$j, $red);
                     $img->line(100+$width_per_group*($i+$_[4]-3/2), $max*$weight_per_nt+100-$m5top3minus{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[4]-1/2), $max*$weight_per_nt+100-$m5top3minus{$i}*$weight_per_nt+$j, $blue);
               }
         }

         binmode output_result;
         print output_result GD::Image::png($img);
         close (output_result);
    }

    if ($type=~/m3top3/) {
         $output_file=join "", ($title, "_overlappingcheck m3top3 pairs +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt.png");
         open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

         $max=0;
         for ($i=(2-$_[4]-$_[6]);$i<=0;$i++) {
               if ($m3top3{$i}>$max) {
               $max=$m3top3{$i};
               }
         }

         $scale=0;
         if ($max>=1) {
              $scale=10**(int(log($max)/log(10)));
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         } elsif ($max>0) {
              $scale=10**(int(log($max)/log(10))-1);
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         }
         if ($scale>0) {
              $max=(int($max/$scale)+1)*$scale;
         }

         $weight_per_nt=0;
         if ($max!=0) {
              $weight_per_nt=500/$max;
         }

         $img=GD::Image->new(($_[4]+$_[6]-1)*$width_per_group+110,$max*$weight_per_nt+200);
         local $white=$img->colorAllocate(255,255,255);
         local $red=$img->colorAllocate(255,0,0);
         local $black=$img->colorAllocate(0,0,0);
         local $blue=$img->colorAllocate(0,0,255);
         local $green=$img->colorAllocate(0,255,0);
         local $yellow=$img->colorAllocate(255,255,0);
         local $purple=$img->colorAllocate(255,0,255);
         $img->fill(1,1,$white);

         $img->line(90,$max*$weight_per_nt+100,($_[4]+$_[6]-1)*$width_per_group+100,$max*$weight_per_nt+100,$black);
         $img->line(100,50,100,$max*$weight_per_nt+150,$black);
         if ($scale>0) {
              for ($i=$scale;$i<=$max;$i=$i+$scale) {
                    $img->line(95,100+$max*$weight_per_nt-$i*$weight_per_nt,100,100+$max*$weight_per_nt-$i*$weight_per_nt,$black);
                    $img->string(gdSmallFont,20,92+$max*$weight_per_nt-$i*$weight_per_nt,"$i",$black);
              }
         }
         $img->string(gdMediumBoldFont,30,50,"Pairs",$black);
         $temp=join "", ("overlappingcheck m3top3 pairs +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt");
         $img->string(gdMediumBoldFont,30,160+$max*$weight_per_nt,"$temp",$black);

         for ($i=(2-$_[4]-$_[6]);$i<=0;$i++) {
               $img->string(gdSmallFont,97+($i+$_[4]+$_[6]-3/2)*$width_per_group,70,$i,$black);
               $img->string(gdSmallFont,97+($i+$_[4]+$_[6]-3/2)*$width_per_group,130+$max*$weight_per_nt,$i,$black);
         }
         for ($i=(3-$_[4]-$_[6]);$i<=0;$i++) {
               for ($j=-$pensize;$j<=$pensize;$j++) {
                     $img->line(100+$width_per_group*($i+$_[4]+$_[6]-5/2), $max*$weight_per_nt+100-$m3top3{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[4]+$_[6]-3/2), $max*$weight_per_nt+100-$m3top3{$i}*$weight_per_nt+$j, $black);
               }
         }

         binmode output_result;
         print output_result GD::Image::png($img);
         close (output_result);

         $output_file=join "", ($title, "_overlappingcheck m3top3 reads +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt.png");
         open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

         $max=0;
         for ($i=(2-$_[4]-$_[6]);$i<=0;$i++) {
               if ($m3top3total{$i}>$max) {
               $max=$m3top3total{$i};
               }
         }

         $scale=0;
         if ($max>=1) {
              $scale=10**(int(log($max)/log(10)));
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         } elsif ($max>0) {
              $scale=10**(int(log($max)/log(10))-1);
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         }
         if ($scale>0) {
              $max=(int($max/$scale)+1)*$scale;
         }

         $weight_per_nt=0;
         if ($max!=0) {
              $weight_per_nt=500/$max;
         }

         $img=GD::Image->new(($_[4]+$_[6]-1)*$width_per_group+110,$max*$weight_per_nt+200);
         local $white=$img->colorAllocate(255,255,255);
         local $red=$img->colorAllocate(255,0,0);
         local $black=$img->colorAllocate(0,0,0);
         local $blue=$img->colorAllocate(0,0,255);
         local $green=$img->colorAllocate(0,255,0);
         local $yellow=$img->colorAllocate(255,255,0);
         local $purple=$img->colorAllocate(255,0,255);
         $img->fill(1,1,$white);

         $img->line(90,$max*$weight_per_nt+100,($_[4]+$_[6]-1)*$width_per_group+100,$max*$weight_per_nt+100,$black);
         $img->line(100,50,100,$max*$weight_per_nt+150,$black);
         if ($scale>0) {
              for ($i=$scale;$i<=$max;$i=$i+$scale) {
                    $img->line(95,100+$max*$weight_per_nt-$i*$weight_per_nt,100,100+$max*$weight_per_nt-$i*$weight_per_nt,$black);
                    $img->string(gdSmallFont,20,92+$max*$weight_per_nt-$i*$weight_per_nt,"$i",$black);
              }
         }
         $img->string(gdMediumBoldFont,30,50,"Reads",$black);
         $temp=join "", ("overlappingcheck m3top3 reads +", $_[3], "~", $_[4], "nt -", $_[5], "~", $_[6], "nt");
         $img->string(gdMediumBoldFont,30,160+$max*$weight_per_nt,"$temp",$black);

         for ($i=(2-$_[4]-$_[6]);$i<=0;$i++) {
               $img->string(gdSmallFont,97+($i+$_[4]+$_[6]-3/2)*$width_per_group,70,$i,$black);
               $img->string(gdSmallFont,97+($i+$_[4]+$_[6]-3/2)*$width_per_group,130+$max*$weight_per_nt,$i,$black);
         }
         for ($i=(3-$_[4]-$_[6]);$i<=0;$i++) {
               for ($j=-$pensize;$j<=$pensize;$j++) {
                     $img->line(100+$width_per_group*($i+$_[4]+$_[6]-5/2), $max*$weight_per_nt+100-$m3top3total{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[4]+$_[6]-3/2), $max*$weight_per_nt+100-$m3top3total{$i}*$weight_per_nt+$j, $purple);
                     $img->line(100+$width_per_group*($i+$_[4]+$_[6]-5/2), $max*$weight_per_nt+100-$m3top3plus{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[4]+$_[6]-3/2), $max*$weight_per_nt+100-$m3top3plus{$i}*$weight_per_nt+$j, $red);
                     $img->line(100+$width_per_group*($i+$_[4]+$_[6]-5/2), $max*$weight_per_nt+100-$m3top3minus{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[4]+$_[6]-3/2), $max*$weight_per_nt+100-$m3top3minus{$i}*$weight_per_nt+$j, $blue);
               }
         }

         binmode output_result;
         print output_result GD::Image::png($img);
         close (output_result);
    }

    if ($type=~/m3top5/) {
         $output_file=join "", ("Duplex_pattern_showpairs_+", $_[3], "to", $_[4], "_-", $_[5], "to", $_[6], "nt.png");
         open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

         $max=0;
         for ($i=(1-$_[6]);$i<=($_[4]-1);$i++) {
               if ($m3top5{$i}>$max) {
               $max=$m3top5{$i};
               }
         }

         $scale=0;
         if ($max>=1) {
              $scale=10**(int(log($max)/log(10)));
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         } elsif ($max>0) {
              $scale=10**(int(log($max)/log(10))-1);
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         }
         if ($scale>0) {
              $max=(int($max/$scale)+1)*$scale;
         }

         $weight_per_nt=0;
         if ($max!=0) {
              $weight_per_nt=500/$max;
         }

         $img=GD::Image->new(($_[4]+$_[6]-1)*$width_per_group+110,$max*$weight_per_nt+200);
         local $white=$img->colorAllocate(255,255,255);
         local $red=$img->colorAllocate(255,0,0);
         local $black=$img->colorAllocate(0,0,0);
         local $blue=$img->colorAllocate(0,0,255);
         local $green=$img->colorAllocate(0,255,0);
         local $yellow=$img->colorAllocate(255,255,0);
         local $purple=$img->colorAllocate(255,0,255);
         $img->fill(1,1,$white);

         $img->line(90,$max*$weight_per_nt+100,($_[4]+$_[6]-1)*$width_per_group+100,$max*$weight_per_nt+100,$black);
         $img->line(100,50,100,$max*$weight_per_nt+150,$black);
         if ($scale>0) {
              for ($i=$scale;$i<=$max;$i=$i+$scale) {
                    $img->line(95,100+$max*$weight_per_nt-$i*$weight_per_nt,100,100+$max*$weight_per_nt-$i*$weight_per_nt,$black);
                    $img->string(gdSmallFont,20,92+$max*$weight_per_nt-$i*$weight_per_nt,"$i",$black);
              }
         }
         $img->string(gdMediumBoldFont,30,50,"Pairs",$black);

         for ($i=(1-$_[6]);$i<=($_[4]-1);$i++) {
               $img->string(gdSmallFont,97+($i+$_[6]-1/2)*$width_per_group,70,$i,$black);
               $img->string(gdSmallFont,97+($i+$_[6]-1/2)*$width_per_group,130+$max*$weight_per_nt,$i,$black);
         }
         for ($i=(2-$_[6]);$i<=($_[4]-1);$i++) {
               for ($j=-$pensize;$j<=$pensize;$j++) {
                     $img->line(100+$width_per_group*($i+$_[6]-3/2), $max*$weight_per_nt+100-$m3top5{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[6]-1/2), $max*$weight_per_nt+100-$m3top5{$i}*$weight_per_nt+$j, $black);
               }
         }

         binmode output_result;
         print output_result GD::Image::png($img);
         close (output_result);

         $output_file=join "", ("Duplex_pattern_showreads_+", $_[3], "to", $_[4], "_-", $_[5], "to", $_[6], "nt.png");
         open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

         $max=0;
         for ($i=(1-$_[6]);$i<=($_[4]-1);$i++) {
               if ($m3top5total{$i}>$max) {
               $max=$m3top5total{$i};
               }
         }

         $scale=0;
         if ($max>=1) {
              $scale=10**(int(log($max)/log(10)));
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         } elsif ($max>0) {
              $scale=10**(int(log($max)/log(10))-1);
              if ($max<=(2*$scale)) {
                   $scale=$scale/5;
              } elsif ($max>=(8*$scale)) {
                   $scale=2*$scale;
              }
         }
         if ($scale>0) {
              $max=(int($max/$scale)+1)*$scale;
         }

         $weight_per_nt=0;
         if ($max!=0) {
              $weight_per_nt=500/$max;
         }

         $img=GD::Image->new(($_[4]+$_[6]-1)*$width_per_group+110,$max*$weight_per_nt+200);
         local $white=$img->colorAllocate(255,255,255);
         local $red=$img->colorAllocate(255,0,0);
         local $black=$img->colorAllocate(0,0,0);
         local $blue=$img->colorAllocate(0,0,255);
         local $green=$img->colorAllocate(0,255,0);
         local $yellow=$img->colorAllocate(255,255,0);
         local $purple=$img->colorAllocate(255,0,255);
         $img->fill(1,1,$white);

         $img->line(90,$max*$weight_per_nt+100,($_[4]+$_[6]-1)*$width_per_group+100,$max*$weight_per_nt+100,$black);
         $img->line(100,50,100,$max*$weight_per_nt+150,$black);
         if ($scale>0) {
              for ($i=$scale;$i<=$max;$i=$i+$scale) {
                    $img->line(95,100+$max*$weight_per_nt-$i*$weight_per_nt,100,100+$max*$weight_per_nt-$i*$weight_per_nt,$black);
                    $img->string(gdSmallFont,20,92+$max*$weight_per_nt-$i*$weight_per_nt,"$i",$black);
              }
         }
         $img->string(gdMediumBoldFont,30,50,"Reads",$black);

         for ($i=(1-$_[6]);$i<=($_[4]-1);$i++) {
               $img->string(gdSmallFont,97+($i+$_[6]-1/2)*$width_per_group,70,$i,$black);
               $img->string(gdSmallFont,97+($i+$_[6]-1/2)*$width_per_group,130+$max*$weight_per_nt,$i,$black);
         }
         for ($i=(2-$_[6]);$i<=($_[4]-1);$i++) {
               for ($j=-$pensize;$j<=$pensize;$j++) {
                     $img->line(100+$width_per_group*($i+$_[6]-3/2), $max*$weight_per_nt+100-$m3top5total{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[6]-1/2), $max*$weight_per_nt+100-$m3top5total{$i}*$weight_per_nt+$j, $purple);
                     $img->line(100+$width_per_group*($i+$_[6]-3/2), $max*$weight_per_nt+100-$m3top5plus{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[6]-1/2), $max*$weight_per_nt+100-$m3top5plus{$i}*$weight_per_nt+$j, $red);
                     $img->line(100+$width_per_group*($i+$_[6]-3/2), $max*$weight_per_nt+100-$m3top5minus{$i-1}*$weight_per_nt+$j, 100+$width_per_group*($i+$_[6]-1/2), $max*$weight_per_nt+100-$m3top5minus{$i}*$weight_per_nt+$j, $blue);
               }
         }

         binmode output_result;
         print output_result GD::Image::png($img);
         close (output_result);
    }
}

sub basebiascurve_oneside {

    local *store=$_[1];
    my $strand=$_[2];
    my $i=0;
    my $j=0;
    my %plusorminus=();
    my $temp="";
    my $hundredpercent=300;
    my $max_plusorminus=0;
    my $comple="";
    my $pensize=2;

    my $output_file=join "", ("Base_bias_", $strand, "_", $_[3], "nt.png");
    open (output_result, ">./result/$output_file") or die "Couldn't open: $!";;

    for ($i=1;$i<=$_[3];$i++) {
          foreach ("G","C","U","A") {
                   $plusorminus[$i]{$_}=0;
          }
    }

    my $total=0;
    for ($i=0;$i<=$#store;$i++) {
          if (($store[$i][1] eq "+") and (length($store[$i][4])==$_[3]) and ($strand ne "minus")) {
                $total=$total+1/$store[$i][5];
                for ($j=1;$j<=$_[3];$j++) {
                      $temp=substr ($store[$i][4], $j-1, 1);
                      if ($temp eq "T") {
                           $temp="U";
                      }
                      $plusorminus[$j]{$temp}=$plusorminus[$j]{$temp}+1/$store[$i][5];
                }
          }
          if (($store[$i][1] eq "-") and (length($store[$i][4])==$_[3]) and ($strand ne "plus")) {
                $total=$total+1/$store[$i][5];
                $comple=&complementary ($store[$i][4]);
                for ($j=1;$j<=$_[3];$j++) {
                      $temp=substr ($comple, $j-1, 1);
                      if ($temp eq "T") {
                           $temp="U";
                      }
                      $plusorminus[$j]{$temp}=$plusorminus[$j]{$temp}+1/$store[$i][5];
                }
          }
    }

    for ($i=1;$i<=$_[3];$i++) {
          if ($total!=0) {
               foreach ("G","C","U","A") {
                        $plusorminus[$i]{$_}=$plusorminus[$i]{$_}/$total*$hundredpercent;
                        if ($plusorminus[$i]{$_}>$max_plusorminus) {
                             $max_plusorminus=$plusorminus[$i]{$_};
                        }
               }
          }
    }

    use GD;
    my $width_per_group=20;
    local $img=GD::Image->new($_[3]*$width_per_group+60,$max_plusorminus+200);

    local $white=$img->colorAllocate(255,255,255);
    local $red=$img->colorAllocate(255,0,0);
    local $black=$img->colorAllocate(0,0,0);
    local $blue=$img->colorAllocate(0,0,255);
    local $green=$img->colorAllocate(0,255,0);
    local $yellow=$img->colorAllocate(255,255,0);
    $img->fill(1,1,$white);

    $img->line(40,$max_plusorminus+100,$_[3]*$width_per_group+50,$max_plusorminus+100,$black);
    $img->line(50,50,50,$max_plusorminus+150,$black);

    my $temptemp=0;
    for ($i=0;$i<=0.9;$i=$i+0.1) {
          if ($max_plusorminus>$i*$hundredpercent) {
               $img->line(45,100+$max_plusorminus-($i+0.1)*$hundredpercent,50,100+$max_plusorminus-($i+0.1)*$hundredpercent,$black);
               $temptemp=($i+0.1)*100;
               $img->string(gdSmallFont,25,95+$max_plusorminus-($i+0.1)*$hundredpercent,"$temptemp%",$black);
          }
    }

    $img->string(gdSmallFont,5,40,int($total),$black);

    my $x=50;
    my $y=15;
    $img->filledRectangle($x, $y, $x+10, $y+10, $red);
    $img->filledRectangle($x+10, $y, $x+20, $y+10, $black);
    $img->filledRectangle($x+20, $y, $x+30, $y+10, $green);
    $img->filledRectangle($x+30, $y, $x+40, $y+10, $blue);
    $img->string(gdSmallFont,$x+2,$y-15,"A",$black);
    $img->string(gdSmallFont,$x+12,$y-15,"U",$black);
    $img->string(gdSmallFont,$x+22,$y-15,"C",$black);
    $img->string(gdSmallFont,$x+32,$y-15,"G",$black);

    for ($i=1;$i<=$_[3];$i++) {
          $img->string(gdSmallFont,47+($i-1+1/2)*$width_per_group,30,$i,$black);
          $img->string(gdSmallFont,47+($i-1+1/2)*$width_per_group,160+$max_plusorminus,$i,$black);
          $img->line(50+$width_per_group*($i-1/2),100+$max_plusorminus,50+$width_per_group*($i-1/2),105+$max_plusorminus,$black);
    }

    for ($i=2;$i<=$_[3];$i++) {
          for ($j=-$pensize;$j<=$pensize;$j++) {
                $img->line(50+$width_per_group*($i-3/2), $max_plusorminus+100-$plusorminus[$i-1]{"A"}+$j, 50+$width_per_group*($i-1/2), $max_plusorminus+100-$plusorminus[$i]{"A"}+$j, $red);
          }
          for ($j=-$pensize;$j<=$pensize;$j++) {
                $img->line(50+$width_per_group*($i-3/2), $max_plusorminus+100-$plusorminus[$i-1]{"U"}+$j, 50+$width_per_group*($i-1/2), $max_plusorminus+100-$plusorminus[$i]{"U"}+$j, $black);
          }
          for ($j=-$pensize;$j<=$pensize;$j++) {
                $img->line(50+$width_per_group*($i-3/2), $max_plusorminus+100-$plusorminus[$i-1]{"C"}+$j, 50+$width_per_group*($i-1/2), $max_plusorminus+100-$plusorminus[$i]{"C"}+$j, $green);
          }
          for ($j=-$pensize;$j<=$pensize;$j++) {
                $img->line(50+$width_per_group*($i-3/2), $max_plusorminus+100-$plusorminus[$i-1]{"G"}+$j, 50+$width_per_group*($i-1/2), $max_plusorminus+100-$plusorminus[$i]{"G"}+$j, $blue);
          }
    }

    my $temppercentage=0;
    print ("$output_file\n");
    for ($i=1;$i<=$_[3];$i++) {
          print ("$i");
          print (" A");
          $temppercentage=$plusorminus[$i]{"A"}/$hundredpercent*100;
          printf "%.2f", $temppercentage;
          print ("%");
          print (" U");
          $temppercentage=$plusorminus[$i]{"U"}/$hundredpercent*100;
          printf "%.2f", $temppercentage;
          print ("%");
          print (" C");
          $temppercentage=$plusorminus[$i]{"C"}/$hundredpercent*100;
          printf "%.2f", $temppercentage;
          print ("%");
          print (" G");
          $temppercentage=$plusorminus[$i]{"G"}/$hundredpercent*100;
          printf "%.2f", $temppercentage;
          print ("%");
          print ("\n");
    }

    binmode output_result;
    print output_result GD::Image::png($img);
    close (output_result);

}

sub export_siRNA_list_excel {
    local *store=$_[1];
    my $strand=$_[2];
    my $i=0;
    my $j=0;
    my $temp="";
    my $TtoU="";
    my $filenum=1;
    my $morethan=$_[5];
    my $k=0;

    use Spreadsheet::WriteExcel;
    my $excel=Spreadsheet::WriteExcel->new("./result/Small_RNA_details_range$_[3]to$_[4]nt.xls");
    my $sheet=$excel->add_worksheet();

    $sheet->write(0,0,"segment order");
    $sheet->write(0,1,"strand");
    $sheet->write(0,2,"position on +");
    $sheet->write(0,3,"5'head allseg");
    $sheet->write(0,4,"sequence");
    $sheet->write(0,5,"reads");
    $sheet->write(0,6,"size");
    $sheet->write(0,7,"summary");
    $sheet->write(0,8,"from terminal");

    $sheet->write(0,10,"segment order");
    $sheet->write(0,11,"strand");
    $sheet->write(0,12,"position on +");
    $sheet->write(0,13,"5'head allseg");
    $sheet->write(0,14,"sequence");
    $sheet->write(0,15,"reads");
    $sheet->write(0,16,"size");
    $sheet->write(0,17,"summary");
    $sheet->write(0,18,"from terminal");


    my $blank=$excel->add_format();
    $blank->set_color('white');
    $j=0;
    my $segmentorder=0;
    my $end5=0;
    my $end3=0;
    for ($i=0;$i<=$#store;$i++) {
          if ((1/$store[$i][5]>=$_[5]) and ($store[$i][1] eq "+") and ($strand ne "minus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $j++;
               for ($k=1;$k<=$segmentnumber;$k++) {
                     if ($segmentID[$k]=~/\Q$store[$i][2]\E/) {
                          $segmentorder=$k;
                          last;
                     }
               }
               $sheet->write($j,10,"$segmentorder");
               $sheet->write($j,11,"+");

               $temp=$store[$i][3];
               $temp=$temp-$jointsize[$segmentorder-1];
               $temp=join "", ($temp, "~", $temp+length($store[$i][4])-1);

               $sheet->write($j,12,"$temp");
               $sheet->write($j,13,"$store[$i][3]");
               $TtoU=$store[$i][4];
               $TtoU=~s/T/U/g;
               $sheet->write($j,14,"$TtoU");
               $temp=1/$store[$i][5];
               $sheet->write($j,15,"$temp");
               $temp=length($store[$i][4]);
               $sheet->write($j,16,"$temp nt");
               $temp=$store[$i][3]+length($store[$i][4])-1;
               $sheet->write($j,3,"$temp",$blank);

               $end5=$store[$i][3]-$jointsize[$segmentorder-1];
               $end3=$jointsize[$segmentorder-1]+$truesize[$segmentorder]-$store[$i][3]-length($store[$i][4])+2;
               $temp=join "",("(", int(1/$store[$i][5]), ")");
               if ($end5<=$end3) {
                     $sheet->write($j,18,"$end5");
                     $temp=join " ",($end5, $TtoU, $end5+length($store[$i][4])-1, $temp);
               } else {
                     $sheet->write($j,18,"$end3");
                     $temp=join " ",($end3+length($store[$i][4])-1, $TtoU, $end3, $temp);
               }
               $sheet->write($j,17,"$temp");
          }
          if ((1/$store[$i][5]>=$_[5]) and ($store[$i][1] eq "-") and ($strand ne "plus") and (length($store[$i][4])>=$_[3]) and (length($store[$i][4])<=$_[4])) {
               $j++;
               for ($k=1;$k<=$segmentnumber;$k++) {
                     if ($segmentID[$k]=~/\Q$store[$i][2]\E/) {
                          $segmentorder=$k;
                          last;
                     }
               }
               $sheet->write($j,0,"$segmentorder");
               $sheet->write($j,1,"-");

               $temp=$store[$i][3];
               $temp=$temp-$jointsize[$segmentorder-1];
               $temp=join "", ($temp, "~", $temp+length($store[$i][4])-1);

               $sheet->write($j,2,"$temp");
               $temp=$store[$i][3]+length($store[$i][4])-1;
               $sheet->write($j,3,"$temp");
               $TtoU=$store[$i][4];
               $TtoU=&complementary($TtoU);
               $TtoU=~s/T/U/g;
               $sheet->write($j,4,"$TtoU");
               $temp=1/$store[$i][5];
               $sheet->write($j,5,"$temp");
               $temp=length($store[$i][4]);
               $sheet->write($j,6,"$temp nt");
               $sheet->write($j,13,"$store[$i][3]",$blank);

               $temp="";
               for ($k=length($TtoU)-1;$k>=0;$k--) {
                     $temp=join "",($temp,substr($TtoU,$k,1));
               }
               $TtoU=$temp;

               $end5=$store[$i][3]-$jointsize[$segmentorder-1];
               $end3=$jointsize[$segmentorder-1]+$truesize[$segmentorder]-$store[$i][3]-length($store[$i][4])+2;
               $temp=join "",("(", int(1/$store[$i][5]), ")");
               if ($end5<=$end3) {
                     $sheet->write($j,8,"$end5");
                     $temp=join " ",($end5, $TtoU, $end5+length($store[$i][4])-1, $temp);
               } else {
                     $sheet->write($j,8,"$end3");
                     $temp=join " ",($end3+length($store[$i][4])-1, $TtoU, $end3, $temp);
               }
               $sheet->write($j,7,"$temp");
               $TtoU=$store[$i][4];
               $TtoU=&complementary($TtoU);
               $TtoU=~s/T/U/g;
          }
          if ($j==60000) {
               $j=0;
               $filenum++;
               $excel=Spreadsheet::WriteExcel->new("./result/Small_RNA_details_range$_[3]to$_[4]nt_file$filenum.xls");
               $sheet=$excel->add_worksheet();

               $sheet->write(0,0,"segment order");
               $sheet->write(0,1,"strand");
               $sheet->write(0,2,"position on +");
               $sheet->write(0,3,"5'head allseg");
               $sheet->write(0,4,"sequence");
               $sheet->write(0,5,"reads");
               $sheet->write(0,6,"size");
               $sheet->write(0,7,"summary");
               $sheet->write(0,8,"from terminal");

               $sheet->write(0,10,"segment order");
               $sheet->write(0,11,"strand");
               $sheet->write(0,12,"position on +");
               $sheet->write(0,13,"5'head allseg");
               $sheet->write(0,14,"sequence");
               $sheet->write(0,15,"reads");
               $sheet->write(0,16,"size");
               $sheet->write(0,17,"summary");
               $sheet->write(0,18,"from terminal");

               $blank=$excel->add_format();
               $blank->set_color('white');
          }
    }
}

sub data_process {

    open (input_data, "<./$_[0]") or die "Couldn't open: $!";

    my %multiple_hit=();
    @virstore=();
    my $temp="";
    my $i=0;
    my $j=0;
    my $segmentorder=0;

    my @tempdata=();
    if (-s "./$_[0]") {
         while (1) {
                chomp ($temp=<input_data>);
                @tempdata=();
                @tempdata=split /\s+/, $temp;
                $virstore[$i]=[@tempdata];
                $multiple_hit{$virstore[$i][0]}++;
                if (eof) {
                    last;
                }
                $i++;
         }
    }

    for ($i=0;$i<=$#virstore;$i++) {
          $virstore[$i][0]=~/_/;
          $virstore[$i][5]=$multiple_hit{$virstore[$i][0]}/$'; ## now (1/[i][5]) is the value for calculating
          $virstore[$i][3]++; ##because bowtie returns result starting from [0]
          for ($j=1;$j<=$segmentnumber;$j++) {
                if ($segmentID[$j]=~/\Q$virstore[$i][2]\E/) {
                     $segmentorder=$j;
                     last;
                }
          }
          $virstore[$i][3]=$virstore[$i][3]+$jointsize[$segmentorder-1];
    }

    close (input_data);

     &totalstore_generation ($_[0],$_[1]);
     &length_distribution_percentage ($_[0], *totalstore, "plus & minus", 18, 28, $_[1]);
     &length_distribution_percentage ($_[0], *totalstore, "plus & minus", 18, 35, $_[1]);
     &count_reads_number_figure ($_[0], *virstore, "plus & minus", 18, 28, $_[1]);
     &export_miRNA_list_excel ($_[0], *mirstore, "plus & minus", 18, 32, 0, $_[1]);
     &lengthdistribution_show_5basebias ($_[0], *virstore, "plus & minus", 18, 28, "y", $_[1]);
     &lengthdistribution_show_5basebias ($_[0], *virstore, "plus & minus", 18, 28, "n", $_[1]);
     &lengthdistribution_show_5basebias ($_[0], *virstore, "plus & minus", 18, 35, "y", $_[1]);
     &lengthdistribution_show_5basebias ($_[0], *virstore, "plus & minus", 18, 35, "n", $_[1]);
     &hit_distribution ($_[0], *virstore, "plus & minus", 18, 28);
     &hit_distribution ($_[0], *virstore, "plus & minus", 21, 23);
     &overlappingcheck ($_[0], *virstore, "m3top5", 22, 22, 22, 22, $_[1]);
     &export_siRNA_list_excel ($_[0], *virstore, "plus & minus", 18, 28, 0, $_[1]);
}


if (!(-d "result")) {
    system("mkdir ./result");
    system("mkdir ./result1");
    system("mkdir ./result2");
}
use 5.010;
@databasename_files=glob"*.fastq";
     foreach $databasename_files(@databasename_files) {
          if($databasename_files=~/withrepeats.fastq/){
               $database_filename=join "",$`,$&;
               print ("$database_filename is the input fastq file\n") ;
               last;}
          if(eof){last}}

if($database_filename eq"")
     {foreach $databasename_files(@databasename_files) {
          if($databasename_files=~/.fastq/){
               $database_filename=join "",$`,$&;
               print ("$database_filename is the input fastq file\n\n") ;
               last;}
          if(eof){last}}
     }

$reply="";
$output_file_virus_trunc="";
$output_file_virus_trunc=$database_filename;
$description="withrepeats";
if ($output_file_virus_trunc=~/adapter_removal_datapurify_withrepeats.fastq$/) 
     {print("trimmed file detected, skipping trimminng programs \n" );
     }
else{$reply="y";}

if ($reply eq "y") {
     print ("The file is not trimmed, autotrimming running\n");
     print ("5' none and 3' AGATCGGA is default setting, see line 2388\n");
     print ("Trimming, please wait\n");
     $ada5="none";
     $ada3="AGATCGGA";
     &truncate($output_file_virus_trunc,$ada5,$ada3); 
     $output_file_virus_trunc=&datapurify ($output_file_virus_trunc);       
}
print ("The mapping will be done allowing no mismatches.\n\n");
@reference_files=glob"*.exo.fa";
foreach $reference_files(@reference_files) {
     if($reference_files=~/.exo.fa/){
          $reference=join "",$`,$&;
          print ("$reference is the input reference file of exogenouse dsRNA\n") ;
          last;}
     if(eof){last}}
$genomesize=&genome_size($reference);
@microRNA_file=glob"*.miRNA.fa";
     foreach $microRNA_test(@microRNA_file) {
     if($microRNA_test=~/.miRNA.fa/){
     $microRNA_file=join "",$`,$&;
          print ("$microRNA_file is the input reference file of miRNA \n\n") ;
          last;}
     if(eof){last}}

##generating result files
print("processing no normed data\n");
     &bowtie_program;
     &data_process ($output_file_virus_map);
     print ("Done.\n");
     use File::Copy ;
     rename ("./result","./No_norm");
     rename ("./result1","./result");

print("processing data normalized by total 18_28nt\n");
     $output_file_virus_trunc=&normalization_total ($output_file_virus_trunc, 18, 28);
     &bowtie_program;
     &data_process ($output_file_virus_map);
     print ("Done.\n");
     use File::Copy ;
     rename ("./result","./Norm_to_total18_28nt");
     rename ("./result2","./result");

print("processing data normalized by miRNA\n");
     $output_file_virus_trunc=&normalization_microRNA ($output_file_virus_trunc, 18, 32);
     &bowtie_program;
     &data_process ($output_file_virus_map);
     print ("Done.\n");
     use File::Copy ;
     rename ("./result","./Norm_to_miRNA");

unlink glob "*.ebwt";
unlink glob "*_by_total.fastq";
unlink glob "*_by_microRNA.fastq";
unlink glob "*adapter_removal.fastq";
unlink glob "*.txt_withrepeats";
print ("\nProgram done. Please find all results in 'No_Norm' 'Norm_to_total' and 'Norm_to_miRNA' folder.\n");

exit;