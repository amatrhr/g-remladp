#!/usr/bin/perl -w

$dummy = 0;

use Getopt::Std;

getopts('c:h:m:s:');

#Read in map file to convert genomic position identifiers to rs number identifiers
open(MAP, "$opt_m")||die "Can't find file $opt_m\n";
while (<MAP>) {
  chomp;
  @temp = split;
  $rsid{$temp[0]}{$temp[3]} = $temp[1];
}
close MAP;
print "Read Map file\n";

#Read in sample file

$count = 0;
open(SAMPLE,"$opt_s")||die "Can't find file $opt_s\n";
<SAMPLE>;
<SAMPLE>;
while (<SAMPLE>) {
  chomp;
  @temp = split;
  $famid = substr($temp[0],0,5);
  $qlet = substr($temp[0],5,1);
  if($qlet eq "M") {
    push @famids, $famid;                       #Just based on mothers because if no mother then pattern of transmission cannot be determined
    $mother_id{$famid} = $temp[0];
    $motherpos{$famid}= $count;
  } else {
    $qlet{$famid} = $qlet;
    $child_id{$famid} = $temp[0];
    $childpos{$famid}= $count;
#    $sex{$famid} = $temp[5];
  }
  $count++
}
close SAMPLE;
print "Read sample file\n";

foreach $famid (@famids) {
  if(defined($mother_id{$famid}) && defined($child_id{$famid})) {
    push @selected_fam, $famid;                 #By definition this will be one mother and one child (even in the case of twins)
  }
}


#Read in haplotypes
$count = 0;
open HAP, '-|', '/bin/zcat', "$opt_h"||die "Can't find file $opt_h\n";
while (<HAP>) {
  chomp;
  ($dummy, $snp, $pos, $a1, $a2, @temp) = split;
#  if(!defined($selected{$snp})) {              #Select only SNPs that have been genotyped
#    next;
#  }
  unless(!defined($rsid{$opt_c}{$pos})) {
    $snp = $rsid{$opt_c}{$pos};
  } else {
  print "Warning $snp does not appear to have an rs identifier\n"
  }
  push @snps, $snp;
  $pos{$snp} = $pos;
  $allele{$snp}{"0"} = $a1;
  $allele{$snp}{"1"} = $a2;
  $allele{$snp}{"5"} = 0;
  $allele{$snp}{"9"} = 0;
  $genotypes{$snp} = join "", @temp;
  $sum = 0;                                     #Calculate allele frequencies for each SNP
  for (@temp) {
    $sum += $_;
  }
  $num = @temp;
  $freq{$snp} = 1 - ($sum/$num);                #Freq of A1
  print "Read in $snp. Count = $count\n";
  $snp_number{$snp} = $count;
  $count++;
}
close HAP;
$nsnp = $count-1;
#Now convert into useful format
#should just do for selected individuals***
foreach $famid (@selected_fam) {

  $temp = $motherpos{$famid}*2;
  $template1 = "x$temp"."A1";
  $temp++;
  $template2 = "x$temp"."A1";

  @hap1 = ();
  @hap2 = ();

  foreach $snp (@snps) {

    ($value1) = unpack($template1,$genotypes{$snp});
    ($value2) = unpack($template2,$genotypes{$snp});

    push @hap1, $value1;
    push @hap2, $value2;

  }

  $hap1_mother{$famid} = join "", @hap1;
  $hap2_mother{$famid} = join "", @hap2;

  $temp = $childpos{$famid}*2;
  $template1 = "x$temp"."A1";
  $temp++;
  $template2 = "x$temp"."A1";

  @hap1 = ();
  @hap2 = ();

  foreach $snp (@snps) {

    ($value1) = unpack($template1,$genotypes{$snp});
    ($value2) = unpack($template2,$genotypes{$snp});

    push @hap1, $value1;
    push @hap2, $value2;

  }

  $hap1_child{$famid} = join "", @hap1;
  $hap2_child{$famid} = join "", @hap2;

  print "Resolved family $famid\n";

}

%genotypes=();




#Determine transmissions
#Resolved documents maternal transmission 1 = hap1; 2 = hap2; 0 = not known, but doesn't matter; 9 = unknown and needs resolving; 5 = incompatible genotypes

foreach $famid (@selected_fam) {

  foreach $snp (@snps) {

  $temp = $snp_number{$snp};
  $template1 = "x$temp"."A1";

    ($hap1mum) = unpack($template1,$hap1_mother{$famid});
    ($hap2mum) = unpack($template1,$hap2_mother{$famid});
    ($hap1child) = unpack($template1,$hap1_child{$famid});
    ($hap2child) = unpack($template1,$hap2_child{$famid});

    if($hap1mum == 0 && $hap2mum == 0) {

      if($hap1child == 0 && $hap2child == 0) {
        $maternal = 0;
        $paternal = 0;
        $resolved = 0;
      } elsif($hap1child == 0 && $hap2child == 1 || $hap1child == 1 && $hap2child == 0) {
        $maternal = 0;
        $paternal = 1;
        $resolved = 0;
      } elsif($hap1child == 1 && $hap2child == 1) {
        $maternal = 5;
        $paternal = 5;
        $resolved = 5;
      } else {
        die "Unknown genotype for individual $famid at SNP $snp\n";
      }

    } elsif($hap1mum == 1 && $hap2mum == 1) {

      if($hap1child == 0 && $hap2child == 0) {
        $maternal = 5;
        $paternal = 5;
        $resolved = 5;
      } elsif($hap1child == 0 && $hap2child == 1 || $hap1child == 1 && $hap2child == 0) {
        $maternal = 1;
        $paternal = 0;
        $resolved = 0;
      } elsif($hap1child == 1 && $hap2child == 1) {
        $maternal = 1;
        $paternal = 1;
        $resolved = 0;
      } else {
        die "Unknown genotype for individual $famid at SNP $snp\n";
      }

    } elsif($hap1mum == 0 && $hap2mum == 1) {

      if($hap1child == 0 && $hap2child == 0) {
        $maternal = 0;
        $paternal = 0;
        $resolved = 1;
      } elsif($hap1child == 0 && $hap2child == 1 || $hap1child == 1 && $hap2child == 0) {
        $maternal = 9;
        $paternal = 9;
        $resolved = 9;
      } elsif($hap1child == 1 && $hap2child == 1) {
        $maternal = 1;
        $paternal = 1;
        $resolved = 2;
      } else {
        die "Unknown genotype for individual $famid at SNP $snp\n";
      }

    } elsif($hap1mum == 1 && $hap2mum == 0) {

      if($hap1child == 0 && $hap2child == 0) {
        $maternal = 0;
        $paternal = 0;
        $resolved = 2;
      } elsif($hap1child == 0 && $hap2child == 1 || $hap1child == 1 && $hap2child == 0) {
        $maternal = 9;
        $paternal = 9;
        $resolved = 9;
      } elsif($hap1child == 1 && $hap2child == 1) {
        $maternal = 1;
        $paternal = 1;
        $resolved = 1;
      } else {
        die "Unknown genotype for individual $famid at SNP $snp\n";
      }

    } else {
        die "Unknown genotype for individual $famid at SNP $snp\n";
    }

    push @maternal, $maternal;
    push @paternal, $paternal;
    push @resolved, $resolved;

  }

  $maternal_hap{$famid} = join "", @maternal;
  $paternal_hap{$famid} = join "", @paternal;
  $resolved{$famid} = join "", @resolved;
  @maternal = ();
  @paternal = ();
  @resolved = ();
  print "Transmission determined $famid\n";

}


#Print map file
open(MAP,">$opt_c.map");
foreach $snp (@snps) {
  print MAP "$opt_c $snp 0 $pos{$snp}\n";
}
close MAP;




print "Filling in gaps\n";

foreach $famid (@selected_fam) {


  @hap1_mother = split('', $hap1_mother{$famid});
  @hap2_mother = split('', $hap2_mother{$famid});
  @maternal = split('', $maternal_hap{$famid});
  @paternal = split('', $paternal_hap{$famid});
  @resolved = split('', $resolved{$famid});

  #Come in from left
  $left = $snps[0];
  $direction = 1;
  $pos{$left} = -999999999;
  for($i=0;$i<=$nsnp;$i++) {

    if($resolved[$i] == 1 || $resolved[$i] == 2) {
      $left = $snps[$i];
      $direction = $resolved[$i];
    } elsif($resolved[$i] == 9) {
      $dist_left{$snps[$i]} = $pos{$snps[$i]} - $pos{$left};
      $hap_left{$snps[$i]} = $direction;
    }

  }

  #Come in from right
  $right = $snps[$nsnp];
  $direction = 1;
  $pos{$right} = 999999999;
  for($i=$nsnp;$i>=0;$i--) {

    if($resolved[$i] == 1 || $resolved[$i] == 2) {
      $right = $snps[$i];
      $direction = $resolved[$i];
    } elsif($resolved[$i] == 9) {
      $dist_right{$snps[$i]} = $pos{$right} - $pos{$snps[$i]};
      $hap_right{$snps[$i]} = $direction;
    }

  }

  #Resolve unknown transmissions. Remember these are always 1/2 1/2 transmissions
  for($i=0;$i<=$nsnp;$i++) {

    if($resolved[$i] == 9) {
      if($dist_left{$snps[$i]} <= $dist_right{$snps[$i]}) {
        if($hap_left{$snps[$i]} == 1) {
          $maternal[$i] = $hap1_mother[$i];
          $paternal[$i] = $hap2_mother[$i]; #This formulation OK because heterozygous haplotypes
          $resolved[$i] = 1;
        } else {
          $maternal[$i] = $hap2_mother[$i];
          $paternal[$i] = $hap1_mother[$i];
          $resolved[$i] = 2;
        }
      } else {
        if($hap_right{$snps[$i]} == 1) {
          $maternal[$i] = $hap1_mother[$i];
          $paternal[$i] = $hap2_mother[$i];
          $resolved[$i] = 1;
        } else {
          $maternal[$i] = $hap2_mother[$i];
          $paternal[$i] = $hap1_mother[$i];
          $resolved[$i] = 2;
        }
      }
    }

  }

  $maternal_hap{$famid} = join "", @maternal;
  $paternal_hap{$famid} = join "", @paternal;
  $resolved{$famid} = join "", @resolved;

  print "Resolved family $famid\n";

}


#Print mlfino file
open(MLINFO,">$opt_c.mlinfo");
print MLINFO "SNP Al1 Al2 Freq1 MAF Quality Rsq\n";
foreach $snp (@snps) {
  if($freq{$snp} <= 0.5) {
    print MLINFO "$snp $allele{$snp}{0} $allele{$snp}{1} $freq{$snp} $freq{$snp} 1 1\n";
  } else {
    $maf = 1-$freq{$snp};
    print MLINFO "$snp $allele{$snp}{0} $allele{$snp}{1} $freq{$snp} $maf 1 1\n";
  }
}
close MLINFO;

#Print mldose files
open(PED,">$opt_c.mldose");
open(MATERNAL, ">maternal_$opt_c.mldose");
open(PATERNAL, ">paternal_$opt_c.mldose");
#open(RESOLVED, ">resolved.txt"); #Check to show algorithm working correctly
foreach $famid (@selected_fam) {

  @maternal = split('', $maternal_hap{$famid});
  @paternal = split('', $paternal_hap{$famid});
#  @resolved = split('', $resolved{$famid});

  print PED "$famid->$qlet{$famid} MLDOSE ";
  print MATERNAL "$famid->$qlet{$famid} MLDOSE ";
  print PATERNAL "$famid->$qlet{$famid} MLDOSE ";
#  print RESOLVED "$famid->$qlet{$famid} MLDOSE ";

  for($i=0;$i<=$nsnp;$i++) {

#    print RESOLVED "$resolved[$i] ";

    if($maternal[$i] == 0 && $paternal[$i] == 0) {
      print PED "0 ";
      print MATERNAL "0 ";
      print PATERNAL "0 ";
    } elsif ($maternal[$i] == 0 && $paternal[$i] == 1){
      print PED "-1 ";
      print MATERNAL "0 ";
      print PATERNAL "2 ";
    } elsif ($maternal[$i] == 1 && $paternal[$i] == 0) {
      print PED "1 ";
      print MATERNAL "2 ";
      print PATERNAL "0 ";
    } elsif ($maternal[$i] == 1 && $paternal[$i] == 1) {
      print PED "0 ";
      print MATERNAL "2 ";
      print PATERNAL "2 ";
    } elsif ($maternal[$i] == 5 && $paternal[$i] == 5) {
      print PED "NA ";
      print MATERNAL "NA ";
      print PATERNAL "NA ";
    }

#    unless($maternal[$i] == 5 || $paternal[$i] == 5) {
#      print PED "$allele{$snps[$i]}{$maternal[$i]} $allele{$snps[$i]}{$paternal[$i]} ";
#      print MATERNAL "$allele{$snps[$i]}{$maternal[$i]} $allele{$snps[$i]}{$maternal[$i]} ";
#      print PATERNAL "$allele{$snps[$i]}{$paternal[$i]} $allele{$snps[$i]}{$paternal[$i]} ";
#    } else {

  }

  print PED "\n";
  print MATERNAL "\n";
  print PATERNAL "\n";
#  print RESOLVED "\n";

  print "Printed individual $famid\n";

}

close PED;
close MATERNAL;
close PATERNAL;
#close RESOLVED;
