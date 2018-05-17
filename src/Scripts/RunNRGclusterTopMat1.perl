#!/usr/bin/perl
##
## Script to start the NRG code in the UT/Winnbago clusters 
##
##using getopt
##
##use Getopt::Std;
##%options=();
##getopts("",\%options);

sub Print_Help(){
  print "Usage $0 -i=d1(--dir0=d1)  -f=d2(--dirF=d2)  (-C,--command=\"Command -m Option\") \n";
  print "Options : \n";
  print "\t--help\t\t: prints this message \n";
  print "\t--addtoname=string \t: adds \"string\" to the output file names extension\n";
  print "\t--queue=queuename \t: sends to queue \"queuename\"\n";
  print "\t--special \t: adds special features in the script (edit it here) \n";
  print "\t--pbsresources=\"nodes=wilson01:ppn=2\" , etc\t: text for -l option in qsub (advanced). Useful to request a single node (default: nodes=1:ppn=2) \n";
  print "\t--CFS\t: DM_NRG -C option: uses Complete Fock Space calculation. \n";
  print "\t--Nw=nw \t: DM_NRG -n nw option: Number of omegas per shell (default: 1) \n";
  print "\t--SubGap\t: DM_NRG -G option: get the sub-gap spectrum \n";
  exit;

}
use Getopt::Long;

my $Nw=1;

my $DMNRGExecName="DM_NRG ";
GetOptions("queue=s"=>\$QueueName,
           "i|dir0=i" =>\$Dir0,
           "f|dirF=i" =>\$DirF,
           "C|command=s" => \$Command,
           "addtoname=s" => \$AddToName,
           "pbsresources=s" => \$PBSresources,
           "Nw=i" => \$Nw,
           "SubGap" => \$SubGap,
           "CFS" => \$CFS,
           "special" => \$Special,
           "h|help" => \$HelpOpt);

Print_Help() if defined $HelpOpt;
Print_Help() if ( (!defined($Dir0))||(!defined($DirF)) );
if (defined($SubGap)){$DMNRGExecName.=" -G";}
if (defined($CFS)){
  $UseCFS=1;  
  $DMNRGExecName.=" -C -n ${Nw}";
  print "Using CFS option: ".$DMNRGExecName."\n";
}


##
## Check hostname
##
                                                                              
   chomp($HostName=`/bin/hostname --long`);
   for ($HostName)
    {
      if (/winnebago/){print "Winnebago cluster \n"; }
      elsif (/correlated/){print "Correlated cluster. \n";}
      elsif (/glyph/){print "Glyph cluster. \n";}
      elsif (/wilson/){print "Wilson cluster. \n";}
      elsif (/aguia000/){print "Aguia cluster (PBS). \n";}
      elsif (/WilsonLaptop1/){print "Laptop (testing). \n";}
      elsif (/TOPMAT1.cs1f5cloud.internal/){print "TOPMAT1 VM. \n";}
      else{print "Cant recognize cluster. Exiting... \n"; exit(0);}
    }
                                                                              

$test1="Running NRG ";
$test1.="in cluster $HostName ";
$test1.="in queue $QueueName " if defined $QueueName;
$test1.="in dirs $Dir0 through $DirF ";
##$test1.=" and avoiding $AvoidNode" if defined $AvoidNode;
$test1.=" NOW! \n";
print $test1;

##
##  Choose code:
##

my $UseAllCodes=0;
if (!defined($Command)){
  $ChoiceCode=3;
  print "Choose code: \n";
  print "      1     - OneChQS \n";
  print "      2     - TwoChQS \n";
  print "      3     - NRG_main (default) \n";
  print "        3.1 - NRG_main with ./DM_NRG and ./Conductance \n";
  print "      4     - DM_NRG \n";
  print "        4.1 - DM_NRG with ./Conductance (allows for G vs T) \n";
  print "      5     - Conductance \n \n";
  print "  choice : ";
##  chomp($ChoiceCode = <STDIN>);
  chomp($test1 = <STDIN>);
  if ( ($test1 ne '')&&($test1>0)&&($test1<=5) ){
    $ChoiceCode=$test1;
  }else{print "Invalid choice. Using default (=$ChoiceCode) \n";}

  for ($ChoiceCode){
   if (/^1/){
    $CodeName="OneChQS";
   }
   elsif (/^2/){
    $CodeName="TwoChQS";
   }
   elsif (/^3/){
    $CodeName="NRG_main";
    if ($ChoiceCode=~ m/^3.1/){
        $ChoiceCode=3;
	$UseAllCodes=1;
	print "Using ALL codes...\n";
##        $UseCFS=0;
##        print "Use Complete Fock Space? (y/n) "; chomp($ChoiceYN = <STDIN>);
##        if ($ChoiceYN eq "y"){$UseCFS=1;}
    }
   }
   elsif (/^4/){
    $CodeName="DM_NRG";
    if ($ChoiceCode=~ m/^4.1/){
        $ChoiceCode=4;
	$UseAllCodes=1;
	print "Adding Conductance code...\n";
    }
   }
   elsif (/^5/){
    $CodeName="Conductance";
   }
   else{print "Code not valid. Exiting... \n";exit; }
  }
  print " --- Code = $CodeName \n";
##
##  Choose model:
##
  if (($ChoiceCode==4)||($ChoiceCode==5)){
    $ChoiceModel=7;
  }
  else{
    $ChoiceModel=0;
    print "Choose model : \n";
    print "      0     - Anderson QS basis (default) \n";
    print "        0.1 - Anderson Q basis (NRG_main only) \n";
    print "        0.2 - Anderson (Q Sz) basis (NRG_main only) \n";
    print "        0.3 - Anderson S basis with SC leads (NRG_main only) \n";
    print "        0.4 - Anderson Sz basis with SC leads (NRG_main only) \n";
    print "      1     - Kondo \n";
    print "      2     - Phonon (1ch) \n";
    print "      3     - Chain \n";
    print "      4     - CM phonon (2ch) \n";
    print "        4.1 - CM phonon w/ parity (QSP basis, NRG_main only) \n";
    print "      5     - SMM \n";
    print "      6     - Double Quantum Dot (DQD) 1chQS \n";
    print "        6.1 - DQD with Zeeman in both dots (1chQSz) (NRG_main only) \n";
    print "      7     - DM_NRG/Conductance calculation (no model) \n";
    print "      8     - Majorana model (Nup Pdn basis)  \n";
    print "      9     - Double Majorana model (Pup Pdn basis)  \n";
    print "  choice : ";
#    chomp($ChoiceModel = <STDIN>);
    chomp($test1 = <STDIN>);
    if ($test1 ne ''){
      $ChoiceModel=$test1;
    }else{print "Using default (=$ChoiceModel) \n";}
  }
  for ($ChoiceModel){
   if (/^0/){
    $ModelName="Anderson";
    if ($ChoiceModel =~ m/^0.1/){$ModelName="1chQ_Anderson";}
    if ($ChoiceModel =~ m/^0.2/){$ModelName="1chQSz_Anderson";}
    if ($ChoiceModel =~ m/^0.3/){$ModelName="1chS_AndersonSC";}
    if ($ChoiceModel =~ m/^0.4/){$ModelName="1chSz_AndersonSC";}
   }
   elsif (/^1/){
    $ModelName="Kondo";
   }
   elsif (/^2/){
    $ModelName="Phonon";
   }
   elsif (/^3/){
    $ModelName="Chain";
   }
   elsif (/^4/){
    $ModelName="CMphonon";
    if ($ChoiceModel =~ m/^4.1/){$ModelName="2chQSP_CMphonon";}
   }
   elsif (/^5/){
    $ModelName="SMM";
   }
   elsif (/^6/){
    $ModelName="1chQS_DQD";
    if ($ChoiceModel =~ m/^6.1/){$ModelName="1chQSz_DQD";}
   }
   elsif (/^7/){
    $ModelName="";
   }
   elsif (/^8/){
    $ModelName="1chNupPdn_Majorana";
   }
   elsif (/^9/){
    $ModelName="1chPupPdn_Majorana";
   }
     else{print "Model not valid. Exiting... \n";exit; }
  }
 
  print " --- Model: $ModelName \n";
##
##  Choose band :
##
  if ($ChoiceCode==4){ 
    $ChoiceBand=0;
  } ## end if ChoiceCode=DM_NRG
  elsif ($ChoiceCode==5){
    $ChoiceBand=0;
  }## ChoiceCode=Conductance 
  else{
    $ChoiceBand=0;
    print "Choose host band type : \n";
    print "      0     - Constant (square) DoS [rho(e)=rho_0 -D < e < D] (default) \n";
    print "      1     - Side Dot (needs lanc.in) \n";
    print "      2     - Const (to use z-trick, for instance) \n";
    print "      3     - PowerLaw (needs lanc.in) \n";
    print "      4     - FromFile (needs HybFunc.dat) \n";
    print "      5     - Cavity (needs lanc.in) \n";
    print "  choice : ";
#    chomp($ChoiceBand = <STDIN>);
    chomp($test1 = <STDIN>);
    if ($test1 ne ''){
      $ChoiceBand=$test1;
    }else{print "Using default (=$ChoiceBand) \n";}
  }
  for ($ChoiceBand){
   if (/^0/){
    $BandName="SquareWilson";
   }
   elsif (/^1/){
    $BandName="SideDot";
   }
   elsif (/^2/){
    $BandName="Const";
   }
   elsif (/^3/){
    $BandName="PowerLaw";
   }
   elsif (/^4/){
    $BandName="FromFile";
   }
   elsif (/^5/){
    $BandName="Cavity";
   }
   else{print "Band not implemented. Using Default \n";
     $ChoiceBand=0;
     $BandName="SquareWilson";
   }
  }
# end for ChoiceBand

##
## Temperature
##
  $Mtemp=1000;
  $betabar=0.727;
  $FiniteT=0;
  $NewBBar=0;
  $LoopMtemp=0;
  print "Mtemp = 1000 (T=0) and betabar=0.727. Change? (y/n) ";
  chomp($ChoiceYN = <STDIN>);
  if ($ChoiceYN eq "y"){
    print "Mtemp = "; chomp($Input = <STDIN> );
    if ( ($Input ne '') && ($Input>0) && ($Input<1000) ){
      $Mtemp=$Input;
      $FiniteT=1;
      print "Mtemp = $Mtemp \n";
##      print "Loop Mtemp in each dir? (Mtemp0=$Mtemp ) (y/n) ";
      print "Loop Mtemp (in same dir)? (Mtemp0=$Mtemp ) (y/n) ";
      chomp($ChoiceLoopMtemp = <STDIN> );
      if ($ChoiceLoopMtemp eq "y"){
        $LoopMtemp=1;
        $Mtemp0=$Mtemp;
        print "MtempFinal = "; chomp($MtempFinal = <STDIN> );
        print "MtempStep = "; chomp($MtempStep = <STDIN> );
      } else {
        $Mtemp0=$Mtemp;$MtempFinal=$Mtemp;$MtempStep=2;
      }## end if LoopMtemp
    } else {print "Not valid. Keeping Mtemp = $Mtemp;"}
    print "(default value: 0.727) betabar = "; chomp($Input = <STDIN> );
    if ( ($Input ne '') && ($Input>0) && ($Input<1.0) ){
      $betabar=$Input;
      $NewBBar=1;
      print "betabar = $betabar \n";
    } else {print "Not valid. Keeping betabar = $betabar\n";}
  }## end if choice=y

##
##  Ztrick: command line 
##
 $z0=0.6;
 $zF=1.41;
 $zstep=0.2;
 print "Use Ztrick (z=".$z0."...".$zF." step ".$zstep." -> Nz=".(int(($zF-$z0)/$zstep)+1).") ? (y/n) :";
 chomp($UseZtrickYN = <STDIN>);
 if ($UseZtrickYN eq "y"){
   $UseZtrick=1;
   if ($BandName eq "SquareWilson"){
     print "z-trick: Changing Band choice to \"Const\"... \n";
     $ChoiceBand=2;$BandName="Const";
   } ## end if BandName=SquareWilson
 }else{$UseZtrick=0;}
##$UseZtrick=0;
}
# end if !defined($Command)
else{
## Strip "./" from Command
 $Command =~ s/\.//;
 $Command =~ s/\///;
 print "Code/Option is $Command \n";
}

##
## Loop on dirs 
##
##
## Executing
##

## Setting Dir-independent variables
   if(defined($Command)){$ExecName=$Command;}
   else{
     $ExecName=$CodeName;
     if ($NewBBar==1){$DMNRGExecName.=" -b ".$betabar;}
     if ($FiniteT==1){
       $DMNRGExecNameTeq0=$DMNRGExecName;
       $DMNRGExecName.=" -M ".$Mtemp;
     }
     if ($ChoiceCode==4) {
       $ExecName=$DMNRGExecName;
     } elsif ($ChoiceCode!=5) {
       if ($ModelName ne "") {
	 $ExecName.=" -m ".$ModelName;
       }
       if ($ChoiceBand != 0) {
	 $ExecName.=" -b ".$BandName;
       }
     }# end if ChoiceCode=4
   }# end if defined($Command)
 
for ($Dir=$Dir0;
    $Dir<=$DirF;
    $Dir+=1){

   $TwoChDir=$ENV{'TWOCHDIR'};
   print "Root Dir: $TwoChDir\n"; 
   $User=$ENV{'USER'};
   $Extension="NRG_Code_".$CodeName;
   if ($UseCFS==1){$Extension.="_CFS";}
   if ($ModelName ne ""){$Extension.="_Model_".$ModelName;}
   if ($ChoiceBand != 0){$Extension.="_Band_".$BandName;}




   $LocalDir=$TwoChDir."/Runs/run$Dir";
   print "Running from : $LocalDir \n";
   $Extension.="_Dir".$Dir;
   if(defined($AddToName)){$Extension.=$AddToName;}

###
### Preparing submission
###

# Build Script file
#
   $ScriptFileName="sub_script".$Extension;
   $String="\#!/bin/sh \n";

   if ($HostName =~ m/winnebago/)
     {
      $ScriptFileName="SLURM_".$ScriptFileName;
      $String.="\#SLURM -s\n\#\n\#SLURM -J job".$Extension."\n\#\n";
      $String.="\#SLURM -e $LocalDir/job".$Extension.".error\n\#\n";
      $String.="\#SLURM -o $LocalDir/job".$Extension.".output\n\#\n";
##      $String.="\#SLURM -D ".$LocalDir."\n\#\n";
##      $SubCommand="srun -b ";
      $SubCommand="sbatch ";
      $SubCommand.="-p $QueueName " if defined $QueueName;
      $SubCommand.=" $LocalDir/$ScriptFileName \n";
     }
   elsif ($HostName =~ m/correlated/)
     {
      $ScriptFileName="LAVA_".$ScriptFileName;
      $String.="\#BSUB -J job".$Extension."\n\#\n";
      $String.="\#BSUB -e job".$Extension.".error\n\#\n";
      $String.="\#BSUB -o job".$Extension.".output\n\#\n";
##   $String.="\#BSUB -D ".$LocalDir."\n\#\n";
      $SubCommand="bsub ";
      $SubCommand.="-q $QueueName " if defined $QueueName;
      $SubCommand.=" < $LocalDir/$ScriptFileName \n";
     }
   elsif ( ($HostName =~ m/glyph/)||($HostName =~ m/wilson/)||($HostName =~ m/aguia000/)||($HostName =~ m/WilsonLaptop1/)||($HostName =~ m/TOPMAT1.cs1f5cloud.internal/) )
     {
      $ScriptFileName="PBS_".$ScriptFileName;
      $String.="\#PBS -N job".$Extension."\n\#\n";
      $String.="\#PBS -e $LocalDir/job".$Extension.".error\n\#\n";
      $String.="\#PBS -o $LocalDir/job".$Extension.".output\n\#\n";
      if (defined($PBSresources)){$String.="\#PBS -l $PBSresources \n";}
      else{$String.="\#PBS -l nodes=1:ppn=2 \n";}
      $SubCommand="qsub ";
      $SubCommand.="-q $QueueName " if defined $QueueName;
      $SubCommand.=" $LocalDir/$ScriptFileName \n";
## ## -d path : Defines the working directory path and stores it in $PBS_O_INITDIR
###      $SubCommand.=" $ScriptFileName -d $LocalDir \n";
     }
   else {print "Which cluster is this? \n"; exit(0);}
   $String.="cd ".$LocalDir."\n";
## Adding z-trick
   if ($UseZtrick==1){
######
## Using Z-trick
######
##     $z0=0.6;
##     $zF=1.41;
##     $zstep=0.2;
     for ($zz=$z0;$zz<=$zF;$zz+=$zstep){
       $ExecNameZ=$ExecName." -z ".$zz;
       $ExtensionZ=$Extension."_zEQ".$zz;
       $String.="nice ./$ExecNameZ > output_$ExtensionZ.txt \n";
       ##$String.="sed '3 c21' ThisCodePars.dat > a.dat \n";
       ##$String.="cat a.dat > ThisCodePars.dat \n";
       if (defined($Special)){
         $String.="perl -I../ ../ChangeNRGParams.perl -i $Dir -f $Dir --choiceparam=9 --p0=2 --pstep=0 --indir  \n";
         $String.="nice ./$ExecNameZ > output_$ExtensionZ.txt \n";
         $String.="perl -I../ ../ChangeNRGParams.perl -i $Dir -f $Dir --choiceparam=9 --p0=3 --pstep=0 --indir  \n";
       }
## end special
       if ($UseAllCodes==1){
	 if ($FiniteT!=1) { ## T=0: onde file for each Z.
	     if ($UseCFS==1) {
	     $String.="sed '3 c\21 ' ThisCodePars.dat > a.dat \n";
	     $String.="cat a.dat > ThisCodePars.dat \n";
	     $String.="nice ./".$DMNRGExecName."  > output_DMNRG_CFS_zEQ${zz}_Dir${Dir}.txt \n";
	   } else {
	     $String.="nice ./".$DMNRGExecName." -n ${Nw} > output_DMNRG_zEQ${zz}_Dir${Dir}.txt \n";
	   }			## end if UseCFS==1
        ## Finite-T: loop over Mtemp for each Z.
        } else { ## else $FiniteT==1
           for ($Mtemp=$Mtemp0;$Mtemp<=$MtempFinal;$Mtemp+=$MtempStep){
             ##$DMNRGExecName=$DMNRGExecNameTeq0." -M ".$Mtemp;  ## Need to add -n Nw !!!
             $DMNRGExecName=$DMNRGExecNameTeq0." -M ".$Mtemp." -n ".$Nw;  ## Need to add a CFS flag for finite-T here!
             $ExtensionMtemp=$Extension."_MtempEQ".$Mtemp."_zEQ".$zz;
             $String.="nice ./$DMNRGExecName > output_$ExtensionMtemp.txt \n";
## No conductance in this mode. I would need to add another Mtemp loop;
##             $String.="nice ./Conductance > output_Conductance";
##             $String.="nice ./Conductance -R -G > output_Resistivity";
##             if ($UseCFS==1){$String.="_CFS";}
##             $String.="_Mtemp".$Mtemp."_zEQ".$zz."_Dir".$Dir.".txt \n";
             $String.="ls rho_\*_OmegaRhow\* \| sed \"s\/rho_\\\(\.\*\\)\.dat\/mv \& rhoMtemp".$Mtemp."_\\1\.dat/\" \| sh \n";
             $String.="rm -f rhoDM\*.bin \n";
           }
           # end loop in Mtemp
           $Mtemp=$Mtemp0;
        } ## end if FInite_T==1

	## Clean up NRG bin files
        $String.="rm -f \*.bin \n";
       }
       # end if UseAllCodes=1 (using Z-trick with DM-NRG)
     }
     # end loop in Z

     if ( ($UseAllCodes==1)&&($ChoiceCode != 5) ){
       if ($FiniteT != 1) {
	 $String.="ls rho_\*_OmegaRhow.dat \| sed \"s\/rho_\\(\.\*\\)_OmegaRhow\.dat\/cp \-p \& rho_\\1_OmegaRhow_zEQ1\.00\.dat/\" \| sh \n";
	 ## Add a script here to average rho files
	 $String.="./Zavg.perl --type=\"rho_0_0_OmegaRhow\" --noask \n";
       } else {
	 for ($Mtemp=$Mtemp0;$Mtemp<=$MtempFinal;$Mtemp+=$MtempStep){
	   $String.="ls rhoMtemp".$Mtemp."_\*_OmegaRhow.dat \| sed \"s\/rhoMtemp".$Mtemp."_\\(\.\*\\)_OmegaRhow\.dat\/cp \-p \& rhoMtemp".$Mtemp."_\\1_OmegaRhow_zEQ1\.00\.dat/\" \| sh \n";
	 ## Add a script here to average rho files and copy to usual rhoAVG
	   $String.="./Zavg.perl --type=\"rhoMtemp".$Mtemp."_0_0_OmegaRhow\" --noask \n";
	   $String.="cp -p rhoMtemp".$Mtemp."_0_0_OmegaRhow_zEQAVG.00.dat rho_0_0_OmegaRhow_zEQAVG.00.dat \n";
         } ## end loop in Mtemp
       } ##end if Finite-T==1
## No conductance in this mode. I would need to add another Mtemp loop;
##       $String.="nice ./Conductance > output_Conductance_dir$Dir.txt \n";
     } ## end if UseAllCodes after Ztrick loop

     if (defined($Special)){
       $String.="./Zavg.perl --noask \n";
       $String.="ln -s SuscepImp1ChQS_Anderson_zEQAVG.dat SuscepImp_25_727.dat\n";
     }
     ## end if Special after Ztrick loop
   } else{ ## end if Z-trick
######
## no Z-trick
######
    print " ExecName is ".$ExecName."\n";
    $String.="nice ./$ExecName > output_$Extension.txt \n";
## Add stuff to script
     if (defined($Special)){
      $String.="perl -I../ ../ChangeNRGParams.perl -i $Dir -f $Dir --choiceparam=9 --p0=2 --pstep=0 --indir  \n";
      $String.="nice ./$ExecName > output_$Extension.txt \n";
##      $String.="../GrepNRGOpAvg.perl --indir --chiloc \n"; ## ChiLoc (not Suscep)
##      $String.="cp -p SuscepImp1ChQS_Anderson_zEQ1.dat SuscepImp_25_727.dat \n";
##      $String.="./CalcTK 2 $Dir \n";
##      $String.="../GrepNRGOpAvg.perl --indir --grep=\"Ndot\" -p 4 --Tzero -i $Dir -f $Dir \n";
     }
##   end special
     print "UseAllCodes = $UseAllCodes \n" ;
     if ($UseAllCodes==1){
        print "Adding DM_NRG and/or Conductance \n";
        if ( ($ChoiceCode==4)&&($FiniteT==1) ){
##          $String.="nice ./Conductance > output_Conductance"; 
          $String.="nice ./Conductance -R -G > output_Resistivity";
          if ($UseCFS==1){$String.="_CFS";}
          $String.="_Mtemp".$Mtemp0."_Dir".$Dir.".txt \n";
          $String.="ls rho_\*_OmegaRhow\* \| sed \"s\/rho_\\\(\.\*\\)\.dat\/mv \& rhoMtemp".$Mtemp."_\\1\.dat/\" \| sh \n";
          $String.="rm -f rhoDM\*.bin \n";
          if($LoopMtemp==1){
           for ($Mtemp=$Mtemp0+$MtempStep;$Mtemp<=$MtempFinal;$Mtemp+=$MtempStep){
             $DMNRGExecName=$DMNRGExecNameTeq0." -M ".$Mtemp;
             $ExtensionMtemp=$Extension."_MtempEQ".$Mtemp;
             $String.="nice ./$DMNRGExecName > output_$ExtensionMtemp.txt \n";
##             $String.="nice ./Conductance > output_Conductance";
             $String.="nice ./Conductance -R -G > output_Resistivity";
             if ($UseCFS==1){$String.="_CFS";}
             $String.="_Mtemp".$Mtemp."_Dir".$Dir.".txt \n";
             $String.="ls rho_\*_OmegaRhow\* \| sed \"s\/rho_\\\(\.\*\\)\.dat\/mv \& rhoMtemp".$Mtemp."_\\1\.dat/\" \| sh \n";
             $String.="rm -f rhoDM\*.bin \n"; 
           } 
           # end loop in Mtemp
           $Mtemp=$Mtemp0;
         } 
         #end if LoopMtemp==1
        } else {
          if($ChoiceCode != 4){ 
           if ($UseCFS==1){$String.="nice ./".$DMNRGExecName." > output_DMNRG_CFS_Dir$Dir.txt \n";}
           else {$String.="nice ./".$DMNRGExecName." -n ${Nw} > output_DMNRG_Dir$Dir.txt \n";}
          } 
          if($ChoiceCode != 5){ $String.="nice ./Conductance > output_Conductance_dir$Dir.txt \n";}
        }
     }
##   end if UseAllCodes==1
   }
## end if UseZtrick

  $String.="echo Finished > $LocalDir/Finish_$Extension \n";
  print "Building $ScriptFileName ... ";

  open (OUT,"> $LocalDir/".$ScriptFileName);
    print OUT $String;
  close(OUT);

  print "... Done.\n";
#
# Run PB/SLURM...
#
  print "Executing $SubCommand ...\n";
  system($SubCommand);
  print "... Done!\n";

}
## End Loop in Dirs
