#!/usr/bin/perl
#This program was created by Bo Segerman, National veterinary Institute, Uppsala, Sweden, bo.segerman@sva.se
######################################################################################################
######################################################################################################
######################################################################################################


#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
#####################################  COMPONENT: GLOBAL INITIATION ##################################

# version 2022.09.16

use IO::Zlib;
use Cwd;

$global_version = "2022.09.16";


############################### END COMPONENT: GLOBAL INITIATION  ####################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================





#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\






#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
######################### COMPONENT: DECLARATION OF TOOLS ############################################

# version 2022.03.04

%global_toolhash = ();
@global_tooldesc = ();

$global_toolhash{'version'} = \&sub_version; 
$global_toolhash{'-version'} = \&sub_version; 
$global_toolhash{'--version'} = \&sub_version; 
$global_toolhash{'-v'} = \&sub_version; 
push(@global_tooldesc,"version (--verbose)");


$global_toolhash{'help'} = \&sub_help; 
$global_toolhash{'-help'} = \&sub_help; 
$global_toolhash{'--help'} = \&sub_help; 
push(@global_tooldesc,"help");


$global_toolhash{'fileinfo'} = \&sub_fileinfo; 
push(@global_tooldesc,"fileinfo filename (--verbose --tabular --header)");

$global_toolhash{'downsample'} = \&sub_downsample; 
push(@global_tooldesc,"downsample R1.fastq.gz (R2.fastq.gz) coverage=xx genomesize=xx|organism=xx OR(bases=xx)");

$global_toolhash{'filtercontigs'} = \&sub_filtercontigs; 
push(@global_tooldesc,"filtercontigs infile.fasta (coveragethreshold=x% sizethreshold=x)");

$global_toolhash{'kmeroverlapp'} = \&sub_kmeroverlapp; 
push(@global_tooldesc,"kmeroverlapp assembly1.fasta assembly2.fasta (kmersize=x)");

$global_toolhash{'kmerhisto'} = \&sub_kmerhisto; 
push(@global_tooldesc,"kmerhisto R1.fastq.gz (R2.fasta.gz) (reference=reference.fasta)");



############################### END COMPONENT: DECLARATION OF TOOLS  #################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================





#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\






#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
######################### COMPONENT: CONFIGURATION DATA ##############################################

# version 2022.09.16


#### FIRST LOAD BUILT-IN DEFAULT VALUES 

#################
### settings

$setting_filereadbuffersize = 1000000;   ## size of chunks in bytes that are read by filereaders .... 
 					  ## small size may affect speed negatively, large size increases RAM usage
$setting_fastqinitialanalysis = 100000;
$setting_16staglength = 40;
$setting_16s_threshold_default = 0.1;

#################################  Reference sequences ##############################

@constants_species = ();
@constants_species_16Sseq = ();
@constants_species_reference = ();
@constants_species_16Sthreshold = ();

#Campylobacter jejuni 
push(@constants_species,"Campylobacter jejuni/coli");
push(@constants_species_16Sseq,"TTTTTATGGAGAGTTTGATCCTGGCTCAGAGTGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAACGATGAAGCTTTTAGCTTGCTAGAAGTGGATTAGTGGCGCACGGGTGAGTAAGGTATAGTTAATCTGCCCTACACAAGAGGACAACAGTTGGAAACGACTGCTAATACTCTATACTCCTGCTTAACACAAGTTGAGTAGGGAAAGTTTTTCGGTGTAGGATGAGACTATATAGTATCAGCTAGTTGGTAAGGTAATGGCTTACCAAGGCTATGACGCTTAACTGGTCTGAGAGGATGATCAGTCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGAATATTGCGCAATGGGGGAAACCCTGACGCAGCAACGCCGCGTGGAGGATGACACTTTTCGGAGCGTAAACTCCTTTTCTTAGGGAAGAATTCTGACGGTACCTAAGGAATAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTACTCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGATTATCAAGTCTCTTGTGAAATCTAATGGCTTAACCATTAAACTGCTTGGGAAACTGATAGTCTAGAGTGAGGGAGAGGCAGATGGAATTGGTGGTGTAGGGGTAAAATCCGTAGATATCACCAAGAATACCCATTGCGAAGGCGATCTGCTGGAACTCAACTGACGCTAAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTACACTAGTTGTTGGGGTGCTAGTCATCTCAGTAATGCAGCTAACGCATTAAGTGTACCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATAGACGGGGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGATACGCGAAGAACCTTACCTGGGCTTGATATCCTAAGAACCTTATAGAGATATGAGGGTGCTAGCTTGCTAGAACTTAGAGACAGGTGCTGCACGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCACGTATTTAGTTGCTAACGGTTCGGCCGAGCACTCTAAATAGACTGCCTTCGTAAGGAGGAGGAAGGTGTGGACGACGTCAAGTCATCATGGCCCTTATGCCCAGGGCGACACACGTGCTACAATGGCATATACAATGAGACGCAATACCGCGAGGTGGAGCAAATCTATAAAATATGTCCCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGCCGGAATCGCTAGTAATCGTAGATCAGCCATGCTACGGTGAATACGTTCCCGGGTCTTGTACTCACCGCCCGTCACACCATGGGAGTTGATTTCACTCGAAGCCGGAATACTAAACTAGTTACCGTCCACAGTGGAATCAGCGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGAGAACCTGCGGTTGGATCACCTCCT");
push(@constants_species_reference, "GCF_000009085.1");
push(@constants_species_16Sthreshold,0.1);

## Campylobacter coli and C. jejuni have nearly identical 16S
##Campylobacter coli 
#push(@constants_species,"Campylobacter coli");
#push(@constants_species_16Sseq,"TTTATGGAGAGTTTGATCCTGGCTCAGAGTGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAACGATGAAGCTTCTAGCTTGCTAGAAGTGGATTAGTGGCGCACGGGTGAGTAAGGTATAGTTAATCTGCCCTACACAAGAGGACAACAGTTGGAAACGACTGCTAATACTCTATACTCCTGCTTAACACAAGTTGAGTAGGGAAAGTTTTTCGGTGTAGGATGAGACTATATAGTATCAGCTAGTTGGTAAGGTAATGGCTTACCAAGGCTATGACGCTTAACTGGTCTGAGAGGATGATCAGTCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGAATATTGCGCAATGGGGGAAACCCTGACGCAGCAACGCCGCGTGGAGGATGACACTTTTCGGAGCGTAAACTCCTTTTCTTAGGGAAGAATTCTGACGGTACCTAAGGAATAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTACTCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGATTATCAAGTCTCTTGTGAAATCTAATGGCTTAACCATTAAACTGCTTGGGAAACTGATAGTCTAGAGTGAGGGAGAGGCAGATGGAATTGGTGGTGTAGGGGTAAAATCCGTAGATATCACCAAGAATACCCATTGCGAAGGCGATCTGCTGGAACTCAACTGACGCTAAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTACACTAGTTGTTGGGGTGCTAGTCATCTCAGTAATGCAGCTAACGCATTAAGTGTACCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATAGACGGGGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGATACGCGAAGAACCTTACCTGGGCTTGATATCCTAAGAACCTTTTAGAGATAAGAGGGTGCTAGCTTGCTAGAACTTAGAGACAGGTGCTGCACGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCACGTATTTAGTTGCTAACGGTTCGGCCGAGCACTCTAAATAGACTGCCTTCGTAAGGAGGAGGAAGGTGTGGACGACGTCAAGTCATCATGGCCCTTATGCCCAGGGCGACACACGTGCTACAATGGCATATACAATGAGACGCAATACCGCGAGGTGGAGCAAATCTATAAAATATGTCCCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGCCGGAATCGCTAGTAATCGTAGATCAGCCATGCTACGGTGAATACGTTCCCGGGTCTTGTACTCACCGCCCGTCACACCATGGGAGTTGATTTCACTCGAAGCCGGAATACTAAACTAGTTACCGTCCACAGTGGAATCAGCGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGAGAACCTGCGGTTGGATCACCTCCTTT");
#push(@constants_species_reference, "GCF_009730395.1");
#push(@constants_species_16Sthreshold,0.1);

#Salmonella enterica
push(@constants_species,"Salmonella enterica"); 
push(@constants_species_16Sseq,"TTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTGGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCTACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTAGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTTAGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTACCTTAAAGAAGC");
push(@constants_species_reference, "GCF_000006945.2");
push(@constants_species_16Sthreshold,0.1);

#Listeria monocytogenes
push(@constants_species,"Listeria monocytogenes"); 
push(@constants_species_16Sseq,"TATTTTAAAGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAACGAACGGAGGAAGAGCTTGCTCTTCCAATGTTAGTGGCGGACGGGTGAGTAACACGTGGGCAACCTGCCTGTAAGTTGGGGATAACTCCGGGAAACCGGGGCTAATACCGAATGATAAGATGTGGCGCATGCCACGCCTTTGAAAGATGGTTTCGGCTATCGCTTACAGATGGGCCCGCGGTGCATTAGCTAGTTGGTAGGGTAATGGCCTACCAAGGCAACGATGCATAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGTATGAAGAAGGTTTTCGGATCGTAAAGTACTGTTGTTAGAGAAGAACAAGGATAAGAGTAACTGCTTGTCCCTTGACGGTATCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGCGCGCAGGCGGTCTTTTAAGTCTGATGTGAAAGCCCCCGGCTTAACCGGGGAGGGTCATTGGAAACTGGAAGACTGGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGATATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGACCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTTTGACCACTCTGGAGACAGAGCTTTCCCTTCGGGGACAAAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGATTTTAGTTGCCAGCATTTAGTTGGGCACTCTAAAGTGACTGCCGGTGCAAGCCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGATAGTACAAAGGGTCGCGAAGCCGCGAGGTGGAGCTAATCCCATAAAACTATTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGCCGGAATCGCTAGTAATCGTGGATCAGCATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTAGGGTAACCTTTATGGAGCCAGCCGCCGAAGGTGGGACAGATAATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTCT");
push(@constants_species_reference, "GCF_000196035.1");
push(@constants_species_16Sthreshold,0.1);

#Escherichia coli
push(@constants_species,"Escherichia coli"); 
push(@constants_species_16Sseq,"TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGAGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT");
push(@constants_species_reference, "GCF_000008865.2");
push(@constants_species_16Sthreshold,0.1);

#Mycobacterium tuberculosis
push(@constants_species,"Mycobacterium tuberculosis"); 
push(@constants_species_16Sseq,"TTTTGTTTGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGGTCTCTTCGGAGATACTCGAGTGGCGAACGGGTGAGTAACACGTGGGTGATCTGCCCTGCACTTCGGGATAAGCCTGGGAAACTGGGTCTAATACCGGATAGGACCACGGGATGCATGTCTTGTGGTGGAAAGCGCTTTAGCGGTGTGGGATGAGCCCGCGGCCTATCAGCTTGTTGGTGGGGTGACGGCCTACCAAGGCGACGACGGGTAGCCGGCCTGAGAGGGTGTCCGGCCACACTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACCTCTTTCACCATCGACGAAGGTCCGGGTTCTCTCGGATTGACGGTAGGTGGAGAAGAAGCACCGGCCAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCGCGTTGTTCGTGAAATCTCACGGCTTAACTGTGAGCGTGCGGGCGATACGGGCAGACTAGAGTACTGCAGGGGAGACTGGAATTCCTGGTGTAGCGGTGGAATGCGCAGATATCAGGAGGAACACCGGTGGCGAAGGCGGGTCTCTGGGCAGTAACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGGTGGGTACTAGGTGTGGGTTTCCTTCCTTGGGATCCGTGCCGTAGCTAACGCATTAAGTACCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGCGGAGCATGTGGATTAATTCGATGCAACGCGAAGAACCTTACCTGGGTTTGACATGCACAGGACGCGTCTAGAGATAGGCGTTCCCTTGTGGCCTGTGTGCAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCTCATGTTGCCAGCACGTAATGGTGGGGACTCGTGAGAGACTGCCGGGGTCAACTCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGCCCCTTATGTCCAGGGCTTCACACATGCTACAATGGCCGGTACAAAGGGCTGCGATGCCGCGAGGTTAAGCGAATCCTTAAAAGCCGGTCTCAGTTCGGATCGGGGTCTGCAACTCGACCCCGTGAAGTCGGAGTCGCTAGTAATCGCAGATCAGCAACGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGTCATGAAAGTCGGTAACACCCGAAGCCAGTGGCCTAACCCTCGGGAGGGAGCTGTCGAAGGTGGGATCGGCGATTGGGACGAAGTCGTAACAAGGTAGCCGTACCGGAAGGTGCGGCTGGATCACCTCCTTTCT");
push(@constants_species_reference, "GCF_000195955.2");
push(@constants_species_16Sthreshold,0.1);

#Klebsiella pneumoniae
push(@constants_species,"Klebsiella pneumoniae");
push(@constants_species_16Sseq,"AGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGTAGCACAGAGAGCTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGTGGGGGACCTTCGGGCCTCATGCCATCAGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGCGGGGAGGAAGGCGTAAGGTTAATAACCTCTCGATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGATTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAATCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCATATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTATGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTAGATCAGAATGCTACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCT");
push(@constants_species_reference, "GCF_000240185.1");
push(@constants_species_16Sthreshold,0.1);

#Streptococcus pneumoniae
push(@constants_species,"Streptococcus pneumoniae");
push(@constants_species_16Sseq,"TTTAATGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTAGAACGCTGAAGGAGGAGCTTGCTTCTCTGGATGAGTTGCGAACGGGTGAGTAACGCGTAGGTAACCTGCCTGGTAGCGGGGGATAACTATTGGAAACGATAGCTAATACCGCATAAGAGTGGATGTTGCATGACATTTGCTTAAAAGGTGCACTTGCATCACTACCAGATGGACCTGCGTTGTATTAGCTAGTTGGTGGGGTAACGGCTCACCAAGGCGACGATACATAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTAGGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGAGTGCTAGGTGTTAGACCCTTTCCGGGGTTTAGTGCCGTAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGACCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCTCTGACCGCTCTAGAGATAGAGCTTTCCTTCGGGACAGAGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTATTGTTAGTTGCCATCATTTAGTTGGGCACTCTAGCGAGACTGCCGGTAATAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGCTGGTACAACGAGTCGCAAGCCGGTGACGGCAAGCTAATCTCTTAAAGCCAGTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGTCGGAATCGCTAGTAATCGCGGATCAGCACGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCGTAAGGAGCCAGCCGCCTAAGGTGGGATAGATGATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTT");
push(@constants_species_reference, "GCF_002076835.1");
push(@constants_species_16Sthreshold,0.1);

#Staphylococcus aureus
push(@constants_species,"Staphylococcus aureus");
push(@constants_species_16Sseq,"TTTTATGGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAACGGACGAGAAGCTTGCTTCTCTGATGTTAGCGGCGGACGGGTGAGTAACACGTGGATAACCTACCTATAAGACTGGGATAACTTCGGGAAACCGTAGCTAATACCGGATAATATTTTGAACCGCATGGTTCAAAAGTGAAAGACGGTCTTGCTGTCACTTATAGATGGATCCGCGCTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAACGATGCATAGCCGACCTGAGAGGGTGATCGGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACATATGTGTAAGTAACTGTGCACATCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAAAACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGCAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGACCGCAAGGTTGAAACTCAAAGGAATTGACGGGGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTTTGACAACTCTAGAGATAGAGCCTTCCCCTTCGGGGGACAAAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTAAGCTTAGTTGCCATCATTAAGTTGGGCACTCTAAGTTGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGATTTGGGCTACACACGTGCTACAATGGACAATACAAAGGGCAGCGAAACCGCGAGGTCAAGCAAATCCCATAAAGTTGTTCTCAGTTCGGATTGTAGTCTGCAACTCGACTACATGAAGCTGGAATCGCTAGTAATCGTAGATCAGCATGCTACGGTGAATACGTTCCCGGGTCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGCCGGTGGAGTAACCTTTTAGGAGCTAGCCGTCGAAGGTGGGACAAATGATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTCT");
push(@constants_species_reference, "GCF_000013425.1");
push(@constants_species_16Sthreshold,0.1);

#######################
### reverse complements

@constants_species_16Sseq_rc = ();

for($i=0;$i<@constants_species_16Sseq;$i++){
  push(@constants_species_16Sseq_rc, &private_reverscomplement($constants_species_16Sseq[$i]));
}

#######################
### typical genomesizes

%constants_typicalgenomesize = ();

$constants_typicalgenomesize{"campylobacter_jejuni"}       = 1700000;
$constants_typicalgenomesize{"campylobacter_jejuni/coli"}  = 1700000;
$constants_typicalgenomesize{"campylobacter_coli"}         = 1700000;
$constants_typicalgenomesize{"campylobacter"}              = 1700000;
$constants_typicalgenomesize{"salmonella_enterica"}        = 5000000;
$constants_typicalgenomesize{"salmonella"}                 = 5000000;
$constants_typicalgenomesize{"listeria_monocytogenes"}     = 3000000;
$constants_typicalgenomesize{"listeria"}                   = 3000000;
$constants_typicalgenomesize{"escherichia_coli"}           = 5000000;
$constants_typicalgenomesize{"ecoli"}                      = 5000000;
$constants_typicalgenomesize{"mycobacterium_tuberculosis"} = 4400000;
$constants_typicalgenomesize{"klebsiella_pneumoniae"}      = 5700000;
$constants_typicalgenomesize{"streptococcus_pneumoniae"}   = 2150000;
$constants_typicalgenomesize{"staphylococcus_aureus"}      = 2800000;

############################  Adapter sequences #############################

@constants_adapters = ();
@constants_adapters_seq = ();


push(@constants_adapters,"Illumina Tagmentation: CTGTCTCTTATACACATCT");
push(@constants_adapters_seq,"CTGTCTCTTATACACATCT");

push(@constants_adapters,"Illumina TruSeq R1: AGATCGGAAGAGCACACGTCTGAAC");
push(@constants_adapters_seq,"AGATCGGAAGAGCACACGTCTGAAC");

push(@constants_adapters,"Illumina TruSeq R2: AGATCGGAAGAGCGTCGTGTAGGGA");
push(@constants_adapters_seq,"AGATCGGAAGAGCGTCGTGTAGGGA");


push(@constants_adapters,"Illumina Tagmentation alt: ATGTGTATAAGAGACA");
push(@constants_adapters_seq,"ATGTGTATAAGAGACA");



############################ extensions ####################################


@constants_R1_extensions = ();
@constants_R2_extensions = ();

push(@constants_R1_extensions,"_1.fastq");
push(@constants_R2_extensions,"_2.fastq");

push(@constants_R1_extensions,"_1.fastq.gz");
push(@constants_R2_extensions,"_2.fastq.gz");

push(@constants_R1_extensions,"_1.fq");
push(@constants_R2_extensions,"_2.fq");

push(@constants_R1_extensions,"_1.fq.gz");
push(@constants_R2_extensions,"_2.fq.gz");

push(@constants_R1_extensions,"_r1.fastq");
push(@constants_R2_extensions,"_r2.fastq");

push(@constants_R1_extensions,"_r1.fastq.gz");
push(@constants_R2_extensions,"_r2.fastq.gz");

push(@constants_R1_extensions,"_r1.t.fastq.gz");
push(@constants_R2_extensions,"_r2.t.fastq.gz");

push(@constants_R1_extensions,"_r1.fq");
push(@constants_R2_extensions,"_r2.fq");

push(@constants_R1_extensions,"_r1.fq.gz");
push(@constants_R2_extensions,"_r2.fq.gz");

push(@constants_R1_extensions,"_r1_001.fastq");
push(@constants_R2_extensions,"_r2_001.fastq");

push(@constants_R1_extensions,"_r1_001.fastq.gz");
push(@constants_R2_extensions,"_r2_001.fastq.gz");

push(@constants_R1_extensions,"_r1_001.fq");
push(@constants_R2_extensions,"_r2_001.fq");

push(@constants_R1_extensions,"_r1_001.fq.gz");
push(@constants_R2_extensions,"_r2_001.fq.gz");

push(@constants_R1_extensions,"-1.fastq");
push(@constants_R2_extensions,"-2.fastq");

push(@constants_R1_extensions,"-1.fastq.gz");
push(@constants_R2_extensions,"-2.fastq.gz");

push(@constants_R1_extensions,"-1.fq");
push(@constants_R2_extensions,"-2.fq");

push(@constants_R1_extensions,"-1.fq.gz");
push(@constants_R2_extensions,"-2.fq.gz");

push(@constants_R1_extensions,"-r1.fastq");
push(@constants_R2_extensions,"-r2.fastq");

push(@constants_R1_extensions,"-r1.fastq.gz");
push(@constants_R2_extensions,"-r2.fastq.gz");

push(@constants_R1_extensions,"-r1.t.fastq.gz");
push(@constants_R2_extensions,"-r2.t.fastq.gz");

push(@constants_R1_extensions,"-r1.fq");
push(@constants_R2_extensions,"-r2.fq");

push(@constants_R1_extensions,"-r1.fq.gz");
push(@constants_R2_extensions,"-r2.fq.gz");

push(@constants_R1_extensions,"-r1_001.fastq");
push(@constants_R2_extensions,"-r2_001.fastq");

push(@constants_R1_extensions,"-r1_001.fastq.gz");
push(@constants_R2_extensions,"-r2_001.fastq.gz");

push(@constants_R1_extensions,"-r1_001.fq");
push(@constants_R2_extensions,"-r2_001.fq");

push(@constants_R1_extensions,"-r1_001.fq.gz");
push(@constants_R2_extensions,"-r2_001.fq.gz");



###################### Overide with data from config file  #########################

$program_dir = substr($0,0,rindex($0,"/")) . "/";

$path_configfile = $program_dir . "tngs.config";

if(-f $path_configfile){
  #NOT IMPLEMENTED YET
}

############################### END COMPONENT: CONFIGURATION DATA  #################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================






#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\








#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
#####################################  COMPONENT: ARGUMENT HANDLING ##################################

# version = 2022.03.04

$tool= lc(shift(@ARGV));   #tool to run is always first argument

%arg_parameters = ();
%arg_flags = ();
@arg_remaining = ();

%arg_flags_accessed = ();
%arg_parameters_accessed = ();
@arg_remaining_accessed = ();
 
if(@ARGV > 0){
 do{
   $arg= shift(@ARGV); 

   if(substr($arg,0,2) eq "--"){
      $arg = lc(substr($arg,2));
      $arg_flags{$arg} = 1;
      $arg_flags_accessed{$arg} = 0;
   }
   elsif(index($arg,"=") >=0 ){
      @columns = split(/=/,$arg);
      $arg_parameters{lc(@columns[0])} = @columns[1];
      $arg_parameters_accessed{lc(@columns[0])} = 0;
   }
   elsif(length($arg) > 0){
     push(@arg_remaining,$arg);
     push(@arg_remaining_accessed,0);
   }

 }while($arg);
}

############################################################################################

# version = 2022.03.04

sub arg_hasparameter(){

  my $param = lc(shift(@_));
  
  if(exists( $arg_parameters{$param}) ) { 
     return 1;
  }
  else{
    return 0;
  }

  
}

#############################################################################################

# version = 2022.03.04

sub arg_getparameter(){

  my $param = lc(shift(@_));

  if(exists( $arg_parameters{$param}) ) { 
     $arg_parameters_accessed{$param} = 1;
     return $arg_parameters{$param};
  }
  else{
    return 0;
  }

  
}

#############################################################################################

# version = 2022.03.04

sub arg_hasflag(){
  my $flag = lc(shift(@_));
  
  if(substr($flag,0,2) eq "--"){
      $flag = substr($flag,2);
  }
  
  if(exists( $arg_flags{$flag}) ) { 
     $arg_flags_accessed{$flag} = 1;
     return 1;
  }
  else{
    return 0;
  }

}

#############################################################################################

# version = 2022.03.04

sub arg_getremaining(){

    my $remaining_nr = shift(@_);



    $arg_remaining_accessed[$remaining_nr]=1;
    
    return $arg_remaining[$remaining_nr];
}

#############################################################################################

# version = 2022.03.04

sub arg_warn_if_unused_argument(){

 foreach $key (keys %arg_flags){
    if( $arg_flags_accessed{$key} == 0) {
       print STDERR "WARNING: unused flag --$key\n";
    }
 }
 foreach $key (keys %arg_parameters){
    $value = $arg_parameters{$key};
    if( $arg_parameters_accessed{$key} == 0) {
       print STDERR "WARNING: unused parameter $key=$value\n";
    }
 } 
 for(my $i=0;$i<@arg_remaining;$i++){
    $value = $arg_remaining[$i];
    if($arg_remaining_accessed[$i]==0){
       print STDERR "WARNING: unused argument $value\n";
    }
 }
}



############################### END COMPONENT: ARGUMENT HANDLING  ####################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================





#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\






#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
##################################### COMPONENT: RUN TOOL ############################################

# version = 2022.03.04


if(exists( $global_toolhash{$tool}) ) { 

  $global_toolhash{$tool}->();

}
else{

  &sub_help();

}

############################### END COMPONENT: RUN TOOL ##############################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================





#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\







#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
####################################### sub_version ##################################################

# version = 2022.03.04

sub sub_version(){

 my $sub_name = "version";
 my $sub_version = "2022.03.04";
 
 my $subfunction_arg = shift(@_);
 if($subfunction_arg eq "version") { return; } # prevent recursive loop   
 
 print "tngs.pl version $global_version\n";
 
 if( &arg_hasflag("verbose") ){
   print "\n";
   print "sub $sub_name version=$sub_version\n";
   foreach $key (keys %global_toolhash)
   {
      if(substr($key,0,1) ne "-"){
         $global_toolhash{$key}->("version");
      }
   }
 }

}

######################################### END: sub_version ###########################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================





#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\







#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
####################################### sub_help #####################################################

# version = 2022.03.04

sub sub_help(){

 my $sub_name = "help";
 my $sub_version = "2022.03.04";
 
 my $subfunction_arg = shift(@_);
 if($subfunction_arg eq "version") { print "sub $sub_name version=$sub_version\n"; return; } 

  
 print   "tngs.pl version $global_version\n\n";
 print   "Avaliable tools/commands:\n";

       foreach my $description (@global_tooldesc){
         print "$description\n";
       }
} 

######################################### END: sub_help ##############################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================






#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\






#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
######################################### sub_fileinfo ###############################################

# version = 2022.09.16

sub sub_fileinfo(){

 my $sub_name = "fileinfo";
 my $sub_version = "2022.09.16";
 
 my $subfunction_arg = shift(@_);
 if($subfunction_arg eq "version") { 
     print "sub $sub_name:\n"; 
     &sub_fileinfo_fastq("version");
     &sub_fileinfo_fasta("version");
     return; 

 } 


  my $fasta = 0; 
   if( &arg_hasflag("fasta") ){
     $fasta = 1; 
  }   
 
   my $fastq = 0; 
   if( &arg_hasflag("fastq") ){
     $fastq = 1; 
  }    
   
  my $path = &arg_getremaining(0);
  
  if(-f $path ){
  
    if(private_isfastq($path)){
      &sub_fileinfo_fastq();
    }
    elsif(private_isfasta($path)){
      &sub_fileinfo_fasta();
    }
  }
  
  elsif(-d $path ){
  
    if($fastq==1){
      &sub_fileinfo_fastq();
    }
    elsif($fasta==1){
      &sub_fileinfo_fasta();
    }
    else{
      print STDERR "ERROR: If input is directory, filetype to scann must be given (--fasta --fastq) [Input dir=$path]\n";
      exit 0;
    }
  }
  else{
    print STDERR "ERROR: input path is not a file or a directory [path=$path]\n";
    exit 0;
  }
  
}

######################################### END: sub_fileinfo ##########################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================





#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\








#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
######################################### sub_fileinfo_fastq #########################################

# version = 2022.09.16

sub sub_fileinfo_fastq(){

  ##############
  #versioninfo
  ##############
  
  my $sub_name = "fileinfo_fastq";
  my $sub_version = "2022.09.16";
 
  my $subfunction_arg = shift(@_);
  if($subfunction_arg eq "version") { print "    component $sub_name version=$sub_version\n";  return; } 


  #############
  #arguments
  #############
  
  my $path = &arg_getremaining(0);  # this is the filename or the directoryname
  

  
  my $tabular = 0; # default = non tabular
  
  if(-d $path){
     $tabular = 1; # if input is directory, tabular is default
  }
  
  if( &arg_hasflag("tabular") ){
     $tabular = 1; 
  }
  if( &arg_hasflag("notabular") ){
     $tabular = 0; 
  }     
  
  my $printheader = 1; 
  if( &arg_hasflag("noheader") ){
     $printheader = 0;   # do not print headerrow in tabular mode
  }
  my $verbose = 0;
  if( &arg_hasflag("verbose") ){
     $verbose = 1;  # print more information
  }
 
  my $analyzepairs = 1;
  if( &arg_hasflag("unpaired") ){
     $analyzepairs = 0;  # do not co-analyze R1 and R2 pairs
  }
    
  my $recursive = 0; # recursively scan subdirs as well
  if( &arg_hasflag("recursive") ){
     $recursive = 1; 
  }
  
  &arg_warn_if_unused_argument();

  ###################################
  # prepare list of files to analyze
  ###################################
  
 
  #################
  # set input mode
  #################
  
  my $inputmode = 0;  # 1=single file 2=directory 3=recursive scan of directory
  
  if(-f $path){
     $inputmode = 1;
  }
 
  if(-d $path){
  
     if($recursive==1) {  
       $inputmode = 3;
     }
     else{
       $inputmode = 2;
     }
  }
  
  if($inputmode==0){
    print STDERR "ERROR: input path is not a file or a directory [path=$path]\n";
    exit 0;
  }
  

   
  ########################
  # prepare list of files
  ########################
  
  my @listoffiles = ();
    
  if($inputmode == 1){  #single file input
  
     if($analyzepairs==1){  # try to find the paired R2 (or R1 file)
     
       if(&private_hasR1Extension($path)){  ## is R1 - look for R2

          $R1_path = $path;
          $R2_path = private_toggleR1R2Extension($R1_path);
  
          if( -f $R1_path && -f $R2_path ){
             $path = "$R1_path\t$R2_path";
          }
       }
       elsif(&private_hasR2Extension($path)){ ## is R2 - look for R1
       
          $R2_path = $path;
          $R1_path = private_toggleR1R2Extension($R2_path);
          
          if( -f $R1_path && -f $R2_path ){
             $path = "$R1_path\t$R2_path";
          }
       }
     }
  
     push(@listoffiles, $path);
     
  } # end single file input
  
  if($inputmode == 2){  # directory input (but not recursive)
  
      &private_pushfilesindir($path,"fastq",\@listoffiles);
      if($analyzepairs==1){
        &private_mergeR1R2(\@listoffiles);
      }
  
  } 
  if($inputmode == 3){  # directory input recursively
      &private_pushfilesindir($path,"fastq",\@listoffiles,"recursive");
      if($analyzepairs==1){
        &private_mergeR1R2(\@listoffiles);
      }
  }  
 
 
  if(@listoffiles == 0){
    print STDERR "ERROR: No files to analyze [path=$path]\n";
    exit 0;
  }
  
  ###########################################
  # run anlysis
  ###########################################   
  
  my $filecounter = 0;

  ### Foreach file/filepair to analyze
  foreach my $analysispath (@listoffiles) {

  my $analysispath_R1 = "";
  my $analysispath_R2 = "";  
  my $pair = 0;

  if(index($analysispath,"\t")>=0){  # if the anlysispath variable contains a R1/R2 pair of files
  
     my @columns = split(/\t/,$analysispath);
     
     $analysispath_R1 = $columns[0];
     $analysispath_R2 = $columns[1];
     $pair =1;
     

  }

   if($pair==0){
     &private_filereader_open($analysispath,"FQ1");
   }
   if($pair==1){
     &private_filereader_open($analysispath_R1,"FQ1");
     &private_filereader_open($analysispath_R2,"FQ2");
   }

   $filecounter++;
   
   ## if tabular - make sure header only is printed first time
   if( $tabular==1  && $filecounter>1){
     $printheader = 0;
   } 
   
   ######################
   ### per-file analysis variables
      
   my $count_lines = 0;
   
   #### variables for singleton fastq or the R1 fastq if pair
   
   my @fastqseq = ();   
   
   my $count_seq = 0;
   my $count_bases = 0;
   my $count_initiallyanalysis = 0;
     
   my $count_A = 0;
   my $count_C = 0;
   my $count_G = 0;
   my $count_T = 0;
   my $count_N = 0;
   my $count_X = 0;
   my $count_Q30 = 0;  
   my $count_Q20 = 0;  
                      
   my $maxreadlength = -1;
   my $minreadlength = -1;

   my @specieshits = ();
   my @specieshits_discr = ();
   my @adaptorhits_12bp = ();
   my @adaptorhits_15bp = ();
   my @adaptorhits_full = ();
   

   #### variables for the R2 fastq if it pair
   
   my @fastqseq_R2 = ();  
   
   my $count_seq_R2 = 0;  
   my $count_bases_R2 = 0;        
   my $count_initiallyanalysis_R2 = 0;
   
   my $count_A_R2 = 0;
   my $count_C_R2 = 0;
   my $count_G_R2 = 0;
   my $count_T_R2 = 0;
   my $count_N_R2 = 0;
   my $count_X_R2 = 0;   
   my $count_Q30_R2 = 0;  
   my $count_Q20_R2 = 0; 
   
   my $maxreadlength_R2 = -1;
   my $minreadlength_R2 = -1;
   
   my @specieshits_R2 = ();
   my @specieshits_discr_R2 = ();
   my @adaptorhits_12bp_R2 = ();
   my @adaptorhits_15bp_R2 = ();
   my @adaptorhits_full_R2 = ();
   
           
   ###Fill array variables with zeroes 
   
   for(my $i=0;$i<@constants_species_16Sseq;$i++){
    push(@specieshits ,0);  
    push(@specieshits_nodiscr, 0);  
    push(@specieshits_R2 ,0);  
    push(@specieshits_nodiscr_R2, 0);    
   }      
   for(my $i=0;$i<@constants_adapters;$i++){
    push(@adaptorhits_12bp ,0);  
    push(@adaptorhits_15bp, 0);  
    push(@adaptorhits_full, 0);  
    push(@adaptorhits_12bp_R2 ,0);  
    push(@adaptorhits_15bp_R2, 0);  
    push(@adaptorhits_full_R2, 0);    
   }
      
    
   ###################################
   ## Here starts the loop through the current analysisfile/filepair
         
   do{

     ## read a sequence from the singleton/R1 fastq file
     @fastqseq = &private_filereader_getlines(4,"FQ1");
   
     ## Analyze singleton/R1 file sequence --but make sure it is not empty (end of file)
     if(@fastqseq==4){      
     
       $count_seq++;

       
       my $seq = $fastqseq[1];  
       my $qstring = $fastqseq[3];
       chop($seq); # remove newline character
       chop($qstring); # remove newline character
       
       my $seq_len = length($seq);

       $count_bases = $count_bases + $seq_len; # count total singleton/R1 bases
       
       
       ## keep track of longest and shortest read in singleton/R1 file
       if($maxreadlength == -1 || $seq_len>$maxreadlength){
	  $maxreadlength=$seq_len;
       }
       if($minreadlength == -1 || $seq_len<$minreadlength){
	  $minreadlength=$seq_len;
       }
	
       ## count base types	
       my $nrA  = $seq =~ tr/A/A/; 
       $count_A = $count_A + $nrA;
       my $nrC  = $seq =~ tr/C/C/; 
       $count_C = $count_C+ $nrC;
       my $nrG  = $seq =~ tr/G/G/; 
       $count_G = $count_G + $nrG;
       my $nrT  = $seq =~ tr/T/T/; 
       $count_T = $count_T+ $nrT;
       my $nrN  = $seq =~ tr/N/N/; 
       $count_N = $count_N+ $nrN;
       
       ## uncounted letters
       $count_X = $count_X + $seq_len - $nrA - $nrC - $nrG - $nrT - $nrN;
       
       ## count Q30+
       my $nrQ30plus  = $qstring =~ tr/?@ABCDEFGHI/?@ABCDEFGHI/; 
       $count_Q30 = $count_Q30 + $nrQ30plus;
       
       ## count Q20-
       my $nrQ20minus  = $qstring =~ tr|!"#$%&'()*+,-./012345|!"#$%&'()*+,-./012345|; 
       $count_Q20 = $count_Q20 + $nrQ20minus;
       
       ##############################################################################################
       ############################ Analysis made only on first inital part (default=100K) sequences
       if($count_seq<=$setting_fastqinitialanalysis){
       

          $count_initiallyanalysis++;  # this is the actual size of the initial analysis -- nesessary to keep track of e.g. if less than 100K reads are present

          ########## 16 S analysis
       
          my $seqtag = ""; 
          
          if(length($seq)>=$setting_16staglength){
             $seqtag = substr($seq,0,$setting_16staglength); # use first part (default 40 bases) as tag to scan 16S references
          }
          
          my $hit_index = -1;  # signals if there was a hit 
         
          ## search for 16S hits
          for(my $i=0;$i<@constants_species_16Sseq;$i++){  # for all 16S ref-sequences
          
            if($seqtag ne "" && index($constants_species_16Sseq[$i],$seqtag)>=0){  # if tag hits (and tag is suffecient long)
            
              $specieshits_nodiscr[$i] =  $specieshits_nodiscr[$i] +1;  # count hits - all hits
              if($hit_index==-1){  #first hit (index is still -1)
                 $hit_index = $i;
              }
              else{   # not first hit ... set hit index to -2 to signal non-discriminating hit
                 $hit_index = -2;
              }
            }
            elsif($seqtag ne "" && index($constants_species_16Sseq_rc[$i],$seqtag)>=0){  # repeat on the reverse complement of the 16S ref sequences

              $specieshits_nodiscr[$i] =  $specieshits_nodiscr[$i] +1;  # count hits - all hits
              if($hit_index==-1){  #first hit (index is still -1)
                 $hit_index = $i;
              }
              else{  # not first hit ... set hit index to -2 to signal non-discriminating hit
                 $hit_index = -2;
              }
            }
         }
         if( $hit_index>=0){  #there was a uniqe hit - count it
             $specieshits[$hit_index] =  $specieshits[$hit_index] +1;
         }
         
         
         ################ adapter analysis
         
         for(my $i=0;$i<@constants_adapters;$i++){
            if(index($seq, substr($constants_adapters_seq[$i],0,12))>=0){
              $adaptorhits_12bp[$i] =  $adaptorhits_12bp[$i] +1; 
            }
            
            if(index($seq, substr($constants_adapters_seq[$i],0,15))>=0){
              $adaptorhits_15bp[$i] =  $adaptorhits_15bp[$i] +1; 
            }
            
            if(index($seq, $constants_adapters_seq[$i])>=0){
              $adaptorhits_full[$i] =  $adaptorhits_full[$i] +1; 
            }
            
         }
         
         ##########
       
       } # end initial part (default first 100k) seq analysis
       ###################################################################################
 
       
     } #end if fastq read has 4 rows - should be 4 => a full read or 0 => end of file
     else{
       if(@fastqseq==3) { print STDERR "WARNING: incompleat sequence read (only 3 lines) in $analysispath\n";}  ## warn truncated seq
       if(@fastqseq==2) { print STDERR "WARNING: incompleat sequence read (only 2 lines) in $analysispath\n";}  ## warn truncated seq
       if(@fastqseq==1) { print STDERR "WARNING: incompleat sequence read (only 1 line) in $analysispath\n";}  ## warn truncated seq
     }
     ##############################################################
     ##############################################################
     ## this marks the end of the singleton/R1 sequence analysis
     ##############################################################
     
     
     #########################################################################
     ############################ R2 read analysis ###########################
     if($pair==1){

       ## read a sequence from the R2 fastq file     
       @fastqseq_R2 = &private_filereader_getlines(4,"FQ2");
     
       ## Analyze R2 file sequence --but make sure it is not empty (end of file)
       if(@fastqseq_R2==4){      

       my $seq = $fastqseq_R2[1];  # remove newline character
       my $qstring = $fastqseq_R2[3];  # remove newline character
       chop($seq); # remove newline character
       chop($qstring); # remove newline character
       
       my $seq_len = length($seq);
       
       $count_seq_R2++; 
     
       $count_bases_R2 = $count_bases_R2 + $seq_len; # count total R2 bases
       
       ## keep track of longest and shortest read in R2 file
       if($maxreadlength_R2 == -1 || $seq_len>$maxreadlength_R2){
	  $maxreadlength_R2=$seq_len;
       }
       if($minreadlength_R2 == -1 || $seq_len<$minreadlength_R2){
	  $minreadlength_R2=$seq_len;
       }

       ## count base types	
       my $nrA  = $seq =~ tr/A/A/; 
       $count_A_R2 = $count_A_R2 + $nrA;
       my $nrC  = $seq =~ tr/C/C/; 
       $count_C_R2 = $count_C_R2+ $nrC;
       my $nrG  = $seq =~ tr/G/G/; 
       $count_G_R2 = $count_G_R2 + $nrG;
       my $nrT  = $seq =~ tr/T/T/; 
       $count_T_R2 = $count_T_R2+ $nrT;
       my $nrN  = $seq =~ tr/N/N/; 
       $count_N_R2 = $count_N_R2+ $nrN;

       ## uncounted letters
       $count_X_R2 = $count_X_R2 + $seq_len - $nrA - $nrC - $nrG - $nrT - $nrN;

       ## count Q30+          
       my $nrQ30plus  = $qstring =~ tr/?@ABCDEFGHI/?@ABCDEFGHI/; 
       $count_Q30_R2 = $count_Q30_R2 + $nrQ30plus;

       ## count Q20-       
       my $nrQ20minus  = $qstring =~ tr|!"#$%&'()*+,-./012345|!"#$%&'()*+,-./012345|; 
       $count_Q20_R2 = $count_Q20_R2 + $nrQ20minus;
       
       
       ##############################################################################################
       ############################ Analysis made only on first inital part (default=100K) sequences  For R2 reads
       if($count_seq_R2<=$setting_fastqinitialanalysis){
       
       $count_initiallyanalysis_R2++; # this is the actual size of the initial analysis -- nesessary to keep track of e.g. if less than 100K reads are present
       
       ########## 16 S analysis

         my $seqtag = ""; 
          
         if(length($seq)>=$setting_16staglength){
             $seqtag = substr($seq,0,$setting_16staglength); # use first part (default 40 bases) as tag to scan 16S references
         }
                 
         my $hit_index = -1; 
         
         for(my $i=0;$i<@constants_species_16Sseq;$i++){  # for all 16S ref-sequences
            if($seqtag ne "" && index($constants_species_16Sseq[$i],$seqtag)>=0){  # if tag hits
              $specieshits_nodiscr_R2[$i] =  $specieshits_nodiscr_R2[$i] +1;  # count hits - all hits
              if($hit_index==-1){  #first hit (index is still -1)
                 $hit_index = $i;
              }
              else{   # not first hit ... set hit index to -2 to signal non-discriminating hit
                 $hit_index = -2;
              }
            }
            elsif($seqtag ne "" && index($constants_species_16Sseq_rc[$i],$seqtag)>=0){  # repeat on the reverse complement of the 16S ref sequences
              $specieshits_nodiscr_R2[$i] =  $specieshits_nodiscr_R2[$i] +1;  # count hits - all hits
              if($hit_index==-1){  #first hit (index is still -1)
                 $hit_index = $i;
              }
              else{  # not first hit ... set hit index to -2 to signal non-discriminating hit
                 $hit_index = -2;
              }
            }
         }
         if( $hit_index>=0){  #there was a uniqe hit - count it
             $specieshits_R2[$hit_index] =  $specieshits_R2[$hit_index] +1;
         }
         
         
         ################ adapter analysis
         
         for(my $i=0;$i<@constants_adapters;$i++){
            if(index($seq, substr($constants_adapters_seq[$i],0,12))>=0){
              $adaptorhits_12bp_R2[$i] =  $adaptorhits_12bp_R2[$i] +1; 
            }
            
            if(index($seq, substr($constants_adapters_seq[$i],0,15))>=0){
              $adaptorhits_15bp_R2[$i] =  $adaptorhits_15bp_R2[$i] +1; 
            }
            
            if(index($seq, $constants_adapters_seq[$i])>=0){
              $adaptorhits_full_R2[$i] =  $adaptorhits_full_R2[$i] +1; 
            }
            
         }
         
         ##########
       
       } # end initial part (default first 100k) seq analysis
       ###################################################################################
       
     
       
     } #end if R2 fastq read has 4 rows - should be 4 => a full read or 0 => end of file
     else{
       if(@fastqseq_R2==3) { print STDERR "WARNING: incompleat R2 sequence read (only 3 lines) in $analysispath\n";}
       if(@fastqseq_R2==2) { print STDERR "WARNING: incompleat R2 sequence read (only 2 lines) in $analysispath\n";}
       if(@fastqseq_R2==1) { print STDERR "WARNING: incompleat R2 sequence read (only 1 line) in $analysispath\n";}
     }
     
     } # end if paired = if R2 is to be analyzed

     ##############################################################
     ##############################################################
     ## this marks the end of the R2 sequence analysis
     ##############################################################     

     


  }while(@fastqseq==4); #while sequences left to analyze in the file

  &private_filereader_close("FQ1");

  if($pair==1){
   &private_filereader_close("FQ2");
  }
  
  ##############################################################
  ##############################################################
  ## this marks the end of the sequence analysis loop (both singleton/R1 and R2)
  ############################################################## 


  ##################################################################
  ##################################################################
  ## Calculate some metrix based on the collected data
  ################################################################

  ## singleton/R1 metrix
  
  my $thepath = $analysispath;
  my $thefile = &private_getlastpartofpath($thepath);
         
  my $avlen=0;

  my $percent_Q30 = 0;
  my $percent_Q20 = 0;

  my $percent_GC = 0;
  my $percent_A = 0;
  my $percent_C = 0;
  my $percent_G = 0;
  my $percent_T = 0;
  my $percent_N = 0;
  
  if($count_seq>0){
    $avlen = $count_bases/$count_seq;
  }
  
  if($count_bases>0){
     $percent_Q30 = $count_Q30/$count_bases*100;
     $percent_Q20 = $count_Q20/$count_bases*100;
  }
  
  if(($count_G + $count_C + $count_T + $count_A) > 0){
    $percent_GC = ($count_G + $count_C)/ ($count_G + $count_C + $count_T + $count_A) *100.0;
    $percent_A = ($count_A)/ ($count_G + $count_C + $count_T + $count_A ) *100.0;
    $percent_C = ($count_C)/ ($count_G + $count_C + $count_T + $count_A ) *100.0;
    $percent_G = ($count_G)/ ($count_G + $count_C + $count_T + $count_A ) *100.0;
    $percent_T = ($count_T)/ ($count_G + $count_C + $count_T + $count_A ) *100.0;
  }
  if(($count_G + $count_C + $count_T + $count_A + $count_N) > 0){
    $percent_N = ($count_N)/ ($count_G + $count_C + $count_T + $count_A + $count_N) *100.0;  
  }
  
  
  
  ################################
  ##### 16S hits

  my @specieshits_formated = ();
  my @specieshits_organism = ();
  my @specieshits_value = ();
  my @specieshits_percent = ();
  my @specieshits_thresholds = ();  
  my $specieshits_count = 0;
  
  my %tobeprinted_hash = ();
  my %tobeprinted_hash_nodiscr = ();
    
  ## pass 1 -- format data and assocate with their counts
    for(my $i=0;$i<@constants_species_16Sseq;$i++){ 
    
      #discriminating (unique) hits
      if($specieshits[$i]>0) {
      
        my $tresholdtouse = $constants_species_16Sthreshold[$i];
        
        my $percenthit = $specieshits[$i]/$count_initiallyanalysis*100.0;
        my $hitoutput = "";
          if($percenthit>$tresholdtouse ) { 
            $hitoutput = $specieshits[$i] . " (" . sprintf("%.3f", $percenthit ) . " %)  [OK]"; }
          else { 
            $hitoutput = $specieshits[$i] . " (" . sprintf("%.3f", $percenthit ) . " %)  [WARNING - LOW LEVELS]"; 
          }
 
          my $formated_output =  $constants_species[$i] . " " . $hitoutput;
          $tobeprinted_hash{ $formated_output . "\t" . $constants_species[$i]  . "\t" . $constants_species_16Sthreshold[$i]  } =  $specieshits[$i]  ;
          
      }
      
      #non-discriminating (hitting multiple species)
      if($specieshits_nodiscr[$i]>0) {
          my $percenthit = $specieshits_nodiscr[$i]/$count_initiallyanalysis*100.0;
          $tobeprinted_hash_nodiscr{ $constants_species[$i] . " " . $specieshits_nodiscr[$i] . " (" .  sprintf("%.3f", $percenthit ) . " %)\n" } =  $specieshits_nodiscr[$i]  ;
      }
      
    }
    ## pass 2 -- sort and make printout arrays

    my $all16shits_formated = "";
    my $all16shits_formated_lines = "";
    my $all16shits_nodiscr_formated_lines = "";
    
    foreach my $name (sort { $tobeprinted_hash{$b} <=> $tobeprinted_hash{$a} } keys %tobeprinted_hash) {
  
         my @columns = split(/\t/,$name); 
  
         my $formated = $columns[0];
         my $organism = $columns[1];
         my $threshold_16s = $columns[1];
         my $hitcount = $tobeprinted_hash{$name};
         my $percenthit =   $hitcount/$count_initiallyanalysis*100.0;
      
          push(@specieshits_organism,$organism);
          push(@specieshits_formated,$formated);
          push(@specieshits_value,$hitcount);
          push(@specieshits_percent,$percenthit);
          push(@specieshits_thresholds,$threshold_16s);
          $specieshits_count++;
          
          $all16shits_formated_lines = $all16shits_formated_lines . $formated . "\n";
          
          if($specieshits_count==1){
            $all16shits_formated = $formated;
          }
          else{
            $all16shits_formated = $all16shits_formated . ":" . $formated;
          }
    }



    if($all16shits_formated eq ""){
       $all16shits_formated = "none";
       $all16shits_formated_lines = "none\n";
    }
    
    my $best16Shit =  "none";
    my $best16Shit_formated =  "none";
    my $best16Shit_percent  =  "0";
    my $best16Shit_threshold =  0;
    my $best16Shit_cov =  "NA";
    my $best16Shit_covQ30 = "NA";


    if(@specieshits_organism>=1){
       $best16Shit =  $specieshits_organism[0];
       $best16Shit_formated =  $specieshits_formated[0];
       $best16Shit_percent  =  $specieshits_percent[0];
       $best16Shit_threshold = $specieshits_thresholds[0];
    }


    
    ### non-disriminating hits formated output
    foreach my $name (sort { $tobeprinted_hash_nodiscr{$b} <=> $tobeprinted_hash_nodiscr{$a} } keys %tobeprinted_hash_nodiscr) {
      $all16shits_nodiscr_formated_lines = $all16shits_nodiscr_formated_lines . $name;
    }
    
    if($all16shits_nodiscr_formated_lines eq ""){
       $all16shits_nodiscr_formated_lines = "none\n";
    }
    
    #### END 16S hits
    ############################################
    
    ############################################
    ## ADAPTER hits
 
    my @adaptorhits_formated = ();
    my @adaptorhits_value_12 = ();
    my @adaptorhits_value_15 = ();
    my @adaptorhits_value_full = ();        
    my @adaptorhits_percent_12 = ();
    my @adaptorhits_percent_15 = ();
    my @adaptorhits_percent_full = ();
    my $adaptorhits_count = 0;
  
    my %tobeprinted_hash_adaptors = ();
    
    ## pass 1 -- format printouts and assocate with their counts
    for(my $i=0;$i<@constants_adapters_seq;$i++){

        if($adaptorhits_12bp[$i]>0) {
        
          $percenthit_12 = $adaptorhits_12bp[$i]/$count_initiallyanalysis*100.0;
          $percenthit_15 = $adaptorhits_15bp[$i]/$count_initiallyanalysis*100.0;
          $percenthit_full = $adaptorhits_full[$i]/$count_initiallyanalysis*100.0;
          
          $formated = $constants_adapters[$i] . " 12bp=" . $adaptorhits_12bp[$i] . " ($percenthit_12%)" . " 15bp=" . $adaptorhits_15bp[$i] . " ($percenthit_15%)" . " full_seq=" . $adaptorhits_full[$i]  . " ($percenthit_full%)" ;
          
          $tobeprinted_hash_adaptors{ $formated . "\t" . $adaptorhits_12bp[$i] . "\t" . $adaptorhits_15bp[$i] . "\t" . $adaptorhits_full[$i] . "\t" . $percenthit_12 . "\t" . $percenthit_15 . "\t" . $percenthit_full } =  $adaptorhits_12bp[$i]  ;

        }
    }
  
    ## pass 2 -- sort and printout in order
    
    my $allAdapterhits_formated = "";
    my $allAdapterhits_formated_lines = "";
    
    foreach my $name (sort { $tobeprinted_hash_adaptors{$b} <=> $tobeprinted_hash_adaptors{$a} } keys %tobeprinted_hash_adaptors) {

        my @columns = split(/\t/,$name); 

        my $formated = $columns[0];
        my $hits_12 = $columns[1];
        my $hits_15 = $columns[2];
        my $hits_full = $columns[3];
        my $hits_12_percent = $columns[4];
        my $hits_15_percent = $columns[5];
        my $hits_full_percent= $columns[6];

        push(@adaptorhits_formated,$formated);
        push(@adaptorhits_value_12,$hits_12);
        push(@adaptorhits_value_15,$hits_15);
        push(@adaptorhits_value_full,$hits_full);
        push(@adaptorhits_percent_12,$hits_12_percent);
        push(@adaptorhits_percent_15,$hits_15_percent);
        push(@adaptorhits_percent_full,$hits_full_percent);

	$adaptorhits_count++;

	$allAdapterhits_formated_lines = $allAdapterhits_formated_lines . $formated . "\n";
		                
       if($adaptorhits_count==1) {
          $allAdapterhits_formated =  $formated;
       } 
       else{ 
          $allAdapterhits_formated = $allAdapterhits_formated . ":" . $formated;
       } 
   }
  
   if($allAdapterhits_formated eq ""){
       $allAdapterhits_formated = "none";
       $allAdapterhits_formated_lines = "none\n";
    }
    
    my $bestAdapterhit_formated =  "none";
    my $bestAdapterhit_percent  =  "0";


    if(@adaptorhits_formated>=1){
       $bestAdapterhit_formated =  $adaptorhits_formated[0];
       $bestAdapterhit_percent  =  $adaptorhits_value_12[0]/$count_initiallyanalysis*100.0;
    }


  ############################################
    ## ADAPTER hits R2
 
    my @adaptorhits_formated_R2 = ();
    my @adaptorhits_value_12_R2 = ();
    my @adaptorhits_value_15_R2 = ();
    my @adaptorhits_value_full_R2 = ();        
    my @adaptorhits_percent_12_R2 = ();
    my @adaptorhits_percent_15_R2 = ();
    my @adaptorhits_percent_full_R2 = ();
    my $adaptorhits_count_R2 = 0;
  
    my %tobeprinted_hash_adaptors_R2 = ();
    
    ## pass 1 -- format printouts and assocate with their counts
    for(my $i=0;$i<@constants_adapters_seq;$i++){

        if($adaptorhits_12bp_R2[$i]>0) {
        
          $percenthit_12_R2 = $adaptorhits_12bp_R2[$i]/$count_initiallyanalysis_R2*100.0;
          $percenthit_15_R2 = $adaptorhits_15bp_R2[$i]/$count_initiallyanalysis_R2*100.0;
          $percenthit_full_R2 = $adaptorhits_full_R2[$i]/$count_initiallyanalysis_R2*100.0;
          
          $formated_R2 = $constants_adapters[$i] . " 12bp=" . $adaptorhits_12bp_R2[$i] . " ($percenthit_12_R2%)" . " 15bp=" . $adaptorhits_15bp_R2[$i] . " ($percenthit_15_R2%)" . " full_seq=" . $adaptorhits_full_R2[$i]  . " ($percenthit_full_R2%)" ;
          
          $tobeprinted_hash_adaptors_R2{ $formated_R2 . "\t" . $adaptorhits_12bp_R2[$i] . "\t" . $adaptorhits_15bp_R2[$i] . "\t" . $adaptorhits_full_R2[$i] . "\t" . $percenthit_12_R2 . "\t" . $percenthit_15_R2 . "\t" . $percenthit_full_R2 } =  $adaptorhits_12bp_R2[$i]  ;

        }
    }
  
    ## pass 2 -- sort and printout in order

    my $allAdapterhits_formated_R2 = "";
    my $allAdapterhits_formated_lines_R2 = "";
    
    foreach my $name (sort { $tobeprinted_hash_adaptors_R2{$b} <=> $tobeprinted_hash_adaptors_R2{$a} } keys %tobeprinted_hash_adaptors_R2) {

        my @columns = split(/\t/,$name); 
  
        my $formated = $columns[0];
        my $hits_12 = $columns[1];
        my $hits_15 = $columns[2];
        my $hits_full = $columns[3];
        my $hits_12_percent = $columns[4];
        my $hits_15_percent = $columns[5];
        my $hits_full_percent= $columns[6];

        push(@adaptorhits_formated_R2,$formated);
        push(@adaptorhits_value_12_R2,$hits_12);
        push(@adaptorhits_value_15_R2,$hits_15);
        push(@adaptorhits_value_full_R2,$hits_full);
        push(@adaptorhits_percent_12_R2,$hits_12_percent);
        push(@adaptorhits_percent_15_R2,$hits_15_percent);
        push(@adaptorhits_percent_full_R2,$hits_full_percent);

	$adaptorhits_count_R2++;

	$allAdapterhits_formated_lines_R2 = $allAdapterhits_formated_lines_R2 . $formated . "\n";
		                
       if($adaptorhits_count_R2==1) {
          $allAdapterhits_formated_R2 =  $formated;
       } 
       else{ 
          $allAdapterhits_formated_R2 = $allAdapterhits_formated_R2 . ":" . $formated;
       } 
   }
  
   if($allAdapterhits_formated_R2 eq ""){
       $allAdapterhits_formated_R2 = "none";
       $allAdapterhits_formated_lines_R2 = "none\n";
    }
    
    my $bestAdapterhit_formated_R2 =  "none";
    my $bestAdapterhit_percent_R2  =  "0";
    my $bestAdapterhit_percent_R1R2  =  "0";

    if(@adaptorhits_formated_R2>=1){
       $bestAdapterhit_formated_R2 =  $adaptorhits_formated_R2[0];
       $bestAdapterhit_percent_R2  =  $adaptorhits_value_12_R2[0]/$count_initiallyanalysis_R2*100.0;
       
       $bestAdapterhit_percent_R1R2 = $bestAdapterhit_percent;
       if( $bestAdapterhit_percent_R2 > $bestAdapterhit_percent_R1R2){
          $bestAdapterhit_percent_R1R2 = $bestAdapterhit_percent_R2;
       }
    }
    
    ##################### END ADAPTER HITS    

  ####
  
  ## R2 and pair metrix

    my $thepath_R1 = "";
    my $thepath_R2 = "";
    my $thefile_R1  = "";
    my $thefile_R2 = "";
    
    my $avlen_R1=0;
    my $avlen_R2=0;
    my $avlen_R1R2 = 0;
          
    my $count_seq_R1R2 = 0;
    my $count_bases_R1R2 = 0;
    my $count_Q30_R1R2  = 0;
    my $count_N_R1R2 = 0;
    my $percent_Q30_R1 = 0;
    my $percent_Q30_R2 = 0;
    my $percent_Q30_R1R2 = 0;
    my $count_Q20_R1R2 = 0;
    my $percent_Q20_R1 = 0;
    my $percent_Q20_R2 = 0;
    my $percent_Q20_R1R2 = 0;

    my $percent_GC_R2 = 0;
    my $percent_A_R2 = 0;
    my $percent_C_R2 = 0;
    my $percent_G_R2 = 0;
    my $percent_T_R2 = 0;
    my $percent_N_R2 = 0;

    my $percent_GC_R1R2 = 0;
    my $percent_A_R1R2 = 0;
    my $percent_C_R1R2 = 0;
    my $percent_G_R1R2 = 0;
    my $percent_T_R1R2 = 0;
    my $percent_N_R1R2 = 0;
    
    my $minreadlength_R1R2 = 0;
    my $maxreadlength_R1R2 = 0;
          
    if($pair==1) { 
    
        $thepath_R1 = $analysispath_R1; 
        $thepath_R2 = $analysispath_R2;
     
        $thefile_R1 = &private_getlastpartofpath($thepath_R1);
        $thefile_R2 = &private_getlastpartofpath($thepath_R2);
     
        $count_seq_R1R2 = $count_seq + $count_seq_R2;
        
        $count_bases_R1R2 = $count_bases + $count_bases_R2;
        $count_Q30_R1R2 = $count_Q30 + $count_Q30_R2;
        
        if($count_seq_R2>0){
           $avlen_R2 = $count_bases_R2/$count_seq_R2;
        }
        if($count_seq>0 && $count_seq_R2>0){
           $avlen_R1R2 = ($count_bases+$count_bases_R2)/($count_seq+$count_seq_R2);
        }  
        
        if($count_bases>0){
          $percent_Q30_R1 = $count_Q30/$count_bases*100;
          $percent_Q20_R1 = $count_Q20/$count_bases*100;
          
          if($count_bases_R2>0){
            $percent_Q30_R2 = $count_Q30_R2/$count_bases_R2*100;
          }  
          $percent_Q30_R1R2 = ($count_Q30 + $count_Q30_R2)/($count_bases + $count_bases_R2)*100;

          $count_Q20_R1R2 = $count_Q20 + $count_Q20_R2;
          $percent_Q20_R1 = $count_Q20/$count_bases*100;
          if($count_bases_R2>0){
            $percent_Q20_R2 = $count_Q20_R2/$count_bases_R2*100;
          }  
          $percent_Q20_R1R2 = ($count_Q20 + $count_Q20_R2)/($count_bases + $count_bases_R2)*100;
        }
  
        if(($count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2) > 0){
          $percent_GC_R2 = ($count_G_R2 + $count_C_R2)/ ($count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2) *100.0;
          $percent_A_R2 = ($count_A_R2)/ ($count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 ) *100.0;
          $percent_C_R2 = ($count_C_R2)/ ($count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 ) *100.0;
          $percent_G_R2 = ($count_G_R2)/ ($count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 ) *100.0;
          $percent_T_R2 = ($count_T_R2)/ ($count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 ) *100.0;
        }
        if(($count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 + $count_N_R2) > 0){
          $percent_N_R2 = ($count_N_R2)/ ($count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 + $count_N_R2) *100.0;           
        } 
        if(($count_G + $count_C + $count_T + $count_A + $count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2) > 0){
          $percent_GC_R1R2 = ($count_G + $count_C + $count_G_R2 + $count_C_R2)/ ($count_G + $count_C + $count_T + $count_A + $count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2) *100.0;
          $percent_A_R1R2 = ($count_A + $count_A_R2)/ ($count_G + $count_C + $count_T + $count_A + $count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 ) *100.0;
          $percent_C_R1R2 = ($count_C + $count_C_R2)/ ($count_G + $count_C + $count_T + $count_A + $count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 ) *100.0;
          $percent_G_R1R2 = ($count_G + $count_G_R2)/ ($count_G + $count_C + $count_T + $count_A + $count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 ) *100.0;
          $percent_T_R1R2 = ($count_T + $count_T_R2)/ ($count_G + $count_C + $count_T + $count_A + $count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 ) *100.0;
        }
        if(($count_G + $count_C + $count_T + $count_A + $count_N + $count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 + $count_N_R2) > 0){
          $percent_N_R1R2 = ($count_N + $count_N_R2)/ ($count_G + $count_C + $count_T + $count_A + $count_N + $count_G_R2 + $count_C_R2 + $count_T_R2 + $count_A_R2 + $count_N_R2) *100.0;          
        }
        
        $count_N_R1R2 = $count_N + $count_N_R2;
          
        $minreadlength_R1R2 = $minreadlength;
        $maxreadlength_R1R2 = $maxreadlength;
  
        if($minreadlength_R2<$minreadlength){
            $minreadlength_R1R2 = $minreadlength_R2;
        }
        if($maxreadlength_R2>$maxreadlength){
            $maxreadlength_R1R2 = $maxreadlength_R2;
        }
  

    } # end is paired analysis (calc metrix)


    if($best16Shit_percent>$best16Shit_threshold){

        $tmp=$best16Shit;
        $tmp =~ s/\s/_/g;  #remove spaces
        my $genomesize=  &private_gettypicalgenomesize($tmp);

       if($genomesize>0 && $pair==1){
         $best16Shit_cov = $count_bases_R1R2/$genomesize ;  
         $best16Shit_covQ30 =  $count_Q30_R1R2/$genomesize  ; 
       }
       else{
         $best16Shit_cov = $count_bases/$genomesize ; 
         $best16Shit_covQ30 =  $count_Q30/$genomesize  ;    
       }
    }
    

  #####################################################################
  #####################################################################
  ## OUTPUT the results
  #####################################################################
  #####################################################################



  #####################################
  # TABULAR UNPAIRED OUTPUT
 
  if($tabular==1 && $pair==0){  

    if(&arg_hasflag("verbose")){  
      if($printheader==1){ #  print header
        print "File\tPath\tNrReads\tNrBases\tQ30plus_bases\tQ30plus_percent\tQ20less_bases\tQ20plus_percent" 
        . "\tGC\tA\tC\tG\tT\tN\tMinReadlength\tMaxReadlength\tAverageReadlength"
        . "\tBest_16SHits\tBest_16SHits_percent\tCoverage\tCoverageQ30+\tAll_16SHits\tAdapters_percent_12bphit\tAdapterHits\n";
      }
 
      print "$thefile\t$thepath\t$count_seq\t$count_bases\t$count_Q30\t$percent_Q30\t$count_Q20\t$percent_Q20"
      . "\t$percent_GC\t$percent_A\t$percent_C\t$percent_G\t$percent_T\t$percent_N\t$minreadlength\t$maxreadlength\t$avlen\t"
      . "$best16Shit\t$best16Shit_percent\t$best16Shit_cov\t$best16Shit_covQ30\t$all16shits_formated\t"
      . "$bestAdapterhit_percent\t$allAdapterhits_formated\n";  
    }
    else{
      if($printheader==1){ #  print header
        print "File\tPath\tNrReads\tNrBases\tQ30plus_bases\tQ30plus_percent\tQ20less_bases\tQ20plus_percent" 
        . "\tN\tMinReadlength\tMaxReadlength\tAverageReadlength"
        . "\tBest_16SHits\tBest_16SHits_percent\tCoverage\tCoverageQ30+\tAll_16SHits\tAdapters_percent_12bphit\tAdapterHits\n";
      }
 
      print "$thefile\t$thepath\t$count_seq\t$count_bases\t$count_Q30\t$percent_Q30\t$count_Q20\t$percent_Q20"
      . "\t$percent_N\t$minreadlength\t$maxreadlength\t$avlen\t"
      . "$best16Shit\t$best16Shit_percent\t$best16Shit_cov\t$best16Shit_covQ30\t$all16shits_formated\t"
      . "$bestAdapterhit_percent\t$allAdapterhits_formated\n";  
    }
    

  
 } ## END TABULAR UNPAIRED OUTPUT 
###########################################


#################################################
#### TABULAR AND PAIRED OUTPUT (R1+R2)

  if($tabular==1 && $pair==1){  
  
    if(&arg_hasflag("verbose")){             
      if($printheader==1){ # print header
        print "R1_File\tR2_File\tR1_path\tR2_path\tNrReads_R1\tNrReads_R2\tNrReads_Total\tNrBases_R1\tNrBases_R2\tNrBases_Total\tQ30plus_bases_R1\tQ30plus_bases_R1_percent\tQ30plus_bases_R2\tQ30plus_bases_R2_percent\tQ30plus_bases_Total\tQ30plus_bases_Total_percent" 
        . "\tQ20less_bases_R1\tQ20less_bases_R1_percent\tQ20less_bases_R2\tQ20less_bases_R2_percent\tQ20less_bases_Total\tQ20less_bases_Total_percent"
        . "\tGC_R1\tGC_R2\tGC_Total\tA_R1\tC_R1\tG_R1\tT_R1\tN_R1\tA_R2\tC_R2\tG_R2\tT_R2\tN_R2\tA_Total\tC_Total\tG_Total\tT_Total\tN_Total"
        . "\tMinReadlength_R1\tMaxReadlength_R1\tAverageReadlength_R1\tMinReadlength_R2\tMaxReadlength_R2\tAverageReadlength_R2\tMinReadlength_Total\tMaxReadlength_Total\tAverageReadlength_Total"
        . "\tBest_16SHits\tBest_16SHits_percent\tCoverage\tCoverageQ30+\tAll_16SHits\tAdapters_percent_12bphit\tAdapterHitsR1\tAdapterHitsR2\n";
      }
 
      print "$thefile_R1\t$thefile_R2\t$thepath_R1\t$thepath_R2\t$count_seq\t$count_seq_R2\t$count_seq_R1R2\t$count_bases\t$count_bases_R2\t$count_bases_R1R2\t$count_Q30\t$percent_Q30_R1\t$count_Q30_R2\t$percent_Q30_R2\t$count_Q30_R1R2\t$percent_Q30_R1R2"
        . "\t$count_Q20\t$percent_Q20_R1\t$count_Q20_R2\t$percent_Q20_R2\t$count_Q20_R1R2\t$percent_Q20_R1R2"
        . "\t$percent_GC\t$percent_GC_R2\t$percent_GC_R1R2\t$percent_A\t$percent_C\t$percent_G\t$percent_T\t$percent_N\t$percent_A_R2\t$percent_C_R2\t$percent_G_R2\t$percent_T_R2\t$percent_N_R2\t$percent_A_R1R2\t$percent_C_R1R2\t$percent_G_R1R2\t$percent_T_R1R2\t$percent_N_R1R2"
        . "\t$minreadlength\t$maxreadlength\t$avlen\t$minreadlength_R2\t$maxreadlength_R2\t$avlen_R2\t$minreadlength_R1R2\t$maxreadlength_R1R2\t$avlen_R1R2\t"
        . "$best16Shit\t$best16Shit_percent\t$best16Shit_cov\t$best16Shit_covQ30\t$all16shits_formated\t"
        . "$bestAdapterhit_percent_R1R2\t$allAdapterhits_formated\t$allAdapterhits_formated_R2\n";           

  } #end verbose
  else{
   if($printheader==1){ #  print header
    print "R1_File\tR2_File\tR1_path\tR2_path\tNrReads_R1\tNrReads_R2\tNrReads_Total\tNrBases_R1\tNrBases_R2\tNrBases_Total\tQ30plus_bases_R1\tQ30plus_bases_R1_percent\tQ30plus_bases_R2\tQ30plus_bases_R2_percent\tQ30plus_bases_Total\tQ30plus_bases_Total_percent" 
    . "\tQ20less_bases_R1\tQ20less_bases_R1_percent\tQ20less_bases_R2\tQ20less_bases_R2_percent\tQ20less_bases_Total\tQ20less_bases_Total_percent"
    . "\tGC_Total\tN_Total"
    . "\tMinReadlength_Total\tMaxReadlength_Total\tAverageReadlength_Total"
    . "\tBest_16SHits\tBest_16SHits_percent\tCoverage\tCoverageQ30+\tAll_16SHits\tAdapters_percent_12bphit\tAdapterHitsR1\tAdapterHitsR2\n";
  }
 
  print "$thefile_R1\t$thefile_R2\t$thepath_R1\t$thepath_R2\t$count_seq\t$count_seq_R2\t$count_seq_R1R2\t$count_bases\t$count_bases_R2\t$count_bases_R1R2\t$count_Q30\t$percent_Q30_R1\t$count_Q30_R2\t$percent_Q30_R2\t$count_Q30_R1R2\t$percent_Q30_R1R2"
       . "\t$count_Q20\t$percent_Q20_R1\t$count_Q20_R2\t$percent_Q20_R2\t$count_Q20_R1R2\t$percent_Q20_R1R2"
       . "\t$percent_GC_R1R2\t$percent_N_R1R2"
       . "\t$minreadlength_R1R2\t$maxreadlength_R1R2\t$avlen_R1R2\t"
       . "$best16Shit\t$best16Shit_percent\t$best16Shit_cov\t$best16Shit_covQ30\t$all16shits_formated\t"
       . "$bestAdapterhit_percent_R1R2\t$allAdapterhits_formated\t$allAdapterhits_formated_R2\n";   
  }
  
  

  
 } ## END TABULAR AND PAIRED OUTPUT (R1+R2)
 #################################################
 ################################################  
 
 
 ##############################################################
 ######## NON-TABULAR UNPAIRED OUTPUT ############################# 

 
  if($tabular==0 && $pair==0){  
  
    #### format data for printout
  
  my $format_count_seq      = private_formatstring_readcount($count_seq);
  my $format_count_bases    = private_formatstring_basecount($count_bases);
  my $format_count_Q30 = private_formatstring_basecount($count_Q30);
  my $format_count_Q20 = private_formatstring_basecount($count_Q20);
 
  my $format_avlen = sprintf("%.2f", $avlen ) . " bases";

  my $format_percent_Q30 = sprintf("%.2f", $percent_Q30 ) . "%";
  my $format_percent_Q20 = sprintf("%.2f", $percent_Q20 ) . "%";
  
  my $format_percent_GC = sprintf("%.2f", $percent_GC ) . "%";
  my $format_percent_A = sprintf("%.2f", $percent_A ) . "%";
  my $format_percent_C = sprintf("%.2f", $percent_C ) . "%";
  my $format_percent_G = sprintf("%.2f", $percent_G ) . "%";
  my $format_percent_T = sprintf("%.2f", $percent_T ) . "%";
  my $format_percent_N = sprintf("%.4f", $percent_N ) . "%";
  my $format_adaptors = sprintf("%.3f", $bestAdapterhit_percent ) . "%";  
  
  print "File:          $analysispath\n";
  print "Nr reads:      $count_seq  ($format_count_seq)\n";
  print "Nr bases:      $count_bases ($format_count_bases)\n"; 
  print "bases >= Q30:  $count_Q30 ($format_count_Q30,  $format_percent_Q30)\n";
  print "bases <= Q20:  $count_Q20 ($format_count_Q20,  $format_percent_Q20)\n";
  print "GC content:    $format_percent_GC\n";
  if(&arg_hasflag("verbose")){
  print "bases:         A = $format_percent_A, C = $format_percent_C, G = $format_percent_G, T = $format_percent_T\n";
  }
  print "'N' bases:     $count_N ($format_percent_N) \n"; 
  print "Seqlength:     $minreadlength-$maxreadlength (average = $format_avlen)\n";

 if($best16Shit_percent>$best16Shit_threshold){
   print "Species:       $best16Shit\n";

   if($best16Shit_cov ne "NA"){
    $cov = sprintf("%.1f", $best16Shit_cov ) . "X"; 
    $covQ30 = sprintf("%.1f", $best16Shit_covQ30 ) . "X"; 
    print "Coverage:      $cov\n";
    print "Cov. (Q30+):   $covQ30\n";  
   } 
  }
  print "Adaptors:      $format_adaptors\n";
  
  ###############
  print "\nDiscriminating 16S sequence-tags detected in $count_initiallyanalysis reads: \n";
  print "$all16shits_formated_lines\n";
  
  #### VERBOSE -NON-TABULAR OUTPUT - also print non-unique 16S hits
  
   if( &arg_hasflag("verbose") ){
    
     print "\nNon-discriminating 16S sequence-tags detected in $count_initiallyanalysis reads: \n";
     print "$all16shits_nodiscr_formated_lines\n";
      
   } # end verbose
   ################
  
  print "\nAdaptor sequences detected in $count_initiallyanalysis  reads: \n";
  
  print "$allAdapterhits_formated_lines";
  
  
  print "-------------------\n";
  }  # end if pair==0
  #############################################
  
  
  
 ##############################################################
 ######## NON-TABULAR PAIRED OUTPUT ############################# 
 
  if($tabular==0 && $pair==1){  

  my $format_count_seq_R1      = private_formatstring_readcount($count_seq);
  my $format_count_seq_R2      = private_formatstring_readcount($count_seq_R2);
  my $format_count_seq_R1R2      = private_formatstring_readcount($count_seq_R1R2);
  my $format_count_seq_pairs   = private_formatstring_readpaircount($count_seq);
    
  my $format_count_bases_R1    = private_formatstring_basecount($count_bases);
  my $format_count_bases_R2    = private_formatstring_basecount($count_bases_R2);
  my $format_count_bases_R1R2    = private_formatstring_basecount($count_bases_R1R2);  
  
  my $format_count_Q30bases_R1 = private_formatstring_basecount($count_Q30);
  my $format_count_Q30bases_R2 = private_formatstring_basecount($count_Q30_R2);
  my $format_count_Q20bases_R1 = private_formatstring_basecount($count_Q20);  
  my $format_count_Q20bases_R2 = private_formatstring_basecount($count_Q20_R2);
   
  my $format_count_Q30bases_R1R2 = private_formatstring_basecount($count_Q30_R1R2);
  my $format_count_Q20bases_R1R2 = private_formatstring_basecount($count_Q20_R1R2);
 
  my $format_avlen = sprintf("%.2f", $avlen ) . " bases";  
  my $format_avlen_R2 = sprintf("%.2f", $avlen_R2 ) . " bases";
  my $format_avlen_R1R2 = sprintf("%.2f", $avlen_R1R2 ) . " bases";

  my $format_percent_Q30_R1 = sprintf("%.2f", $percent_Q30_R1 ) . "%";
  my $format_percent_Q20_R1 = sprintf("%.2f", $percent_Q20 ) . "%";
  my $format_percent_Q30_R2 = sprintf("%.2f", $percent_Q30_R2 ) . "%";
  my $format_percent_Q20_R2 = sprintf("%.2f", $percent_Q20_R2 ) . "%";
  my $format_percent_Q30_R1R2 = sprintf("%.2f", $percent_Q30_R1R2 ) . "%";
  my $format_percent_Q20_R1R2 = sprintf("%.2f", $percent_Q20_R1R2 ) . "%";
  
  my $format_percent_GC_R1 = sprintf("%.2f", $percent_GC ) . "%";
  my $format_percent_A_R1 = sprintf("%.2f", $percent_A ) . "%";
  my $format_percent_C_R1 = sprintf("%.2f", $percent_C ) . "%";
  my $format_percent_G_R1 = sprintf("%.2f", $percent_G ) . "%";
  my $format_percent_T_R1 = sprintf("%.2f", $percent_T ) . "%";
  my $format_percent_N_R1 = sprintf("%.4f", $percent_N ) . "%";
  
  my $format_percent_GC_R2 = sprintf("%.2f", $percent_GC_R2 ) . "%";
  my $format_percent_A_R2 = sprintf("%.2f", $percent_A_R2 ) . "%";
  my $format_percent_C_R2 = sprintf("%.2f", $percent_C_R2 ) . "%";
  my $format_percent_G_R2 = sprintf("%.2f", $percent_G_R2 ) . "%";
  my $format_percent_T_R2 = sprintf("%.2f", $percent_T_R2 ) . "%";
  my $format_percent_N_R2 = sprintf("%.4f", $percent_N_R2 ) . "%";
 
  my $format_percent_GC_R1R2 = sprintf("%.2f", $percent_GC_R1R2 ) . "%";
  my $format_percent_A_R1R2 = sprintf("%.2f", $percent_A_R1R2 ) . "%";
  my $format_percent_C_R1R2 = sprintf("%.2f", $percent_C_R1R2 ) . "%";
  my $format_percent_G_R1R2 = sprintf("%.2f", $percent_G_R1R2 ) . "%";
  my $format_percent_T_R1R2 = sprintf("%.2f", $percent_T_R1R2 ) . "%";
  my $format_percent_N_R1R2 = sprintf("%.4f", $percent_N_R1R2 ) . "%";

  my $format_adaptors = sprintf("%.3f", $bestAdapterhit_percent_R1R2 ) . "%";  
    
  print "File (R1):     $analysispath_R1\n";
  print "File (R2):     $analysispath_R2\n";
  
  if($count_seq == $count_seq_R2){
    print "Nr read pairs: $count_seq ($format_count_seq_pairs)\n";
  }
  else {
    print "Nr read pairs: UNMATCHED COUNT! R1=$count_seq ($format_count_seq_R1) R2=$count_seq_R2 ($format_count_seq_R2) \n";
  }

  print "Nr reads:      $count_seq + $count_seq_R2 = $count_seq_R1R2 ($format_count_seq_R1R2)\n";
  print "Nr bases:      $count_bases ($format_count_bases_R1) + $count_bases_R2 ($format_count_bases_R2) = $count_bases_R1R2 ($format_count_bases_R1R2)\n"; 
  
  print "bases >= Q30:  $count_Q30 ($format_count_Q30bases_R1,  $format_percent_Q30_R1) + $count_Q30_R2 ($format_count_Q30bases_R2,  $format_percent_Q30_R2) = $count_Q30_R1R2 ($format_count_Q30bases_R1R2,  $format_percent_Q30_R1R2)\n";
  print "bases <= Q20:  $count_Q20 ($format_count_Q20bases_R1,  $format_percent_Q20_R1) + $count_Q20_R2 ($format_count_Q20bases_R2,  $format_percent_Q20_R2) = $count_Q20_R1R2 ($format_count_Q20bases_R1R2,  $format_percent_Q20_R1R2)\n";

  print "GC content:    R1: $format_percent_GC_R1 R2: $format_percent_GC_R2 Total: $format_percent_GC_R1R2\n";

  if(&arg_hasflag("verbose")){
   print "bases R1:      A = $format_percent_A_R1, C = $format_percent_C_R1, G = $format_percent_G_R1, T = $format_percent_T_R1\n";
   print "bases R2:      A = $format_percent_A_R2, C = $format_percent_C_R2, G = $format_percent_G_R2, T = $format_percent_T_R2\n";
   print "bases Total:   A = $format_percent_A_R1R2, C = $format_percent_C_R1R2, G = $format_percent_G_R1R2, T = $format_percent_T_R1R2\n";
  }
  
  print "'N' bases:     R1: $count_N ($format_percent_N_R1) R2: $count_N_R2 ($format_percent_N_R2) Total: $count_N_R1R2 ($format_percent_N_R1R2)\n"; 
  
  if(&arg_hasflag("verbose") ){
  print "Seqlength R1:  $minreadlength-$maxreadlength (average = $format_avlen)\n";
  print "Seqlength R2:  $minreadlength_R2-$maxreadlength_R2 (average = $format_avlen_R2)\n";
  print "Seqlength Tot: $minreadlength_R1R2-$maxreadlength_R1R2 (average = $format_avlen_R1R2)\n\n";
  }
  else{
  print "Seqlength:     R1: $minreadlength-$maxreadlength (av. = $format_avlen) R2: $minreadlength_R2-$maxreadlength_R2 (av. = $format_avlen_R2) Total average = $format_avlen_R1R2\n";
  }

  if($best16Shit_percent>$best16Shit_threshold){
  print "Species:       $best16Shit\n";

  if($best16Shit_cov ne "NA"){
    $cov = sprintf("%.1f", $best16Shit_cov ) . "X"; 
    $covQ30 = sprintf("%.1f", $best16Shit_covQ30 ) . "X"; 
    print "Coverage:      $cov\n";
    print "Cov. (Q30+):   $covQ30\n";  
   }
  }
  print "Adaptors:      $format_adaptors\n";
    
  ###############
  print "\nDiscriminating 16S sequence-tags detected in $count_initiallyanalysis R1-reads: \n";
  print "$all16shits_formated_lines";
  
  #### VERBOSE -NON-TABULAR OUTPUT - also print non-unique 16S hits
  
   if( &arg_hasflag("verbose") ){
    
     print "\nNon-discriminating 16S sequence-tags detected in $count_initiallyanalysis R1-reads: \n";
     print "$all16shits_nodiscr_formated_lines";
      
   } # end verbose
   ################
  
  print "\nAdaptor sequences detected in $count_initiallyanalysis  R1-reads: \n"; 
  print "$allAdapterhits_formated_lines";
 
  print "\nAdaptor sequences detected in $count_initiallyanalysis_R2  R2-reads: \n";
  print "$allAdapterhits_formated_lines_R2"; 
  
  print "-------------------\n";
  }  # end if pair==0
  
  

  
 } ### END of foreach analysisfile   
} # end sub_fileinfo


######################################### END: sub_fileinfo_fastq ####################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================





#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\








#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
######################################### sub_fileinfo_fasta #########################################

sub sub_fileinfo_fasta(){

  ##############
  #versioninfo
  ##############
  
  my $sub_name = "fileinfo_fasta";
  my $sub_version = "2022.09.16";
 
  my $subfunction_arg = shift(@_);
  if($subfunction_arg eq "version") { print "    component $sub_name version=$sub_version\n";  return; } 
  
  #############
  #arguments
  #############
  
  my $path = &arg_getremaining(0);  # this is the filename or the directoryname
       
  my $verbose = 0;  # print out more info 
  if( &arg_hasflag("verbose") ){
     $verbose = 1; 
  }
   
  my $tabular = 0; # print out info as a tab delimited table
  
  if(-d $path){
     $tabular = 1; # if input is directory, tabular is default
  }
  
  if( &arg_hasflag("tabular") ){
     $tabular = 1; 
  }
  if( &arg_hasflag("notabular") ){
     $tabular = 0; 
  }   
  my $printheader = 1; # do not print headerrow (only used in tabular mode)
  if( &arg_hasflag("noheader") ){
     $printheader = 0; 
  }
  
  my $recursive = 0; # recursively scan subdirs  -NOt implemented yet
  if( &arg_hasflag("recursive") ){
     $recursive = 1; 
  }
 
  my $overwrite = 0; # overwrite outputfile if it already exists   -NOt implemented yet
  if( &arg_hasflag("overwrite") ){
     $overwrite = 1; 
  }

  my $append = 0; # append to outputfile if it already exists   -NOt implemented yet
  if( &arg_hasflag("append") ){
     $append = 1; 
  }
    
  my $outputfile = ""; # output filename ... default is to prit to STDOUT   -NOt implemented yet
  if( &arg_hasparameter("output") ){
      $outputfile = &arg_getparameter("output");
  }
  
  if(-e $outputfile && $overwrite==0 && $append ==0 ){
    die("Output filename already exists. Use another name OR use --overwrite OR use --append");
  }
  if(-d $outputfile  ){
    die("Output filename is used by a directory. Use another name");
  }      
  

    
  &arg_warn_if_unused_argument();
 
  ###################################
  # prepare list of files to analyze
  ###################################
  
 
  #################
  # set input mode
  #################
  
  my $inputmode = 0;  # 1=single file 2=directory 3=recursive scan of directory
  
  if(-f $path){
    $inputmode = 1;
  }
 
  if(-d $path){
  
    if($recursive==1) {  
      $inputmode = 3;
    }
    else{
      $inputmode = 2;
    }
  }
  
  if($inputmode==0){
    die( "input path is not a file or directory");
  }
  

   
  ########################
  # prepare list of files
  ########################
  
  my @listoffiles = ();
    
  if($inputmode == 1){
     push(@listoffiles, $path);
  }
  
  if($inputmode == 2){
      &private_pushfilesindir($path,"fasta",\@listoffiles);
  } 
  if($inputmode == 3){
      &private_pushfilesindir($path,"fasta",\@listoffiles,"recursive");
  }  
 
  
  ########################
  # run aanlysis
  ########################   
  
  $filecounter = 0;
   
  foreach $analysispath (@listoffiles) {
   
   
   &private_filereader_open($analysispath,"FA1");
   
   $filecounter++;
   
   ## only print header first time
   if( $tabular==1  && $filecounter>1){
     $printheader = 0;
     $append = 1;
   } 
            
   my $count_seq = 0;   
   my $count_bases = 0;
   my $count_A = 0;
   my $count_C = 0;
   my $count_G = 0;
   my $count_T = 0;
   my $count_N = 0;
           
   my @specieshits = ();
   for(my $i=0;$i<@constants_species_16Sseq;$i++){
    push(@specieshits ,0);    
   }

   my $line = "";  
   my $docontinue = "\n";  
   my $part = "";
   my $currenttag = "";
   my $currentseq = "";    
   
   my @seqlengths = ();   
   my @seqcoverages = ();   
      
   while(defined $docontinue){


     $line = &private_filereader_Getline("FA1");
     $docontinue = $line;
     chop($line);
     
     if( substr($line,0,1) eq ">" ||  !defined $line){
     
         if(defined($line)){
           $count_seq++;
         }
         
         ## Spades style coverageinfo
         if(index($line,"_cov_")>=0){
            $cov = substr($line,index($line,"_cov_")+5);
            push(@seqcoverages,$cov);
         } 
            
             
         #print last sub-sequence, if exists
         if($count_seq>1  || ($count_seq==1 && !defined $line)){
                      
            my $seqlen =  length($currentseq);
            push(@seqlengths,$seqlen);

            


            $count_bases = $count_bases + $seqlen;
            my $nrA  = $currentseq =~ tr/aA/aA/; 
            $count_A = $count_A + $nrA;
            my $nrC  = $currentseq =~ tr/cC/cC/; 
            $count_C = $count_C+ $nrC;
            my $nrG  = $currentseq =~ tr/gG/gG/; 
            $count_G = $count_G + $nrG;
            my $nrT  = $currentseq =~ tr/tT/tT/; 
            $count_T = $count_T+ $nrT;
            my $nrN  = $currentseq =~ tr/nN/nN/; 
            $count_N = $count_N+ $nrN;
                     
                     
            ####### 16S analysis ##########3
            #loop species
            for(my $i=0;$i<@constants_species_16Sseq;$i++){
              my $seq16s = $constants_species_16Sseq[$i];
              
              #make tags  200 bases- step 10 bases per tag
              for(my $j=0;$j<length($seq16s)-200;$j=$j+10){
                 
                $seqtag = substr($seq16s,$j,200);
                $seqtag_rc = private_reverscomplement($seqtag);
                 
                my $hit_index = -1; # marks no-hit

		 for(my $k=0;$k<@constants_species_16Sseq;$k++){  # step through all species and make sure the seqtag is unique (discriminatory for its species)
		   if($k != $i && ( index($seqtag,$constants_species_16Sseq[$k])>=0 || index($seqtag_rc,$constants_species_16Sseq[$k])>=0 )){
		     $hit_index = -2; # seqtag exists in other species => mark as non-discriminatory seqtag
		   }
		 }

                 if(index($currentseq,$seqtag)>=0 && $hit_index==-1 ){
                       $hit_index = $i;
                 }
                 elsif(index($currentseq,$seqtag_rc)>=0 && $hit_index==-1){
                       $hit_index = $i;
                 }
                 
                 if( $hit_index>=0){
                     $specieshits[$hit_index] =  $specieshits[$hit_index] +1;
                 }    
                
              } #end for each seqtag
              
            } # end foreach species
            #### END 16S analysis

         
            ### Adaptor analysis
            for(my $i=0;$i<@constants_adapters;$i++){
         
              if(index($currentseq, $constants_adapters_seq[$i])>=0){
                 $adaptorhits_full[$i] =  $adaptorhits_full[$i] +1; 
              }
            
            } # end for all adaptors
            #### END adaptor analysis
            
	 }
	 
         $currentseq = "";  # clear current seq
       }
       else{

          $line =~ s/\s+//g;  #remove spaces
          $line =~ s/[0-9]//g;  #remove digits
       
          $currentseq = $currentseq . $line ;
         
       }

   } ## END While defined $line

   
          
  &private_filereader_close("FA1");

 
  ############################
  #Process data for printout
  ############################
  
  
  #N50 calculation
  
  my @seqlengths_sorted = sort  { $b <=> $a } @seqlengths;

  my $cumlen = 0;
  my $N50 = 0;
  my $L50 = 0;
  my $N90 = 0;
  my $L90 = 0;
  
  my $contigs_over_500 = 0;
  my $contigs_over_1000 = 0;
  my $contigs_over_500_size = 0;
  my $contigs_over_1000_size = 0;
        
  for(my $i =0; $i<@seqlengths_sorted ; $i++){
     $contiglen = $seqlengths_sorted[$i];
     $cumlen = $cumlen  + $contiglen;
     
     if($contiglen>=500){
       $contigs_over_500++;
       $contigs_over_500_size =  $contigs_over_500_size + $contiglen;
     }

     if($contiglen>=1000){
       $contigs_over_1000++;
       $contigs_over_1000_size =  $contigs_over_1000_size + $contiglen;       
     }
          
     if($cumlen>=$count_bases*0.50 && $N50 ==0){
       $N50 =  $contiglen;
       $L50 = $i+1;
     }
     if($cumlen>=$count_bases*0.90 && $N90 ==0){
       $N90 =  $contiglen;
       $L90 = $i+1;
     }
  }

  #Max-min
  my $minlength = $seqlengths_sorted[-1];
  my $maxlength = $seqlengths_sorted[0];
  
  my $format_count_seq   = &private_formatstring_readcount($count_seq);
  my $format_count_bases = &private_formatstring_basecount($count_bases);

  my $avlen = 0;
  if($count_seq>0){
     $avlen = $count_bases/$count_seq;
  }
  my $format_avlen = sprintf("%.2f", $avlen );

  ##coverages

  $lowcoverage_count_format = "NA";
  $lowcoverage_cumsize_format = "NA";
  $lowcoverage_cumsize_format_p = "NA";
  $filtered_genomesize = "NA";
  $filtered_genomesize_format = "NA";
  $averagekmercoverage_format = "NA";
    
  $lowcoverage_count = 0;
  $lowcoverage_cumsize = 0; 
  $average_cov = 0;
  
  $tmp =@seqlengths;
  $tmp2=@seqcoverages;
    
   if( @seqlengths ==  @seqcoverages){   # coverageinfo exists
      

      my $cumcov = 0;
      my $cumlen = 0;   
  
      for(my $i=0; $i<@seqlengths;$i++){
         $cumcov = $cumcov + $seqlengths[$i]*$seqcoverages[$i];
         $cumlen = $cumlen + $seqlengths[$i];
      }
  
      $average_cov = 0;
      if($cumlen>0){  
         $average_cov = $cumcov/$cumlen;
      }
       $averagekmercoverage_format = sprintf("%.2f", $average_cov);
      @sequence_coverage_percent = ();
  
      for(my $i=0; $i<@seqlengths;$i++){
        if($average_cov>0){
          my $percent = $seqcoverages[$i]/$average_cov *100;
          push( @sequence_coverage_percent, $percent);
          if($percent<10){
             $lowcoverage_count++;
             $lowcoverage_cumsize = $lowcoverage_cumsize  + $seqlengths[$i];
          }
        }
        else{  
          push( @sequence_coverage_percent, 0); 
        }
      }

     $lowcoverage_count_format = $lowcoverage_count;
     $ass_perc = sprintf("%.4f", $lowcoverage_cumsize / $count_bases*100);
     $lowcoverage_cumsize_format_p = $ass_perc;
     $lowcoverage_cumsize_format = $lowcoverage_cumsize;
     
     $filtered_genomesize =  $count_bases-$lowcoverage_cumsize;
     $filtered_genomesize_format = &private_formatstring_basecount($filtered_genomesize);
   } # end if coveragedata exists
   
  ##GC content
  my $gc = 0;
  my $pa = 0;
  my $pc = 0;
  my $pg = 0;
  my $pt = 0;
  my $pn = 0;
  
  if(($count_G + $count_C + $count_T + $count_A) != 0){
    $gc = ($count_G + $count_C)/ ($count_G + $count_C + $count_T + $count_A) *100.0;
    $pa = ($count_A)/ ($count_G + $count_C + $count_T + $count_A ) *100.0;
    $pc = ($count_C)/ ($count_G + $count_C + $count_T + $count_A ) *100.0;
    $pg = ($count_G)/ ($count_G + $count_C + $count_T + $count_A ) *100.0;
    $pt = ($count_T)/ ($count_G + $count_C + $count_T + $count_A ) *100.0;
    $pn = ($count_N)/ ($count_G + $count_C + $count_T + $count_A + $count_N) *100.0;  
  }

  my $format_gc = sprintf("%.2f", $gc ) . "%";
  my $format_pa = sprintf("%.2f", $pa ) . "%";
  my $format_pc = sprintf("%.2f", $pc ) . "%";
  my $format_pg = sprintf("%.2f", $pg ) . "%";
  my $format_pt = sprintf("%.2f", $pt ) . "%";
  my $format_pn = sprintf("%.4f", $pn ) . "%";

  my $abspath = Cwd::abs_path($analysispath);
  my $filename = &private_getlastpartofpath($analysispath);
  
  ###16S list

  my %tobeprinted_hash = ();
  my %tobeprinted_hash_tab = ();
    
  my $bestspecies = "";
  my $score = "";
  my $allspecies_tab = "";
  my $allspecies_rows = "";
  
     for(my $i=0;$i<@constants_species_16Sseq;$i++){
       if($specieshits[$i]>0) {
          $tobeprinted_hash{ $constants_species[$i] . " (score=" . $specieshits[$i] . ")\n" } =  $specieshits[$i]  ;
          $tobeprinted_hash_tab { $constants_species[$i] } =  $specieshits[$i]  ;
       }
   }
  
   foreach my $name (sort { $tobeprinted_hash{$b} <=> $tobeprinted_hash{$a} } keys %tobeprinted_hash) {
    $allspecies_rows= $allspecies_rows . $name;
   }

   my $placeinlist = 1;
   foreach my $name (sort { $tobeprinted_hash_tab{$b} <=> $tobeprinted_hash_tab{$a} } keys %tobeprinted_hash_tab) {
     my $thescore =$tobeprinted_hash_tab{$name};
     if($placeinlist==1){
         $allspecies_tab=  $name . "\t" .  $thescore . "\t" ;
     }
     $allspecies_tab=  $allspecies_tab . "$name (score=$thescore):"; ;
     
   }
   
   my @keys = keys %tobeprinted_hash;
   
   if(@keys==0){
     $allspecies_rows = "none\n";
     $allspecies_tab = "none\tNA\t";

   }
   
   ## Adapters
   
   %tobeprinted_hash = ();
   
   my $adaptersrow = "";
   my $adaptersdetected = "No";
   
    
   for(my $i=0;$i<@constants_adapters_seq;$i++){
        if($adaptorhits_full[$i]>0) {
          $tobeprinted_hash{ $constants_adapters[$i] . " hits=" . $adaptorhits_full[$i]  . "\n" } =  $adaptorhits_full[$i]  ;
        }
   }

   my @keys = keys %tobeprinted_hash;
   
   if(@keys>0){ 
     $adaptersrow =  "\nAdaptor sequences detected: \n";  
     $adaptersdetected = "yes";   
       
     foreach my $name (sort { $tobeprinted_hash{$b} <=> $tobeprinted_hash{$a} } keys %tobeprinted_hash) {
       $adaptersrow = $adaptersrow .  $name;
     }
   }
   
  ######## TABULAR OUTPUT ######################################
  if($tabular){
   if($printheader==1){  #  print header
     if($verbose==1){
       print "file\tpath\tcontigs\tbases\tgc_percent\ta_percent\tc_percent\tg_percent\tt_percent\tn_percent\tshortest_contig\tlongest_contig\taveragelen\tN50\tL50\tN90\tL90\tcontigs_larger500\tassemblysize_larger500\tcontigs_larger1000\tassemblysize_larger1000\tspecies\tspeciesscore\tall_speciecas_detected\tadaptors_detected\tlowcoverage_contigs\tlowcoverage_size\tlowcoverage_percent\tFiltered_genomesize\taveragecoverage\n";
     }
     else{
       print "file\tpath\tcontigs\tbases\tgc_percent\tn_percent\tshortest_contig\tlongest_contig\taveragelen\tN50\tL50\tN90\tL90\tcontigs_larger500\tassemblysize_larger500\tcontigs_larger1000\tassemblysize_larger1000\tspecies\tspeciesscore\tall_speciecas_detected\tadaptors_detected\tlowcoverage_contigs\tlowcoverage_size\tlowcoverage_percent\tFiltered_genomesize\taveragecoverage\n";
     }
   }
   if($verbose==1){
      print "$filename\t$abspath\t$count_seq\t$count_bases\t$gc\t$pa\t$pc\t$pg\t$pt\t$pn\t$minlength\t$maxlength\t$format_avlen\t$N50\t$L50\t$N90\t$L90\t$contigs_over_500\t$contigs_over_500_size\t$contigs_over_1000\t$contigs_over_1000_size\t$allspecies_tab\t$adaptersdetected\t$lowcoverage_count_format\t$lowcoverage_cumsize_format\t$lowcoverage_cumsize_format_p\t$filtered_genomesize\t$averagekmercoverage_format\n";
   }
   else{
      print "$filename\t$abspath\t$count_seq\t$count_bases\t$gc\t$pn\t$minlength\t$maxlength\t$format_avlen\t$N50\t$L50\t$N90\t$L90\t$contigs_over_500\t$contigs_over_500_size\t$contigs_over_1000\t$contigs_over_1000_size\t$allspecies_tab\t$adaptersdetected\t$lowcoverage_count_format\t$lowcoverage_cumsize_format\t$lowcoverage_cumsize_format_p\t$filtered_genomesize\t$averagekmercoverage_format\n";
   }   
  }
  
  
  ######### NON-TABULAR OUTPUT ################################
  else{
   print "File:         $filename\n";
   print "Path:         $abspath\n";
   print "Nr sequences: $count_seq\n";
   print "Nr bases:     $count_bases ($format_count_bases)\n"; 
   print "GC content:   $format_gc\n";
   print "bases:        A = $format_pa, C = $format_pc, G = $format_pg, T = $format_pa\n";
   print "'N' bases:    $count_N ($format_pn) \n"; 
   print "Seqlength:    $minlength-$maxlength (average = $format_avlen)\n";
   print "N50:          $N50\n";
   print "L50:          $L50\n";
   print "N90:          $N90\n";
   print "L90:          $L90\n";
   print "Contigs 500+  $contigs_over_500\n";
   print "Assembly 500+ $contigs_over_500_size\n";
   print "Contigs 1kb+  $contigs_over_1000\n";
   print "Assembly 1kb+ $contigs_over_1000_size\n";

  if($lowcoverage_count_format ne "NA"){   
   print "Low-cov cont  $lowcoverage_count_format\n";
   print "Low-cov size  $lowcoverage_cumsize_format\n";
   print "Low-cov perc  $lowcoverage_cumsize_format_p\n";
   print "Filt. size    $filtered_genomesize ($filtered_genomesize_format)\n";
   print "Av kmer cov   $averagekmercoverage_format\n";
  }
  
   ###############
   print "\nBest 16S sequence matches: \n";

   print $allspecies_rows;
  
  
  print $adaptersrow;
  
  
   ################

  

   
  } # end else (not tabular output)
  #### END NON-TABULAR
  
 } ### END of foreach analysisfile 
} # end sub_fileinfo



######################################### END: sub_fileinfo_fasta ####################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================





#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\








#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
######################################### sub_downsample #############################################


sub sub_downsample(){

 my $sub_name = "downsample";
 my $sub_version = "2022.03.04";
 
 my $subfunction_arg = shift(@_);
 if($subfunction_arg eq "version") { print "sub $sub_name version=$sub_version\n"; return; } 
 
   my $path1 = &arg_getremaining(0);
   my $path2 = &arg_getremaining(1);

   ### Determin wanted bases
   my $wanted_bases = 0;
   my $coverage = 0;
   my $genomesize = 0;
   
   if( &arg_hasparameter("bases") ){
      $wanted_bases = &arg_getparameter("bases");
   }
   elsif( &arg_hasparameter("coverage") && &arg_hasparameter("genomesize")){
      $coverage = &arg_getparameter("coverage");
      $genomesize = &arg_getparameter("genomesize");
      $wanted_bases = int($genomesize*$coverage);
   }
   elsif( &arg_hasparameter("coverage") && &arg_hasparameter("organism")){
      $coverage = &arg_getparameter("coverage");
      $genomesize=&private_gettypicalgenomesize(&arg_getparameter("organism"));    
      $wanted_bases = int($genomesize*$coverage);
   }

   if($wanted_bases==0){
      die("Error: Not enough information. bases or coverage and genomesize/organism must be given");
   }
  
   &arg_warn_if_unused_argument();
     
   #### R1 output
   ## default is to ad a d-tag
   my $outpath1 = &addtagtopath($path1,".d");
   
   ## overide default if  user has specified two output files
   if( &arg_hasparameter("out1") && &arg_hasparameter("out2")  ){
      $outpath1 = &arg_getparameter("out1");
   }
   ## overide default if  user has specified one output files and only given one inputfile
   if( &arg_hasparameter("out") && $path2 eq ""){
      $outpath1 = &arg_getparameter("out");
   }   
   ## overide default if  user has specified one output files and only given one inputfile
   if( &arg_hasparameter("out1") && $path2 eq ""){
      $outpath1 = &arg_getparameter("out1");
   }   

   &private_filereader_open($path1,"FQIN1");
   &private_filewriter_open($outpath1,"FQOUT1");
   
   #### R2 output
   my $outpath2 = "";
   if(-f $path2){

      ## default is to ad a d-tag
      $outpath2 = &addtagtopath($path2,".d");
      
      ## overide default if  user has specified two output files
      if( &arg_hasparameter("out1") && &arg_hasparameter("out2") ){
         $outpath2 = &arg_getparameter("out2");
      }
      
      &private_filereader_open($path2,"FQIN2");
      &private_filewriter_open($outpath2,"FQOUT2");
   }
   
   ### Mode
   $mode = 1; # count all bases
   if( &arg_hasflag("q30bases") ){
     $mode = 2; # count Q30+ bases
   }
   
   
   my $count_seq = 0;
   my $count_bases = 0;
   my $count_Q30 = 0;
   my @fastqseq1 = ();
 
   do{
     
     @fastqseq1 = &private_filereader_getlines(4,"FQIN1");
     
     if(@fastqseq1==4){      
       $count_seq++;

       $count_bases = $count_bases + length($fastqseq1[1]);
       
       my $nrQ30plus  = $fastqseq1[3] =~ tr/?@ABCDEFGHI/?@ABCDEFGHI/; 
       $count_Q30 = $count_Q30 + $nrQ30plus;
       &private_filewriter_writeline_array(\@fastqseq1,"FQOUT1");
       

       if($outpath2 ne ""){
         
         my @fastqseq2 = &private_filereader_getlines(4,"FQIN2");
       
         if(@fastqseq2==4){      

           $count_bases = $count_bases + length($fastqseq2[1]);
       
           my $nrQ30plus  = $fastqseq2[3] =~ tr/?@ABCDEFGHI/?@ABCDEFGHI/; 
           $count_Q30 = $count_Q30 + $nrQ30plus;
           &private_filewriter_writeline_array(\@fastqseq2,"FQOUT2");
         }
       }


       if($mode ==1 && $count_bases>=$wanted_bases){
         @fastqseq1 = ();  #break the loop
       }
       elsif($mode ==2 && $count_Q30>=$wanted_bases){
         @fastqseq1 = ();  #break the loop
       }
     }
  }while(@fastqseq1==4);

  &private_filereader_close("FQIN1");
  &private_filewriter_close("FQOUT1");
  
  if(-f $path2){
    &private_filereader_close("FQIN2");
    &private_filewriter_close("FQOUT2");
  }

 
 ## WARNING message if wanted bases was not obtained 
 if($mode ==1 && $count_bases<$wanted_bases){ 
    my $percent = $count_bases/$wanted_bases*100;
    my $format_percent = sprintf("%.2f", $percent ) . "%";
    
    my $format_reached_coverage = "";
    if($coverage>0){
      my $reached_coverage = $coverage*$count_bases/$wanted_bases;
      $format_reached_coverage =  " (" . sprintf("%.2f", $reached_coverage ) . "X)";
    }
    print "WARNING: not enough bases. Reached percent of goal: $percent" . "$format_reached_coverage\n";
 }
 elsif($mode ==2 && $count_Q30<$wanted_bases){ 
    my $percent = $count_bases/$wanted_bases*100;
    my $format_percent = sprintf("%.2f", $percent ) . "%";    
    
    my $format_reached_coverage = "";
    if($coverage>0){
      my $reached_coverage = $coverage*$count_bases/$wanted_bases;
      $format_reached_coverage =  " (" . sprintf("%.2f", $reached_coverage ) . "X)";
    }
    print "WARNING: not enough bases. Reached percent of goal: $percent" . "$format_reached_coverage\n";
 }


}
######################################### END: sub_downsample ####################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================





#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\









#
#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
######################################### sub_filtercontigs ##########################################

sub sub_filtercontigs{

 my $sub_name = "filtercontigs";
 my $sub_version = "2022.03.04";
 
 my $subfunction_arg = shift(@_);
 if($subfunction_arg eq "version") { print "sub $sub_name version=$sub_version\n"; return; } 
 

   my $path1 = &arg_getremaining(0);

   # default values
   my $coveragecutoff = 10;
   my $sizecutof =200;
   
   ###################
   # override defaults
   if( &arg_hasparameter("coveragethreshold") ){
      my $value =&arg_getparameter("coveragethreshold");
      if(substr($value,-1) eq "%"){
         $value = substr($value,0,-1);
      }
      if($value>0){
         $coveragecutoff = $value;
      }
   }
   if( &arg_hasparameter("sizethreshold") ){
      my $value =&arg_getparameter("sizethreshold");
      if($value>=0){
         $sizecutof = $value;
      }
   }  
   
   # flags
   my $verbose = 0; 
   if( &arg_hasflag("verbose") ){
     $verbose = 1; 
   }
   
   &arg_warn_if_unused_argument();
      
      
   #################3
   ## outfile
   my $outpath1 = &addtagtopath($path1,".f");

   &private_filewriter_open($outpath1,"FAOUT1");
   &private_filereader_open($path1,"FA1");

   
   print "filtering contigs in $path1 with coveragethreshold=$coveragecutoff% and sizethreshold=$sizecutof" . "bp\n";
   
   my $count_seq = 0;
   my $count_bases = 0;
   
   my @sequence_names =  ();
   my @sequence_lengths =  ();
   my @sequence_coverage =  ();
   my $current_index = -1;

   my $warnings = 0;
   
   ################################
   #### Determin average coverage
   do{

       $line = &private_filereader_getline("FA1");
       $docontinue = $line;
       chop($line);
       
       #Nametag
       if(substr($line,0,1) eq ">"){

          $count_seq++;     
          $cov = 0;
          
          #spades-style coverage info      
          if(index($line,"_cov_")>=0){
             $cov = substr($line,index($line,"_cov_")+5);
          } 
          
          # No coverage info found
          if($cov==0){
            $warnings++;
          } 
                     
          push(@sequence_names, $line);
          push(@sequence_lengths, 0);
          push(@sequence_coverage,$cov);


          # a previous sequence exists (lengthdata collected)
          if($count_seq>1){                   
              my $seqlen =  length($currentseq);        
              $sequence_lengths[$current_index] = $seqlen;
	  }
	  
	  
	  $current_index++;  #index indicating where the nextcomming seq-length info should be stored in the array (one below sequence count)
          $currentseq = "";
       }
       else{

          $line =~ s/\s+//g;  #remove spaces
          $line =~ s/[0-9]//g;  #remove digits
          $currentseq = $currentseq . $line ;
          
       }
       
 
  }while($docontinue ne undef);

            
  my $seqlen =  length($currentseq);       
  $sequence_lengths[$current_index] = $seqlen;
               
  &private_filereader_close("FA1"); 

  if($warnings>0){
     print "WARNING: no coverage info found in $warnings outof $count_seq contig  No filtering based on coverage will be made\n"; 
  }
  
  my $cumcov = 0;
  my $cumlen = 0;   
  for(my $i=0; $i<@sequence_lengths;$i++){
    $cumcov = $cumcov + $sequence_lengths[$i]*$sequence_coverage[$i];
    $cumlen = $cumlen + $sequence_lengths[$i];
  }
  
  my $average_cov = 0;
  if($cumlen>0){  
    $average_cov = $cumcov/$cumlen;
  }
  
  @sequence_coverage_percent = ();
  
  for(my $i=0; $i<@sequence_lengths;$i++){
     if($average_cov>0){
       push( @sequence_coverage_percent, $sequence_coverage[$i]/$average_cov *100);
     }
     else{  push( @sequence_coverage_percent, 0); }
  }
    
  
  ##########################################
  ### read file again -- pass2 = do filtering 
   
  &private_filereader_open($path1,"FA1");

   $count_seq=0;
   my $written_seq = 0;
   my $written_bases = 0;
   my $keepseq = 0;
   
   do{

     $line = &private_filereader_getline("FA1");
     $docontinue = $line;
     chop($line);
     
       if(substr($line,0,1) eq ">"){

          my $cov_perc = $sequence_coverage_percent[$count_seq];
          my $cov = $sequence_coverage[$count_seq];
          my $seq_len = $sequence_lengths[$count_seq];
          my $name = $sequence_names[$i];

          $keepseq = 0;
          
          ##KEEP
          if( ($cov_perc>=$coveragecutoff || $cov_perc==0 ) && $seq_len >= $sizecutof){
            $keepseq =1;
            $written_seq++;
            $written_bases =  $written_bases+$seq_len;
            if($verbose){
              print "#$count_seq\tKEEP\tcum_len=$written_bases\tcontig_len=$seq_len\tcontig_cov=$cov ($cov_perc%)\n";
            }
          }
          ## FILTER
          else
          {
           if($verbose){
             print "#$count_seq\tFILTER\tcum_len=$written_bases\tcontig_len=$seq_len\tcontig_cov=$cov ($cov_perc%)\n";
           }
          }

          $count_seq++;  #counter
          
          ## WRITE IF KEEP
          if($keepseq==1){
             &private_filewriter_writeline($line,"FAOUT1");
          }


       }
       else{
       
          ## WRITE IF KEEP
          if($keepseq==1){
             &private_filewriter_writeline($line,"FAOUT1");
          }
         
       }
       
 
     
  }while($docontinue ne undef);

  
          
          
  &private_filereader_close("FA1");
  &private_filewriter_close("FAOUT1");

  my $filtered_seq = @sequence_lengths - $written_seq;
  my $filtered_bases = $cumlen -  $written_bases;
  print "Size before filtering = $cumlen\nSize after filtering = $written_bases\n";
  if($average_cov>0){
  print "Average coverage = $average_cov\n";
  }
  print "Filtered sequences = $filtered_seq\nFiltered bases = $filtered_bases\n";

}
######################################### END: sub_filtercontigs ####################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================







#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\









#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
#####################################  sub_kmeroverlapp ##############################################

sub sub_kmeroverlapp{

 my $sub_name = "sub_kmeroverlapp";
 my $sub_version = "2022.09.08";
 
 my $subfunction_arg = shift(@_);
 if($subfunction_arg eq "version") { print "sub $sub_name version=$sub_version\n"; return; } 
 
 
 my $path1 = &arg_getremaining(0); # path file1
 my $path2 = &arg_getremaining(1);  #path file2
 
 my $kmersize = 20; # default

 if( &arg_hasparameter("kmersize") ){
      my $value =&arg_getparameter("kmersize");
      if($value>1){
         $kmersize = $value;
      }
      else{
         die("Error: Illegal kmersize given to sub_kmeroverlapp:$value\n");
      }
 }

 
  my $mode = 1;  #default mode   Mode 1= quantify using all kmers, mode 2= quantify using all UNIQUE kmers
  
  if( &arg_hasflag("unique") ){
     $mode = 2; 
  }
  
 &arg_warn_if_unused_argument();
    
    
 if(!( -f  $path1)  ){
   die("Error: Not a file:$path1\n");
 }
 if(!( -f  $path2)  ){
   die("Error: Not a file:$path2\n");
 }   
  

  
  ############################
  ### read kmers in assembly 1
  
  &private_filereader_open($path1,"FA1");

  my $currentseq = "";
  my $count_seq = 0;
  my %seq1_Kmers = ();
 
  do{

     $line = &private_filereader_Getline("FA1");
     $docontinue = $line;
     chop($line);

     ### Tag
     if(substr($line,0,1) eq ">"){

         $count_seq++;

         if($count_seq>1){
           
           for(my $i=0;$i< length($currentseq)-$kmersize;$i++){
              my $kmer = substr($currentseq,$i,$kmersize);
              my $rc_kmer = &private_reverscomplement($kmer);
              
              if(!exists( $seq1_Kmers{$kmer}) ) {  
                $seq1_Kmers{$kmer} = 1; 
              }
              else{
                $seq1_Kmers{$kmer} = $seq1_Kmers{$kmer} + 1; 
              }
              if(!exists( $seq1_Kmers{$rc_kmer}) ) {  
                $seq1_Kmers{$rc_kmer} = 1; 
              }
              else{
                $seq1_Kmers{$rc_kmer} = $seq1_Kmers{$rc_kmer} + 1; 
              }
           
           } # end for all kmers
           
         } # end if seq>1 (seq data avaliabe from previous seq)

         $currentseq = "";


      }  # end tag
      ### sequencepart 
      else{
   
          $line =~ s/\s+//g;  #remove spaces
          $line =~ s/[0-9]//g;  #remove digits
          $currentseq = $currentseq . uc($line) ;
       }
       
     

  }while($docontinue ne undef);

  ## process last subseq
   for(my $i=0;$i< length($currentseq)-$kmersize;$i++){
      my $kmer = substr($currentseq,$i,$kmersize);
      my $rc_kmer = &private_reverscomplement($kmer);
              
      if(!exists( $seq1_Kmers{$kmer}) ) {  
         $seq1_Kmers{$kmer} = 1; 
      }
      else{
        $seq1_Kmers{$kmer} = $seq1_Kmers{$kmer} + 1; 
      }
      if(!exists( $seq1_Kmers{$rc_kmer}) ) {  
        $seq1_Kmers{$rc_kmer} = 1; 
      }
      else{
        $seq1_Kmers{$rc_kmer} = $seq1_Kmers{$rc_kmer} + 1; 
      }
          
    } # end for all kmers  
            
  &private_filereader_close("FA1");


  ############################
  ### read kmers in assembly 2
  
  &private_filereader_open($path2,"FA2");

  $currentseq = "";
  $count_seq = 0;
  my %seq2_Kmers = ();
 
  do{

     $line = &private_filereader_Getline("FA2");
     $docontinue = $line;
     chop($line);
     
     ### Tag
     if(substr($line,0,1) eq ">"){

         $count_seq++;
         
         if($count_seq>1){
           
           for(my $i=0;$i< length($currentseq)-$kmersize;$i++){
              my $kmer = substr($currentseq,$i,$kmersize);
              my $rc_kmer = &private_reverscomplement($kmer);
              
              if(!exists( $seq2_Kmers{$kmer}) ) {  
                $seq2_Kmers{$kmer} = 1; 
              }
              else{
                $seq2_Kmers{$kmer} = $seq2_Kmers{$kmer} + 1; 
              }
              if(!exists( $seq2_Kmers{$rc_kmer}) ) {  
                $seq2_Kmers{$rc_kmer} = 1; 
              }
              else{
                $seq2_Kmers{$rc_kmer} = $seq2_Kmers{$rc_kmer} + 1; 
              }
           
           } # end for all kmers
           
         } # end if seq>1 (seq data avaliabe from previous seq)

         $currentseq = "";


      }  # end tag
      ### sequencepart 
      else{
   
          $line =~ s/\s+//g;  #remove spaces
          $line =~ s/[0-9]//g;  #remove digits
          $currentseq = $currentseq . uc($line) ;
       }
       

  }while($docontinue ne undef);

  ## process last subseq
   for(my $i=0;$i< length($currentseq)-$kmersize;$i++){
      $kmer = substr($currentseq,$i,$kmersize);
      $rc_kmer = &private_reverscomplement($kmer);
              
      if(!exists( $seq2_Kmers{$kmer}) ) {  
         $seq2_Kmers{$kmer} = 1; 
      }
      else{
        $seq2_Kmers{$kmer} = $seq2_Kmers{$kmer} + 1; 
      }
      if(!exists( $seq2_Kmers{$rc_kmer}) ) {  
        $seq2_Kmers{$rc_kmer} = 1; 
      }
      else{
        $seq2_Kmers{$rc_kmer} = $seq2_Kmers{$rc_kmer} + 1; 
      }
          
    } # end for all kmers  
            
  &private_filereader_close("FA2");  
   
  ####################
  ## Do comparisons

    
  my $Seq1_in_seq2 = 0;
  my $Seq1_notin_seq2 = 0;
   
  foreach $key (keys %seq1_Kmers){

    $count = $seq1_Kmers{$key};
    
    if(exists( $seq2_Kmers{$key}) ){
      if($mode==2){
         $Seq1_in_seq2++;
       }
       if($mode==1){
         $Seq1_in_seq2 = $Seq1_in_seq2  + $count ;
       }
    }
    else{
       if($mode==2){
         $Seq1_notin_seq2++;
       }
       if($mode==1){
         $Seq1_notin_seq2 = $Seq1_notin_seq2  + $count ;
       }
    }

  }  ## end foreach key

  
  my $Seq2_in_seq1 = 0;
  my $Seq2_notin_seq1 = 0;
  

  
   
  foreach $key (keys %seq2_Kmers){

    $count = $seq2_Kmers{$key};
    
    if(exists( $seq1_Kmers{$key}) ){
       if($mode==2){
         $Seq2_in_seq1++;
       }
       if($mode==1){
         $Seq2_in_seq1 = $Seq2_in_seq1  + $count ;
       }
    }
    else{
       if($mode==2){
         $Seq2_notin_seq1++;
       }
       if($mode==1){
         $Seq2_notin_seq1 = $Seq2_notin_seq1  + $count ;
       }
    }

  }  ## end foreach key


 my $percent_seq1_in_seq2 = $Seq1_in_seq2  / ($Seq1_in_seq2 + $Seq1_notin_seq2) *100.0;
 
 my $percent_seq2_in_seq1 = $Seq2_in_seq1  / ($Seq2_in_seq1 + $Seq2_notin_seq1) *100.0;


print "$path1\t$path2\tSeq1-kmers_in_Seq2:\t$percent_seq1_in_seq2\tSeq2-kmers_in_Seq1:\t$percent_seq2_in_seq1\n";

}







######################################### END: sub_kmeroverlapp ####################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================








#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\









#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
#####################################  sub_kmerhisto ##############################################

sub sub_kmerhisto{

 my $sub_name = "sub_kmerhisto";
 my $sub_version = "2022.03.04";
 
 my $subfunction_arg = shift(@_);
 if($subfunction_arg eq "version") { print "sub $sub_name version=$sub_version\n"; return; } 
 
 
 my $path1 = &arg_getremaining(0); # path file fastq R1
 my $path2 = &arg_getremaining(1);  #path file fastq R2
 
 my $path_ref = "";
  
 my $kmersize = 20; # default
 my $histo_bins = 500; # default
 my $maxbases = -1;

 if( &arg_hasparameter("kmersize") ){
      my $value =&arg_getparameter("kmersize");
      if($value>1){
         $kmersize = $value;
      }
      else{
         die("Error: Illegal kmersize given to sub_kmeroverlapp:$value\n");
      }
 }
 
 
  if( &arg_hasparameter("reference") ){
      $path_ref =&arg_getparameter("reference");

  }
 
 
 if( &arg_hasparameter("bins") ){
      my $value =&arg_getparameter("bins");
      if($value>1){
         $histo_bins = $value;
      }
      else{
         die("Error: Illegal number of bins given to sub_kmeroverlapp:$value\n");
      }
 }
 
  if( &arg_hasparameter("maxbases") ){
      my $value =&arg_getparameter("maxbases");
      if($value>1){
         $maxbases = $value;
      }
      else{
         die("Error: Illegal maxbases number given to sub_kmeroverlapp:$value\n");
      }
 }
 
 
  $gcbiasplot = 0;
  if( &arg_hasflag("gcbiasplot") ){
     $gcbiasplot = 1;
  }
  
  
 &arg_warn_if_unused_argument();
    
 my $mode = 1; # mode 1 = use all avaliable kmers in fastq files.. mode 2= use only kmers from reference fasta
 
 if($path_ref ne ""){
    $mode = 2;
 }
 ###############################################33


  
  ################################
  ## Read fastq k-mers
  
   &private_filereader_open($path1,"FQIN1");

   my $hasR2 = 0;
   if(-f $path2){
      &private_filereader_open($path2,"FQIN2");
      $hasR2 =1;
   }
   

   
   my $count_seq = 0;
   my $count_bases = 0;
   
   my @fastqseq1 = ();
   my @fastqseq2 = ();
   
   my $addedkmers_fastq_count = 0;
   my $addedkmers_fastq = 0;
   
   my %fastq_Kmers = ();
   
   do{
     
     @fastqseq1 = &private_filereader_getlines(4,"FQIN1");
     
     if(@fastqseq1==4){   
        
       $count_seq++;
    
       $count_bases = $count_bases + length($fastqseq1[1]);
       $seqlen = length($fastqseq1[1]);
       
       for($i=0;$i<$seqlen-$kmersize;$i++){
       
          $kmer = substr($fastqseq1[1],$i,$kmersize);
          
          if(exists( $fastq_Kmers{$kmer}) ) {  
            $fastq_Kmers{$kmer} = $fastq_Kmers{$kmer} +1;
            $addedkmers_fastq_count++;
          }
          else{
            $fastq_Kmers{$kmer} = 1; 
            $addedkmers_fastq++;
            $addedkmers_fastq_count++;
          }
            
       }  ## end for k-mers in seq
           
           
       if($hasR2 ==1){
         
         my @fastqseq2 = &private_filereader_getlines(4,"FQIN2");
       
         if(@fastqseq2==4){      

           $count_bases = $count_bases + length($fastqseq2[1]);
           $seqlen = length($fastqseq2[1]);

           for($i=0;$i< $seqlen-$kmersize;$i++){
       
              $kmer = substr($fastqseq2[1],$i,$kmersize);
              
              if(exists( $fastq_Kmers{$kmer}) ) {  
                $fastq_Kmers{$kmer} = $fastq_Kmers{$kmer} +1;
                $addedkmers_fastq_count++;
              }
              else{
                $fastq_Kmers{$kmer} = 1; 
                $addedkmers_fastq++;
                $addedkmers_fastq_count++;
              }
            
           }  ## end for k-mers in seq
           
      
         }
       } # end R2


       
     }
     
     if($maxbases>0 && $count_bases>$maxbases){
        @fastqseq1 = ();
     }
     
  }while(@fastqseq1==4);

  &private_filereader_close("FQIN1");

  if(-f $path2){
    &private_filereader_close("FQIN2");
  }


 ###### IF MODE 2 = REFERENCE BASED
 
 if($mode==2) {
 
   @histo = ();
   @histo_GCs = ();
   @histo_ATs = ();
   my $total_kmers = 0;

   for($i=0;$i<$histo_bins+1;$i++){
     push(@histo,0);
     push(@histo_GCs,0);
     push(@histo_ATs,0);
   }

   my $count_seq = 0;
   my $refsize = 0;
   my $count_refkmers = 0;
   my $count_refkmers_unique = 0;
   
   &private_filereader_open($path_ref,"FA1");

   my $currentseq = "";

   #read the reference file
   do{  

     $line = &private_filereader_getline("FA1");
     $docontinue = $line;
     chop($line);
     
     ### Tag
     if(substr($line,0,1) eq ">"){

         $count_seq++;
         
         if($count_seq>1){
           
           $refsize = $refsize + length($currentseq);
           $seqlen = length($currentseq);
           
           #extract k-mers
           for(my $i=0;$i<$seqlen-$kmersize;$i++){
           
              my $kmer = substr($currentseq,$i,$kmersize);
              my $rc_kmer = &private_reverscomplement($kmer);
              
              my $nrAT_k  = $kmer =~ tr/aAtT/aAtT/; 
              my $nrGC_k  = $kmer =~ tr/cCgG/cCgG/; 
                  
              $k_cov = 0;
              
              if(exists( $fastq_Kmers{$kmer}) ) {  
                 $k_cov = $k_cov + $fastq_Kmers{$kmer};
              }
              if(exists( $fastq_Kmers{$rc_kmer}) ) {  
                 $k_cov = $k_cov + $fastq_Kmers{$rc_kmer};
              }              
              
              if($k_cov>$histo_bins){
                $k_cov = $histo_bins;
              }
              
              $histo[$k_cov] = $histo[$k_cov] + 1;
              
              $histo_GCs[$k_cov] = $histo_GCs[$k_cov] + $nrGC_k;
              $histo_ATs[$k_cov] = $histo_ATs[$k_cov] + $nrAT_k;
              
              $total_kmers++;
             
           } # end for all kmers
           
         } # end if seq>1 (seq data avaliabe from previous seq)

         $currentseq = "";


      }  # end tag
      ### sequencepart 
      else{
   
          $line =~ s/\s+//g;  #remove spaces
          $line =~ s/[0-9]//g;  #remove digits
          
          $currentseq = $currentseq . $line ;
       }
       

  }while($docontinue ne undef);

  ## process last subseq
  
  $refsize = $refsize + length($currentseq);
  $seqlen = length($currentseq);
           
  #extract k-mers
  for(my $i=0;$i<$seqlen-$kmersize;$i++){
           
       my $kmer = substr($currentseq,$i,$kmersize);
       my $rc_kmer = &private_reverscomplement($kmer);
              
       my $nrAT_k  = $kmer =~ tr/aAtT/aAtT/; 
       my $nrGC_k  = $kmer =~ tr/cCgG/cCgG/; 
                  
       $k_cov = 0;
              
       if(exists( $fastq_Kmers{$kmer}) ) {  
          $k_cov = $k_cov + $fastq_Kmers{$kmer};
       }
       if(exists( $fastq_Kmers{$rc_kmer}) ) {  
          $k_cov = $k_cov + $fastq_Kmers{$rc_kmer};
       }              
             
       if($k_cov>$histo_bins){
         $k_cov = $histo_bins;
       }
              
       $histo[$k_cov] = $histo[$k_cov] + 1;
             
       $histo_GCs[$k_cov] = $histo_GCs[$k_cov] + $nrGC_k;
       $histo_ATs[$k_cov] = $histo_ATs[$k_cov] + $nrAT_k;
              
       $total_kmers++;
             
  } # end for all kmers
            
  &private_filereader_close("FA1");


   print "Reference: $path_ref\nContigs: $count_seq\nSize $refsize \nkmers: $total_kmers\n\n";  
   print "Fastq1: $path1\nFastq2: $path2\nkmers total: $addedkmers_fastq_count\nkmers unique: $addedkmers_fastq\n\n";
   
   print "kmer coverage histogram:\n\nkmer_coverage\tcount\n";
   
   for($i=0;$i<$histo_bins+1;$i++){
      $val = $histo[$i];
      print "$i\t$val\n";
   }   
   
   
   if($gcbiasplot ==1){

   print "kmer-count\tkmer-GC_content\n";

   for($i=0;$i<$histo_bins+1;$i++){
     $gc = $histo_GCs[$i];
     $at = $histo_ATs[$i];
     if(($gc+$at)>0){
        $percentgc = $gc/($gc+$at )*100.0;
     }
     else{
       $percentgc=0;
     }
     print "$i\t$percentgc\n";
   }

 }
 
   
 }  ## end if reference exists
 ################
 
  ###### IF MODE 1 = NO REFERENCE 
 
 if($mode==1) {
 
  
   print "Fastq1: $path1\nFastq2: $path2\nkmers total: $addedkmers_fastq_count\nkmers unique: $addedkmers_fastq\n\n";
   
   print "kmer coverage histogram:\n\nkmer_coverage\tcount\n";
 
   @histogram = ();
   @histo_GCs = ();
   @histo_ATs = ();
   
  my %RCkmers_used = ();
 
   for($i=0;$i<$histo_bins;$i++){
      push(@histogram,0);
      push(@histo_GCs,0);
      push(@histo_ATs,0);
   }
 
   foreach $key (keys %fastq_Kmers)
   {
    $kmer = $key;
    

       
    #if not already used (as a reverse complement of another kmer)
    if(!exists $RCkmers_used{$kmer}){
    
       $rc_kmer = &private_reverscomplement($kmer);
    
       $kmer_count = $fastq_Kmers{$kmer};
    
       my $nrAT_k  = $kmer =~ tr/aAtT/aAtT/; 
       my $nrGC_k  = $kmer =~ tr/cCgG/cCgG/; 
    
       # if the reverse complement exists
       if(exists $fastq_Kmers{$rc_kmer} ){
         $kmer_count = $kmer_count + $fastq_Kmers{$rc_kmer};
         $RCkmers_used{$rc_kmer} = 1;
             
         my $tmp_nrAT_k  = $rc_kmer =~ tr/aAtT/aAtT/; 
         my $tmp_nrGC_k  = $rc_kmer =~ tr/cCgG/cCgG/; 
         
         $nrAT_k  = $nrAT_k + $tmp_nrAT_k;
         $nrGC_k  = $nrGC_k + $tmp_nrGC_k; 
       
       }
       
       if($kmer_count>$histo_bins){
         $kmer_count=$histo_bins;
       }
       $histogram[$kmer_count] = $histogram[$kmer_count] +1;
       $histo_GCs[$kmer_count] = $histo_GCs[$kmer_count] + $nrGC_k;
       $histo_ATs[$kmer_count] = $histo_ATs[$kmer_count] + $nrAT_k;
           
    } # end if not already used
    

  
  }

  print "kmers in fastq: \nuniqe = $addedkmers_fastq\ntotcount = $addedkmers_fastq_count\n\n";

  my $offset = 1;
  my $offset_value = $histogram[1];
  my $continue=1;
  
  for($i=2;$i<$histo_bins;$i++){
     $value = $histogram[$i];
     if($continue==1 && $value<$offset_value){
       $offset = $i;
       $offset_value = $value;
     }
     else{
        $continue=0;
     }
  }
  
  

  for($i=0;$i<$histo_bins;$i++){
     $value = $histogram[$i];
     if($i>=$offset && $value>$maxval){
       $maxval = $value;
       $peakpos = $i;
     }
  }
  
  my $balance_below = 0;
  my $balance_above = 0;
  
  for($i=0;$i<$histo_bins;$i++){
     $value = $histogram[$i];
     if($i>=$offset && $i<$peakpos){
        $balance_below  = $balance_below + $value;
     }
     if($i>=$offset && $i>$peakpos){
        $balance_above  = $balance_above + $value;
     }
     if($i>=$offset && $i==$peakpos){
        $balance_above  = $balance_above + $value/2;
        $balance_below  = $balance_below + $value/2;
     }
  }
  
   if($peakpos>0){
   
   if($peakpos<10){
      print "WARNING peak position is at low value. Genomesize prediction may be uncertain\n";
   }
   
   my $balance = 0;
   if($balance_below+$balance_above>0){
     $balance = $balance_above/($balance_below+$balance_above);
   }
   $format_balance = int($balance*100);
   print "peak position (kmer-coverage): $peakpos";
   print "peak symetry (% kmers above peak. Should be approx 50): $balance\n";
       
   if($balance >0.6 || $balance < 0.4){
      print "WARNING peak is not symetric--estimated genome size may be wrong!!\nThis may be caused by the use of biased library prep kits such as Nextera XT\n";
   }

   $genomesizeest = int($addedkmers_fastq_count / $peakpos);
   
   print "Estimated genomesize = $genomesizeest\n(calculated as [total kmers]/[peak kmer-coverage] = $addedkmers_fastq_count /  $peakpos )\n";
 }
  
  print "kmer-count\tfrequency\n";
  
  my $peakpos = 0;
  my $maxval = 0; 

  for($i=0;$i<$histo_bins;$i++){
     $value = $histogram[$i];
     print "$i\t$value\n";
  }
 

 
 if($gcbiasplot ==1){

   print "kmer-count\tkmer-GC_content\n";

   for($i=0;$i<$histo_bins+1;$i++){
     $gc = $histo_GCs[$i];
     $at = $histo_ATs[$i];
     if(($gc+$at)>0){
        $percentgc = $gc/($gc+$at )*100.0;
     }
     else{
       $percentgc=0;
     }
     print "$i\t$percentgc\n";
   }
  } # gcbiasplot
 } # end if no reference mode mode=1
 
 }

######################################### END: sub_kmerhisto ####################################
######################################################################################################
######################################################################################################
######################################################################################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#=====================================================================================================
#=====================================================================================================








#/////////////////////////////////////////////////////////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\








 
#=====================================================================================================
#=====================================================================================================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
######################################################################################################
######################################################################################################
#####################################  COMPONENT: TOOLS ##############################################

# version 2022.03.04

######################################################################################################

sub private_assert_filepath(){

  # version 2022.03.04
  
  my $path = shift(@_);
  
  if($path eq ""){
    die("Fatal Error: No filename\n");
  }
  if(! -f $path ){
    die("Fatal Error. Not a file: $path\n");
  }

}

######################################################################################################

sub private_assert_dirpath(){

  # version 2022.03.04
  
  my $path = shift(@_);
  
  if($path eq ""){
    die("Fatal Error. No filename\n");
  }
  if(! -d $path ){
    die("Fatal Error. Not a directory: $path\n");
  }

}

######################################################################################################

sub private_assert_fileordirpath(){

  # version 2022.03.04
  
  my $path = shift(@_);
  
  if($path eq ""){
    die("Fatal Error: No filename\n");
  }
  if(! -d $path ){
    if(! -f $path ){
       die("Fatal Error. Not a file or directory: $path\n");
    }   
  }

}

######################################################################################################


sub private_isgz(){

  # version 2022.03.04
  
  my $path = lc(shift(@_));
  
  if(substr($path,-3) eq ".gz"){

    return 1;
  }
  else{
    return 0;
  }

}

######################################################################################################


sub private_isfastq(){

  # version 2022.03.04

  my $path = lc(shift(@_));

  if(substr($path,-6) eq ".fastq"){
    return 1;
  }
  elsif(substr($path,-3) eq ".fq"){      
    return 1;
  }
  elsif(substr($path,-9) eq ".fastq.gz"){
    return 1;
  }
  elsif(substr($path,-6) eq ".fq.gz"){
    return 1;
  }
  else{
    return 0;
  }

}

######################################################################################################


sub private_isfasta(){

  # version 2022.03.04
  
  my $path = lc(shift(@_));
  
  if(substr($path,-6) eq ".fasta"){
    return 1;
  }
  elsif(substr($path,-3) eq ".fa"){
    return 1;
  }
  elsif(substr($path,-4) eq ".fna"){
    return 1;
  }
  elsif(substr($path,-4) eq ".faa"){
    return 1;
  }
  elsif(substr($path,-9) eq ".fasta.gz"){
    return 1;
  }
  elsif(substr($path,-6) eq ".fa.gz"){
    return 1;
  }
  elsif(substr($path,-7) eq ".fna.gz"){
    return 1;
  }
  elsif(substr($path,-7) eq ".faa.gz"){
    return 1;
  }
  else{
    return 0;
  }

}

######################################################################################################
######################################################################################################
######################################################################################################

## FILEREADER
# version 2022.03.04

#file_reader_data
%global_filehandles_read = ();
%global_filehandles_read_isopen = ();
%global_filehandles_read_incomleteline = ();
%global_filehandles_read_lines = ();
%global_filehandles_read_lines_countcomplete = ();
%global_filehandles_read_lines_current = ();
%global_filehandles_read_islastchunk = ();

#file_writer_data
%global_filehandles_write = ();
%global_filehandles_write_isopen = ();

######################################################################################################

sub private_filereader_open(){

   # version 2022.03.04
   my $path = shift(@_);
   my $fileid = shift(@_);
   
   &private_assert_filepath($path);
   
   
   if( exists( $global_filehandles_read_isopen{$fileid})){
     if($global_filehandles_read_isopen{$fileid} ==1){
       close($global_filehandles_read{$fileid});
     }
   }
  
   if(&private_isgz($path)){
     $global_filehandles_read{$fileid} = IO::Zlib->new("$path","rb") ;
   }
   else{
     open( $global_filehandles_read{$fileid}, '<',  $path);
   }
   
   $global_filehandles_read_isopen{$fileid} = 1;
   
   $global_filehandles_read_incomleteline{$fileid} = "";  # first chunk => no incomplete line from previous chunk
   
   &private_filereader_readchunk($fileid);
   
}

######################################################################################################

sub private_filewriter_open(){

   # version 2022.03.04

   my $path = shift(@_);
   my $fileid = shift(@_);
   
   if( exists( $global_filehandles_write_isopen{$fileid})){
     if($global_filehandles_write_isopen{$fileid} ==1){
       close($global_filehandles_write{$fileid});
     }
   }

   if(&private_isgz($path)){
     $global_filehandles_write{$fileid} = IO::Zlib->new("$path","wb") ;
   }
   else{
     open( $global_filehandles_write{$fileid}, '>',  $path);
   }
   
   $global_filehandles_write_isopen{$fileid} = 1;
   

}



######################################################################################################

sub private_filereader_close(){

  # version 2022.03.04
   
  my $fileid = shift(@_);
  
   close($global_filehandles_read{$fileid});
   $global_filehandles_read_isopen{$fileid} = 0;

}

######################################################################################################

sub private_filewriter_close(){

  # version 2022.03.04
    
  my $fileid = shift(@_);
  
   close($global_filehandles_write{$fileid});
   $global_filehandles_write_isopen{$fileid} = 0;

}

######################################################################################################

sub private_filereader_readchunk(){

  # version 2022.03.04
  
   my $fileid = shift(@_);

   my $buffer = "";
   my $lines_count = 0;
   
   ## load chunk
   my $readbytecount = $global_filehandles_read{$fileid}->read($buffer, $setting_filereadbuffersize );
 #  print "$readbytecount\n";
   ### no more data
   if($readbytecount==0 && $global_filehandles_read_incomleteline{$fileid} eq "" ){
      $lines_count = 0; 
      $global_filehandles_read_lines_countcomplete{$fileid} = 0;
      $global_filehandles_read_lines_current{$fileid} = 0; 
      $global_filehandles_read_lines{$fileid} = [];
      return;
    }
   
   ### only last line remaining (could theoretically hapen if the last buffer-read happened to read all the remaining data in the file and the file did not end with a newline)
   if($readbytecount==0 && $global_filehandles_read_incomleteline{$fileid} ne "" ){
      $lines_count = 1; 
      $global_filehandles_read_lines_countcomplete{$fileid} = 1;
      $global_filehandles_read_lines_current{$fileid} = 0; 
      $global_filehandles_read_lines{$fileid} = [];
      push(@{ $global_filehandles_read_lines{$fileid} },$global_filehandles_read_incomleteline{$fileid});
      return;
    }
    

   #if the chunk is smaller than the buffer (did not fill the buffer) => no more data to read
   if($readbytecount<$setting_filereadbuffersize){ 
      $global_filehandles_read_islastchunk{$fileid} = 1; 
   }
   else{
      $global_filehandles_read_islastchunk{$fileid} = 0; 
   }
   
   # if last char is a newline...the last line is complete... otherwise it is not complete untill next chunk has been read
   my $chunk_lastchar = substr($buffer,-1);
   my $lastlinecomplete;
   
   if($chunk_lastchar eq "\n"){
     $lastlinecomplete = 1;
   }
   else{
     $lastlinecomplete = 0;
   }
   
   # split the buffer into lines   
   $global_filehandles_read_lines{$fileid} = [ split(/\n/, $buffer) ];
   $global_filehandles_read_lines_current{$fileid} = 0; 
   $lines_count = @{ $global_filehandles_read_lines{$fileid} }; 

   # add incomplete line to first line .. (first line may be some data or only an empty string if only a newline was in the beginning of the buffer)
   if($global_filehandles_read_incomleteline{$fileid} ne ""){
     $global_filehandles_read_lines{$fileid}->[0] = $global_filehandles_read_incomleteline{$fileid} . $global_filehandles_read_lines{$fileid}->[0];
   #  print "added" . $global_filehandles_read_incomleteline{$fileid}  . "\n";
   }
   # save the last line if it is not complete
   if($lastlinecomplete == 0 && $global_filehandles_read_islastchunk{$fileid}==0){
     $global_filehandles_read_lines_countcomplete{$fileid} =  $lines_count -1;
     $global_filehandles_read_incomleteline{$fileid} = $global_filehandles_read_lines{$fileid}->[$lines_count-1];
    # print "not complete\n";
   }
   else{
     $global_filehandles_read_lines_countcomplete{$fileid} =  $lines_count;
     $global_filehandles_read_incomleteline{$fileid} = "";
   #  print "complete\n";
   }
   
 # print "count:$lines_count\n";
  
}
######################################################################################################

sub private_filereader_getlines(){

  # version 2022.03.04
     
   my $nrlines = shift(@_);
   my $fileid = shift(@_);
   
   my @returnvalue = (); 

   for(my $i =0; $i<$nrlines ;$i++){
      my $line = &private_filereader_Getline($fileid);

      if(!defined($line)) {  # this is the end of file signal
         return @returnvalue;
      }
      push(@returnvalue,$line)
   }
   
   return @returnvalue;
}

######################################################################################################
   
sub private_filereader_Getline(){

  # version 2022.09.16
   
  my $fileid = shift(@_);
   
  if($global_filehandles_read_lines_current{$fileid}>=$global_filehandles_read_lines_countcomplete{$fileid} && $global_filehandles_read_islastchunk{$fileid} ==0){
     &private_filereader_readchunk($fileid);
  }
 
 if($global_filehandles_read_lines_current{$fileid}<$global_filehandles_read_lines_countcomplete{$fileid}){
    my $ret =  $global_filehandles_read_lines{$fileid}->[$global_filehandles_read_lines_current{$fileid}];
    $global_filehandles_read_lines_current{$fileid}= $global_filehandles_read_lines_current{$fileid} +1;

    return $ret . "\n";  #add a newline so that string is never empty (can be interpreted as undef in loops)
 }

 return undef;
 
}
   
######################################################################################################
sub private_filewriter_writeline(){

  # version 2022.03.04
  
  my $line = shift(@_);
  my $fileid = shift(@_);

  print { $global_filehandles_write{$fileid} } $line ."\n" ;

}

######################################################################################################
sub private_filewriter_writeline_array(){

  # version 2022.03.04
  
  my @lines = @{shift(@_)};
  my $fileid = shift(@_);
 
  foreach $line (@lines) {
     &private_filewriter_writeline($line,$fileid);
  }
  #"FQOUT1"
}


######################################################################################################   


sub private_formatstring_readcount(){

  # version 2022.03.04
  
  my $count_seq = shift(@_);
     
  my $format_count_seq = "$count_seq reads";

  if($count_seq>=1000 && $count_seq < 1000000) {
    $tmp = $count_seq/1000;
    $format_count_seq = sprintf("%.2f", $tmp ) . " thousand reads";
  }
  if($count_seq>=1000000 && $count_seq < 1000000000) {
    $tmp = $count_seq/1000000;
    $format_count_seq = sprintf("%.2f", $tmp ) . " million reads";
  }
  if($count_seq>=1000000000 ) {
    $tmp = $count_seq/1000000000;
    $format_count_seq = sprintf("%.2f", $tmp ) . " billion reads";
  }
  
  return $format_count_seq;
} 

######################################################################################################   


sub private_formatstring_readpaircount(){

  # version 2022.03.04
  
  my $count_seq = shift(@_);
     
  my $format_count_seq = "$count_seq reads";

  if($count_seq>=1000 && $count_seq < 1000000) {
    $tmp = $count_seq/1000;
    $format_count_seq = sprintf("%.2f", $tmp ) . " thousand read pairs";
  }
  if($count_seq>=1000000 && $count_seq < 1000000000) {
    $tmp = $count_seq/1000000;
    $format_count_seq = sprintf("%.2f", $tmp ) . " million read pairs";
  }
  if($count_seq>=1000000000 ) {
    $tmp = $count_seq/1000000000;
    $format_count_seq = sprintf("%.2f", $tmp ) . " billion read pairs";
  }
  
  return $format_count_seq;
} 

######################################################################################################   


sub private_formatstring_basecount(){

  # version 2022.03.04
  
  my $count_bases = shift(@_);
     
  my $format_count_bases = "$count_bases bases";

  if($count_bases>=1000 && $count_bases < 1000000) {
    $tmp = $count_bases/1000;
    $format_count_bases = sprintf("%.2f", $tmp ) . " thousand bases";
  }
  if($count_bases>=1000000 && $count_bases < 1000000000) {
    $tmp = $count_bases/1000000;
    $format_count_bases = sprintf("%.2f", $tmp ) . " million bases";
  }
  if($count_bases>=1000000000 ) {
    $tmp = $count_bases/1000000000;
    $format_count_bases = sprintf("%.2f", $tmp ) . " billion bases";
  }
  
  return $format_count_bases;
} 

######################################################################################################   
   
sub private_reverscomplement(){

  # version 2022.03.04
  
  my $seq = shift(@_);

  my $rc = reverse($seq);
  
  $rc =~ tr/acgtACGT/tgcaTGCA/; 

  return $rc;
}

######################################################################################################   

sub private_gettypicalgenomesize(){

  # version 2022.03.04
  
  my $organism = lc($_[0]);
  
   if(exists( $constants_typicalgenomesize{$organism}) ) { 
     return $constants_typicalgenomesize{$organism};
  }
  else{
    return 0;
  }
}

###################################################################################################### 

sub private_pushfilesindir(){

  my $dir = $_[0]; 
  my $filetype = lc($_[1]);

  my $recursive = lc($_[3]); 
  
  $dir= &private_uniformdirpath($dir);
  

  
  opendir(my $dirhandle, $dir) or CORE::say( "Cannot open $dir");
  
  while(my $file = readdir($dirhandle)){  
     

     if($filetype eq "fasta" && &private_isfasta($file)){
        push(@{$_[2]}, $dir . $file);
     } 
     
     if($filetype eq "fastq" && &private_isfastq($file)){
        push(@{$_[2]}, $dir . $file);
     } 
   
     if($recursive  eq "recursive" && -d $dir . $file  &&  $file ne "." &&  $file ne ".."){
        my $subdir= $dir . $file;

        $subdir= &private_uniformdirpath($subdir);
        
        &private_pushfilesindir($subdir,$filetype,$_[2],$recursive);
     } 
  }
  closedir( $dirhandle);
  
   
} 

###################################################################################################### 


sub private_mergeR1R2(){


  my $count = @{$_[0]};
 
  my @R1s = ();
  my @R2s = ();
  my @singletons = ();
  
  ## split R1, R2 and singletons
  foreach $file  (@{$_[0]}) {
  
  
    if(  &private_hasR1Extension($file)==1 ){
       push(@R1s,$file);
    }
    elsif(  &private_hasR2Extension($file)==1 ){
       push(@R2s,$file);
    }
    else{
      push(@singletons,$file);
    }
   
  }

   #empty base array
   @{$_[0]} = ();

   # re-populate with merge
   %merged_R2s = ();
   
   ## go through R1s
   for(my $i=0; $i< @R1s;$i++){
   
     $R1read = $R1s[$i];
     $expected_R2 = &private_toggleR1R2Extension($R1read);
     $R2read = "";
     
     for(my $j=0; $j< @R2s;$j++){
       $tmp_R2 = $R2s[$j];
       if($tmp_R2 eq $expected_R2){
         $R2read = $tmp_R2;
         $merged_R2s{$tmp_R2} = 1; 
       }
     }
     
     if($R2read ne ""){
        push(@{$_[0]},"$R1read\t$R2read");
     }
     else{
        push(@{$_[0]},"$R1read");
     }
       
   }
   
   ## go through R2s
   for(my $i=0; $i< @R2s;$i++){
   
     $R2read = $R2s[$i];
     $added = 0;
     
     if( ! exists( $merged_R2s{$R2read}) ) { 
       push(@{$_[0]},"$R2read"); 
     }
       
   }

   ## go through singletons
   for(my $i=0; $i< @singletons;$i++){
   
     $read = $singletons[$i];
     push(@{$_[0]},"$read"); 
      
   }

}

###################################################################################################### 

sub private_uniformdirpath(){

  my $path = shift(@_);

  if(substr($path,-1) ne "/"){
      $path = $path . "/";
    }
    
  return $path;
    
}

###################################################################################################### 

sub private_getlastpartofpath(){

  my $path = shift(@_);

  if(index($path, "/")){
    $path = substr($path, index($path, "/")+1);
  }
    
  return $path;
    
}

###################################################################################################### 


sub private_hasR1Extension(){

    my $path = lc(shift(@_));
    
    if(substr($path,-1) eq "/" || substr($path,-1) eq "\\"){
      $path = substr($path,0,-1); 
    }
    
    foreach $R1_ext (@constants_R1_extensions) {
      
         my $len = length($R1_ext);
         
         if(substr($path,-1*$len) eq  $R1_ext){
         
            return 1;
         
         }
      
    }
    
    return 0;


}

###################################################################################################### 


sub private_hasR2Extension(){

    my $path = lc(shift(@_));
    
    if(substr($path,-1) eq "/" || substr($path,-1) eq "\\"){
      $path = substr($path,0,-1); 
    }
    
    foreach $R2_ext (@constants_R2_extensions) {
      
         my $len = length($R2_ext);
         
         if(substr($path,-1*$len) eq  $R2_ext){
         
            return 1;
         
         }
      
    }
    
    return 0;


}

###################################################################################################### 


sub private_toggleR1R2Extension(){

    my $path_casesens = shift(@_);
    my $path = lc($path_casesens);


    if(substr($path,-1) eq "/" || substr($path,-1) eq "\\"){
      $path = substr($path,0,-1); 
    }
    
    for($i=0;$i<@constants_R1_extensions;$i++) {
      
       my $R1_ext = $constants_R1_extensions[$i];
       my $R2_ext = $constants_R2_extensions[$i];
       
       $value = substr($path,-1*length($R1_ext));

       #### R1 pattern  
       if(substr($path,-1*length($R1_ext)) eq  $R1_ext){
         
          
         my $len = length($R1_ext);
         
         my $caseSens_R1_ext = substr($path_casesens,-1*$len); 
         my $caseSens_R2_ext = "";
         
         
         for($j=0;$j<length($caseSens_R1_ext);$j++){
            
             if(substr($caseSens_R1_ext,$j,1) eq lc(substr($caseSens_R1_ext,$j,1)) ){
             
                $caseSens_R2_ext = $caseSens_R2_ext . lc(substr($R2_ext,$j,1)); 
             }
             else{
                 $caseSens_R2_ext = $caseSens_R2_ext . uc(substr($R2_ext,$j,1)); 
             }
         
         }
         
         my $R2path = substr($path_casesens,0,-1*$len) . $caseSens_R2_ext;
         return $R2path;
         
        }  ## end R1 pattern

       #### R2 pattern  
       if(substr($path,-1*length($R2_ext)) eq  $R2_ext){
         
         $len = length($R2_ext);
         
         $caseSens_R2_ext = substr($path_casesens,-1*$len); 
         $caseSens_R1_ext = "";
         
         for($j=0;$j<length($caseSens_R2_ext);$j++){
            
             if(substr($caseSens_R2_ext,$j,1) eq lc(substr($caseSens_R2_ext,$j,1)) ){
             
                $caseSens_R1_ext = $caseSens_R1_ext . lc(substr($R1_ext,$j,1)); 
             }
             else{
                 $caseSens_R1_ext = $caseSens_R1_ext . uc(substr($R1_ext,$j,1)); 
             }
         
         }
         
         my $R1path = substr($path_casesens,0,-1*$len) . $caseSens_R1_ext;
         return $R1path;
         
        }  ## end R1 pattern
      
 
      
    }
    
       #no R1 or R2 pattern match
      return "";



}

######################################################################################################   

sub addtagtopath(){

  # version 2022.03.04
  
  my $path = shift(@_);
  my $tag = shift(@_);

  my $newpath=$path . $tag;
  my $extension = "";

  $extension = ".fastq.gz";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  }
  
  $extension = ".fastq";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  }


  $extension = ".fq.gz";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  } 
  $extension = ".fq";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  }
  
  $extension = ".fasta.gz";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  }
  
  $extension = ".fasta";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
    
     return $newpath;
  }
  $extension = ".fa.gz";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  } 
  $extension = ".fa";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  }
  
  $extension = ".fna.gz";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  } 
  $extension = ".fna";
  if(substr(lc($path),-1*length($extension)) eq $extension   ){
     $newpath = substr($path,0,-1*length($extension)) . $tag . substr($path,-1*length($extension));
     return $newpath;
  }
  
  return $newpath;
}

