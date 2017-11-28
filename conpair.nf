#!/usr/bin/env nextflow

// requires (in path):
// conpair (https://github.com/nygenome/Conpair)
// python 2.7 : www.python.org
// numpy 1.7.0 or higher : www.numpy.org
// scipy 0.14.0 or higher : www.scipy.org
// GATK 2.3 or higher : www.broadinstitute.org/gatk/download
// java : http://java.com

params.help = null
params.ref = null
params.markers = null
params.tumor_bam_folder = null
params.normal_bam_folder = null
params.bam_folder = null
params.conpair_dir = null
params.gatk_jar = null
params.concordance_out = 'concordance_summary.txt'
params.contamination_out = 'contamination_summary.txt'
if (params.markers){
  markers_tag = "--markers"
} else {
  markers_tag = ""
}

if (params.help) {
    log.info ''
    log.info '------------------------------------------------------------'
    log.info ' NEXTFLOW for Conpair (https://github.com/nygenome/Conpair) '
    log.info '------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run conpair.nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --ref ref.fasta' 
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '   When using Tumor/Normal pairs:'     
    log.info '    --tumor_bam_folder   FOLDER                  Folder containing tumor BAM files.'
    log.info '    --normal_bam_folder  FOLDER                  Folder containing matched normal BAM files.'
    log.info '   In other cases:'     
    log.info '    --bam_folder         FOLDER                  Folder containing BAM files.'
    log.info '   In all cases:'                
    log.info '    --ref                FILE (with index)       Reference fasta file indexed by bwa.'
    log.info '    --conpair_dir        DIR                     Conpair directory.'
    log.info '    --gatk_jar           FILE                    GATK.jar explicit path.'   
    log.info 'Optional arguments:'
    log.info '    --suffix_tumor       STRING                  Suffix identifying tumor bam (default: "_T").'
    log.info '    --suffix_normal      STRING                  Suffix identifying normal bam (default: "_N").'
    log.info '    --concordance_out    FILE                    Output file name for concorfance summary (default: concordance_summary.txt).'
    log.info '    --contamination_out  FILE                    Output file name for concorfance summary (default: contamination_summary.txt).'
    log.info ''
    exit 1
}

assert (params.ref != true) && (params.ref != null) : "please specify --ref option (--ref reference.fasta(.gz))"

if(params.bam_folder != null) {
    assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"
} else {
    assert (params.normal_bam_folder != true) && (params.normal_bam_folder != null) : "please specify --normal_bam_folder option (--normal_bam_folder bamfolder)"
    assert (params.tumor_bam_folder != true) && (params.tumor_bam_folder != null) : "please specify --tumor_bam_folder option (--tumor_bam_folder bamfolder)"
}

assert (params.conpair_dir != true) && (params.conpair_dir != null) : "please specify --conpair_dir option (--conpair_dir /path/to/conpair/)"
assert (params.gatk_jar != true) && (params.gatk_jar != null) : "please specify --gatk_jar option (--gatk_jar /path/to/gatk.jar)"

fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )
fasta_ref_gzi = file( params.ref+'.gzi' )
fasta_ref_dict = file( params.ref.replace(".fasta",".dict").replace(".fa",".dict") )

params.suffix_tumor = "_T"
params.suffix_normal = "_N"
params.out_folder = "."

try { assert fasta_ref.exists() : "\n WARNING : fasta reference not located in execution directory. Make sure reference index is in the same folder as fasta reference" } catch (AssertionError e) { println e.getMessage() }
if (fasta_ref.exists()) {assert fasta_ref_fai.exists() : "input fasta reference does not seem to have a .fai index (use samtools faidx)"}
if (fasta_ref.exists()) {assert fasta_ref_dict.exists() : "input fasta reference does not seem to have a .dict index (use GATK to index)"}
if (fasta_ref.exists() && params.ref.tokenize('.')[-1] == 'gz') {assert fasta_ref_gzi.exists() : "input gz fasta reference does not seem to have a .gzi index (use samtools faidx)"}

if(params.bam_folder) {

    try { assert file(params.bam_folder).exists() : "\n WARNING : input BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
    assert file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0 : "BAM folder contains no BAM"

    // recovering of bam files
    bams = Channel.fromPath( params.bam_folder+'/*.bam' )
              .ifEmpty { error "Cannot find any bam file in: ${params.bam_folder}" }
              .map {  path -> [ path.name.replace(".bam",""), path ] }

    // recovering of bai files
    bais = Channel.fromPath( params.bam_folder+'/*.bam.bai' )
              .ifEmpty { error "Cannot find any bai file in: ${params.bam_folder}" }
              .map {  path -> [ path.name.replace(".bam.bai",""), path ] }

    // building bam-bai pairs
    bam_bai = bams
              .phase(bais)
              .map { bam, bai -> [ bam[1], bai[1] ] }
              
	process GATK_pileup {

	   tag { bam_tag }

	   input:
	   file bam_bai
	   file fasta_ref
	   file fasta_ref_fai
	   file fasta_ref_gzi
	   file fasta_ref_dict 

	   output:
	   set val(bam_tag), file("${bam_tag}.pileup") into pileup

	   shell:
	   bam_tag = bam_bai[0].baseName
	   '''
	   set +u
	   export CONPAIR_DIR=!{params.conpair_dir}
	   export GATK_JAR=!{params.gatk_jar}
	   export PYTHONPATH=${PYTHONPATH}:!{params.conpair_dir}

	   !{params.conpair_dir}/scripts/run_gatk_pileup_for_sample.py -B !{bam_tag}.bam -O !{bam_tag}.pileup --reference !{fasta_ref} --temp_dir_java tmp --remove_chr_prefix
	   '''
	}
             
	process contamination {
   
	   tag { bam_tag }

	   publishDir params.out_folder+'/all_contamination/', mode: 'move'

	   input:
	   set val(bam_tag), file(pileup_file) from pileup
   
	   output:
	   file("${bam_tag}_contamination.txt")
	   stdout contamination_summary

	   shell:
	   '''
	   set +u
	   export CONPAIR_DIR=!{params.conpair_dir}
	   export PYTHONPATH=${PYTHONPATH}:!{params.conpair_dir}/modules/
 
	   !{params.conpair_dir}/scripts/estimate_tumor_normal_contamination.py -T !{pileup_file} -N !{pileup_file} --outfile !{bam_tag}_contamination.txt
	   echo -ne !{bam_tag}"\t"
	   grep Normal !{bam_tag}_contamination.txt | sed 's/.*: //'
	   '''
	}

	contamination_summary.collectFile(name: params.contamination_out, storeDir: params.out_folder, seed: 'SAMPLE\tCONTAMINATION\n')             
                
}  else {

	try { assert file(params.tumor_bam_folder).exists() : "\n WARNING : input tumor BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
	assert file(params.tumor_bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0 : "tumor BAM folder contains no BAM"
	try { assert file(params.normal_bam_folder).exists() : "\n WARNING : input normal BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
	assert file(params.normal_bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0 : "normal BAM folder contains no BAM"

	// FOR TUMOR 
	// recovering of bam files
	tumor_bams = Channel.fromPath( params.tumor_bam_folder+'/*'+params.suffix_tumor+'.bam' )
		    .ifEmpty { error "Cannot find any bam file in: ${params.tumor_bam_folder}" }
		    .map {  path -> [ path.name.replace("${params.suffix_tumor}.bam",""), path ] }

	// recovering of bai files
	tumor_bais = Channel.fromPath( params.tumor_bam_folder+'/*'+params.suffix_tumor+'.bam.bai' )
		    .ifEmpty { error "Cannot find any bai file in: ${params.tumor_bam_folder}" }
		    .map {  path -> [ path.name.replace("${params.suffix_tumor}.bam.bai",""), path ] }

	// building bam-bai pairs
	tumor_bam_bai = tumor_bams
		    .phase(tumor_bais)
		    .map { tumor_bam, tumor_bai -> [ tumor_bam[0], tumor_bam[1], tumor_bai[1] ] }

	// FOR NORMAL 
	// recovering of bam files
	normal_bams = Channel.fromPath( params.normal_bam_folder+'/*'+params.suffix_normal+'.bam' )
		    .ifEmpty { error "Cannot find any bam file in: ${params.normal_bam_folder}" }
		    .map {  path -> [ path.name.replace("${params.suffix_normal}.bam",""), path ] }

	// recovering of bai files
	normal_bais = Channel.fromPath( params.normal_bam_folder+'/*'+params.suffix_normal+'.bam.bai' )
		    .ifEmpty { error "Cannot find any bai file in: ${params.normal_bam_folder}" }
		    .map {  path -> [ path.name.replace("${params.suffix_normal}.bam.bai",""), path ] }

	// building bam-bai pairs
	normal_bam_bai = normal_bams
		    .phase(normal_bais)
		    .map { normal_bam, normal_bai -> [ normal_bam[0], normal_bam[1], normal_bai[1] ] }

	// building 4-uplets corresponding to {tumor_bam, tumor_bai, normal_bam, normal_bai}
	tn_bambai = tumor_bam_bai
		.phase(normal_bam_bai)
		.map {tumor_bb, normal_bb -> [ tumor_bb[1], tumor_bb[2], normal_bb[1], normal_bb[2] ] }    
	// here each element X of tn_bambai channel is a 4-uplet. X[0] is the tumor bam, X[1] the tumor bai, X[2] the normal bam and X[3] the normal bai.

	process GATK_pileup {

	   tag { tumor_normal_tag }

	   input:
	   file tn from tn_bambai
	   file fasta_ref
	   file fasta_ref_fai
	   file fasta_ref_gzi
	   file fasta_ref_dict 

	   output:
	   set val(tumor_normal_tag), file("${tumor_normal_tag}${params.suffix_normal}.pileup"), file("${tumor_normal_tag}${params.suffix_tumor}.pileup") into pileup_pairs1, pileup_pairs2

	   shell:
	   tumor_normal_tag = tn[0].baseName.replace(params.suffix_tumor,"")
	   '''
	   set +u
	   export CONPAIR_DIR=!{params.conpair_dir}
	   export GATK_JAR=!{params.gatk_jar}
	   export PYTHONPATH=${PYTHONPATH}:!{params.conpair_dir}/modules/

	   !{params.conpair_dir}/scripts/run_gatk_pileup_for_sample.py -B !{tumor_normal_tag}!{params.suffix_tumor}.bam -O !{tumor_normal_tag}!{params.suffix_tumor}.pileup --reference !{fasta_ref} !{markers_tag} !{params.markers} --temp_dir_java tmp --remove_chr_prefix
	   !{params.conpair_dir}/scripts/run_gatk_pileup_for_sample.py -B !{tumor_normal_tag}!{params.suffix_normal}.bam -O !{tumor_normal_tag}!{params.suffix_normal}.pileup --reference !{fasta_ref} !{markers_tag} !{params.markers} --temp_dir_java tmp --remove_chr_prefix
	   '''
	}

	process concordance {
   
	   tag { tumor_normal_tag }

	   publishDir params.out_folder+'/all_concordance/', mode: 'move'

	   input:
	   set val(tumor_normal_tag), file(normal), file(tumor) from pileup_pairs1
   
	   output:
	   file("${tumor_normal_tag}_concordance.txt")
	   stdout concordance_summary

	   shell:
	   '''
	   set +u
	   export CONPAIR_DIR=!{params.conpair_dir}
	   export PYTHONPATH=${PYTHONPATH}:!{params.conpair_dir}/modules/
	 
	   !{params.conpair_dir}/scripts/verify_concordance.py -T !{tumor} -N !{normal} --outfile !{tumor_normal_tag}_concordance.txt
	   echo -ne !{tumor_normal_tag}"\t"
	   grep Concordance !{tumor_normal_tag}_concordance.txt | sed 's/.*: //' 
	   '''
	}

	concordance_summary.collectFile(name: params.concordance_out, storeDir: params.out_folder, seed: 'SAMPLE\tCONCORDANCE\n')


	process contamination {
   
	   tag { tumor_normal_tag }

	   publishDir params.out_folder+'/all_contamination/', mode: 'move'

	   input:
	   set val(tumor_normal_tag), file(normal), file(tumor) from pileup_pairs2
   
	   output:
	   file("${tumor_normal_tag}_contamination.txt")
	   stdout contamination_summary

	   shell:
	   '''
	   set +u
	   export CONPAIR_DIR=!{params.conpair_dir}
	   export PYTHONPATH=${PYTHONPATH}:!{params.conpair_dir}/modules/
 
	   !{params.conpair_dir}/scripts/estimate_tumor_normal_contamination.py -T !{tumor} -N !{normal} --outfile !{tumor_normal_tag}_contamination.txt
	   echo -ne !{tumor_normal_tag}"\t"
	   grep Normal !{tumor_normal_tag}_contamination.txt | sed 's/.*: //'  | tr '\n' '\t'
	   grep Tumor !{tumor_normal_tag}_contamination.txt | sed 's/.*: //'
	   '''
	}

	contamination_summary.collectFile(name: params.contamination_out, storeDir: params.out_folder, seed: 'SAMPLE\tNORMAL\tTUMOR\n')
}
