//! Command-line tools for processing HiFiRe3 data.

// imports
use std::env;
use std::error::Error;

// modules
mod formats;
mod inserts;
mod junctions;
mod snvs;
mod sites;
mod tools;
mod benchmark;

// constants
const TOOLS_NAME: &str = "hf3_tools";

// load and process data
fn main() -> Result<(), Box<dyn Error>> {

    // read command line arguments
    let mut args: Vec<String> = env::args().collect();
    args = args[1..].to_vec(); // drop executable name 
    if args.len() == 0 { // check for something to do, i.e., a tool to run
        eprintln!("{}: missing tool or command", TOOLS_NAME);
        Err(format!("usage: {} <tool> [additional arguments]", TOOLS_NAME))?
    }
    let tool = args[0].clone(); // drop tool name
    // args = args[1..].to_vec(); // uncomment if passing additional arguments to tools

    // dispatch to tool or command
    match tool.as_str() {

        /*--------------------------------------------------------------
        basecall ONT
        ------------------------------------------------------------- */
        // RE site-aware trimming for ONT reads
        "trim_ont" => tools::trim_ont::stream(),

        // compression and reformatting of legacy ONT POD5 files
        "reformat_ont" => tools::reformat_ont::stream(),

        /*--------------------------------------------------------------
        basecall PacBio
        ------------------------------------------------------------- */
        // merge fwd and rev PacBio reads into pbFree "unleaded" reads
        "basecall_pacbio" => tools::basecall_pacbio::main(),

        /*--------------------------------------------------------------
        analyze fragments
        ------------------------------------------------------------- */

        // perform pre- and post-alignment processing steps on Ultima aligned reads
        "post_process_ultima" => tools::post_process_ultima::main(),

        // sort (per read) and examine read alignments in the output stream from minimap2
        "analyze_alignments" => tools::analyze_alignments::stream(),

        // assess inserts based on RE site matching and insert sizing/stem lengths
        "analyze_inserts" => tools::analyze_inserts::stream(),

        /*--------------------------------------------------------------
        analyze SVs
        ------------------------------------------------------------- */
        // prepare reads for SV analysis by splitting BAM to chrom-level files
        "split_bam_by_chrom_sv" => tools::split_bam_by_chrom_sv::main(),

        // index fragments to collect read paths and junctions
        "analyze_svs" => tools::analyze_svs::main(),

        /*--------------------------------------------------------------
        analyze SNVs
        ------------------------------------------------------------- */
        // prepare reads for SNV analysis by splitting BAM to chrom-level files
        "split_bam_by_chrom_snv" => tools::split_bam_by_chrom_snv::main(),

        // analyze RE fragment haplotypes for clonal and subclonal variants
        "analyze_snvs" => tools::analyze_snvs::main(),

        /*--------------------------------------------------------------
        merge SVs
        ------------------------------------------------------------- */
        // fuzzy merge final junctions across related sample of different library types
        "merge_svs" => tools::merge_svs::main(),

        /*--------------------------------------------------------------
        benchmark SVs
        ------------------------------------------------------------- */
        // digest allelic consensus sequences to determine RE fragment sizes
        "digest_sv_alleles" => tools::digest_sv_alleles::main(),
        
        // run SV benchmark comparison
        "compare_sv_calls" => tools::compare_sv_calls::main(),

        /*--------------------------------------------------------------
        unrecognized pipeline action tool
        ------------------------------------------------------------- */
        _ => Err(format!("{}: unknown tool or command: {}", TOOLS_NAME, tool))?
    }
}
