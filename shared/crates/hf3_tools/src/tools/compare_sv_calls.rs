//! Compare SVs as called by:
//!     - HiFiRe3 (FinalJunctions) as the test set
//!     - Sniffles2 and CuteSV (VCF) at two truth sets
//! Options from `library-properties`, `restriction-enzyme`, and `targets` are
//! applied to the HiFiRe3 test set to establish operating parameters, i.e., to
//! restrict attention to SV junctions from the truth sets that could have been 
//! called by the test set. The truth sets are used to establish the external
//! list of expected SVs and their RE fragment behavior. Thus, `compare_sv_calls` 
//! supports comparisons of:
//!     - HiFiRe3 and Sniffles2/CuteSV run on the same data, as well as
//!     - HiFiRe3 run on HiFiRe3 data and Sniffles2/CuteSV run on external truth data.
//! to benchmark both the HiFiRe3 pipeline alone and the library+pipeline 
//! combination. Along the way, Sniffles2 and CuteSV are also implicitly compared,
//! although exact comparisons are confounded by different encoding between the
//! SV callers while neither reports a detailed assessment of junction properties.

// dependencies
use std::error::Error;
use std::fs;
use std::path::PathBuf;
use rustc_hash::FxHashMap;
use rust_htslib::bcf;
use rust_htslib::bcf::Read as VcfRead;
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use mdi::{InputCsv, OutputCsv};
use genomex::genome::{Chroms, TargetRegions, Exclusions};
use crate::junctions::FinalJunction;
use crate::benchmark::*;

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "compare_sv_calls";
pub_key_constants!(
    // from environment variables
    N_CPU
    SEQUENCING_PLATFORM
    GENOME_FASTA
    IS_COMPOSITE_GENOME
    ENZYME_NAME
    MIN_N_OBSERVED
    MIN_SV_MAPQ
    MIN_SV_LEN
    MATCH_DISTANCE
    TRUTHSET_DIR
    TRUTHSET_IS_COMPOSITE_GENOME
    PRIMARY_GENOME_EXCLUSIONS
    SV_FINAL_JUNCTIONS_FILE_1
    BENCHMARK_ADMIXED_CALLS
    BENCHMARK_AGGREGATED_CALLS
    // counter keys
    N_INPUT_JXNS_HF3
    N_INPUT_JXNS_SNF
    N_INPUT_JXNS_CSV
    N_OUTPUT_JXNS_HF3
    N_OUTPUT_JXNS_SNF
    N_OUTPUT_JXNS_CSV
);

/// Tool for processing VCF SV calls.
pub struct ComparisonTool {
    min_n_observed: u16,
    _min_sv_mapq: u8,
    min_sv_len: u32,
    match_distance: u32,
    pub chroms:  Chroms,
    pub targets: TargetRegions,
    pub is_composite_genome: bool,
    truthset_is_composite_genome: bool,
    pub primary_genome_exclusions: Exclusions,
    pub re: ReSite,
    jxn_calls: Vec<JunctionCall>,
}

// function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[
        SEQUENCING_PLATFORM, GENOME_FASTA, ENZYME_NAME, 
        SV_FINAL_JUNCTIONS_FILE_1, TRUTHSET_DIR, PRIMARY_GENOME_EXCLUSIONS,
        BENCHMARK_ADMIXED_CALLS, BENCHMARK_AGGREGATED_CALLS
    ]);
    cfg.set_usize_env( &[N_CPU, MIN_N_OBSERVED, MIN_SV_MAPQ, MIN_SV_LEN, MATCH_DISTANCE]);
    cfg.set_bool_env( &[IS_COMPOSITE_GENOME, TRUTHSET_IS_COMPOSITE_GENOME]);

    // initialize counters
    let ctrs = Counters::new(TOOL, &[
        (N_INPUT_JXNS_HF3,  "number of input lines in HiFiRe3 final junctions"),
        (N_INPUT_JXNS_CSV,  "number of input lines in CuteSV truth-set VCF"),
        (N_INPUT_JXNS_SNF,  "number of input lines in Sniffles truth-set VCF"),
        (N_OUTPUT_JXNS_HF3, "number of HF3 output junctions"),
        (N_OUTPUT_JXNS_CSV, "number of CSV output junctions"),
        (N_OUTPUT_JXNS_SNF, "number of SNF output junctions"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // collect the working chromosomes and targets
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::from_env(&mut w, false);

    // create the digestion tool
    let mut tool = ComparisonTool {
        min_n_observed: *w.cfg.get_usize(MIN_N_OBSERVED) as u16,
        _min_sv_mapq:   *w.cfg.get_usize(MIN_SV_MAPQ) as u8,
        min_sv_len:     *w.cfg.get_usize(MIN_SV_LEN) as u32,
        match_distance: *w.cfg.get_usize(MATCH_DISTANCE) as u32,
        chroms:  chroms,
        targets: targets,
        is_composite_genome: *w.cfg.get_bool(IS_COMPOSITE_GENOME),
        truthset_is_composite_genome: *w.cfg.get_bool(TRUTHSET_IS_COMPOSITE_GENOME),
        primary_genome_exclusions: Exclusions::new(w.cfg.get_string(PRIMARY_GENOME_EXCLUSIONS), false),
        re: ReSite::from_re_name(w.cfg.get_string(ENZYME_NAME), None),
        jxn_calls: Vec::new(),
    };

    // load junctions into a single table identified by source
    load_hf3_jxns(&mut w, &mut tool)?;
    load_vcf_jxns(&mut w, &mut tool, "sniffles", N_INPUT_JXNS_SNF, N_OUTPUT_JXNS_SNF)?;
    load_vcf_jxns(&mut w, &mut tool, "cutesv", N_INPUT_JXNS_CSV, N_OUTPUT_JXNS_CSV)?;
    let jxn_id = fuzzy_match_junctions(&mut w, &mut tool);
    rescue_cutesv_dels(&mut w, &mut tool, jxn_id);
    tool.jxn_calls.retain(|jxn| jxn.jxn_id > 0); // drop rejected junctions after group assessment
    write_admixed_calls(&mut w, &mut tool);
    write_aggregated_calls(&mut w, &mut tool);

    // print counts
    w.ctrs.print_grouped(&[
        &[N_INPUT_JXNS_HF3,  N_INPUT_JXNS_SNF,  N_INPUT_JXNS_CSV],
        &[N_OUTPUT_JXNS_HF3, N_OUTPUT_JXNS_SNF, N_OUTPUT_JXNS_CSV],
    ]);
    Ok(())
}

// -----------------------------------------------------------------------------
// load junctions from all data sets
// -----------------------------------------------------------------------------
/// Load HiFiRe3 Final Junctions.
fn load_hf3_jxns(
    w: &mut Workflow, 
    tool: &mut ComparisonTool
) -> Result<(), Box<dyn Error>>{
    w.log.print("loading HiFiRe3 junctions");

    // open the HiFiRe3 final junctions file
    let filepath = w.cfg.get_string(SV_FINAL_JUNCTIONS_FILE_1);
    let mut csv = InputCsv::open(filepath);

    // process junction lines
    for result in csv.deserialize::<FinalJunction>() {
        let jxn = result?;
        w.ctrs.increment(N_INPUT_JXNS_HF3);

        // only continue with primary species canonical chromosomes
        let Some(mut chrom_1) = hf3_chrom_index_to_name(
            jxn.chrom_index1_1, tool
        ) else { continue };
        let Some(mut chrom_2) = hf3_chrom_index_to_name(
            jxn.chrom_index1_2, tool
        ) else { continue };

        // adjust output to be consistent with Sniffles VCF encoding
        // includes strict lexical chromosome ordering
        let mut pos1_1 = jxn.ref_pos1_1;
        let mut pos1_2 = jxn.ref_pos1_2; 
        let mut strands = jxn.get_strand_encoding();
        let pos1_2_csv = jxn.adjust_pos1_2_to_cutesv(pos1_2);
        let sv_size = jxn.adjust_to_vcf(
            &mut chrom_1,
            &mut pos1_1,
            &mut chrom_2,
            &mut pos1_2,
            &mut strands,
        );

        // filter to retain only on-target junctions
        // redundant with SV calling, but ensure consistency with load_snf_jxns()
        if !sv_is_on_target(
            &chrom_1, pos1_1, &chrom_1, pos1_1, &tool
        ) { continue  }

        // reject SVs in excluded regions
        if sv_is_excluded(
            &chrom_1, pos1_1, &chrom_1, pos1_1, &tool
        ) { continue  }

        // commit the junction call
        tool.jxn_calls.push(JunctionCall { 
            jxn_id: 0,

            sv_type: SvType::from_final_junction(&jxn), 
            strands: strands,

            source:     SvSource::HiFiRe3,

            chrom_1, 
            pos1_1, 
            chrom_2, 
            pos1_2, 
            pos1_2_csv,

            sv_size:     sv_size, 
            frag_size:   jxn.get_median_insert_size(),
            stem_length: jxn.min_stem_length,

            offset:  jxn.offset,
            jxn_seq: jxn.jxn_seq.clone(),
            
            n_observed: jxn.n_instances_dedup, //n_reads_dedup,
            mapq:       jxn.max_min_mapq,
            divergence: jxn.min_max_divergence,

            jxn_fail:     jxn.jxn_failure_flag,
            aln_fail:     jxn.aln_failure_flag,
            all_chimeric: (jxn.any_was_not_chimeric == 0) as u8,

            mosaic: false,
            vaf: 0.0,

            qnames: get_sorted_qnames_str(&jxn.qnames),
        });
        w.ctrs.increment(N_OUTPUT_JXNS_HF3);
    }
    Ok(())
}

/// Load truth-set VCF as generated by benchmark::align+digest.
fn load_vcf_jxns(
    w: &mut Workflow, 
    tool: &mut ComparisonTool, 
    sv_caller: &str,
    n_input_jxns: &str,
    n_output_jxns: &str,
) -> Result<(), Box<dyn Error>>{
    w.log.print(&format!("loading {} truthset junctions", sv_caller));
    let sv_source = match sv_caller {
        "sniffles" => SvSource::Sniffles,
        _ =>          SvSource::CuteSV,
    };

    // find and open the truth set VCF file
    let dir = w.cfg.get_string(TRUTHSET_DIR);
    let mut filepath = None::<PathBuf>;
    for entry in fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_file() && path
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with(&format!("{}.frag_sizes.vcf.gz", sv_caller))) {
            filepath = Some(path);
            break;
        }
    }
    let filepath = filepath.expect(&format!("Missing {} truth set VCF.", sv_caller));
    let mut vcf = bcf::Reader::from_path(filepath)?;

    // process junction lines, as returned by `digest_sv_alleles`
    for result in vcf.records() {
        let jxn = result?;
        w.ctrs.increment(n_input_jxns);

        // unpack the initial encoding
        let alt = get_vcf_alt(&jxn)?;
        let (chrom_1, mut pos1_1) = get_vcf_left_node(&jxn)?;
        let (sv_type, sv_type_vcf) = SvType::from_vcf(&jxn, &chrom_1, alt);
        let (chrom_2, mut pos1_2) = get_vcf_right_node(&jxn, &chrom_1, alt, &sv_type_vcf)?;

        // only continue with primary species canonical chromosomes
        let Some(mut chrom_1) = simplify_chrom_name(
            &chrom_1, tool.truthset_is_composite_genome
        ) else { continue };
        let Some(mut chrom_2) = simplify_chrom_name(
            &chrom_2, tool.truthset_is_composite_genome
        ) else { continue };

        // orient translocation junctions using lexical-sorted chroms
        let mut strands = get_vcf_strand_str(
            &jxn, 
            &sv_source,
            &sv_type_vcf,
            alt,
        )?;
        if chrom_1 > chrom_2 {
            (chrom_1, pos1_1, chrom_2, pos1_2) = 
            (chrom_2, pos1_2, chrom_1, pos1_1);
            strands = BndType::invert_jxn_strands(strands);
        }

        // filter to retain only on-target junctions
        if !sv_is_on_target(
            &chrom_1, pos1_1, &chrom_1, pos1_1, &tool
        ) { continue  }

        // reject SVs in excluded regions
        if sv_is_excluded(
            &chrom_1, pos1_1, &chrom_1, pos1_1, &tool
        ) { continue  }

        // collect additional values
        let pos1_2_csv = if sv_source == SvSource::CuteSV { pos1_2 } else { 0 };
        let sv_size = get_vcf_sv_size(&chrom_1, pos1_1, &chrom_2, pos1_2);
        let qnames = get_vcf_qnames(&jxn)?.expect("Missing INFO/RNAMES");

        // commit junction call
        tool.jxn_calls.push(JunctionCall { 
            jxn_id: 0,

            sv_type, 
            strands: strands,

            source: sv_source.clone(),

            chrom_1, 
            pos1_1, 
            chrom_2, 
            pos1_2, 
            pos1_2_csv,

            sv_size, 
            frag_size:   get_frag_size(&jxn, &tool)?,
            stem_length: get_stem_length(&jxn, &tool)?,

            offset: -999,
            jxn_seq: "?".to_string(), // microhomology not available
            
            n_observed: get_vcf_n_observed(&jxn, &sv_source)?,
            mapq: jxn.qual() as u8,
            divergence: 0.0,

            jxn_fail: 0,
            aln_fail: 0,
            all_chimeric: 0,

            mosaic: get_vcf_mosaic(&jxn, &sv_source)?,
            vaf:    get_vcf_vaf(&jxn, &sv_source)?,

            qnames: get_sorted_qnames_bytes(qnames),
        });
        w.ctrs.increment(n_output_jxns);
    }
    Ok(())
}

// -----------------------------------------------------------------------------
// write a table of admixed SV/junction calls
// -----------------------------------------------------------------------------
/// Group, i.e., fuzzy match, junctions across algorithms by setting their `jxn_id`.
pub fn fuzzy_match_junctions(
    w: &mut Workflow, 
    tool: &mut ComparisonTool
) -> u32 {
    w.log.print("grouping junctions");
    let mut jxn_id = 0_u32;
    let mut jxns: Vec<_> = tool.jxn_calls.iter_mut().map(|x| x).collect();

    // fuzzy group by sv_type, strands, chrom_1, pos_1
    jxns.sort_unstable_by_key(|jxn| (
        jxn.sv_type.clone(), jxn.strands.clone(), jxn.chrom_1.clone(), jxn.pos1_1
    ));
    for jxns in jxns.chunk_by_mut(|a, b|{
        a.sv_type == b.sv_type && a.strands == b.strands && a.chrom_1 == b.chrom_1 && 
        b.pos1_1 - a.pos1_1 <= tool.match_distance
    }){

        // fuzzy group by chrom_2, pos_2
        jxns.sort_unstable_by_key(|jxn| (
            jxn.chrom_2.clone(), jxn.pos1_2
        ));
        for jxns in jxns.chunk_by_mut(|a, b|{
            a.chrom_2 == b.chrom_2 && 
            b.pos1_2 - a.pos1_2 <= tool.match_distance
        }){

            // keep a junction if any algorithm claims at least min_n_observed reads
            let max_n_observed = jxns.iter()
                .map(|x| x.n_observed)
                .max()
                .unwrap();

            // keep a junction if all algorithms said it was larger than min_sv_len
            let min_sv_size = jxns.iter()
                .map(|x| x.sv_size)
                .min()
                .unwrap();

            // set a shared id on all fuzzy-matched junctions in a group
            if min_sv_size >= tool.min_sv_len &&
               max_n_observed >= tool.min_n_observed {                
                jxn_id += 1;
                for jxn in jxns {
                    jxn.jxn_id = jxn_id;
                }                
            }
        }
    }
    jxn_id
}

/// Continue fuzzy grouping by matching HiFiRe3 to CuteSV DEL events, where
/// CuteSV encoding differences lead to match failure at large junction offsets.
pub fn rescue_cutesv_dels(
    w: &mut Workflow, 
    tool: &mut ComparisonTool,
    mut jxn_id: u32,
){
    w.log.print("rescuing CuteSV DEL junctions");
    let mut jxns: Vec<_> = tool.jxn_calls.iter_mut()
        .filter(|x| {
            x.source != SvSource::Sniffles &&
            x.sv_type == SvType::Del
        })
        .collect();

    // fuzzy group by sv_type, strands, chrom_1, pos_1
    jxns.sort_unstable_by_key(|jxn| (
        jxn.chrom_1.clone(), jxn.pos1_1
    ));
    for jxns in jxns.chunk_by_mut(|a, b|{
        a.chrom_1 == b.chrom_1 && 
        b.pos1_1 - a.pos1_1 <= tool.match_distance
    }){

        // fuzzy group by chrom_2, pos_2
        jxns.sort_unstable_by_key(|jxn| (
            jxn.chrom_2.clone(), jxn.pos1_2_csv
        ));
        for jxns in jxns.chunk_by_mut(|a, b|{
            a.chrom_2 == b.chrom_2 && 
            b.pos1_2_csv - a.pos1_2_csv <= tool.match_distance
        }){
            if jxns.len() > 1 {
                let mut hf3_jxn_ids: Vec<_> = jxns.iter()
                    .filter_map(|jxn| if jxn.source == SvSource::HiFiRe3 {
                        Some(jxn.jxn_id)
                    } else {
                        None
                    })
                    .collect();
                if hf3_jxn_ids.len() > 0 {
                    if hf3_jxn_ids[0] == 0 {
                        jxn_id += 1;
                        hf3_jxn_ids[0] = jxn_id;
                    }
                    for jxn in jxns {
                        jxn.jxn_id = hf3_jxn_ids[0];
                    }  
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// write a table of admixed SV/junction calls
// -----------------------------------------------------------------------------
fn write_admixed_calls(w: &mut Workflow, tool: &mut ComparisonTool){
    w.log.print("sorting and writing admixed calls");
    // as per write_aggregated_calls, an algorithm can have >1 call
    // report all in the admixed table, in order of decreasing evidence by algorithm
    tool.jxn_calls.sort_unstable_by_key(|jxn| 
        (jxn.jxn_id, jxn.source.clone(), -(jxn.n_observed as i16))
    );
    let csv = OutputCsv::open_env(&mut w.cfg, BENCHMARK_ADMIXED_CALLS, None);
    csv.serialize_all::<JunctionCall>(&tool.jxn_calls);
}

// -----------------------------------------------------------------------------
// write a table of SV/junction calls aggregated by junction group
// -----------------------------------------------------------------------------
fn write_aggregated_calls(w: &mut Workflow, tool: &mut ComparisonTool){
    w.log.print("sorting and writing aggregated calls");
    let mut jxn_calls: Vec<JunctionCallCollapsed> = Vec::new();
    for jxns in tool.jxn_calls.chunk_by(|a, b| {
        a.jxn_id == b.jxn_id
    }) {

        // determine which algorithms called the grouped junction
        // due to grouping differences between the initial callers and this script,
        // it is possible for an algorithm to have more than one call per output group
        let mut by_algo: FxHashMap<SvSource, Vec<&JunctionCall>> = FxHashMap::default();
        for jxn in jxns.into_iter() {
            by_algo
                .entry(jxn.source.clone())
                .or_insert_with(Vec::new)
                .push(jxn);
        }

        // for the single line table, keep the algorithm call with max n_observed
        for jxns in by_algo.values_mut() {
            if jxns.len() > 1 {
                if let Some((arg_max_i, _)) = jxns
                    .iter()
                    .enumerate()
                    .max_by_key(|(_, jxn)| jxn.n_observed)
                {
                    jxns.swap(0, arg_max_i);
                    jxns.truncate(1);
                }
            }
        }

        // assign a reference junction in priority HF3 > Sniffles > CuteSV
        let ref_jxn = if let Some(jxns) = by_algo.get(&SvSource::HiFiRe3) {
            jxns[0]
        } else if let Some(jxns) = by_algo.get(&SvSource::Sniffles) {
            jxns[0]
        } else if let Some(jxns) = by_algo.get(&SvSource::CuteSV) {
            jxns[0]
        } else { continue };

        // collect additional values in an algorithmic fashion
        let mut pos1_2_csv    = 0_u32;
        let mut frag_size_hf3 = 0_u32;
        let mut frag_size_vcf = 0_u32;
        let mut mosaic = false;
        let mut vaf = 0.0;
        let mut algorithms = 0_u8;
        let hf3 = if let Some(jxns) = by_algo.get(&SvSource::HiFiRe3) {
            pos1_2_csv = jxns[0].pos1_2_csv;
            frag_size_hf3 = jxns[0].frag_size;
            algorithms |= AlgorithmOutcomes::HiFiRe3 as u8;
            jxns[0].n_observed
        } else { 0 };
        let csv = if let Some(jxns) = by_algo.get(&SvSource::CuteSV) {
            pos1_2_csv = jxns[0].pos1_2_csv;
            frag_size_vcf = jxns[0].frag_size;
            vaf = jxns[0].vaf; // unless overridden by Sniffles
            algorithms |= AlgorithmOutcomes::CuteSV as u8;
            jxns[0].n_observed
        } else { 0 };
        let snf = if let Some(jxns) = by_algo.get(&SvSource::Sniffles) {
            frag_size_vcf = jxns[0].frag_size; 
            mosaic = jxns[0].mosaic;
            vaf = jxns[0].vaf;
            algorithms |= AlgorithmOutcomes::Sniffles as u8;
            jxns[0].n_observed
        } else { 0 };

        jxn_calls.push(JunctionCallCollapsed { 
            jxn_id:         ref_jxn.jxn_id, 

            sv_type:        ref_jxn.sv_type.clone(), 
            strands:        ref_jxn.strands.clone(), 

            n_algorithms:   by_algo.len() as u8, 
            algorithms,
            hf3, 
            snf, 
            csv, 

            chrom_1:        ref_jxn.chrom_1.clone(), 
            pos1_1:         ref_jxn.pos1_1, 
            chrom_2:        ref_jxn.chrom_2.clone(), 
            pos1_2:         ref_jxn.pos1_2, 
            pos1_2_csv, 

            sv_size:        ref_jxn.sv_size, 
            frag_size_hf3, 
            frag_size_vcf, 

            offset:     ref_jxn.offset, 
            jxn_seq:    ref_jxn.jxn_seq.clone(), 

            mapq:       ref_jxn.mapq, 
            divergence: ref_jxn.divergence, 

            jxn_fail:     ref_jxn.jxn_fail, 
            aln_fail:     ref_jxn.aln_fail, 
            all_chimeric: ref_jxn.all_chimeric, 

            mosaic, 
            vaf, 

            qnames: ref_jxn.qnames.clone(), 
        });
    }

    jxn_calls.sort_unstable_by_key(|jxn| jxn.jxn_id);
    let csv = OutputCsv::open_env(&mut w.cfg, BENCHMARK_AGGREGATED_CALLS, None);
    csv.serialize_all::<JunctionCallCollapsed>(&jxn_calls);
}
