//! Handle the input BAM stream; dispatch strand pairs for merging.

// modules
pub mod usable;

// dependencies
use std::error::Error;
use std::path::Path;
use rustc_hash::FxHashMap;
use std::str::from_utf8_unchecked;
use crossbeam::channel::Sender;
use rust_htslib::{bam::{Read, Reader, Record as BamRecord}};
use mdi::pub_key_constants;
use mdi::workflow::Config;
use genomex::bam::tags;
use crate::formats::hf3_tags::*;
use super::{BufferedStrand, StrandPair, MergeResult};
use usable::{is_usable, UnusableReason, StrandSource};

/// StrandBuffer structure for caching PacBio strands
/// while waiting for their matching strand.
struct StrandBuffer {
    movies: Vec<String>,
    pub strand_buffer: FxHashMap<usize, BufferedStrand>,  
}
impl StrandBuffer {
    /// Create a new, empty StrandBuffer structure.
    pub fn new() -> Self {
        StrandBuffer {
            movies: Vec::new(),
            strand_buffer: FxHashMap::default(),
        }
    }
    /// Convert a movie name and ZMW hole number into a unique usize strand buffer key.
    pub fn get_key(&mut self, movie: &str, zmw: usize) -> usize {
        let movie_idx = if let Some(idx) = self.movies.iter().position(|m| m == movie) {
            idx
        } else {
            self.movies.push(movie.to_string());
            self.movies.len() - 1
        };
        movie_idx << 32 | zmw
    }
    /// Push a first encountered strand into the buffer for later merging.
    pub fn insert(
        &mut self, 
        key:    usize, 
        source: StrandSource,
        usable: bool, 
        reason: UnusableReason,
        ff:     u8, 
        ec:     f32, 
        strand: &BamRecord
    ) {
        self.strand_buffer.insert(key, BufferedStrand {
            source: source,
            usable,
            reason,
            ff,
            ec,
            seq:      strand.seq().as_bytes().iter().map(|&c| c as char).collect(),
            qual:     strand.qual().to_vec(),
            ip:       tags::get_tag_u8_vec_opt(strand, INTER_PULSE_DURATION),
            pw:       tags::get_tag_u8_vec_opt(strand, PULSE_WIDTH),
        });
    }
    /// Attempt to remove a first strand from the buffer for merging to a second strand.
    pub fn remove(&mut self, key: &usize) -> Option<BufferedStrand> {
        self.strand_buffer.remove(&key)
    }
}

// constants
pub_key_constants!(
    // from environment variables
    FAIL_BAM_FILE
    HIFI_BAM_FILE
    // read outcomes
    ONE_USABLE_STRAND
    BOTH_STRANDS_UNUSABLE
    ORPHAN_STRAND
);

// called by crossbeam scope to create the single input thread
pub fn stream_bam_files(
    tx_strand_pair:  Sender<StrandPair>,
    tx_merge_result: Sender<MergeResult>,
) -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[FAIL_BAM_FILE, HIFI_BAM_FILE]);

    // instantiate buffering
    let mut this_strand = BamRecord::new();
    let mut strand_buffer = StrandBuffer::new();

    // scan the fail and hifi BAM files in order to build a cache of by-strand reads
    // merge strands when both are encountered and usable
    for bam_file_key in &[FAIL_BAM_FILE, HIFI_BAM_FILE] {
    // for bam_file_key in &[HIFI_BAM_FILE] {
        let bam_file = cfg.get_string(bam_file_key);
        if bam_file == "NA" { continue } // make FAIL_BAM_FILE optional
        let bam_file_name = Path::new(&bam_file).file_name().unwrap().to_str().unwrap();
        eprintln!("processing {}", bam_file_name);
        let mut bam = Reader::from_path(&bam_file)?;
        let source = if *bam_file_key == FAIL_BAM_FILE {
            StrandSource::Fail
        } else {
            StrandSource::Hifi
        };

        //////////////////////////////////////
        // let mut counter: usize = 0;

        while let Some(r) = bam.read(&mut this_strand) {
            r.expect("Failed to parse BAM record");
            process_strand(
                source,
                &this_strand, 
                &mut strand_buffer, 
                &tx_strand_pair, 
                &tx_merge_result,

                // &mut counter,

            )?;
            
            // ////////////////////////////////////////////////////////
            // if counter >= 10000 {
            //     break; // TODO: remove this break after testing
            // }
        }
    }
    
    // report any residual orphan strands; hopefully there are none
    for (_key, strand) in strand_buffer.strand_buffer.into_iter() {
        tx_merge_result.send(MergeResult{
            sources: strand.source,
            reason:  strand.reason.to_str(),
            outcome: ORPHAN_STRAND,
            qname: Vec::new(), // can't parse qname without read data from 'this' strand
            seq:   String::new(),
            qual:  Vec::new(),
            ff: 0,
            ec: 0.0,
            dt: None,
            dd: None,
            sk: None,
        })?;
    }
    Ok(())
}

// dispatch one strand of a PacBio by-strand read to either the StrandBuffer
// or for processing with its partner strand
fn process_strand(
    this_source:       StrandSource,
    this_strand:       &BamRecord,
    strand_buffer:     &mut StrandBuffer,
    tx_strand_pair:    &Sender<StrandPair>,
    tx_merge_result:   &Sender<MergeResult>,
    // counter:         &mut usize,
) -> Result<(), Box<dyn Error>> {

    // parse ZMW from read name: m64011_190714_120746/14/ccs/rev
    // fail if the QNAME or hole number can't be parsed
    let qname = unsafe { from_utf8_unchecked(this_strand.qname()) };
    let parts: Vec<&str> = qname.split('/').collect();
    if parts.len() < 4 { 
        return Err(format!("Unexpected name format: {}", qname).into());
    }
    let zmw = parts[1].parse::<u32>().map_err(|e|
        format!("Failed to parse ZMW from read name {}: {}", qname, e)
    )?;

    // assemble a buffer key for this strand's ZMW that is unique across movies
    let buffer_key = strand_buffer.get_key(parts[0], zmw as usize);

    // assess whether this strand is usable for merging
    let (this_usable, this_reason, this_ff, this_ec) = is_usable(this_strand);

    // process the second strand of a matching by-strand read pair
    // remove the cached data from the map for memory management
    if let Some(prev_strand) = strand_buffer.remove(&buffer_key){
        let qname = format!("{}/{}/ccs", parts[0], zmw).into_bytes();
        let ff = prev_strand.ff | this_ff; // read reports both strands' fail flags
        let this_is_hifi = this_source == StrandSource::Hifi;

        // merge two usable strands into a single read if this strand is a hifi strand
        // used downstream for both SNV and SV detection
        if prev_strand.usable && this_usable && this_is_hifi { // a fail strand can be used to error correct a hifi strand but not another fail strand
            tx_strand_pair.send(StrandPair{
                sources: this_source.get_sources(&prev_strand.source),
                qname,
                ff,
                this: BufferedStrand {
                    source: this_source,
                    usable: this_usable,
                    reason: this_reason,
                    ff,
                    ec:   this_ec,
                    seq:  this_strand.seq().as_bytes().iter().map(|&c| c as char).collect(),
                    qual: this_strand.qual().to_vec(),
                    ip:   tags::get_tag_u8_vec_opt(this_strand, INTER_PULSE_DURATION),
                    pw:   tags::get_tag_u8_vec_opt(this_strand, PULSE_WIDTH),
                },
                prev: prev_strand,
            })?;
            // *counter += 1;

        // transmit reads with only one usable strand directly to writer without 
        // the two-strand tag used downstream for SV detection only
        } else if this_usable && this_is_hifi { // single strands reported for SV analysis must be hifi
            tx_merge_result.send(MergeResult{
                sources: prev_strand.source, // the source and reason for the other unusable strand
                reason:  prev_strand.reason.to_str(),
                outcome: ONE_USABLE_STRAND,
                qname,
                seq:  this_strand.seq().as_bytes().iter().map(|&c| c as char).collect(),
                qual: this_strand.qual().to_vec(),
                ff,
                ec: this_ec,
                dt: None,
                dd: None,
                sk: None,
            })?;
        } else if prev_strand.usable && prev_strand.source == StrandSource::Hifi {
            tx_merge_result.send(MergeResult{
                sources: this_source,
                reason:  this_reason.to_str(),
                outcome: ONE_USABLE_STRAND,
                qname,
                seq:  prev_strand.seq,
                qual: prev_strand.qual,
                ff,
                ec:   prev_strand.ec,
                dt:   None,
                dd:   None,
                sk:   None,
            })?;

        // reject reads with no usable strands due to adapter parsing failures, low coverage, etc.
        } else {
            tx_merge_result.send(MergeResult{
                sources: this_source.get_sources(&prev_strand.source),
                reason:  this_reason.get_reasons(&prev_strand.reason),
                outcome: BOTH_STRANDS_UNUSABLE,
                qname,
                seq: String::new(),
                qual: Vec::new(),
                ff,
                ec: 0.0,
                dt: None,
                dd: None,
                sk: None,
            })?;
        }

    // cache the needed bits of the first encountered strand of a ZMW while waiting for the second
    } else {
        strand_buffer.insert(
            buffer_key, 
            this_source,
            this_usable, 
            this_reason,
            this_ff, 
            this_ec, 
            this_strand
        );
    }
    Ok(())
}
