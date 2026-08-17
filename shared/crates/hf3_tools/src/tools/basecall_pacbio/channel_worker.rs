//! Handling strand merging, including decisions not to merge.
//! Generate outcomes for writing and tallying. Operates in parallel 
//! worker threads.

// modules
mod strand_merger;

// dependencies
use std::error::Error;
use crossbeam::{channel::{Receiver, Sender}};
use minimap2::{Aligner as Minimap2, Built};
use faimm::IndexedFasta;
use mdi::pub_key_constants;
use super::{StrandPair, KineticsInstance, MergeResult, MAX_READ_LEN};
use strand_merger::StrandMerger;

// constants
pub_key_constants!(
    // usability reasons
    BOTH_STRANDS_USABLE // all StrandPairs sent via tx_strand_pair have two usable strand 
    //--------------------------
    // merging outcomes
    FAILED_UNMAPPABLE
    TWO_STRANDS_MERGED_WITH_REF
    TWO_STRANDS_MERGED_NO_REF
);

// called by crossbeam scope to create parallel strand-merging worker threads
pub fn merge_strand_pairs(
    rx_strand_pair:  Receiver<StrandPair>,
    tx_kinetics:     Sender<KineticsInstance>,
    tx_corr_to_ref:  Sender<String>,
    tx_merge_result: Sender<MergeResult>,
    minimap2:        &Minimap2<Built>,
    fa:              &IndexedFasta,
) -> Result<(), Box<dyn Error>> {

    // initialize tools for this worker thread
    let mut strand_merger = StrandMerger::new();

    // process strand pairs as they arrive
    for strand_pair in rx_strand_pair.iter() {
        let sources = strand_pair.sources;

        // create ref_on_this map for assigning is_ref to heteroduplex strands during three-stranded error correction
        // along the way, swap this and prev strands to make this == forward whenever possible
        // continue with basecalling even if ref_on_this is not used due to low MAPQ of initial this_on_ref alignment
        //     two-strand error correction will still be performed, with distinct operations in the dd tag
        if let Some(strand_pair) = strand_merger.set_ref_on_this(strand_pair, minimap2, fa, false) {

            // create prev_on_this map for discovering heteroduplex strands
            // abort and report just this strand for SV analysis if strands don't align sufficiently well
            if let Some(failure) = strand_merger.set_prev_on_this(&strand_pair) {
                tx_merge_result.send(MergeResult {
                    sources,
                    reason:   BOTH_STRANDS_USABLE,
                    outcome:  failure,
                    qname:    strand_pair.qname,
                    seq:      strand_pair.this.seq,
                    qual:     strand_pair.this.qual,
                    ff:       strand_pair.ff,
                    ec:       strand_pair.this.ec,
                    dt:       None,
                    dd:       None,
                    sk:       None,
                })?;
            } else {

                // merge the two strands into a single read
                let (seq, qual, dt_tag, dd_tag, sk_tag) = 
                    strand_merger.merge_strands(&strand_pair, &tx_kinetics, &tx_corr_to_ref);

                // transmit the merged read with the two-strand tag set
                tx_merge_result.send(MergeResult {
                    sources,
                    reason: BOTH_STRANDS_USABLE,
                    outcome:  if strand_merger.has_ref_on_this { 
                        TWO_STRANDS_MERGED_WITH_REF 
                    } else { 
                        TWO_STRANDS_MERGED_NO_REF 
                    },
                    qname:    strand_pair.qname,
                    seq:      seq,
                    qual:     qual,
                    ff:       strand_pair.ff,
                    ec:       strand_pair.this.ec + strand_pair.prev.ec,
                    dt:       Some(dt_tag),
                    dd:       Some(dd_tag),
                    sk:       Some(sk_tag),
                })?;
            }
        
        // short-circuit unmapped reads, they are dropped downstream as unusable
        } else {
            tx_merge_result.send(MergeResult {
                sources,
                reason:   BOTH_STRANDS_USABLE,
                outcome:  FAILED_UNMAPPABLE,
                qname:    Vec::new(),
                seq:      String::new(),
                qual:     Vec::new(),
                ff:       0,
                ec:       0.0,
                dt:       None,
                dd:       None,
                sk:       None,
            })?;
        }
    }
    Ok(())
}
