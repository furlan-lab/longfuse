use anyhow::{Context, Result};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::{HashMap, HashSet};

/// A simple structure to hold information about one alignment
/// from one long read (primary or supplementary).
#[derive(Debug, Clone)]
struct SplitAlignment {
    contig: String,
    start: i64,
    end: i64,
    strand: bool, // false -> negative, true -> positive
}

/// A fusion event might contain two or more splits from distinct genomic loci.
#[derive(Debug)]
struct FusionEvent {
    read_id: String,
    splits: Vec<SplitAlignment>,
}

/// Parse a single BAM record, returning all its alignments (primary + supplementary).
fn parse_record(record: &bam::Record, tid2name: &[String]) -> Vec<(String, SplitAlignment)> {
    let mut result = Vec::new();
    
    // If the record is unmapped or has an invalid TID, skip
    if record.tid() < 0 {
        return result;
    }

    let read_id = String::from_utf8_lossy(record.qname()).to_string();
    let tid = record.tid() as usize;
    let contig = tid2name[tid].clone();
    let start = record.pos();
    let end = record.calculate_end() as i64;
    let strand = !record.is_reverse();

    // Primary or supplementary alignment itself
    if record.is_primary() || record.is_supplementary() {
        let split_alignment = SplitAlignment {
            contig: contig.clone(),
            start,
            end,
            strand,
        };
        result.push((read_id.clone(), split_alignment));
    }

    // Parse the SA tag for supplementary alignments
    // SA tag format: contig, pos, strand, CIGAR, mapQ, NM;
    if let Ok(aux) = record.aux(b"SA") {
        // If the tag is not a string, skip
        let sa_bytes = match aux.string() {
            Ok(bytes) => bytes,
            Err(_) => return result,
        };
        let sa_str = String::from_utf8_lossy(sa_bytes);
        let sa_entries: Vec<&str> = sa_str.split(';').filter(|x| !x.is_empty()).collect();
        for sa in sa_entries {
            let fields: Vec<&str> = sa.split(',').collect();
            if fields.len() < 6 {
                continue;
            }
            let sa_contig = fields[0].to_string();
            let sa_pos = fields[1].parse::<i64>().unwrap_or(-1);
            let sa_strand = fields[2] == "+";
            // For a more accurate end, you'd parse the CIGAR string in fields[3].
            // Here, we only record the start and a placeholder end.
            let sa_split = SplitAlignment {
                contig: sa_contig,
                start: sa_pos - 1, // SA pos is 1-based
                end: sa_pos,       // naive placeholder; real code should parse CIGAR
                strand: sa_strand,
            };
            result.push((read_id.clone(), sa_split));
        }
    }
    result
}

/// Determine if a set of splits for a given read constitutes a potential fusion.
fn detect_fusion(read_id: &str, splits: &[SplitAlignment]) -> Option<FusionEvent> {
    if splits.len() < 2 {
        return None;
    }

    // Check distinct contigs
    let distinct_contigs: HashSet<&str> = splits.iter().map(|s| s.contig.as_str()).collect();
    if distinct_contigs.len() > 1 {
        // Potential fusion if it spans multiple contigs
        return Some(FusionEvent {
            read_id: read_id.to_string(),
            splits: splits.to_vec(),
        });
    } else {
        // If on the same contig, check for large distances (example threshold = 1Mb)
        let mut suspicious = false;
        for i in 0..splits.len() {
            for j in i + 1..splits.len() {
                let dist = (splits[i].start - splits[j].start).abs();
                if dist > 1_000_000 {
                    suspicious = true;
                    break;
                }
            }
            if suspicious {
                break;
            }
        }
        if suspicious {
            return Some(FusionEvent {
                read_id: read_id.to_string(),
                splits: splits.to_vec(),
            });
        }
    }
    None
}

fn main() -> Result<()> {
    // Path to your BAM file
    let bam_path = "input_long_reads.bam";

    // Open the BAM file
    let mut bam_reader = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;

    // Build a vec of contig names to map tid -> contig
    let header = bam_reader.header().clone();
    let tid2name: Vec<String> = header
        .target_names()
        .iter()
        .map(|bytes| String::from_utf8_lossy(bytes).to_string())
        .collect();

    // ----------------------------------------------------------------------
    // 1) Read all records into memory (serial). For large files, consider streaming.
    // ----------------------------------------------------------------------
    let mut all_records = Vec::new();
    for rec_res in bam_reader.records() {
        let record = rec_res?;
        all_records.push(record);
    }
    println!("Loaded {} records from BAM.", all_records.len());

    // ----------------------------------------------------------------------
    // 2) Extract primary + supplementary alignments in parallel
    //    Each record can yield multiple (read_id, SplitAlignment) pairs.
    // ----------------------------------------------------------------------
    let mut all_splits: Vec<(String, SplitAlignment)> = all_records
        .par_iter()
        .flat_map(|record| parse_record(record, &tid2name))
        .collect();

    // ----------------------------------------------------------------------
    // 3) Group by read_id
    //    (Serial grouping for simplicity; you could also parallelize with something like
    //    a DashMap or a parallel grouping approach, but this is usually fast enough.)
    // ----------------------------------------------------------------------
    let mut split_map: HashMap<String, Vec<SplitAlignment>> = HashMap::new();
    for (read_id, split_align) in all_splits.drain(..) {
        split_map
            .entry(read_id)
            .or_insert_with(Vec::new)
            .push(split_align);
    }

    // ----------------------------------------------------------------------
    // 4) Detect potential fusion events in parallel.
    //    We parallelize over the entries of the HashMap.
    // ----------------------------------------------------------------------
    let fusion_events: Vec<FusionEvent> = split_map
        .par_iter()
        .filter_map(|(read_id, splits)| detect_fusion(read_id, splits))
        .collect();

    // ----------------------------------------------------------------------
    // 5) Print out detected events
    // ----------------------------------------------------------------------
    for event in &fusion_events {
        println!("Potential Fusion Event for read: {}", event.read_id);
        for (i, split) in event.splits.iter().enumerate() {
            println!(
                "  Split {}: {}:{}-{} strand={}",
                i + 1,
                split.contig,
                split.start,
                split.end,
                if split.strand { "+" } else { "-" }
            );
        }
    }

    println!("Total fusion events detected: {}", fusion_events.len());
    Ok(())
}
