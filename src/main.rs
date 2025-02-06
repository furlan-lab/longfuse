use anyhow::{Context, Result};
use clap::{Parser, ArgAction, command};
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, Read};
use rust_htslib::bam::record::Cigar;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write};

/// A single alignment segment (primary or supplementary) from one read.
#[derive(Debug, Clone)]
struct SplitAlignment {
    contig: String,
    start: i64, // 0-based
    end: i64,   // 1-based
    strand: bool, // true = forward, false = reverse
}

/// A naive "fusion event": a read with multiple distinct splits.
#[derive(Debug)]
struct FusionEvent {
    read_id: String,
    splits: Vec<SplitAlignment>,
}

/// Command-line arguments. 
#[derive(Parser, Debug)]
#[command(
    name = "longfuse_strict_filters",
    version = "0.1.0",
    about = "Rust CLI for detecting naive fusions with optional strict filters."
)]
struct Cli {
    /// Path to the input BAM file
    #[arg(long = "bam", short = 'b')]
    bam_path: String,

    /// Number of threads for Rayon
    #[arg(long = "threads", short = 't', default_value_t = 1)]
    threads: usize,

    /// Output file (if omitted, prints to STDOUT)
    #[arg(long = "output", short = 'o')]
    output: Option<String>,

    /// Distance threshold for intra-contig splits to be called suspicious
    #[arg(long = "threshold", default_value_t = 1_000_000)]
    threshold: i64,

    /// Minimum MAPQ below which an alignment is discarded
    #[arg(long = "min-mapq", default_value_t = 20)]
    min_mapq: u8,

    /// Minimum aligned length below which an alignment is discarded
    #[arg(long = "min-length", default_value_t = 100)]
    min_length: u32,

    /// Whether to allow random/unplaced contigs (default false = exclude)
    #[arg(long = "allow-random", action = ArgAction::SetFalse, default_value_t = false)]
    allow_random: bool,

    /// Whether to skip secondary alignments (default true)
    #[arg(long = "skip-secondary", action = ArgAction::SetTrue, default_value_t = true)]
    skip_secondary: bool,

    /// Whether to deduplicate repeated (contig, start, end, strand) splits (default true)
    #[arg(long = "deduplicate", action = ArgAction::SetTrue, default_value_t = true)]
    deduplicate: bool,
}

/// Return the total aligned length (in bp) of a `bam::Record` from its CIGAR.
fn aligned_length_from_cigar(record: &bam::Record) -> u32 {
    record
        .cigar()
        .iter()
        .filter_map(|op| match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => Some(*len),
            _ => None,
        })
        .sum()
}

/// Parse the CIGAR string from an SA tag entry (like "50M2S") to get its aligned length.
fn aligned_length_from_cigar_str(cigar_str: &str) -> u32 {
    // Very naive parser: sum up the digits for M/EQ/X.
    // Real code might parse carefully with regex or a small parser.
    // E.g. "100M10S" => 100 aligned
    // This is a quick approach:
    use regex::Regex;
    lazy_static::lazy_static! {
        static ref RE: Regex = Regex::new(r"(\d+)([MIDNSHP=X])").unwrap();
    }
    let mut sum = 0;
    for cap in RE.captures_iter(cigar_str) {
        let length: u32 = cap[1].parse().unwrap_or(0);
        let op = &cap[2];
        match op {
            "M" | "=" | "X" => sum += length,
            // Diff, etc.
            _ => {}
        }
    }
    sum
}

/// Check if a contig is considered "random"/"unplaced".
fn is_random_contig(contig: &str) -> bool {
    contig.contains("_random")
        || contig.starts_with("chrUn_")
        || contig.starts_with("chrUn")
        || contig.contains("_KI270")
        || contig.starts_with("chrUn-")
}

/// Gathers primary + supplementary alignments (optionally filtering out
/// unwanted alignments based on user CLI options).
fn parse_record(
    record: &bam::Record,
    tid2name: &[String],
    cli: &Cli,
) -> Vec<(String, SplitAlignment)> {
    let mut segments = Vec::new();

    // If unmapped or TID invalid
    if record.tid() < 0 {
        return segments;
    }

    // If skipping secondary alignments
    if cli.skip_secondary && record.is_secondary() {
        return segments;
    }

    // Check record's MAPQ
    let mapq = record.mapq();
    if mapq < cli.min_mapq {
        return segments;
    }

    // Determine contig name
    let contig = &tid2name[record.tid() as usize];

    // Exclude random/unplaced if allow_random is false
    if !cli.allow_random && is_random_contig(contig) {
        return segments;
    }

    // Compute aligned length from the CIGAR
    let aln_len = aligned_length_from_cigar(record);
    if aln_len < cli.min_length {
        return segments;
    }

    let read_id = String::from_utf8_lossy(record.qname()).to_string();
    let start_0 = record.pos() as i64;
    let end_1 = record.cigar().end_pos() as i64; // 1-based
    let strand = !record.is_reverse();

    // "Primary" means not secondary nor supplementary
    let is_primary = !record.is_secondary() && !record.is_supplementary();
    // If it's either primary or supplementary, store it
    if is_primary || record.is_supplementary() {
        let split = SplitAlignment {
            contig: contig.clone(),
            start: start_0,
            end: end_1,
            strand,
        };
        segments.push((read_id.clone(), split));
    }

    // Parse SA tag if present
    if let Ok(Aux::String(sa_str)) = record.aux(b"SA") {
        // Multiple alignments separated by ';'
        // Each alignment: contig, pos, strand, CIGAR, mapQ, NM
        for sa_entry in sa_str.split(';').filter(|s| !s.is_empty()) {
            let fields: Vec<&str> = sa_entry.split(',').collect();
            if fields.len() < 6 {
                continue;
            }
            let sa_contig = fields[0].to_string();
            let sa_pos = fields[1].parse::<i64>().unwrap_or(-1);
            let sa_strand = fields[2] == "+";
            let sa_cigar = fields[3];
            let sa_mapq = fields[4].parse::<u8>().unwrap_or(0);

            // Filter by MAPQ
            if sa_mapq < cli.min_mapq {
                continue;
            }

            // Exclude random contigs if not allowed
            if !cli.allow_random && is_random_contig(&sa_contig) {
                continue;
            }

            // Compute aligned length from the SA CIGAR
            let sa_aln_len = aligned_length_from_cigar_str(sa_cigar);
            if sa_aln_len < cli.min_length as u32 {
                continue;
            }

            // 1-based => 0-based start, naive end
            // If you want more precise end, parse the CIGAR carefully.
            let sa_start = sa_pos - 1;
            let sa_end = sa_pos; // naive

            let sa_split = SplitAlignment {
                contig: sa_contig,
                start: sa_start,
                end: sa_end,
                strand: sa_strand,
            };
            segments.push((read_id.clone(), sa_split));
        }
    }

    segments
}

/// Optionally deduplicate repeated `(read_id, contig, start, end, strand)`.
fn deduplicate_splits(pairs: Vec<(String, SplitAlignment)>) -> Vec<(String, SplitAlignment)> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for (rid, split) in pairs {
        let key = (
            rid.clone(),
            split.contig.clone(),
            split.start,
            split.end,
            split.strand,
        );
        if seen.insert(key) {
            out.push((rid, split));
        }
    }
    out
}

/// Determine if a read's splits indicate a potential fusion:
/// - Multiple distinct contigs => fusion.
/// - Same contig but separated by > threshold => fusion.
fn detect_fusion(
    read_id: &str,
    splits: &[SplitAlignment],
    threshold: i64,
) -> Option<FusionEvent> {
    if splits.len() < 2 {
        return None;
    }
    let distinct_contigs: HashSet<_> = splits.iter().map(|s| s.contig.as_str()).collect();
    if distinct_contigs.len() > 1 {
        // multi-chrom => call it a fusion
        return Some(FusionEvent {
            read_id: read_id.to_string(),
            splits: splits.to_vec(),
        });
    }
    // If same contig, check distance
    for i in 0..splits.len() {
        for j in i + 1..splits.len() {
            let dist = (splits[i].start - splits[j].start).abs();
            if dist > threshold {
                return Some(FusionEvent {
                    read_id: read_id.to_string(),
                    splits: splits.to_vec(),
                });
            }
        }
    }
    None
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Set Rayon thread pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .expect("Failed to configure Rayon");

    // Open the BAM file
    let mut reader = bam::Reader::from_path(&cli.bam_path)
        .with_context(|| format!("Could not open BAM: {}", &cli.bam_path))?;

    // (Optional) htslib multi-threaded BGZF reading:
    // reader.set_threads(cli.threads).unwrap_or(());

    // Contig names
    let header = reader.header().clone();
    let tid2name: Vec<String> = header
        .target_names()
        .iter()
        .map(|bs| String::from_utf8_lossy(bs).to_string())
        .collect();

    // Read everything (serial)
    let mut records = Vec::new();
    for rec_res in reader.records() {
        records.push(rec_res?);
    }
    eprintln!("Loaded {} records.", records.len());

    // Parallel-extract splits
    let mut splitted: Vec<(String, SplitAlignment)> = records
        .par_iter()
        .flat_map(|r| parse_record(r, &tid2name, &cli))
        .collect();

    // Deduplicate if enabled
    if cli.deduplicate {
        splitted = deduplicate_splits(splitted);
    }

    // Group by read_id
    let mut split_map: HashMap<String, Vec<SplitAlignment>> = HashMap::new();
    for (rid, split) in splitted {
        split_map.entry(rid).or_insert_with(Vec::new).push(split);
    }

    // Detect fusions in parallel
    let fusions: Vec<FusionEvent> = split_map
        .par_iter()
        .filter_map(|(rid, splits)| detect_fusion(rid, splits, cli.threshold))
        .collect();

    // Output
    let mut writer: Box<dyn Write> = if let Some(path) = cli.output.as_ref() {
        Box::new(File::create(path).with_context(|| format!("Cannot create {}", path))?)
    } else {
        Box::new(io::stdout())
    };

    for event in &fusions {
        writeln!(writer, "Fusion Event in read: {}", event.read_id)?;
        for (i, split) in event.splits.iter().enumerate() {
            writeln!(
                writer,
                "  Split {} => {}:{}-{} strand={}",
                i + 1,
                split.contig,
                split.start,
                split.end,
                if split.strand { "+" } else { "-" }
            )?;
        }
        writeln!(writer)?;
    }
    writeln!(writer, "Total fusions detected: {}", fusions.len())?;
    Ok(())
}
