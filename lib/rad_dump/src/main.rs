//! Step 13 RAD readers (port of alevin-fry-0183 `examples/rad_mol_dump.rs` to libradicl 0.21).
//!
//! Derived from alevin-fry (https://github.com/COMBINE-lab/alevin-fry), BSD 3-Clause License,
//! Copyright (c) 2020-, University of Maryland. All rights reserved. The upstream terms apply
//! to the derived portions.
//!
//! Changes from the original: libradicl 0.21 (RAD v2 headers, per-chunk codec taken from the
//! file tags so an lz4-collated RAD is read transparently), sample identity is the raw b0
//! barcode (the caller maps it to a sample from the sample barcode list instead of voting via
//! featureDump), and a second mode reads an uncollated map RAD.
//!
//! mol <collated.rad> <cells16.txt> <out_prefix>
//!     one record per (b0, cell, UMI) molecule of a cell in cells16.txt, gene by the alevin-fry
//!     cr-like rule (per UMI, +1 per distinct gene per read; unique arg-max wins; a tie is AMBIG):
//!         key u64 = (bidx << 56) | (cell << 24) | umi     gene u32 (u32::MAX = AMBIG)
//!     writes <p>.bin <p>.genes.txt <p>.bmap.tsv (b0raw -> bidx, records) <p>.stats.json
//! map <map.rad> <out_prefix>
//!     one (cell u64, gene u32) pair per distinct gene of every mapped record (12 bytes LE)
//!     writes <p>.bin <p>.genes.txt <p>.stats.json
//!
//! Barcodes and UMIs are the 2-bit codes piscem writes (A=0 C=1 G=2 T=3, first base most
//! significant).  Gene ids are the reference name up to the first '|'.

use anyhow::{bail, Context};
use libradicl::header::RadPrelude;
use libradicl::readers::ParallelChunkReader;
use libradicl::record::{MultiBarcodeReadRecord, MultiBarcodeRecordContext};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::num::NonZeroUsize;

const AMBIG: u32 = u32::MAX;

fn code_of(seq: &[u8]) -> anyhow::Result<u64> {
    let mut c = 0u64;
    for &b in seq {
        let v = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => bail!("non-ACGT base in barcode {:?}", String::from_utf8_lossy(seq)),
        };
        c = (c << 2) | v;
    }
    Ok(c)
}

struct Opened {
    reader: BufReader<File>,
    prelude: RadPrelude,
    codec: libradicl::codec::ChunkCodec,
    tid2gid: Vec<u32>,
    gnames: Vec<String>,
}

fn open_rad(path: &str) -> anyhow::Result<Opened> {
    let file = File::open(path).with_context(|| format!("open {path}"))?;
    let mut reader = BufReader::with_capacity(1 << 22, file);
    let prelude = RadPrelude::from_bytes(&mut reader)?;
    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut reader)?;
    let codec = libradicl::codec::chunk_codec_from_tag_map(&file_tag_map)?;
    prelude.get_record_context_prefer_roles::<MultiBarcodeRecordContext>()?;
    let mut gnames = Vec::new();
    let mut gidx: HashMap<String, u32> = HashMap::new();
    let mut tid2gid = Vec::with_capacity(prelude.hdr.ref_count as usize);
    for rn in &prelude.hdr.ref_names {
        let g = rn.split('|').next().unwrap_or("").to_string();
        let id = *gidx.entry(g.clone()).or_insert_with(|| {
            gnames.push(g.clone());
            (gnames.len() - 1) as u32
        });
        tid2gid.push(id);
    }
    eprintln!("{path}: refs {} genes {} codec {}", tid2gid.len(), gnames.len(), codec.as_str());
    Ok(Opened { reader, prelude, codec, tid2gid, gnames })
}

fn write_genes(prefix: &str, g: &[String]) -> anyhow::Result<()> {
    let mut w = BufWriter::new(File::create(format!("{prefix}.genes.txt"))?);
    for n in g {
        writeln!(w, "{n}")?;
    }
    Ok(w.flush()?)
}

fn mol(rad: &str, cellf: &str, prefix: &str) -> anyhow::Result<()> {
    let mut cells: HashSet<u64> = HashSet::new();
    for line in BufReader::new(File::open(cellf)?).lines() {
        let l = line?;
        let t = l.trim();
        if t.is_empty() {
            continue;
        }
        if t.len() != 16 {
            bail!("cells16 entry is not 16 nt: {t:?}");
        }
        cells.insert(code_of(t.as_bytes())?);
    }
    if cells.is_empty() {
        bail!("empty cell set");
    }
    let Opened { mut reader, prelude, codec, tid2gid, gnames } = open_rad(rad)?;
    write_genes(prefix, &gnames)?;
    let mut w = BufWriter::with_capacity(1 << 24, File::create(format!("{prefix}.bin"))?);
    let mut chunk_reader = ParallelChunkReader::<MultiBarcodeReadRecord>::new(&prelude, NonZeroUsize::new(1).unwrap())
        .with_chunk_codec(codec);

    type Res = (u64, u64, u64, u64, u64, Vec<u64>, Vec<u64>);
    let (nrec, nchunk, nmol, nemit, nambig, b0s, b0recs) = std::thread::scope(|scope| -> anyhow::Result<Res> {
        let chunks = chunk_reader.chunk_iter();
        let h = scope.spawn(|| -> anyhow::Result<Res> {
            let (mut nrec, mut nchunk, mut nmol, mut nemit, mut nambig) = (0u64, 0u64, 0u64, 0u64, 0u64);
            let mut bmap: HashMap<u64, u8> = HashMap::new();
            let mut b0s: Vec<u64> = Vec::new();
            let mut b0recs: Vec<u64> = Vec::new();
            let mut votes: HashMap<u64, HashMap<u32, u32>> = HashMap::new();
            let mut gset: Vec<u32> = Vec::new();
            let mut buf: Vec<u8> = Vec::with_capacity(1 << 24);
            for meta_chunk in chunks {
                for chunk in meta_chunk.iter() {
                    nchunk += 1;
                    votes.clear();
                    let mut key_cell: Option<(u64, u64)> = None;
                    for record in &chunk.reads {
                        nrec += 1;
                        if record.barcodes.len() != 2 {
                            bail!("record has {} barcodes, expected 2 (sample, cell)", record.barcodes.len());
                        }
                        let bc = (record.barcodes[0], record.barcodes[1]);
                        match key_cell {
                            None => key_cell = Some(bc),
                            Some(k) if k != bc => bail!("collated chunk mixes (b0, cell) {k:?} and {bc:?}"),
                            _ => {}
                        }
                        if record.refs.is_empty() {
                            continue;
                        }
                        gset.clear();
                        for &t in &record.refs {
                            gset.push(tid2gid[t as usize]);
                        }
                        gset.sort_unstable();
                        gset.dedup();
                        let e = votes.entry(record.umi).or_default();
                        for &g in &gset {
                            *e.entry(g).or_insert(0) += 1;
                        }
                    }
                    let Some((b0, cell)) = key_cell else { continue };
                    let bidx = *bmap.entry(b0).or_insert_with(|| {
                        b0s.push(b0);
                        b0recs.push(0);
                        (b0s.len() - 1) as u8
                    });
                    b0recs[bidx as usize] += chunk.reads.len() as u64;
                    if b0s.len() > 255 {
                        bail!("more than 255 distinct b0 values");
                    }
                    let emit = cells.contains(&cell);
                    for (&umi, gm) in votes.iter() {
                        nmol += 1;
                        let (mut best_g, mut best_c, mut tie) = (0u32, 0u32, false);
                        for (&g, &c) in gm.iter() {
                            if c > best_c {
                                best_c = c;
                                best_g = g;
                                tie = false;
                            } else if c == best_c {
                                tie = true;
                            }
                        }
                        if tie {
                            nambig += 1;
                        }
                        if !emit {
                            continue;
                        }
                        let key = ((bidx as u64) << 56) | (cell << 24) | umi;
                        buf.extend_from_slice(&key.to_le_bytes());
                        buf.extend_from_slice(&(if tie { AMBIG } else { best_g }).to_le_bytes());
                        nemit += 1;
                        if buf.len() >= (1 << 23) {
                            w.write_all(&buf)?;
                            buf.clear();
                        }
                    }
                }
            }
            w.write_all(&buf)?;
            w.flush()?;
            Ok((nrec, nchunk, nmol, nemit, nambig, b0s, b0recs))
        });
        chunk_reader.start(&mut reader, None::<fn(u64, u64)>)?;
        h.join().map_err(|_| anyhow::anyhow!("consumer thread panicked"))?
    })?;
    let mut bw = BufWriter::new(File::create(format!("{prefix}.bmap.tsv"))?);
    writeln!(bw, "b0raw\tbidx\trecords")?;
    for (i, b0) in b0s.iter().enumerate() {
        writeln!(bw, "{b0}\t{i}\t{}", b0recs[i])?;
    }
    bw.flush()?;
    let mut sw = File::create(format!("{prefix}.stats.json"))?;
    writeln!(sw, "{{\"num_records\": {nrec}, \"num_chunks\": {nchunk}, \"num_molecules_all_cells\": {nmol}, \"num_emitted_in_cellset\": {nemit}, \"num_ambiguous_all_cells\": {nambig}, \"num_genes\": {}}}", gnames.len())?;
    eprintln!("mol: records {nrec} chunks {nchunk} molecules {nmol} emitted {nemit} ambiguous {nambig}");
    Ok(())
}

fn map(rad: &str, prefix: &str) -> anyhow::Result<()> {
    let Opened { mut reader, prelude, codec, tid2gid, gnames } = open_rad(rad)?;
    write_genes(prefix, &gnames)?;
    let mut w = BufWriter::with_capacity(1 << 24, File::create(format!("{prefix}.bin"))?);
    let mut chunk_reader = ParallelChunkReader::<MultiBarcodeReadRecord>::new(&prelude, NonZeroUsize::new(1).unwrap())
        .with_chunk_codec(codec);
    let (nrec, npairs) = std::thread::scope(|scope| -> anyhow::Result<(u64, u64)> {
        let chunks = chunk_reader.chunk_iter();
        let h = scope.spawn(|| -> anyhow::Result<(u64, u64)> {
            let (mut nrec, mut npairs) = (0u64, 0u64);
            let mut gset: Vec<u32> = Vec::new();
            let mut buf: Vec<u8> = Vec::with_capacity(1 << 24);
            for meta_chunk in chunks {
                for chunk in meta_chunk.iter() {
                    for record in &chunk.reads {
                        nrec += 1;
                        let cell = *record.barcodes.last().context("record without barcodes")?;
                        gset.clear();
                        for &t in &record.refs {
                            gset.push(tid2gid[t as usize]);
                        }
                        gset.sort_unstable();
                        gset.dedup();
                        if gset.is_empty() {
                            bail!("mapped record with no references");
                        }
                        for &g in &gset {
                            buf.extend_from_slice(&cell.to_le_bytes());
                            buf.extend_from_slice(&g.to_le_bytes());
                            npairs += 1;
                        }
                        if buf.len() >= (1 << 23) {
                            w.write_all(&buf)?;
                            buf.clear();
                        }
                    }
                }
            }
            w.write_all(&buf)?;
            w.flush()?;
            Ok((nrec, npairs))
        });
        chunk_reader.start(&mut reader, None::<fn(u64, u64)>)?;
        h.join().map_err(|_| anyhow::anyhow!("consumer thread panicked"))?
    })?;
    let mut sw = File::create(format!("{prefix}.stats.json"))?;
    writeln!(sw, "{{\"num_records\": {nrec}, \"num_pairs\": {npairs}, \"num_genes\": {}}}", gnames.len())?;
    eprintln!("map: records {nrec} (read, gene) pairs {npairs}");
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    match a.first().map(String::as_str) {
        Some("mol") if a.len() == 4 => mol(&a[1], &a[2], &a[3]),
        Some("map") if a.len() == 3 => map(&a[1], &a[2]),
        _ => bail!("usage: rad_dump mol <collated.rad> <cells16.txt> <out_prefix> | map <map.rad> <out_prefix>"),
    }
}
