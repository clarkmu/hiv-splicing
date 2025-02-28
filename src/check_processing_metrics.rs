use std::fs::{ read_to_string, File };
use std::io::{ self, BufReader, Read };
use std::time::Instant;
use peak_alloc::PeakAlloc;
use bio::io::fastq::Reader;
use bio::io::fasta::Record;
use rayon::prelude::*;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn process_record(line: &str) {
    // Simulate processing a line
    // println!("Processing record: {}", line);
}

pub fn check_processing_metrics(file_path: &str) -> io::Result<()> {
    // Measure time for reading with chunking
    let start = Instant::now();
    let mut file = File::open(file_path)?;

    let mut chunk = [0; 1024];
    let current_mem = PEAK_ALLOC.current_usage_as_mb();
    println!("RAM chunking {:.2}MB", current_mem);
    while file.read(&mut chunk)? != 0 {
        let lines = String::from_utf8_lossy(&chunk);
        for line in lines.lines() {
            process_record(line);
        }

        // let records = Reader::from_bufread(&chunk);
        // for record in records.records() {
        //     println!("{:?}", record.id());
        // }
    }
    let duration = start.elapsed();
    println!("Time chunking: {:.2}", duration.as_secs_f32());
    drop(file);

    let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    println!("The max amount that has been used {}GB", peak_mem);
    // return Ok(());

    // Measure time for reading into memory
    let start = Instant::now();
    let lines = read_to_string(file_path)?;
    let current_mem = PEAK_ALLOC.current_usage_as_mb();
    println!("RAM read_to_string {:.2}MB", current_mem);
    for line in lines.lines() {
        process_record(line);
    }
    let duration = start.elapsed();
    println!("Time read_to_string: {:.2}", duration.as_secs_f64());
    drop(lines);

    let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    println!("The max amount that has been used {}GB", peak_mem);

    // Measure time for reading with rust-bio
    let start = Instant::now();
    let fasta_reader = Reader::from_file(file_path).unwrap();
    let current_mem = PEAK_ALLOC.current_usage_as_mb();
    println!("RAM rust-bio {:.2}MB", current_mem);

    // let records = fasta_reader.records().collect::<Vec<_>>();
    // println!("Number of records: {}", records.count());

    for result in fasta_reader.records() {
        match result {
            Ok(result_data) => {
                // println!("Record: {:?}", result_data.id());
                process_record(result_data.id());
            }
            Err(e) => {
                println!("Error reading record {:?}", e);
            }
        }
    }
    let duration = start.elapsed();
    println!("Time rust-bio: {:.2}", duration.as_secs_f32());
    // drop(fasta_reader);

    let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    println!("The max amount that has been used {}GB", peak_mem);

    // Measure time for reading with rust-bio
    let start = Instant::now();
    let fasta_reader = Reader::from_file(file_path).unwrap();
    let current_mem = PEAK_ALLOC.current_usage_as_mb();
    println!("RAM rust-bio into_par_iter {:.2}MB", current_mem);
    fasta_reader
        .records()
        .collect::<Vec<_>>()
        .into_par_iter()
        .for_each(|result| {
            if let Ok(result_data) = result {
                process_record(result_data.id());
            } else {
                println!("Error reading record");
            }
        });
    let duration = start.elapsed();
    println!("Time rust-bio into_par_iter: {:.2}", duration.as_secs_f32());
    // drop(fasta_reader);

    let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    println!("The max amount that was used {}GB", peak_mem);

    Ok(())
}
