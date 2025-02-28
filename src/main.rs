#![allow(unused)]

mod check_processing_metrics;

use std::{ fs::File, io::Write };
use bio::io::fasta::{ Record, Reader };
use rayon::prelude::*;
use check_processing_metrics::check_processing_metrics;

pub struct SpliceConfig {
    pub subtype: String,
    pub strain: String,
    pub splice_form_name: String,
    pub donor: String,
    pub donor_sequence: String,
    pub receptor: String,
    pub receptor_sequence: String,
}

pub struct SpliceResults {
    pub splice_form: String,
    pub id_sequence: String,
}

fn main() {
    check_processing_metrics(
        &format!("/app/input/Simonetti/A_468d23_CAGATCA_S1_L001_R1_001.fastq")
    ).unwrap();

    return;

    // set thread count if a maximum is needed
    // rayon::ThreadPoolBuilder::new().num_threads(cores).build_global().unwrap();

    let user_input = SpliceConfig {
        subtype: "B".to_string(),
        strain: "HXB2".to_string(),
        splice_form_name: "A1".to_string(),
        donor: "5'SS".to_string(),
        receptor: "3'SS".to_string(),
        donor_sequence: "GTAGTGTG".to_string(),
        receptor_sequence: "AGTGTG".to_string(),
    };

    let filename: &str = "/app/src/seqs.fasta";

    // TODO: decompress file if not

    let reader = Reader::from_file(filename).unwrap();

    let mut results: Vec<SpliceResults> = vec![];

    // let mut buffer = Vec::new();
    // buffer.read_to_end(&mut buffer).unwrap();

    // let f = File::open(filename).expect("Unable to open file");
    // let f = BufReader::new(f);

    // for line in f.lines() {
    //     let line = line.expect("Unable to read line");
    //     println!("Line: {}", line);
    // }

    // let chunk_size = 0x1000 /* or anything proper (e.g. by row) */;
    // for nth_chunk in 0..=buffer.len() / chunk_size {
    //     let start = nth_chunk * chunk_size;
    //     let chunk = &mut buffer[start..(start + chunk_size).min(buffer.len())];
    //     // parallel_operation(chunk);
    //     // output(&buffer);
    // }

    // reader
    //     .records()
    //     .into_par_iter()
    //     .for_each(|result| {
    //         if result.is_err() {
    //             println!("Error reading record");
    //             return;
    //         }
    //         let record = result.unwrap();
    //         let splice_results = identify_hiv_spliceform(&user_input, record);
    //         results.push(splice_results);
    //     });

    for result in reader.records() {
        if result.is_err() {
            println!("Error reading record");
            continue;
        }
        let record = result.unwrap();
        let splice_results = identify_hiv_spliceform(&user_input, record);
        results.push(splice_results);
    }

    // convert results to csv file
    let mut writer = File::create("/app/src/results.csv").expect("");
    writer.write(b"splice_form,id_sequence\n").expect("");
    for result in results {
        writer
            .write(format!("{},{}\n", result.splice_form, result.id_sequence).as_bytes())
            .expect("");
    }
}

fn identify_hiv_spliceform(user_input: &SpliceConfig, record: Record) -> SpliceResults {
    let splice_form = user_input.splice_form_name.clone();
    let id_sequence = user_input.donor_sequence.clone();

    println!("{}: {}", splice_form, id_sequence);

    SpliceResults {
        splice_form,
        id_sequence,
    }
}
