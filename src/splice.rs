#![allow(dead_code, unused_variables)]

use bio::alignment::Alignment;
use bio::alphabets;
use bio::data_structures::bwt::{ bwt, less, Occ };
use bio::data_structures::fmindex::{ BackwardSearchResult, FMIndex, FMIndexable };
use bio::data_structures::suffix_array::suffix_array;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use std::collections::HashMap;
use std::io;

fn main_vars() {
    //ref is reference sequence from post D1 cut to end of transcript
    //this is the NL4-3 sequence
    let reference_sequence: &[
        u8;
        8883
    ] = b"GTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCGGTATTAAGCGGGGGAGAATTAGATAAATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAACAATATAAACTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTTTTAGAGACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAATAGCAGTCCTCTATTGTGTGCATCAAAGGATAGATGTAAAAGACACCAAGGAAGCCTTAGATAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAGGCACAGCAAGCAGCAGCTGACACAGGAAACAACAGCCAGGTCAGCCAAAATTACCCTATAGTGCAGAACCTCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTAATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAATACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGATTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACACATAATCCACCTATCCCAGTAGGAGAAATCTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGATTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAAGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGGAGCGACACTAGAAGAAATGATGACAGCATGTCAGGGAGTGGGGGGACCCGGCCATAAAGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATCCAGCTACCATAATGATACAGAAAGGCAATTTTAGGAACCAAAGAAAGACTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACATAGCCAAAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCCACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTTTGGGGAAGAGACAACAACTCCCTCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAGCTTCCCTCAGATCACTCTTTGGCAGCGACCCCTCGTCACAATAAAGATAGGGGGGCAATTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGCGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGCTGCACTTTAAATTTTCCCATTAGTCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGAAGGAAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGATTTCTGGGAAGTTCAATTAGGAATACCACATCCTGCAGGGTTAAAACAGAAAAAATCAGTAACAGTACTGGATGTGGGCGATGCATATTTTTCAGTTCCCTTAGATAAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAGTGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTCATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAACTGAGACAACATCTGTTGAGGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAGGACAGCTGGACTGTCAATGACATACAGAAATTAGTGGGAAAATTGAATTGGGCAAGTCAGATTTATGCAGGGATTAAAGTAAGGCAATTATGTAAACTTCTTAGGGGAACCAAAGCACTAACAGAAGTAGTACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGGGAGATTCTAAAAGAACCGGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAAGGGTGCCCACACTAATGATGTGAAACAATTAACAGAGGCAGTACAAAAAATAGCCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAATTACCCATACAAAAGGAAACATGGGAAGCATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTCAATACCCCTCCCTTAGTGAAGTTATGGTACCAGTTAGAGAAAGAACCCATAATAGGAGCAGAAACTTTCTATGTAGATGGGGCAGCCAATAGGGAAACTAAATTAGGAAAAGCAGGATATGTAACTGACAGAGGAAGACAAAAAGTTGTCCCCCTAACGGACACAACAAATCAGAAGACTGAGTTACAAGCAATTCATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTGACAGACTCACAATATGCATTGGGAATCATTCAAGCACAACCAGATAAGAGTGAATCAGAGTTAGTCAGTCAAATAATAGAGCAGTTAATAAAAAAGGAAAAAGTCTACCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATGGGTTGGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGAAGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTACCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGGGAAGCCATGCATGGACAAGTAGACTGTAGCCCAGGAATATGGCAGCTAGATTGTACACATTTAGAAGGAAAAGTTATCTTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTAATTCCAGCAGAGACAGGGCAAGAAACAGCATACTTCCTCTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAGTACATACAGACAATGGCAGCAATTTCACCAGTACTACAGTTAAGGCCGCCTGTTGGTGGGCGGGGATCAAGCAGGAATTTGGCATTCCCTACAATCCCCAAAGTCAAGGAGTAATAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAGATCCAGTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATCAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAACACATGGAAAAGATTAGTAAAACACCATATGTATATTTCAAGGAAAGCTAAGGACTGGTTTTATAGACATCACTATGAAAGTACTAATCCAAAAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAAATTAGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGACCTAGCAGACCAACTAATTCATCTGCACTATTTTGATTGTTTTTCAGAATCTGCTATAAGAAATACCATATTAGGACGTATAGTTAGTCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAGTACTTGGCACTAGCAGCATTAATAAAACCAAAACAGATAAAGCCACCTTTGCCTAGTGTTAGGAAACTGACAGAGGACAGATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCATACAATGAATGGACACTAGAGCTTTTAGAGGAACTTAAGAGTGAAGCTGTTAGACATTTTCCTAGGATATGGCTCCATAACTTAGGACAACATATCTATGAAACTTACGGGGATACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATCCATTTCAGAATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCAAGTTTGTTTCATGACAAAAGCCTTAGGCATCTCCTATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCAGTAAGTAGTACATGTAATGCAACCTATAATAGTAGCAATAGTAGCATTAGTAGTAGCAATAATAATAGCAATAGTTGTGTGGTCCATAGTAATCATAGAATATAGGAAAATATTAAGACAAAGAAAAATAGACAGGTTAATTGATAGACTAATAGAAAGAGCAGAAGACAGTGGCAATGAGAGTGAAGGAGAAGTATCAGCACTTGTGGAGATGGGGGTGGAAATGGGGCACCATGCTCCTTGGGATATTGATGATCTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTTAACATGTGGAAAAATGACATGGTAGAACAGATGCATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGAAGAATGATACTAATACCAATAGTAGTAGCGGGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAGCATAAGAGATAAGGTGCAGAAAGAATATGCATTCTTTTATAAACTTGATATAGTACCAATAGATAATACCAGCTATAGGTTGATAAGTTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATCAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGATGTAGTAATTAGATCTGCCAATTTCACAGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAATTAATTGTACAAGACCCAACAACAATACAAGAAAAAGTATCCGTATCCAGAGGGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATGCCACTTTAAAACAGATAGCTAGCAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACACTCCCATGCAGAATAAAACAATTTATAAACATGTGGCAGGAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAATATTACTGGGCTGCTATTAACAAGAGATGGTGGTAATAACAACAATGGGTCCGAGATCTTCAGACCTGGAGGAGGCGATATGAGGGACAATTGGAGAAGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCTGCACGTCAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGATATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAACAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAGATTTGGAATAACATGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAATCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTAGCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAACTTGCTCAATGCCACAGCCATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTATTACAAGCAGCTTATAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAAGATGGGTGGCAAGTGGTCAAAAAGTAGTGTGATTGGATGGCCTGCTGTAAGGGAAAGAATGAGACGAGCTGAGCCAGCAGCAGATGGGGTGGGAGCAGTATCTCGAGACCTAGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTAACAATGCTGCTTGTGCCTGGCTAGAAGCACAAGAGGAGGAAGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAAAGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTGGCAGAACTACACACCAGGGCCAGGGGTCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGATAAGGTAGAAGAGGCCAATAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCTGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACGTGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATGCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCA";

    //this is the RC of $ref, from the end of the transcript to D1
    let rc: &[
        u8;
        8883
    ] = b"TGAAGCACTCAAGGCAAGCTTTATTGAGGCTTAAGCAGTGGGTTCCCTAGTTAGCCAGAGAGCTCCCAGGCTCAGATCTGGTCTAACCAGAGAGACCCAGTACAGGCAAAAAGCAGCTGCTTATATGCAGCATCTGAGGGCTCGCCACTCCCCAGTCCCGCCCAGGCCACGCCTCCCTGGAAAGTCCCCAGCGGAAAGTCCCTTGTAGCAAGCTCGATGTCAGCAGTTCTTGAAGTACTCCGGATGCAGCTCTCGGGCCACGTGATGAAATGCTAGGCGGCTGTCAAACCTCCACTCTAACACTTCTCTCTCAGGGTCATCCATTCCATGCAGGCTCACAGGGTGTAACAAGCTGGTGTTCTCTCCTTTATTGGCCTCTTCTACCTTATCTGGCTCAACTGGTACTAGCTTGTAGCACCATCCAAAGGTCAGTGGATATCTGACCCCTGGCCCTGGTGTGTAGTTCTGCCAATCAGGGAAGTAGCCTTGTGTGTGGTAGATCCACAGATCAAGGATATCTTGTCTTCTTTGGGAGTGAATTAGCCCTTCCAGTCCCCCCTTTTCTTTTAAAAAGTGGCTAAGATCTACAGCTGCCTTGTAAGTCATTGGTCTTAAAGGTACCTGAGGTGTGACTGGAAAACCCACCTCTTCCTCCTCTTGTGCTTCTAGCCAGGCACAAGCAGCATTGTTAGCTGCTGTATTGCTACTTGTGATTGCTCCATGTTTTTCTAGGTCTCGAGATACTGCTCCCACCCCATCTGCTGCTGGCTCAGCTCGTCTCATTCTTTCCCTTACAGCAGGCCATCCAATCACACTACTTTTTGACCACTTGCCACCCATCTTATAGCAAAATCCTTTCCAAGCCCTGTCTTATTCTTCTAGGTATGTGGCGAATAGCTCTATAAGCTGCTTGTAATACTTCTATAACCCTATCTGTCCCCTCAGCTACTGCTATGGCTGTGGCATTGAGCAAGTTAACAGCACTATTCTTTAGTTCCTGACTCCAATACTGTAGGAGATTCCACCAATATTTGAGGGCTTCCCACCCCCTGCGTCCCAGAAGTTCCACAATCCTCGTTACAATCAAGAGTAAGTCTCTCAAGCGGTGGTAGCTGAAGAGGCACAGGCTCCGCAGATCGTCCCAGATAAGTGCTAAGGATCCGTTCACTAATCGAATGGATCTGTCTCTGTCTCTCTCTCCACCTTCTTCTTCTATTCCTTCGGGCCTGTCGGGTCCCCTCGGGATTGGGAGGTGGGTCTGAAACGATAATGGTGAATATCCCTGCCTAACTCTATTCACTATAGAAAGTACAGCAAAAACTATTCTTAAACCTACCAAGCCTCCTACTATCATTATGAATAATTTTATATACCACAGCCAATTTGTTATGTTAAACCAATTCCACAAACTTGCCCATTTATCTAATTCCAATAATTCTTGTTCATTCTTTTCTTGCTGGTTTTGCGATTCTTCAATTAAGGAGTGTATTAAGCTTGTGTAATTGTTAATTTCTCTGTCCCACTCCATCCAGGTCATGTTATTCCAAATCTGTTCCAGAGATTTATTACTCCAACTAGCATTCCAAGGCACAGCAGTGGTGCAAATGAGTTTTCCAGAGCAACCCCAAATCCCCAGGAGCTGTTGATCCTTTAGGTATCTTTCCACAGCCAGGATTCTTGCCTGGAGCTGTTTGATGCCCCAGACTGTGAGTTGCAACAGATGCTGTTGCGCCTCAATAGCCCTCAGCAAATTGTTCTGCTGCTGCACTATATCAGACAATAATTGTCTGGCCTGTACCGTCAGCGTCATTGACGTGCAGCCCATAGTGCTTCCTGCTGCTCCCAAGAACCCAAGGAACAAAGCTCCTATTCCCACTGCTCTTTTTTCTCTCTGCACCACTCTTCTCTTTGCCTTGGTGGGTGCTACTCCTAATGGTTCAATTTTTACTACTTTATATTTATATAATTCACTTCTCCAATTGTCCCTCATATCGCCTCCTCCAGGTCTGAAGATCTCGGACCCATTGTTGTTATTACCACCATCTCTTGTTAATAGCAGCCCAGTAATATTTGATGAACATCTAATTTGTCCACTGATGGGAGGGGCATACATTGCTTTTCCTACTTCCTGCCACATGTTTATAAATTGTTTTATTCTGCATGGGAGTGTGATTGTGTCACTTCCTTCAGTGTTATTTGACCCTTCAGTACTCCAAGTACTATTAAACCAAGTACTATTAAACAGTTGTGTTGAATTACAGTAGAAAAATTCCCCTCCACAATTAAAACTGTGCGTTACAATTTCTGGGTCCCCTCCTGAGGATTGCTTAAAGATTATTGTTTTATTATTTCCAAATTGTTCTCTTAATTTGCTAGCTATCTGTTTTAAAGTGGCATTCCATTTTGCTCTACTAATGTTACAATGTGCTTGTCTCATATTTCCTATTTTTCCTATTGTAACAAATGCTCTCCCTGGTCCCCTCTGGATACGGATACTTTTTCTTGTATTGTTGTTGGGTCTTGTACAATTAATTTCTACAGATGTGTTCAGCTGTACTATTATGGTTTTAGCATTGTCTGTGAAATTGGCAGATCTAATTACTACATCTTCTTCTGCTAGACTGCCATTTAACAGCAGTTGAGTTGATACTACTGGCCTGATTCCATGTGTACATTGTACTGTGCTGACATTTGTACATGGTCCTGTTCCATTGAACGTCTTATTATTACATTTTAGAATCGCAAAACCAGCCGGGGCACAATAATGTATGGGAATTGGCTCAAAGGATACCTTTGGACAGGCCTGTGTAATGACTGAGGTGTTACAACTTATCAACCTATAGCTGGTATTATCTATTGGTACTATATCAAGTTTATAAAAGAATGCATATTCTTTCTGCACCTTATCTCTTATGCTTGTGCTGATATTGAAAGAGCAGTTTTTTATCTCTCCTTTCTCCATTATCATTCTCCCGCTACTACTATTGGTATTAGTATCATTCTTCAAATCAGTGCACTTTAAACTAACACAGAGTGGGGTTAATTTTACACATGGCTTTAGGCTTTGATCCCATAAACTGATTATATCCTCATGCATCTGTTCTACCATGTCATTTTTCCACATGTTAAAATTTTCTGTCACATTTACCAATACTACTTCTTGTGGGTTGGGGTCTGTGGGTACACAGGCATGTGTGGCCCAAACATTATGTACCTCTGTATCATATGCTTTAGCATCTGATGCACAAAATAGAGTGGTGGTTGCTTCCTTCCACACAGGTACCCCATAATAGACTGTGACCCACAATTTTTCTGTAGCACTACAGATCATCAATATCCCAAGGAGCATGGTGCCCCATTTCCACCCCCATCTCCACAAGTGCTGATACTTCTCCTTCACTCTCATTGCCACTGTCTTCTGCTCTTTCTATTAGTCTATCAATTAACCTGTCTATTTTTCTTTGTCTTAATATTTTCCTATATTCTATGATTACTATGGACCACACAACTATTGCTATTATTATTGCTACTACTAATGCTACTATTGCTACTATTATAGGTTGCATTACATGTACTACTTACTGCTTTGATAGAGAAGCTTGATGAGTCTGACTGTTCTGATGAGCTCTTCGTCGCTGTCTCCGCTTCTTCCTGCCATAGGAGATGCCTAAGGCTTTTGTCATGAAACAAACTTGGCAATGAAAGCAACACTTTTTACAATAGCAATTGGTACAAGCAGTTTTAGGCTGACTTCCTGGATGCTTCCAGGGCTCTAGTCTAGGATCTACTGGCTCCATTTCTTGCTCTCCTCTGTCGAGTAACGCCTATTCTGCTATGTCGACACCCAATTCTGAAATGGATAAACAGCAGTTGTTGCAGAATTCTTATTATGGCTTCCACTCCTGCCCAAGTATCCCCGTAAGTTTCATAGATATGTTGTCCTAAGTTATGGAGCCATATCCTAGGAAAATGTCTAACAGCTTCACTCTTAAGTTCCTCTAAAAGCTCTAGTGTCCATTCATTGTATGGCTCCCTCTGTGGCCCTTGGTCTTCTGGGGCTTGTTCCATCTGTCCTCTGTCAGTTTCCTAACACTAGGCAAAGGTGGCTTTATCTGTTTTGGTTTTATTAATGCTGCTAGTGCCAAGTACTGTAGAGATCCTACCTTGTTATGTCCTGCTTGATATTCACACCTAGGACTAACTATACGTCCTAATATGGTATTTCTTATAGCAGATTCTGAAAAACAATCAAAATAGTGCAGATGAATTAGTTGGTCTGCTAGGTCAGGGTCTACTTGTGTGCTATATCTCTTTTTCCTCCATTCTATGGAGACTCCCTGACCCAAATGCCAGTCTCTTTCTCCTGTATGCAGACCCCAATATGTTGTTATTACTAATTTAGCATCCCCTAGTGGGATGTGTACTTCTGAACTTATTTTTGGATTAGTACTTTCATAGTGATGTCTATAAAACCAGTCCTTAGCTTTCCTTGAAATATACATATGGTGTTTTACTAATCTTTTCCATGTGTTAATCCTCATCCTGTCTACTTGCCACACAATCATCACCTGCCATCTGTTTTCCATAATCCCTGATGATCTTTGCTTTTCTTCTTGGCACTACTTTTATGTCACTATTATCTTGTATTACTACTGCCCCTTCACCTTTCCAGAGGAGCTTTGCTGGTCCTTTCCAAACTGGATCTCTGCTGTCCCTGTAATAAACCCGAAAATTTTGAATTTTTGTAATTTGTTTTTGTAATTCTTTAGTTTGTATGTCTGTTGCTATTATGTCTACTATTCTTTCCCCTGCACTGTACCCCCCAATCCCCCCTTTTCTTTTAAAATTGTGGATGAATACTGCCATTTGTACTGCTGTCTTAAGATGTTCAGCCTGATCTCTTACCTGTCCTATAATTTTCTTTAATTCTTTATTCATAGATTCTATTACTCCTTGACTTTGGGGATTGTAGGGAATGCCAAATTCCTGCTTGATCCCCGCCCACCAACAGGCGGCCTTAACTGTAGTACTGGTGAAATTGCTGCCATTGTCTGTATGTACTGTTTTTACTGGCCATCTTCCTGCTAATTTTAAGAGGAAGTATGCTGTTTCTTGCCCTGTCTCTGCTGGAATTACTTCTGCTTCTATATATCCACTGGCTACATGAACTGCTACCAAGATAACTTTTCCTTCTAAATGTGTACAATCTAGCTGCCATATTCCTGGGCTACAGTCTACTTGTCCATGCATGGCTTCCCCTTTTAGCTGACATTTATCACAGCTGGCTACTATTTCTTTTGCTACTACAGGTGGTAGGTTAAAATCACTAGCCATTGCTCTCCAATTACTGTGATATTTCTCATGTTCTTCTTGGGCCTTATCTATTCCATCTAAAAATAGTACTTTCCTGATTCCAGCACTGACCAACCCATCTACTTGTTCATTTCCTCCAATTCCTTTGTGTGCTGGTACCCATGCCAGGTAGACTTTTTCCTTTTTTATTAACTGCTCTATTATTTGACTGACTAACTCTGATTCACTCTTATCTGGTTGTGCTTGAATGATTCCCAATGCATATTGTGAGTCTGTCACTATGTTTACTTCTAATCCCGAATCCTGCAAAGCTAGATGAATTGCTTGTAACTCAGTCTTCTGATTTGTTGTGTCCGTTAGGGGGACAACTTTTTGTCTTCCTCTGTCAGTTACATATCCTGCTTTTCCTAATTTAGTTTCCCTATTGGCTGCCCCATCTACATAGAAAGTTTCTGCTCCTATTATGGGTTCTTTCTCTAACTGGTACCATAACTTCACTAAGGGAGGGGTATTGACAAACTCCCACTCAGGAATCCAGGTGGCTTGCCAATACTCTGTCCACCATGCTTCCCATGTTTCCTTTTGTATGGGTAATTTAAATTTAGGAGTCTTTCCCCATATTACTATGCTTTCTGTGGCTATTTTTTGTACTGCCTCTGTTAATTGTTTCACATCATTAGTGTGGGCACCCTTCATTCTTGCATATTTTCCTGTTTTCAGATTTTTAAATGGCTCTTGATAAATTTGATATGTCCATTGGCCTTGCCCCTGCTTCTGTATTTCTGCTATTAAGTCTTTTGATGGGTCATAATACACTCCATGTACCGGTTCTTTTAGAATCTCCCTGTTTTCTGCCAGTTCTAGCTCTGCTTCTTCTGTTAGTGGTACTACTTCTGTTAGTGCTTTGGTTCCCCTAAGAAGTTTACATAATTGCCTTACTTTAATCCCTGCATAAATCTGACTTGCCCAATTCAATTTTCCCACTAATTTCTGTATGTCATTGACAGTCCAGCTGTCCTTTTCTGGCAGCACTATAGGCTGTACTGTCCATTTATCAGGATGGAGTTCATAACCCATCCAAAGGAATGGAGGTTCTTTCTGATGTTTTTTGTCTGGTGTGGTAAATCCCCACCTCAACAGATGTTGTCTCAGTTCCTCTATTTTTGTTCTATGCTGCCCTATTTCTAAGTCAGATCCTACATACAAATCATCCATGTATTGATAGATGACTATGTCTGGATTTTGTTTTCTAAAAGGCTCTAAGATTTTTGTCATGCTACACTGGAATATTGCTGGTGATCCTTTCCATCCCTGTGGAAGCACATTGTACTGATATCTAATCCCTGGTGTCTCATTGTTTATACTAGGTATGGTAAATGCAGTATACTTCCTGAAGTCTTTATCTAAGGGAACTGAAAAATATGCATCGCCCACATCCAGTACTGTTACTGATTTTTTCTGTTTTAACCCTGCAGGATGTGGTATTCCTAATTGAACTTCCCAGAAATCTTGAGTTCTCTTATTAAGTTCTCTGAAATCTACTAATTTTCTCCATTTAGTACTGTCTTTTTTCTTTATGGCAAATACTGGAGTATTGTATGGATTTTCAGGCCCAATTTTTGAAATTTTTCCTTCCTTTTCCATTTCTGTACAAATTTCTACTAATGCTTTTATTTTTTCTTCTGTCAATGGCCATTGTTTAACTTTTGGGCCATCCATTCCTGGCTTTAATTTTACTGGTACAGTCTCAATAGGACTAATGGGAAAATTTAAAGTGCAGCCAATCTGAGTCAACAGATTTCTTCCAATTATGTTGACAGGTGTAGGTCCTACTAATACTGTACCTATAGCTTTATGTCCGCAGATTTCTATGAGTATCTGATCATACTGTCTTACTTTGATAAAACCTCCAATTCCCCCTATCATTTTTGGTTTCCATCTTCCTGGCAAATTCATTTCTTCTAATACTGTATCATCTGCTCCTGTATCTAATAGAGCTTCCTTTAATTGCCCCCCTATCTTTATTGTGACGAGGGGTCGCTGCCAAAGAGTGATCTGAGGGAAGCTAAAGGATACAGTTCCTTGTCTATCGGCTCCTGCTTCTGAGAGGGAGTTGTTGTCTCTTCCCCAAACCTGAAGCTCTCTTCTGGTGGGGCTGTTGGCTCTGGTCTGCTCTGAAGAAAATTCCCTGGCCTTCCCTTGTGGGAAGGCCAGATCTTCCCTAAAAAATTAGCCTGTCTCTCAGTACAATCTTTCATTTGGTGTCCTTCCTTTCCACATTTCCAACAGCCCTTTTTCCTAGGGGCCCTGCAATTTTTGGCTATGTGCCCTTCTTTGCCACAATTGAAACACTTAACAGTCTTTCTTTGGTTCCTAAAATTGCCTTTCTGTATCATTATGGTAGCTGGATTTGTTACTTGGCTCATTGCTTCAGCCAAAACTCTTGCTTTATGGCCGGGTCCCCCCACTCCCTGACATGCTGTCATCATTTCTTCTAGTGTCGCTCCTGGTCCCAATGCTTTTAAAATAGTCTTACAATCTGGGTTCGCATTTTGGACCAACAAGGTTTCTGTCATCCAATTTTTTACCTCTTGTGAAGCTTGCTCGGCTCTTAGAGTTTTATAGAATCGGTCTACATAGTCTCTAAAGGGTTCCTTTGGTCCTTGTCTTATGTCCAGAATGCTGGTAGGGCTATACATTCTTACTATTTTATTTAATCCCAGGATTATCCATCTTTTATAGATTTCTCCTACTGGGATAGGTGGATTATGTGTCATCCATCCTATTTGTTCCTGAAGGGTACTAGTAGTTCCTGCTATGTCACTTCCCCTTGGTTCTCTCATCTGGCCTGGTGCAATAGGCCCTGCATGCACTGGATGCAATCTATCCCATTCTGCAGCTTCCTCATTGATGGTCTCTTTTAACATTTGCATGGCTGCTTGATGTCCCCCCACTGTGTTTAGCATGGTATTTAAATCTTGTGGGGTGGCTCCTTCTGATAATGCTGAAAACATGGGTATTACTTCTGGGCTGAAAGCCTTCTCTTCTACTACTTTTACCCATGCATTTAAAGTTCTAGGTGATATGGCCTGATGTACCATTTGCCCCTGGAGGTTCTGCACTATAGGGTAATTTTGGCTGACCTGGCTGTTGTTTCCTGTGTCAGCTGCTGCTTGCTGTGCCTTTTTCTTACTTTTGTTTTGCTCTTCCTCTATCTTATCTAAGGCTTCCTTGGTGTCTTTTACATCTATCCTTTGATGCACACAATAGAGGACTGCTATTGTATTATATAATGATCTAAGTTCTTCTGATCCTGTCTGAAGGGATGGTTGTAGCTGTCCCAGTATTTGTCTACAGCCTTCTGATGTCTCTAAAAGGCCAGGATTAACTGCGAATCGTTCTAGCTCCCTGCTTGCCCATACTATATGTTTTAGTTTATATTGTTTCTTTCCCCCTGGCCTTAACCGAATTTTTTCCCATTTATCTAATTCTCCCCCGCTTAATACCGACGCTCTCGCACCCATCTCTCTCCTTCTAGCCTCCGCTAGTCAAAATTTTTGGCGTACTCAC";

    //this is the RC of the 600 bases after D1
    let rc_ref: &[
        u8;
        600
    ] = b"TTTAAATCTTGTGGGGTGGCTCCTTCTGATAATGCTGAAAACATGGGTATTACTTCTGGGCTGAAAGCCTTCTCTTCTACTACTTTTACCCATGCATTTAAAGTTCTAGGTGATATGGCCTGATGTACCATTTGCCCCTGGAGGTTCTGCACTATAGGGTAATTTTGGCTGACCTGGCTGTTGTTTCCTGTGTCAGCTGCTGCTTGCTGTGCCTTTTTCTTACTTTTGTTTTGCTCTTCCTCTATCTTATCTAAGGCTTCCTTGGTGTCTTTTACATCTATCCTTTGATGCACACAATAGAGGACTGCTATTGTATTATATAATGATCTAAGTTCTTCTGATCCTGTCTGAAGGGATGGTTGTAGCTGTCCCAGTATTTGTCTACAGCCTTCTGATGTCTCTAAAAGGCCAGGATTAACTGCGAATCGTTCTAGCTCCCTGCTTGCCCATACTATATGTTTTAGTTTATATTGTTTCTTTCCCCCTGGCCTTAACCGAATTTTTTCCCATTTATCTAATTCTCCCCCGCTTAATACCGACGCTCTCGCACCCATCTCTCTCCTTCTAGCCTCCGCTAGTCAAAATTTTTGGCGTACTCAC";

    //this is the RC of the 600 bases after A7
    let rc_complete: &[
        u8;
        600
    ] = b"TGTGCTTCTAGCCAGGCACAAGCAGCATTGTTAGCTGCTGTATTGCTACTTGTGATTGCTCCATGTTTTTCTAGGTCTCGAGATACTGCTCCCACCCCATCTGCTGCTGGCTCAGCTCGTCTCATTCTTTCCCTTACAGCAGGCCATCCAATCACACTACTTTTTGACCACTTGCCACCCATCTTATAGCAAAATCCTTTCCAAGCCCTGTCTTATTCTTCTAGGTATGTGGCGAATAGCTCTATAAGCTGCTTGTAATACTTCTATAACCCTATCTGTCCCCTCAGCTACTGCTATGGCTGTGGCATTGAGCAAGTTAACAGCACTATTCTTTAGTTCCTGACTCCAATACTGTAGGAGATTCCACCAATATTTGAGGGCTTCCCACCCCCTGCGTCCCAGAAGTTCCACAATCCTCGTTACAATCAAGAGTAAGTCTCTCAAGCGGTGGTAGCTGAAGAGGCACAGGCTCCGCAGATCGTCCCAGATAAGTGCTAAGGATCCGTTCACTAATCGAATGGATCTGTCTCTGTCTCTCTCTCCACCTTCTTCTTCTATTCCTTCGGGCCTGTCGGGTCCCCTCGGGATTGGGAGGTGGGT";

    //this is the RC of the 600 bases after D4
    let rc_incomp: &[
        u8;
        600
    ] = b"TATTGGTATTAGTATCATTCTTCAAATCAGTGCACTTTAAACTAACACAGAGTGGGGTTAATTTTACACATGGCTTTAGGCTTTGATCCCATAAACTGATTATATCCTCATGCATCTGTTCTACCATGTCATTTTTCCACATGTTAAAATTTTCTGTCACATTTACCAATACTACTTCTTGTGGGTTGGGGTCTGTGGGTACACAGGCATGTGTGGCCCAAACATTATGTACCTCTGTATCATATGCTTTAGCATCTGATGCACAAAATAGAGTGGTGGTTGCTTCCTTCCACACAGGTACCCCATAATAGACTGTGACCCACAATTTTTCTGTAGCACTACAGATCATCAATATCCCAAGGAGCATGGTGCCCCATTTCCACCCCCATCTCCACAAGTGCTGATACTTCTCCTTCACTCTCATTGCCACTGTCTTCTGCTCTTTCTATTAGTCTATCAATTAACCTGTCTATTTTTCTTTGTCTTAATATTTTCCTATATTCTATGATTACTATGGACCACACAACTATTGCTATTATTATTGCTACTACTAATGCTACTATTGCTACTATTATAGGTTGCATTACATGTACTACTTAC";
}

fn splice() {
    // a given text
    let text = b"ACAGCTCGATCGGTA$";
    let pattern = b"ATCG";

    // Create an FM-Index for the given text.

    // instantiate an alphabet
    let alphabet = alphabets::dna::iupac_alphabet();
    // calculate a suffix array
    let sa = suffix_array(text);
    // calculate the Burrows-Wheeler-transform
    let bwt = bwt(text, &sa);
    // calculate the vectors less and Occ (occurrences)
    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, 3, &alphabet);
    // set up FMIndex
    let fmindex = FMIndex::new(&bwt, &less, &occ);
    // do a backwards search for the pattern
    let interval = fmindex.backward_search(pattern.iter());
    let mut partial_match_len = 0;
    // get the locations where the pattern matched (completely in this case).
    let positions = match interval {
        BackwardSearchResult::Complete(saint) => saint.occ(&sa),
        BackwardSearchResult::Partial(saint, l) => {
            partial_match_len = l;
            saint.occ(&sa)
        }
        BackwardSearchResult::Absent => Vec::new(),
    };
    // Iterate over a FASTQ file, use the alphabet to validate read
    // sequences and search for exact matches in the FM-Index.

    // create FASTQ reader
    let mut reader = fastq::Reader::new(io::stdin());
    let mut record = fastq::Record::new();
    let mut partial_match_len = 0;
    reader.read(&mut record).expect("Failed to parse record");
    while !record.is_empty() {
        let check = record.check();
        if check.is_err() {
            panic!("I got a rubbish record!");
        }
        // obtain sequence
        let seq = record.seq();
        // check, whether seq is in the expected alphabet

        if alphabet.is_word(seq) {
            let interval = fmindex.backward_search(seq.iter());
            // get the positions where seq matched completely
            // or where the maximal matching suffix of seq occurred.
            let positions = match interval {
                BackwardSearchResult::Complete(saint) => saint.occ(&sa),
                BackwardSearchResult::Partial(saint, l) => {
                    partial_match_len = l;
                    saint.occ(&sa)
                }
                BackwardSearchResult::Absent => Vec::new(),
            };
        }
        reader.read(&mut record).expect("Failed to parse record");
    }

    println!("Hello, world!");
}

// returns the reverse complement of a DNA sequence
fn return_rc(dna_sequence: &str) -> String {
    // Reverse the sequence
    let rev = dna_sequence.chars().rev();
    let mut reverse_complement_strand = String::new();
    let complements = ['T', 'G', 'C', 'A'];
    let bases = ['A', 'C', 'G', 'T'];

    // Complement each base ( A <=> T and C <=> G )
    for ch in rev {
        if let Some(pos) = bases.iter().position(|&b| b == ch.to_ascii_uppercase()) {
            reverse_complement_strand.push(complements[pos]);
        }
    }

    return reverse_complement_strand;
}

// finds the index of an unknown acceptor sequence, adapted to be the FL genome index
fn identify_acceptor(unknown_acceptor_sequence: &str, ref_sequence: &str) -> Option<usize> {
    // the unknownAcceptor Sequence is 10 bases long
    // ref is reference sequence from post D1 cut to end of transcript
    // this is the NL4-3 sequence
    let search_index = ref_sequence.find(unknown_acceptor_sequence);
    if let Some(idx) = search_index {
        let actual_acceptor_index = idx + 743 + 1; // index of post D1, +1 because Ruby indexes at 0
        return Some(actual_acceptor_index);
    }
    return None;
}

// returns t/f if edit distance is less than maxAllowedEditDistance
fn close_enough(
    sequence: Option<&[u8]>,
    reference_sequence: &[u8],
    max_allowed_edit_distance: usize
) -> bool {
    let mut score_table = vec![];

    if sequence.is_none() {
        return false;
    }

    let seq = sequence.unwrap();

    let reference_length = reference_sequence.len();
    let sequence_length = seq.len();

    for xidx in 0..=sequence_length {
        score_table.push(vec![]);
    }

    // reference is across the top
    // test sequence is on the left column
    // initialize score table
    // initialize top row
    for xidx in 0..=reference_length {
        score_table[0].push(xidx);
    }

    // initialize first column
    for xidx in 1..=sequence_length {
        score_table[xidx].push(xidx);
    }

    // TODO: replace this with bio::alignment
    // https://docs.rs/bio/latest/bio/alignment/struct.Alignment.html

    // let alignment: Alignment = Alignment::new(sequence, reference_sequence, 1, 1, 1, false);
    // let edit_distance = alignment.edit_distance();
    // return edit_distance <= max_allowed_edit_distance;

    // and now we come to the actual alignment
    for sequence_idx in 0..sequence_length {
        for reference_idx in 0..reference_length {
            // score table reference is across the top, test sequence is down the left side

            let diag = score_table[sequence_idx][reference_idx];
            let up = score_table[sequence_idx][reference_idx + 1];
            let left = score_table[sequence_idx + 1][reference_idx];

            let row_array = [diag, left + 1, up + 1];

            if reference_sequence[reference_idx] == seq[sequence_idx] {
                // it's a perfect match
                // append the values to the reference index row array
                let min_score = row_array.iter().min().unwrap();
                score_table[sequence_idx + 1].push(*min_score);
            } else {
                // they don't match, cost to mutate = 1, cost to insert or delete = 1
                let min_score = row_array.iter().min().unwrap();
                score_table[sequence_idx + 1].push(*min_score);
            }
        }
    }

    let edit_distance = score_table[sequence_length][reference_length];

    return edit_distance <= max_allowed_edit_distance;
}

fn sort_forward_sequence(sequence: &str, unknown_acceptors: &mut Vec<[&str; 2]>) -> Vec<&str> {
    let mut splice_map = vec![];
    let l = sequence.len();

    // look for search sequences of forward primer idx_fp
    if sequence[4..8].contains("TGCTGA") {
        let fp_idx = sequence[4..8].find("TGCTGA");
        let fp_index = fp_idx.unwrap() + 4; //start of fp
    } else if close_enough(Some(&sequence[4..9]), "TGCTGAAGC", 1) {
        let fp_index = 4;
    } else {
        // not finding forward primer
        return splice_map;
    }

    // remove forward primer dimers
    if sequence[10..50].contains("CTGAACT") {
        splice_map.push("PrimerDimer");
        return splice_map;
    }

    // look here for D1 sequence
    let mut idx = None;
    if sequence[fp_index + 33..50].contains("CGACTG") {
        idx = sequence[fp_index + 33..50].find("CGACTG");
    } else {
        for d in fp_index + 30..50 {
            if close_enough(Some(&sequence[d..d + 9]), "CGGCGACTG", 1) {
                idx = Some(d + 3);
                break;
            }
        }
    }

    if idx.is_none() {
        return splice_map; // failed to find preD1 sequence
    }

    let acceptors: HashMap<&str, &str> = [
        ("GTGAGT", "unspliced"),
        ("GGACAG", "A1"),
        ("AATCTG", "A2"),
        ("AATTGG", "A3"),
        ("TGTTGC", "A4d"),
        ("TTTGTT", "A4c"),
        ("CCTTAG", "A4a"),
        ("GCATCT", "A4b"),
        ("GAAGAA", "A5a"),
        ("AAGCGG", "A5b"),
        ("ACCCAC", "A7"),
    ]
        .iter()
        .cloned()
        .collect();

    let mutated_acceptors: HashMap<&str, &str> = [
        ("GGACAGCAG", "A1"),
        ("AATCTGCTA", "A2"),
        ("AATTGGGTG", "A3"),
        ("TGTTGCTTT", "A4d"),
        ("TTTGTTTCA", "A4c"),
        ("CCTTAGGCA", "A4a"),
        ("GCATCTCCT", "A4b"),
        ("GAAGAAGCG", "A5a"),
        ("AAGCGGAGA", "A5b"),
        ("ACCCACCTC", "A7"),
    ];

    let acceptor_sequence = &sequence[idx.unwrap() + 6..idx.unwrap() + 12];

    if acceptors.contains_key(acceptor_sequence) {
        splice_map.push("D1");
        splice_map.push(acceptors[acceptor_sequence]);
        return splice_map;
    } else if mutated_acceptors.contains_key(acceptor_sequence) {
        splice_map.push("D1");
        splice_map.push(mutated_acceptors[acceptor_sequence]);
        return splice_map;
    } else {
        let unknown_acceptor = &sequence[idx.unwrap() + 6..idx.unwrap() + 16];
        let unknown_acceptor_idx = identify_acceptor(unknown_acceptor, rc_ref);
        unknown_acceptors.push([unknown_acceptor_idx || "no index", unknown_acceptor]);
        return splice_map;
    }
}
