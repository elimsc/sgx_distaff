use distaff::{ Program, ProgramInputs, ProofOptions };
use std::prelude::v1::*;
mod utils;

pub mod collatz;
pub mod comparison;
pub mod conditional;
pub mod fibonacci;
pub mod merkle;
pub mod range;

pub struct Example {
    pub program         : Program,
    pub inputs          : ProgramInputs,
    pub num_outputs     : usize,
    pub options         : ProofOptions,
    pub expected_result : Vec<u128>
}