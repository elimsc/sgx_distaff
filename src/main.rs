#![cfg_attr(not(target_env = "sgx"), no_std)]
#![cfg_attr(
    all(target_env = "sgx", target_vendor = "mesalock"),
    feature(rustc_private)
)]

#[cfg(not(target_env = "sgx"))]
#[macro_use]
extern crate sgx_tstd as std;
use std::prelude::v1::*;

extern crate distaff;
extern crate tiny_keccak;

use distaff::{assembly, ProgramInputs, ProofOptions};
use tiny_keccak::{Hasher, Sha3};

fn main() {
    // // this is our program, we compile it from assembly code
    // let program = assembly::compile("begin push.3 push.5 add end").unwrap();

    // // let's execute it
    // let (outputs, proof) = distaff::execute(
    //     &program,
    //     &ProgramInputs::none(), // we won't provide any inputs
    //     1,                      // we'll return one item from the stack
    //     &ProofOptions::default(),
    // ); // we'll be using default options

    // // the output should be 8
    // assert_eq!(vec![8], outputs);
    let mut sha3 = Sha3::v256();
    let mut output = [0u8; 32];
    let expected = b"\
        \x64\x4b\xcc\x7e\x56\x43\x73\x04\x09\x99\xaa\xc8\x9e\x76\x22\xf3\
        \xca\x71\xfb\xa1\xd9\x72\xfd\x94\xa3\x1c\x3b\xfb\xf2\x4e\x39\x38\
    ";

    sha3.update(b"hello");
    sha3.update(b" ");
    sha3.update(b"world");
    sha3.finalize(&mut output);

    assert_eq!(expected, &output);
}
