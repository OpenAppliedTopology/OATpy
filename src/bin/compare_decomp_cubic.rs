use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind, Indexing};
use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::cubical::{Cube, CubicalComplex};

use exhact::decomp_row_use_pairs::decomp_row_use_pairs;
use exhact::decomp_col_use_pairs::decomp_col_use_pairs;
use exhact::decomp_row::decomp_row;
use exhact::decomp_col::decomp_col;

use exhact::csm::CSM;

extern crate csv;

#[macro_use]
extern crate npy_derive;
extern crate npy;

use ndarray::{Array2, Array3};
use ndarray_npy::{read_npy, write_npy};
use tuple_conv::RepeatedTuple;

use std::collections::{HashMap, HashSet};

use std::io::{Read, Write, BufReader, BufRead};
use npy::NpyData;

use csv::ReaderBuilder;
use std::error::Error;
use std::fs::File;
use std::env;

use ordered_float::OrderedFloat;
use num::rational::Ratio;

use std::time::Instant;

fn main() -> Result<(), Box<dyn Error>> {
///////////////////////////////////////////////////////////////////////////////////////////////////
// Accessing commond line arguments

	let args: Vec<String> = env::args().collect();
	let mut data_file_name = args[1].clone();
	let dim: usize = args[2].trim().parse().expect("Please type a number!");
	let field: usize = args[3].trim().parse().expect("Please type a number!");
	let print_or_not: bool = args[4].trim().parse().expect("Please type a booliean!");

	data_file_name.push_str("/pixel_births.npy");
	println!("Array data is in file: {:?}", data_file_name);
	let arr: Array3<f64> = read_npy(data_file_name).unwrap();
/////////////////// Read array data from .npy file

    let mut entry_array = Vec::new();
    let mut max_value = OrderedFloat(0.0);
    for entry in arr.iter() {
        entry_array.push(OrderedFloat(*entry));
        if max_value < OrderedFloat(*entry) { max_value = OrderedFloat(*entry); }
    }
	println!("The shape of matrix data is: {:?}", arr.dim());

///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
    let ringmetadata = RingMetadata{
        ringspec: RingSpec::Modulus(field),
        identity_additive: 0,
        identity_multiplicative: 1,
    };

    let chx = CubicalComplex {
		ringmetadata,
		shape: arr.dim().to_vec(),
		entry_array,
		max_value,
		safe_homology_degrees_to_build_boundaries: (1..(dim+1)).collect(),
		major_dimension: MajorDimension::Row
	};

///////////////////////////////////////////////////////////////////////////////////////////////////

	use serde::{Serialize, Deserialize};
    #[derive(Serialize, Deserialize)]
    struct PairedKeys(f64, Vec<u32>, f64, Vec<u32>);

	let mut pairs_file_name = args[1].clone();
	pairs_file_name.push_str("/pairs_dim");
	pairs_file_name.push_str(&dim.to_string());
	pairs_file_name.push_str(".csv");
	println!("Reading pairs from file: {}", pairs_file_name);

	let mut paired_major_keys = Vec::new();
	let mut paired_minor_keys = Vec::new();

	let paired_keys_file = File::open(pairs_file_name)?;
	let buffered = BufReader::new(paired_keys_file);

	for line in buffered.lines() {
		let record: PairedKeys = serde_json::from_str(&line?).unwrap();
		let major = Cube {
			filvalue: OrderedFloat(record.0),
			coordinates: record.1
		};

		let minor = Cube {
			filvalue: OrderedFloat(record.2),
			coordinates: record.3
		};

		paired_major_keys.push(major);
		paired_minor_keys.push(minor);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

	let matrix_row = chx.get_smoracle(
		MajorDimension::Row,
		ChxTransformKind::Boundary
	);

	let matrix_col = chx.get_smoracle(
		MajorDimension::Col,
		ChxTransformKind::Boundary
	);


	let mut maj_to_reduce = Vec::new();
	let mut min_to_reduce = HashSet::new();
	while let Some(maj_key) = paired_major_keys.pop() {
		if let Some(min_key) = paired_minor_keys.pop() {
			maj_to_reduce.push(maj_key);
			min_to_reduce.insert(min_key);
		}
	}

	let mut keys_dim1 = chx.keys_ordered(1);
	let mut keys_dim0 = chx.keys_ordered(0);

	println!("Num of all rows: {}", keys_dim0.len());
	println!("Num of paired rows: {}", maj_to_reduce.len());

	let now = Instant::now();
	let (rowoper, indexing) = decomp_row(&matrix_row, &mut keys_dim0);
	println!("Runs {} nanos", now.elapsed().as_nanos());
/*
	let now = Instant::now();
	let (rowoper2, indexing2) = decomp_row_use_pairs(&matrix, &mut maj_to_reduce, &mut min_to_reduce);
	println!("Runs {} nanos", now.elapsed().as_nanos());
*/
	//let now = Instant::now();
	//let (rowoper3, indexing3) = decomp_col_use_pairs(&matrix, &mut min_to_reduce, &mut maj_to_reduce);
	//println!("Runs {} nanos", now.elapsed().as_nanos());
///////////////////////////////////////////////////////////////////////////////////////////////////

	Ok(())
}
