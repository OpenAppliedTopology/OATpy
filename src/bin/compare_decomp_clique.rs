use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind, Indexing};
use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::clique::{Simplex, CliqueComplex};

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

	data_file_name.push_str("/dismat.npy");
	println!("Dissimilarity matrix data is in file: {:?}", data_file_name);
	let arr: Array2<f64> = read_npy(data_file_name).unwrap();
/////////////////// Read dissimilarity matrix data into a Vec<Vec<FilVal>>

	let mut dis_mat = Vec::new();
	let mut min_max = OrderedFloat(0.0);
	for row in arr.outer_iter() {
		let mut vector = Vec::new();
		let mut max = OrderedFloat(0.0);
		for entry in row.iter() {
			vector.push(OrderedFloat(*entry));
			if OrderedFloat(*entry) > max { max = OrderedFloat(*entry); }
		}
		dis_mat.push(vector);
		if max < min_max || min_max == OrderedFloat(0.0) { min_max = max; }
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

    let ringmetadata = RingMetadata{
        ringspec: RingSpec::Modulus(field),
        identity_additive: 0,
        identity_multiplicative: 1,
    };

	let chx = CliqueComplex {
		dissimilarity_matrix: dis_mat,
		dissimilarity_value_max: min_max,
		safe_homology_degrees_to_build_boundaries: (1..(dim+1)).collect(),
		major_dimension: MajorDimension::Row,
		ringmetadata,
		simplex_count: Vec::new()
	};

///////////////////////////////////////////////////////////////////////////////////////////////////

	use serde::{Serialize, Deserialize};
    #[derive(Serialize, Deserialize)]
    struct PairedKeys(f64, Vec<u16>, f64, Vec<u16>);

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
		let major = Simplex {
			filvalue: OrderedFloat(record.0),
			vertices: record.1
		};

		let minor = Simplex {
			filvalue: OrderedFloat(record.2),
			vertices: record.3
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

	let mut maj_to_reduce_vec = Vec::new();
	let mut min_to_reduce_vec = Vec::new();
	let mut maj_to_reduce_set = HashSet::new();
	let mut min_to_reduce_set = HashSet::new();
	while let Some(maj_key) = paired_major_keys.pop() {
		if let Some(min_key) = paired_minor_keys.pop() {
			maj_to_reduce_vec.push(maj_key.clone());
			min_to_reduce_vec.push(min_key.clone());
			maj_to_reduce_set.insert(maj_key);
			min_to_reduce_set.insert(min_key);
		}
	}

	maj_to_reduce_vec.sort();
	min_to_reduce_vec.sort_by(|b, a| a.partial_cmp(b).unwrap());

	let mut keys_dim1 = chx.keys_ordered(1);
	let mut keys_dim0 = chx.keys_ordered(0);
	
	keys_dim1.sort_by(|b, a| a.partial_cmp(b).unwrap());

	println!("Num of all rows: {}", keys_dim0.len());
	println!("Num of all cols: {}", keys_dim1.len());
	println!("Num of paired rows: {}", maj_to_reduce_vec.len());

	let now = Instant::now();
	let (rowoper, indexing_row) = decomp_row(&matrix_row, &mut keys_dim0);
	println!("decomp_row with all rows runs {} nanos and outputs {} pairs.", now.elapsed().as_nanos(), indexing_row.index_2_majkey.len());

	let now = Instant::now();
	let (coloper, indexing_col) = decomp_col(&matrix_col, &mut keys_dim1);
	println!("decomp_col runs {} nanos and outputs {} pairs.", now.elapsed().as_nanos(), indexing_col.index_2_majkey.len());

	let now = Instant::now();
	let (rowoper2, indexing_row_2) = decomp_row_use_pairs(&matrix_row, &mut maj_to_reduce_vec, &mut min_to_reduce_set);
	println!("decomp_row_use_pairs runs {} nanos and outputs {} pairs.", now.elapsed().as_nanos(), indexing_row_2.index_2_majkey.len());

// all columns
// hashset
// resolve the generater bug

	let now = Instant::now();
	let (coloper2, indexing_col_2) = decomp_col_use_pairs(&matrix_col, &mut min_to_reduce_vec, &mut maj_to_reduce_set);
	println!("decomp_col_use_pairs runs {} nanos and outputs {} pairs.", now.elapsed().as_nanos(), indexing_col_2.index_2_majkey.len());

///////////////////////////////////////////////////////////////////////////////////////////////////

	Ok(())
}
