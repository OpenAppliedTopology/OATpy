use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::clique::{Simplex, CliqueComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind, Indexing};

extern crate csv;

#[macro_use]
extern crate npy_derive;
extern crate npy;

use ndarray::{Array2, Array3};
use ndarray_npy::{read_npy, write_npy};
use tuple_conv::RepeatedTuple;

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
	//let mut maxdis: OrderedFloat<f32> = args[4].trim().parse().expect("Please type a number!");
	let print_or_not: bool = args[4].trim().parse().expect("Please type a booliean!");

	data_file_name.push_str("/dismat.npy");
	println!("Dissimilarity matrix data is in file: {:?}", data_file_name);
	let arr: Array2<f64> = read_npy(&data_file_name).unwrap();
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

	println!("Testing ...");
    let now = Instant::now();
    let blocks = factor_chain_complex(&chx, dim+1);
    println!("Runs {} seconds", now.elapsed().as_secs());
    println!("Done!");

    let mut barcode = blocks.barcode(dim);
    barcode.sort();
    println!("Barcode contains {} bars.", barcode.len());
    if print_or_not == true {
        for (start, end) in barcode.iter() {
            println!("{},{}", start, end);
        }
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

	use serde::{Serialize, Deserialize};
    #[derive(Serialize, Deserialize)]
    struct PairedKeys(f64, Vec<u16>, f64, Vec<u16>);

    let maj_to_reduce = blocks.dim_indexing[dim].index_2_majkey.clone();
    let min_to_reduce = blocks.dim_indexing[dim].index_2_minkey.clone();

	let mut save_as = args[1].clone();
	save_as.push_str("/pairs_dim");
	save_as.push_str(&dim.to_string());
	save_as.push_str(".csv");
	println!("Pairs are saved in file: {}", save_as);
    let mut file = File::create(save_as).unwrap();

    for ind in 0..maj_to_reduce.len(){
        let OrderedFloat(maj_val) = maj_to_reduce[ind].filvalue;
		let OrderedFloat(min_val) = min_to_reduce[ind].filvalue;
        let tem_keys = PairedKeys(
			maj_val,
			maj_to_reduce[ind].vertices.clone(),
			min_val,
			min_to_reduce[ind].vertices.clone(),
		);

        let mut tmp = serde_json::to_string(&tem_keys).unwrap();
        tmp += "\n";
        file.write(tmp.to_string().as_bytes()).unwrap();
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

    Ok(())
}
