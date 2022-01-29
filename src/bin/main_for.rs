use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::cubical::{Cube, CubicalComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind};

extern crate csv;
extern crate npy_derive;
extern crate npy;

use ndarray::{Array2, Array3};
use ndarray_npy::read_npy;
use std::io::Read;
use npy::NpyData;

use csv::ReaderBuilder;
use std::error::Error;
use std::fs::File;
use std::env;

use ordered_float::OrderedFloat;
use num::rational::Ratio;

use std::time::Instant;

//fn main() {


fn main() -> Result<(), Box<dyn Error>> {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Accessing commond line arguments

	let args: Vec<String> = env::args().collect();
	let dim_file_name = &args[1];
	let array_file_name = &args[2];
	let dim: usize = args[3].trim().parse().expect("Please type a number!");
	let field: usize = args[4].trim().parse().expect("Please type a number!");
	//let mut maxdis: OrderedFloat<f32> = args[5].trim().parse().expect("Please type a number!");
	let print_or_not: bool = args[5].trim().parse().expect("Please type a booliean!");

///////////////////////////////////////////////////////////////////////////////////////////////////
// Read dissimilarity matrix data into a Vec<Vec<FilVal>>
	/*
	let dim_file = File::open(dim_file_name)?;
	let mut dim_reader = ReaderBuilder::new().has_headers(false).from_reader(dim_file);
	let mut dim_iter = dim_reader.deserialize();

	//let array_file = File::open(array_file_name)?;
	//let mut array_reader = ReaderBuilder::new().has_headers(false).from_reader(array_file);
	//let mut array_iter = array_reader.deserialize();

	let mut shape = Vec::new();
	while let Some(result) = dim_iter.next() {
		let record: u16 = result?;
		shape.push(record);
	}

	let mut min_value = OrderedFloat(0.0);
	let mut max_value = OrderedFloat(0.0);
	let mut entry_array:Vec<f64> = Vec::new();
	*/


	/////////////////// Read array data from .csv file
	/*
	if let Some(result) = array_iter.next() {
		let record: f64 = result?;
		min_value = OrderedFloat(record);
		max_value = OrderedFloat(record);
		entry_array.push(OrderedFloat(record));
	}
	while let Some(result) = array_iter.next() {
		let record: f64 = result?;
		entry_array.push(OrderedFloat(record));
		if min_value > OrderedFloat(record) { min_value = OrderedFloat(record); }
		if max_value < OrderedFloat(record) { max_value = OrderedFloat(record); }
	}
	*/
	/////////////////// Read array data from .npy file
	let arr: Array2<f64> = read_npy("grf_anisotropic2d_1000pixperaxis/pixel_births.npy").unwrap();
	//let arr: Array2<f64> = read_npy("pixel_births.npy").unwrap();
	println!("{:?}", arr.dim());

	let mut vec = arr.into_raw_vec();
	println!("{:?}", vec.len());
/*
	let mut buf = vec![];
    std::fs::File::open(array_file_name).unwrap().read_to_end(&mut buf).unwrap();

    //let data: NpyData<f64> = NpyData::from_bytes(&buf).unwrap();

	println!("{:?}", shape);
	println!("{:?}", buf.len());
	for record in buf {
        entry_array.push(record);
    }
	for ii in 0..10 {
		println!("{:?}", entry_array[ii]);
	}
*/
	//println!("{:?}", shape);
///////////////////////////////////////////////////////////////////////////////////
// Read dissimilarity matrix data into a Vec<Vec<FilVal>>
/*
	let file = File::open(input_file)?;
	let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);
	let mut iter = reader.deserialize();

	let mut dismat_sample = Vec::new();
	let mut radius = OrderedFloat(0.0);
	let mut max = OrderedFloat(0.0);
	if let Some(result) = iter.next() {
		let record: Vec<f32> = result?;
		let mut ordered_float_vector = Vec::new();
		for item in record.iter() {
			ordered_float_vector.push(OrderedFloat(item.clone()));
		}
		if let Some(max_value) = ordered_float_vector.iter().max() {
			 max = *max_value;
			 radius = *max_value;
		}
		dismat_sample.push(ordered_float_vector);
	}

	while let Some(result) = iter.next() {
		let record: Vec<f32> = result?;
		let mut ordered_float_vector = Vec::new();
		for item in record.iter() {
			ordered_float_vector.push(OrderedFloat(item.clone()));
		}
		if let Some(max_value) = ordered_float_vector.iter().max() {
			if *max_value > max { max = *max_value; }
			if *max_value < radius { radius = *max_value; }
		}
		dismat_sample.push(ordered_float_vector);
	}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////
// Setup ring data

	let ringmetadata = RingMetadata{
		ringspec: RingSpec::Modulus(field),
		identity_additive: 0,
		identity_multiplicative: 1,
	};

	//println!("Ring: {:?}", ringmetadata.identity_additive);
///////////////////////////////////////////////////////////////////////////////////////////////////
// Compute persistence of clique complex
/*
	if maxdis == OrderedFloat(0.0) { maxdis = radius; }
	let chx = CliqueComplex {
		dissimilarity_matrix: dismat_sample,
		dissimilarity_value_max: maxdis,
		safe_homology_degrees_to_build_boundaries: (1..(dim+1)).collect(),
		major_dimension: MajorDimension::Row,
		ringmetadata: ringdata,
		simplex_count: Vec::new()
	};
*/
///////////////////////////////////////////////////////////////////////////////////////////////////
// Compute persistence of clique complex
	//let entry_array = vec![1,1,1,1,0,1,1,1,1];
	//let mut entry_array = vec![0;27];
	//let array = vec![
	//0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
	//0.1, 10.1, 2.1, 30.1, 4.1, 50.1, 6.1, 70.1, 8.1,
	//0.2, 1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2
	//0.3, 10.1, 2.2, 30.3, 4.4, 50.5, 6.6, 70.7, 8.8,
	//0.4, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8,
	//0.5, 10.1, 2.2, 30.3, 4.4, 50.5, 6.6, 70.7, 8.8,
	//0.6, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8,
	//0.7, 10.1, 2.2, 30.3, 4.4, 50.5, 6.6, 70.7, 8.8,
	//0.8, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8
	//];
	//let array = vec![
	//-1.86,	-1.46,	-1.57,	-1.33,	-1.75,	-2.10,	-2.46,	-3.42,	-3.34,	-2.72,
	//-1.28,	-1.64,	-1.78,	-1.80,	-3.21,	-2.77,	-3.61,	-4.34,	-3.56,	-2.94,
	//-3.20,	-2.80,	-3.49,	-2.29,	-2.46,	-3.21,	-4.11,	-3.82,	-3.44,	-3.06,
	//-4.17,	-4.28,	-3.56,	-2.72,	-3.10,	-3.46,	-3.73,	-3.79,	-3.94,	-4.24,
	//-3.63,	-3.08,	-3.14,	-4.06,	-4.37,	-4.17,	-3.73,	-3.52,	-4.32,	-4.73,
	//-3.23,	-2.88,	-3.30,	-3.74,	-3.76,	-4.20,	-3.72,	-4.35,	-3.26,	-4.54,
	//-3.85,	-3.58,	-3.89,	-3.36,	-3.52,	-4.15,	-3.68,	-4.11,	-3.51,	-4.38,
	//-3.76,	-2.96,	-2.54,	-3.39,	-3.21,	-3.56,	-3.39,	-3.15,	-4.25,	-4.28,
	//-4.17,	-2.63,	-2.52,	-3.31,	-3.51,	-3.89,	-3.45,	-3.98,	-4.37,	-4.17,
	//-3.82,	-3.10,	-2.80,	-3.44,	-2.71,	-3.74,	-3.92,	-4.19,	-3.48,	-4.44,
	//];

	//let mut entry_array = Vec::new();
	//for val in array.iter() {
	//	entry_array.push(OrderedFloat(val.clone()));
	//}
	//entry_array = vec![70;27]-entry_array;
	//entry_array[10] = 1; entry_array[12] = 2; entry_array[14] = 3; entry_array[16] = 4;
/*
	let chx = CubicalComplex {
		ringmetadata,
		shape: vec![20,20],
		entry_array,
		max_value: OrderedFloat(70.0),
		safe_homology_degrees_to_build_boundaries: (1..(dim+1)).collect(),
		major_dimension: MajorDimension::Row
	};
*/



/*
	let chx = CubicalComplex {
		ringmetadata,
		shape,
		entry_array,
		max_value,
		safe_homology_degrees_to_build_boundaries: (1..(dim+1)).collect(),
		major_dimension: MajorDimension::Row
	};

///////////////////////////////////////////////////////////////////////////////////////////////////
// Print out the outputs

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

*/


///////////////////////////////////////////////////////////////////////////////////////////////////
// Write pivot pairs to a file




///////////////////////////////////////////////////////////////////////////////////////////////////
/*
	let adjacency_matrix = matrix.adjacency_matrix();
	println!("adjacency_matrix:");
	for ii in 1..50{
		println!("{:?}", adjacency_matrix[ii]);
	}
*/
/*
	let mut iter = chx.keys_unordered_itr(1);
	println!("Here:");
	for ii in 0..100 {
		if let Some(cube) = iter.next() {
			println!("{:?}", cube);
		}
	}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////
// Test for cofaces and faces
	//let matrix = chx.get_smoracle(MajorDimension::Row, ChxTransformKind::Boundary);
	//let majkey = &Cube { filvalue: OrderedFloat(-4.240518), coordinates: vec![0, 16] };

	//let majkey = &Cube{filvalue: matrix.disvec[106].clone(), vertices: vec![2,6]};
	//let majkey = &Simplex{filvalue: OrderedFloat(0.0), vertices: vec![ 6, 7, 8]};
	//println!("cofaces:");
	//for (key, val) in matrix.maj_itr(&majkey) {
	//	println!("key={:?}, val={:?}", key, val);
	//}
///////////////////////////////////////////////////////////////////////////////////////////////////
	//let indexing = &blocks.dim_indexing[2];
	//println!("matched pairs:");
	//for ind in 0..indexing.index_2_majkey.len(){
	//	majkey = &indexing.index_2_majkey[ind];
	//	let minkey = &indexing.index_2_minkey[ind];
	//	println!("{:?}	{:?}", majkey.clone(), minkey);
	//}

	//println!("indexing.minkey_2_index: {:?}", indexing.minkey_2_index);
	//println!("indexing.majkey_2_index: {:?}", indexing.majkey_2_index);

	//println!("max = {}", max);
	//println!("rad = {}", radius);
	//let simp = Simplex{filvalue: OrderedFloat(1.0), vertices: vec![1, 2]};
	//let matched_basis = blocks.get_matched_basis_vector(1, majkey);
	//println!("matched basis: {:?}", matched_basis);

	// Ensure that we got the original array back
	//assert_eq!(array_read, array);

/*
	let barfile = File::open(barcode_file)?;
	let mut rdr = ReaderBuilder::new().has_headers(false).from_reader(barfile);
	let mut iter = rdr.deserialize();

	let mut barcode_1 = Vec::new();
	while let Some(result) = iter.next() {
		let record: (f64, f64) = result?;
		barcode_1.push((OrderedFloat(record.0), OrderedFloat(record.1)));
	}
	barcode_1.sort();
	barcode.sort();
	println!("Contains {} bars.", barcode_1.len());
	assert_eq!(barcode, barcode_1);
*/
/*
	for (OrderedFloat(a), OrderedFloat(b)) in barcode_1.iter() {
	println!("{:?},{:?}", a, b);
	}
*/
/*
	println!("Paired bars:");
	for ii in 0..barcode.len() {
		let (OrderedFloat(a), OrderedFloat(b)) = barcode[ii];
		let (OrderedFloat(c), OrderedFloat(d)) = barcode_1[ii];
		println!("ii = {}", ii);
		println!("{:?}", (a,b));
		println!("{:?}", (c,d));
		//println!("{:?} vs. {:?}", (a,b), (c,d));
/*
		if barcode[ii] != barcode_1[ii] {
			println!("{:?}", barcode[ii]);
			println!("{:?}", barcode_1[ii]);
		}
*/
	}
*/
	Ok(())
}


/*

// -------------------------------------------------------------------------


    use persistence::clique::Simplex;
    use std;


    // ----------------------------------------------------------------------------------
    // Set maximum threshold values for homology dimension and dissimilarity
    // ----------------------------------------------------------------------------------
    let dim = 3;
    let maxdis = 2;


    // ----------------------------------------------------------------------------------
    // Define an object to represent the ring Z/3Z
    // ----------------------------------------------------------------------------------
    let ringmetadata = persistence::matrix::RingMetadata{
    	ringspec: RingSpec::Modulus(3),
    	identity_additive: 0,
    	identity_multiplicative: 1,
    };


    // ----------------------------------------------------------------------------------
    // Build a "dissimilarity matrix" as a vector of vectors
    // ----------------------------------------------------------------------------------
    let dismat = vec![  vec![0,  1,  2,  1],
                        vec![1,  0,  1,  2],
                        vec![2,  1,  0,  1],
                        vec![1,  2,  1,  0]  ];


    // ----------------------------------------------------------------------------------
    // Construct the corresponding filtered clique complex
    // ----------------------------------------------------------------------------------
    let chx = persistence::clique::CliqueComplex {
        // the distance/dissimilarity matrix
        dissimilarity_matrix: dismat,
        // threshold to stop the filtration
        dissimilarity_value_max: maxdis,
        // sets "safeguards" on dimension; we'll get warnings if we try to
        // get boundary matrices in dimension higher than dim+1
        safe_homology_degrees_to_build_boundaries: (1..dim+1).collect(),
        // set the default major dimension (for sparse matrices) to be row
        major_dimension: MajorDimension::Row,
        // indicates we want Z/3Z coefficients
        ringmetadata: ringmetadata,
        // don't worry about this
        simplex_count: Vec::new()
    };


    // ----------------------------------------------------------------------------------
    // Get a (row-major) sparse matrix oracle for the boundary operator
    // ----------------------------------------------------------------------------------
    let D =  chx.get_smoracle(persistence::matrix::MajorDimension::Row,
                              persistence::chx::ChxTransformKind::Boundary);


    // ----------------------------------------------------------------------------------
    // Define a weighted 0-simplex and 1-simplex
    // ----------------------------------------------------------------------------------
    let simplex_d1 = persistence::clique::Simplex{
        filvalue: 1,
        vertices : vec![0, 1]
    };
    let simplex_d0 = persistence::clique::Simplex{
        filvalue: 0,
        vertices : vec![0]
    };


    // ----------------------------------------------------------------------------------
    // Compute / check the filtration values
    // ----------------------------------------------------------------------------------
    std::assert_eq!(chx.key_2_filtration( &simplex_d0 ), 0);
    std::assert_eq!(chx.key_2_filtration( &simplex_d1 ), 1);


    // ----------------------------------------------------------------------------------
    // Access a row of the boundary matrix
    //  - the boundary matrix is "row-major" so we call rows "major fields"
    // ----------------------------------------------------------------------------------
    // create an iterator that runs over the structural nonzero entries of the row
    // each item returned by the iterator is a tuple of form (weighted_simplex, coefficient)
    let major_field = D.maj_itr( &simplex_d0 );


    // the following for-loop should print the following:
    // ```text
    // Structural nonzero entries of a row corresponding to a 0-simplex:
    // (Simplex { filvalue: 1, vertices: [0, 1] }, -1)
    // (Simplex { filvalue: 2, vertices: [0, 2] }, -1)
    // (Simplex { filvalue: 1, vertices: [0, 3] }, -1)
    // ```
    println!("Structural nonzero entries of a row corresponding to a 0-simplex:");
    for item in major_field  {
        println!("{:?}", item);
    }

    // check to ensure that the output is correct:
    let mut correct_val : Vec< (Simplex<i64>, i16) >  =
                      vec![ (Simplex{ filvalue: 1, vertices: vec![0, 1] }, -1),
                            (Simplex{ filvalue: 2, vertices: vec![0, 2] }, -1),
                            (Simplex{ filvalue: 1, vertices: vec![0, 3] }, -1) ];

    let major_field2 = D.maj_itr( &simplex_d0 ); // this re-creates the iterator
    std::assert_eq!( major_field2.eq(correct_val.iter().map( |x| x.clone() ) ), true);


    // ----------------------------------------------------------------------------------
    // Access a column of the boundary matrix
    //  - the boundary matrix is "row-major" so we call columns "column fields"
    // ----------------------------------------------------------------------------------
    // create an iterator that runs over the structural nonzero entries of the row
    // each item returned by the iterator is a tuple of form (weighted_simplex, coefficient)
    let minor_field = D.min_itr( &simplex_d1 );

    // the following for-loop should print the following:
    // ```text
    // Structural nonzero entries of a column corresponding to a 1-simplex:
    // (Simplex { filvalue: 0, vertices: [1] }, 1)
    // (Simplex { filvalue: 0, vertices: [0] }, -1)
    // ```
    println!("Structural nonzero entries of a column corresponding to a 1-simplex:");
    for item in minor_field  {
        println!("{:?}", item);
    }

    // check to ensure that the output is correct:
    let mut correct_val : Vec< (Simplex<i64>, i16) >  =
                      vec![  (Simplex{ filvalue: 0, vertices: vec![1] },  1),
                             (Simplex{ filvalue: 0, vertices: vec![0] }, -1)  ];

    let minor_field2 = D.min_itr( &simplex_d1 ); // this re-creates the iterator
    std::assert_eq!( minor_field2.eq(correct_val.iter().map( |x| x.clone() ) ), true);


    // ----------------------------------------------------------------------------------
    // Factor the complex
    // ----------------------------------------------------------------------------------
    let factored_complex = persistence::chx::factor_chain_complex(&chx, dim+1);


    // ----------------------------------------------------------------------------------
    // Read barcodes / check for correctness
    // ----------------------------------------------------------------------------------
    // predefine the set of correct solutions
    let correct_barcodes = vec![    vec![ (0, 2), (0, 1), (0, 1), (0, 1) ], // dimension 0
                                    vec![ (1, 2) ],                         // dimension 1
                                    vec![],                                 // dimension 2
                                    vec![]                                  // dimension 3
                                ];

    // confirm that correct solutions and actual solutions are the same
    for i in 0..dim+1 {
        assert_eq!( correct_barcodes[i], factored_complex.barcode(i) );
    }


    // ----------------------------------------------------------------------------------
    // Get the cycle representative associated with a (weighted) simplex
    //    - this may corresopnd to a "length-0 bar", depending on choice
    // ----------------------------------------------------------------------------------
    // get the vector *represented as a hashmap*
    let basis_vec = factored_complex.get_matched_basis_vector(1, &simplex_d1);

    // this for-loop should print the following:
    // ```text
    // Basis vector corresponding to a 1-simplex:
    // (Simplex { filvalue: 1, vertices: [0, 1] }, 1)
    // ```
    println!("Basis vector corresponding to a 1-simplex:");
    for item in basis_vec.iter()  {
        println!("{:?}", item);
    }

    // check to ensure that the output is correct:
        // write a vector with correct values
    let mut correct_val : Vec< (Simplex<i64>, i16) >  =
                      vec![  (Simplex{ filvalue: 1, vertices: vec![0, 1] },  1)  ];

        // convert the hashmap to an iterator

        // remark: if we had tried to create this iterator via a command
        // > factored_complex.get_matched_basis_vector(1, &simplex_d1).iter().map( etc. )
        // then we would have gotten an error; the reason is that by declaring the `basis_vec`
        // variable name we signaled rust that it should have a certain lifetime; see
        // https://stackoverflow.com/questions/54056268/temporary-value-is-freed-at-the-end-of-this-statement/54056716#54056716
        // for further discussion
    let basis_vec_iter = basis_vec
                            .iter() // get iterator for the hashmap
                            .map( |x| ( x.0.clone(), x.1.clone() ) ); // replace ref's with clones

        // check equality
    std::assert_eq!( basis_vec_iter.eq( correct_val.iter().cloned() ) , true);


    // ----------------------------------------------------------------------------------
    // Get the birth and death filtraitons of a chain
    //    - here "chain" is formalized as a hashmap mapping keys to coefficients
    // ----------------------------------------------------------------------------------
    // this part is not implemented yet.
}
*/
