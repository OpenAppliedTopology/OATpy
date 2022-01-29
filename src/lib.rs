
//! A library for applied algebra and topology.

//! # Summary
//!
//!  This package is designed for work with chain complexes in topological data analysis.
//!
//! These are **big** objects -- they have matrices with billions or even trillions of rows and columns.  Most computers can't hold that much information in memory.
//! Luckily, we don't often need to store that much information at one time.  Many jobs can be accomplished using just a few rows or columns -- or by processing a sequence of rows/columns, one at a time.  Once we're done with each piece we can clear it from memory.  What we really need is a  **matrix oracle** -- something that provides just the information we need, when we need it.
//!
//!
//! This package provides a streamlined system for working with oracles in algebraic topology.
//!
//  There are many types of chain complexes (simplicial, cubical, etc.) and different types of oracles work better for each one.  However, it's best if all the different oracles obey the same, standardized set of commands (you don't want to have to learn a new command to get the row of a sparse matrix every time someone creates a new oracle).  In rust, you can standardize a common set of commands so that they work the same way on multiple different objects by using a [trait](https://doc.rust-lang.org/book/ch10-02-traits.html). In fact, that's all a trait really is.

//!
//! * reference information for [sparse matrix oracles](matrix)
//! * reference information for [chain complex oracles](chx)
//!

//! # Sparse matrix oracles
//!
//!

//! A *sparse matrix oracle* is an object that allows users to access rows, columns, and other information about a sparse matrix.  This library provides a standardized set of commands that you can use to access information from any sparse matrix oracle (you don't need different commands for different oracles).  This standard set of commands is formalized as a [Rust trait](https://doc.rust-lang.org/book/ch10-02-traits.html).

//! * a [Rust trait](https://doc.rust-lang.org/book/ch10-02-traits.html) is a set of commands that can be applied in the same way to several different objects
//! * the [`SmOracle`](crate::matrix::SmOracle) trait is a set of commands for accessing information about a sparse matrix
//! * a sparse matrix oracle is defined as any [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that implements the commands in the [`SmOracle`](crate::matrix::SmOracle) trait
//!

// The [`SmOracle`](matrix/trait.SmOracle.html) trait is a standard set of commands for interacting with sparse matrix oracles.  It gives you access to rows, columns, and other information.
//
// * Definition and examples
// 	* A sparse matrix oracle is a  [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that implements the [`SmOracle`](matrix/trait.SmOracle.html) trait.
// 	* A list of sparse matrix oracles can be found at the bottom of the [`SmOracle`](matrix/trait.SmOracle.html) documentation page.
// 	* **You can also build your own**!  In concrete terms, this means making a [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that [implements](https://doc.rust-lang.org/book/ch10-02-traits.html) the [`SmOracle`](crate::matrix::SmOracle) trait.
// * Operations
// 	* There are serval matrix operations you can perform with matrix oracles; for example, you can multiply a matrix with a sparse vector, obtain an *R = DV* decomposition, and solve *Ax = b* for *x*.  See the [`solver`](solver) module for more details.




//!
//!
//!	```
//!	Example 1: (to fill in later)
//! 	- create a row-major csm
//! 	- access a row of the csm
//!		- access a column of the csm
//! 	- perform triangular solve
//! 	- check the triangular solve by confirming that Ax = b
//! 	- perform R=DV factorization
//! ```
//!


//! # Chain complex oracles
//!
//! A *chain complex oracle* is an object that allows users to access information about a chain complex:
//! * matrices: boundary matrices, change of basis matrices, etc.
//! * indices: row/column indices, pivot indices, etc.
//! * (co)cycle representatives
//! * barcodes, persistence diagrams, birth/death times, etc.
//!
//! This library provides a universal set of commands to work with chain complex oracles.  Just like the commands for matrix oracles, this set of commands is formalized as a Rust trait.
//!
//! * The [`ChainComplex`](crate::chx::ChainComplex) trait is a standard set of commands for  interacting with filtered chain complexes (simplicial, cubical, etc.).
//! * A chain complex oracle is defined as any [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that implements the commands in the [`ChainComplex`](crate::chx::ChainComplex) trait
//!



// * Definition and examples
// 	* A chain complex oracle is a [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that implements the [`ChainComplex`](crate::chx::ChainComplex) trait.
// 	* A list of existing chain complex oracles, including [clique](clique) and [cubical](cubical) complexes, can be found at the bottom of the documentation page for the [`ChainComplex`](crate::chx::ChainComplex) trait.
// 	* **You can also build your own**!  In concrete terms, this means making a [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that [implements](https://doc.rust-lang.org/book/ch10-02-traits.html) the [`ChainComplex`](crate::chx::ChainComplex) trait.
//

// * Information provided by chain complex oracles includes:
//	* matrices: boundary matrices, change of basis matrices, etc.
// 	* indices: row/column indices, pivot indices, etc.
// 	* (co)cycle representatives
//	* barcodes, persistence diagrams, birth/death times, etc.
//


//! ## Examples
//!
//! Work with a filtered clique complex
//!
// * build a complex
// * access rows/columns of the boundary matrix
// * compute barcodes
// * compute cycle representatives


//! ## WARNING IF YOU COPY/PASTE THIS TEXT FROM A WEB BROWSER IT MAY CHANGE ALL CAPS TO LOWER CASE,
//! CAUSING THE CODE TO FAIL; INSTEAD COPY THIS EXAMPLE FROM THE lib.rs FILE AND REMOVE THE
//! COMMENTS MARKERS //!
//! ```
//! use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
//! use exhact::cubical::{Cube, CubicalComplex};
//! use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind};
//! use exhact::clique::Simplex;
//! use std;
//!
//!
//! fn main() {
//!
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Set maximum threshold values for homology dimension and dissimilarity
//!     // ----------------------------------------------------------------------------------
//!     let dim = 3;
//!     let maxdis = 2;
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Define an object to represent the ring Z/3Z
//!     // ----------------------------------------------------------------------------------
//!     let ringmetadata = exhact::matrix::RingMetadata{
//!     	ringspec: RingSpec::Modulus(3),
//!     	identity_additive: 0,
//!     	identity_multiplicative: 1,
//!     };
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Build a "dissimilarity matrix" as a vector of vectors
//!     // ----------------------------------------------------------------------------------
//!     let dismat = vec![  vec![0,  1,  2,  1],
//!                         vec![1,  0,  1,  2],
//!                         vec![2,  1,  0,  1],
//!                         vec![1,  2,  1,  0]  ];
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Construct the corresponding filtered clique complex
//!     // ----------------------------------------------------------------------------------
//!     let chx = exhact::clique::CliqueComplex {
//!         // the distance/dissimilarity matrix
//!         dissimilarity_matrix: dismat,
//!         // threshold to stop the filtration
//!         dissimilarity_value_max: maxdis,
//!         // sets "safeguards" on dimension; we'll get warnings if we try to
//!         // get boundary matrices in dimension higher than dim+1
//!         safe_homology_degrees_to_build_boundaries: (1..dim+1).collect(),
//!         // set the default major dimension (for sparse matrices) to be row
//!         major_dimension: MajorDimension::Row,
//!         // indicates we want Z/3Z coefficients
//!         ringmetadata: ringmetadata,
//!         // don't worry about this
//!         simplex_count: Vec::new()
//!     };
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Get a (row-major) sparse matrix oracle for the boundary operator
//!     // ----------------------------------------------------------------------------------
//!     let D =  chx.get_smoracle(exhact::matrix::MajorDimension::Row,
//!                               exhact::chx::ChxTransformKind::Boundary);
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Define a weighted 0-simplex and 1-simplex
//!     // ----------------------------------------------------------------------------------
//!     let simplex_d1 = exhact::clique::Simplex{
//!         filvalue: 1,
//!         vertices : vec![0, 1]
//!     };
//!     let simplex_d0 = exhact::clique::Simplex{
//!         filvalue: 0,
//!         vertices : vec![0]
//!     };
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Compute / check the filtration values
//!     // ----------------------------------------------------------------------------------
//!     std::assert_eq!(chx.key_2_filtration( &simplex_d0 ), 0);
//!     std::assert_eq!(chx.key_2_filtration( &simplex_d1 ), 1);
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Access a row of the boundary matrix
//!     //  - the boundary matrix is "row-major" so we call rows "major fields"
//!     // ----------------------------------------------------------------------------------
//!     // create an iterator that runs over the structural nonzero entries of the row
//!     // each item returned by the iterator is a tuple of form (weighted_simplex, coefficient)
//!     let major_field = D.maj_itr( &simplex_d0 );
//!
//!
//!     // the following for-loop should print the following:
//!     // ```text
//!     // Structural nonzero entries of a row corresponding to a 0-simplex:
//!     // (Simplex { filvalue: 1, vertices: [0, 1] }, -1)
//!     // (Simplex { filvalue: 2, vertices: [0, 2] }, -1)
//!     // (Simplex { filvalue: 1, vertices: [0, 3] }, -1)
//!     // ```
//!     println!("Structural nonzero entries of a row corresponding to a 0-simplex:");
//!     for item in major_field  {
//!         println!("{:?}", item);
//!     }
//!
//!     // check to ensure that the output is correct:
//!     let mut correct_val : Vec< (Simplex<i64>, i16) >  =
//!                       vec![ (Simplex{ filvalue: 1, vertices: vec![0, 1] }, -1),
//!                             (Simplex{ filvalue: 2, vertices: vec![0, 2] }, -1),
//!                             (Simplex{ filvalue: 1, vertices: vec![0, 3] }, -1) ];
//!
//!     let major_field2 = D.maj_itr( &simplex_d0 ); // this re-creates the iterator
//!     std::assert_eq!( major_field2.eq(correct_val.iter().map( |x| x.clone() ) ), true);
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Access a column of the boundary matrix
//!     //  - the boundary matrix is "row-major" so we call columns "column fields"
//!     // ----------------------------------------------------------------------------------
//!     // create an iterator that runs over the structural nonzero entries of the row
//!     // each item returned by the iterator is a tuple of form (weighted_simplex, coefficient)
//!     let minor_field = D.min_itr( &simplex_d1 );
//!
//!     // the following for-loop should print the following:
//!     // ```text
//!     // Structural nonzero entries of a column corresponding to a 1-simplex:
//!     // (Simplex { filvalue: 0, vertices: [1] }, 1)
//!     // (Simplex { filvalue: 0, vertices: [0] }, -1)
//!     // ```
//!     println!("Structural nonzero entries of a column corresponding to a 1-simplex:");
//!     for item in minor_field  {
//!         println!("{:?}", item);
//!     }
//!
//!     // check to ensure that the output is correct:
//!     let mut correct_val : Vec< (Simplex<i64>, i16) >  =
//!                       vec![  (Simplex{ filvalue: 0, vertices: vec![1] },  1),
//!                              (Simplex{ filvalue: 0, vertices: vec![0] }, -1)  ];
//!
//!     let minor_field2 = D.min_itr( &simplex_d1 ); // this re-creates the iterator
//!     std::assert_eq!( minor_field2.eq(correct_val.iter().map( |x| x.clone() ) ), true);
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Factor the complex
//!     // ----------------------------------------------------------------------------------
//!     let factored_complex = exhact::chx::factor_chain_complex(&chx, dim+1);
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Read barcodes / check for correctness
//!     // ----------------------------------------------------------------------------------
//!     // predefine the set of correct solutions
//!     let correct_barcodes = vec![    vec![ (0, 2), (0, 1), (0, 1), (0, 1) ], // dimension 0
//!                                     vec![ (1, 2) ],                         // dimension 1
//!                                     vec![],                                 // dimension 2
//!                                     vec![]                                  // dimension 3
//!                                 ];
//!
//!     // confirm that correct solutions and actual solutions are the same
//!     for i in 0..dim+1 {
//!         assert_eq!( correct_barcodes[i], factored_complex.barcode(i) );
//!     }
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Get the cycle representative associated with a (weighted) simplex
//!     //    - this may corresopnd to a "length-0 bar", depending on choice
//!     // ----------------------------------------------------------------------------------
//!     // get the vector *represented as a hashmap*
//!     let basis_vec = factored_complex.get_matched_basis_vector(1, &simplex_d1);
//!
//!     // this for-loop should print the following:
//!     // ```text
//!     // Basis vector corresponding to a 1-simplex:
//!     // (Simplex { filvalue: 1, vertices: [0, 1] }, 1)
//!     // ```
//!     println!("Basis vector corresponding to a 1-simplex:");
//!     for item in basis_vec.iter()  {
//!         println!("{:?}", item);
//!     }
//!
//!     // check to ensure that the output is correct:
//!         // write a vector with correct values
//!     let mut correct_val : Vec< (Simplex<i64>, i16) >  =
//!                       vec![  (Simplex{ filvalue: 1, vertices: vec![0, 1] },  1)  ];
//!
//!         // convert the hashmap to an iterator
//!
//!         // remark: if we had tried to create this iterator via a command
//!         // > factored_complex.get_matched_basis_vector(1, &simplex_d1).iter().map( etc. )
//!         // then we would have gotten an error; the reason is that by declaring the `basis_vec`
//!         // variable name we signaled rust that it should have a certain lifetime; see
//!         // https://stackoverflow.com/questions/54056268/temporary-value-is-freed-at-the-end-of-this-statement/54056716#54056716
//!         // for further discussion
//!     let basis_vec_iter = basis_vec
//!                             .iter() // get iterator for the hashmap
//!                             .map( |x| ( x.0.clone(), x.1.clone() ) ); // replace ref's with clones
//!
//!         // check equality
//!     std::assert_eq!( basis_vec_iter.eq( correct_val.iter().cloned() ) , true);
//!
//!
//!     // ----------------------------------------------------------------------------------
//!     // Get the birth and death filtraitons of a chain
//!     //    - here "chain" is formalized as a hashmap mapping keys to coefficients
//!     // ----------------------------------------------------------------------------------
//!     // this part is not implemented yet.
//! }
//! ```





pub mod csm;
pub mod matrix;

pub mod chx;
pub mod solver;

pub mod clique;
pub mod cubical;

pub mod decomp_row;
pub mod decomp_col;
pub mod decomp_row_use_pairs;
pub mod decomp_col_use_pairs;
