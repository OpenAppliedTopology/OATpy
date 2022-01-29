



/*!

Chain complex oracles (for filtered, based chain complexes)

# Chain complex oracles


A *chain complex oracle* is an object that allows users to access boundary matrices, birth times, and other information about a chain complex.  This library provides a standardized set of commands that you can use to access information from any chain complex oracle (you don't need different commands for different oracles).  This standard set of commands is formalized as a [Rust trait](https://doc.rust-lang.org/book/ch10-02-traits.html).

* a [Rust trait](https://doc.rust-lang.org/book/ch10-02-traits.html) is a set of commands that can be applied in the same way to several different objects
* the [`ChainComplex`](crate::chx::ChainComplex) trait is set of commands for accessing information about a chain complex
* a chain complex oracleis defined as any [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that implements the commands in the [`ChainComplex`](crate::matrix::SmOracle) trait

# Factoring
Some information (e.g. persistence diagrams, cycle representatives) won't be available right away.  For this you'll need to *factor* the complex.  You can do this yourself, or you can use the [`factor_chain_complex`](factor_chain_complex) function.  This will produce a new chain complex with lots more data
* change of basis matrices
* indices: row/column indices, pivot indices, etc.
* (co)cycle representatives
* barcodes, persistence diagrams, birth/death times, etc.

 # Quick facts

* List of chain complex oracles
    * A list of existing chain complex oracles, including [clique](crate::clique) and [cubical](crate::cubical) complexes, can be found at the bottom of the documentation page for the [`ChainComplex`](crate::chx::ChainComplex) trait.
    * **You can also build your own**!  In concrete terms, this means defining a new [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that [implements](https://doc.rust-lang.org/book/ch10-02-traits.html) the [`ChainComplex`](crate::chx::ChainComplex) trait.



 ```
 Example 1: (to fill in later)
     - create a clique complex
     - access a boundary row
     - access a boundary column
     - get the birth parameter of a simplex
     - factor the complex
     - access the barcode
     - access a cycle representative
     - compute the birth/death time of a chain
 ```

 ```
 Example 2: (to fill in later)
     - create a cubical complex
     - access a boundary row
     - access a boundary column
     - get the birth parameter of a cube
     - factor the complex
     - access the barcode
     - access a cycle representative
     - compute the birth/death time of a chain
 ```



*/


// # Examples

// The intended workflow for a chain complex oracle is as follows

// * **Build a chain complex oracle**
//     * A chain complex oracle is a [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that implements the [`ChainComplex`](crate::chx::ChainComplex) trait.


// * **Explore the chain complex**
//     * Use the commands in the [`ChainComplex`](ChainComplex) trait to access boundary matrices, birth times, and more.

// * **Factor the complex**
//     *   Some information (e.g. persistence diagrams, cycle representatives) won't be available right away.  For this you'll need to *factor* the complex.  You can do this yourself, or you can use the [`factor_chain_complex`](factor_chain_complex) function.  This will produce a new chain complex with lots more data
//         * change of basis matrices
//         * indices: row/column indices, pivot indices, etc.
//         * (co)cycle representatives
//         * barcodes, persistence diagrams, birth/death times, etc.


use std::hash::Hash;
use std::collections::HashMap;
use std::marker::PhantomData;
use std::fmt::Debug;
use std::ops::{Add, Neg, AddAssign, Mul};

use crate::matrix::{MajorDimension, SmOracle, InvMod};
use crate::csm::CSM;
use crate::solver::{multiply_hash_smoracle, triangular_solver_version_2};
use crate::decomp_row::decomp_row;

// ================================================================================================
// FILE DESCRIPTION
// This file contains some outlines for the ChainComplex trait, and its associated objects, functions, etc.
//
// FILE HISTORY
// 20201119: created by Greg Henselman-Petrusek
// 20201119-present: edited by Greg Henselman-Petrusek and Haibin Hang
// ================================================================================================









// -----------------------------------------------------------------------------------------------
// ABBREVIATIONS
// -----------------------------------------------------------------------------------------------


// Chx - chain complex
// Gbf - graded bifiltration
// Subquo - a subquotient (that is, a quotient of a subspace by another subspace)


// -----------------------------------------------------------------------------------------------
// SCRATCHWORK
// -----------------------------------------------------------------------------------------------


// PROPOSED PRINCIPLE: AVOID USING UUID'S WHERE POSSIBLE, ESPECIALLY IF YOU NEED TO BE ABLE TO TELL THAT TWO OBJECTS ARE THE SAME (FOR EXAMPLE, TO SUBQUOTIENTS OF THE SAME VECTOR SPACE)

// CHX SHOULD SPECIFY WHETHER YOU WANT ROW OR COLUMN MAJOR

// SUBQUOTIENTS SHOULD HAVE 2 FILTRATIN PARAMETERS AND 1 DEGREE PARAMETER

// BASIS <-> BASIS INDEX SET <-> MATRIX INDEX SET <-> MATRIX

// SUBQUOID <-> (VECSPACEID, SUBSPACESPEC) // SUBSPACE SPEC IS A LATTICE ELEMENT

// VECSPACE ID <-> UUID



// ------------------------------------------------------------------------------------------------
// FILTRATION PARAMETERS
// ------------------------------------------------------------------------------------------------



// This struct records a filtration parameter for a subspace.  The idea is that
// bot < realnum(s) < realnum(t) < top
// whenever s < t are real numbers -- even when those real numbers are floats representing +/- infinity.
/// unused
pub enum ExtendedLine<Filtration>{
    Top,
    Bot,
    Realnum(Filtration)
}

// This is two number lines side by side, so that
// lower(x) < upper(y)
// for all x and y
/// unused
pub enum ExtendedLineTwice<Filtration>{
    Upper(ExtendedLine<Filtration>), // indexes subspaces in the INVERSE image of the boundary operator
    Lower(ExtendedLine<Filtration>) // indexes subspaces in the image of the boundary operator
}

// -----------------------------------------------------------------------------------------------
// HOW TO SPECIFY A PAIR OF FUNCTIONS OF FORM:
// (INDEX SET) --> (BASIS) --> SUBQUOTIENT
// -----------------------------------------------------------------------------------------------

// MOTIVATION
// ----------
// USUALLY WE DON'T INDEX A MATRIX BY THE VECTORS IN A BASIS B.  INSTEAD, WE FIX A BIJECTION f: I -> B (F0R EXAMPLE, I = [0, ..., N]) AND INDEX THE MATRIX BY I.  IT'S IMPORTANT TO BE PRECISE ABOUT NOT JUST I, BUT THE FUNCTION f ITSELF.

// STRUCTS
// -------
// HOW TO INTERPRET THIS STRUCT
// - Gbf STANDS FOR 'GRADED BIFILTRATION'
// - Subquo STANDS FOR SUBQUOTIENT
// - (m, t, lower(s)) REPRESENTS {m-chains born by time t} \cap \partial(m+1 chains born by s)
// - (m, t, upper(s)) REPRESENTS {m-chains born by time t} \cap \partial^{-1}(m-1 chains born by s)
// - [(a0,b0,c0),...,(an,bn,cn)] REPRESENTS THE SUM OF ALL THE SUBSPACES (ai,bi,ci)
// - **THIS STRUCT REPRESENTS THE QUOTIENT SPACE OF THE numerator BY THE denominator**
// homology in dimension 1 at time 0.5
//
// ```
// let x = GbfSubquoSpec{
//     numerator: (1, ExtendedLine::realnum(0.5), ExtendedLineTwice::upper(ExtendedLine::bot)), //  {space of 1-cycles born by time 0.5}
//     denominator: (1, ExtendedLine::realnum(0.5), ExtendedLineTwice::lower(ExtendedLine::realnum(0.5)))
//     //  {space of 1-boundaries at time 0.5}
// }
// ```
/// unused
pub struct GbfSubquoSpec<Filtration>{
    numerator: Vec<(usize, ExtendedLine<Filtration>, ExtendedLineTwice<Filtration>)>,
    denominator: Vec<(usize, ExtendedLine<Filtration>, ExtendedLineTwice<Filtration>)>,
}

// EXAMPLE: IF WE'RE COMPUTING THE PERSISTENT HOMOLOGY OF A CLIQUE COMPLEX, THEN
// `standard` MEANS THE BASIS OF SIMPLICES
// `matching` MEANS THE MATCHING BASIS COMPUTED BY THE PH COMPUTATION
// `other`    MEANS A CUSTOMIZED BASIS DESIGNED BY THE USER
/// unused
pub enum ChxSubquoBasisSpec{
    Standard,
    Matching,
    Other(String)
}

// THIS STRUCT DOESN'T CONTAIN MUCH REAL DATA; IT'S JUST A `NAMING CONVENTION` THAT LETS THE USER SPECIFY A PARTICULAR BIJECTIVE FUNCTION OF FORM f:I ->B.  The `KeyKind` PARAMETER SPECIFIES THE TYPES OF THE KEYS, SINCE THESE ARE NOT SPECIFIED BY THE UNDERLYING MATHEMATICAL OBJECTS.
// ++++++++++++++++ REMARK: Add "PhantomData" such that "ChxSubquoBasisBijectionSpec<IndexKind, Filtration>" in the future. ++++++++++++++++++
/// unused
pub struct ChxSubquoBasisBijectionSpec<Filtration>{
    subquo: GbfSubquoSpec<Filtration>, // the subquotient
    basis_spec: ChxSubquoBasisSpec, // this has 3 options (standard, matching, other)
    index_spec: String  // this can have lots of variety: simplicial, integer_filtlex, etc.
}

// -----------------------------------------------------------------------------------------------
// THE CHAIN COMPLEX TRAIT: HOW TO OBTAIN
//      MATRICES
//      BASIS INDICES
//      BIRTH/DEATH TIMES
// -----------------------------------------------------------------------------------------------

/// A formal symbol that represents either the identity map or the boundary map.
pub enum ChxTransformKind{
    Identity,
    Boundary
}

/// Stores data about pivot rows/columns and birth times.
pub struct BasisIndexData<Filtration, IndexKind>{
    pub index:                  IndexKind,
    pub birth:                  Option<Filtration>,
    pub homology_degree:        Option<usize>,
    pub matched_index_upper:    Option<IndexKind>,
    pub matched_index_lower:    Option<IndexKind>,
    pub matched_birth_upper:    Option<Filtration>,
    pub matched_birth_lower:    Option<Filtration>
}

//are returned as objects that implement the [SmOracle] trait.  All the boundary matrices have the same type, and use the same type of elements for  major and minor keys.

/**
Methods for based, filtered chain complexes

# Overview

A *based filtered chain complex* is a filtered chain complex *C* together with a basis *E(n)* for the chains in each dimension *n*.  We assume that *E(n)* is indexed by a bijective function *f:I(n) -> E(n)*, and refer to basis vectors by their indices.

There are several types of information you might want from a based chain complex:

* **Boundary matrices** The *n*-dimensional boundary matrix has rows and columns indexed by *I(n-1)* and *I(n)*, respectively.

* **Row and column indices** The index sets *I(n)* and the pivot pairs of each matrix factorization.  We rely on the chain complex oracle to supply this information because [SmOracle](SmOracle) trait objects don't always contain complete information about index sets.  **Sets *I(n)* and *I(m)* should not intersect when *m* and *n* are different.**

* **Safe dimensions** Sometimes boundary matrices are so big they can't fit in memory.  It's usually best to avoid building these.

* **Birth parameters** The "birth time" of each basis vector.

* **Persistence information** meaning persistence diagrams, (co)cycle representatives, etc.
    * You generally need to "factor" [INSERT REFERENCE TO MATCHING BASES] a chain complex in order to get this type of data.  There is a [`factorization function`](chx/function.factor_chain_complex) to perform this function; the output is stored in a  [`FactoredComplexBlockCsm`](FactoredComplexBlockCsm) object.


The [ChainComplex] trait gives a general set of conventions for extracting all this and more.


 ```
 Example 1: (to fill in later)
     - create a clique complex
     - access a boundary row
     - access a boundary column
     - get the birth parameter of a simplex
     - factor the complex
     - access the barcode
     - access a cycle representative
     - compute the birth/death time of a chain
 ```

 ```
 Example 2: (to fill in later)
     - create a cubical complex
     - access a boundary row
     - access a boundary column
     - get the birth parameter of a cube
     - factor the complex
     - access the barcode
     - access a cycle representative
     - compute the birth/death time of a chain
 ```


# Type parameters (different from associated types)

* `MatrixIndexKey` the type of the row and column indices.
* `SnzVal` the type of the structural nonzero elements of the coefficient ring (different rings can have elements of the same type).
* `Filtration` the type of the filtration values; typically these will be either floats or integers.






*/
pub trait ChainComplex<MatrixIndexKey, SnzVal, Filtration> where
MatrixIndexKey: PartialEq + Eq + Hash + Clone,
SnzVal: Clone
{

    /// Every matrix returned by the [ChainComplex] object will have this type.
    type Matrix: SmOracle<MatrixIndexKey, MatrixIndexKey, SnzVal>;

    /// Returns the sparse matrix oracle which representing the boundary operators
    fn get_smoracle(
        &self,
        major_dim       :   MajorDimension, // row or column
        transform       :   ChxTransformKind
        //basis_index_spec_maj: ChxSubquoBasisBijectionSpec<Filtration>,
        //basis_index_spec_min: ChxSubquoBasisBijectionSpec<Filtration>
    ) -> Self::Matrix;

    /// Return iterator of keys (unordered)
    fn keys_unordered_itr(&self, dim: usize) -> Box<dyn Iterator<Item=MatrixIndexKey> + '_>;

    /// Return all keys(simplices) of given dimension with increasing filtration value order
    fn keys_ordered(&self, h_degree: usize) -> Vec<MatrixIndexKey>;

    /// Return the filtration value of given simplex(key)
    fn key_2_filtration(&self, key: &MatrixIndexKey) -> Filtration;

    /// Return the maximum filtration value (cutoff value)
    fn max_filtration(&self) -> Filtration;
/*
    // This is the principle means of associating filtration parameters and grades (ie dimensions) to each (element that indexes a) basis vector
    // This should return an error message if
    //   - the user asks to return death times of standard basis vectors
    //   - the user tries to get a basis of a "death" subspace (or a subspace of a death subspace)
    fn get_basis_index_iterator<IndexKind>(
        basisindexspec: ChxSubquoBasisBijectionSpec<Filtration>,
        return_birthtime: bool,
        return_deathtime: bool,
        return_grade:     bool
    ) -> BasisIndexData<Filtration, IndexKind>;

    /// Returns the dimension of the subquotient vector space.  This function may return a None value if the dimensions of one or more chain spaces have not been computed (e.g., if we have not yet counted all the simplices of a simplicial complex).  In addition, it may be expensive or slow to calculate the rank (for example, if this requires enumerating a large number of simplices).  Setting the `not_safe` keywork argument to false should prevent this calculation from running and pring a message explaining why the calculation did not run.
    fn rank( &self, subquo_spec: GbfSubquoSpec<Filtration>, not_safe: bool) -> Option<usize>;

    /// Returns a sorted vector containing the chain degrees where the complex will return non-None values.
    fn safe_homology_degrees_to_build_boundary( &self ) -> Vec<usize>;

    /// Returns true if birth times have been computed for the matching basis in this degree.
    fn filtration_parameters_available( &self,
        matched_codimension: usize,
        basisindexspec: ChxSubquoBasisBijectionSpec<Filtration>
    ) -> bool;

    fn get_filtration_parameter<IndexKeyKind>( &self,
        basisindexspec: ChxSubquoBasisBijectionSpec<Filtration>,
        index: MatrixIndexKey
    ) -> ExtendedLine<Filtration>;
*/
}

/// Records some bijections between intergers and keys. Pivot pairs are indexed by the same usize integer. Particularly, major keys are indexed in decreasing filtration value order
pub struct Indexing<MinKey, MajKey> {
    pub minkey_2_index: HashMap<MinKey, usize>,
    pub majkey_2_index: HashMap<MajKey, usize>,
    pub index_2_minkey: Vec<MinKey>,
    pub index_2_majkey: Vec<MajKey>,
    pub ordered_minind: Vec<usize>    //increasing filtration value order
}
/// Methods of Indexing struct
impl<MinKey: Hash + Eq, MajKey: Hash + Eq> Indexing<MinKey, MajKey>{
    /// Construct a trivial Indexing instance
    pub fn new() -> Indexing<MinKey, MajKey> {
        Indexing {
            majkey_2_index: HashMap::new(),
            minkey_2_index: HashMap::new(),
            index_2_majkey: Vec::new(),
            index_2_minkey: Vec::new(),
            ordered_minind: Vec::new()   // increasing order
        }
    }
    /// Construct a trivial Indexing instance with given capacity
    ///
    /// # Parameters
    /// - `capacity`: the given capacity
    pub fn with_capacity( capacity: usize ) -> Indexing<MinKey, MajKey> {
        Indexing {
            majkey_2_index: HashMap::with_capacity(capacity),
            minkey_2_index: HashMap::with_capacity(capacity),
            index_2_majkey: Vec::with_capacity(capacity),
            index_2_minkey: Vec::with_capacity(capacity),
            ordered_minind: Vec::with_capacity(capacity),   // increasing order
        }
    }

    /// Shrink the capacity of the Indexing to save memory
    pub fn shrink_to_fit( &mut self ) {
        self.majkey_2_index.shrink_to_fit();
        self.minkey_2_index.shrink_to_fit();
        self.index_2_majkey.shrink_to_fit();
        self.index_2_minkey.shrink_to_fit();
        self.ordered_minind.shrink_to_fit();
    }

}

/// A factored chain complex with change-of-basis matrices
///
/// See also [Factoring](crate::chx) in the documentation for the `chx` module.
pub struct FactoredComplexBlockCsm<'a, MatrixIndexKey, SnzVal, Filtration, OriginalChx> where
MatrixIndexKey: PartialEq + Eq + Hash + Clone,
SnzVal: Clone,
OriginalChx: ChainComplex<MatrixIndexKey, SnzVal, Filtration>
{
    pub phantom: PhantomData<Filtration>,
    pub original_complex: &'a OriginalChx,   // Reference to the original complex
    pub dim_rowoper: Vec<CSM<usize, SnzVal>>,
    pub dim_indexing: Vec<Indexing<MatrixIndexKey, MatrixIndexKey>>,
}
/// Methods of FactoredComplexBlockCsm struct
impl<'a, MatrixIndexKey, SnzVal, Filtration, OriginalChx> FactoredComplexBlockCsm<'a, MatrixIndexKey, SnzVal, Filtration, OriginalChx> where
Filtration: PartialOrd + Clone,
MatrixIndexKey: PartialEq + Eq + Ord + Hash + Clone + Debug,
SnzVal: Clone + PartialEq + Neg<Output=SnzVal> + AddAssign + Mul<Output = SnzVal> + InvMod<Output = SnzVal> + Debug,
OriginalChx: ChainComplex<MatrixIndexKey, SnzVal, Filtration>
{
    /// Returns the persistence barcode of given dimension
    pub fn barcode(&self, h_degree: usize) -> Vec<(Filtration, Filtration)> {

        let mut barcode = Vec::new();
        if h_degree >= self.dim_rowoper.len() { return barcode; }
        //let keys = self.original_complex.keys_ordered(h_degree);
        let indexing = &self.dim_indexing[h_degree];
        let mut indexing_upper: &Indexing<MatrixIndexKey, MatrixIndexKey> = &Indexing::new();
        if h_degree+1 < self.dim_rowoper.len() {
            indexing_upper = &self.dim_indexing[h_degree+1];
        }

        for key in self.original_complex.keys_unordered_itr(h_degree) {
            if indexing.minkey_2_index.contains_key(&key) { continue; }
            else if indexing_upper.majkey_2_index.contains_key(&key) {
                let ind = indexing_upper.majkey_2_index[&key];
                let matched_key = &indexing_upper.index_2_minkey[ind];
                let diam1 = self.original_complex.key_2_filtration(&key);
                let diam2 = self.original_complex.key_2_filtration(matched_key);
                if diam1 == diam2 { continue; }
                barcode.push((diam1, diam2));
            } else {
                let diam1 = self.original_complex.key_2_filtration(&key);
                let diam2 = self.original_complex.max_filtration();
                if diam1 == diam2 { continue; }
                barcode.push((diam1, diam2));
            }
        }
        return barcode;
    }

    /// Return matched basis vector of given original_basis at given dimension
    pub fn get_matched_basis_vector(
        &self,
        h_degree:           usize,
        original_basis:     &MatrixIndexKey
    ) -> HashMap<MatrixIndexKey, SnzVal>
    {

        let matrix = self.original_complex.get_smoracle(
            MajorDimension::Row,
            ChxTransformKind::Boundary
        );

        if self.dim_indexing[h_degree+1].majkey_2_index.contains_key(original_basis) {
            let indexing = &self.dim_indexing[h_degree+1];
            let index = indexing.majkey_2_index[original_basis];
            let minkey = &indexing.index_2_minkey[index];

            let order = &(0..indexing.index_2_majkey.len()).collect();
            let rowoper = self.dim_rowoper[h_degree+1].change_major_dim(order);
            let mut reduced = CSM::new(MajorDimension::Col, matrix.ring().clone());
            let mut boundary = CSM::new(MajorDimension::Col, matrix.ring().clone());

            let mut col_boundary = HashMap::new();
            let mut col_reduced = HashMap::new();
            let mut col_inverse = HashMap::new();
            let mut bb = HashMap::new();
            for minind in indexing.ordered_minind.iter() {
                col_boundary.clear();
                let minor_key = &indexing.index_2_minkey[*minind];
                for (key, val) in matrix.min_itr(&minor_key) {
                    boundary.push_snzval(key.clone(), val.clone());
                    if indexing.majkey_2_index.contains_key(&key) {
                        col_boundary.insert(indexing.majkey_2_index[&key], val);
                    }
                }
                boundary.majptr.push(boundary.minind.len());
                boundary.nummaj += 1;

                col_reduced.clear();
                col_reduced = multiply_hash_smoracle(&col_boundary, &rowoper);
                reduced.append_maj(&mut col_reduced);

                if minor_key == minkey {
                    bb.insert(*minind, matrix.ring().identity_multiplicative.clone());
                    break;
                }
            }

            col_inverse.clear();
            col_inverse = triangular_solver_version_2(&reduced, &indexing.ordered_minind, &mut bb);
            return multiply_hash_smoracle(&col_inverse, &boundary);
        } else if self.dim_indexing[h_degree].minkey_2_index.contains_key(original_basis){
            let indexing = &self.dim_indexing[h_degree];
            let order = &(0..indexing.index_2_majkey.len()).collect();
            let rowoper = self.dim_rowoper[h_degree].change_major_dim(order);
            let mut reduced = CSM::new(MajorDimension::Col, matrix.ring().clone());

            let mut col_boundary = HashMap::new();
            let mut col_reduced = HashMap::new();
            let mut col_inverse = HashMap::new();
            let mut bb = HashMap::new();
            for minind in indexing.ordered_minind.iter() {
                col_boundary.clear();
                let minor_key = &indexing.index_2_minkey[*minind];
                for (key, val) in matrix.min_itr(minor_key) {
                    if indexing.majkey_2_index.contains_key(&key) {
                        col_boundary.insert(indexing.majkey_2_index[&key], val);
                    }
                }
                col_reduced.clear();
                col_reduced = multiply_hash_smoracle(&col_boundary, &rowoper);
                reduced.append_maj(&mut col_reduced);

                if minor_key == original_basis {
                    bb.insert(*minind, matrix.ring().identity_multiplicative.clone());
                    break;
                }
            }

            col_inverse.clear();
            col_inverse = triangular_solver_version_2(&reduced, &indexing.ordered_minind, &mut bb);

            let mut new_basis = HashMap::new();
            for (ind, val) in col_inverse.drain(){
                let min_ind = indexing.ordered_minind[ind];
                new_basis.insert(indexing.index_2_minkey[min_ind].clone(), val);
            }
            return new_basis;
        } else {
            let indexing = &self.dim_indexing[h_degree];
            let order = &(0..indexing.index_2_majkey.len()).collect();
            let rowoper = self.dim_rowoper[h_degree].change_major_dim(order);
            let mut reduced = CSM::new(MajorDimension::Col, matrix.ring().clone());

            let mut col_boundary = HashMap::new();
            let mut col_reduced = HashMap::new();
            for minind in indexing.ordered_minind.iter() {
                col_boundary.clear();
                let minor_key = &indexing.index_2_minkey[*minind];
                for (key, val) in matrix.min_itr(minor_key) {
                    if indexing.majkey_2_index.contains_key(&key) {
                        col_boundary.insert(indexing.majkey_2_index[&key], val);
                    }
                }
                col_reduced.clear();
                col_reduced = multiply_hash_smoracle(&col_boundary, &rowoper);
                reduced.append_maj(&mut col_reduced);
            }

            col_boundary.clear();
            for (key, val) in matrix.min_itr(original_basis) {
                if indexing.majkey_2_index.contains_key(&key) {
                    col_boundary.insert(indexing.majkey_2_index[&key], val);
                }
            }

            let mut col_inverse = triangular_solver_version_2(&reduced, &indexing.ordered_minind, &mut col_boundary);
            let mut new_basis = HashMap::new();
            for (ind, val) in col_inverse.drain(){
                let min_ind = indexing.ordered_minind[ind];
                new_basis.insert(indexing.index_2_minkey[min_ind].clone(), val);
            }
            let mut new_basis = HashMap::new();
            new_basis.insert(original_basis.clone(), matrix.ring().identity_multiplicative.clone());
            return new_basis;
        }
    }

}


/// Returns a [factored chain complex](FactoredComplexBlockCsm), giving access to barcodes, (co)cycle representatives, etc.
///
/// # Parameters
/// - `original_complex`: The chaincomplex instance based on which we want to compute barcodes
/// - `max_homology_degree`: The max degree of which we want to compute barcodes
/// # Returns
/// an FactoredComplexBlockCsm struct which records decomposition information of each dimension
/// # See also
/// [Factoring](crate::chx) in the documentation for the `chx` module.
pub fn factor_chain_complex<'a, MatrixIndexKey, SnzVal, Filtration, OriginalChx>(
    original_complex:       &'a OriginalChx,
    max_homology_degree:    usize
) -> FactoredComplexBlockCsm<'a, MatrixIndexKey, SnzVal, Filtration, OriginalChx> where
MatrixIndexKey: PartialEq + Eq + Ord + Hash + Clone + Debug,
SnzVal: Clone + PartialEq + Neg<Output=SnzVal> + AddAssign + Mul<Output = SnzVal> + InvMod<Output = SnzVal> + Debug + Add,
Filtration: Debug + PartialOrd,
OriginalChx: ChainComplex<MatrixIndexKey, SnzVal, Filtration>
{
    let matrix = original_complex.get_smoracle(
        MajorDimension::Row,
        ChxTransformKind::Boundary
    );

    let mut blocks = FactoredComplexBlockCsm{
		phantom: PhantomData,
		original_complex: original_complex,
		dim_rowoper: vec![CSM::new(MajorDimension::Row, matrix.ring().clone())],
		dim_indexing: vec![Indexing::new()]
	};

	let mut maj_to_reduce = Vec::new();
	for ii in 1..(max_homology_degree+1){
		maj_to_reduce.clear();
		let major_keys = original_complex.keys_ordered(ii-1);
        for key in major_keys.iter() {
            if !blocks.dim_indexing[ii-1].minkey_2_index.contains_key(key) {
				maj_to_reduce.push(key.clone());
			}
        }
        //println!("length of maj_to_reduce {}", maj_to_reduce.len());
		let (rowoper, indexing) = decomp_row(&matrix, &mut maj_to_reduce);
        //println!("{}", indexing.index_2_majkey.len());
		blocks.dim_rowoper.push(rowoper);
		blocks.dim_indexing.push(indexing);
	}

    return blocks;
}
