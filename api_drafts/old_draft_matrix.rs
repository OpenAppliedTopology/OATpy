// ================================================================================================
// FILE DESCRIPTION
// This file contains some outlines matrix (oracle) objects and traits.
//
// FILE HISTORY
// 20201120: created by Greg Henselman-Petrusek
// 20201120-present: edited by Greg Henselman-Petrusek and Haibin Hang
// ================================================================================================


// ------------------------------------------------------------------------------------------------
// OVERVIEW OF CONTENTS
// ------------------------------------------------------------------------------------------------


// This code is a general outline for how to interact with two kinds of objects: matrix oracles and row oracles.
//
//
// Matrix oracles
// --------------
//
// Matrix oracles are integral to the draft API.  There are a lot of details to consider vis-a-vis matrix oracle API's, including:
// - Indexing
//      - In general, a matrix is a function M: I x J -> CoefficientRing.  It's sometimes handy to let I and J be something other than sets of integers.  So we need some general system of dealing with different types of indices.
// - Matrix representations
//      - For M to represent a lienar map T, you generally need a bijections
//              i: I -> (basis for the codomain of T)
//              j: J -> (basis for the domain of T)
//        Thus, for based matrices, we need some way to specify the bases AND the bijections
//
//
// Row oracles
// -----------
//
// A row oracle is like a matrix oracle, except that
// - there are no 'column indices' associate with a row oracle (it just retruns a collection of values for each row)
// - you can ask a row oracle for multiple types of data
// So, for example, you might define a row oracle that eats a simplex S and returns a list of any one fo the following things: 1. the cofacets of S, 2. the diameters fo the cofacets of S, 3. pairs of form (cofact, diameter of cofacet)
//


// -----------------------------------------------------------------------------------------------
// (IGNORE IF YOU ARE A FIRST TIME READER) THINGS TO CONSIDER, VIS-A-VIS MATRICES
// -----------------------------------------------------------------------------------------------


// things to consider:
//     - possible formats for each oracle (resp, each row): array, iterator, hashmap, function
//         - if hashmap or function
//             - what are the keys?
//         - if array or iterator
//             - questions about order
//                 - how is order determined?
//                 - do we require sorted (with respect to filtration) order?
//                     - in the priority queue implementation, one never actually places a sequence of entries in sorted order (just puts everything in a priority queue)
//                     - however, the process that generates the priority queue *does* follow a linear order
//                     - however, parallelized implementations may not preserve this order
//                 - do we allow user to choose the order?
//                 - do we require elements to appear in filtration order?
//                 - what if someone wants elements in opposite of filtration order? (which is common)
//                     - should we havea filtration object with an attribute to flag whether its ordered?
//             - what if someone wants an iterator that generates filtration values in the form of (key, value) pairs, eg (simplex, simplex_diameter)?
//     - should there be a filtration object?
//         - in part, this addresses/is addressed by the BasisIndexSpec struct


// -----------------------------------------------------------------------------------------------
// ROW AND MATRIX ORACLES
// -----------------------------------------------------------------------------------------------


/// This trait formalizes what we want from a "row oracle," something that takes a key as input, and returns a row (not necessarily of a matrix ... possibly of a sparse table) as output.  Organizing data into a family of rows is quite common in PH computation (CSR matrix format, boundary matrix rows stored as heaps, etc.), and this just formalizes this organization.  NOTA BENE: this trait is distinguished by having a well-defined type for the keys that index rows, but no defined type for column indices.  Thus it makes no sense to multiply two of these objects together, as in matrix multiplication.
pub enum RowTypes{
    Array,
    Hashmap,
    Vector
}

pub enum KeyTypes{
    Usize,
    Simplex
}

pub trait RowOracle{

    /// Return the type of the keys that index rows.  Note that there is no corresponding colkey_type, since we do not assume anything about the columns of the table.
    pub fn type_rowkey( &self ) -> KeyTypes;

    /// Retreive a row
    ///
    /// # Arguments
    /// * key: an instance of self.rowkey_type()
    /// * dataID: specifies which of several different rows to return; for example, one might ask EITHER for an array of column pointers OR an array of structural nonzero values
    /// * row_type: determines the type of the output `row`.  Allowable types include: Array (user required to implement this), Iterator (user required to implement this, but defaults to iterating across the array), Hash (optional), Function (optional)
    pub fn row( &self, key: String, dataID: String, rowtype: String ) -> RowTypes;

    /// !!! I'm not sure how necessary/important this function is.  We might want to delete it.
    /// Return a row-slice (either occupying new memory or just masking) of the row oracle
    ///
    /// # Attributes
    /// * rows: the set of rows to slice
    /// * mask: if True, then return a masked array

    //pub fn rowslice( &self, rows: [rowkey_type], mask: bool ) -> Self;

}

/// Some useful extensions to the `RowOracle` trait, mostly to do with defining/searching the set of all valid row index keys.
pub trait RowOracle_exhaustive: RowOracle{

    /// Returns an iterator that runs over all rows in the RowOracle, ehaustively.
    ///
    /// # Arguments
    /// * dataID: specifies which data to produce
    /// * rowtype: specifies the type of each row
    /// * orderspec: specifies the order in which to run over rows
    ///
    /// # Returns
    /// * An iterator `I` that returns a tuple of form (rowindex_key, row) on each iteration. NOTA BENE: (i) Even if `type(row_index_key) == usize`, we do not require that `row_index_key = p` on the `p`th iteration. (ii) The function {iteration's} -> {rowindex_key's} should be injective. (iii) In general, there may be several different order specs that may be better or worse suited to specific applications (eg lexicographic versus filtration-with-ties-broken-by-lexicographic).
    pub fn iterator_keyval( &self, dataID: String, rowtype: Type, orderspec: Option<String> ) -> Iterator;

    /// Same as iterator_keyval, but returns only rowindex_key's
    pub fn iterator_key( &self, dataID: String, orderspec: Option<String> ) -> Iterator ;

    /// Returns True if and only if key is a valid row index key.
    pub fn is_row_index_key<T>( &self, key: T ) -> bool;

    /// Returns [a,b] if and only if the set of rowindex_key's is [a, .., b]
    pub fn rowindex_minmax( &self ) -> Option<[usize]>;

    /// Returns True if and only if the row index of the `p`th iteration is `p`
    pub fn iterator_index_matches_row_index( & self, dataID: String, orderspec: Option<String> ) -> bool;
}

// A general RowOracle format that imitates CSR format for numerical matrices.  Similar to traditional CSR, this oracle lets one extract column pointers and entry coefficients for a specific row of the matrix.  The distinguishing features of this type, as compared with the RowOracle, are: (i) it has a well-defined type for column index keys, (ii) it has a well-defined matrix entry type, (iii) it has some special calls for specific rows.  NOTA BENE: we do note require this to implement the `RowOracle_exhaustive` trait.
pub trait RowOracle_PseudoCSR: RowOracle{

    type ColKey;
    type SnzVal;
    /// Returns the type of the column index keys. (MAYBE WE SHOULD CHANGE THIS NAME TO MATCH ROWS)
    pub fn type_colkey( &self ) -> Self::ColKey;

    /// Returns the type of the structural nonzero entries.
    pub fn type_snzval( &self ) -> Self::SnzVal;

    /// Return the column pointers for row `rowindex_key`.  Format this collection of pointers as `rowtype`.
    pub fn colind<S, T where T==self.type_rowkey()>( &self, rowindex_key: T, rowtype: S ) -> T;

    /// Return the row of structural nonzeros indexed by `rowindex_key` formatted as `rowtype`
    pub fn snzval<S, T where T==self.type_rowkey()>( &self, rowindex_key: T, rowtype: S ) -> T;

    /// Return the sequence of (colindex_key, snzval) tuples for this row, formatted as `rowtype`.  This command can be useful for PH computations, where it takes essentially the same amount of work to generate the column indices, the snz values, or both (simultaneously).
    pub fn index_value_pairs<S, T where T==self.rowkey_type()>( &self, rowindex_key: T, rowtype: S ) -> T;

    /// Return entry [i,j]
    pub fn entry<S,T,U>( &self, rowind: S, colind: T) -> U;
}


/// A unique id for a specific bijection between a basis and a set of index keys.  In general, the bijection itself can't be derived from the data in this object; it's just a 'name.'
struct BasisIndexSpec{
    basisID: BasisID,
    indexspec: String
}

/// This trait gives complete data about a matrix representation of linear map stored on computer.  In particular, suppose
///     T: a linear map
///     B: a basis for the domain
///     C: a basis for the codomain
///     X: a bijection I -> B
///     Y: a biJection J -> C
///     M: an I x J matrix such that M[i,j] is the proper coefficient of the matrix rep of T with respect to B and C
/// This trait provides 'names' for X and Y (in the forrm of `BasisIndexSpec`s) and an oracle to evaluate the rows of M.

struct Based_CSR{
    csr: CmpRowFmt,
    ind_of_rows: BasisIndexSpec,
    ind_of_cols: BasisIndexSpec,
}

trait BasedMatrix_RowOracle_PseudoCSR{

    /// Returns the underlying matrix oracle.
    fn oracle( &self ) -> &impl RowOracle_PseudoCSR;

    /// Specifies the function that indexes matrix rows.
    fn basisindexspec_row( &self ) -> BasisIndexSpec;

    /// Specifies the function that indexes matrix cols.
    fn basisindexspec_col( &self ) -> BasisIndexSpec;

    // FUNCTIONS TO IMPLEMENT FOR THIS TYPE:
    // - product, sum, scale, rowslice, colslice
    // - each oper should be implemented with (i) lazy/masked and/or (ii) internal CSR output format
    //     - if internal CSR, then lefthand oracle must have exhasutive rows
}
