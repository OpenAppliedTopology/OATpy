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


/// This struct records a filtration parameter for a subspace.  The idea is that
/// bot < realnum(s) < realnum(t) < top
/// whenever s < t are real numbers -- even when those real numbers are floats representing +/- infinity.
enum ExtendedLine<T>{
    top,
    bot,
    realnum(T)
}

/// This is two number lines side by side, so that
/// lower(x) < upper(y)
/// for all x and y
enum ExtendedLineTwice<T>{
    upper(ExtendedLine<T>), // indexes subspaces in the INVERSE image of the boundary operator
    lower(ExtendedLine<T>) // indexes subspaces in the image of the boundary operator
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
struct GbfSubquoSpec<Filtration>{
    numerator: Vec<(usize, ExtendedLine<Filtration>, ExtendedLineTwice<Filtration>)>,
    denominator: Vec<(usize, ExtendedLine<Filtration>, ExtendedLineTwice<Filtration>)>,
}

// EXAMPLE: IF WE'RE COMPUTING THE PERSISTENT HOMOLOGY OF A CLIQUE COMPLEX, THEN
// `standard` MEANS THE BASIS OF SIMPLICES
// `matching` MEANS THE MATCHING BASIS COMPUTED BY THE PH COMPUTATION
// `other`    MEANS A CUSTOMIZED BASIS DESIGNED BY THE USER
enum ChxSubquoBasisSpec{
    standard,
    matching,
    other(String)
}

// THIS STRUCT DOESN'T CONTAIN MUCH REAL DATA; IT'S JUST A `NAMING CONVENTION` THAT LETS THE USER SPECIFY A PARTICULAR BIJECTIVE FUNCTION OF FORM f:I ->B.  The `KeyKind` PARAMETER SPECIFIES THE TYPES OF THE KEYS, SINCE THESE ARE NOT SPECIFIED BY THE UNDERLYING MATHEMATICAL OBJECTS.
struct ChxSubquoBasisBijectionSpec<Filtration, IndexKind>{
    subquo: GbfSubquoSpec<Filtration>, // the subquotient
    basis_spec: ChxSubquoBasisSpec, // this has 3 options (standard, matching, other)
    index_spec: String,  // this can have lots of variety: simplicial, integer_filtlex, etc.
    phantom_indexkind: PhantomData<IndexKind>
}





// -----------------------------------------------------------------------------------------------
// THE CHAIN COMPLEX TRAIT: HOW TO OBTAIN
//      MATRICES
//      BASIS INDICES
//      BIRTH/DEATH TIMES
// -----------------------------------------------------------------------------------------------

// THIS ENUM LETS US SPECIFY WHETHER WE WANT A ROW-MAJOR MATRIX OR A COL-MAJOR MATRIX
enum MajorDimension{
    row,
    col
}

// DOES THE MATRIX REPRESENT THE IDENTITY FUNCTION OR THE BOUNDARY OPERATOR?
enum ChxTransformKind{
    identity,
    boundary
}


struct BasisIndexData<Filtration, IndexKind>{
    index:                  IndexKind,
    birth:                  Option(Filtration),
    homology_degree:        Option(HomologyGrade)
    matched_index_upper:    Option(IndexKind),
    matched_index_lower:    Option(IndexKind),
    matched_birth_upper:    Option(Filtration),
    matched_birth_lower:    Option(Filtration)
}


/// A trait that summarizes everything we want to do with a filtered chain complex
pub trait ChainComplex<StandardBasisIndexKey, SnzVal, Filtration>{
    fn get_smoracle<MinKey, MajKey, SnzVal, T: impl smoracle< MajKey, MinKey, SnzVal >>(
        &self,
        major_dim           :   MajorDimension, // row or column
        transform           :   ChxTransformKind,
        basis_index_spec_maj:   ChxSubquoBasisBijectionSpec< Filtration, MajKey >,
        basis_index_spec_min:   ChxSubquoBasisBijectionSpec< Filtration, MinKey >
    ) -> T;

    // This is the principle means of associating filtration parameters and grades (ie dimensions) to each (element that indexes a) basis vector
    // This should return an error message if
    //   - the user asks to return death times of standard basis vectors
    //   - the user tries to get a basis of a "death" subspace (or a subspace of a death subspace)
    fn get_basis_index_iterator<Filtration,
                                IndexKind,
                                BasisIndexIterator: Iterator<Item = BasisIndexData<Filtration, IndexKind>>
                                >
                                (
        basisindexspec:                ChxSubquoBasisBijectionSpec< Filtration, IndexKind >,
        return_birth:                  bool,
        return_homology_degree:        bool,
        return_matched_index_upper:    bool,
        return_matched_index_lower:    bool,
        return_matched_birth_upper:    bool,
        return_matched_birth_lower:    bool
    ) -> Option(BasisIndexIterator);

    /// Returns the dimension of the subquotient vector space.  This function may return a None value if the dimensions of one or more chain spaces have not been computed (e.g., if we have not yet counted all the simplices of a simplicial complex).  In addition, it may be expensive or slow to calculate the rank (for example, if this requires enumerating a large number of simplices).  Setting the `unsafe` keywork argument to false should prevent this calculation from running and pring a message explaining why the calculation did not run.
    fn rank( &self, subquo_spec: GbfSubquoSpec, unsafe: bool) -> Option<usize>;

    /// Returns a sorted vector containing the chain degrees where the complex will return non-None values.
    fn safe_homology_degrees_to_build_boundary( &self ) -> Vec<usize>;

    /// Returns true if birth times have been computed for the matching basis in this degree.
    fn filtration_parameters_available<IndexKey>(   &self,
                                                    matched_codimension: usize,
                                                    basis_index_spec:   ChxSubquoBasisBijectionSpec
                                                                          <Filtration, IndexKey>
                                ) -> bool;

    fn get_filtration_parameter<IndexKeyKind>(  basisindexspec: ChxSubquoBasisBijectionSpec
                                                                    < Filtration, IndexKeyKind>,
                                                index:          IndexKeyKind
                                            ) -> InternalFiltration{
    };
}

// -----------------------------------
// IMPORTANT CHANGES
// -----------------------------------
//
// - Indices for boundary matrices are TUPLES (simpliex, diameter)

// -----------------------------------
// THINGS WE MAY NEED
// -----------------------------------
//
// - TRANSPOSE FUNCTION
// - TRIANGULAR SOLVE (OR AT LEAST A MODIFICATION OF THE REDUCTION ALGORITHM THAT ALLOWS ONE TO REDUCE A SPARSE VECTOR THAT'S NOT IN THE MATRIX)
// - SPARSE VECTOR - SPARSE VECTOR MULTIPLICATION




pub struct CliqueComplex<MatrixIndexKey, SnzVal: Clone, FiltrationParamKind> {
    pub dissimilarity_matrix: Vec<FiltrationParamKind>,
    pub dissimilarity_value_max: FiltrationParamKind,
    pub safe_homology_degrees_to_build_boundaries: Vec<usize>,
    pub major_dimension: MajorDimension,
    pub ringmetadata: RingMetadata<SnzVal>,
    pub simplex_count: Vec<(usize,usize)>
}


pub struct FactoredChainComplexBlockCsm<    OriginalComplex : ChainComplex< MatrixIndexKey,
                                                                            SnzVal,
                                                                            Filtration>,
                                            MatrixIndexKey,
                                            SnzVal,
                                            Filtration>{
    //  Reference to the original complex
    pub original_complex: OriginalComplex,
    //  Maps converting integers to indices
    pub integer_to_piv_pair: Vec<(usize, Vec<(MatrixIndexKey, MatrixIndexKey)>)>,
    pub unpaired_indices: Vec<(usize, Vec<MatrixIndexKey>)>, // a list of the unpaired simplices
    //  Maps converting indices to integers
    pub hashmap_pivmin_to_int: Vec<(usize, Hash<MatrixIndexKey, usize>)>,
    pub hashmap_pivmaj_to_int: Vec<(usize, Hash<MatrixIndexKey, usize>)>,
    //  Square invertible matrices
    pub invertible_blocks_low: Vec<(usize,CSM<usize, SnzVal>)>,
        // Interpretation:
        // invertible_blocks_hig = RowOprMatrix[pivrows, pivrows]
    pub invertible_blocks_low_inverse: Vec<(usize,CSM<usize, SnzVal>)>,
    pub invertible_blocks_hig: Vec<(usize,CSM<usize, SnzVal>)>,
    pub invertible_blocks_hig_inverse: Vec<(usize,CSM<usize, SnzVal>)>
        // Interpretation:
        // invertible_blocks_hig_inverse = RowOprMatrix[pivrows, pivrows] * Dn[pivrows, pivcols]
}

impl FactoredComplexBlockCsm{

fn factor_chain_complex(&self, max_homology_degree:usize)) -> FactoredComplexBlockCsm{
    // initialize FactoredComplexBlockCsm with empty vectors of pivot pairs and change of basis blocks
    // for dim = 0 ... m, do
    //    factor matrix dn
    //    store the change of basis matrices and pivot pairs in the FactoredChainComplexBlockCsm (we probably want to store the transposes of those matrices, in fact)
    //    before storing these matrices, we delete everything except for the pivot row and cols
}
//  (LD=PR)
//   D(R^{-1}) R (basis)
fn get_matching_basis(&self, homology_degree: usize) -> CSM< MatrixIndexKey, MatrixIndexKey, SnzVal >{
    // IF THE SYSTEM HAS ALREADY CALCULATED THE FOLLOWING MATRICES
    //      Dn = boundary operator (n chains) -> (n-1 chains)
    //      Ln = RowOperMatrix // indexed by (integer PivotRows) x (integer PivotRows)
    //         = X, where (n, X) is in invertible_blocks_low
    //      Rn = Ln * Dn[PivotRows, PivotCols]
    //         = Y, where (n, Y) is in invertible_blocks_hig_inverse
    // THEN THIS IS ONE WAY TO GET THE VECTORS OF THE MATCHING BASIS:
    //      INPUT: an n-simplex named <simplex>:
    //         [[[[[[ if <simplex> is a pivot col of Dn:
    //              -convert <simplex> to an integer m via the hashmap named hashmap_pivmin_to_int
    //              -let <col>:= (Rn^{-1})[:,m]
    //              -convert the (integer) indices of <col> to simplices
    //              -return this sparse vector ]]]]]]
    //          elseif <simplex> is a pivot row of D(n+1):
    //              -convert <simplex> to an integer m via the hashmap named hashmap_pivmaj_to_int
    //              -let <col> := (L(n+1))^{-1}[:, m]
    //              -convert the integer indices of <col> to n-simplices, then multiply the
    //                  vector with Dn to get a sparse vector <vec> representing a linear
    //                  combination of (n-1)-simplices
    //              -if this vector has a nonzero coefficient in a row labeled <n-1 simplex> such //                  that <n-1 simplex> is not a pivot row of Dn, then delete this entry from
    //                  the sparse vector
    //              - solve Ux = <vec> for x, where U is the square invertible matrix equal to
    //                  RowOperMatrix[pivotrows,pivotrows] * Dn[pivotrow, pivotcols]
    //              - return <vec> + x
    //          otherwise:
    //              - follow the instructions for the preceding elseif statement, except delete
    //                  the first two steps and instead let <vec> = Dn[:, <simplex>]
    )

}
