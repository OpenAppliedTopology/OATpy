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


// Matrix oracles are integral to the draft API.  There are a lot of details to consider vis-a-vis matrix oracle API's.  These include:

// Indexing
// --------
// In general, a matrix is a function M: I x J -> SnzValRing.  It's sometimes handy to let I and J be something other than intervals of form [0, 1, ..., n].  So we need some general system of dealing with different types of indices.  We also want to accomodate situation where we have functions f: I' -> I and g: J' -> J, and someone has asked for the matrix N such that N[i',j'] = M[f(i'), g(j')].  N = M * (fxg).

// Sparse matrix data structures
// -----------------------------
// There are many different (reasonable) data structures for a sparse matrix.  The philosophy adopted here is to aske the person who writes the data structure to do the work of figuring out how best to export their data to all the other structures.

// Matrix representations
// ----------------------
// For M to represent a lienar map T, you generally need bijections
//      i: I -> (basis for the codomain of T)
//      j: J -> (basis for the domain of T)
// Thus, for based matrices, we need some way to specify the bases AND the bijections.


// The following sketch API is an attempt to deal with some of these issues.


// -----------------------------------------------------------------------------------------------
// RINGS
// -----------------------------------------------------------------------------------------------


// QUESTIONS TO ANSWER
// -------------------

// (1) In principle it might be possible to to write a function of form additive_identity(x: RingSpec) = IdentityEnum, where IdentityEnum is an enum that holds all the different identity elements.  This would allow the programmer to get access to an identity element, but I'm not sure how much work it would really save, since they could also get an identity element via an if-then or match statement.

// STRUCTS AND FUNCTIONS
// -----------------------

// Specifies the ring that the matrix entries belong to.  All the options below should implement the RingSpec trait (so that, for example, they can identify if an element is 1 or 0)
use core::ops::Range;
use std::collections::HashMap;
use std::collections::HashSet;
//use std::hash::Hash;
use crate::csm::CSM;

//use crate::draft_chx::BasisID;
struct BasisID{}


pub enum RingSpec{
    Integer,
    Modulus(usize),
    Rational,
    Float
}

pub struct RingMetadata<ElementType>{
    ringspec: RingSpec,
    identity_additive: ElementType,
    identity_multiplicative: ElementType
}


trait Ring<ElementType>{
    fn get_0( &self ) -> ElementType;
    fn get_1( &self ) -> ElementType;    
    fn invert( &self, ElementType ) -> ElementType;        
    fn add( &self, ElementType ) -> ElementType;        
    fn multiply( &self, ElementType ) -> ElementType;        
}




// A function that eats a RingSpec and returns a function that decides whether or not x is 0.  Each ring should have an associated type, and T should match the type of the `ring` parameter passed to the function.  For example, if ring = RingSpec::modulus(5), then `T` should be `usize`
fn additive_identity_fun<T>( ring: RingSpec) -> Box<dyn Fn(T) -> bool > {
    Box::new(|x| true)
}

// Determine whether x is 0 in this ring.
fn is_0<T>( ring: RingSpec, x: T ) -> bool {
    true
}

// Determine whether x is 1 in this ring.
fn is_1<T>( ring: RingSpec, x: T ) -> bool {
    true
}


// -----------------------------------------------------------------------------------------------
// STRUCTS AND ENUMS USED BY MATRIX ORACLES
// -----------------------------------------------------------------------------------------------


// Similar to an enum used in the Rust sprs crate, this enum is used to specify whether a matrix is stored in row-major or col-major format.
pub enum MajorDimension{
    Row,
    Col
}

// A type for sparse vectors.
pub struct SparseVector<Key, Val>{
    ind: Vec<Key>,
    snz: Vec<Val>
}

// Given an index function f: I -> X, this enum offers several different data formats for representing the set I (that is, the domain of the index function f).  If the user chooses a vector format, then all pairs of elements in the vector should be distinct.
enum IndexDomain<IndexType>{
    Range(Range<IndexType>),  // example: 0..4 means that the index set is [0, 1, 2, 3]
    Vec(Vec<IndexType>),      // example: [0, 1, 2, 3]
    Hashset(HashSet<IndexType>),
    Iterator(Box<dyn Iterator<Item=IndexType>>)   // a custom option
}

// This enum specifies an index function f: I -> I', where I is an index set with elements of type Key, and I' is another set, with elements of type I'.
enum IndexFunction<Key, Val> {
    Identity, // represents the identity function
    Funciton(Box<dyn Fn(Key) -> Val>), // an actual Rust function
    Hashmap(HashMap<Key, Val>), // a function recorded as a hash map
    Vec(Vec<Val>) // a function of form i |-> vec[i]
}

// Suppose we have a matrix M: IxJ -> Field.  We have a funciton f: J' -> J, and we'd like to create a sparse matrix representation of the matrix N: IxJ -> Field defined by N[i,j] = M[i, f(j)].  If the majs of our matrix are stored as vectors of tuples [(col_index, entry), ..., (col_index, entry)] then we need to create a pair (col_index', entry) for every col_index' in J' such that f(col_index') = col_index.  In essence, we need to compute the inverse image of col_index under f.
enum IndexInvImg<Key, Val> {
    Identity, // represents the identity function
    Funciton(Box<dyn Fn(Key) -> Vec<Val>>), // an actual Rust function
    Hashmap(HashMap<Key, Vec<Val>>), // a function recorded as a hash map
    //vec(Vec<Val>) // a function of form i |-> vec[i]
}


// -----------------------------------------------------------------------------------------------
// THE SPARSE VECTOR ORACLE (SVOracle) TRAIT
// -----------------------------------------------------------------------------------------------


// NOTES
// -----
// I wrote the following enum in the past, but I now think that a sparse vector trait is unnecessary (in particular, I the SmOracle trait provides users with everything they need).
// enum SmOracleTypes{
//     CSR: SmOracle_CSR,
//     lazyflag: SmOracle_lazyflag // ripser style matrix
// }


// -----------------------------------------------------------------------------------------------
// THE SPARSE MATRIX ORACLE (SmOracle) TRAIT
// -----------------------------------------------------------------------------------------------


// IDEA
// ----
//
// ExHACT needs sparse matrices to represent linear maps.  There are many classical data structures for a sparse matrix, and the best choice usually depends on context. Moreover, we expect that ExHACT users will regularly develop THEIR OWN data structures.  To promote interoperability, we're wrapping these in a trait that outputs data in each of four different data formats.

// KEY ASSERTIONS
// --------------
//
// WE WANT TO
// - GIVE DEVELOPERS MAXIMUM FREEDOM TO ADD FUNCTIONALITY TO THEIR ORACLES
// - MAKE AS FEW ASSUMPTIONS ABOUT THESE ORACLES AS POSSIBLE
// CONSEQUENTLY, WE WILL NOT REQUIRE THE USER TO PROVIDE MINOR SLICES (eg, A COLUMN OF A CSR MATRIX, OR A ROW OF A CSC MATRIX), BUT WE **WILL** PROVIDE THE COMMANDS THAT THEY SHOULD USE IF THE WANT TO.

// QUESTIONS + ANSWERS (SEE ALSO THE NOTES ON MATRIX OPERATIONS, BELOW)
// -------------------

//  - QUESTION: SHOULD ORACLES EVER RETURN REFERENCES?
//  - ANSWER: WE CAN ASK TO GET SLICES / REFERENCES ONLY FROM CSR STRUCTS; THIS SHOULD BE A METHOD ON THE CSR STRUCT, NOT A GENERAL TRAIT

//  - QUESTION: HOW TO DEFINE THE SET OF LEGAL ROW AND COLUMN INDICES
//  - ANSWER:
//      - WITH A WRAPPER OBJECT, NAMELY, A IndexedSmo

//  - QUESTION: HOW TO HANDLE THE TRANSPOSE OPERATION
//  - ANSWER: TRANSPOSE SHOULD BE PART OF THE SmOracle TRAIT (BECAUSE DIFFERENT DATA STRUCTURES MIGHT IMPLEMENT THE TRANSPOSE OPERATION DIFFERENTLY); HOWEVER
//      - WE SHOULDN'T REQUIRE THE TRANSPOSE OPERATION ON A STRUCT OF TYPE X RETURN ANOTHER STRUCT OF TYPE X (SINCE THERE MIGHT NOT BE A NATURAL WAY TO ENCODE THE TRANSPOSED MATRIX WITH THE SAME DATA STRUCTURE THAT X USES), AND
//      - THERE SHOULD BE A DEFAULT IMPLEMENTATION WHICH EXPORTS TO A CSR-TYPE MATRIX ORACLE FROM A IndexedSmo (WE CAN'T NECESSARILY TRANSPOSE A SMORACLE BECUASE, FOR EXAMPLE, IT MIGHT HAVE INFINITELY MANY ROWS)

//  - QUESTION: HOW TO HANDLE SCALING
//  - DISCUSSION: GREG ISN'T SURE WHETHER IT'S A GOOD IDEA TO ALLOW SOMEONE TO SCALE THE MATRIX ORACLE ITSELF; IT MIGHT MAKE MORE SENSE TO JUST SCALE THINGS VIA THE `PORT` OBJECTS.  PERHAPS WE COULD WAIT TO IMPLEMENT A SCALING METHOD UNTIL LATER, IF IT TURNS OUT WE REALLY NEED ONE.

//  - QUESTION: HOW TO ORGANIZE BINARY OPERATIONS (PRODUCT, SUM)
//      - QUESION: SHOULD WE INCLUDE SAFETY CHECKS TO MAKE SURE THAT THE INDICES FOR THE ROWS + COLUMNS OF THE TWO FACTORS MATCH ONE ANOTHER?
//      - QUESTION: SHOULD WE OFFER LAZY VERSIONS?

// - QUESTION: SHOULD WE ALLOW SCALAR MATRICES (IE SCALAR MULTIPLES OF IDENTITY MATRIX)?
// - DISCUSSION:
//      - THESE CAN BE SUPER HELPFUL, BUT THEY DO ADD SOME COMPLEXITY WHEN IMPLIMENTING MATRIX OPERATIONS LIKE PRODUCT
//      - PERHAPS WE COULD SKIRT THIS ISSUE VIA IndexedSmoS?

// - QUESTION: SHOULD WE HAVE A STRUCT TO REPRESENT AN IDENTITY MATRIX AND/OR SCALAR MATRIX?
// - DISCUSSION: THESE CAN BE QUITE HANDY, AND WOULD BE FAIRLY EASY TO IMPLEMENT.  IF WE DO THIS, THEN ONE OF THE MAIN TRICKS WILL BE FIGURING OUT HOW TO PREVENT THE SYSTEM FROM MULTIPLYING THINGS BY 1 UNNECESARILY.


// ORACLE TYPES
// ------------

// This enum lists all the structs that implement the SmOracle trait.
// The lazy_product struct doesn't contain much data, just a sequence of references to matrix SmOracles [M0, ..., Mn].  Suppose these are row-major oracles, and the user wants row k.  Then the struct will get row k from SmOracle M0, then multiply this row by M2, ..., Mn, and return the resulting vector to the user.
// The lazy_sum struct behaves similarly.
enum OracleKind</*MajKey,*/ MinKey, SnzVal>{
    Csm(CSM<MinKey, SnzVal>), // classic csr format
    //lazy_flag(Oracle_LazyCliqe<MajKey, MinKey, SnzVal>), // ripser style matrix oracle
    //lazy_product(Oracle_Product<MajKey, MinKey, SnzVal>), // this implements lazy evaluation of product
    //lazy_sum(Oracle_Sum<MajKey, MinKey, SnzVal>) // this implements lazy evaluation of sum
}

// ORACLE TRAIT
// ------------

// REMAINING QUESTIONS
// - EVERY FUNCTION (EXCEPT FOR HASH MAPS AND FUNCTIONS) RETURNS A COLLECTION OF VALUES IN A CERTAIN ORDER.  SHOULD WE REQUIRE THE ORDER TO BE THE SAME IN ALL CASES?
//   - EXAMPLE: SHOULD WE REQUIRE THAT maj_minind(k) = maj_sprsvec(k).ind?
pub trait SmOracle< MajKey: PartialEq, MinKey: PartialEq, SnzVal > {

    // Tells you what the coefficient ring is; this can't be inferreed from the SnzVal type, since the SnzVal type can't hold the modulus of a finite field (because there are infinitely many moduli, and we'll have to compile different versions of every function for each element type)
    fn ring( &self ) -> RingSpec;

    // Tells you whether the underlying struct is row-major or col-major
    fn maj_dim( &self ) -> MajorDimension;

    // GET THE VECTORS OF SNZ VALUES AND CORRESPONDING MINOR INDICES
    fn maj_indmin( &self, majkey: MajKey ) -> Option<Vec<MinKey>>;
    fn maj_snzval( &self, majkey: MajKey ) -> Option<Vec<SnzVal>>;

    // DEVELOPERS MAY RETURN NONE HERE
    fn min_indmaj( &self, minkey: MinKey ) -> Option<Vec<MajKey>>;
    fn min_snzval( &self, minkey: MinKey ) -> Option<Vec<SnzVal>>;

    // DEFAULT IMPLEMENTATION: return (self.maj_minind(..), self.maj_snzval(..))
    fn maj_sprsvec( &self, majkey: MajKey ) -> Option<SparseVector<MinKey, SnzVal>>;

    // DEFAULT IMPLEMENTATION: return (self.maj_minind(..), self.maj_snzval(..))
    // UNLESS maj_snzval OR min_snzval RETURN None; IN THAT CASE RETURN None
    fn min_sprsvec( &self, minkey: MinKey ) -> Option<SparseVector<MajKey, SnzVal>>;

    // DEFAULT IMPLEMENTATION: (SEE https://stackoverflow.com/questions/29669287/how-can-i-zip-more-than-two-iterators FOR DETAIL)
    //      RETURN izip!(self.maj_minind(majkey).itr(), self.maj_minind(majkey).itr())
    fn maj_itr( &self, majkey: MajKey ) -> Box<dyn Iterator<Item=(MinKey, SnzVal)> + '_>;

    // DEFAULT IMPLEMENTATION: (SEE https://stackoverflow.com/questions/29669287/how-can-i-zip-more-than-two-iterators FOR DETAIL)
    //      RETURN izip!(self.min_majind(minkey).itr(), self.min_majind(minkey).itr())
    //      UNLESS EITHER OF THE CALLED FUNCTIONS RETURNS None; IN THAT CASE RETURN None
    fn min_itr( &self, minkey: MinKey ) -> Box<dyn Iterator<Item=(MajKey, SnzVal)> + '_>;

    // DEFAULT IMPLEMENTATION: PUSH VALUES FROM maj_itr/min_itr INTO A HASHMAP
    fn maj_hash( &self, majkey: MajKey ) -> HashMap<MinKey, SnzVal>;
    fn min_hash( &self, minkey: MinKey ) -> HashMap<MajKey, SnzVal>;

    // FORMAT: f(minumn_index) = matrix_coefficient
    // DEFAULT: WRAP A HASHMAP IN A FUNCTION
    // For documentation on using functions as outputs, see https://doc.rust-lang.org/book/ch19-05-advanced-functions-and-closures.html
    fn maj_fn( &self, majkey: MajKey ) -> Box<dyn Fn(MinKey) -> Option<SnzVal> + '_>;
    fn min_fn( &self, majkey: MinKey ) -> Box<dyn Fn(MajKey) -> Option<SnzVal> + '_>;

    // DEFAULT IMPLEMENTATION: SEARCH FOR THE CORRECT MINOR INDEX BY RUNNING OVER AN ITERABLE PRODUCED BY self.maj_itr(majkey)
    fn entry( &self, majkey: MajKey, minkey: MinKey ) -> Option<SnzVal> {
        for (key, val) in self.maj_itr(majkey) {
            if key == minkey { return Some(val); }
        }
        return None;
    }

    // This method has a default implementation.
    fn maj_length( &self, majkey: MajKey ) -> usize {
        self.maj_itr( majkey ).count()
    }

    // This method has a default implementation. (Should return None if self.min_itr returns None)
    fn min_length( &self, minkey: MinKey ) -> usize {
        self.min_itr( minkey ).count()
    }

    // Total # of nonzeros in the matrix (if this number can actually be calculated)
    fn countsnz( &self ) -> Option<usize>;

    fn finiteminors( &self ) -> Option<bool>; // does not mean finite total matrix

    fn finitemajors( &self ) -> Option<bool>; // does not mean finite total matrix
}


// -----------------------------------------------------------------------------------------------
// PORTALS TO SPARSE MATRIX ORACLES (IndexedSmos)
// -----------------------------------------------------------------------------------------------


M[[1,2,3], [1,3,1]]

// A struct that indexes into a sparse matrix oracle (SMO).
// This struct should implement the SmOracle trait.
pub struct IndexedSmo<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal, IndexType>{
    oracle: &'a OracleKind</*MajKeyInner,*/ MinKeyInner, SnzVal>,   // a reference to a SmOracle
    indexmaj_dom: Option<Box<IndexDomain<IndexType>>>,  // the range of legal maj indices
    indexmaj_fun: Option<IndexFunction<MajKeyOuter, MajKeyInner>>, // an object implementing the maj index function
    indexmin_invimg: Option<IndexInvImg<MinKeyOuter, MinKeyInner>>,  // we need to be able to compute the inverse image of the minumn index function.  See the comment above the IndexInvImg enum for further details.

    // As a group we decided to exclude the scale factor for now. We can look for another place to let people do scaling.
    // scalefactor: Option<SnzVal>, // records whether we should scale the matrix, and if so, by what sclalar

    // I added this attributes earlier but now I think they should not be part of this struct.  The user would have to set these values themselves, and if they got it wrong then there could be major consequences.
    //If a user really wants to record the number of majs, then they can implement the ExactIterator trait on their iterator (or perhaps we could add an attribute to IndexDomain to record the length of the iterator).
    // nummaj: Option<usize>, // cardinality of the major index set (optional)
    // nummin: Option<usize>, // cardinality of the putative minor index set (optional)
    // numsnz: Option<usize>  // number of structural nonzeros (optional)
}

impl<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal, IndexType>
    IndexedSmo<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal, IndexType>{

    // Returns an iterator that runs over all the major key values we want to index.  If we think of a IndexedSmo as a submatrix M[I,J], then `majkey_iterator` is like the set I.
    fn majkey_iterator( &self ) {} //-> Box<dyn Iterator<MajKeyOuter>>

    // This function exports the IndexedSmo to a CSR SmOracle.
    fn IntoCSM( &self ) {} /* -> CSM {

        // Calculate the number of majs
        nummaj = //

        // Initialize the majptr
        majptr = Vec::with_capacity( nummaj + 1 ) // this might be a bad way to initialize ... i don't know

        // Get an iterator that will run over all the major indices.
        majkey_iterator = self.majkey_iterator()

        // Create a hashmap from keys to integers.
        majkey_hash = HashMap<MakKey_outer, usize>

        // Use the iterator in a for-loop to calculate the majptr vector
        for (slice_number, maj_key) in enumerate(majkey_iterator){
            majkey_hash.push(maj_key, slice_number)
            numsnz_thismaj = // this might be tricky to calculate ... unless this Port uses all the columns that the oracle provides, we'll have to count just the ones that the Port asks for (possibly multiple times, depending on indexrel_min)
        }

        // Initialize the minind and snzval vectors
        numsnz = majkey_iterator.last()
        minind = Vec::with_capacity(numsnz)
        snzval = Vec::with_capacity(numsnz) // not sure how to make this initialize as a vector with entries of the proper type

        // Fill in the minind and snzval vectors; this process will probably look a lot like the algorithm for transposing a CSR matrix
        ..

        // Wrap all the calculated vectors into a new CSR struct!
    }*/
}



// Add later
/*
impl<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal, IndexType>
    SmOracle<MajKeyOuter, MinKeyOuter, SnzVal> for
    IndexedSmo<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal, IndexType>{

    // We have to remember to scale the entries by `scalefactor`.

    fn maj_indmin( &self, majkey: MajKeyOuter ) -> Option<Vec<MinKeyOuter>> {
        None
    }
    fn maj_snzval( &self, majkey: MajKeyOuter ) -> Option<Vec<SnzVal>> {
        None
    }
    fn min_indmaj( &self, majkey: MinKeyOuter ) -> Option<Vec<MajKeyOuter>> {
        None
    }
    fn min_snzval( &self, majkey: MinKeyOuter ) -> Option<Vec<SnzVal>> {
        None
    }

    // ... etc.

}
*/




// -----------------------------------------------------------------------------------------------
// SUMMARY OF MATRIX OPERATIONS
// -----------------------------------------------------------------------------------------------


// UNARY
// -----------------------------------------------------------------------------------------------


// SLICE (THAT IS, TURNING M INTO M[I,J])
// --------------------------------------
// IndexedSmos do this job.  If you want to take a slice of an oracle, just make a IndexedSmo around it.  If you want to take a slice of a IndexedSmo, then you just have to (manually) change the index functions.

// TRANSPOSE FOR SMORacles
// -----------------------
// The author of each matrix oracle is responsible for writing their own transpose function. (In many cases, they can simply export to CSR format, then use the built-in function to transpose CSR's).

// TRANSPOSE FOR IndexedSmos
// ----------------------
// I am exhausted.  I don't think it's essential to implement the transpose operation for IndexedSmos right now.  Let's leave it for later.

// This function should return an error if IndexedSmo.indexmaj_dom does not provide a finite set (if IndexedSmo.indexmaj_domain = Everything).
fn transpose<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal, IndexType>
    (S: IndexedSmo<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal, IndexType>) {}// -> IndexedSmo<MinKey, usize, MajKey, MajKey, Element, >;
    // IndexedSmo<MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal>


// SCALING FOR SMORacles
// ---------------------
// The easiest thing is to wrap your SMORacle in a IndexedSmo with an appropriate scalefactor.  If you don't want to specify an index domain / function / inverse image you don't have to ... just use the default `identity_function` values for the relevant enums.  If a user doesn't want to do that, then they can write their own function for modifying a particular implemenation of SmOracle.

// SCALING FOR IndexedSmos
// --------------------
// Just change the scale factor.


// BINARY
// -----------------------------------------------------------------------------------------------


// STRATEGY
// --------
// The general strategy for binary matrix operations (product and sum) is to implement `lazy` realizations on SMORacles.  For further details, see the description of Oracle_Product and Oracle_Sum above.
// If users want, they can wrap the resulting lazy object in a IndexedSmo and export to CSR format.
// If users need something with better performance, then they can write it themselves!

// EXAMPLES
// --------

// This function returns a lazy `product` SmOracle, but the function doesn't have much work to do.  It just places `v` into a SmOracle_Product struct.
// **All SmOracles should have the same major dimension.**
fn lazy_prod<SmOracleTypes>( sequence: Vec<Option<SmOracleTypes>> ){}// -> Oracle_Product;

// This function returns a lazy `sum` SmOracle, but the function doesn't have much work to do.  It just places `v` into a SmOracle_Sum struct.
// **All SmOracles should have the same major dimension.**
fn lazy_sum<SmOracleTypes> ( sequence: Vec<Option<SmOracleTypes>> ){}// -> Oracle_Sum; //  returns a lazy 'sum' SmOracle

// NOTES
// -----
// WE SHOULD PROBABLY HAVE SOME WAY TO CHECK THAT INDICES ALIGN .. AT LEAST, THE USER SHOULD HAVE THE OPTION TO DO THIS IN CERTAIN SPECIAL CASES, LIKE THE CASE WHERE EVERY IndexedSmo HAS AN INDEX DOMAIN OF FORM [0, .., N]


// -----------------------------------------------------------------------------------------------
// MATRIX VECTOR OPERATIONS
// -----------------------------------------------------------------------------------------------


// There are essentially three different sparse vector formats compatible with sparse matrix / vector multiplication, and all three can be converted to iterators at very low cost.  Therefore I'm proposing that we format the `vector` argument as an iterator.

fn prod_major<T, MajKey, MinKey, SnzVal>(matrix: Box<dyn SmOracle<MajKey, MinKey, SnzVal>>, vector: T, majdim: MajorDimension)-> SparseVector<MinKey, SnzVal>
where
    T: Iterator<Item=(MajKey, SnzVal)>,
{
    SparseVector {
        ind: Vec::new(),
        snz: Vec::new()
    }
}// -> SparseVector // the major dimension decides whether we multiply on the left or the right.

fn prod_minor<T, MajKey, MinKey, SnzVal>(matrix: Box<dyn SmOracle<MajKey, MinKey, SnzVal>>, vector: Vec<T>, majdim: MajorDimension) -> SparseVector<MajKey, SnzVal>
where
    T: Iterator<Item=(MajKey, SnzVal)>,
{
    SparseVector {
        ind: Vec::new(),
        snz: Vec::new()
    }
}// -> SparseVector // the major dimension decides whether we multiply on the left or the right.  Note that we have to use a IndexedSmo in this conext, since we need to know the IndexDomain of the major dimension.


// -----------------------------------------------------------------------------------------------
// VECTOR VECTOR OPERATIONS
// -----------------------------------------------------------------------------------------------

fn sum_vector<T, MajKey, SnzVal>( vector1: T, vector2: T ) -> SparseVector<MajKey, SnzVal>
where
    T: Iterator<Item=(MajKey, SnzVal)>,
{
    SparseVector {
        ind: Vec::new(),
        snz: Vec::new()
    }
}

// -----------------------------------------------------------------------------------------------
// ASSOCIATING ORACLES WITH BASES (STILL UNDER CONSTRUCTION)
// -----------------------------------------------------------------------------------------------


// GREG PROPOSES TABLING THE FOLLOWING IDEAS UNTIL LATER; WE MAY NOT NEED THEM
// ---------------------------------------------------------------------------
//
// REASONS FOR ABANDONING: IT SEEMS THAT KEEPING TRACK OF BASES AND THEIR INDEXING SYSTEMS IS HARD IN GENERAL.  THE SPECIFIC PROBLEM THAT MOTIVATED MUCH OF OUR CAN BE SOLVED MORE EFFECTILY WITH CUSTOMIZED SOLUTION (SPECIFIC TO CHAIN COMPLEXES) THAN WITH A GENERAL ONE.  WE CAN FOCUS ON THE CASE OF CHAIN COMPLEXES FOR NOW, AND IF MORE GENERAL PROBLEMS ARISE IN THE FUTURE, WE CAN REVISIT THE PROBLEM OF TRACKING BASES AT THAT TIME.
//
// struct BasisBijectionSpec{
//     subquo_spec : SubquoSpec,
//
// }
//
// enum BijectiveFunction{
//     Identity,
//     Funciton(Box<dyn Fn()->()>)
// }
//
// enum BijectionDomain<T>{
//     Range,
//     Vec(),
// }
//
// /// A unique id for a specific bijection between a basis and a set of index keys.  In general, the bijection itself can't be derived from the data in this object; it's just a 'name.'
// struct BasisBijectionSpec{
//     basisID: BasisID,
//     bijectiondomain: BijectionDomain,
//     bijectivefunction: BijectiveFunction,
// }
//
// /// This struct gives complete data about a matrix representation of linear map stored on computer.  In particular, suppose
// ///     T: a linear map
// ///     B: a basis for the domain
// ///     C: a basis for the codomain
// ///     X: a bijection I -> B
// ///     Y: a biJection J -> C
// ///     M: an I x J matrix such that M[i,j] is the proper coefficient of the matrix rep of T with respect to B and C
// /// This struct provides 'names' for X and Y (in the forrm of `BasisIndexSpec`s) and an oracle to evaluate the majs of M.
// struct BasedCsm<MinKey, SnzVal>{
//     csm: CSM<MinKey, SnzVal>,
//     ind_of_majs: BasisBijectionSpec,
//     ind_of_mins: BasisBijectionSpec
// }
