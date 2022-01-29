// ================================================================================================
// FILE DESCRIPTION
// This file contains some outlines matrix (oracle) objects and traits.
//
// FILE HISTORY
// 20201120: created by Greg Henselman-Petrusek
// 20201120-present: edited by Greg Henselman-Petrusek and Haibin Hang
// ================================================================================================


/*!
 Sparse matrix oracles

 # Sparse matrix oracles

A *sparse matrix oracle* is an object that allows users to access rows, columns, and other information about a sparse matrix.  This library provides a standardized set of commands that you can use to access information from any sparse matrix oracle (you don't need different commands for different oracles).  This standard set of commands is formalized as a [Rust trait](https://doc.rust-lang.org/book/ch10-02-traits.html).

* a [Rust trait](https://doc.rust-lang.org/book/ch10-02-traits.html) is a set of commands that can be applied in the same way to several different objects
* the [`SmOracle`](crate::matrix::SmOracle) trait is set of commands for accessing information about a sparse matrix
* a sparse matrix oracle is defined as any [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that implements the commands in the [`SmOracle`](crate::matrix::SmOracle) trait

# Quick facts

* List of oracles
     * A list of sparse matrix oracles can be found at the bottom of the [`SmOracle`](crate::matrix::SmOracle) documentation page.
     * **You can  build your own oracle**!  In concrete terms, this means defining a [structure](https://doc.rust-lang.org/book/ch05-01-defining-structs.html) that [implements](https://doc.rust-lang.org/book/ch10-02-traits.html) the [`SmOracle`](crate::matrix::SmOracle) trait.
 * Matrix operations
     * There are serval matrix operations you can perform with matrix oracles; for example, you can multiply a matrix with a sparse vector, obtain an *R = DV* decomposition, and solve *Ax = b* for *x*.  See the [`solver`](crate::solver) module for more details.

# Examples


 ```
 Example 1: (to fill in later)
     - create a row-major csm with usize snz vals and columns indexed by floats
     - access a row as an iterator
     - access a column as a sparse vector

 Example 2: (to fill in later)
     - create the same csm, but with columns indexed by usize
     - access the same row as in Example 1

 Example 3: (to fill in later)
     - create a csm with snz values in Z/3Z
     - multiply the matrix with a sparse vector

 ```

*/





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

use core::ops::Range;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;
use std::fmt::Debug;
use crate::csm::CSM;
use num::rational::Ratio;


struct BasisID{}


// -----------------------------------------------------------------------------------------------
// RINGS
// -----------------------------------------------------------------------------------------------


// QUESTIONS TO ANSWER
// -------------------

// (1) In principle it might be possible to to write a function of form additive_identity(x: RingSpec) = IdentityEnum, where IdentityEnum is an enum that holds all the different identity elements.  This would allow the programmer to get access to an identity element, but I'm not sure how much work it would really save, since they could also get an identity element via an if-then or match statement.

// STRUCTS AND FUNCTIONS
// -----------------------

/// Represents the "name" of a given rings.  Unlike a [`RingMetadata`](RingMetadata) object, it does not contain sufficient information for Rust to infer how to execute addition, subtraction, multiplication, and (sometimes) division.
#[derive(Clone)]
pub enum RingSpec{
    Integer,
    Modulus(usize),
    Rational,
    Float
}

/// Stores the data needed for rust to perform basic ring operations:
/// * addition
/// * subtraction
/// * multiplication
/// * (sometimes) division
#[derive(Clone)]
pub struct RingMetadata<ElementType: Clone>{
    pub ringspec: RingSpec,                     // the type of ring
    pub identity_additive: ElementType,         // additive identity
    pub identity_multiplicative: ElementType    // multiplicative identiry
}


type IntegerType = i16;
/// The Euclidean algorithm
///
/// # Parameters
/// - `prime`: an integer, most of the case is a prime number
/// - `integer`: an interger
///
/// # Returns
/// - an interger named 'inverse' such that '(inverse*integer)%prime = gcd(prime,integer)'
/// the greatest common divisor gcd(prime,integer)
pub fn euclidean(prime: IntegerType, integer:IntegerType) -> (IntegerType,IntegerType) {
    let mut integer = integer%prime;
    if integer < 0 { integer += prime; }
    let mut x;
    let mut y0 = prime;
    let mut y1 = integer;
    let mut a0 = 0;
    let mut a1 = 1;
    while y1>0 && y0%y1 != 0 {
        x = y0/y1;
        let temp = y0;
        y0 = y1;
        y1 = temp%y0;
        let temp = a0;
        a0 = a1;
        a1 = (temp - x*a0)%prime;
    }
    if a1 < 0 {
        a1 += prime;
    }
    return (a1,y1);
}

/// A trait which makes sure we can inverse coefficients of different types. Even though each type has their own inverse method, to work on them uniformly, we need to define this trait to unify their inverse operations
pub trait InvMod {
    type Output; // the output type

    /// Find the inverse module given (prime) interger
    ///
    /// # Parameters
    /// -`module`: given (prime) interger
    ///
    /// # Returns
    /// If inverse exists return Some(inverse); otherwise return None
    fn inverse_mod(&self, module: usize) -> Option<Self::Output>;

    /// module by given interger
    ///
    /// # Parameters
    /// -`module`: given (prime) interger
    fn module_by(&self, module: usize) -> Self::Output;
}

/// Implimentation of trait InvMod for rational type
impl InvMod for Ratio<IntegerType>{
    type Output = Ratio<IntegerType>;
    fn inverse_mod(&self, _prime: usize) -> Option<Self::Output>{
        if *self.numer() == 0 || *self.denom() == 0 { return None; }
        else {
            return Some( Ratio::new(*self.denom(), *self.numer()) );
        }
    }

    fn module_by(&self, prime: usize) -> Self::Output {
        return self.clone();
    }
}

/// Implimentation of trait InvMod for interger type
impl InvMod for IntegerType{
    type Output = IntegerType;

    fn inverse_mod(&self, prime: usize) -> Option<Self::Output>{
        if euclidean(prime as IntegerType, *self).1 != 1 { return None; }
        else { return Some(euclidean(prime as IntegerType, *self).0); }
    }

    fn module_by(&self, prime: usize) -> Self::Output {
        let mut simple_form = *self%(prime as IntegerType);
        if simple_form < 0 { simple_form = simple_form + prime as IntegerType; }
        return simple_form;
    }
}

/// Methods of RingMetadata struct
impl<ElementType: Clone + PartialEq + InvMod<Output=ElementType>> RingMetadata<ElementType> {

    /// Find the inverse of x and returns None if does not exist
    pub fn inverse( &self, x: &ElementType ) -> Option<ElementType> {
        match self.ringspec {
            RingSpec::Integer => {
                return Some(x.clone());
            }
            RingSpec::Modulus(prime) => {
                x.inverse_mod(prime)
            }
            RingSpec::Rational => {
                x.inverse_mod(0)
            }
            RingSpec::Float => { return Some(x.clone()); }
        }
    }

    /// Simplify the representation of x. For example, 6 in Modulus(5) is simplified to 1; 6/3 in rational is simplified to 2/1.
    pub fn simplify(&self, x: &ElementType) -> ElementType {
        match self.ringspec {
            RingSpec::Integer => {
                return x.clone();
            }
            RingSpec::Modulus(prime) => {
                x.module_by(prime)
            }
            RingSpec::Rational => { return x.clone(); }
            RingSpec::Float => { return x.clone(); }
        }
    }

    /// Determine whether x is 0 in this ring.
    pub fn is_0( &self, x: &ElementType ) -> bool {
        match self.inverse(x) {
            Some(_) => { return false; }
            None => { return true; }
        }
    }

}




// -----------------------------------------------------------------------------------------------
// STRUCTS AND ENUMS USED BY MATRIX ORACLES
// -----------------------------------------------------------------------------------------------


/// Similar to an enum used in the Rust sprs crate, this enum is used to specify whether a matrix is stored in row-major or col-major format.
#[derive(Clone, PartialEq)]
pub enum MajorDimension{
    Row,
    Col
}

/// A struct for sparse vectors.
pub struct SparseVector<Key, Val>{
    pub ind: Vec<Key>, // vector of indices
    pub snz: Vec<Val>  // vector of structural non zero entries
}

/// Methods of SparseVector
impl<Key: Debug, Val: Debug> SparseVector<Key, Val>{

    /// Construct a trivial SparseVector with no entries
    pub fn new() -> SparseVector<Key, Val> {
        SparseVector{
            ind: Vec::new(),
            snz: Vec::new()
        }
    }

    /// Print out the SparseVector
    pub fn print_sparse_vector(&self) {
        println!("  ind: {:?}", &self.ind);
        println!("  snz: {:?}", &self.snz);
    }
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
    Function(Box<dyn Fn(Key) -> Val>), // an actual Rust function
    Hashmap(HashMap<Key, Val>), // a function recorded as a hash map
    Vec(Vec<Val>) // a function of form i |-> vec[i]
}
// Original : M(I,J), f: I'-> I, g: J'-> J
//
// Suppose we have a matrix M: IxJ -> Field.  We have a function f: J' -> J, and we'd like to create a sparse matrix representation of the matrix N: IxJ -> Field defined by N[i,j] = M[i, f(j)].  If the majs of our matrix are stored as vectors of tuples [(col_index, entry), ..., (col_index, entry)] then we need to create a pair (col_index', entry) for every col_index' in J' such that f(col_index') = col_index.  In essence, we need to compute the inverse image of col_index under f.
enum IndexInvImg<Key, Val> {
    Identity, // represents the identity function
    Function(Box<dyn Fn(Key) -> Vec<Val>>), // an actual Rust function
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


//  - QUESTION: SHOULD WE HAVE SO MANY OUTPUT FORMATS FOR MAJOR/MINOR FIELDS?
//  - ANSWER: KEEPING THE REQUIRED IMPLEMENTATIONS TO A SMALL NUMBER WHILE HAVING LOTS OF PROVIDED IMPLEMENTATIONS SEEMS LIKE A REASONABLE WAY TO LIFT BURDEN FROM DEVELOPERS WHILE PROVIDING A CONSISTENT / MODULAR FRAMEWORK FOR DIFFERENT DATA TYPES

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
//      - PERHAPS WE COULD SKIRT THIS ISSUE VIA SMOPORTS?

// - QUESTION: SHOULD WE HAVE A STRUCT TO REPRESENT AN IDENTITY MATRIX AND/OR SCALAR MATRIX?
// - DISCUSSION: THESE CAN BE QUITE HANDY, AND WOULD BE FAIRLY EASY TO IMPLEMENT.  IF WE DO THIS, THEN ONE OF THE MAIN TRICKS WILL BE FIGURING OUT HOW TO PREVENT THE SYSTEM FROM MULTIPLYING THINGS BY 1 UNNECESARILY.


// ORACLE TYPES
// ------------

// This enum lists all the structs that implement the SmOracle trait.
// The lazy_product struct doesn't contain much data, just a sequence of references to matrix SmOracles [M0, ..., Mn].  Suppose these are row-major oracles, and the user wants row k.  Then the struct will get row k from SmOracle M0, then multiply this row by M2, ..., Mn, and return the resulting vector to the user.
// The lazy_sum struct behaves similarly.
enum OracleKind</*MajKey,*/ MinKey, SnzVal: Clone>{
    Csm(CSM<MinKey, SnzVal>), // classic csr format
    //Clique(CliqueBoundaryMatrix<Simplex, Simplex, i32>), // ripser style matrix oracle
    //lazy_product(Oracle_Product<MajKey, MinKey, SnzVal>), // this implements lazy evaluation of product
    //lazy_sum(Oracle_Sum<MajKey, MinKey, SnzVal>) // this implements lazy evaluation of sum
}

// ORACLE TRAIT
// ------------



// REMAINING QUESTIONS
// - EVERY FUNCTION (EXCEPT FOR HASH MAPS AND FUNCTIONS) RETURNS A COLLECTION OF VALUES IN A CERTAIN ORDER.  SHOULD WE REQUIRE THE ORDER TO BE THE SAME IN ALL CASES?
//   - EXAMPLE: SHOULD WE REQUIRE THAT maj_minind(k) = maj_sprsvec(k).ind?
/**
Methods for sparse matrix oracles

# Overview

There are several types of information you might want from a sparse matrix:

* **Index sets** A *matrix* is a rectangle of coefficients *M[i,j]*, where *i* runs over some index set *I* and *j* runs over some index set *J*.  We don't assume that *I* or *J* are sets of integers (they could be sets of formal symbols, for example).

* **Coefficient ring** The coefficients *M[i,j]* lie in some ring, *R*.  You need to know the ring in order to add, subtract, multiply, and (sometimes) divide.  You also need to know the additive and multiplicative identities.

* **Major dimension** A matrix is **column-major** if it is easier to access columns than rows, and **row-major** if the converse holds.  (Traditional computer science literature uses this term a little differently; there *row-major* and *column-major*  refer to the way that entries are stored in memory.  Here we don't assume anything about the way information is stored, so major dimension is technically just a "label" that someone has assigned to the matrix.)

* **Rows and columns = major and minor fields** We call the rows and columns of a row-major matrix its "major fields" and "minor fields," respectively.  The reverse holds for column-major matrices.

The `SmOracle` trait gives a general set of conventions for extracting all this and more.

# Examples

A list of sparse matrix oracles (that is, implementors of the `SmOracle` trait) can be found at the bottom of this page.

 ```
 Example 1: (to fill in later)
     - create a row-major csm with columns indexed by floats
     - access a row of the csm
     - access a column of the csm

 Example 2: (to fill in later)
     - create the same row-major csm, but with columns indexed by usize
     - access a row of the csm
     - access a column of the csm

Example 3: (to fill in later) -- maybe this is too many examples?
     - create the same row-major csm, but with columns indexed by simplices
     - access a row of the csm
     - access a column of the csm

Example 4: (to fill in later) -- maybe this is too many examples?
     - create the same row-major csm, but with columns indexed by cubes
     - access a row of the csm
     - access a column of the csm
 ```

# Common abbreviations

The type parameters for this trait can be interpreted as follows:

* `MajKey` the type of the elements that index the major dimension
* `MinKey` the type of the elements that index the minor dimension
* `SnzVal` the type of the the structural nonzero entries (this can be the same for several different rings)

*/
pub trait SmOracle<MajKey, MinKey, SnzVal> where
MajKey: PartialEq + Eq + Hash + Clone,
MinKey: PartialEq + Eq + Hash + Clone,
SnzVal: Clone
{

    /// The coefficient ring of the matrix, represented by a [`RingMetadata`](RingMetadata) object.
    fn ring(&self) -> &RingMetadata<SnzVal>;

    /// The major dimension of the oracle (row or column)
    fn maj_dim(&self) -> MajorDimension { MajorDimension::Row }

    // DEFAULT IMPLEMENTATION: (SEE https://stackoverflow.com/questions/29669287/how-can-i-zip-more-than-two-iterators FOR DETAIL)
    //      RETURN izip!(self.maj_minind(majkey).itr(), self.maj_minind(majkey).itr())
    /// The major field indexed by `majkey`, represented by an iterator of tuples of form `(minor_index, coefficient)`.
    fn maj_itr(&self, majkey: &MajKey) -> Box<dyn Iterator<Item=(MinKey, SnzVal)> + '_>;

    // DEFAULT IMPLEMENTATION: (SEE https://stackoverflow.com/questions/29669287/how-can-i-zip-more-than-two-iterators FOR DETAIL)
    //      RETURN izip!(self.min_majind(minkey).itr(), self.min_majind(minkey).itr())
    //      UNLESS EITHER OF THE CALLED FUNCTIONS RETURNS None; IN THAT CASE RETURN None
    /// The minor field indexed by `minkey`, represented by an iterator of tuples of form `(major_index, coefficient)`.
    fn min_itr(&self, minkey: &MinKey) -> Box<dyn Iterator<Item=(MajKey, SnzVal)> + '_>;

    // GET THE VECTORS OF SNZ VALUES AND CORRESPONDING MINOR INDICES
    /// The minor indices of the structural nonzero entries in the major field indexed by `majkey`, represented as a vector.
    fn maj_indmin(&self, majkey: &MajKey) -> Option<Vec<MinKey>> {
        let mut vec = Vec::new();
        for (key, _) in self.maj_itr(majkey) {
            vec.push(key.clone());
        }
        if vec.len() == 0 { return None; }
        else { return Some(vec); }
    }

    /// The structural nonzero entries in the major field indexed by `majkey`, represented as a vector.
    fn maj_snzval(&self, majkey: &MajKey) -> Option<Vec<SnzVal>> {
        let mut vec = Vec::new();
        for (_, val) in self.maj_itr(majkey) {
            vec.push(val);
        }
        if vec.len() == 0 { return None; }
        else { return Some(vec); }
    }

    // DEVELOPERS MAY RETURN NONE HERE
    /// The major indices of the structural nonzeros in field `minkey`, represented by a vector.
    fn min_indmaj(&self, minkey: &MinKey) -> Option<Vec<MajKey>> {
        let mut vec = Vec::new();
        for (key, _) in self.min_itr(minkey) {
            vec.push(key.clone());
        }
        if vec.len() == 0 { return None; }
        else { return Some(vec); }
    }

    /// The structural nonzeros in field `minkey`, represented by a vector.
    fn min_snzval(&self, minkey: &MinKey) -> Option<Vec<SnzVal>> {
        let mut vec = Vec::new();
        for (_, val) in self.min_itr(minkey) {
            vec.push(val.clone());
        }
        if vec.len() == 0 { return None; }
        else { return Some(vec); }
    }

    /// The major field indexed by `majkey`, represented by a sparse vector.
    fn maj_sprsvec(&self, majkey: &MajKey) -> Option<SparseVector<MinKey, SnzVal>> {
        let mut ind = Vec::new();
        let mut snz = Vec::new();
        for (key, val) in self.maj_itr(majkey) {
            ind.push(key);
            snz.push(val);
        }
        if ind.len() == 0 { return None; }
        return Some(SparseVector{
            ind: ind,
            snz: snz
        });
    }

    // DEFAULT IMPLEMENTATION: return (self.maj_minind(..), self.maj_snzval(..))
    // UNLESS maj_snzval OR min_snzval RETURN None; IN THAT CASE RETURN None
    /// The minor field indexed by `minkey`, represented by a sparse vector.
    fn min_sprsvec(&self, minkey: &MinKey) -> Option<SparseVector<MajKey, SnzVal>> {
        let mut ind = Vec::new();
        let mut snz = Vec::new();
        for (key, val) in self.min_itr(minkey) {
            ind.push(key);
            snz.push(val);
        }
        if ind.len() == 0 { return None; }
        return Some(SparseVector{
            ind: ind,
            snz: snz
        });
    }

    // DEFAULT IMPLEMENTATION: PUSH VALUES FROM maj_itr/min_itr INTO A HASHMAP
    /// The major field indexed by `majkey`, represented by a hashmap.
    fn maj_hash(&self, majkey: &MajKey) -> HashMap<MinKey, SnzVal> {
        let mut res = HashMap::new();
        for (key, val) in self.maj_itr(majkey) {
            res.insert(key,val);
        }
        return res;
    }

    /// The minor field indexed by `minkey`, represented by a hash map
    fn min_hash(&self, minkey: &MinKey) -> HashMap<MajKey, SnzVal> {
        let mut res = HashMap::new();
        for (key, val) in self.min_itr(minkey) {
            res.insert(key,val);
        }
        return res;
    }

    // FORMAT: f(minumn_index) = matrix_coefficient
    // DEFAULT: WRAP A HASHMAP IN A FUNCTION
    // For documentation on using functions as outputs, see https://doc.rust-lang.org/book/ch19-05-advanced-functions-and-closures.html
    /// The major field indexed by `majkey`, represented by a function.
    fn maj_fn<'a>(&'a self, majkey: &'a MajKey) -> Box<dyn Fn(MinKey) -> Option<SnzVal> + 'a> {
        Box::new( move |x| {
            for (key, val) in self.maj_itr(majkey) {
                if key == x { return Some(val); }
            }
            return None;
        })
	}

    /// The minor field indexed by `minkey`, represented by a function.
	fn min_fn<'a>(&'a self, minkey: &'a MinKey) -> Box<dyn Fn(MajKey) -> Option<SnzVal> + 'a> {
        Box::new( move |x| {
            for (key, val) in self.min_itr(minkey) {
                if key == x { return Some(val); }
            }
            return None;
        })
	}

    // DEFAULT IMPLEMENTATION: SEARCH FOR THE CORRECT MINOR INDEX BY RUNNING OVER AN ITERABLE PRODUCED BY self.maj_itr(majkey)
    /// Coefficient in minor field `majkey`, minor field `minkey`.
    fn entry(&self, majkey: &MajKey, minkey: &MinKey) -> Option<SnzVal> {
        for (key, val) in self.maj_itr(majkey) {
            if key == *minkey { return Some(val); }
        }
        return None;
    }

    /// The number of structural nonzero entries in major field `majkey`.
    fn maj_length(&self, majkey: &MajKey) -> usize {
        self.maj_itr(majkey).count()
    }

    // This method has a default implementation. (Should return None if self.min_itr returns None)
    /// The number of structural nonzero entries in major field `majkey`.
    fn min_length(&self, minkey: &MinKey) -> usize {
        self.min_itr(minkey).count()
    }

    /// Total number of nonzeros in the matrix (if this number can actually be calculated)
	fn countsnz(&self) -> Option<usize> {
		None
	}

    /// True if every minor field has finitely many structural nonzero entries
	fn finiteminors(&self) -> Option<bool> {
		Some(true)
	}

    /// True if every major field has finitely many structural nonzero entries
	fn finitemajors(&self) -> Option<bool> {
		Some(true)
	}

    /// Function that tells if there is a way to judge this location if pivot or not before decomposition. This is especially desighed for Rips complex.
    fn is_pivot(&self, majkey: &MajKey, minkey: &MinKey) -> Option<bool> {
        return None;
    }

    ///
    fn is_apparent(&self, minkey: &MinKey) -> Option<MajKey> {
        return None;
    }
}


// -----------------------------------------------------------------------------------------------
// INDEXED SPARSE MATRIX ORACLES (IndexedSmo)
// -----------------------------------------------------------------------------------------------


// A struct that indexes into a sparse matrix oracle (SMO).
// This struct should implement the SmOracle trait.
/// Sparse matrix oracle that indexes into another sparse matrix oracle
pub struct IndexedSmo<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal: Clone, IndexType>{
    oracle: &'a OracleKind</*MajKeyInner,*/ MinKeyInner, SnzVal>,   // a reference to a SmOracle
    indexmaj_dom: Option<Box<IndexDomain<IndexType>>>,  // the range of legal maj indices
    indexmaj_fun: Option<IndexFunction<MajKeyOuter, MajKeyInner>>, // an object implementing the maj index function
    indexmin_invimg: Option<IndexInvImg<MinKeyOuter, MinKeyInner>>,  // we need to be able to compute the inverse image of the minumn index function.  See the comment above the IndexInvImg enum for further details.
    scalefactor: Option<SnzVal>, // records whether we should scale the matrix, and if so, by what sclalar

    // I added this attributes earlier but now I think they should not be part of this struct.  The user would have to set these values themselves, and if they got it wrong then there could be major consequences.
    //If a user really wants to record the number of majs, then they can implement the ExactIterator trait on their iterator (or perhaps we could add an attribute to IndexDomain to record the length of the iterator).
    // nummaj: Option<usize>, // cardinality of the major index set (optional)
    // nummin: Option<usize>, // cardinality of the putative minor index set (optional)
    // numsnz: Option<usize>  // number of structural nonzeros (optional)
}

impl<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal: Clone, IndexType>
    IndexedSmo<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal, IndexType>{

    /// Returns an iterator that runs over all the major key values we want to index.  If we think of a IndexedSmo as a submatrix M[I,J], then `majkey_iterator` is like the set I.
    fn majkey_iterator( &self ) {} //-> Box<dyn Iterator<MajKeyOuter>>

    /// This function exports the IndexedSmo to a CSR SmOracle.
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




// ------------------------------------------------------------------------------------------------
// SUMMARY OF MATRIX OPERATIONS
// ------------------------------------------------------------------------------------------------


// UNARY
// ------------------------------------------------------------------------------------------------


// SLICE (THAT IS, TURNING M INTO M[I,J])
// --------------------------------------
// SmoPorts do this job.  If you want to take a slice of an oracle, just make a IndexedSmo around it.  If you want to take a slice of a IndexedSmo, then you just have to (manually) change the index functions.

// TRANSPOSE FOR SMORacles
// -----------------------
// The author of each matrix oracle is responsible for writing their own transpose function. (In many cases, they can simply export to CSR format, then use the built-in function to transpose CSR's).

// TRANSPOSE FOR SmoPorts
// ----------------------
// I am exhausted.  I don't think it's essential to implement the transpose operation for SmoPorts right now.  Let's leave it for later.

/// This function should return an error if IndexedSmo.indexmaj_dom does not provide a finite set (if IndexedSmo.indexmaj_domain = Everything).
fn transpose <'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal: Clone, IndexType>
    (S: IndexedSmo<'a, MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal, IndexType>) {

    }// -> IndexedSmo<MinKey, usize, MajKey, MajKey, Element, >;
    // IndexedSmo<MajKeyOuter, MajKeyInner, MinKeyOuter, MinKeyInner, SnzVal>


// SCALING FOR SMORacles
// ---------------------
// The easiest thing is to wrap your SMORacle in a IndexedSmo with an appropriate scalefactor.  If you don't want to specify an index domain / function / inverse image you don't have to ... just use the default `identity_function` values for the relevant enums.  If a user doesn't want to do that, then they can write their own function for modifying a particular implemenation of SmOracle.

// SCALING FOR SmoPorts
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

/// This function returns a lazy `product` SmOracle, but the function doesn't have much work to do.  It just places `v` into a SmOracle_Product struct.
/// **All SmOracles should have the same major dimension.**
fn lazy_prod<SmOracleTypes>( sequence: Vec<Option<SmOracleTypes>> ){}// -> Oracle_Product;

/// This function returns a lazy `sum` SmOracle, but the function doesn't have much work to do.  It just places `v` into a SmOracle_Sum struct.
/// **All SmOracles should have the same major dimension.**
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

enum BijectiveFunction{
    Identity,
    Funciton(Box<dyn Fn()->()>)
}

enum BijectionDomain{
    Range,
    Array
}

/// A unique id for a specific bijection between a basis and a set of index keys.  In general, the bijection itself can't be derived from the data in this object; it's just a 'name.'
struct BasisBijectionSpec{
    basisID: BasisID,
    bijectiondomain: BijectionDomain,
    bijectivefunction: BijectiveFunction,
}

/// This struct gives complete data about a matrix representation of linear map stored on computer.  In particular, suppose
///     T: a linear map
///     B: a basis for the domain
///     C: a basis for the codomain
///     X: a bijection I -> B
///     Y: a biJection J -> C
///     M: an I x J matrix such that M[i,j] is the proper coefficient of the matrix rep of T with respect to B and C
/// This struct provides 'names' for X and Y (in the forrm of `BasisIndexSpec`s) and an oracle to evaluate the majs of M.
struct BasedCsm<MinKey, SnzVal: Clone>{
    csm: CSM<MinKey, SnzVal>,
    ind_of_majs: BasisBijectionSpec,
    ind_of_mins: BasisBijectionSpec
}
