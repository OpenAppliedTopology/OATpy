/*!

Algebraic operations for sparse matrix oracles

# Algegraic operations for sparse matrix oracles

```
Example:
     - perform triangular solve for x such that Ax = b

Example:
     - check the solution by confirming that Ax = b

```


*/

use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::ops::{Neg, AddAssign, Mul};
use std::cmp::Reverse;
use std::hash::Hash;
use std::fmt::Debug;
use crate::csm::CSM;
use crate::decomp_row::update_heap_hash; // Haibin: maybe we should just put update_heap_hash() in this solver.rs file
use crate::matrix::{SparseVector, InvMod, SmOracle, RingMetadata};
use crate::chx::Indexing;

/// Solve upper*xx = bb
///
/// # Parameters
/// - `uppper`: an upper triangular matrix in CSM format with diagnal entries to be 1
/// - `bb`: constant column sparse vector
/// # Returns
/// xx the unknow sparse vector
pub fn triangular_solver<SnzVal: Neg + Neg<Output=SnzVal> + AddAssign + Mul + Mul<Output=SnzVal> + Clone + InvMod + InvMod<Output = SnzVal> + PartialEq + Debug>(
    upper:  &CSM<usize, SnzVal>,
    bb:     &SparseVector<usize, SnzVal>
) -> SparseVector<usize, SnzVal>{

    let mut hash = HashMap::new();
    let mut heap = BinaryHeap::new();
    for ii in 0..bb.ind.len() {
        hash.insert(bb.ind[ii].clone(), bb.snz[ii].clone());
        heap.push(Reverse(bb.ind[ii].clone()));
    }

    let mut xx: SparseVector<usize, SnzVal> = SparseVector::new();

    while let Some(Reverse(index)) = heap.pop() {
        if let Some(value) = hash.remove(&index) {
            if !upper.ringmetadata.is_0(&value) {
                xx.ind.push(index.clone());
                xx.snz.push(value.clone());
                let mut row = upper.maj_hash(&index);
                row.remove(&index);
                let scale = -value;
                update_heap_hash(&upper.ringmetadata, &mut heap, &mut hash, &mut row, &scale);
            }
        }
    }
    return xx;
}

/// Solve upper*xx = bb, solution is indexed by new key.
///
/// # Parameters
/// - `uppper`: a non-singular upper triangular matrix in CSM format
/// - `bb`: constant column sparse vector
/// - `majind_pivot_minind`: a vector indicates pivot entry locations
/// # Returns
/// xx the unknow sparse vector
pub fn triangular_solver_version_2<SnzVal>(
    upper:                  &CSM<usize, SnzVal>, // column major csm
    majind_pivot_minind:    &Vec<usize>,
    mut bb:                 &mut HashMap<usize, SnzVal>
) -> HashMap<usize, SnzVal> where
SnzVal: Neg + Neg<Output=SnzVal> + AddAssign + Mul + Mul<Output=SnzVal> + Clone + InvMod + InvMod<Output = SnzVal> + PartialEq + Debug
{
    let mut xx: HashMap<usize, SnzVal> = HashMap::new();
    let mut ii = upper.nummaj;
    while ii > 0 {
        if let Some(value) = bb.remove(&majind_pivot_minind[ii-1]) {
            let mut hash = upper.maj_hash(&(ii-1));
            if let Some(dominator) = hash.remove(&majind_pivot_minind[ii-1]) {
                if let Some(inverse) = upper.ring().inverse(&dominator){
                    let mut scale = value*inverse;
                    scale = upper.ring().simplify(&scale);
                    let neg_scale = -scale.clone();
                    add_assign_hash(&upper.ring(), &mut bb, &mut hash, &neg_scale);
                    xx.insert(ii-1, scale); // solution is indexed by new key after ordering.
                }
            }
        }
        ii -= 1;
    }
    return xx;
}

/// Update the hash by adding scaled row to it
///
/// # Parameters
/// -`ringmetadata`: ring data of entries
/// -`hash`: the row to be updated
/// -`row`: a row represented as a hash map
/// - `scale`: the scale value
pub fn add_assign_hash<MinKey, SnzVal>(
    ringmetadata:   &RingMetadata<SnzVal>,
    hash:           &mut HashMap<MinKey, SnzVal>,
    row:            &mut HashMap<MinKey, SnzVal>,
    scale:          &SnzVal
) where
MinKey: Eq+Hash+Ord+Clone,
SnzVal: Clone + AddAssign + Mul<Output = SnzVal> + PartialEq + InvMod<Output = SnzVal>
{
    for (key, val) in row.drain() {
        let value = scale.clone()*val;
        if let Some(x) = hash.get_mut(&key) {
            *x += value;
            if ringmetadata.is_0(x) {
                hash.remove(&key);
            }
        }
        else if !ringmetadata.is_0(&value){
            hash.insert(key, value);
        }
    }
}

/// Multiply a row vector (represented as a hash map) with a smoracle matrix
///
/// # Parameters
/// -`hash`: row vector represented as a hash map
/// -`matrix`: the sparse matrix oracle
/// # Returns
/// the product vector represented as a hash map
pub fn multiply_hash_smoracle<MajKey, MinKey, SnzVal, Matrix>(
    hash:     &HashMap<MajKey, SnzVal>,
    matrix:   &Matrix
) -> HashMap<MinKey, SnzVal>
where
MinKey: PartialEq + Eq + Hash + Clone + Ord,
MajKey: PartialEq + Eq + Hash + Clone + Ord,
SnzVal: Clone + AddAssign + Mul<Output = SnzVal> + PartialEq + InvMod<Output = SnzVal>,
Matrix: SmOracle<MajKey, MinKey, SnzVal>
{
    let mut product = HashMap::new();
    for (key, val) in hash.iter() {
        let mut row = matrix.maj_hash(key);
        add_assign_hash(&matrix.ring(), &mut product, &mut row, &val);
    }
    return product;
}

/// Multiply a row vector (represented as a hash map) with a sparse matrix oracle. In this version, hash key is usize.
///
/// # Parameters
/// -`hash`: row vector represented as a hash map
/// -`matrix`: the sparse matrix oracle
/// -`index_pivot_majkey`: a vector transform the each key (which is usize type) of hash to major key type of the matrix
/// # Returns
/// the product vector represented as a hash map
pub fn multiply_hash_smoracle_version2<MajKey, MinKey, SnzVal, Matrix>(
    hash:                   &HashMap<usize, SnzVal>,
    index_pivot_majkey:     &Vec<MajKey>,
    matrix:                 &Matrix
) -> HashMap<MinKey, SnzVal>
where
MinKey: PartialEq + Eq + Hash + Clone + Ord,
MajKey: PartialEq + Eq + Hash + Clone + Ord,
SnzVal: Clone + AddAssign + Mul<Output = SnzVal> + PartialEq + InvMod<Output = SnzVal>,
Matrix: SmOracle<MajKey, MinKey, SnzVal>
{
    let mut product = HashMap::new();
    for (majind, snzval) in hash.iter() {
        let mut row = matrix.maj_hash(&index_pivot_majkey[*majind]);
        add_assign_hash(&matrix.ring(), &mut product, &mut row, &snzval);
    }
    return product;
}

/// Solve matrix*xx = hash
///
/// # Parameters
/// - `matrix`: a non-singular pivot sparse matrix oracle
/// - `hash`: constant column sparse vector represented as a hash map
/// - `minkey_pivot_majkey`:  A hash map indicates the pivots locations
/// # Returns
/// xx the unknow sparse vector represented as a hash map
pub fn solver<MinKey, MajKey, SnzVal, Matrix>(
    matrix:                 &Matrix,
    minkey_pivot_majkey:    &HashMap<MinKey, MajKey>,
    mut hash:               &mut HashMap<MinKey, SnzVal>
) -> HashMap<MajKey, SnzVal> where
MinKey: PartialEq + Eq + Hash + Clone + Ord,
MajKey: PartialEq + Eq + Hash + Clone + Ord,
SnzVal: Neg + Neg<Output=SnzVal> + AddAssign + Mul + Mul<Output=SnzVal> + Clone + InvMod + InvMod<Output = SnzVal> + PartialEq,
Matrix: SmOracle<MajKey, MinKey, SnzVal>
{
    let mut heap: BinaryHeap<Reverse<MinKey>> = BinaryHeap::new();
    for (key, _val) in hash.iter() {
        heap.push(Reverse(key.clone()));
    }

    let mut xx: HashMap<MajKey, SnzVal> = HashMap::new();

    while let Some(Reverse(minkey)) = heap.pop() {
        //if minkey_pivot_majkey.contains_key(&minkey) {
            let majkey = &minkey_pivot_majkey[&minkey];
            if let Some(value) = hash.remove(&minkey) {
                let mut row: HashMap<MinKey, SnzVal> = matrix.maj_hash(majkey);
                if let Some(dominator) = row.remove(&minkey) {
                    if let Some(inverse) = matrix.ring().inverse(&dominator){
                        let mut scale = -value.clone()*inverse.clone();
                        scale = matrix.ring().simplify(&scale);
                        update_heap_hash(&matrix.ring(), &mut heap, &mut hash, &mut row, &scale);
                        xx.insert(majkey.clone(), value*inverse);
                    }
                }
            }
        //}
    }
    return xx;
}


//pub fn submatrix_of_smoracle<MajKey, MinKey, SnzVal, Matrix>(matrix: Matrix, maj_table: Vec<MajKey>, min_table: Vec<MinKey>) -> CSM<usize, SnzVal> { }
