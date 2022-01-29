

/*!

UU factorization of a sparse matrix oracle

*/


use std::collections::BinaryHeap;
use std::cmp::Reverse;
use std::hash::Hash;
use std::fmt::Debug;
use std::collections::{HashMap, HashSet};
use std::ops::{Add, Neg, AddAssign, Mul};

use crate::matrix::{SmOracle, InvMod, MajorDimension, RingMetadata};
use crate::csm::CSM;
use crate::chx::Indexing;

use crate::solver::{add_assign_hash, multiply_hash_smoracle_version2};

/// Scale a given row and add it to a sparse vector (represented as a hash table  associated with a binary heap which record the priorities of each entry)
///
/// # Parameters
/// - `ringmetadata`: The coefficient ring information
/// - `hash`: A hash map representing a sparse row we want to update
/// - `heap`: A heap that indicates the priority of each entry in the sparse row we want to update
/// - `row`: The row we want to scale and add to spare vector (hash and heap)
/// - `scale`: The scale value
fn update_heap_hash<MinKey, SnzVal>(
    ringmetadata:   & RingMetadata<SnzVal>,
    heap:           &mut BinaryHeap<MinKey>,
    hash:           &mut HashMap<MinKey, SnzVal>,
    row:            &mut HashMap<MinKey, SnzVal>,
    scale: &SnzVal) where
MinKey: Eq + Hash + Ord + Clone,
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
            heap.push(key.clone());
            hash.insert(key, value);
        }
    }
}

/// UU row(maj) decomposition of a sparse matrix with selected rows and columns
/// # Parameters
/// - `matrix`: a sparse matrix oracle
/// - `maj_to_reduce`: an vector indicates the rows(majs) and the order to perform decomposition
/// - `min_to_reduce`: an hashset indicates the cols(mins) to perform decomposition
///
/// # Returns
/// - A sparse matrix in compressed row format representing the row operation.
/// - A Indexing recording pivot matrix, order of pivot major indices, order of pivot minor indices.
pub fn decomp_col_use_pairs<MajKey, MinKey, SnzVal, Matrix>(
    matrix:             &Matrix,
    maj_to_reduce:      &mut Vec<MajKey>,
    min_to_reduce:      &mut HashSet<MinKey>
) -> (CSM<usize, SnzVal>, Indexing<MinKey, MajKey>) where
MajKey: PartialEq + Eq + Hash + Clone + Ord + Debug,
MinKey: PartialEq + Eq + Hash + Clone + Ord + Debug,
SnzVal: Add + Neg<Output=SnzVal> +Clone + PartialEq + InvMod<Output=SnzVal> + Mul<Output=SnzVal> + AddAssign + Debug,
Matrix: SmOracle<MajKey, MinKey, SnzVal>
{
    //initialize "rowoper" and "indexing"
    let length:f64 = maj_to_reduce.len() as f64;
    let capacity:usize = (length*1.2) as usize;
    //let capacity: usize = maj_to_reduce.len();
	let mut rowoper: CSM<usize, SnzVal> = CSM::with_capacity(capacity, MajorDimension::Row, matrix.ring().clone());

    let mut indexing = Indexing::with_capacity(capacity);

    let mut heap_reduced = BinaryHeap::new();   //with_capacity(capacity);
    let mut heap_rowoper = BinaryHeap::new();   //with_capacity(capacity);
    let mut hash_reduced = HashMap::new();      //with_capacity(capacity);
    let mut hash_rowoper = HashMap::new();      //with_capacity(capacity);

    let mut row_reduced = HashMap::new();       //with_capacity(capacity);
    let mut row_rowoper = HashMap::new();       //with_capacity(capacity);

    let one: SnzVal = matrix.ring().identity_multiplicative.clone();

    //eliminate cols
	while let Some(majkey) = maj_to_reduce.pop() {

        heap_reduced.clear();
        heap_rowoper.clear();
        hash_reduced.clear();
        hash_rowoper.clear();

        for (key, val) in matrix.maj_itr(&majkey) {
            heap_reduced.push(key.clone());
            hash_reduced.insert(key,val);
        }

		while let Some(minkey) = heap_reduced.pop(){
            if let Some(leading_entry) = hash_reduced.remove(&minkey){
                if matrix.ring().is_0(&leading_entry) { continue; }
                else if indexing.minkey_2_index.contains_key(&minkey){
                    let index = indexing.minkey_2_index[&minkey];
                    row_rowoper.clear();
                    row_reduced.clear();
                    row_rowoper = rowoper.maj_hash(&index);

                    //row_reduced = multiply_hash_smoracle_version2(&row_rowoper, &indexing.index_2_majkey, matrix);

                    let mut row = HashMap::new();
                    for (majind, snzval) in row_rowoper.iter() {
                        row.clear();
                        for (key, val) in matrix.maj_itr(&indexing.index_2_majkey[*majind]){
                            if min_to_reduce.contains(&key) {
                                row.insert(key, val);
                            }
                        }
                        add_assign_hash(&matrix.ring(), &mut row_reduced, &mut row, &snzval);
                    }

                    if let Some(dominator) = row_reduced.remove(&minkey) {
                        if let Some(inverse) = matrix.ring().inverse(&dominator){
                            let mut scale = -leading_entry*inverse;
                            scale = matrix.ring().simplify(&scale);
                            update_heap_hash(&matrix.ring(), &mut heap_reduced, &mut hash_reduced, &mut row_reduced, &scale);
                            update_heap_hash(&matrix.ring(), &mut heap_rowoper, &mut hash_rowoper, &mut row_rowoper, &scale);
                        }
                    }
                } else {
                    indexing.minkey_2_index.insert(minkey.clone(), rowoper.nummaj);
                    indexing.majkey_2_index.insert(majkey.clone(), rowoper.nummaj);
                    indexing.index_2_majkey.push(majkey);
                    indexing.index_2_minkey.push(minkey);

                    rowoper.push_snzval(rowoper.nummaj, one.clone());
                    rowoper.append_maj(&mut hash_rowoper);
                    break;
                }
            }
        }
	}

    heap_reduced.clear();
    for minkey in indexing.index_2_minkey.iter() {
        heap_reduced.push(minkey.clone());
    }
    while let Some(minkey) = heap_reduced.pop(){
        indexing.ordered_minind.push(indexing.minkey_2_index[&minkey]);
    }

    rowoper.shrink_to_fit();
    indexing.shrink_to_fit();
    return (rowoper, indexing);
}
