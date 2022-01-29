use std::hash::{Hash, Hasher};
use priority_queue::PriorityQueue;
use std::collections::HashMap;
use crate::matrix::SmOracle;
use crate::csm::CSM;

/// The key for priority queue
struct PqKey<SnzVal> {
    id: usize,
    value: SnzVal,
}

impl<D> Hash for PqKey<D> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl<D> PartialEq for PqKey<D> {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for PqKey {}

/// Append the priority queue by scaled indicated row of a matrix
///
/// # Parameters
/// - `pq`: A priority queue
/// - `matrix`: A sparse matrix
/// - `row_ind`: The index of indicated row
/// - `scale`: The multiplication scale
///
/// # Returns
/// - An appended priority queue
fn pq_append<MinKey, MajKey, SnzVal, M: SmOracle<MajKey, MinKey, SnzVal>>(pq:&mut PriorityQueue<PqKey<SnzVal>, usize>, matrix:&M, row_ind:usize, majind: &Vec<MajKey>, scale:SnzVal, add_lead:bool) {
    let (ind_vec,val_vec) = matrix.maj_sprsvec(majind[row_ind]);
    let mut aa:usize = 0;
    if add_lead == false { aa=1; }
    for jj in aa..val_vec.len() {
        let pqkey = PqKey{id:matrix.rowptr[row_ind]+jj, value:scale*val_vec[jj]};
        let pqval = matrix.min_length()-ind_vec[jj];
        pq.push(pqkey, pqval);
    }
}

/// Append the priority queue by scaled given index vector and value vector
///
/// # Parameters
/// - `pq`: A priority queue
/// - `col_vec`: The vector of colum indices
/// - `val_vec`: The vector of entry values
/// - `scale`: The scale
/// # Returns
/// - An appended priority queue
fn pq_append_toclear<M: SmOracle<MajKey, MinKey, SnzVal>>(pq:&mut PriorityQueue<PqKey, usize>, matrix: &M, reduced:&M, row_ind:usize) {
    let (col_vec, val_vec) = matrix.row(row_ind);
    for jj in 0..val_vec.len() {
        let pqkey = PqKey{id:reduced.colind.len()+jj, value:val_vec[jj]};
        let pqval = reduced.min_length()-col_vec[jj];
        pq.push(pqkey, pqval);
    }
}

/// UU decompose of a sparse matrix
/// # Parameters
/// - `matrix`: a sparse matrix in compressed row format.
///
/// # Returns
/// - A sparse matrix in compressed row format representing the reduced matrix.
/// - A sparse matrix in compressed row format representing the row operation.
/// - A hash table recording pivot row index of each pivot colum index.
pub fn cohomology<MajKey, MinKey, SnzVal, Matrix: SmOracle<MajKey, MinKey, SnzVal>>(
    matrix:             &Matrix,
    maj_to_reduce:      &Vec<MajKey>
) -> (CSM<usize, SnzVal>, Indexing<MinKey, MajKey>) where
MajKey: PartialEq + Eq + Hash + Clone + Ord + Debug,
MinKey: PartialEq + Eq + Hash + Clone + Ord + Debug,
SnzVal: Add + Neg<Output=SnzVal> +Clone + PartialEq + InvMod<Output=SnzVal> + Mul<Output=SnzVal> + AddAssign + Debug,
Matrix: SmOracle<MajKey, MinKey, SnzVal>
{

	//initialize "reduced", "row_oper", "row_oper_inv" and "permut"
    let nummaj = matrix.maj_length();
    let nummin = matrix.min_length();
	let mut reduced = CSM::new(nummaj, nummin);
	let mut maj_oper = CSM::new(numamj, nummaj);
	let mut maj_oper_inv = CSM::new(numamj, nummaj);
	let mut pivot_hash = HashMap::new();

	//row elimination
    let mut pq_for_reduced = PriorityQueue::new(); // add some comments about pq
    let mut pq_for_row_oper = PriorityQueue::new();
	for ii in 0..nummaj {
		let row_index = nummaj-1-ii;
		row_oper.push(row_index, 1.0);             // push_snzval
		row_oper_inv.push(row_index, 1.0);

        if matrix.maj_itr(majind[row_index]) == None {
            reduced.rowptr.push(reduced.minind.len());
            row_oper.rowptr.push(row_oper.minind.len());
            row_oper_inv.rowptr.push(row_oper_inv.minind.len());
            continue;
        }

		pq_for_reduced.clear();
		pq_for_reduced.push(PqKey{id:reduced.minind.len()+nummaj, value:0.0}, 0);
        pq_append_toclear(&mut pq_for_reduced, &matrix, &reduced, row_index);

		pq_for_row_oper.clear();
		pq_for_row_oper.push(PqKey{id:row_oper.colind.len()+nummaj, value:0.0}, 0);

		let mut lead_entry = 0.0; let mut priority = 0;
		if let Some((key, pval)) = pq_for_reduced.pop(){
			lead_entry = key.value; priority = pval;
		}

		let mut look_up_table = true;
		while let Some((pqkey,pqval)) = pq_for_reduced.pop() {
			if pqval == priority {
				lead_entry += pqkey.value;
			} else {
				if lead_entry != 0.0 && look_up_table == true {
                    let col_piv_ind = nummaj-priority;
					match permut.get(&col_piv_ind){
						Some(&row_piv_ind) => {
                            let scale = -lead_entry/reduced.snzval[nummaj-1-row_piv_ind];
                            pq_append(&mut pq_for_reduced, &reduced, nummaj-1-row_piv_ind, scale, false);
							row_oper_inv.push(row_piv_ind, scale);
                            pq_append(&mut pq_for_row_oper, &row_oper, nummaj-1-row_piv_ind, scale, true);
						}
						None => {
                            permut.insert(nummin-priority, row_index);
                            look_up_table = false;
                            reduced.push(nummin-priority, lead_entry);
                        }
					}
				} else if lead_entry != 0.0 && look_up_table == false {
                    reduced.push(nummin-priority, lead_entry);
                }
				lead_entry = pqkey.value; priority = pqval;
			}
        }
        reduced.rowptr.push(reduced.colind.len());
        row_oper_inv.rowptr.push(row_oper_inv.colind.len());

        let mut lead_entry = 0.0;
        let mut priority = 0;
        while let Some((key,pval))=pq_for_row_oper.pop() {
        		if pval==priority { lead_entry = lead_entry+key.value; }
        		else {
            		if lead_entry != 0.0 { row_oper.push(nummaj-priority,lead_entry); }
            		lead_entry = key.value; priority=pval;
        		}
        }
        row_oper.rowptr.push(row_oper.colind.len());
	}
    let reduced = reduced.reverse_row_order();
    let row_oper = row_oper.reverse_row_order();
	return (reduced,row_oper,permut);
}
