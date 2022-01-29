/*!

Filtered clique complexes (also known as a *Vietoris-Rips* complexes)


# Clique chain complex oracles


Mathematically, a filtered clique complex is determined by two pieces of data:

* a **dissimilarity matrix**
    * we represent this by a symmetric matrix *S* that has been flattened into a vector *v*
* a **threshold**  that determines where we stop growing the filtration
    * we represent this as a real number *t*


The boundary matrices of this chain complex oracle are indexed by `Simplex` objects which are
(essentially) pairs of form `(simplex, filtration_value)`.
```
/*
Example:
    - construct a dissimilarity matrix
    - construct the associated clique complex
    - access a row + column of the boundary matrix
    - compute the dimension 1 barcode

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
*/

```
*/



use std::ops::Neg;
use std::convert::TryInto;
use std::hash::Hash;
use std::fmt::Debug;

use crate::matrix::{SmOracle, RingMetadata, MajorDimension};
use crate::chx::{ChainComplex, ChxTransformKind};

use std::cmp::Ordering;

// The specific type we use to store vertices
type Vertex = u16;

/// A simplex associated with a filtration value
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct Simplex<FilVal: Clone + Debug>{
    pub filvalue: FilVal,     // the filtration value of simplex
    pub vertices: Vec<Vertex> // a sorted vector in strictly descending order recording vertices
}

impl <FilVal>
    PartialOrd for Simplex<FilVal> 
    where FilVal: Eq + PartialEq + Ord + PartialOrd + Clone + Debug

{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl <FilVal>
    Ord for Simplex<FilVal> 
    where FilVal: Eq + PartialEq + Ord + PartialOrd + Clone + Debug
{
    fn cmp(&self, other: &Self) -> Ordering {

        // first compare filtration values
        let mut comp = self.filvalue.cmp( & other.filvalue );
        if comp != Ordering::Equal { return comp }

        // next compare simplex dimensions
        comp = self.vertices.len().cmp( & other.vertices.len() );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        return self.vertices.cmp( & other.vertices )
    }
}


/// Boundary matrix represented as an [`SmOracle`] trait object.
pub struct CliqueBoundaryMatrix<FilVal, SnzVal: Clone>{
    pub ringmetadata: RingMetadata<SnzVal>, // ring meta data
    /// The "maxdis" value represents the maximum of filtration value
    pub maxdis: FilVal,
    /// A vector representing the dissimilarity matrix by listing all its rows
	pub dismat: Vec<Vec<FilVal>>,
    /// A vector representing the neighborhood within "maxdis" of each vertex
    pub cutoff_matrix: Vec<Vec<Vertex>>
}

/// A function that produces "cutoff_matrix"
///
/// # Parameters
/// - `dismat`: A dissimilarity/distance matrix
/// - `maxdis`: The radius of neighborhood
fn cutoff_matrix<FilVal: PartialOrd + Debug>(dismat: &Vec<Vec<FilVal>>, maxdis: FilVal) -> Vec<Vec<Vertex>> {
    let mut output: Vec<Vec<Vertex>> = Vec::new();
    for row in dismat.iter() {
        let mut index_vec = Vec::new();
        for jj in 0..row.len() {
            if row[jj] <= maxdis {
                index_vec.push(jj.try_into().unwrap());
            }
        }
        output.push(index_vec);
    }
    return output;
}

/// Methods of CliqueBoundaryMatrix struct
impl<FilVal: Clone + Debug + PartialOrd + Ord, SnzVal: Clone> CliqueBoundaryMatrix<FilVal, SnzVal> {
    /// Output the diameter of a simplex
    ///
    /// # Parameters
    /// -`vertices`: vertices of the simplex
    /// # Returns
    /// Return None if some edge of simplex is out of cutoff value; return Some(diam) otherwise
    fn diam(&self, vertices: &Vec<Vertex>) -> Option<FilVal> {
        if vertices.is_empty() { return None };

        let mut a = usize::from( vertices[0] );     // first vertex
        let mut b = a.clone();                      // second vertex

        let mut diam_bound = &self.dismat[ a ][ b ];   // lower bound on diameter
        let mut diam = diam_bound.clone(); 
        for ii in 0..vertices.len() {
            a = usize::from( vertices[ii] ); 
            for jj in 0..ii {
                b = usize::from( vertices[jj] ); 
                diam_bound = &self.dismat[ a ][ b ];
                if diam_bound > &diam { diam = diam_bound.clone(); }
            }
        }
        return Some(diam);
    }

}

/// Find the largest minimum index of a vector
fn find_min<FilVal: PartialOrd>(vector: Vec<FilVal>) -> usize {
    let mut index = 0;
    for ii in 0..vector.len() {
        if vector[ii] < vector[index] {
            index = ii;
        }
    }
    return index;
}

/// An iterator that iterates all cofacets of a given simplex
struct CofacetIter<'a, FilVal: Clone, SnzVal: Clone> {
	pub clique: &'a CliqueBoundaryMatrix<FilVal, SnzVal>,
	pub next_cofacet_vertices: Vec<Vertex>,
    pub simplex_filvalue: FilVal,
	pub insertion_location: usize,
    pub candidate_location: usize,
    pub first_vert: Vertex,
	pub coeff: SnzVal,
}

/// implement standard methods of Iterator for CofacetIter struct
impl<'a, FilVal, SnzVal> Iterator for CofacetIter<'a, FilVal, SnzVal> where
FilVal: Clone + Debug + PartialOrd,
SnzVal: Clone + Neg<Output = SnzVal>
{
	type Item = (Simplex<FilVal>, SnzVal);

    fn next( &mut self ) -> Option<Self::Item> {

        let candidacies = &self.clique.cutoff_matrix[usize::from(self.first_vert)];
        let vertices = &mut self.next_cofacet_vertices;

        loop {
            if self.candidate_location >= candidacies.len(){
                return None;
            }
            let new_vertex = candidacies[self.candidate_location];
            let mut flag = false;
            vertices[self.insertion_location] = new_vertex;
            let mut new_filval = self.simplex_filvalue.clone();
            for vert in vertices.iter() {
                if self.clique.dismat[usize::from(new_vertex)][usize::from(*vert)] > self.clique.maxdis {
                    flag = true;
                    break;
                } else if self.clique.dismat[usize::from(new_vertex)][usize::from(*vert)] > new_filval {
                    new_filval = self.clique.dismat[usize::from(new_vertex)][usize::from(*vert)].clone();
                }
            }
            if flag { self.candidate_location += 1; continue; }

            while self.insertion_location < vertices.len()-1 &&
                                new_vertex >= vertices[self.insertion_location+1] {
                if new_vertex == vertices[self.insertion_location+1] {
                    flag = true;
                    break;
                }
                vertices[self.insertion_location] = vertices[self.insertion_location+1];
                self.insertion_location += 1;
                self.coeff = -self.coeff.clone();
            }
            if flag { self.candidate_location += 1; continue; }

            if self.insertion_location >= vertices.len() {
                return None;
            }
            vertices[self.insertion_location] = new_vertex;
        //println!("new_filval={:?}, simplex_filvalue={:?}", new_filval, self.simplex_filvalue);
            let simp = Simplex{ vertices: vertices.clone(), filvalue: new_filval };
            self.candidate_location += 1;
            return Some((simp, self.coeff.clone()));
        }
        //println!("candidate_location={:?}", self.candidate_location);
	}
}

/// A iterator that iterates all facets of a given simplex
struct FacetIter<'a, FilVal: Clone + Debug, SnzVal: Clone> {
    pub clique: &'a CliqueBoundaryMatrix<FilVal, SnzVal>,
	pub simp: Simplex<FilVal>,
	pub removal_location: usize,
	pub coeff: SnzVal
}

/// implement standard methods of Iterator for FacetIter struct
impl<FilVal, SnzVal> Iterator for FacetIter<'_, FilVal, SnzVal> where
FilVal: Clone + PartialOrd + Debug + Ord,
SnzVal: Clone + Neg<Output = SnzVal>
{
	type Item = (Simplex<FilVal>, SnzVal);

	fn next( &mut self ) -> Option<Self::Item> {
        if self.removal_location == self.simp.vertices.len() { return None; }
        let mut simplex = self.simp.clone();
        simplex.vertices.remove(self.removal_location);
        if let Some(diam) = self.clique.diam(&simplex.vertices) {
            simplex.filvalue = diam;
            self.coeff = -self.coeff.clone();
            self.removal_location += 1;
            return Some((simplex, self.coeff.clone()));
        } else {
            return None;
        }
    }
}

// See matrix.rs file for specific definition of SmOracle trait
impl<FilVal, SnzVal> SmOracle<Simplex<FilVal>, Simplex<FilVal>, SnzVal> for CliqueBoundaryMatrix<FilVal, SnzVal> where
FilVal: PartialOrd + Clone + Debug + Eq + Hash + Ord,
SnzVal: Neg<Output = SnzVal> + Clone
{

    fn ring( &self ) -> &RingMetadata<SnzVal> {
        &self.ringmetadata
    }

	fn maj_dim( &self ) -> MajorDimension { MajorDimension::Row }

	fn maj_itr(&self, majkey: &Simplex<FilVal>) -> Box<dyn Iterator<Item=(Simplex<FilVal>, SnzVal)> + '_> {
        let mut vertices = majkey.vertices.clone();
        let first_vert = vertices[0];
        vertices.insert(0,0);
        vertices.shrink_to_fit();
        Box::new(CofacetIter{
        	clique: &self,
        	next_cofacet_vertices: vertices,
            simplex_filvalue: majkey.filvalue.clone(),
        	insertion_location: 0,
            candidate_location: 0,
        	coeff: self.ringmetadata.identity_multiplicative.clone(),
            first_vert,
        })
	}

	fn min_itr( &self, minkey: &Simplex<FilVal> ) -> Box<dyn Iterator<Item=(Simplex<FilVal>, SnzVal)> + '_> {
        Box::new( FacetIter {
            clique: &self,
            simp: minkey.clone(),
            removal_location: 0,
            coeff: -self.ringmetadata.identity_multiplicative.clone()
        })
	}

	// Total # of nonzeros in the matrix (if this number can actually be calculated)
    // Fill later
	fn countsnz( &self ) -> Option<usize> {
        //for
        //for iter in self.maj_itr()
		None
	}

	fn finiteminors( &self ) -> Option<bool> {
		Some(true)
	}

	fn finitemajors( &self ) -> Option<bool> {
		Some(true)
	}

    fn is_pivot(&self, majkey: &Simplex<FilVal>, minkey: &Simplex<FilVal>) -> Option<bool> {
        if majkey.filvalue == minkey.filvalue && majkey.vertices[0] < minkey.vertices[0]{
            return Some(true);
        }
        else {
            return Some(false);
        }
    }

}


/// Based, filtered chain complex implementing the [`ChainComplex`](ChainComplex) trait
pub struct CliqueComplex<SnzVal: Clone, FilVal> {
    pub ringmetadata: RingMetadata<SnzVal>,
    pub dissimilarity_matrix: Vec<Vec<FilVal>>,
    pub dissimilarity_value_max: FilVal,
    pub safe_homology_degrees_to_build_boundaries: Vec<usize>,
    pub major_dimension: MajorDimension,
    pub simplex_count: Vec<(usize,usize)>
}

/// An iterator that iterates all simplices of given dimension
struct SimplexIter<'a, FilVal: Clone + Debug, SnzVal: Clone> {
    pub chx: &'a CliqueComplex<SnzVal, FilVal>,
    pub filvec: Vec<FilVal>,
    pub vec: Vec<Vertex>,
    pub val: Vertex,
    pub loc: usize
}
/// implement standard methods of Iterator for SimplexIter struct
impl<'a, FilVal, SnzVal> Iterator for SimplexIter<'a, FilVal, SnzVal> where
FilVal: Clone + PartialOrd + Debug,
SnzVal: Clone
{
    type Item = Simplex<FilVal>;

    fn next(&mut self) -> Option<Self::Item> {
        let size = self.chx.dissimilarity_matrix.len();
        if self.vec.len() > size { return None; }
        loop {
            while usize::from(self.val) <= size-self.vec.len()+self.loc {
                self.vec[self.loc] = self.val;
                self.filvec[self.loc+1] = self.filvec[self.loc].clone();
                for ii in 0..self.loc {
                    if self.filvec[self.loc+1] < self.chx.dissimilarity_matrix[usize::from(self.val)][usize::from(self.vec[ii])] {
                        self.filvec[self.loc+1] = self.chx.dissimilarity_matrix[usize::from(self.val)][usize::from(self.vec[ii])].clone();
                    }
                }
                if self.filvec[self.loc+1] <= self.chx.dissimilarity_value_max {
                    if self.loc == self.vec.len()-1 {
                        self.val += 1;
                        return Some(Simplex{
                                vertices: self.vec.clone(),
                                filvalue: self.filvec[self.loc+1].clone()
                            });
                    } else {
                        self.loc += 1;
                    }
                }
                self.val += 1;
            }
            if self.loc == 0 { return None; }
            self.loc = 	self.loc - 1;
            self.val =	self.vec[self.loc] + 1;

        }
    }
}

// Implimentation of ChainComplex trait. See chx.rs for definition of ChainComplex trait.
impl<SnzVal, FilVal> ChainComplex<Simplex<FilVal>, SnzVal, FilVal> for CliqueComplex<SnzVal, FilVal> where
SnzVal: Clone + Neg<Output = SnzVal> + Debug,
FilVal: PartialOrd + Copy + Debug + Eq + Hash + Ord
{
    type Matrix = CliqueBoundaryMatrix<FilVal, SnzVal>;

    fn get_smoracle(
        &self,
        major_dim       :   MajorDimension, // row or column
        transform       :   ChxTransformKind
    ) -> CliqueBoundaryMatrix<FilVal, SnzVal> {
        CliqueBoundaryMatrix {
            ringmetadata: self.ringmetadata.clone(),
            maxdis: self.dissimilarity_value_max,
        	dismat: self.dissimilarity_matrix.clone(),
            cutoff_matrix: cutoff_matrix(&self.dissimilarity_matrix, self.dissimilarity_value_max)
        }
    }

    /// Return an iterator of all simplices for a given dimension
    ///
    /// # Parameters
    /// -`dim`: The given dimension
    fn keys_unordered_itr(&self, dim:usize) -> Box<dyn Iterator<Item=Simplex<FilVal>> + '_> {
        Box::new( SimplexIter {
            chx: &self,
            filvec: vec![self.dissimilarity_matrix[0][0].clone(); dim+2],
            vec: vec![0; dim+1],
            val: 0,
            loc: 0
        })
    }

    fn keys_ordered(&self, dim: usize) -> Vec<Simplex<FilVal>> {
        let mut vector: Vec<Simplex<FilVal>> = Vec::new();

        for simp in self.keys_unordered_itr(dim) {
            vector.push(simp.clone());
        }

        vector.sort();
        return vector;
    }

    fn key_2_filtration(&self, key: &Simplex<FilVal>) -> FilVal {
        return key.filvalue;
    }

    fn max_filtration(&self) -> FilVal {
        return self.dissimilarity_value_max;
    }
}
