
/*!

Filtered cubical complexes

# Cubical chain complex oracles

We use the [T-construction](https://arxiv.org/pdf/2005.04597.pdf) to turn any n-dimensional array into a filtered cubical complex.

```
Example:
    - use Haibin's excellent test example // it was something like this
           0  5  2  4  3
           0  6  2  7  3
           0  5  2  4  3

    - construct the associated clique complex
    - access a row + column of the boundary
    - compute + print the dimension 1 barcode
```

*/


use std::ops::Neg;
use std::convert::TryInto;
use std::hash::Hash;
use std::fmt::Debug;

use crate::matrix::{SmOracle, RingMetadata, MajorDimension};
use crate::chx::{ChainComplex, ChxTransformKind};

type Coordi = u32;

#[derive(Debug, PartialEq, Clone, PartialOrd, Eq, Ord, Hash)]
pub struct Cube<FilVal: Clone + Debug>{
    pub filvalue: FilVal,           // the filtration value of the cube
    pub coordinates: Vec<Coordi>    // a vector represents the location of the cube
}

// shape [3,4,2]
// index order in entry_array: [0,0,0] < [0,0,1] < [0,1,0] < [0,1,1] < [0,2,0] < ...
pub struct CubicalBoundaryMatrix<FilVal, SnzVal: Clone> {
    pub ringmetadata: RingMetadata<SnzVal>, // ring meta data
    pub majordimension: MajorDimension,
    pub shape: Vec<usize>,            // shape of the n-dimensional matrix
    pub entry_array: Vec<FilVal>, // entries from n-dimensional matrix in lexicographical order
    pub max_value: FilVal
}

impl<FilVal: PartialOrd + Clone, SnzVal: Clone> CubicalBoundaryMatrix<FilVal, SnzVal> {
    /// Transform the top dim cube coordinates to index in 1-d vector "entry_array".
    pub fn coordinates_to_index(&self, coordinates: &Vec<Coordi>) -> usize {
        // For example, if shape= [3,2,5], then product = [30,10,5,1].
        let mut loc = self.shape.len();
        let mut product = vec![1; loc+1];

        while loc > 0 {
            loc -= 1;
            product[loc] = product[loc+1]*(self.shape[loc] as Coordi);
        }

        let mut index = 0;
        for loc in 0..coordinates.len() {
            index += (coordinates[loc]-1)/2*product[loc+1];
        }
        return index as usize;
    }

    pub fn filvalue(&self, coordinates: &Vec<Coordi>) -> Option<FilVal> {
        let mut upper = Vec::new();
        let mut lower = Vec::new();
        let mut locations = Vec::new();

        for loc in 0..coordinates.len() {
            let coordinate = coordinates[loc];
            if coordinate > 2*self.shape[loc] as Coordi{
                return None;
            } else if coordinate == 2*self.shape[loc] as Coordi{
                upper.push(coordinate-1);
                lower.push(coordinate-1);
            } else if coordinate == 0 {
                upper.push(coordinate+1);
                lower.push(coordinate+1);
            } else if coordinate%2 == 0 {
                upper.push(coordinate+1);
                lower.push(coordinate-1);
                locations.push(loc);
            } else {
                upper.push(coordinate);
                lower.push(coordinate);
            }
        }

        let mut current = lower.clone();
        let mut filvalue = self.entry_array[self.coordinates_to_index(&upper)].clone();
        while current < upper {
            let index = self.coordinates_to_index(&current);
            if filvalue > self.entry_array[index] {
                filvalue = self.entry_array[index].clone();
            }

            for ii in 0..locations.len() {
                let ind = locations.len()-ii-1;
                if current[locations[ind]] == upper[locations[ind]] {
                    current[locations[ind]] = lower[locations[ind]];
                } else {
                    current[locations[ind]] = upper[locations[ind]];
                    break;
                }
            }
        }
        return Some(filvalue);
    }
}

/// An iterator that iterates all cofacets of a given simplex
struct CofacetIter<'a, FilVal: Clone, SnzVal: Clone> {
	pub cubical: &'a CubicalBoundaryMatrix<FilVal, SnzVal>,
    pub current_cube_coordinates: Vec<Coordi>,
    pub edit_location: usize,
    pub just_arrived: bool,
	pub one: SnzVal,
    pub is_last: bool
}

/// implement standard methods of Iterator for CofacetIter struct
impl<'a, FilVal, SnzVal> Iterator for CofacetIter<'a, FilVal, SnzVal> where
FilVal: Clone + Debug + PartialOrd,
SnzVal: Clone + Neg<Output = SnzVal>
{
	type Item = (Cube<FilVal>, SnzVal);

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_last { return None; }

        while self.current_cube_coordinates[self.edit_location]%2 == 1 {
            if self.edit_location == 0 { return None; }
            else {
                self.edit_location -= 1;
                self.just_arrived = true;
            }
        }

        let mut coordinates = self.current_cube_coordinates.clone();
        let mut coeff = -self.one.clone();
        if self.just_arrived {
            if coordinates[self.edit_location] == 0 {
                coordinates[self.edit_location] += 1;
                coeff = self.one.clone();
                if self.edit_location == 0 { self.is_last = true; }
                else {
                    self.edit_location -= 1;
                    self.just_arrived = true;
                }
            } else if coordinates[self.edit_location] == 2*self.cubical.shape[self.edit_location] as Coordi{
                coordinates[self.edit_location] -= 1;
                if self.edit_location == 0 { self.is_last = true; }
                else {
                    self.edit_location -= 1;
                    self.just_arrived = true;
                }
            } else {
                coordinates[self.edit_location] -= 1;
                self.just_arrived = false;
            }
        } else {
            coordinates[self.edit_location] += 1;
            coeff = self.one.clone();
            if self.edit_location == 0 { self.is_last = true; }
            else {
                self.edit_location -= 1;
                self.just_arrived = true;
            }
        }

        if let Some(filvalue) = self.cubical.filvalue(&coordinates) {
            let cube = Cube{ filvalue, coordinates };
            return Some((cube, coeff));
        } else { return None; }

	}
}


/// A iterator that iterates all facets of a given simplex
struct FacetIter<'a, FilVal: Clone + Debug, SnzVal: Clone> {
    pub cubical: &'a CubicalBoundaryMatrix<FilVal, SnzVal>,
    pub current_cube_coordinates: Vec<Coordi>,
    pub edit_location: usize,
    pub just_arrived: bool,
	pub one: SnzVal,
    pub is_last: bool
}


/// implement standard methods of Iterator for FacetIter struct
impl<FilVal, SnzVal> Iterator for FacetIter<'_, FilVal, SnzVal> where
FilVal: Clone + PartialOrd + Debug + Ord,
SnzVal: Clone + Neg<Output = SnzVal>
{
	type Item = (Cube<FilVal>, SnzVal);

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_last { return None; }

        while self.current_cube_coordinates[self.edit_location]%2 == 0 {
            if self.edit_location == 0 { return None; }
            else {
                self.edit_location -= 1;
                self.just_arrived = true;
            }
        }

        let mut coordinates = self.current_cube_coordinates.clone();
        let mut coeff = -self.one.clone();
        if self.just_arrived {
                coordinates[self.edit_location] -= 1;
                self.just_arrived = false;
        } else {
            coordinates[self.edit_location] += 1;
            coeff = self.one.clone();
            if self.edit_location == 0 { self.is_last = true; }
            else {
                self.edit_location -= 1;
                self.just_arrived = true;
            }
        }

        if let Some(filvalue) = self.cubical.filvalue(&coordinates) {
            let cube = Cube{ filvalue, coordinates };
            return Some((cube, coeff));
        } else { return None; }

	}
}


// See matrix.rs file for specific definition of SmOracle trait
impl<FilVal, SnzVal> SmOracle<Cube<FilVal>, Cube<FilVal>, SnzVal> for CubicalBoundaryMatrix<FilVal, SnzVal> where
FilVal: PartialOrd + Clone + Debug + Eq + Hash + Ord,
SnzVal: Neg<Output = SnzVal> + Clone
{

    fn ring(&self) -> &RingMetadata<SnzVal> {
        &self.ringmetadata
    }

	fn maj_itr(&self, majkey: &Cube<FilVal>) -> Box<dyn Iterator<Item=(Cube<FilVal>, SnzVal)> + '_> {

        if self.majordimension == MajorDimension::Row {
            let mut current_cube_coordinates = majkey.coordinates.clone();
            current_cube_coordinates.shrink_to_fit();

            return Box::new(CofacetIter{
                cubical: &self,
                just_arrived: true,
                is_last: false,
                edit_location: current_cube_coordinates.len()-1,
                current_cube_coordinates,
                one: self.ringmetadata.identity_multiplicative.clone()
            });
        } else {
            let mut current_cube_coordinates = majkey.coordinates.clone();
            current_cube_coordinates.shrink_to_fit();

            return Box::new(FacetIter{
                cubical: &self,
                just_arrived: true,
                is_last: false,
                edit_location: current_cube_coordinates.len()-1,
                current_cube_coordinates,
                one: self.ringmetadata.identity_multiplicative.clone()
            });
        }

	}

	fn min_itr( &self, minkey: &Cube<FilVal> ) -> Box<dyn Iterator<Item=(Cube<FilVal>, SnzVal)> + '_> {

        if self.majordimension == MajorDimension::Row {
            let mut current_cube_coordinates = minkey.coordinates.clone();
            current_cube_coordinates.shrink_to_fit();

            return Box::new(FacetIter{
                cubical: &self,
                just_arrived: true,
                is_last: false,
                edit_location: current_cube_coordinates.len()-1,
                current_cube_coordinates,
                one: self.ringmetadata.identity_multiplicative.clone()
            });
        } else {
            let mut current_cube_coordinates = minkey.coordinates.clone();
            current_cube_coordinates.shrink_to_fit();

            return Box::new(CofacetIter{
                cubical: &self,
                just_arrived: true,
                is_last: false,
                edit_location: current_cube_coordinates.len()-1,
                current_cube_coordinates,
                one: self.ringmetadata.identity_multiplicative.clone()
            });
        }
	}
}


/// The CubicalComplex struct which is an instance of the ChainComplex trait
pub struct CubicalComplex<SnzVal: Clone, FilVal> {
    pub ringmetadata: RingMetadata<SnzVal>,
    pub entry_array: Vec<FilVal>,
    pub shape: Vec<usize>,
    pub max_value: FilVal,
    pub safe_homology_degrees_to_build_boundaries: Vec<usize>,
    pub major_dimension: MajorDimension
}

impl<SnzVal: Clone, FilVal: Clone + PartialOrd> CubicalComplex<SnzVal, FilVal> {
    /// Transform the top dim cube coordinates to index in 1-d vector "entry_array".
    pub fn coordinates_to_index(&self, coordinates: &Vec<Coordi>) -> usize {
        // For example, if shape= [3,2,5], then product = [30,10,5,1].
        let mut loc = self.shape.len();
        let mut product = vec![1; loc+1];

        while loc > 0 {
            loc -= 1;
            product[loc] = product[loc+1]*(self.shape[loc] as Coordi);
        }

        let mut index = 0;
        for loc in 0..coordinates.len() {
            index += (coordinates[loc]-1)/2*product[loc+1];
        }
        return index as usize;
    }

    pub fn filvalue(&self, coordinates: &Vec<Coordi>) -> Option<FilVal> {
        let mut upper = Vec::new();
        let mut lower = Vec::new();
        let mut locations = Vec::new();

        for loc in 0..coordinates.len() {
            let coordinate = coordinates[loc];
            if coordinate > 2*self.shape[loc] as Coordi {
                return None;
            } else if coordinate == 2*self.shape[loc] as Coordi {
                upper.push(coordinate-1);
                lower.push(coordinate-1);
            } else if coordinate == 0 {
                upper.push(coordinate+1);
                lower.push(coordinate+1);
            } else if coordinate%2 == 0 {
                upper.push(coordinate+1);
                lower.push(coordinate-1);
                locations.push(loc);
            } else {
                upper.push(coordinate);
                lower.push(coordinate);
            }
        }

        let mut current = lower.clone();
        let mut filvalue = self.entry_array[self.coordinates_to_index(&upper)].clone();
        while current < upper {
            let index = self.coordinates_to_index(&current);
            if filvalue > self.entry_array[index] {
                filvalue = self.entry_array[index].clone();
            }

            for ii in 0..locations.len() {
                let ind = locations.len()-ii-1;
                if current[locations[ind]] == upper[locations[ind]] {
                    current[locations[ind]] = lower[locations[ind]];
                } else {
                    current[locations[ind]] = upper[locations[ind]];
                    break;
                }
            }
        }
        return Some(filvalue);
    }

}


/// An iterator that iterates all cubes of given dimension
struct CubeIter<'a, FilVal: Clone + Debug, SnzVal: Clone> {
    pub chx: &'a CubicalComplex<SnzVal, FilVal>,
    pub coordinates: Vec<Coordi>,
    pub dim: usize,
    pub sum: usize,
    pub first_cube: bool
}

/// implement standard methods of Iterator for CubeIter struct
impl<'a, FilVal, SnzVal> Iterator for CubeIter<'a, FilVal, SnzVal> where
FilVal: Clone + PartialOrd + Debug,
SnzVal: Clone
{
    type Item = Cube<FilVal>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.dim > self.coordinates.len() { return None; }
        let coordinates = &mut self.coordinates;
        let end = coordinates.len()-1;
        if self.first_cube {
            self.first_cube = false;
            if let Some(filvalue) = self.chx.filvalue(coordinates) {
                let cube = Cube{ coordinates:coordinates.clone(), filvalue };
                coordinates[end] += 2;
                return Some(cube);
            } else {
                return None;
            }
        }

        let mut loc = end;
        while coordinates[loc] > 2*self.chx.shape[loc] as Coordi {
            if loc == 0 { return None; }
            if coordinates[loc]%2 == 1 { self.sum -= 1; }
            coordinates[loc] = 0;
            loc -= 1;
            if coordinates[loc]%2 == 1 { self.sum -= 1; }
            else { self.sum += 1; }
            coordinates[loc] += 1;
            if self.sum > self.dim {
                coordinates[loc] += 1;
            }
        }

        loc = end;
        while self.dim > self.sum {
            coordinates[loc] += 1;
            self.sum += 1;
            if loc > 0 { loc -= 1; }
        }

        if let Some(filvalue) = self.chx.filvalue(coordinates){
            let cube = Cube{ coordinates:coordinates.clone(), filvalue };
            coordinates[end] += 2;
            return Some(cube);
        } else {
            return None;
        }

    }
}


// Implimentation of ChainComplex trait. See chx.rs for definition of ChainComplex trait.
impl<SnzVal, FilVal> ChainComplex<Cube<FilVal>, SnzVal, FilVal> for CubicalComplex<SnzVal, FilVal> where
SnzVal: Clone + Neg<Output = SnzVal> + Debug,
FilVal: PartialOrd + Copy + Debug + Eq + Hash + Ord
{
    type Matrix = CubicalBoundaryMatrix<FilVal, SnzVal>;

    fn get_smoracle(
        &self,
        major_dim       :   MajorDimension, // row or column
        transform       :   ChxTransformKind
    ) -> CubicalBoundaryMatrix<FilVal, SnzVal> {
        CubicalBoundaryMatrix {
            ringmetadata: self.ringmetadata.clone(),
            majordimension: MajorDimension::Row,
            shape: self.shape.clone(),
            entry_array: self.entry_array.clone(),
            max_value: self.max_value
        }
    }

    /// Return an iterator of all cubes for a given dimension
    ///
    /// # Parameters
    /// -`dim`: The given dimension
    fn keys_unordered_itr(&self, dim:usize) -> Box<dyn Iterator<Item=Cube<FilVal>> + '_> {
        let mut coordinates = vec![0; self.shape.len()];
        if dim <= self.shape.len() {
            for ii in self.shape.len()-dim..self.shape.len(){
                coordinates[ii] = 1;
            }
        }

        Box::new( CubeIter {
            chx: &self,
            coordinates,
            dim,
            sum: dim,
            first_cube: true
        })
    }

    fn keys_ordered(&self, dim: usize) -> Vec<Cube<FilVal>> {
        let mut vector: Vec<Cube<FilVal>> = Vec::new();

        for cube in self.keys_unordered_itr(dim) {
            vector.push(cube.clone());
        }

        vector.sort();
        return vector;
    }

    fn key_2_filtration(&self, key: &Cube<FilVal>) -> FilVal {
        return key.filvalue;
    }

    fn max_filtration(&self) -> FilVal {
        return self.max_value;
    }
}
