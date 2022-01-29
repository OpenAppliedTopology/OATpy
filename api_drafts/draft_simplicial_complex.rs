
// ================================================================================================
// FILE DESCRIPTION
// This file contains some outlines for the SimplicialComplex trait, and its associated objects, functions, etc.
//
// FILE HISTORY
// 20201220: created by Greg Henselman-Petrusek
// 20201220-present: edited by Greg Henselman-Petrusek and Haibin Hang
// ================================================================================================


pub trait SimplicialComplex<VertexType>{


	// determines the order on the vertices
	fn le( v0: VertexType, v1: VertexType) -> bool

	// 
	fn cup_product(cochain0: SparseVector<VertexType>, cochain1: SparseVector<VertexType>) -> SparseVector<VertexType> 
}