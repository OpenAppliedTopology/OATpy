


use persistence::matrix::{SmOracle, MajorDimension};
use persistence::csm::CSM;





/// Provides access to the upper triangular matrices (and their inverses) in an U-match
/// decomposiion.
///
/// The commands associated with this struct (e.g. those allowing one to invert or change major
/// dimension of a row/column operation matrix) can be much more efficient than the generic
/// options.
struct UMatch<MajKey, MinKey, SnzVal>{
    smoracle: &SmOracle<MajKey, MinKey, SnzVal>,    // the matrix to factor
    factor_data: CSM<usize, SnzVal>,                // partial change of basis matrix
    pivot_bijections: persistence::chx::Indexing    // indexing information
}


impl UMatch<MajKey, MinKey, SnzVal> {
   

    /// Set the major dimension of the CSM stored in the `factor_data` field.
    ///
    /// Changing the major dimension may transpose the `factor_data` CSM stored internally to the struct.  This can dramatically increase the efficiency of the change of basis oracles.
    /// 
    fn set_major_dim( &self, MajorDimension ) {
        self.seed_data = transpose(self.cob);
    }


    /// Get the major dimension of the `factor_data` attribute. 
    fn get_major_dim( &self ) -> &MajorDimension {
        &self.factor_data.majdim
    }
    
    /// Return an upper-unitriangular row/column-operation matrix.
    ///
    /// The major dimension of the output matrix will be the major dimension of the `factor_data`
    /// attribute (not necessarily the major dimension of the matrix being factored).
    ///
    /// # Parameters
    /// 
    /// * `which_side` -- row versus column operation matrix
    /// * `invert` -- if true, then return the inverse of the operation matrix
    ///
    fn cob( &self, which_side : MajorDimension, invert : bool ) -> 
        UMatchCob < MajKey, Minkey, SnazVal > {
    }


}

/// A matrix oracle for the change of basis matrices in a U-match decomposition.
struct UMatchCob<MajKey, MinKey, SnzVal> {
    umatch: &UMatch<MajKey, MinKey, SnzVal>,
    which_side: MajorDimension
}

impl SmOracle < MajKey, MinKey, SnzVal > for UMatchCob {
    
    /// to fill in
    
}

