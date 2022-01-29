use crate::GetRow;

/// Sparse matrix in compressed row form
//#[derive(PartialEq,Eq,Debug)]
pub struct CmpRowFmt {
    //pub maxnumnz: usize,               // maximum number of non-zero entries
    pub numrow: usize,                 // number of rows
    pub numcol: usize,                 // number of columns
    pub rowptr: Vec<usize>,            // row end pointer
    pub colind: Vec<usize>,            // colum indices, size maxnumnz
    pub snzval: Vec<f64>,              // structural non-zero entry values
    //reverse_row_order: bool,           // reverse rows or not
}

impl CmpRowFmt {
    /// Construct a sparse matrix in compressed row format with all entries being zero
    ///
    /// # Parameters
    /// - `row_num`: Number of rows
    /// - `col_num`: Number of columns
    ///
    /// # Returns
    /// - A sparse matrix in compressed row format with all entries being zero
	pub fn initialize_csr(row_num:usize, col_num:usize /*, reverse:bool*/)->CmpRowFmt{
    	return CmpRowFmt{
            //maxnumnz: 0,
    		numrow:row_num,
    		numcol:col_num,
    		rowptr: vec![0],
    		colind: Vec::new(),
    		snzval: Vec::new(),
            //reverse_row_order: reverse,
    	};
	}

    /// Print the sparse matrix in compressed row format
    pub	fn print_csr(&self){
        println!("Matrix in compressed row format:");
    	println!("numrow={}, numcol={};\nrowptr={:?}; \ncolind={:?};\nsnzval={:?}.",self.numrow,self.numcol,self.rowptr,self.colind,self.snzval);
	}

    /*
    /// Print the sparse matrix in dense format
    pub	fn print_dense(&self){
        for ii in 0..self.numrow {
            for jj in 0..self.numcol {
                if self.rowptr[]
                let mut kk:usize = self.colind[self.rowptr[ii]];
                if jj == self.colind[kk] {print!("{}, ", self.snzval[kk]); kk+=1; }
                else {print!("0, ");}
            }
            println!(" ");
        }
    }
    */

    pub fn reverse_row_order(self)-> CmpRowFmt {
        let mut reversed = CmpRowFmt::initialize_csr(self.numrow,self.numcol);
        for ii in 0..self.numrow {
            let row_ind = self.numrow-1-ii;
            for jj in self.rowptr[row_ind]..self.rowptr[row_ind+1] {
                reversed.push(self.colind[jj],self.snzval[jj]);
            }
            reversed.rowptr.push(reversed.colind.len());
        }
        reversed
    }


    /// Add a new entry to the sparse matrix with out updating the rowptr
    ///
    /// # Parameters
    /// - `index`: The colum index of the entry
    /// - `value`: The value of the entry
	pub fn push(&mut self, index:usize, value:f64 /*, newrow:bool*/) {
		self.colind.push(index);
		self.snzval.push(value);
        /*
        if newrow==true {
            self.rowptr.push(self.maxnumnz);
        } else {
            let aa = self.rowptr.len()-1;
            self.rowptr[aa] = self.maxnumnz;
        }
        */
	}

}

impl GetRow for CmpRowFmt{
    /// Get a specific row of a spare matrix
    ///
    /// # Parameters
    /// - `row_index`: The row index of the indicated row
    ///
    /// # Returns
    /// - A vector of colum indices of the indicated row
    /// - A vector of entry values of the indicated row
    fn row(&self,row_index:usize) -> (&[usize],&[f64]) {
    	let ind_vec = &self.colind[self.rowptr[row_index]..self.rowptr[row_index+1]];
    	let val_vec = &self.snzval[self.rowptr[row_index]..self.rowptr[row_index+1]];
    	return (ind_vec,val_vec);
    }
}
