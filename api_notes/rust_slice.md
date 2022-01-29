* slices with vectors:

	```
	fn main() {
	    let x = vec![0; 8];         // create vector of 8 zeros
	    let y = &x[ 1..3 ];         // create slice by indexing into the vector
	    let z = y.to_vec();         // create a new vector (not a reference)
		let w = z.iter();
	    for element in w {
	    println!("{}", element);
	    }
	}
	```    
