* ndarray crate
	* [a nice intro for people used to numpy](https://docs.rs/ndarray/0.12.1/ndarray/doc/ndarray_for_numpy_users/index.html)

* the [use](https://doc.rust-lang.org/rust-by-example/mod/use.html) method for extracting names from deeply nested lists	

* lifetimes
	* [link](https://doc.rust-lang.org/1.9.0/book/lifetimes.html)
	* intimately related to passing references as outputs
		* some examples:
		
			```
			fn print(s: &str); // elided
			fn print<'a>(s: &'a str); // expanded
				
			fn debug(lvl: u32, s: &str); // elided
			fn debug<'a>(lvl: u32, s: &'a str); // expanded
				
			// In the preceding example, `lvl` doesn’t need a lifetime because it’s not a
			// reference (`&`). Only things relating to references (such as a `struct`
			// which contains a reference) need lifetimes.
				
			fn substr(s: &str, until: u32) -> &str; // elided
			fn substr<'a>(s: &'a str, until: u32) -> &'a str; // expanded
				
			fn get_str() -> &str; // ILLEGAL, no inputs
				
			fn frob(s: &str, t: &str) -> &str; // ILLEGAL, two inputs
			fn frob<'a, 'b>(s: &'a str, t: &'b str) -> &str; // Expanded: Output lifetime is ambiguous
				
			fn get_mut(&mut self) -> &mut T; // elided
			fn get_mut<'a>(&'a mut self) -> &'a mut T; // expanded
				
			fn args<T: ToCStr>(&mut self, args: &[T]) -> &mut Command; // elided
			fn args<'a, 'b, T: ToCStr>(&'a mut self, args: &'b [T]) -> &'a mut Command; // expanded
				
			fn new(buf: &mut [u8]) -> BufWriter; // elided
			fn new<'a>(buf: &'a mut [u8]) -> BufWriter<'a>; // expanded
			```