## Printing

From [stack exchange](https://stackoverflow.com/questions/21747136/how-do-i-print-the-type-of-a-variable-in-rust)

```
#![feature(core_intrinsics)]

fn print_type_of<T>(_: &T) {
    println!("{}", unsafe { std::intrinsics::type_name::<T>() });
}

fn main() {
    print_type_of(&32.90);          // prints "f64"
    print_type_of(&vec![1, 2, 4]);  // prints "std::vec::Vec<i32>"
    print_type_of(&"foo");          // prints "&str"
}
```

## Type parameters

### Use

Rust lets you use type parameters in three places:

* traits: [see "generic traits" here](https://doc.rust-lang.org/reference/items/traits.html)

	```
	trait Seq<T> {
	    fn len(&self) -> u32;
	    fn elt_at(&self, n: u32) -> T;
	    fn iter<F>(&self, f: F) where F: Fn(T);
	}
	```
	(Note that this is different from associated types)
* structs: [see "In Struct Definitions"](https://doc.rust-lang.org/book/ch10-01-syntax.html)

	```
	struct Point<T> {
	    x: T,
	    y: T,
	}
	```
* methods and functions:

	```
	fn largest<T>(list: &[T]) -> &T {
    let mut largest = &list[0];

    for item in list {
        if item > largest {
            largest = item;
        }
    }

    largest
	}
	```


### Performance

There seems to be no loss of performance for using type parameters, in general, since the compiler makes separate versions of code for every possible type.


## Trait objects

### References

[This](https://brson.github.io/rust-anthology/1/all-about-trait-objects.html) seems to have some very good information.

### Use

A trait object is something expressed as `impl trait`.

### Object safety // associated types

You can only turn a trait into a trait object if the trait is "safe": [link](https://doc.rust-lang.org/reference/items/traits.html)

> Object safe traits can be the base trait of a trait object. A trait is object safe if it has the following qualities (defined in RFC 255):
> 
> * It must not require Self: Sized
> * All associated functions must either have a where Self: Sized bound, or
> * Not have any type parameters (although lifetime parameters are allowed), and
> * Be a method that does not use Self except in the type of the receiver.
> * It must not have any associated constants.
> * All supertraits must also be object safe.
> 
> When there isn't a Self: Sized bound on a method, the type of a method receiver must be one of the following types:
> 
> * &Self
> * &mut Self
> * Box<Self>
> * Rc<Self>
> * Arc<Self>
> * Pin<P> 
> 
> where P is one of the types above


The general trend seems to be that trait objects shouldn't have type parameters. **HOWEVER** you can use trait objects with associated types [link](https://depth-first.com/articles/2020/06/22/returning-rust-iterators/):

```
struct Wrapper {
    value: u8
}

struct Container {
    items: Vec<Wrapper>
}

impl Container {
    fn values(&self) -> impl Iterator<Item=u8> + '_ {
        self.items.iter().map(|wrapper| wrapper.value)
    }
}
```


