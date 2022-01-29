## References as attribute values

It seems that an attribute can be a reference, but you have to specify the lifetime of the reference?  See [this example](https://doc.rust-lang.org/1.9.0/book/lifetimes.html).


```
struct Foo<'a> {
    x: &'a i32,
}

impl<'a> Foo<'a> {
    fn x(&self) -> &'a i32 { self.x }
}

fn main() {
    let y = &5; // this is the same as `let _y = 5; let y = &_y;`
    let f = Foo { x: y };

    println!("x is: {}", f.x());
}
```
