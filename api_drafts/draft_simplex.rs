/// The simplex struct and its implementations

/// The struct for simplex
pub struct Simplex<T>{
    pub vertices: Vec<T>
}


// INTUITION FOR THE COFACE ITERATOR:
// start with: [2, 4, 7]
//
// let L = an empty list
//
// v = [*, 2, 4, 7]
// v[0] = 0
// append v.clone() to L
// v[0] = 1
// append v.clone() to L
//
// v[0] = 2 // so v = [2, *, 4, 7]
// v[1] = 3
// append v.clone() to L
//
// v[1] = 4 // so v = [2, 4, *, 7]
// v[2] = 5
// append v.clone() to L
// v[2] = 6
// append v.clone() to L
//
// v[3] = 7 // so v = [2, 4, 7, *]

struct CofaceItr<Coeff, Vertex, Dist>{
    distmat: &,
    face: Vec<Vertex>, // THIS VECTOR SHOULD BE SORTED
    coface_template: Vec<Vertex>, // THE LENGTH OF THIS VECTOR SHOULD EQUAL ITS CAPACITY
    next_vertex_to_check: Vertex, // THIS WILL BE 0 IF THE DISTMAT HAS 0 ROWS/COLUMNS
    num_facevertices_below: usize,
    coefficient: Coeff,
    radius_max: Dist
}

impl<RingKind> Iterator for CofaceItr {
    type Item = (Vec<u16>, RingKind)

    fn next( &self ) -> Item{
        // PULLING OUT VARIABLE NAMES AS SHOWN BELOW MAY OR MAY NOT IMPROVE PERFORMANCE ... I'LL ASK BRYN ABOUT IT
        let dmat = self.distmat
        let numrow_dmat = dmat.len_of(Axis(0))
        let face = self.face
        let numvert_face = face.len()
        let coface_template = self.coface_template

        // these conditions imply that no cofaces exist
        if numrow_dmat == 0 | face[numvert_face-1] == numrow_dmat {
            return None
        }

        // we have a face of form [a_1 < ... < a_n] and a starting vertex a_{m-1} < v < a_m.  We will search over all the vertices in search_range:= [v, ..., a_m-1] to see if we can insert one.
        while true{
            let inserting_top_vrt : bool = num_facevertices_below == numvert_face
            let search_range = match {
                inserting_top_vrt => self.next_vertex_to_check .. numrow_dmat,
                _ => self.next_vertex_to_check .. self.face[self.num_facevertices_below]
            }
            for vertex_to_insert in search_range{
                diameter = maximum( dmat.slice( [ s![ vertex_to_insert, .. ] )[ face ] ) // I DON'T KNOW IF THIS INDEXING WILL ACTUALLY WORK
                if diameter <= self.radius_max{
                    coface_template[self.num_facevertices_below] = vertex_to_insert
                    return coface_template.clone()
                }
            }

            // if we reach this point, it means we have searched across the whole search_range:= [v, ..., a_m-1] without finding a single vertex to insert

            if inserting_top_vrt{
                // if we make it to hear, then we have generated all the cofaces and we are done

                // reset everything to initial condition (this means 'sliding' all the vertices in the template vector over to the righthand side)
                for ind in 0 .. numvert_face{
                    template_vector[numvert_face-ind] = template_vector[numvert_face-ind-1]
                }

                // find the next starting location (ie, the first vertex not in the face).  Note that there will always be at least one vertex outside the face, since an if-statement at the top of this function catches the case where numvrt_fact = numrow_dmat.
                for candidate_vertex in 0..numvert_face+1{
                    if face[candidate_vertex] != candidate_vertex {
                        self.next_vertex_to_check = candidate_vertex
                        self.num_facevertices_below = candidate_vertex
                        self.coefficient = (-1)^candidate_vertex
                        break
                    }
                }

                // return none, since we have generated all the cofaces
                return None
            } else{
                // if we make it hear, then we must define a new search range
                self.next_vertex_to_check =

                // THIS MUST BE FILLED IN
            }

        }


    }
}





impl CofaceItr{

}
