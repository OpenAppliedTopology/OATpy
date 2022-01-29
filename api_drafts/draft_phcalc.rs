
// ================================================================================================
// FILE DESCRIPTION
// This file is intended as a blueprint for running the persistence calculation using the FilteredChainComplex api.
//
// FILE HISTORY
// 20201202: created by Greg Henselman-Petrusek
// 20201202-present: edited by Greg Henselman-Petrusek and Haibin Hang
// ================================================================================================


// ------------------------------------------------------------------------------------------------
// OVERVIEW OF CONTENTS
// ------------------------------------------------------------------------------------------------

// This is just a rough overview of the persistence computation.


// ------------------------------------------------------------------------------------------------
// PERSISTENCE COMPUTATION
// ------------------------------------------------------------------------------------------------

fn workflow.execute() -> OutputObject{

    // Initialize the output object
    let OutObj : OutputObject = Deafault::default();

    // Load `OutObj` with all the user input data from the workflow obj
    // ..

    // Construct the filtered chain complex
    let OutObj.Chx : ChainComplex{ .. };

    // Get some pointers to variables we will use often.
    let Chx = OutObj.Chx;
    let PHspec = OutObj.PHspec; // PHspec is an attribute of OutObj containing the specifications/instructions for the PH computation
    let PH_chaindegrees = PHspec.chaindegrees; // the homology degrees in which to compute PH

    // Compute PH
    //------------------------------------------------------------------

    for chaindeg_PHloop in PH_chaindegrees{

        // Construct the boundary matrix
        let bdmat = Chx.matrix(chaindeg = chaindeg_PHloop, ... )

        // Determine which rows to reduce.
        let rows_to_reduce == [an array of row index keys; keys could be simplices, integers, or something else]

        // Initialize the L and R matrices
        let L = initialize L;
        let R = initizlize R;

        // Compute the LR factorization
        for row_to_reduce in rows_to_reduce{
            // do stuff
        }

        // Store the L and R matrices in the Chx object
        // ..

    }
}
