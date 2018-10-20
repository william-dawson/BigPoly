# Fake CheSS Conversion

In this program, I demonstrate how one might convert back and forth between
the CheSS segmented CSR matrix format to NTPoly.

(THIS PROJECT STILL ISN'T DONE!)

## How To Compile
Go into the Build folder and type:

> cmake .. -DCMAKE_PREFIX_PATH="/path/to/ntpoly/install;/path/to/bigdft/install"

Note the semicolon between paths.

## Fake CheSS Approach

In this first version, we will use a fake version of the CheSS library. Like the
real deal, the `FakeChesModule` implements a segmented csr matrix. After this
program's first version is completed, we can hopefully drop the real CheSS in and
use the same logic.

## Data Distribution Differences

I make the following assumptions about how CheSS data structures are distributed
in BigDFT. First, the Hamiltonian is copied across all processes. Second, the
density matrix is distributed evenly by columns.

It is not important to know how NTPoly distributes its data. Instead, each process
can ask an NTPoly matrix for its starting rows and columns. We can then extract
the necessary elements from the CheSS data structure, put them in a triplet list,
and call `FillMatrixFromTripletList`.

On the way back, we get the local data elements

## Redistribution Strategy
