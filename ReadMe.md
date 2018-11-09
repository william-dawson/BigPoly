# BigPoly

In this program, I provide some utilities that demonstrate the interoperability
of BigDFT and NTPoly.

## How To Compile
Go into the Build folder and type:

> cmake .. -DCMAKE_PREFIX_PATH="/path/to/ntpoly/install;/path/to/bigdft/install"

Note the semicolon between paths.

## How To Run

This program generates a single executable called `BigPoly`. To see a list
of input options, type:

> /path/to/BigPoly --help
