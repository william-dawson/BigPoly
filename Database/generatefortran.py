"""This script takes a yaml file, and turns it into a string contained inside
a fortran module.

Usage:
  python generatefortran.py input_file.yaml output_file.f90

  The name of the module will be the same as the ``output_file``.
"""
from __future__ import print_function
from sys import argv
from yaml import load, dump
from os.path import basename

if __name__ == "__main__":
    # Read input
    try:
        infile = argv[1]
        outfile = argv[2]
    except:
        print(__doc__)
        exit(-1)

    # Read in database file
    with open(infile) as ifile:
        database = load(ifile)
    data_str = dump(database, default_flow_style=False)
    data_str.replace('"','\\"')

    # Setup the module
    modname = str(basename(outfile).replace(".f90","Module"))
    ostring = ""
    ostring += "MODULE "+modname+"\n"
    ostring += "\tIMPLICIT NONE\n"

    # Add the string
    concat = "//"
    cr = "new_line('A')"
    ostring += "\tCHARACTER(LEN=*), PARAMETER, PUBLIC :: datastream =&\n"
    for i, line in enumerate(data_str.splitlines()):
        ostring+='"'+line+'"'+concat+cr
        if i < len(data_str.splitlines()) - 1:
            ostring += concat+"&\n"
    ostring +='\n'

    # Cleanup
    ostring += "END MODULE\n"

    # Write to file
    with open(outfile, "w") as ofile:
        ofile.write(ostring)
