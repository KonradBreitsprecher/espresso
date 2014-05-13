from __future__ import print_function
import featuredefs, sys

if len(sys.argv) != 3:
    print("Usage: %s DEFFILE CPPFILE" % sys.argv[0], file=sys.stderr)
    exit(2)

deffilename, cfilename = sys.argv[1:3]
    
print("Reading definitions from " + deffilename + "...")
defs = featuredefs.defs(deffilename)
print("Done.")

# generate cpp-file
print("Writing " + cfilename + "...")
cfile = open(cfilename, 'w');

cfile.write("""
#include "config.hpp"
#include <iostream>
using namespace std;
int main() {

cout << "# This file was autogenerated." << endl
     << "# Do not modify it or your changes will be overwritten!" << endl;

""")

template="""
#ifdef {0}
cout << "DEF {0} = 1" << endl;
#else
cout << "DEF {0} = 0" << endl;
#endif
"""

for feature in defs.allfeatures:
    cfile.write(template.format(feature))

cfile.write("""
}
""")
    
cfile.close()
print("Done.")

