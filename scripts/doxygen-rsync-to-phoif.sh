#!/bin/bash
cd ../doc/doxygen
./doxygen  
rsync -ave ssh html/ phoif:~/onzafdtd.org/