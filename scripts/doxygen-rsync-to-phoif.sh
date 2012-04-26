#!/bin/bash
cd ../doc/doxygen
./go
rsync -ave ssh html/ phoif:~/onzafdtd.org/