#!/bin/bash
FILES=`find ../src -type f -name \*.cc -o -name \*.h`
for file in $FILES
do
 echo >> result-cpplint.log    
 ./cpplint.py $file 2>>result-cpplint.log
done
cat result-cpplint.log |grep ^../src
echo
cat result-cpplint.log |grep errors
rm result-cpplint.log
