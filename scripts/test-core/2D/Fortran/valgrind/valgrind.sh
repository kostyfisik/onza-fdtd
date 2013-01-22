 #!/bin/bash
FILE=$1
valgrind --tool=cachegrind --L2=8388608,8,64 ./$FILE
