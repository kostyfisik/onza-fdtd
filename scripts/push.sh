#!/bin/bash
hg -v push ssh://onza-dev@onzafdtd.org//home/nfs-shared/onza-dev/onza-fdtd
echo Updating dir on onzafdtd.org
ssh onza-dev@onzafdtd.org "cd onza-fdtd && hg update"
