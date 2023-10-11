#!/bin/bash
#
cp ../src/*.java .
#cp ../src/pangenie.sh .
docker build --progress=plain -t fcunial/sv-merging .
docker push fcunial/sv-merging
rm -f *.java 
#rm -f pangenie.sh
