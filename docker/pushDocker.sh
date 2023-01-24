#!/bin/bash
#
cp ../src/*.java .
docker build --progress=plain -t fcunial/sv-merging .
docker push fcunial/sv-merging
rm -f *.java
