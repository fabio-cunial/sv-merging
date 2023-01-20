#!/bin/bash
#
docker build --progress=plain -t fcunial/sv-merging .
docker push fcunial/sv-merging
