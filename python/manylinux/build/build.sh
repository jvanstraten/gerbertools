#!/bin/bash
set -e

rm -rf ../dist
mkdir -p ../dist
docker build --pull -t gerbertools-manylinux ../../.. -f Dockerfile
docker run -v `pwd`/../dist:/io/dist gerbertools-manylinux python/manylinux/build/helper.sh
