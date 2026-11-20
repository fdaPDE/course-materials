#!/bin/bash
set -e
export DEBIAN_FRONTEND=noninteractive
apt-get -y update
apt-get -y install libcurl4-openssl-dev libssl-dev libjpeg-dev libpng-dev libicu-dev cmake make libuv1-dev pandoc libabsl-dev libgdal-dev gdal-bin libgeos-dev libproj-dev libsqlite3-dev libudunits2-dev zlib1g-dev libfontconfig1-dev libfreetype6-dev libfribidi-dev libharfbuzz-dev libx11-dev git libgit2-dev libtiff-dev libwebp-dev libxml2-dev
