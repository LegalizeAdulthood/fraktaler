#!/bin/bash
for prefix in ${HOME}/win/posix/*
do
  cat OpenEXR.pc.in |
  sed "s|PREFIX|$prefix|g" > "$prefix/lib/pkgconfig/OpenEXR.pc"
done
