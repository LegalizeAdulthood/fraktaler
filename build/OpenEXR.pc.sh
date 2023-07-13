#!/bin/bash
for prefix in ${HOME}/opt/windows/posix/*
do
  sed "s|PREFIX|$prefix|g" < OpenEXR.pc.in > "$prefix/lib/pkgconfig/OpenEXR.pc"
done
