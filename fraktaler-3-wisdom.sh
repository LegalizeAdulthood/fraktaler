#!/bin/bash
cli="$(pwd)/fraktaler-3-cli.gcc"
cl="$(pwd)/fraktaler-3-cl.gcc"
tmp="$(mktemp -d)"
echo > "${tmp}/wisdom.log" -e "# platform\tdevice\tnumbertype\tmantissa\texponent\terror\tseconds"
{
  printf "%s: without OpenCL\n-1.0: system CPU\n" "-1"
  clinfo --list --raw
} |
while read platformdevice name
do
  pattern='^([-0-9]+)\.([-0-9]+):$'
  if [[ "${platformdevice}" =~ ${pattern} ]]
  then
    platform="${BASH_REMATCH[1]}"
    device="${BASH_REMATCH[2]}"
    devicename="${name}"
    if (( platform < 0 ))
    then
      numbertypes=("float:24:8" "double:53:11" "long double:64:15" "floatexp:24:24" "softfloat:32:31" "float128:113:15")
      fraktaler="${cli}"
    else
      numbertypes=("float:24:8" "double:53:11" "floatexp:24:24" "softfloat:32:31")
      fraktaler="${cl}"
    fi
    for numbertypespec in "${numbertypes[@]}"
    do
      echo "${numbertypespec}" |
      {
        IFS=":" read numbertype mantissa exponent
        cat > "${tmp}/wisdom.f3.toml" <<EOF
program = "fraktaler-3"
version = "0-338-g27a93b9"
location.real = "0.352465816656845732"
location.imag = "0.392188990843255425"
location.zoom = "4.1943021e6"
image.subframes = 1
render.filename = "${tmp}/wisdom-${platform}-${device}-${numbertype}"
algorithm.number_types = ["${numbertype}"]
opencl.platform = ${platform}
opencl.device = ${device}
opencl.tile_width = 512
opencl.tile_height = 288

[[formula]]
power = 2
neg_y = false
neg_x = false
abs_y = false
abs_x = false
[[formula]]
power = 2
neg_y = false
neg_x = false
abs_y = true
abs_x = false
[[formula]]
power = 2
neg_y = false
neg_x = false
abs_y = false
abs_x = false
[[formula]]
power = 2
neg_y = false
neg_x = false
abs_y = true
abs_x = true


EOF
        printf "%s\t%s\t%s\t%s\t%s\n" "${platformname}" "${devicename}" "${numbertype}" "${mantissa}" "${exponent}"
        /usr/bin/time -a -o "${tmp}/wisdom.log" -f "${platform}\t${device}\t${numbertype}\t${mantissa}\t${exponent}\t%x\t%e" "${fraktaler}" "${tmp}/wisdom.f3.toml" >/dev/null 2>/dev/null
      }
    done
  else
    platformname="${name}"
  fi
done
cp -avi "${tmp}/wisdom.log" fraktaler-3.wisdom
