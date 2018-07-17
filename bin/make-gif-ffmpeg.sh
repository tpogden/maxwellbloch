#!/bin/bash

# Takes an MP4 file output from make-fixed-frame-mp4.py and converts it to an
# animated gif file.

FILENAME=""
INSCALE=900 # The width in pixels of the input MP4.
SCALE=900 # The width in pixels of the output gif.
INFPS=30 # The frames-per-second of the input MP4.
FPS=30 # The frames-per-second of the output gif.

# https://stackoverflow.com/a/14203146
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -f|--file)
    FILENAME="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--in-scale)
    INSCALE="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--scale)
    SCALE="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--in-fps)
    INFPS="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--fps)
    FPS="$2"
    shift # past argument
    shift # past value
    ;;
    *) # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo "Converting ${FILENAME} to gif format"
echo "Input Scale: $INSCALE"
echo "Input FPS: $INFPS"

echo "Output Scale: $SCALE"
echo "Output FPS: $FPS"

Y = ${FILENAME%.mp4}
GIF="$Y.gif" # Assumes file extension is .mp4

ffmpeg -y -i $FILENAME -vf \
    fps=$INFPS,scale=$INSCALE:-1:flags=lanczos,palettegen palette.png

echo "Outputting to $GIF"

ffmpeg -y -i $FILENAME -i palette.png -filter_complex \
    "fps=$FPS,scale=$SCALE:-1:flags=lanczos[x];[x][1:v]paletteuse" $GIF
