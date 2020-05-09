# Scripts for Making Movies

Scripts are provided for producing MP4 and GIF movies showing propagation.

### `make-fixed-frame-mp4.py`

This Python script takes an MBSolve problem defined in a JSON file and outputs
an MP4 video showing the propagation.

```sh
optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Path of input file.
  -c SPEED_OF_LIGHT, --speed-of-light SPEED_OF_LIGHT
                        Speed of Light in the system units. (Default: 0.1)
  -m Y_MIN, --y-min Y_MIN
                        Minimum of the y-axis. (Default: 0.0)
  -y Y_MAX, --y-max Y_MAX
                        Maximum of the y-axis maximum. (Default: 1.0)
  -z ZOOM, --zoom ZOOM  To use interpolation on the output data, select the
                        order of interpolation. (e.g. 2, 4). Note this may
                        introduce numerical artefacts. (Default: 1)
```

Only the path of the input file is required.

### `make-ffmpeg-gif.sh`

This bash script takes an MP4 file output from
make-fixed-frame-mp4.py and converts it to an animated gif file.

```sh
optional arguments:
  -f FILENAME, --file FILENAME
                        Path of input MP4 file.
  -n INSCALE, --in-scale INSCALE
                        The width in pixels of the input. (Default: 900)
  -s SCALE, --scale
                        The width in pixels of the output gif. (Default: 900)
  -i INFPS, --in-fps
                        The frames-per-second of the input MP4. (Default: 30)
  -p FPS, --fps FPS
                        The frames-per-second of the output gif. (Default: 30)
```

Only the path of the MP4 file is required. It requires [FFmpeg][ff], which can
be installed in a conda env with
```
conda install -c conda-forge ffmpeg
```

[ff]: https://www.ffmpeg.org/
