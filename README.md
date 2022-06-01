# dotplot
Produce .svg dotplot from .paf alignement


## Install

```bash
cmake .
make
```

## Usage

```
$ ./dotplotter --help
Allowed options:
  --help                 Prints help message
  --paf arg              Input file (.paf format)
  --q arg (=100000)      Minimum query length. Default: 100000
  --t arg (=100000)      Minimum target length. Default: 100000
  --a arg (=10000)       Minimum alignment length. Default: 10000
  --m arg (=0)           Minimum mapping quality. Default: 0
  --s arg (=1)           Stroke width. Default: 1.0
  --i arg (=0.300000012) Minimum alignment identity. Default: 0.3
  --x arg (=4096)        x-size. Default: 4096px
  --y arg (=4096)        y-size. Default: 4096px
  --off arg (=80)        offset of origin. Default: 80px
  --fs arg (=12)         font size. Default: 12pt
```
