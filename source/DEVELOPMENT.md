# Development of TranD

To develop TranD:
* code fork/clone the repository
```bash
$ git clone https://github.com/McIntyre-Lab/TranDi_EA
```
* create a conda environment (conda install mamba first)
```bash
$ mamba create -yp ./conda
$ conda activate ./conda
```
* Install the dependencies and the trand package (in dev mode, so all changes in the code are
  reflected as you rerun the 'trand' script)
```bash
$ mamba install python bedtools
$ cd TranDi_EA/source
$ pip install -e .
```
* Run the 'trand' script.
