# LDEO_BPR

## LDEO Ocean Bottom Pressure Recorder: BPR/APG Logger

### Installation

From the CLI, ensure you are in the root project directory before starting.

If using a conda environment it is safer to install the dependency modules first using conda. Note that if using PIP the installation will install any necessary dependency modules automatically.

To install a new conda environment with the necessary dependencies use the following command:

```sh
conda env create -f environment.yaml
conda activate bpr
```

Or to update an existing conda environment with the necessary dependencies use the following command:

```sh
> conda env update -f environment.yaml
> conda activate bpr
```

Ensure you have the Python environment active that you intend the package to run in. This applies whether using Conda, pyenv, or some other environment manager. Ensure this environment is enabled/active whenever using this package.

Install the package with the following command. (Note: use pip install even if you are using a Conda env.):

```sh
> python -m pip install --editable .
```

The script can now be run/called using one of the following methods:

```sh
> apg_read <parameters>
# or
> python -m ldeo_bpr.apg_read <parameters>
# or
> python "<path>\ldeo_bpr\ldeo_bpr\apg_read.py" <parameters>
```
