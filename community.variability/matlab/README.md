# MATLAB Community Variability

This folder contains dependency-free MATLAB functions for metacommunity
variability metrics using:

```text
X(time, site, taxon)
```

## Install From A Checkout

```matlab
run("community.variability/matlab/install_community_variability.m")
install_community_variability(true)
```

`true` saves the MATLAB path for future sessions. Omit the argument to only add
the path for the current session.

## Fetch Only This Folder From Git

MATLAB does not install subdirectories directly from Git the way Python and R
do. Use sparse checkout when you only need the MATLAB code:

```bash
git clone --filter=blob:none --sparse https://github.com/anhph95/ifcbtools.git
cd ifcbtools
git sparse-checkout set community.variability/matlab
```

Then run the installer above.
