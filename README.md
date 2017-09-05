---
layout: default
title: svt
---

## svt function for MATLAB

`svt` is a MATLAB wrapper function for singular value thresholding.

The function is developed by [Cai Li](http://www4.ncsu.edu/~cli9/) and [Hua Zhou](http://hua-zhou.github.io).

### Compatibility

The code is tested on MATLAB R2013a, but should work on other versions of MATLAB with no or little changes. Current version works on these platforms: Windows 32-bit, Windows 64-bit, Linux 64-bit, and Mac (Intel 64-bit). Type `computer` in MATLAB's command window to determine the platform.

### Download

[ZIP File](https://github.com/Hua-Zhou/svt/zipball/master) or [TAR Ball](https://github.com/Hua-Zhou/svt/tarball/master)

### Installation

1. Download the zip or tar package.
2. Extract the zip or tar package.  
```
unzip PackageName.zip
```
or
```
tar xvzf PackageName.tar.gz
```
3. Rename the folder from *PackageName* to *svt*.  
```
mv PackageName svt
```
4. Add the *svt* folder to MATLAB search path. Start MATLAB, cd to the *svt* directory, and execute the following commands  
`addpath(pwd)	%<-- Add the toolbox to the MATLAB path`  
`savepath		%<-- Save for future MATLAB sessions`
5. Go through following tutorials for the usage. For help of individual functions, type `?` followed by the function name in MATLAB.

### Tutorial

* [Singular value thresholding](http://hua-zhou.github.io/svt/html/demo_svt.html)

### Licensing

svt function for MATLAB is licensed under the [BSD](./html/COPYRIGHT.txt) license. Please use at your own risk.

### How to cite

If you use this function in any way, please cite the software itself along with at least one publication or preprint.

* Software reference  
C Li and H Zhou. MATLAB svt function Version 0.0.1. Available online, July 2014.  
C Li and H Zhou (2017) svt: Singular Value Thresholding in MATLAB. [Journal of Statistical Software](in press)

### Contacts

Cai Li <cli9@ncsu.edu> | Hua Zhou <hua_zhou@ncsu.edu>
