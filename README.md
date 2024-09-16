# Altera FFT Tutorial

8/31/2024 D. W. Hawkins (dwh@caltech.edu)

## Introduction

This repository contains a tutorial on the Intel/Altera FFT.

Directory           | Contents
--------------------|-----------
doc                 | Tutorial documentation
designs             | Quartus designs
matlab              | MATLAB analysis
references          | Reference documentation

## Git LFS Installation

This repository was created using the github web interface, and then checked out using Windows 10 WSL. Update WSL and install Git LFS support using

~~~
$ sudo apt update && sudo apt upgrade -y
$ sudo apt install git-lfs
~~~

Check the repo and install git LFS installed using

~~~
$ git clone git@github.com:d-hawkins/altera_fft_tutorial.git
$ cd altera_fft_tutorial/
$ git lfs install
~~~

The .gitattributes file from another repo was then copied to this repo, and that file checked in.

~~~
$ git add .gitattributes
$ git commit -m "Git LFS" .gitattributes
$ git push
~~~

The .gitattributes file contains file extension patterns for the majority of binary file types that could be checked into the repo (additional patterns can be added as needed).

