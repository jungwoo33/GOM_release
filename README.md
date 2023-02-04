# README by Jungwoo Lee
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## 1. Title: The General Ocean Model (GOM)

## 2. Description: 
A three-dimensional unstructured grid finite-volume model for coastal and estuarine circulation, which I (Jungwoo Lee) named the General Ocean Model (GOM), has been developed. Combining the finite volume and finite difference methods, GOM achieved both the exact conservation and computational efficiency. The propagation term was implemented by a semi-implicit numerical scheme, so-called θ scheme, and the time-explicit Eulerian-Lagrangian Method (ELM) was used to discretize the nonlinear advection term to remove the major simulation limitations of the time step, which appears when solving shallow water equations, by the Courant-Friedrichs-Lewy stability condition. Because the GOM uses orthogonal unstructured computational grids, allowing both triangular and quadrilateral grids, considerable flexibility to resolve complex coastal boundaries is allowed without any transformation of governing equations. More fundamental details of the GOM can be found in the original development paper (Lee et al., 2020; https://doi.org/10.3390/w12102752) or in the model homepage, https://ufgom.org/publications/.

## 3. Installation:
Please make sure to have following programs installed on your computer to use this app:
+ Git
+ gfortran

To use this program, `git clone` the repo down to your local. Now, you are ready to go!

## 4. Usage:
Source codes are pre-compiled both on Linux and Windows with Gfortran. Note: I am releasing just serially compiled executable now, but the source codes are parallelized with OpenMP.

An executable file for both Windows and Linux are provided in `./release`:
  + `run_release_Linux.exe` - for Linux
  + `run_release_Windows.exe` - for Windows

Also, they are located in each project folders, so just navigate to the one of the following example folders:
  + Analytical_test/
    + 10_salt_lock_exchange
  + Project/MB_test/ - Note: this is the Mobile Bay test case.
    + 1_barotropic_test
    + 2_baroclinic_test

Note, if you are trying to test `MB_test` cases, you should unzip `./input/hurricane.zip` file first. I put this file as a zip file since there is a file uploading limit in GitHub.    

Now, the folder structure will look like:
  + your_project_folder
    + input
    + output
    + run_release_***.exe

Now, type `./run_release_***.exe` in the terminal. Then, you will see the simulation progress.

Simulation output files will be located in each example folder's `output` folder.
  + Note: there is an `etc.txt` file in the `output` folder, and ignore this (I just put this to keep `output` folder in GitHub since GitHub ignore an empty folder)

As you can see in the `input` folder, I initially set the output file format with `VTK` file format. So, you will see several `.vtk` files. 

If you have a `Tecplot` license, try it (you will have faster simulation results with the Tecplot options than with the VTK options).

This project has the following directory structure:
+ ./Analytical_test
  + ./10_salt_lock_exchange
    + ./input
    + ./output
    + ./`your executable must be located here`
+ ./Project
  + ./MB_test
    + ./1_barotropic_test
      + ./input
      + ./output
      + ./`your executable must be located here`
    + ./2_baroclinic_test
      + ./input
      + ./output
      + ./`your executable must be located here`
+ ./release
  + run_release_Linux.exe
  + run_release_Windows.exe
+ ./LICENSE: GPL v3 License 
+ ./README.md: readme file
+ ./GOM_Manual_draft_v5.pdf
  + This is the model user's manual.

## 5. License:
### The GNU GPLv3 License
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## 6. Application results:
I will upload some results with animations later.

## 8. Questions?:
If you have any questions, feel free to contact me via information below:<br>
[Email:] jungwoo33@gmail.com

- - -
© 2023 Jungwoo Lee. Confidential and Proprietary. All Rights Reserved.