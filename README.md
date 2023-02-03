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
Source codes are already compiled for either Windows or Linux users and provided with objct files. So, you just need to create an executable file with your gfortran compiler (linker). To do so, follow the instruction below: 
+ Navigate either to `./source_v77/release_Linux` (for a Linux user) or `./source_v77/release_Windows` (for a Windows user).
+ `makefile` is provided in this folder to link all the pre-compiled object files.
+ So, now type: `make all` in the terminal, then the linker will link all object files and produce an executable file:
  + `run_release_Linux.exe` - for Linux
  + `run_release_Windows.exe` - for Windows
+ Copy & paste the executable file to one of the following example folders:
  + Analytical_test/10_salt_lock_exchange
  + Project/MB_test/ - Note: this is the Mobile Bay test case.
    + 1_barotropic_test
    + 2_baroclinic_test
+ Now, type `./run_release_***.exe` in the terminal, then GOM will run and show simulation progress.
  + Note: the executable file should be located along with `input` and `output` folders, i.e., folder structure should look like:
    + your_project_folder
      + input
      + output
      + run_release_***.exe
  + Note: if you are trying to test `MB_test` cases, you should unzip `./input/hurricane.zip` file first. I put this file as a zip file since there is a file uploading limit in GitHub.    
+ Simulation output files will be located in each example folder's `output` folder.
+ As you can see in the `input` folder, I initially set the output file format with `VTK` file format. So, you will see several `.vtk` files. 
+ If you have a `Tecplot` license, try it (you will have faster simulation results with the Tecplot options than with the VTK options).


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
+ ./source_v77
  + ./release_Linux
    + ./`pre-compiled objects files & makefile are here`
  + ./release_Windows
    + ./`pre-compiled objects files & makefile are here`
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