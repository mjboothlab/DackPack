# DackPack - Nucleic Acid Analysis

DackPack is a python-based graphical user interface for the NUPACK<sup>1,2</sup> package, designed to aid nucleic acid sequence design and analysis. As demonstarted below, it allows users to input nucleic acid sequences and parameters, run thermodynamic calculations, and visualise the resulting structures via NucDraw<sup>3</sup>.

Developed by Chelsea Dack (chelsea.dack.24@ucl.ac.uk), PhD student in the [Booth group](https://boothlab.uk/).

## In Developement
- Integration with Multistrand<sup>4,5</sup> for kinetic simulations
- Update to latest NUPACK<sup>1,2</sup> version to enable DNA/RNA hydrids, OMe RNA, and wider salt conditions

## Installation
1. First, make sure the following requirements are installed on your host system:
   - Python 3.10+
   - [NUPACK 4.0.2+](https://www.nupack.org/download/overview) (see NUPACK installation instructions [here](https://docs.nupack.org/start/#installation-requirements))
2. If not already installed, also install the following packages:
    - `pip install pandas`
    - `pip install numpy`
    - `pip install matplotlib`
    - `pip install PySide6`
    - `pip install nucdraw`

## Use
To launch DackPack on MacOS:
   - Open a new Terminal window and type this line, substituting the location of the DackPack.py script:
     - `python3 ~/insert_path_to/DackPack.py`
   - Note: A standalone DackPack application (not included here) has also been developed for MacOS. This allows users to launch the GUI without using Terminal. The NUPACK module still requires separate installation however, due to licence restrictions.

To launch DackPack on Windows (using Ubuntu, assuming NUPACK installation as described [here](https://docs.nupack.org/start/#installation-requirements)):
   - Move the DackPack.py script into Ubuntu's Linux home directory
   - Open a new Ubuntu window and type this line:
     - `python3 DackPack.py`

To analyse nucleic acid sequences:
1. Upon launching the script, this window will appear, where you can input your parameters and nucleic acid sequence(s):
<img src="example_imgs/user_interface.png" width="250" height="350">
2. Click 'Run Analysis' for calculation of equilibrium concentrations (via NUPACK<sup>1,2</sup>). This window will appear:
<img src="example_imgs/example_equilibrium.png" width="600" height="150">
3. Click 'Plot Structure' for structure visualisation of a chosen complex (using NucDraw<sup>3</sup>). This window will appear, where individual bases / strand backbones are coloured by identity and base pairs are coloured by MFE probability:
<img src="example_imgs/example_structure.png" width="400" height="300">
Figure appearance can be customised, for example by changing base size, structure rotation, and figure zoom.

## References:
1. M.E. Fornace, N.J. Porubsky, and N.A. Pierce (2020). A unified dynamic programming framework for the analysis of interacting nucleic acid strands: enhanced models, scalability, and speed. ACS Synth Biol, 9:2665-2678, 2020.
2. R. M. Dirks, J. S. Bois, J. M. Schaeffer, E. Winfree, and N. A. Pierce. Thermodynamic analysis of interacting nucleic acid strands. SIAM Rev, 49:65-88, 2007.
3. 10.5281/zenodo.15352138
4. Joseph M. Schaeffer (2013) Stochastic simulation of the kinetics of multiple interacting nucleic acid strands. Dissertation (Ph.D.), California Institute of Technology. [https://thesis.library.caltech.edu/7460/](https://thesis.library.caltech.edu/7460/)
5. Joseph M. Schaeffer, Chris Thachuk, Erik Winfree (2015) Stochastic simulation of the kinetics of multiple interacting nucleic acid strands. DNA Computing and Molecular Programming (DNA21), Lecture Notes in Computer Science (LNCS) volume 9211, pp 194-211. [https://doi.org/10.1007/978-3-319-21999-8_13](https://doi.org/10.1007/978-3-319-21999-8_13)
