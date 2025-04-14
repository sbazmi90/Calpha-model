# C<sub>α<sub> model: Coarse-Grained Protein Folding and Fold Switching Simulations

This C code implements a one bead-per-amino acid coarse-grained
model for protein folding and protein fold switching, with a structure-based potential.
Conformational sampling is carried out using Langevin dynamics.

This repository developed in the [Wallin Lab](https://www.physics.mun.ca/~jswallin/index.html) under the supervision of **Dr. Stefan Wallin**, this model allows researchers to study **protein folding, fold switching, and particularly macromolecular crowding effects** on proteins.

## How to cite

If you use the Calpha model or build upon it in your own work, please cite the following publications:

Simulations of a protein fold switch reveal crowding-induced population shifts driven by disordered regions
S Bazmi, B Seifi, S Wallin
Communications Chemistry 6, 191 (2023)
[Read on Nature](https://www.nature.com/articles/s42004-023-00995-2)

This study uses the Calpha model to explore how macromolecular crowding impacts protein fold switching, revealing the role of intrinsically disordered regions in driving population shifts between conformational states.


Conformational entropic barriers in topology-dependent protein folding: perspectives from a simple native-centric polymer model
S Wallin, H.S. Chan
Journal of Physics: Condensed Matter 18, S307–S328 (2006)
[Link to paper](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://iopscience.iop.org/article/10.1088/0953-8984/21/32/329801/pdf)

This foundational paper outlines the theoretical framework behind structure-based coarse-grained protein models and insights into entropic barriers and folding pathways using native-centric approaches.


A minimal sequence of commands to compile the code and run a
simulation of a single chain at a fixed temperature is:

make constants

./constants input 1

make fixtemp

./main

Contact: Stefan Wallin, swallin@mun.ca

