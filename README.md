# BHInfTurbHTens
Under construction (!): Need to split the generation of Fig10-11-12-13 which currently are in the same file and Github says is to large. The same applies to Fig16.ipynb that needs to be thiner.

Code and data from the article "Black hole information turbulence and the Hubble tension" Copyright (C) 2025 Juan Luis Cabrera Fernández

This project is licensed under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International - Creative Commons licence (CC BY-NC-SA 4.0). https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en

Composed by jupyter notebook .ipynb files written in Julia. Some notebooks load Julia .jl routines. Certain data outputs are already available and loaded as .txt or .jld files to expedite the process, but they are also generated somewhere else.

To recreate the full set of results as in the paper you may like to run the code  in the following sequence:

1- Fig1-5.ipynb (loads darkcolored2.jl and cascadas_fractal3.jl. Outputs fig 1 and fig 5).

2- Fig2.ipynb (loads darkcolored2.jl and outputs fig 2)

2- Fig3-8-9.ipynb (Loads figXsoloY1.jl with X: HLCDM, EDE or EDEP loads figHLCDMsoloY1.jl , figEDEsoloY1.jl and figEDEPsoloY1.jl. Outputs fig3, fig8 and fig9)

3- Fig6-7.ipynb (loads "dims_ReglaCascada_iter27.txt" with the cascade calculated with eqs. 4-6 using dim_cascade2.jl. Outputs figures 6 and 7)

4- Fig10-11-12-13.ipynb (Paper figures 10,11,12 and 13. Generates .jld files with the posterior data used to plot figures 2 and 16)

5- Fig14.ipynb (Outputs figure 14 and a .jld file with the posterior data used to plot figures 2 and 16)

6- Fig15.ipynb  (Outputs figure 15 and a .jld file with the posterior data used to plot figures 2 and 16)


Code for fig 4 is not provided as it is just a LaTeX rendering. The above .ipynb are also translated as .pdf files for fast review. In some ocassions these versions may be a bit larger as they could include some comments deleted in the uploaded .ipynb files.  



