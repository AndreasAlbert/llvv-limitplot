# llvv-limitplot

For each paramater point, the plotting scripts read the appropriate root file produced by the HiggsCombine tool.
By default, the scripts expect as an argument a top-level directory containing a separate folder for each parameter 
point of the model.
E.g.
```
  ADD/1
  ADD/2
  ADD/3
```
and so on.
Multiple parameters are divided by underscores, e.g.
```
  MV/10_50
  MV/10_100
```  
and so on for vector-mediated DM production.

For details, best check the code, it is quite easy to understand.

To run everything, I usually do:

```
carddir=/somepath/
root -l -b -q plot2D_DMV.cc"(\"$carddir/MA\")"
root -l -b -q plot2D_DMV.cc"(\"$carddir/MV\")"
root -l -b -q plotDM_EWK_1D.cc"(\"$carddir/MEWK\")"
root -l -b -q plotDM_EWK_K1K2.cc"(\"$carddir/MEWK\")"
root -l -b -q plotDM_EWK_K1K2_mu.cc"(\"$carddir/MEWK\")"
root -l -b -q plotUnpart.cc"(\"$carddir/Unpart\")"
root -l -b -q plotADD.cc"(\"$carddir/ADD\")"
root -l -b -q plotS_MChi_1_mu.cc"(\"$carddir/NLOS\")"
root -l -b -q plotP_MChi_1_mu.cc"(\"$carddir/NLOP\")"


root -l -b -q plot_wimpxs_sd.cc;
root -l -b -q plot_wimpxs_si.cc;
````
