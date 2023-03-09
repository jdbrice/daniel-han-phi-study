# Intro to ROOT and study of phi meson 

### get the code

```sh
git clone git@github.com:jdbrice/daniel-han-phi-study.git
cd daniel-han-phi-study
```

### setup the ROOT data I/O library
```sh
root -l -b -q  FemtoPairFormat.h+
```

this will give an error but it isnt a problem
```sh
Processing FemtoPairFormat.h+...
input_line_12:2:3: error: use of undeclared identifier 'FemtoPairFormat'
 (FemtoPairFormat())
```

it should produce a few files, notably `FemtoPairFormat_h.so`

### run the analysis code
```sh
root -l -b -q  ana.C
```

this will run `ana.C` script and produce a plot `plot0.pdf`

### Intro Tasks

- Learn how to open a root file in a TBrowser to view its contents
- Look at the `FemtoPairFormat.h` file to see all of the variables in the input data
- Make plots (histograms) of the data for each variable
- Learn how to reconstruct a parent parcile's 4-vector from daugher particle 4-vectors

