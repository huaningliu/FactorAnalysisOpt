# FactorAnalysisOpt

## To Run the Replicated Simulation
Mean-based method:
- generate data: select all code in `Data-gen-andoBai.R` and run
- benchmark: select all code in `Ando-Bai-Original.R` and run
- our code: select all code in `Ando-Bai-Opt.R` and run
- To compare: Look at the convergence process where the error is printed when running the codes

Quantile-based method:
- benchmark: run `Mul_interactNew.m` on MATLAB softare
- our code: select all code in `mean-based-data-gen.R` and run
- To compare: check the value of `BETA` on both versions and compare that to true beta: (1, 3, 5, 2, 4)

## To Run invented ES factor model and its related Simulation
Two-step framework:
- Select all code in `twoStepFAnew.R` and run, losses will be printed; then select all code in `mean-based-real.R` and run, losses will be printed
- Note: You will need to change `test_period`, `r_tau`, `tau` in order to give each subgroup of data a proper factor parameter.

To compare the 3 factor loadings generated from the two frameworks:
- Select all code in `canonicalAna.R` and run, two tables with p-values will be provided.

Note: data paths need to be changed when running

## Reference
Mean-based Estimation Method: [Bai & Ng 2008](http://www.columbia.edu/~sn2294/pub/eco-002.pdf), [Bai 2005 Estimation 3.2](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.115.5857&rep=rep1&type=pdf)

Mean-based Estimation Benchmark Code (in MATLAB): [Bai 2005 Suppliment code (out-resource)](https://ideas.repec.org/c/boc/bocode/m430011.html)

Quantile-based Estimation Method: [Ando & Bai 2018](https://par.nsf.gov/servlets/purl/10163503)

Qunatile-based Estimation Benchmark Code (in R): [Ando & Bai 2018 Suppliment code](https://www.tandfonline.com/doi/suppl/10.1080/01621459.2018.1543598)

Real world data source: [Fama/French 5 factor regional data](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html)
