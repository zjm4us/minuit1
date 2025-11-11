# minuit1

- expFit.cpp: example using the more general fitting interface in minuit
- expFit.ipynb: equivalent example using lmfit
- SimultaneousExps(lm).ipynb: generation of histograms with correlated signals for simultaneous fit exercise
- rootExample.cpp: just another example of using ROOT classes in a C++ program
- *.root files: various input histograms for fitting exercises

-----

Team member names and computing IDs: Malinda Amarakoon (zjm4us) Taylor Colaizzi (zxk4gs)

-----
Exercise 1 comments:

1A — Comments
- Both fits yield consistent central values, but the NLL fit is statistically more accurate for low-count data.  
- The χ² fit runs faster and is appropriate when counts are large.  
- The NLL fit produces slightly smaller uncertainties, reflecting better modeling of Poisson fluctuations.  
- χ² probability (p-value) was computed to assess fit quality.

---

1B — Comments
- The double Gaussian fit gave a better χ²/ndof and smaller residuals, indicating the data likely contain two overlapping peaks.  
- The Gumbel distribution captured asymmetry but fit quality was worse near the tails.  
- Based on χ² probability, the two-Gaussian model is preferred.  
- Fit results, parameters, and plots were saved in **ex1.pdf**.

Exercise 2 comments:
--
- The simultaneous fit improved parameter precision compared to fitting datasets separately.  
- The shared parameters (mean and σ) were consistent across both histograms.  
- χ²/ndof values indicate good agreement between model and data for both fits.  
- The combined uncertainty on the signal parameters was smaller due to the joint constraint.  
- Plots and fit results were saved in **ex2.pdf**.

-----

Exercise 3 comments:
--
<<<<<<< HEAD
=======
<<<<<<< HEAD
Here are our comments about what we did for exercise 3.
=======
>>>>>>> d436bc1
- The 2D Gaussian fit accurately reproduces the signal peak.  
- Residuals show no significant structure, confirming a good model.  
- Background normalization near unity indicates a consistent background model.  
- Signal yield computed from \(N_\text{sig} = A \pi \sigma_1 \sigma_2\).  
- Results and plots stored in **ex3.pdf**.
<<<<<<< HEAD
=======
>>>>>>> be093c7 (Comments)
>>>>>>> d436bc1

-----
