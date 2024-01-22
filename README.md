# Bayesian Peak Fitting in R
Using Bayesian statistics to fit overlapping peaks for spectroscopy analysis. Determines Gaussian and non-Gaussian fitting parameters with the [quap](https://www.rdocumentation.org/packages/rethinking/versions/2.13/topics/quap) (quadratic approximate posterior distribution) function, part of the [rethinking](https://github.com/rmcelreath/rethinking) package.

Non-Gaussian peaks, e.g. due to target oxidation, are assumed to be fit with Exponentially-Modified Gaussians (EMGs).

## Examples from my [Ph.D. dissertation](https://github.com/Willcarr99/PhDThesis):

</br>

![5deg_9400keVGroup](https://github.com/Willcarr99/BayeSpec-Analysis/assets/55559733/e3e59b3d-784d-49fe-8ae2-5e8b0d4ef19e)
A Gaussian multiplet consisting of 4 peaks. Deuteron counts from the custom focal-plane detector at the TUNL Enge Split-Pole Spectrograph for a <sup>39</sup>K(<sup>3</sup>He, d)<sup>40</sup>Ca nuclear reaction experiment.

</br></br></br>

![EMG_Multiplet](https://github.com/Willcarr99/BayeSpec-Analysis/assets/55559733/53b9dabb-a7d9-4506-9f44-ef377499be6b)
A non-Gaussian multiplet with long tails as a result of target oxidation during the <sup>39</sup>K(<sup>3</sup>He, d)<sup>40</sup>Ca experiment. Different colored peaks correspond to deuterons from different reactions. The orange peak is from the contaminant <sup>13</sup>C(<sup>3</sup>He, d)<sup>14</sup>N. Code currently supports peaks from up to 3 different reactions.