# -Simplification-of-multimodal-histogram

### This is an academic R&D code.
### The multimodal histogram is cut and approximated for each single peak so that the error is below an arbitrary threshold.
### For the algorithm, refer to the following paper.
### [Linear Time Algorithm for Approximating a Curve by a Single-Peaked Curve](https://link.springer.com/chapter/10.1007/978-3-540-24587-2_3)

# Use
## fitting_hist.py
### input  > file name, epsilon

### output > epsilon, num of peak, L2-norm, processing time and figure

### please save the data as csv and work in the same directory.

## gaussian_fitting.ipynb
### This is a program that calculates the Gaussian approximation for each specified interval.
### Use the interval calculated in fitting_hist.py.
