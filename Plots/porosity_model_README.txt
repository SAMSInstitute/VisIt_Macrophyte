Re 10 results in numerical overflow for 4 and 16 towers of length 10/64.
This is due to the algorithm trying to fit a porosity value greater than 3000, resulting in
an overflow when trying to compute the exponential functions in the denomintaor of the D constant.

