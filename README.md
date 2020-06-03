Welcome to the aerosol refractive index confidence explorer (RICE) 1.0. RICE is a tool that helps understand uncertainties associated with aerosol refractive index (m) calculations by calculating the probability that nearby m values could have caused the originally retrieved m value.
To use RICE, you will need:
1. 2D particle size distribution time series wave (dN/dLogD, with size bins in rows, time points in columns, and values as the number concentration)
2. 1D particle size bin location wave, with the midpoint diameter of each size bin used to measure the size distribution. This wave must have the same number of rows as the size distribution.
3.1D single wavelength absorption and scattering coefficient waves. These waves must both have the same number of points as the columns of the size distribution.
4.The relative measurement uncertainties (expressed as fraction of the observation) of the optical measurements, size distribution bin locations**, and size distribution number concentrations**. Ideally, these should be calculated from the specific measurement systems used, not literature- based values.
5. If using the chi-squared method to calculate m, the standard deviations associated with the averaging of the absorption and scattering coefficient waves.
**In RICE 1.0 size distribution and number concentration uncertainties are treated as uniform over the whole size distribution.
To use RICE, follow the following steps:

1.	Calculate the m time series, ideally using the built in “Refractive Index Calculation” function because this performs the calculation using the same methods as RICE. This can be found in the RICE drop down menu.
2.	Open RICE 1.0 from the RICE drop down Menu
3.	Select/input the m times series wave and the other data waves and variables detailed above.
4.	Choose the number of iterations (recommended 10), i*j value (recommended 64, corresponding to 9 discrete n's and 9 discrete k's) for establishing the sampling space.
5.	Choose the number of iterations (recommended 100) and i*j value (recommended 144, corresponding to 13 discrete n's and 13 discrete k's) for calculating the final solutions.
6.	Choose the confidence interval size of the output and whether or not you would like to keep all the raw data from the RICE run (slightly slows down RICE and makes significantly more output files).
Note: RICE runs relatively slow on personal computers, it's not uncommon for RICE to take >10 minutes per time point. Higher i x j and iteration values will exponentially increase RICE's runtime.
After RICE has completed, you'll have the following data:
1.	Confidence Interval Width values of n and k (n_CI_width and k_CI_with, respectively).
2.	Lower bounds and upper bounds of the calculated n and k confidence intervals.
3.	RICE_Output, a 2D wave which contains the sigmoidal fit data allowing for the recalculation of confidence intervals and examination of the goodness of fits.
4.	Flag waves Space_n, Space_k, and Count_. Space_n and Space_k attempt to determine if the established sampling space was inappropriate. A code of 0 suggests the sampling space includes all possible solutions. A code of 3 suggests the sample space is too small to yield an accurate RICE result. Count is a measure of how many of the attempted m calculations yielded solutions with values within the i x j matrix. When count is 0 the sampling space is likely appropriate. A count value of 1 suggests many of the iterations did not yield useable solutions, so the results should be handled with caution.
RICE Reanalysis:
Once a RICE result has been obtained, intervals at different confidence levels can be obtained by using the RICE Reanalysis menu function. This allows the user to calculate new confidence intervals from any RICE_Output wave. Before using the reanalysis function, set the current data folder to the folder which contains the RICE results.

