The examples below both are from the book "Biostatistics, 5th edition, Wiley,   by Wayne W. Daniel (2005)". 

Example 1 below, is an Oneway ANOVA example, where it is not reasonable assume equal variances.

In many applications, the p-values provided by alternative tests are in the same ballpark, leading to the same conclusion, 
but in Example 1 below the approximate Welch test almost leads to the rejection of the null hypotheses at 0.05 level, 
whereas Example 2 below leads to the same conclusion 

# Example 1
x1= c( 1,1,1,1,
       2,2,2,2,
       3,3,3,3,
       4,4,4)
y1= c( 71.8, 66.1, 67.6, 66.4, 42.8, 53.2, 56.1, 56.5, 
       72.5, 62.9, 58.9, 69.3, 47.1, 86.6, 56)
Example1 = data.frame( cbind(y1,x1)  )
out=oneway.anova(y1~x1,Example1)

# The output you get is the following.
                                                                         Method p-value
1                                                          Welch Test (approximate)   0.050
2     Generalized Test Equivalent to Fiducial Test (size-assured, but conservative)   0.100
3 Generalized Test Equivalent to Parametric Bootstratp Test (size close to intended   0.071

# Eaxmple 2
x1= c( 1,1,1,1,1,1,1,1,1,1,1,1,
       2,2,2,2,2,
       3,3,3,3,3,3,3,3)
y1= c( 92, 93, 74, 80.5, 76, 71, 75.5, 88.5, 93, 80.5, 83, 87, 79,
       78, 100, 76.5, 68, 81.5, 75, 76.5, 70.5, 69, 73.8, 74, 80)
Example2= data.frame( cbind(y1,x1)  )
# run as in # out=oneway.anova(y1~x1,data)
out=oneway.anova(y1~x1,Example2)
# The output you get is the following.
                                                                            Method p-value
1                                                          Welch Test (approximate)   0.056
2     Generalized Test Equivalent to Fiducial Test (size-assured, but conservative)   0.073
3 Generalized Test Equivalent to Parametric Bootstratp Test (size close to intended   0.060
