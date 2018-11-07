# Generating the Genetic Relationship Matrix from a SnpArray

The simplest case is equivalent to forming the `m × m` symmetric matrix `X*X'` where
`X` matrix is the centered, scaled and imputed major allele counts in each column.

To break this down, for each column the operation corresponds to:
- Copy the column into a floating-point vector according to, for the additive model,

| UInt8 code | value |
|:----------:| -----:|
| 0x00 | 0.0 |
| 0x01 | missing | 
| 0x02 | 1.0 |
| 0x03 | 2.0 |

- Evaluate the mean, `μ`, and the scaling factor, `σ`, for the column skipping missing values.  For the additive model the scaling factor is `σ = √(μ(1 - μ / 2))`
- Imputation: Replace missing values with `μ`
- Centering: Subtract the mean `μ` from each value.
- Scaling: Divide by `σ`

Reversing the order of imputation and centering helps as imputation becomes replacing missing values by zero.

In a column generated in this way there will be at most four distinct values.
It is more effective to generate the four possible values first, the populate the vector and add its contribution to the `X*X'` product.  These four values are `[-μ, 0, 1 - μ, 2 - μ] ./ σ`.