# BEDFiles.jl
Routines for reading and manipulating GWAS data in .bed files

## Installation

At present this repository is in the form of a Julia project which can only
be used with Julia v0.7 and later.  A beta of v0.7.0 is currently available on the
[download site](https://julialang.org/downloads).  Install that first.

Clone this repository then, in the directory that you have cloned, start julia.
Switch to the package REPL-mode using the key `]` then activate and instantiate
this project.
```julia
(BEDFiles) pkg> activate .

(BEDFiles) pkg> instantiate
  Updating registry at `~/.julia/registries/Uncurated`
  Updating git-repo `https://github.com/JuliaRegistries/Uncurated.git`

```

Use the backspace key to return to the Julia REPL.

## Loading a .bed file

The `BEDFile` struct contains a memory-mapped `.bed` file as a `Matrix{UInt8}`
(unsigned 8-bit integers) along with `m`, the number of virtual rows.  The columns
correspond to SNP positions.  Rows of the internal matrix are packed values from groups of 4 subjects.

```julia
julia> using BEDFiles

julia> const bf = BEDFile("./data/mouse/alldata.bed")
BEDFile(UInt8[0xba 0xba … 0xff 0xff; 0xab 0xab … 0xfe 0xfe; … ; 0xbb 0xbb … 0xff 0xff; 0x2a 0x2a … 0xdf 0xdf], 1940)

julia> size(bf)  # the virtual size of the GWAS data
(1940, 10150)

julia> size(bf.data) # the actual size of the memory-mapped matrix of UInt8s
(485, 10150)
```

In general, data on `m` subjects produces a memory-mapped matrix with `div(m + 3, 4)` rows
```
julia> div(1940 + 3, 4)
485

julia> (1940 + 3) >> 2  # an equivalent, somewhat faster, calculation
485
```

The virtual number of rows can be given as a second argument in the call to `BEDFile`.
If it is omitted, the number of lines in the file with the same name but a `.fam` extension
is used.

Because the file is memory-mapped this operation is fast, even for very large `.bed` files.
```julia
julia> @time BEDFile("./data/mouse/alldata.bed");
  0.000339 seconds (48 allocations: 10.500 KiB)
```

## Raw summaries

Each location in the virtual array is a value between 0 and 3 where, somewhat unintuitively, 1
is the missing value code, 0 indicates homozygous first allele, 2 indicates heterozygous and 3 is
homozygous second allele.  (The information on first and second allele are in the `.bim` file.)
See [`https://www.cog-genomics.org/plink2/formats#bed`]

Counts of each type for each column are returned by `columncounts`.

```julia

julia> columncounts(bf)
4×10150 Array{Int64,2}:
  358   359  252   358    33   359    33  186   360  …    53    56    56    56    56    56    56    56
    2     0    4     3     4     1     4    1     3      171   174   173   173   162   173   174   175
 1003  1004  888  1004   442  1004   481  803  1002      186   242   242   242   242   242   242   242
  577   577  796   575  1461   576  1422  950   575     1530  1468  1469  1469  1480  1469  1468  1467
```

This operation also is reasonably fast
```julia
julia> Threads.nthreads()
1

julia> @time columncounts(bf);
  0.250703 seconds (30.46 k allocations: 1.549 MiB)
```
and is multi-threaded.  That is, when the environment variable `JULIA_NUM_THREADS` has been set before starting Julia the operation is spread over this number of threads.
```julia
julia> Threads.nthreads()  # after restarting with JULIA_NUM_THREADS=4
4

julia> @time columncounts(bf);
  0.071961 seconds (30.24 k allocations: 1.540 MiB)
```

## Operating on columns

The `BEDColumn` type represents a column of the data and provides efficient access to the data, especially when used as an *iterator*.
The magic of [fused vectorized operations](https://docs.julialang.org/en/latest/manual/performance-tips/#More-dots:-Fuse-vectorized-operations-1) and the iterator optimizations in v0.7 generally provide the fastest way to operate on a column.
```julia
julia> bf1221 = BEDColumn(bf, 1221);

 0x02
 0x02
 0x02
 0x02
 0x03
 0x03
 0x02
 0x03
 0x03
    ⋮
 0x03
 0x03
 0x03
 0x03
 0x03
 0x03
 0x03
 0x02
 0x03
```
Although it appears that the column has been expanded to 1940 `UInt8` values,
the column is in fact just a one-dimensional view of the column of the compressed memory-mapped array in `bf`.
```julia
julia> typeof(bf1221.data)
SubArray{UInt8,1,Array{UInt8,2},Tuple{Base.Slice{Base.OneTo{Int64}},Int64},true}

julia> length(bf1221.data)
485
```

In addition there is a fast way of iterating over the 1940 values in the GWAS array.
As a trivial example
```julia
julia> extrema(bf1221)  # the minimum and maximum values in the column
(0x00, 0x03)

julia> @time extrema(bf1221);
  0.000053 seconds (5 allocations: 176 bytes)
```

A more practical example is to find all the positions in this column with missing values.
Recall that `0x01` indicates a missing value.
```julia
julia> findall(isone.(bf1221))
3-element Array{Int64,1}:
  676
  990
 1044

julia> @time findall(isone.(bf1221));
  0.000060 seconds (10 allocations: 4.906 KiB)
```
To break this down, the `isone` function applied to a number returns a `Bool`.
```julia
julia> isone(2)
false

julia> isone(1)
true
```
*Dot vectorization* in Julia means that `isone.` applied to an iterator is itself an iterator formed by applying `isone` to each element in turn.
```julia
julia> show(isone.(1:5))
Bool[true, false, false, false, false]
```
and `findall` applied to an iterator of `Bool` values returns the indices of the `true`
values.

The `countcols` method and other column-oriented operations are performed in parallel by assigning a chunk of columns to each thread.

## Instantiating as a count of the second allele

In some operations on GWAS data the data are converted to counts of the second allele.
This is accomplished by indexing `bedvals` with the column returning a vector of type
`Union{Missing,UInt8}` which is the preferred way in v0.7 of representing possibly data
vectors that may contain missing values.
```julia
julia> bedvals[bf1221]
1940-element Array{Union{Missing, UInt8},1}:
 0x01
 0x01
 0x01
 0x01
 0x02
 0x02
 0x01
 0x02
 0x02
    ⋮
 0x02
 0x02
 0x02
 0x02
 0x02
 0x02
 0x02
 0x01
 0x02

```
The mean for each column in this representation is returned by
```julia
julia> mean(bf, dims=1)
1×10150 Array{Float64,2}:
 1.113  1.11237  1.28099  1.11203  …  1.8009  1.79966  1.79955  1.79943
```
This could be done by applying `mean(skipmissing(col))` to each column.
```julia

## Location of the missing values

Some operations subsetting the rows to only those with complete data.  Discovering which
rows do not have any missing data could be done by iterating across the rows but it is generally faster to iterate over columns.

Recall that the missing value indicator is 1.  The rows with missing values in column `j` are determined as
```julia
julia> findall(isone.(BEDColumn(bf, 1221)))
3-element Array{Int64,1}:
  676
  990
 1044
```
One way to determine the rows with any missing data is convert the row numbers of missing values to a `BitSet`
and take the union over all the columns. 
```julia
julia> BitSet(findall(isone.(bf1221)))
BitSet([676, 990, 1044])

julia> anymsng = mapreduce(j -> BitSet(findall(isone.(BEDColumn(bf, j)))), union, 1:size(bf,2));

julia> @time mapreduce(j -> BitSet(findall(isone.(BEDColumn(bf, j)))), union, 1:10150);
  0.219699 seconds (377.41 k allocations: 65.061 MiB, 3.27% gc time)
```
The bad news here is that only a few row (288, to be exact) don't have any missing data (recall that there are 1940 rows)
```julia
julia> length(anymsng)
1652
```

An alternative is to use `missingpos` to obtain a `SparseMatrixCSC` indicating the positions of the missing values
```julia
julia> msngpos = missingpos(bf)
1940×10150 SparseMatrixCSC{Int8,Int32} with 33922 stored entries:
  [702  ,     1]  =  1
  [949  ,     1]  =  1
  [914  ,     3]  =  1
  [949  ,     3]  =  1
  [1604 ,     3]  =  1
  [1891 ,     3]  =  1
  [81   ,     4]  =  1
  [990  ,     4]  =  1
  [1882 ,     4]  =  1
  ⋮
  [1848 , 10150]  =  1
  [1851 , 10150]  =  1
  [1853 , 10150]  =  1
  [1860 , 10150]  =  1
  [1873 , 10150]  =  1
  [1886 , 10150]  =  1
  [1894 , 10150]  =  1
  [1897 , 10150]  =  1
  [1939 , 10150]  =  1
julia> @time missingpos(bf);
  0.217169 seconds (60.93 k allocations: 48.183 MiB, 1.25% gc time)
```
which can then be reduced using matrix multiplication
```julia
julia> findall(iszero, msngpos * ones(Int, size(msngpos, 2)))'  # rows with no missing data
1×288 LinearAlgebra.Adjoint{Int64,Array{Int64,1}}:
 2  5  11  22  30  37  38  53  56  59  63  65  67  …  1880  1885  1902  1904  1910  1915  1926  1928

julia> @time findall(iszero, msngpos * ones(Int, size(msngpos, 2)))';
  0.000128 seconds (22 allocations: 103.422 KiB)
```
