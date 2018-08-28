var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "BEDFiles.jl",
    "title": "BEDFiles.jl",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#BEDFiles.jl-1",
    "page": "BEDFiles.jl",
    "title": "BEDFiles.jl",
    "category": "section",
    "text": "Routines for reading and manipulating GWAS data in .bed files"
},

{
    "location": "index.html#Background-1",
    "page": "BEDFiles.jl",
    "title": "Background",
    "category": "section",
    "text": "Data from Genome-wide association studies are often saved as a PLINK binary biallelic genotype table or .bed file. To be useful, such files should be accompanied by a .fam file, containing metadata on the rows of the table, and a .bim file, containing metadata on the columns. The .fam and .bim files are in tab-separated format.The table contains the observed allelic type at n single-nucleotide polymorphism (SNP) positions  for m individuals.A SNP corresponds to a nucleotide position on the genome where some degree of variation has been observed in a population, with each individual have one of two possible alleles at that position on each of a pair of chromosomes. The three possible types that can be observed are: homozygous allele 1, coded as 0x00, heterzygous, coded as 0x10, and homozygous allele 2, coded as 0x11. Missing values are coded as 0x01.A single column - one SNP position over all m individuals - is packed into an array of div(m + 3, 4) bytes (UInt8 values)."
},

{
    "location": "index.html#Installation-1",
    "page": "BEDFiles.jl",
    "title": "Installation",
    "category": "section",
    "text": "This package requires Julia v0.7.0 or later, which can be obtained from https://julialang.org/downloads/ or by building Julia from the sources in the https://github.com/JuliaLang/julia repository.The package has not yet been registered and must be installed using the repository location. Start julia and use the ] key to switch to the package manager REPL(v0.7) pkg> add https://github.com/dmbates/BEDFiles.jl.git#master\n  Updating git-repo `https://github.com/dmbates/BEDFiles.jl.git`\n  Updating registry at `~/.julia/registries/Uncurated`\n  Updating git-repo `https://github.com/JuliaRegistries/Uncurated.git`\n Resolving package versions...\n  Updating `~/.julia/environments/v0.7/Project.toml`\n  [6f44c9a6] + BEDFiles v0.1.0 #master (https://github.com/dmbates/BEDFiles.jl.git)\n  Updating `~/.julia/environments/v0.7/Manifest.toml`\n  [6f44c9a6] + BEDFiles v0.1.0 #master (https://github.com/dmbates/BEDFiles.jl.git)\n  [6fe1bfb0] + OffsetArrays v0.6.0\n  [10745b16] + Statistics Use the backspace key to return to the Julia REPL."
},

{
    "location": "index.html#BEDFiles.BEDFile",
    "page": "BEDFiles.jl",
    "title": "BEDFiles.BEDFile",
    "category": "type",
    "text": "BEDFile\n\nRaw .bed file as a shared, memory-mapped Matrix{UInt8}.  The number of rows, m is stored separately because it is not uniquely determined by the size of the data field.\n\n\n\n\n\n"
},

{
    "location": "index.html#Loading-a-.bed-file-1",
    "page": "BEDFiles.jl",
    "title": "Loading a .bed file",
    "category": "section",
    "text": "The BEDFile struct contains the read-only, memory-mapped .bed file as a Matrix{UInt8}, along with m, the number of individuals.BEDFileFor convenience, two Int matrices, columncounts and rowcounts are allocated but not populated until used.The columns correspond to SNP positions. Rows of the internal matrix are packed values from groups of 4 individuals.julia> using BenchmarkTools, BEDFiles\n\njulia> const bf = BEDFile(BEDFiles.datadir(\"mouse.bed\"));\n\njulia> size(bf)      # the virtual size of the GWAS data - 1940 observations at each of 10150 SNP positions\n(1940, 10150)\n\njulia> size(bf.data) # the actual size of the memory-mapped matrix of UInt8s\n(485, 10150)As described above, a column, consisting of m values in the range 0x00 to 0x03, is packed into div(m + 3, 4) bytes.(The calculation div(i + 3, 4) or (i + 3) ÷ 4 occurs in some important loops in the code where it is evaluated as the equivalent, but somewhat faster, expression (i + 3) >> 2 that performs the integer division by 4 via shifting the number two bits to the right.)The virtual number of rows, m, can be given as a second argument in the call to BEDFile. If omitted, m is determined as the number of lines in the .fam file. Because the file is memory-mapped this operation is fast, even for very large .bed files.julia> @benchmark(BEDFile(BEDFiles.datadir(\"mouse.bed\")))\nBenchmarkTools.Trial: \n  memory estimate:  390.42 KiB\n  allocs estimate:  82\n  --------------\n  minimum time:     127.349 μs (0.00% GC)\n  median time:      135.760 μs (0.00% GC)\n  mean time:        150.356 μs (7.05% GC)\n  maximum time:     41.823 ms (99.33% GC)\n  --------------\n  samples:          10000\n  evals/sample:     1This file, from a study published in 2006, is about 5 Mb in size but data from recent studies, which have samples from tens of thousands of individuals at over a million SNP positions, would be in the tens or even hundreds of Gb range."
},

{
    "location": "index.html#Raw-summaries-1",
    "page": "BEDFiles.jl",
    "title": "Raw summaries",
    "category": "section",
    "text": "Counts of each the four possible values for each column are returned by counts.julia> counts(bf, dims=1)\n4×10150 Array{Int64,2}:\n  358   359  252   358    33   359    33  186   360  …    53    56    56    56    56    56    56    56\n    2     0    4     3     4     1     4    1     3      171   174   173   173   162   173   174   175\n 1003  1004  888  1004   442  1004   481  803  1002      186   242   242   242   242   242   242   242\n  577   577  796   575  1461   576  1422  950   575     1530  1468  1469  1469  1480  1469  1468  1467Column 2 has no missing values (code 0x01, the second row in the column-counts table). In that SNP position for this sample, 359 indivduals are homozygous allele 1 (G according to the .bim file), 1004 are heterozygous, and 577 are homozygous allele 2 (A).The counts by column and by row are cached in the BEDFile object. Accesses after the first are extremely fast.julia> @benchmark counts($bf, dims=1)\nBenchmarkTools.Trial: \n  memory estimate:  0 bytes\n  allocs estimate:  0\n  --------------\n  minimum time:     5.115 ns (0.00% GC)\n  median time:      5.298 ns (0.00% GC)\n  mean time:        5.244 ns (0.00% GC)\n  maximum time:     23.231 ns (0.00% GC)\n  --------------\n  samples:          10000\n  evals/sample:     1000"
},

{
    "location": "index.html#BEDFiles.bedvals",
    "page": "BEDFiles.jl",
    "title": "BEDFiles.bedvals",
    "category": "constant",
    "text": "BEDvals\n\nVector{Union{UInt8,Missings.Missing}} of the possible values in a BEDFile\n\n\n\n\n\n"
},

{
    "location": "index.html#Instantiating-as-a-count-of-the-second-allele-1",
    "page": "BEDFiles.jl",
    "title": "Instantiating as a count of the second allele",
    "category": "section",
    "text": "In some operations on GWAS data the data are converted to counts of the second allele, according toBEDFile count\n0x00 0\n0x01 missing\n0x10 1\n0x11 2This can be accomplished by indexing bedvalsbedvalswith the BEDFile or with a view of the BEDFile, producing an array of type Union{Missing,Int8}, which is the preferred way in v0.7 of representing arrays that may contain missing values.julia> bedvals[bf]\n1940×10150 Array{Union{Missing, Int8},2}:\n 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2  …  2         2         2       \n 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       \n 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2     2         2         2       \n 1  1  1  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       \n 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2     1         1         1       \n 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2  …  2         2         2       \n 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2     2         2         2       \n 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       \n 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       \n ⋮              ⋮              ⋮              ⋱                              \n 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2     2         2         2       \n 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2     2         2         2       \n 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2     2         2         2       \n 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2     2         2         2       \n 1  1  1  1  2  1  2  1  1  1  1  2  1  1  2  …  2         2         2       \n 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       \n 1  1  2  1  1  1  1  2  1  1  1  1  2  2  1     2         2         2       \n 1  1  1  1  1  1  1  2  1  1  1  1  2  2  1      missing   missing   missing\n 0  0  0  0  2  0  2  0  0  0  0  2  0  0  2     2         2         2       \n\njulia> sort(unique(ans))\n4-element Array{Union{Missing, Int8},1}:\n 0       \n 1       \n 2       \n  missing"
},

{
    "location": "index.html#Summary-statistics-1",
    "page": "BEDFiles.jl",
    "title": "Summary statistics",
    "category": "section",
    "text": "The package provides methods for the generics mean and var from the Statistics package.julia> mean(bf, dims=1)\n1×10150 Array{Float64,2}:\n 1.113  1.11237  1.28099  1.11203  …  1.8009  1.79966  1.79955  1.79943\n\njulia> var(bf, dims=1)\n1×10150 Array{Float64,2}:\n 0.469929  0.470089  0.462605  0.469365  …  0.223714  0.223818  0.223923These methods make use of the cached column or row counts and thus are very fastjulia> @benchmark mean(bf, dims=1)\nBenchmarkTools.Trial: \n  memory estimate:  79.39 KiB\n  allocs estimate:  2\n  --------------\n  minimum time:     37.777 μs (0.00% GC)\n  median time:      38.186 μs (0.00% GC)\n  mean time:        45.207 μs (14.53% GC)\n  maximum time:     43.815 ms (99.83% GC)\n  --------------\n  samples:          10000\n  evals/sample:     1The column-wise or row-wise standard deviations are returned by std.julia> std(bf, dims=2)\n1940×1 Array{Float64,2}:\n 0.6504997290784408\n 0.6379008244533891\n 0.6558172726141286\n 0.6532675479248437\n 0.6744432174014563\n 0.6519092298111158\n 0.6779881845456428\n 0.6955814098050999\n 0.6437566832989493\n ⋮                 \n 0.6613451651895259\n 0.6659810347614777\n 0.6274577846909379\n 0.6823658517777204\n 0.6695299551061924\n 0.710756592739754 \n 0.6387913736114869\n 0.6736492722732016\n 0.688855476425891 "
},

{
    "location": "index.html#BEDFiles.missingpos",
    "page": "BEDFiles.jl",
    "title": "BEDFiles.missingpos",
    "category": "function",
    "text": "missingpos(f::BEDFile)\n\nReturn a SparseMatrixCSC{Bool,Int32} of the same size as f indicating the positions with missing data\n\n\n\n\n\n"
},

{
    "location": "index.html#Location-of-the-missing-values-1",
    "page": "BEDFiles.jl",
    "title": "Location of the missing values",
    "category": "section",
    "text": "The positions of the missing data are evaluated bymissingposjulia> mp = missingpos(bf)\n1940×10150 SparseArrays.SparseMatrixCSC{Bool,Int32} with 33922 stored entries:\n  [702  ,     1]  =  true\n  [949  ,     1]  =  true\n  [914  ,     3]  =  true\n  [949  ,     3]  =  true\n  [1604 ,     3]  =  true\n  [1891 ,     3]  =  true\n  [81   ,     4]  =  true\n  [990  ,     4]  =  true\n  [1882 ,     4]  =  true\n  ⋮\n  [1848 , 10150]  =  true\n  [1851 , 10150]  =  true\n  [1853 , 10150]  =  true\n  [1860 , 10150]  =  true\n  [1873 , 10150]  =  true\n  [1886 , 10150]  =  true\n  [1894 , 10150]  =  true\n  [1897 , 10150]  =  true\n  [1939 , 10150]  =  true\njulia> @benchmark missingpos($bf)\nBenchmarkTools.Trial: \n  memory estimate:  1.81 MiB\n  allocs estimate:  19273\n  --------------\n  minimum time:     38.009 ms (0.00% GC)\n  median time:      38.161 ms (0.00% GC)\n  mean time:        38.761 ms (1.25% GC)\n  maximum time:     80.250 ms (52.34% GC)\n  --------------\n  samples:          129\n  evals/sample:     1So, for example, the number of missing data values in each column can be evaluated asjulia> sum(mp, dims=1)\n1×10150 Array{Int64,2}:\n 2  0  4  3  4  1  4  1  3  3  0  4  0  …  174  173  173  162  173  174  175although it is faster, but somewhat more obscure, to usejulia> view(counts(bf, dims=1), 2:2, :)\n1×10150 view(::Array{Int64,2}, 2:2, :) with eltype Int64:\n 2  0  4  3  4  1  4  1  3  3  0  4  0  …  174  173  173  162  173  174  175"
},

]}
