# Spatial Clustering of GInteractions

Merges spatially proximal chromatin loops into consensus interactions
using graph-based clustering. Loops are considered overlapping if both
anchors are within `gap` bp of each other. Each resulting cluster is
represented by the **union genomic range (min start to max end)**
spanning all its members.

## Usage

``` r
reduce_ginteractions(gi, gap = 1000)
```

## Arguments

- gi:

  A
  [`GInteractions`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)
  object.

- gap:

  Numeric. Maximum distance (in base pairs) allowed between anchors to
  consider two loops overlapping. Default: 1000.

## Value

A list with two elements:

- `gi`:

  Reduced
  [`GInteractions`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)
  object, one per cluster.

- `membership`:

  Integer vector indicating cluster assignment for each input loop.

Metadata columns include `cluster_id`, `n_members`, and averaged
`score`.

## Examples

``` r
# 1. Load example data (loops that are close to each other)
bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")

if (bedpe_path != "") {
  # Convert BEDPE to GInteractions object
  gi_raw <- bedpe_to_gi(bedpe_path)

  # 2. Run clustering
  # Merge loops if their anchors are within 1000bp
  res <- reduce_ginteractions(gi_raw, gap = 1000)

  # 3. Inspect results
  # The 'gi' element contains the merged consensus loops
  print(res$gi)

  # The 'membership' vector tells which original loop belongs to which cluster
  head(res$membership)

  # Check cluster sizes (how many loops were merged into each cluster)
  table(res$membership)
}
#> StrictGInteractions object with 300 interactions and 4 metadata columns:
#>         seqnames1             ranges1     seqnames2             ranges2 |
#>             <Rle>           <IRanges>         <Rle>           <IRanges> |
#>     [1]      chr1 109492652-109495835 ---      chr1 110644072-110649627 |
#>     [2]      chr1 116687216-116689300 ---      chr1 116844969-116856809 |
#>     [3]      chr1 116876465-116881376 ---      chr1 117819956-117830132 |
#>     [4]      chr1 146374138-146379497 ---      chr1 146442737-146445345 |
#>     [5]      chr1     1370034-1377149 ---      chr1     1469037-1474411 |
#>     ...       ...                 ... ...       ...                 ... .
#>   [296]      chr1     1011273-1015982 ---      chr1     1018256-1022868 |
#>   [297]      chr1     1342739-1364987 ---      chr1     1397119-1403614 |
#>   [298]      chr1 147271364-147274896 ---      chr1 147290756-147294529 |
#>   [299]      chr1 108639596-108642391 ---      chr1 109098569-109102247 |
#>   [300]      chr1   11011474-11014295 ---      chr1   12745038-12747130 |
#>         cluster_id n_members     score    n_reps
#>          <numeric> <integer> <numeric> <integer>
#>     [1]          1         1         1         1
#>     [2]          2         1         1         1
#>     [3]          3         1         2         1
#>     [4]          4         1         1         1
#>     [5]          5         1         3         1
#>     ...        ...       ...       ...       ...
#>   [296]        296         1        12         1
#>   [297]        297         1         7         1
#>   [298]        298         1         4         1
#>   [299]        299         1         1         1
#>   [300]        300         1         1         1
#>   -------
#>   regions: 385 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> 
#>   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#>  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#>  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#>  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#>  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 
#>   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 
```
