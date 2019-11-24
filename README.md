## skyex package

How to install the package?

1. Open Rstudio
2. Run: 

```
install.packages("devtools")
library(devtools)
install_github("suelai/skyex")
library(skyex)
```

And that was it! :smiley:

includes two datasets: pairsManual https://dl.acm.org/citation.cfm?id=3340979, and 
                       restaurants https://www.cs.utexas.edu/users/ml/riddle/data.html
                       
Labeling algorithm based on:

[1]Isaj, Suela, Esteban Zimányi, and Torben Bach Pedersen. "Multi-Source Spatial Entity Linkage." Proceedings of the 16th International Symposium on Spatial and Temporal Databases. ACM, 2019.
APA	

[2]Isaj, Suela, Torben Bach Pedersen, and Esteban Zimányi."Multi-Source Spatial Entity Linkage."  https://arxiv.org/pdf/1911.09016v1.pdf

Includes 16 functions:

|**Module**                |               **Function name** |
| --- | --- |
|Blocking	               |             textual_blocking |
|Blocking	               |             spatial_blockcing |
|Blocking	               |             prefix_blocking |
|Blocking	               |             suffix_blocking |
|Pairwise comparison	   |               text_similarity |
|Pairwise comparison	   |               spatial_similarity |
|Pairwise comparison	   |               semantic_similarity |
|Labeling	               |             skyexf |
|Labeling	               |             skyexd |
|Analysis and visualization	 |         plot.skyexf.cutoffs |
|Analysis and visualization	 |         plot.skyexd.cutoffs |
|Analysis and visualization	 |         plot.skyexd.smooth |
|Analysis and visualization	 |         evaluate.skyex |
|Analysis and visualization	 |         plot.pairs2D |
|Analysis and visualization	 |         plot.pairs3D |
|Analysis and visualization	 |         plot.pairs.interactive.3D |

