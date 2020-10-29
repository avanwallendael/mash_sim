# mash_sim

**Part 1: simulated reads**  
  
We simulated reads in python using [*ipcoal*](https://ipcoal.readthedocs.io/en/latest/) and [*toytree*](https://toytree.readthedocs.io/en/latest/). Scripts were run in a Jupyter notebook: mash_sim_github.ipynb  
We cleaned these reads for analysis in R: process_sim_polyploid.R  
These outputs were analyzed in two pipelines: 1: [*mash*](https://mash.readthedocs.io/en/latest/), and 2: alignment  
1. Simulated haploid and polyploid data analyzed with mash: mash_sim_github.bash  
2. Simulated haploid and polyploid data analyzed with alignment tools: polymiss_github.bash  
Next we removed reads from simulated polyploid data and tested the performance of each pipeline again: make_missing_github.R  

Finally we calculated distances for aligned data and compared results using polymiss_dist_github.R
  
**Part 2: real reads**
  
 We downloaded reads from published datasets for three studies in Panicum, Capsella, and Reynoutria. SRA projects: PRJNA622568, PRJNA299253, PRJNA574173.  
 We performed mash distance estimation ex. capsella_mash.bash
 Then cleaned and visualized results by incorporating published metadata ex. mash_clean_viz.R
