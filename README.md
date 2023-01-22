# scRNA Version 4

## Branch dev/2023

Major refactor of workflow



Now does SCTransform and Integrate

Needs R>=4.x. Run `CMD.Setup_R-4.x` to load it

## Current pipeline

Multiple stages

- `doSeuratV5_01.R`: Initial QC to check filtering paramters and check cell cycle regression plots to see if cell cycle regression is needed.

- Stage 2: Multiple parts now:

    - `doSeuratV5_02a.R` Now does SCTransform and Integrate with Normalize and CC regression
    - `doSeuratV5_02b.R` PCA and clustering


