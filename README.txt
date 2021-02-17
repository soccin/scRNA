# scRNA Version 4

## Branch proj/p11533 [AbdelwaO/ChenS7/Proj_11533]

10X data

## Current pipeline

Multiple stages

- `doSeuratV5_01.R`: Initial QC to check filtering paramters and check cell cycle regression plots to see if cell cycle regression is needed.

- `doSeuratV5_02.R`: Stage 2, redo the QC-filter (but can reset parameters) and now set whether to regress out cell cycle or not.

