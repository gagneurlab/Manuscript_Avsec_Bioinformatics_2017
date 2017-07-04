# Computing the distance features in the clip data

Branstormed with Julien:

Extending the model with GAM:
- scalars could be positions to RNA landmarks or 
- clip signal to the other 30 factors
- signed distance to the nearest clip peak for each of the other factor
  - GAM smooth
- **Any other useful features?**
- Average occupancy in the region
- distance to the nearest peak

## Features to include
Distance to nearest:
  * intron-exon 
  * exon-intron
  * TSS
  * poly-a
  * startcodon 
  * stop-codon

## Open questions
- which isoform to take? - just take the nearest feature


## Plan
- [x] Get the landmarks you're interested in from a GTF
   - use standard hg19 GTF file
   - how was it done in the Beth et al paper?
- [x] For each location, compute the nearest boundaries
- [x] Save the distances
- [x] Visualize the distributions
   - sort by the test statistic
   - is log-maybe the more appropriate transforman
- [x] Run the same processing pipeline also for the eClip data
