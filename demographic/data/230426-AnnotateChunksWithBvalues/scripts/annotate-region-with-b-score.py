"""
Annote region with B-score
The B-score is the mean of the overlapping reiong scores.
"""

import pybedtools

bmap = pybedtools.BedTool('data/caad-bestfit-GRCh38.bmap.bed')
region = pybedtools.BedTool('chunk_283.bed')

# which regions in bmap overlap with region?
overlapping = bmap.intersect(region, wa=True, wb=True)
