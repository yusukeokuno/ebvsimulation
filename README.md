# ebvsimulation
Monte-Carlo simulation for intragenic EBV deletions

1) Prepare deletions.txt, which is tab-delimited data containing sample ID and start and end positions of deletions.
ID1	102471	152317
ID10	143902	147118
ID11	121654	147332
ID11	64161	64895
ID12	135396	139827
ID12	35243	36015
The sample data contains some of the deletions identified in the original publication (Okuno Y et al. XXX XXX)

2) Run the script without any parameters:
php simluation.php
to obtain the simulation result (result.txt).
