sequenceClassifier
==================

Classifies genomic sequences

Needs NLTK package and numpy. 
Inputs are taxonomic hierarchies and training sequences for those hierarchies, 
as well as short unknown segments of DNA to classify.
Output the most likely taxonomic classification of this organsim based on n-gram 
perplexities calculated for the sequence as compared to the reference sequences.

Sample output:
The taxonomy of a human (Mammalia) is estimated as follows: 
The perplexity of this organism's genome when compared to Eukarya is 3.88734593957
The perplexity of this organism's genome when compared to Metazoa is 3.88779645956
The perplexity of this organism's genome when compared to Mammalia is 3.88250879765

The taxonomy of Staphylococcus aureus (Bacteria) is estimated as follows: 
The perplexity of this organism's genome when compared to Bacteria is 3.80761733257

The taxonomy of the common carp (ray-finned fish) is estimated as follows: 
The perplexity of this organism's genome when compared to Eukarya is 3.92766225573
The perplexity of this organism's genome when compared to Metazoa is 3.92960130126
The perplexity of this organism's genome when compared to Actinopterygii is 3.94319357506

The taxonomy of the giant amoeba is estimated as follows: 
The perplexity of this organism's genome when compared to Eukarya is 3.90615895904
The perplexity of this organism's genome when compared to Fungi is 3.90998119253