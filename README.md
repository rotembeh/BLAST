BLAST overview
From Wiki:

“In bioinformatics, BLAST for Basic Local Alignment Search Tool is an algorithm for comparing primary biological sequence information, such as the amino-acid sequences of different proteins or the nucleotides of DNA sequences. A BLAST search enables a researcher to compare a query sequence with a library or database of sequences, and identify library sequences that resemble the query sequence above a certain threshold.”

For this assignment, BLAST will be defined according to this variant:

The input is a Data Base sequence, and one or more query sequences.

Data Base pre-process:

Mapping each w-mer (a substrings of length w) from the db to its index in the db string. Mapping can be done in various ways – you can implement it using a data structure of your choice. This is done once for the db sequence, and used for multiple queries.

Query processing - Finding hits:

For each w-mer in the query ("sliding window"), make a list of neighbour words. Reminder: A neighbour word 𝑠′ of a w-mer 𝑠 is a word of length |s|=w, with score(s,s') ≥  T (using a scoring matrix). T is a chosen threshold.

For each neighbor, find the indexes of its appearances in the text (hits), using the data structures constructed in stage 1. A hit is called an HSP (High-Scoring Segment Pair).

Extension of HSPs:

Extend the hits found in stage 2 to obtain MSPs (Maximal Scoring Pairs). Extension stops when the extended HSP score drops X below the maximal score obtained so far.

Reporting MSPs:

Report top scoring MSPs of the query, after removing duplicate MSPs.
