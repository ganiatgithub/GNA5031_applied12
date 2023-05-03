# GNA5031_applied12

Initial planning:

1. Three isolates grown under antibiotic selection. Give them predicted protein fasta, plasmid nucleotide sequence and 16S rRNA nucleotide file each
2. Align proteins with CARD database via DIAMOND
3. Use results to work out which ARGs are present, and what organism the database says they are from
4. BLAST the 16S rRNA nucleotide file for the isolate to check identity. They match
5. One has no ARGs in the genome. Why? Maybe it's a plasmid-borne ARG. Let's check
6. Align the plasmid sequence against the CARD database. We find an ARG that does not match the taxonomy

Quesion: why doesn't taxonomy match? Because the plasmid has come from something else.
