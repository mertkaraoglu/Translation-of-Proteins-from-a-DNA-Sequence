# Translation-of-Proteins-from-a-DNA-Sequence

Given a DNA sequence, following the transcription process, proteins that would be translated are found.
PS: This code is based on very simple rules that are mostly underfitting the real world results.

**Transcription and Intron Extraction**
1. Find the initial point of the first exon "ATG"
2. Find the Branch Points of the introns
3. Locate the 5' and 3' intron sites of each branch point (GT - AG)
4. Extract and exclude the introns 
5. Search for the last point "TGA" and end the last exon there. 
