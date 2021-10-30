# PerFSeeB
<h1>Periodic full sensitivity spaced seeds</h1>

<h2>Introduction</h2>
We consider sequences of symbols <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt>. Suppose there is a long reference sequence (for a human genome it may have length of 3.2 billion symbols). There is also a set of short sequences (called <i>reads</i>), their size is around 50-300 symbols. We know that the reads are chunks of another long sequence which is in some way is close to the reference sequence. Our goal is to align those reads with respect to the reference sequence. 

Suppose, we have a reference sequence
<tt>ACGACAACCTTGTCGTTGGAGATCGGAAGAGCACACGTCTGAAC</tt>
and a read
<tt>TAGGTGCTCG</tt>. We shift the read with respect to the reference sequence, so

<table>
  <tr><th><tt>ACGACAACCTTGTCGTTGGAGATCGGAAGAGCACACGTCTGAAC</tt></th></tr>
  <tr><th><tt>_______________TAGGTGCTCG___________________</tt></th></tr>
  <tr><th><tt>_______________1011010111___________________</tt></th></tr>
</table>

We use Hamming distance to measure similarity of these two seqeunces. If two symbols at the same position are identical, then the distance is 0, otherwise the distance is 1. The total distance is the sum of elementwise distances. Ideally, we would like to position the second seqeunce in such a way, so the distance attains its minimum value (or is within some ranges of values). 

The standard approach in case of a long seqeunce is to consider its smaller chunks and find their poistions within the reference sequence. After a read is pre-aligned, we calculate the final similarity score based on Hamming or other distances. The goal of this project is to pre-align reads (find candidate positions within the reference sequence).

Reads may have various mutations. Therefore there can be some mismatches (when the reference seqeunce and a read have different symbols at same positions), insertions or deletions of symbols. Here we consider only mismatches. 

If there are several mutations, then using contiguous chunks of symbols to identify candiate positions may not work. Therefore it is better to consider chunks with some symbols to be ignored when sequences are compared. 

For this purpose we used a <i>seed</i>, a sequence of 1 and 0, of length n<sub>S</sub> (the total number of all symbols). By seed's <i>weight</i> we call the total number of 1s in the seed.

<table>
  <tr><th>Seq 1</th><th><tt>TTGGAGATCG</tt></th></tr>
  <tr><th>Seq 2</th><th><tt>TAGGTGCTCG</tt></th></tr>
  <tr><th><tt>101101____</tt></th><th><tt>TTGGAG____</tt></th></tr>
  <tr><th><tt>_101101___</tt></th><th><tt>_TGGAGA___</tt></th></tr>
  <tr><th><tt>__101101__</tt></th><th><tt>__GGAGAT__</tt></th></tr>
  <tr><th><tt>___101101_</tt></th><th><tt>___GAGATC_</tt></th></tr>
  <tr><th><tt>____101101</tt></th><th><tt>____AGATCG</tt></th></tr>
</table>



