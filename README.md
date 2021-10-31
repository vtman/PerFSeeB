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

For this purpose we used a <i>seed</i>, a sequence of 1 and 0, of length n<sub>S</sub> (the total number of all symbols). By seed's <i>weight</i> we call the total number of 1s in the seed. It is assumed that a seed starts and ends with 1-element.

Let us conider a seed <tt>1111</tt> of weight 4. We may seed that this seed cannot be used to find candidate position, since there are no all same symbols defiend by a shifted seed.

<table>
  <tr><th>Seed</th><th>Seq 1</th><th>Seq 2</th><th></th></tr>
  <tr><th></th><th><tt>TTGGAGATCG</tt></th><th><tt>TAGGTGCTCG</tt></th><th></th></tr>
  <tr><th><tt>1111______</tt></th><th><tt>TTGG______</tt></th><th><tt>TAGG______</tt></th><th>&#10060;</th></tr>
  <tr><th><tt>_1111_____</tt></th><th><tt>_TGGA_____</tt></th><th><tt>_AGGT_____</tt></th><th>&#10060;</th></tr>
  <tr><th><tt>__1111____</tt></th><th><tt>__GGAG____</tt></th><th><tt>__GGTG____</tt></th><th>&#10060;</th></tr>
  <tr><th><tt>___1111___</tt></th><th><tt>___GAGA___</tt></th><th><tt>___GTGC___</tt></th><th>&#10060;</th></tr>
  <tr><th><tt>____1111__</tt></th><th><tt>____AGAT__</tt></th><th><tt>____TGCT__</tt></th><th>&#10060;</th></tr> 
  <tr><th><tt>_____1111_</tt></th><th><tt>_____GATC_</tt></th><th><tt>_____GCTC_</tt></th><th>&#10060;</th></tr>
  <tr><th><tt>______1111</tt></th><th><tt>______ATCG</tt></th><th><tt>______CTCG</tt></th><th>&#10060;</th></tr>
</table>


However, if we consider a seed <tt>101101</tt> (it also has weight 4, we call this seed a <i>spaced seed</i>), then there is one position when two sequences are fully matched.

<table>
  <tr><th>Seed</th><th>Seq 1</th><th>Seq 2</th><th></th></tr>
  <tr><th></th><th><tt>TTGGAGATCG</tt></th><th><tt>TAGGTGCTCG</tt></th><th></th></tr>
  <tr><th><tt>101101____</tt></th><th><tt>T_GG_G____</tt></th><th><tt>T_GG_G____</tt></th><th>&#10004;</th></tr>  
  <tr><th><tt>_101101___</tt></th><th><tt>_T_GA_A___</tt></th><th><tt>_A_GT_C___</tt></th><th>&#10060;</th></tr>
  <tr><th><tt>__101101__</tt></th><th><tt>__G_AG_T__</tt></th><th><tt>__G_TG_T__</tt></th><th>&#10060;</th></tr>
  <tr><th><tt>___101101_</tt></th><th><tt>___G_GA_C_</tt></th><th><tt>___G_GC_C_</tt></th><th>&#10060;</th></tr>
  <tr><th><tt>____101101</tt></th><th><tt>____A_AT_G</tt></th><th><tt>____T_CT_G</tt></th><th>&#10060;</th></tr>
</table>

<h2>checkSeedClassic/checkSeed128: Check if a seed is valid</h2>
Suppose we are given a seed of length <tt>L</tt>. We also know the maximum number <tt>m</tt> of mismatches and read's length <tt>r</tt>. We create <tt>T = r-L+1</tt> rows and pad the seed with 0s (just adding extra zero to the left for each new row and removing one zero from the right). A seed is valid if for any arbitrary <tt>m</tt> columns of the matrix there is at least one row such that all corresponding elements are zeros.

<h3>Parameters</h3>

<ol>
  <li>Read length (integer)</li>
  <li>Number of mismatches (integer)</li>
  <li>Seed (string of 1 or 0)</li>
</ol>

<tt>checkSeed128.exe 15 3 10110001</tt>

For the output we get the matrix and information about seed's validity. If the seed is not valid we also get a list of columns when the requirements are not met. For example, we get columns 4, 8, 10 and the matrix (extra separator | is added for convinience).

<table>
  <tr><th><tt>101|1000|10|000000</tt></th></tr>
  <tr><th><tt>010|1100|01|000000</tt></th></tr>
  <tr><th><tt>001|0110|00|100000</tt></th></tr>
  <tr><th><tt>000|1011|00|010000</tt></th></tr>  
  <tr><th><tt>000|0101|10|001000</tt></th></tr>
  <tr><th><tt>000|0010|11|000100</tt></th></tr>
  <tr><th><tt>000|0001|01|100010</tt></th></tr>
  <tr><th><tt>000|0000|10|110001</tt></th></tr>
  <tr><th><tt>___|X___|X_|X_____</tt></th></tr>
</table>


checkSeed128 is a faster version based on SIMD instructions. checkSeedClassic and checkSeed128 have the same input parameters.

<h2>iterSeed: Spaced seeds generated iteratively</h2>

We have a tool to check if a seed is valid. Of course, for a given length of a seed we may generate all possible spaced seeds. For length 1, there is only one seed <tt>1</tt>, for length 2, there is also one seed <tt>11</tt>, for length 3, there are 2 seeds (<tt>101</tt>, <tt>111</tt>), for length 4, there are 4 seeds (<tt>1001</tt>, <tt>1011</tt>, <tt>1101</tt>, <tt>1111</tt>), etc. So, for length k, there are 2<sup>k-2</sup> seeds of length k, therefore there are 2<sup>k-1</sup> seeds of length not more than k. We need to reduce the number of possible seeds to be validated. 
