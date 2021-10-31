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

If there is a valid seed, then all its subseeds are also valid. For example, seed <tt>111001011</tt> is valid for r=15, m=2. Subseeds like <tt>111</tt>, <tt>1011</tt>, <tt>111001</tt>, <tt>1001001</tt> are also valid. Therefore for a step k we may consider all valid spaced seeds of length less than k and pad them with all 0-elements and one end 1-element from the right. So, to find seed <tt>1101001101</tt> (it is also a valid seed for the abovbe parameters) we use seed <tt>11010011</tt> pad it with one <tt>0</tt> and one <tt>1</tt>. We may also perform an additional check: by removing the leftest 1-element of the new seed (and all neighbouring 0-elements) we should get a valid seed. Therefore we need to check if 
<tt>101001101</tt> is in the list of valid seeds generated before. There is no need to check other subseeds of <tt>1101001101</tt>, since either <tt>11010011</tt> or <tt>101001101</tt> contain them.

<h3>Parameters</h3>

<ol>
  <li>Output folder</li>
  <li>Number of mismatches</li>
  <li>Length of reads</li>
	<li>Minimum weight</li>
</ol>

Since there are usually a lot of spaced seeds generated in this way, we try to report only seeds of large weights. Therefore we specify the minimum weight required for a seed to be reported.

<tt>iterSeed.exe C:\MyFolder 15 2 6</tt>

As an output we get the following spaced seeds: <tt>111010011</tt>, <tt>111001011</tt>, <tt>110100111</tt>, <tt>110010111</tt>, <tt>1101001101</tt>, <tt>1011001011</tt>. These seeds have the maximum weight possible (6). Note that if there is a valid seed, then its flipped version should also be in the list, i.e. seeds <tt>110100111</tt>, <tt>110010111</tt> are both in the list. There may also be longer seeds of smaller weight. For example, the maximum length for seeds of weight 6 is 10 (seeds <tt>1101001101</tt> and <tt>1011001011</tt>), however for seeds of weight 5 we get <tt>1001001001001</tt> of length 13.

We have generated spaced seeds for number of mismatches from 2 to 8 and reads' lengths from 10 to 50 (different ranges for different numbers of mismatches). They are in <b>iterSeed</b> folder.

<h2>Seeds of maximum weight</h2>
For practical applications it is better to use seeds of maximum weight. We have four letters <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt>. Let us assume that their chance to be in a seqeunce is the same and does not depend on neighbouring symbols (these assumptions are not completly true). So, if a chance to find a pattern of length k within a reference sequence is P<sub>k</sub>, then the chance to find a pattern of length (k+1) is P<sub>k</sub>/4. Therefore seeds of higher weight allow us to process 4 times less candidate positions within a reference sequence.

For a given length of reads and a number of mismatches, there may be several seeds of maximum weight. For example, for r=45, m=8 and w=6 we get several valid seeds including <tt>1111011</tt>, <tt>1011110001</tt>, <tt>11011000101</tt>, <tt>1001001000001000000001001</tt>, <tt>1100000001000000001100000001</tt>. This means that if we choose one of the shortest seeds (<tt>1111011</tt>) of length 7, then we need to consider 45-7+1=39 chunks of a read (of length 45) and use them to find corresponding candidate positions within a reference seqeunce. However, if we use seed <tt>1100000001000000001100000001</tt> of length 28, the number of chunks becomes 45-28+1=18 (almost 2 times less). So, using longest seeds among seeds of maximum weight can reduce processing times.

Of course, there may several lengths of reads when w is the maximum weight. If a seed is valid for a read of length r, then it is also valid for a read of length (r+1). Therefore, while there may be several seeds valid for various lengths of reads, we pick up only those valid for shortest lengths. For around 80% of seeds obtained using the iterative procedure we may see that the best seeds (longest seeds of maximum weight valid for shortest reads) have a periodic structure: an integer number n<sub>b</sub> of blocks of length <i>T</i> and a "remainder" (first n<sub>d</sub> elements of the block), so the total length is n<sub>s</sub> = n<sub>b</sub> T + n<sub>d</sub>, and the following formula is valid

r = n<sub>s</sub> + T - 1

For example, seed <tt>1110100000000111010000000011101</tt> is valid for m=4, r=43 and can be split up as
<table>
	<tr><th><tt>1110100000000</tt></th><th><tt>1110100000000</tt></th><th><tt>11101</tt></th></tr>
</table>
We seed that T=13, n<sub>b</sub>=2, n<sub>d</sub>=5, so n<sub>s</sub>=31 and 43 = 31 + 13 - 1.

We try to find possible blocks such that we are able to form seeds of the given structure.

<h2>Periodic blocks</h2>

If the formula above is valid for a periodic seed, then to validate the seed it is enough to validate its periodic block. For example, we have seed <tt>110010111001011100101110010111</tt>. We create 7 rows such that

<table>
	<tr><th><tt>110010111001011100101110010111000000</tt></th></tr>
	<tr><th><tt>011001011100101110010111001011100000</tt></th></tr>
	<tr><th><tt>001100101110010111001011100101110000</tt></th></tr>
	<tr><th><tt>000110010111001011100101110010111000</tt></th></tr>
	<tr><th><tt>000011001011100101110010111001011100</tt></th></tr>
	<tr><th><tt>000001100101110010111001011100101110</tt></th></tr>
	<tr><th><tt>000000110010111001011100101110010111</tt></th></tr>
</table>

We split the rows to form
<table>
	<tr><th><tt>1100101</tt></th><th><tt>1100101</tt></th><th><tt>1100101</tt></th><th><tt>1100101</tt></th><th><tt>11000000</tt></th></tr>
	<tr><th><tt>0110010</tt></th><th><tt>1110010</tt></th><th><tt>1110010</tt></th><th><tt>1110010</tt></th><th><tt>11100000</tt></th></tr>
	<tr><th><tt>0011001</tt></th><th><tt>0111001</tt></th><th><tt>0111001</tt></th><th><tt>0111001</tt></th><th><tt>01110000</tt></th></tr>
	<tr><th><tt>0001100</tt></th><th><tt>1011100</tt></th><th><tt>1011100</tt></th><th><tt>1011100</tt></th><th><tt>10111000</tt></th></tr>
	<tr><th><tt>0000110</tt></th><th><tt>0101110</tt></th><th><tt>0101110</tt></th><th><tt>0101110</tt></th><th><tt>01011100</tt></th></tr>
	<tr><th><tt>0000011</tt></th><th><tt>0010111</tt></th><th><tt>0010111</tt></th><th><tt>0010111</tt></th><th><tt>00101110</tt></th></tr>
	<tr><th><tt>0000001</tt></th><th><tt>1001011</tt></th><th><tt>1001011</tt></th><th><tt>1001011</tt></th><th><tt>10010111</tt></th></tr>
</table>
