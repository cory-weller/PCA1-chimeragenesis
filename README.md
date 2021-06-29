# README
This repository includes the code and files used for PCA1 chimeragenesis, including
* codon shuffling of PW5 PCA1 to achieve various degrees of homology to BY PCA1
* removal of homopolymers in the codon-shuffled sequences
* generating mutagenic repair template sequences





## Initialize codon usage table
Yeast codon usage table retrieved from [kazusa.or.jp](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N) and the body of the table was saved as  `codonTable.txt`. 

The contents of `codonTable.txt` was parsed and processed by [`formatCodons.py`](bin/formatCodons.py) into `seqs/codons.txt`:
```
( if [ ! -f "codons.tab" ]; then python3 formatCodons.py > "codons.tab"; fi )
```





## Acquire BY and PW5 PCA1 sequences
BY PCA1 retrieved from [SGD](https://www.yeastgenome.org/locus/S000000499) and renamed from `S288C_YBR295W_PCA1_coding.fsa` to `seqs/BY_PCA1.cds.fasta`.

PW5 PCA1 was confirmed by [Sanger sequencing](https://benchling.com/s/seq-0foOnEFeYKrPPwgw5EJ6), and saved as `seqs/PW5_PCA1.cds.fasta`. Click `alignments` on right sidebar in benchling, then click on the saved alignment to view the chromatogram.






## Generate aligned sequences
```
echo ">PW5_PCA1_pep" > seqs/PW5_PCA1.pep
python3 Xmera/bin/translate.py seqs/PW5_PCA1.cds.fasta >> seqs/PW5_PCA1.pep

echo ">BY_PCA1_pep" > seqs/BY_PCA1.pep
python3 Xmera/bin/translate.py seqs/BY_PCA1.cds.fasta >> seqs/BY_PCA1.pep

cat seqs/BY_PCA1.pep | fold -w 60 > seqs/align.in.pep
cat seqs/PW5_PCA1.pep | fold -w 60 >> seqs/align.in.pep

Xmera/bin/clustalo -i seqs/align.in.pep > seqs/align.pep

lineNo2=$(grep -n ">" seqs/align.pep | tail -n 1 | cut -d ":" -f 1)
let lineNo1="${lineNo2}-1"

# save separate files PW5_PCA1.align.pep and BY_PCA1.align.pep
head -n ${lineNo1} seqs/align.pep > seqs/BY_PCA1.align.pep
tail -n +${lineNo2} seqs/align.pep > seqs/PW5_PCA1.align.pep

echo ">PW5_PCA1" > seqs/PW5_PCA1.align.dna
( cd seqs/ && python3 ../Xmera/bin/aligndna.py PW5_PCA1 | fold -w 60 >> PW5_PCA1.align.dna )

echo ">BY_PCA1" > seqs/BY_PCA1.align.dna
( cd seqs/ && python3 ../Xmera/bin/aligndna.py BY_PCA1 | fold -w 60 >> BY_PCA1.align.dna )
```





## Generate shuffled sequences
```
( cd seqs/ && Rscript ../Xmera/bin/shuffle.R BY_PCA1.align.dna PW5_PCA1.align.dna )
```
The script [`shuffle.R`](shuffle.R) generates five new `fasta` files with various % identity shared with `PCA1.cds.fasta`:
* seqs/PW5_PCA1.align.high.fasta
* seqs/PW5_PCA1.align.low.fasta
* seqs/PW5_PCA1.align.max.fasta
* seqs/PW5_PCA1.align.medium.fasta
* seqs/PW5_PCA1.align.min.fasta


 ## Remove homopolymers

The [`homopolymers.py`](Xmera/bin/homopolymers.py) script removes homopolymers. See `class fasta` within the script for exact replacements.

```
( cd seqs/ && \
python3 ../Xmera/bin/homopolymers.py PW5_PCA1.align.min.fasta BY_PCA1.align.dna > PW5_PCA1.min_homology.fasta
python3 ../Xmera/bin/homopolymers.py PW5_PCA1.align.low.fasta BY_PCA1.align.dna > PW5_PCA1.low_homology.fasta
python3 ../Xmera/bin/homopolymers.py PW5_PCA1.align.medium.fasta BY_PCA1.align.dna > PW5_PCA1.medium_homology.fasta
python3 ../Xmera/bin/homopolymers.py PW5_PCA1.align.high.fasta BY_PCA1.align.dna > PW5_PCA1.high_homology.fasta
python3 ../Xmera/bin/homopolymers.py PW5_PCA1.align.max.fasta BY_PCA1.align.dna > PW5_PCA1.max_homology.fasta

# this should return no matches if there are no stretches of 5+ homopolymers
for file in PW5_PCA1.*homology.fasta; do
    for nucleotide in A C T G; do
        cat ${file} | tr -d "\n" | grep -E "[${nucleotide}]{5,}"
    done
done

for file in PW5_PCA1.*homology.fasta; do
    python3 ../Xmera/bin/formatfasta.py $file | fold -w 70 >> PW5_PCA1_gene_blocks.fasta
done )
```

## Manually edit tandem repeats
The first ~800 nucleotides of PW5 PCA1 contains 5 (plus one partial) tandem repeats. 
I manually edited the region to avoid tandem homology within the shuffled sequence,
then replaced the first 789 nucleotides of ALL shuffled sequences to contain this
manually-shuffled sequence:
```
ATGAAACCTAAGAAGTTTGATTTCTCCCTTGGAACTTCAGACAATGATAGGAAGGCTGGCAACTCTGAAA
ATATTCCTATTGACACGATGTATGGCTCCAATTGGCCTGTGGAGATGCACGCCTCCGATGACCAACGCCT
ATCAAGTCGCAAACAAAAGTGCCTAAACAAATCCGAAGTTATGTGTAACGAAGGAGGTAACCGAAACAAA
GCTGCCTCAGACTCAGACAGCTGTTGCGGTGATGCCCAGCGAGGCGAAAATTACAGTGACGAGTCTTGTG
TAGACAAGTGCTGTGCTGAGAAGGAGAATGAAACGGAAGCCGCGTTTGATTCCGATAGTTGCTGTGGCGA
TGCTCAACGCGGGGAAAACTATAGTGATGAATCATGTGTCAACGAGTGTTGTGCAAAGAAAGAGAACGAG
ACTGAGGCGGCATCGGACTCGGACTCCTGTTGTGGTGACGCTCAACGTGGAGAGAATTACAGTGATGAAA
GCTGTGTAGATAAGTGTTGCGCCGAAAAGGAGAACGAAACCGAAGCAGCTTTCGACTCTGACTCGTGCTG
CGGGGATGCTCAGCGTGGCGAGAACTATTCGGACGAGTCATGTGTTAATGAATGCTGTGCGAAGAAAGAA
AATGAGACCGAGGCCGCTTCAGACAGCGATTCGTGTTGCGGAGACGCGCAGCGCGGTGAAAACTACTCAG
ATGAGTCATGCGTAAATGAGTGCTGCGCTAAGAAGGAAAACGAGACAGAAGCGGCTAGCGGTAGTGACTC
ATGCTGTGGCGATGCACAG
```

```
python3 manuallyEditRepeats.py seqs/PW5_PCA1.min_homology.fasta > seqs/PW5_PCA1.gblocks.fasta
python3 manuallyEditRepeats.py seqs/PW5_PCA1.low_homology.fasta >> seqs/PW5_PCA1.gblocks.fasta
python3 manuallyEditRepeats.py seqs/PW5_PCA1.medium_homology.fasta >> seqs/PW5_PCA1.gblocks.fasta
python3 manuallyEditRepeats.py seqs/PW5_PCA1.high_homology.fasta >> seqs/PW5_PCA1.gblocks.fasta
python3 manuallyEditRepeats.py seqs/PW5_PCA1.max_homology.fasta >> seqs/PW5_PCA1.gblocks.fasta
```

All these gblock sequences translate to the same protein sequence,
which is identical to the translated Sanger-confirmed sequence:
```
MKPKKFDFSLGTSDNDRKAGNSENIPIDTMYGSNWPVEMHASDDQRLSSRKQKCLNKSEV
MCNEGGNRNKAASDSDSCCGDAQRGENYSDESCVDKCCAEKENETEAAFDSDSCCGDAQR
GENYSDESCVNECCAKKENETEAASDSDSCCGDAQRGENYSDESCVDKCCAEKENETEAA
FDSDSCCGDAQRGENYSDESCVNECCAKKENETEAASDSDSCCGDAQRGENYSDESCVNE
CCAKKENETEAASGSDSCCGDAQKDSKFPEKYADKCSIESSSVMIEEIADNYEAECCKGQ
LLPGVKVVSGECGGEQPTCGVREMPHCEPGSSNQQTGRGDFCFESRNSILKKRGFRVGRK
NIEVSGKAECCNISCVERLASRHEKKMFDANANVGVSSSCSSDDLSGKSFSEHYSETYNR
YSSILKNLGCICTYLRSLGKKSCCLPKIRFCSGEDTSIKKKYSHRNSSGRLTTKRAQRDG
KKLSNDTADFACSKSCCRKIMNRAVSSAIYERSSNEIPRSVPIEPIREIDHLNLEAGSTG
NEHVVLSVSGMTCTGCETKLKRSFASLKYVHNLKTSLILSQAEFDLDLAQASVKDIIRHL
SKTTEFKYEQILDHGSTIDVVVPYAAKDFINEEWPQGVTELKIIEKNIVRIYFDAKIIGA
RDLVNKGWNMPVKLAPPSAHPTVEVGRKHLVRVGCTTAISIMLTIPILVMAWAPHLREKV
STMSASMGLATIIQVLIAGPFYSNALKSLIFSRLIEMDLLIVLSTSAAYIFSIVSFGYFV
AGRPLSTEQFFETSSLLVTLIMVGRFVSELARHRAVKSISVRSLQASSAILVDETGNETE
IDIRLLQYGDTFKVLPDSRIPTDGTVISGSSEVDEALITGESMPVPKKCQSIVIAGSVNG
TSTLFVKLIKLPGNNTISTIATMVDEAKLTKPKIQNIADKIASYFVPTIIGITVITFCVW
IGVGISVKKQSRSDAVIQAIIYAITVLIVSCPCAIGLAVPMVFVIASGVAAKRGVIFKSA
ESIEVAHNTSHVVFDKTGTLTEGKLTVVHEIIRDDRLNSRSLLLGLTEGVKHPVSIAIAS
YLKEQSVSAENVFNTKAVTGKGVEGTSQSGLKLQGGNCRWLSYSNDVDVRKALDQGYSVF
CFSVNGSLTAVYALEDSVRADAASTINLLRQRGISLHILSGDDDGAVRSLAVRLGIERSN
VRSHATPAEKGEYIKDIVEGKNFDNSQSKRPVVVFCGDGTNDAVALTQATIGVHINEGSE
IAKLAADVVMLKPKLNNIITMITVSRKAMFRVKLNFIWSFTYNLFAILLAAGVFVDFHIP
PEYAGLGELVSILPVIFVATLLRCASI*
```

The sequences in [`PW5_PCA1.gblocks.fasta`](seqs/PW5_PCA1.gblocks.fasta) were ordered
via gene synthesis.








## Preparing to generate repair templates
Add nucleotide sequences upstream and downstream of the original genome location used for chimeragenesis.  i.e. when chimerizing PCA1 with CAD2 allele, it will be done around the original PCA1 site. This will require two files, `PCA1.upstream` and `PCA1.downstream`. 

I took the `Genomic DNA +/- 1kb` fasta from [SGD](https://www.yeastgenome.org/locus/S000000499#sequence). The first 1000 nucleotides made `PCA1.upstream` and last 1000 nucleotides made `PCA1.downstream`.






## Copy required files to this directory
```
cp ../codon_shuffle/CAD2.min.fasta ./CAD2.min.cds.fasta
cp ../codon_shuffle/CAD2.low.fasta ./CAD2.low.cds.fasta
cp ../codon_shuffle/CAD2.medium.fasta ./CAD2.medium.cds.fasta
cp ../codon_shuffle/CAD2.high.fasta ./CAD2.high.cds.fasta
cp ../codon_shuffle/CAD2.cds.fasta .
```






## Generate Repair Templates that show the method works
Using longest repair templates (80 bp for each, or 160 bp total) and most diverged CAD2 sequence
```
mkdir -p 01_test_method
python3 ./chimera.py BY_PCA1 PW5_PCA1.min --flanking BY_PCA1 --repair-template-length 160 --unique protein > 01_test_method/BY-PW5.min.RT-160.lax.fasta
python3 ./chimera.py PW5_PCA1.min BY_PCA1 --flanking BY_PCA1 --repair-template-length 160 --unique protein > 01_test_method/PW5-BY.min.RT-160.lax.fasta

python3 ./chimera.py BY_PCA1 PW5_PCA1.min --strict --flanking BY_PCA1 --repair-template-length 160 --unique protein > 01_test_method/BY-PW5.min.RT-160.strict.fasta
python3 ./chimera.py PW5_PCA1.min BY_PCA1 --strict --flanking BY_PCA1 --repair-template-length 160 --unique protein > 01_test_method/PW5-BY.min.RT-160.strict.fasta
# 262 rts each * 2 = 524
```






## Look at variety of repair template lengths and sequence homology
This will be a total of 10 transformations (PCA1-CAD2 and CAD2-PCA1 orientations, with 5 levels of sequence homology)
```
mkdir -p 02_RT_length
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.min.RT-all.lax.fasta  ::: BY_PCA1 ::: PW5_PCA1.min :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.low.RT-all.lax.fasta  ::: BY_PCA1 ::: PW5_PCA1.low :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.medium.RT-all.lax.fasta  ::: BY_PCA1 ::: PW5_PCA1.medium :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.high.RT-all.lax.fasta  ::: BY_PCA1 ::: PW5_PCA1.high :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.orig.RT-all.lax.fasta ::: BY_PCA1 ::: PW5_PCA1.max ::: 40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.min_BY_PCA1.RT-all.lax.fasta  ::: PW5_PCA1.min ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.low_BY_PCA1.RT-all.lax.fasta  ::: PW5_PCA1.low ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.medium_BY_PCA1.RT-all.lax.fasta  ::: PW5_PCA1.medium ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.high_BY_PCA1.RT-all.lax.fasta  ::: PW5_PCA1.high ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160

parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.min.RT-all.strict.fasta  ::: BY_PCA1 ::: PW5_PCA1.min :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.low.RT-all.strict.fasta  ::: BY_PCA1 ::: PW5_PCA1.low :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.medium.RT-all.strict.fasta  ::: BY_PCA1 ::: PW5_PCA1.medium :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.high.RT-all.strict.fasta  ::: BY_PCA1 ::: PW5_PCA1.high :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.orig.RT-all.strict.fasta ::: BY_PCA1 ::: PW5_PCA1.max ::: 40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.min_BY_PCA1.RT-all.strict.fasta  ::: PW5_PCA1.min ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.low_BY_PCA1.RT-all.strict.fasta  ::: PW5_PCA1.low ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.medium_BY_PCA1.RT-all.strict.fasta  ::: PW5_PCA1.medium ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.high_BY_PCA1.RT-all.strict.fasta  ::: PW5_PCA1.high ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.orig_BY_PCA1.RT-all.strict.fasta ::: PW5_PCA1.max ::: BY_PCA1 ::: 40 44 50 58 68 80 94 110 128 148 160

```






## Look at the affects of chimerizing at synonymous sites
```
mkdir -p 03_synonymous_RT
python3 ./chimera.py  BY_PCA1 PW5_PCA1.min --flanking BY_PCA1 --repair-template-length 160 --unique dna --oligo-length 190 --repair-template-length 160 > 03_synonymous_RT/BY-PW5.min.RT-160-syn.lax.fasta
python3 ./chimera.py  PW5_PCA1.min BY_PCA1 --flanking BY_PCA1 --repair-template-length 160 --unique dna --oligo-length 190 --repair-template-length 160 > 03_synonymous_RT/PW5-BY.min.RT-160-syn.lax.fasta

python3 ./chimera.py  BY_PCA1 PW5_PCA1.min --strict --flanking BY_PCA1 --repair-template-length 160 --unique dna --oligo-length 190 --repair-template-length 160 > 03_synonymous_RT/BY-PW5.min.RT-160-syn.strict.fasta
python3 ./chimera.py  PW5_PCA1.min BY_PCA1 --strict --flanking BY_PCA1 --repair-template-length 160 --unique dna --oligo-length 190 --repair-template-length 160 > 03_synonymous_RT/PW5-BY.min.RT-160-syn.strict.fasta

# 1198 RTs each * 2 = 2396
```






## Generate chimeras between AcrIIa2 and AcrIIa2b (no codon shuffling needed)
```

```






## Generate all possible chimeras (no alignment needed) for AcrIIa2b and AcrIIa4
```
python3 chimera.py AcrIIa2b AcrIIa4 --all --unique dna --flanking acr --repair-template-length 160 --primer-length 15 --oligo-length 190 > 04_Acrs/AcrIIa2b_AcrIIa4.RT-160.fasta
python3 chimera.py AcrIIa4 AcrIIa2b --all --unique dna --flanking acr --repair-template-length 160 --primer-length 15 --oligo-length 190 > 04_Acrs/AcrIIa4_AcrIIa2b.RT-160.fasta
# 23,056
```






## Add skpp15 primers

Retrieved `skpp15` primers. PCR primer pairs (15-mers) obtained directly from Sri Kosuri (@skosuri). Modified from Elledge barcodes, https://elledge.hms.harvard.edu/?page_id=638

```
wget https://github.com/lasersonlab/pribar/raw/master/skpp15-forward.faa
wget https://github.com/lasersonlab/pribar/raw/master/skpp15-reverse.faa
```

Built reverse complement of of the reverse primers at https://www.bioinformatics.org/sms2/rev_comp.html and removed blank lines with `sed -i '/^$/d' skpp15-reverse-complemented.faa` 

Convert FASTA format to single-line csv

```
for file in $(ls 01_test_method/*.fasta 02_RT_length/*.fasta 03_synonymous_RT/*.fasta); do
    cat ${file} | tr "\n" "@" | sed 's/@>/\n/g' | sed 's/each@/each,/g' | tr -d "@" > ${file%.fasta}.csv
done
```

Concatenate skpp15 primers to RT sequences

```
N=100
for file in $(ls 01_test_method/*.csv 02_RT_length/*.csv 03_synonymous_RT/*.csv); do
    let N=N+1
    let lineNo=N*2
    forward=$(sed -n ${lineNo}p skpp15-forward.faa)
    reverse=$(sed -n ${lineNo}p skpp15-reverse-complemented.faa)
    awk -F "," -v f=${forward} -v r=${reverse} '{print $1,f$2r}' ${file} > ${file%.csv}.skpp${N}.RT.csv
done
```

Zip skpp15 RT oligos

```
find . | grep skpp | grep csv | zip skpp15.RTs.zip -@

```