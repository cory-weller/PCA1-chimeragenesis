# README
This repository includes the code and files used for PCA1 chimeragenesis, including
* codon shuffling of PW5 PCA1 to achieve various degrees of homology to BY PCA1
* removal of homopolymers in the codon-shuffled sequences
* generating mutagenic repair template sequences





## Initialize codon usage table
Yeast codon usage table retrieved from [kazusa.or.jp](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N) and the body of the table was saved as  `codonTable.txt`. 

The contents of `codonTable.txt` was parsed and processed by [`formatCodons.py`](https://github.com/cory-weller/Xmera/blob/b33db0d/bin/formatCodons.py) into `seqs/codons.txt`:
```
( if [ ! -f "codons.txt" ]; then python3 ../Xmera/bin/formatCodons.py > "codons.txt"; fi )
```





## Acquire BY and PW5 PCA1 sequences
BY PCA1 retrieved from [SGD](https://www.yeastgenome.org/locus/S000000499) and renamed from `S288C_YBR295W_PCA1_coding.fsa` to `seqs/BY_PCA1.cds.fasta`.

PW5 PCA1 was confirmed by [Sanger sequencing](https://benchling.com/s/seq-0foOnEFeYKrPPwgw5EJ6), and saved as `seqs/PW5_PCA1.cds.fasta`. Click `alignments` on right sidebar in benchling, then click on the saved alignment to view the chromatogram.






## Generate aligned sequences
```bash
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
```bash
( cd seqs/ && Rscript ../Xmera/bin/shuffle.R BY_PCA1.align.dna PW5_PCA1.align.dna )
```
The script [`shuffle.R`](https://github.com/cory-weller/Xmera/blob/b33db0d/bin/shuffle.R) generates five new `fasta` files with various % identity shared with `PCA1.cds.fasta`:
* seqs/PW5_PCA1.align.high.fasta
* seqs/PW5_PCA1.align.low.fasta
* seqs/PW5_PCA1.align.max.fasta
* seqs/PW5_PCA1.align.medium.fasta
* seqs/PW5_PCA1.align.min.fasta


 ## Remove homopolymers

The [`homopolymers.py`](https://github.com/cory-weller/Xmera/blob/b33db0d/bin/homopolymers.py) script removes homopolymers. See `class fasta` within the script for exact replacements.

```bash
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

```bash
python3 manuallyEditRepeats.py seqs/PW5_PCA1.min_homology.fasta > seqs/PW5_PCA1.min.cds.fasta
python3 manuallyEditRepeats.py seqs/PW5_PCA1.low_homology.fasta > seqs/PW5_PCA1.low.cds.fasta
python3 manuallyEditRepeats.py seqs/PW5_PCA1.medium_homology.fasta > seqs/PW5_PCA1.medium.cds.fasta
python3 manuallyEditRepeats.py seqs/PW5_PCA1.high_homology.fasta > seqs/PW5_PCA1.high.cds.fasta
python3 manuallyEditRepeats.py seqs/PW5_PCA1.max_homology.fasta > seqs/PW5_PCA1.max.cds.fasta
cat seqs/PW5_PCA1.min.cds.fasta \
    seqs/PW5_PCA1.low.cds.fasta \
    seqs/PW5_PCA1.medium.cds.fasta \
    seqs/PW5_PCA1.high.cds.fasta \
    seqs/PW5_PCA1.max.cds.fasta > seqs/PW5_PCA1.gblocks.fasta
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
Add nucleotide sequences upstream and downstream of the original genome location used for chimeragenesis.  i.e. when chimerizing BY and PW5 PCA1, it will be done around the original PCA1 site. This will require two files.

I took the `Genomic DNA +/- 1kb` fasta from [SGD](https://www.yeastgenome.org/locus/S000000499#sequence), saved as [`S288C_YBR295W_PCA1_flanking.fsa`](seqs/S288C_YBR295W_PCA1_flanking.fsa)

I separated it into two files. The first 1000 nucleotides made [`PCA1.upstream.fasta`](seqs/BY_PCA1.upstream.fasta) and last 1000 nucleotides made [`PCA1.downstream`](seqs/BY_PCA1.downstream.fasta):

```bash
echo ">BY_PCA1_upstream" > seqs/BY_PCA1.upstream.fasta
tail -n +2 seqs/S288C_YBR295W_PCA1_flanking.fsa | tr -d "\n" | head -c 1000 | fold >> seqs/BY_PCA1.upstream.fasta

echo ">BY_PCA1_downstream" > seqs/BY_PCA1.downstream.fasta
tail -n +2 seqs/S288C_YBR295W_PCA1_flanking.fsa | tr -d "\n" | tail -c 1000 | fold >> seqs/BY_PCA1.downstream.fasta
```

## Clean up unneeded files
_Note: these files are still available in previous commits._
```bash
rm seqs/align.in.pep
rm seqs/BY_PCA1.align.pep
rm seqs/PW5_PCA1.align.pep
rm seqs/PW5_PCA1.align.high.fasta
rm seqs/PW5_PCA1.align.low.fasta
rm seqs/PW5_PCA1.align.max.fasta
rm seqs/PW5_PCA1.align.medium.fasta
rm seqs/PW5_PCA1.align.min.fasta
rm seqs/PW5_PCA1.min_homology.fasta
rm seqs/PW5_PCA1.low_homology.fasta
rm seqs/PW5_PCA1.medium_homology.fasta
rm seqs/PW5_PCA1.high_homology.fasta
rm seqs/PW5_PCA1.max_homology.fasta
rm seqs/S288C_YBR295W_PCA1_flanking.fsa
```

fileOne="BY_PCA1.cds"
fileTwo="PW5_PCA1.cds"
flanking="BY_PCA1"
mkdir tmp && cd tmp
cp ../seqs/${fileOne}.fasta .
cp ../seqs/${fileTwo}.fasta .
cp ../seqs/${flanking}.upstream.fasta .
cp ../seqs/${flanking}.downstream.fasta .
python3 ../Xmera/bin/buildRTs.py ${fileOne} ${fileTwo} --flanking BY_PCA1 --repair-template-length 160 --unique protein > test_1.fasta
python3 ../Xmera/bin/buildRTs.py ${fileTwo} ${fileOne} --flanking BY_PCA1 --repair-template-length 160 --unique protein > test_2.fasta


## Copy required files to this directory
```
cp ../codon_shuffle/CAD2.min.fasta ./CAD2.min.cds.fasta
cp ../codon_shuffle/CAD2.low.fasta ./CAD2.low.cds.fasta
cp ../codon_shuffle/CAD2.medium.fasta ./CAD2.medium.cds.fasta
cp ../codon_shuffle/CAD2.high.fasta ./CAD2.high.cds.fasta
cp ../codon_shuffle/CAD2.cds.fasta .
```

## Calculate expected number of RT for aligned transitions
In `python`:
```python
BY = """MKPEKLFSGLGTSDGEYGVVNSENISIDAMQDNRGECHRRSIEMHANDNLGLVSQRDCTNRPKITPQECLSETEQICHHGENRTKAGLDV\
----------------------------------------------------------------------------DDAETGGDHTNESRVDECCAEK\
VND------------------------------------TETGLDVDSCCGDAQTGGDHTNESCVDGCCVRD--------------------------\
---------SSVMVEEVTGSCEAVSSKEQLLTSFEVVPSKSEGLQSIHDIRETTRCNT-NSNQHTGKGRLCIESSDSTLKKRSCKVSRQKIEVSSKPE\
CCNISCVERIASRSCEKRTFKGSTNVGISGSSSTDSLSEKFFSEQYSRMYNRYSSILKNLGCICNYLRTLGKESCCLPKVRFCSGEGASKKTKYSYRN\
SSGCLTKKKTHGDKERLSNDNGHADFVCSKSCCTKMKDCAVTSTISGHSSSEISRIVSMEPIE--NHLNLEAGSTGTEHIVLSVSGMSCTGCESKLKK\
SFGALKCVHGLKTSLILSQAEFNLDLAQGSVKDVIKHLSKTTEFKYEQISNHGSTIDVVVPYAAKDFINEEWPQGVTELKIVERNIIRIYFDPKVIGA\
RDLVNEGWSVPVSIAPFSCHPTIEVGRKHLVRVGCTTALSIILTIPILVMAWAPQLREKISTISASMVLATIIQFVIAGPFYLNALKSLIFSRLIEMD\
LLIVLSTSAAYIFSIVSFGYFVVGRPLSTEQFFETSSLLVTLIMVGRFVSELARHRAVKSISVRSLQASSAILVDKTGKETEINIRLLQYGDIFKVLP\
DSRIPTDGTVISGSSEVDEALITGESMPVPKKCQSIVVAGSVNGTGTLFVKLSKLPGNNTISTIATMVDEAKLTKPKIQNIADKIASYFVPTIIGITV\
VTFCVWIAVGIRVEKQSRSDAVIQAIIYAITVLIVSCPCVIGLAVPIVFVIASGVAAKRGVIFKSAESIEVAHNTSHVVFDKTGTLTEGKLTVVHETV\
RGDRHNSQSLLLGLTEGIKHPVSMAIASYLKEKGVSAQNVSNTKAVTGKRVEGTSYSGLKLQGGNCRWLGHNNDPDVRKALEQGYSVFCFSVNGSVTA\
VYALEDSLRADAVSTINLLRQRGISLHILSGDDDGAVRSMAARLGIESSNIRSHATPAEKSEYIKDIVEGRNCDSSSQSKRPVVVFCGDGTNDAIGLT\
QATIGVHINEGSEVAKLAADVVMLKPKLNNILTMITVSQKAMFRVKLNFLWSFTYNLFAILLAAGAFVDFHIPPEYAGLGELVSILPVIFVAILLRYA\
KI*"""

PW5 = """MKPKKFDFSLGTSDNDRKAGNSENIPIDTMYGSN-----WPVEMHASDDQRLSS----------RKQKCLNKSEVMCNEGGNRNKAASD\
SDSCCGDAQRGENYSDESCVDKCCAEKENETEAAFDSDSCCGDAQRGENYSDESCVNECCAKKENETEAASDSDSCCGDAQRGENYSDESCVDKCCAE\
KENETEAAFDSDSCCGDAQRGENYSDESCVNECCAKKENETEAASDSDSCCGDAQRGENYSDESCVNECCAKKENETEAASGSDSCCGDAQKDSKFPE\
KYADKCSIESSSVMIEEIADNYEAECCKGQLLPGVKVVSGECGGEQPTCGVREMPHCEPGSSNQQTGRGDFCFESRNSILKKRGFRVGRKNIEVSGKA\
ECCNISCVERLASRH-EKKMFDANANVGVSSSCSSDDLSGKSFSEHYSETYNRYSSILKNLGCICTYLRSLGKKSCCLPKIRFCSGEDTSIKKKYSHR\
NSSGRLTTKRAQRDGKKLSND--TADFACSKSCCRKIMNRAVSSAIYERSSNEIPRSVPIEPIREIDHLNLEAGSTGNEHVVLSVSGMTCTGCETKLK\
RSFASLKYVHNLKTSLILSQAEFDLDLAQASVKDIIRHLSKTTEFKYEQILDHGSTIDVVVPYAAKDFINEEWPQGVTELKIIEKNIVRIYFDAKIIG\
ARDLVNKGWNMPVKLAPPSAHPTVEVGRKHLVRVGCTTAISIMLTIPILVMAWAPHLREKVSTMSASMGLATIIQVLIAGPFYSNALKSLIFSRLIEM\
DLLIVLSTSAAYIFSIVSFGYFVAGRPLSTEQFFETSSLLVTLIMVGRFVSELARHRAVKSISVRSLQASSAILVDETGNETEIDIRLLQYGDTFKVL\
PDSRIPTDGTVISGSSEVDEALITGESMPVPKKCQSIVIAGSVNGTSTLFVKLIKLPGNNTISTIATMVDEAKLTKPKIQNIADKIASYFVPTIIGIT\
VITFCVWIGVGISVKKQSRSDAVIQAIIYAITVLIVSCPCAIGLAVPMVFVIASGVAAKRGVIFKSAESIEVAHNTSHVVFDKTGTLTEGKLTVVHEI\
IRDDRLNSRSLLLGLTEGVKHPVSIAIASYLKEQSVSAENVFNTKAVTGKGVEGTSQSGLKLQGGNCRWLSYSNDVDVRKALDQGYSVFCFSVNGSLT\
AVYALEDSVRADAASTINLLRQRGISLHILSGDDDGAVRSLAVRLGIERSNVRSHATPAEKGEYIKDIVEGKNFDN-SQSKRPVVVFCGDGTNDAVAL\
TQATIGVHINEGSEIAKLAADVVMLKPKLNNIITMITVSRKAMFRVKLNFIWSFTYNLFAILLAAGVFVDFHIPPEYAGLGELVSILPVIFVATLLRC\
ASI*"""

# begin with n = 2 for pure BY and pure PW5 alleles
n = 2

# iterate through aligned amino acids
for position in zip(BY, PW5):
    BY_aa = position[0]
    PW5_aa = position[1]
    # if not a gap and not identical:
    if BY_aa not in ["-", PW5_aa] and PW5_aa not in ["-", BY_aa]:
        n += 1

print(str(n))
# 261
```





## Generate Repair Templates that show the method works
This experiment uses...
* the longest repair template homology arms (80 bp for each arm, or 160 bp total)
* the most diverged PW5 PCA1 allele (lowest homology to BY PCA1)
* `--mode aligned` and `--unique all` to transition at all aligned codons

```
mkdir -p 01_test_method
python3 Xmera/bin/buildRTs.py \
    seqs/BY_PCA1.cds \
    seqs/PW5_PCA1.min.cds \
    --flanking seqs/BY_PCA1 \
    --repair-template-length 160 \
    --mode aligned \
    --unique protein > 01_test_method/BY-PW5.min.RT-160.fasta

python3 Xmera/bin/buildRTs.py \
    seqs/PW5_PCA1.min.cds \
    seqs/BY_PCA1.cds \
    --flanking seqs/BY_PCA1 \
    --repair-template-length 160 \
    --mode aligned \
    --unique protein > 01_test_method/PW5.min-BY.RT-160.fasta
# 262 rts each * 2 = 524
```






## Look at variety of repair template lengths and sequence homology
```
mkdir -p 02_RT_length_homology
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode aligned --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/BY_PCA1_PW5_PCA1.min.RT-all.fasta  ::: seqs/BY_PCA1.cds ::: seqs/PW5_PCA1.min.cds ::: 40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/BY_PCA1_PW5_PCA1.low.RT-strict.fasta  ::: seqs/BY_PCA1.cds ::: seqs/PW5_PCA1.low.cds ::: 40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/BY_PCA1_PW5_PCA1.medium.RT-strict.fasta  ::: seqs/BY_PCA1.cds ::: seqs/PW5_PCA1.medium.cds ::: 40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/BY_PCA1_PW5_PCA1.high.RT-strict.fasta  ::: seqs/BY_PCA1.cds ::: seqs/PW5_PCA1.high.cds ::: 40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/BY_PCA1_PW5_PCA1.max.RT-strict.fasta  ::: seqs/BY_PCA1.cds ::: seqs/PW5_PCA1.max.cds ::: 40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/BY_PCA1_PW5_PCA1.wt.RT-strict.fasta  ::: seqs/BY_PCA1.cds ::: seqs/PW5_PCA1.cds ::: 40 50 60 70 80 90 100 120 140 160

parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode aligned --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/PW5_PCA1.min_BY_PCA1.RT-all.fasta  ::: seqs/PW5_PCA1.min.cds ::: seqs/BY_PCA1.cds :::  40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/PW5_PCA1.low_BY_PCA1.RT-strict.fasta  ::: seqs/PW5_PCA1.low.cds ::: seqs/BY_PCA1.cds :::  40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/PW5_PCA1.medium_BY_PCA1.RT-strict.fasta  ::: seqs/PW5_PCA1.medium.cds ::: seqs/BY_PCA1.cds :::  40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/PW5_PCA1.high_BY_PCA1.RT-strict.fasta  ::: seqs/PW5_PCA1.high.cds ::: seqs/BY_PCA1.cds :::  40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/PW5_PCA1.max_BY_PCA1.RT-strict.fasta  ::: seqs/PW5_PCA1.max.cds ::: seqs/BY_PCA1.cds :::  40 50 60 70 80 90 100 120 140 160
parallel -j 1 python3 Xmera/bin/buildRTs.py {1} {2} --flanking seqs/BY_PCA1 --repair-template-length {3} --mode strict --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length_homology/PW5_PCA1.wt_BY_PCA1.RT-strict.fasta  ::: seqs/PW5_PCA1.cds ::: seqs/BY_PCA1.cds :::  40 50 60 70 80 90 100 120 140 160

```






## Look at the affects of chimerizing at all codons (including synonymous sites)
```
mkdir -p 03_all_codons
python3 Xmera/bin/buildRTs.py \
    seqs/BY_PCA1.cds \
    seqs/PW5_PCA1.min.cds \
    --flanking seqs/BY_PCA1 \
    --repair-template-length 160 \
    --mode aligned \
    --unique all > 03_all_codons/BY-PW5.min.RT-160.allCodons.fasta

python3 Xmera/bin/buildRTs.py \
    seqs/PW5_PCA1.min.cds \
    seqs/BY_PCA1.cds \
    --flanking seqs/BY_PCA1 \
    --repair-template-length 160 \
    --mode aligned \
    --unique all > 03_all_codons/PW5.min-BY.RT-160.allCodons.fasta
# 1198 RTs each * 2 = 2396
```





## Pick up here
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

Retrieved `skpp15` primers. PCR primer pairs (15-mers) obtained directly from Sri Kosuri (@skosuri). See https://github.com/KosuriLab/DropSynth

```bash
Xmera/bin/addPrimers.py 01_test_method/BY-PW5.min.RT-160.fasta 101 > PCA1.RT.txt
Xmera/bin/addPrimers.py 01_test_method/PW5.min-BY.RT-160.fasta 102 >> PCA1.RT.txt

Xmera/bin/addPrimers.py 02_RT_length_homology/PW5_PCA1.max_BY_PCA1.RT-strict.fasta 103 >> PCA1.RT.txt
Xmera/bin/addPrimers.py 02_RT_length_homology/BY_PCA1_PW5_PCA1.max.RT-strict.fasta 104 >> PCA1.RT.txt


Xmera/bin/addPrimers.py 02_RT_length_homology/BY_PCA1_PW5_PCA1.low.RT-strict.fasta 105 >> PCA1.RT.txt
Xmera/bin/addPrimers.py 02_RT_length_homology/PW5_PCA1.low_BY_PCA1.RT-strict.fasta 106 >> PCA1.RT.txt

Xmera/bin/addPrimers.py 02_RT_length_homology/BY_PCA1_PW5_PCA1.medium.RT-strict.fasta 107 >> PCA1.RT.txt
Xmera/bin/addPrimers.py 02_RT_length_homology/PW5_PCA1.medium_BY_PCA1.RT-strict.fasta 108 >> PCA1.RT.txt

Xmera/bin/addPrimers.py 02_RT_length_homology/BY_PCA1_PW5_PCA1.min.RT-all.fasta 109 >> PCA1.RT.txt
Xmera/bin/addPrimers.py 02_RT_length_homology/PW5_PCA1.min_BY_PCA1.RT-all.fasta 110 >> PCA1.RT.txt

Xmera/bin/addPrimers.py 02_RT_length_homology/BY_PCA1_PW5_PCA1.high.RT-strict.fasta 111 >> PCA1.RT.txt
Xmera/bin/addPrimers.py 02_RT_length_homology/PW5_PCA1.high_BY_PCA1.RT-strict.fasta 112 >> PCA1.RT.txt

Xmera/bin/addPrimers.py 02_RT_length_homology/PW5_PCA1.wt_BY_PCA1.RT-strict.fasta 113 >> PCA1.RT.txt
Xmera/bin/addPrimers.py 02_RT_length_homology/BY_PCA1_PW5_PCA1.wt.RT-strict.fasta 114 >> PCA1.RT.txt

Xmera/bin/addPrimers.py 03_all_codons/BY-PW5.min.RT-160.allCodons.fasta 115 >> PCA1.RT.txt
Xmera/bin/addPrimers.py 03_all_codons/PW5.min-BY.RT-160.allCodons.fasta 116 >> PCA1.RT.txt
```

## Add control oligos

```bash
Xmera/bin/printControlOligos.py chimeragenesis 160 100 > negControl.RT.tmp.txt
Xmera/bin/addPrimers.py negControl.RT.tmp.txt 117 > negControl.RT.txt && \
rm negControl.RT.tmp.txt
```