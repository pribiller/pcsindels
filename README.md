# An Indel Model To Estimate Evolutionary Times For Perfectly Conserved Sequences

## Documentation

> [!NOTE]
> For detailed instructions on how to use and modify the scripts, please refer to:<br>
> [https://pribiller.github.io/pcsindels/index.html](https://pribiller.github.io/pcsindels/index.html)

## Input

The 40-vertebrate dataset files are available for download on the [UCSC website](https://hgdownload.soe.ucsc.edu/downloads.html). A detailed table with direct download links is provided below. The total dataset, including chain and fasta files for each species, requires **45.85 GB** of available space.

| Common name | UCSC name | Div. Time | Chain filename | Chain file size | Fasta filename | Fasta file size |
|:--- | :---: | ---: | :---: | ---: | :---: | ---: |
| Bonobo | panPan3 | 12.1 mya | [panPan3.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/panPan3/vsHg38/panPan3.hg38.all.chain.gz) | 0.11 GB | [panPan3.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/panPan3/bigZips/panPan3.fa.gz) | 0.89 GB |
| Chimp | panTro6 | 12.1 mya | [panTro6.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/panTro6/vsHg38/panTro6.hg38.all.chain.gz) | 0.13 GB | [panTro6.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/panTro6/bigZips/panTro6.fa.gz) | 0.90 GB |
| Gorilla | gorGor6 | 15.1 mya | [gorGor6.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/vsHg38/gorGor6.hg38.all.chain.gz) | 0.10 GB | [gorGor6.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/bigZips/gorGor6.fa.gz) | 0.88 GB |
| Orangutan | ponAbe3 | 15.2 mya | [ponAbe3.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/ponAbe3/vsHg38/ponAbe3.hg38.all.chain.gz) | 0.07 GB | [ponAbe3.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/ponAbe3/bigZips/ponAbe3.fa.gz) | 0.91 GB |
| Baboon | papAnu4 | 28.8 mya | [papAnu4.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/papAnu4/vsHg38/papAnu4.hg38.all.chain.gz) | 0.96 GB | [papAnu4.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/papAnu4/bigZips/papAnu4.fa.gz) | 0.88 GB |
| Crab-eating macaque | macFas5 | 28.8 mya | [macFas5.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/macFas5/vsHg38/macFas5.hg38.all.chain.gz) | 0.47 GB | [macFas5.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/macFas5/bigZips/macFas5.fa.gz) | 0.85 GB |
| Golden snub-nosed monkey | rhiRox1 | 28.8 mya | [rhiRox1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/rhiRox1/vsHg38/rhiRox1.hg38.all.chain.gz) | 0.57 GB | [rhiRox1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/rhiRox1/bigZips/rhiRox1.fa.gz) | 0.87 GB |
| Green monkey | chlSab2 | 28.8 mya | [chlSab2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/chlSab2/vsHg38/chlSab2.hg38.all.chain.gz) | 0.06 GB | [chlSab2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/chlSab2/bigZips/chlSab2.fa.gz) | 0.84 GB |
| Proboscis Monkey | nasLar1 | 28.8 mya | [nasLar1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/nasLar1/vsHg38/nasLar1.hg38.all.chain.gz) | 0.36 GB | [nasLar1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/nasLar1/bigZips/nasLar1.fa.gz) | 0.75 GB |
| Rhesus | rheMac10 | 28.8 mya | [rheMac10.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsHg38/rheMac10.hg38.all.chain.gz) | 0.08 GB | [rheMac10.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz) | 0.88 GB |
| Marmoset | calJac4 | 43.0 mya | [calJac4.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/calJac4/vsHg38/calJac4.hg38.all.chain.gz) | 0.97 GB | [calJac4.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/calJac4/bigZips/calJac4.fa.gz) | 0.87 GB |
| Tarsier | tarSyr2 | 69.0 mya | [tarSyr2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/tarSyr2/vsHg38/tarSyr2.hg38.all.chain.gz) | 1.15 GB | [tarSyr2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/tarSyr2/bigZips/tarSyr2.fa.gz) | 1.04 GB |
| Mouse lemur | micMur2 | 74.0 mya | [micMur2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/micMur2/vsHg38/micMur2.hg38.all.chain.gz) | 0.44 GB | [micMur2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/micMur2/bigZips/micMur2.fa.gz) | 0.73 GB |
| Malayan flying lemur | galVar1 | 79.0 mya | [galVar1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/galVar1/vsHg38/galVar1.hg38.all.chain.gz) | 1.15 GB | [galVar1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/galVar1/bigZips/galVar1.fa.gz) | 0.87 GB |
| Mouse | mm39 | 87.0 mya | [mm39.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsHg38/mm39.hg38.all.chain.gz) | 0.20 GB | [mm39.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz) | 0.81 GB |
| Rabbit | oryCun2 | 87.0 mya | [oryCun2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/oryCun2/vsHg38/oryCun2.hg38.all.chain.gz) | 0.48 GB | [oryCun2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/oryCun2/bigZips/oryCun2.fa.gz) | 0.80 GB |
| Rat | rn7 | 87.0 mya | [rn7.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/rn7/vsHg38/rn7.hg38.all.chain.gz) | 0.63 GB | [rn7.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.fa.gz) | 0.80 GB |
| Alpaca | vicPac2 | 94.0 mya | [vicPac2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/vicPac2/vsHg38/vicPac2.hg38.all.chain.gz) | 0.39 GB | [vicPac2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/vicPac2/bigZips/vicPac2.fa.gz) | 0.65 GB |
| Bison | bisBis1 | 94.0 mya | [bisBis1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/bisBis1/vsHg38/bisBis1.hg38.all.chain.gz) | 0.62 GB | [bisBis1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/bisBis1/bigZips/bisBis1.fa.gz) | 0.85 GB |
| Cat | felCat9 | 94.0 mya | [felCat9.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/felCat9/vsHg38/felCat9.hg38.all.chain.gz) | 0.67 GB | [felCat9.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/felCat9/bigZips/felCat9.fa.gz) | 0.76 GB |
| Chinese pangolin | manPen1 | 94.0 mya | [manPen1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/manPen1/vsHg38/manPen1.hg38.all.chain.gz) | 0.44 GB | [manPen1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/manPen1/bigZips/manPen1.fa.gz) | 0.62 GB |
| Cow | bosTau9 | 94.0 mya | [bosTau9.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/vsHg38/bosTau9.hg38.all.chain.gz) | 0.36 GB | [bosTau9.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.fa.gz) | 0.82 GB |
| Dog | canFam6 | 94.0 mya | [canFam6.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/vsHg38/canFam6.hg38.all.chain.gz) | 0.45 GB | [canFam6.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/bigZips/canFam6.fa.gz) | 0.71 GB |
| Ferret | musFur1 | 94.0 mya | [musFur1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/musFur1/vsHg38/musFur1.hg38.all.chain.gz) | 0.43 GB | [musFur1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/musFur1/bigZips/musFur1.fa.gz) | 0.71 GB |
| Hawaiian monk seal | neoSch1 | 94.0 mya | [neoSch1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/neoSch1/vsHg38/neoSch1.hg38.all.chain.gz) | 0.47 GB | [neoSch1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/neoSch1/bigZips/neoSch1.fa.gz) | 0.72 GB |
| Horse | equCab3 | 94.0 mya | [equCab3.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/equCab3/vsHg38/equCab3.hg38.all.chain.gz) | 0.40 GB | [equCab3.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/equCab3/bigZips/equCab3.fa.gz) | 0.77 GB |
| Little brown bat | myoLuc2 | 94.0 mya | [myoLuc2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/myoLuc2/vsHg38/myoLuc2.hg38.all.chain.gz) | 0.52 GB | [myoLuc2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/myoLuc2/bigZips/myoLuc2.fa.gz) | 0.61 GB |
| Pig | susScr11 | 94.0 mya | [susScr11.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/susScr11/vsHg38/susScr11.hg38.all.chain.gz) | 0.45 GB | [susScr11.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/susScr11/bigZips/susScr11.fa.gz) | 0.76 GB |
| Southern sea otter | enhLutNer1 | 94.0 mya | [enhLutNer1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/enhLutNer1/vsHg38/enhLutNer1.hg38.all.chain.gz) | 0.45 GB | [enhLutNer1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/enhLutNer1/bigZips/enhLutNer1.fa.gz) | 0.74 GB |
| Manatee | triMan1 | 99.0 mya | [triMan1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/triMan1/vsHg38/triMan1.hg38.all.chain.gz) | 0.50 GB | [triMan1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/triMan1/bigZips/triMan1.fa.gz) | 0.87 GB |
| Wallaby | macEug2 | 160.0 mya | [macEug2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/macEug2/vsHg38/macEug2.hg38.all.chain.gz) | 0.52 GB | [macEug2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/macEug2/bigZips/macEug2.fa.gz) | 0.81 GB |
| Platypus | ornAna2 | 180.0 mya | [ornAna2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/ornAna2/vsHg38/ornAna2.hg38.all.chain.gz) | 0.34 GB | [ornAna2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/ornAna2/bigZips/ornAna2.fa.gz) | 0.56 GB |
| Brown kiwi | aptMan1 | 319.0 mya | [aptMan1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/aptMan1/vsHg38/aptMan1.hg38.all.chain.gz) | 0.05 GB | [aptMan1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/aptMan1/bigZips/aptMan1.fa.gz) | 0.40 GB |
| Chicken | galGal6 | 319.0 mya | [galGal6.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/vsHg38/galGal6.hg38.all.chain.gz) | 0.05 GB | [galGal6.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/bigZips/galGal6.fa.gz) | 0.32 GB |
| Garter snake | thaSir1 | 319.0 mya | [thaSir1.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/thaSir1/vsHg38/thaSir1.hg38.all.chain.gz) | 0.13 GB | [thaSir1.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/thaSir1/bigZips/thaSir1.fa.gz) | 0.36 GB |
| Golden eagle | aquChr2 | 319.0 mya | [aquChr2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/aquChr2/vsHg38/aquChr2.hg38.all.chain.gz) | 0.05 GB | [aquChr2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/aquChr2/bigZips/aquChr2.fa.gz) | 0.37 GB |
| Turkey | melGal5 | 319.0 mya | [melGal5.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/melGal5/vsHg38/melGal5.hg38.all.chain.gz) | 0.05 GB | [melGal5.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/melGal5/bigZips/melGal5.fa.gz) | 0.34 GB |
| African clawed frog | xenLae2 | 352.0 mya | [xenLae2.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/vsHg38/xenLae2.hg38.all.chain.gz) | 0.14 GB | [xenLae2.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/bigZips/xenLae2.fa.gz) | 0.77 GB |
| X. tropicalis | xenTro10 | 352.0 mya | [xenTro10.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/xenTro10/vsHg38/xenTro10.hg38.all.chain.gz) | 1.06 GB | [xenTro10.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/xenTro10/bigZips/xenTro10.fa.gz) | 0.42 GB |
| Zebrafish | danRer11 | 429.0 mya | [danRer11.hg38.all.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/vsHg38/danRer11.hg38.all.chain.gz) | 0.38 GB | [danRer11.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz) | 0.51 GB |
| | | | | **16.86 GB** |  | **28.99 GB** |

## Pipeline

To be added later.

### Computational Resources

The following table presents a reference for the computational resources needed for each script in the pipeline.
**Disk space** refers to the total size of the output files, while memory usage, number of cores, and run time offer a general estimate of the computational requirements for running these scripts in a cluster. Links to the scripts and their documentation are also provided.

|  | Step | Script | Memory | Nb. cores | Run time | Disk space | 
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| 0 | Process UCSC data for pipeline use | [breakFasta.py](https://github.com/pribiller/pcsindels/blob/main/scripts/utils/breakFasta.py) |  |  |  |  |
| 1 | [Extract Perfectly Conserved Sequences (PCSs)](https://pribiller.github.io/pcsindels/1_extractPCS.html)  | [1_extractPCS.py](https://github.com/pribiller/pcsindels/blob/main/scripts/1_extractPCS.py) |  |  |  |  |
| 2 | [Split genome into windows](https://pribiller.github.io/pcsindels/2_computeWindows.html) | [2_computeWindows.py](https://github.com/pribiller/pcsindels/blob/main/scripts/2_computeWindows.py) |  |  |  |  |
| 3 | [Setup for estimating evolutionary times](https://pribiller.github.io/pcsindels/3_setupEvolTimes.html) | [3_setupEvolTimes.py](https://github.com/pribiller/pcsindels/blob/main/scripts/3_setupEvolTimes.py) |  |  |  |  |
| ├─3.1 | Parameter α = 1.1 (indels) |  |  |  |  |  |
| └─3.2 | Parameter α = 10 (subs.) |  |  |  |  |  |
| 4 | [Estimate evolutionary times](https://pribiller.github.io/pcsindels/4_estimateEvolTimes.html) | [4_estimateEvolTimes.py](https://github.com/pribiller/pcsindels/blob/main/scripts/4_estimateEvolTimes.py) |  |  |  |  |
| ├─4.1 | Parameter α = 1.1 (indels) |  |  |  |  |  |
| └─4.2 | Parameter α = 10 (subs.) |  |  |  |  |  |
| 5 | Reproduce figures from the paper |  |  |  |  |  |
| ├─5.1 | PCS size distributions comparison (Figure 2) | To be added later |  |  |  |  |
| ├─5.2 | Evolutionary time estimates comparison (Figure 3) | To be added later  |  |  |  |  |
| ├─5.3 | Indel rates comparison (Figure 4) | To be added later  |  |  |  |  |
| └─5.4 | Indel rates in functional classes (Figure 5) | To be added later  |  |  |  |  |
|  |  |  |  |  |  |  |
