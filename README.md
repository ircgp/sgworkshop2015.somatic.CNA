# [Statistical Genetics Workshop 2015](http://www.kogo.or.kr/webapp/event/2015/sgworkshop/1/)
## Session 6. Analysis for CNA and CNV Data
#### Somatic copy number alteration using EXCAVATOR

```bash
git clone https://github.com/ircgp/sgworkshop2015.somatic.CNA.git

cd sgworkshop2015.somatic.CNA

./DO.00.download.EXCAVATOR

env | grep PATH

source ~/.bash_profile

env | grep PATH

./DO.20.ln.WES.sh

./DO.30.ln.hg19.fasta,SureSelect.bed

./DO.40.add.4th.col.in.target.bed.R

```
It takes 4.5 hours to execute the following command.
```
./DO.70.TargetPerla.sh
```
Just remark it and modify the code to copy pre-executed result.
```
vi DO.70.TargetPerla.sh takes

#!/bin/sh
# TargetPerla.pl SourceTarget.txt SS50M.bed SS50M

EXCAVATOR=/DATA/CUK/Somatic_Copy_Number_Alteration/EXCAVATOR_Package_v2.2
rm -rf EXCAVATOR_Package_v2.2
cp -R ${EXCAVATOR} .
```

```
./DO.70.TargetPerla.sh
```

And run the main code.
```
./DO.80.ReadPerla.sh
```
