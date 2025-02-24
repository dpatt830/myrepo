mkdir data
cd data 
mkdir donor1-2dpi
mkdir donor1-6dpi
mkdir donor3-2dpi
mkdir donor3-6dpi

wget -P ./donor1-2dpi/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
wget -P ./donor1-6dpi/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
wget -P ./donor3-2dpi/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
wget -P ./donor3-6dpi/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045

fasterq-dump ./donor1-2dpi/SRR5660030 -O ./donor1-2dpi/
fasterq-dump ./donor1-6dpi/SRR5660033 -O ./donor1-6dpi/
fasterq-dump ./donor3-2dpi/SRR5660044 -O ./donor3-2dpi/
fasterq-dump ./donor3-6dpi/SRR5660045 -O ./donor3-6dpi/
