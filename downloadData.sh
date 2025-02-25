mkdir data

mkdir ./data/donor1-2dpi
mkdir ./data/donor1-6dpi
mkdir ./data/donor3-2dpi
mkdir ./data/donor3-6dpi

wget -P ./data/donor1-2dpi/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
wget -P ./data/donor1-6dpi/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
wget -P ./data/donor3-2dpi/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
wget -P ./data/donor3-6dpi/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045

fasterq-dump ./data/donor1-2dpi/SRR5660030 -O ./data/donor1-2dpi/
fasterq-dump ./data/donor1-6dpi/SRR5660033 -O ./data/donor1-6dpi/
fasterq-dump ./data/donor3-2dpi/SRR5660044 -O ./data/donor3-2dpi/
fasterq-dump ./data/donor3-6dpi/SRR5660045 -O ./data/donor3-6dpi/
