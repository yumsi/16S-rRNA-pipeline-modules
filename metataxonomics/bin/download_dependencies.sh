#!/bin/bash

echo "Welcome!"
echo "Downloading the pipeline dependencies..."

fileURLs=("https://surfdrive.surf.nl/files/index.php/s/tJW0N2Tx5RiAisp" "https://surfdrive.surf.nl/files/index.php/s/RFjX7iT6tbQg6EG" "https://surfdrive.surf.nl/files/index.php/s/OOcV5uW35IfhEZA" "https://surfdrive.surf.nl/files/index.php/s/ctDBSEX4ZsJzTJI")
fileNames=("SILVA-138-SSURef-Full-Seqs.qza" "Silva-v138-full-length-seq-taxonomy.qza" "SILVA-138-Full-Seqs-Classifier.qza" "qiime2-2020.2.yml")

mkdir lib/
cd lib/

for i in ${!fileURLs[@]};
do
	wget -O ${fileNames[i]} ${fileURLs[i]}"/download"
done

echo "Done!"