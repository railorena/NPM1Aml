cut -d':' -f1 samples | tr '\n' ' ' > samples_name.tsv
sed -i -e '$a\' samples_name.tsv
sed -i '1s/^/tag /' samples_name.tsv
sed -i -e '1 e cat samples_name.tsv' sorted_matrix.tsv
rm samples_name.tsv
gzip sorted_matrix.tsv

