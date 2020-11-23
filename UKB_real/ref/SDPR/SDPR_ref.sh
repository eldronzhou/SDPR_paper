
for i in {1..22}; do
SDPR -make_ref -ref_prefix ../1000G/eur_SNPmaf5_nomhc -chr ${i} -ref_dir ./
done
