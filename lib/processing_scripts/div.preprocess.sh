# Usage bash div.preprocess.sh

# To run on my own computer
# source activate qiime

# change directory
cd ../data/processed_food/

# convert food otu tables to biom
biom convert -i food.txt -o food.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
biom convert -i dhydrt.txt -o dhydrt.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
biom convert -i fiber.txt -o fiber.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
biom convert -i fiber.norm.counts.txt -o fiber.norm.counts.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
biom convert -i fiber.norm.counts.no.soy.txt -o fiber.norm.counts.no.soy.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
biom convert -i dhydrt.smry.txt -o dhydrt.smry.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
biom convert -i dhydrt.smry.no.soy.txt -o dhydrt.smry.no.soy.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
biom convert -i dhydrt.norm.counts.txt -o dhydrt.norm.counts.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
biom convert -i dhydrt.pre.txt -o dhydrt.pre.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy
biom convert -i dhydrt.post.txt -o dhydrt.post.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy

# calculate beta diversity tables for food otus
beta_diversity.py -i food.biom -o food_beta -t ../../raw/diet.tree.txt
beta_diversity.py -i dhydrt.biom -o dhydrt_beta -t ../../raw/diet.tree.txt
beta_diversity.py -i fiber.biom -o fiber_beta -t ../../raw/diet.tree.txt
beta_diversity.py -i fiber.norm.counts.biom -o fiber_norm_counts_beta -t ../../raw/diet.tree.txt
beta_diversity.py -i fiber.norm.counts.no.soy.biom -o fiber_norm_counts_no_soy_beta -t ../../raw/diet.tree.txt
beta_diversity.py -i dhydrt.smry.biom -o dhydrt_smry_beta -t ../../raw/diet.tree.txt
beta_diversity.py -i dhydrt.smry.no.soy.biom -o dhydrt_smry_no_soy_beta -t ../../raw/diet.tree.txt
beta_diversity.py -i dhydrt.norm.counts.biom -o dhydrt_norm_counts_beta -t ../../raw/diet.tree.txt
beta_diversity.py -i dhydrt.pre.biom -o dhydrt_pre_beta -t ../../raw/diet.tree.txt
beta_diversity.py -i dhydrt.post.biom -o dhydrt_post_beta -t ../../raw/diet.tree.txt

# calculate alpha diversity tables for food otus
alpha_diversity.py -i dhydrt.biom -o dhydrt_alpha.txt -m PD_whole_tree,shannon,chao1,simpson,observed_otus -t ../../raw/diet.tree.txt
alpha_diversity.py -i food.biom -o food_alpha.txt -m PD_whole_tree,shannon,chao1,simpson,observed_otus -t ../../raw/diet.tree.txt
alpha_diversity.py -i fiber.biom -o fiber_alpha.txt -m PD_whole_tree,shannon,chao1,simpson,observed_otus -t ../../raw/diet.tree.txt
alpha_diversity.py -i dhydrt.pre.biom -o dhydrt_pre_alpha.txt -m PD_whole_tree,shannon,chao1,simpson,observed_otus -t ../../raw/diet.tree.txt
alpha_diversity.py -i dhydrt.post.biom -o dhydrt_post_alpha.txt -m PD_whole_tree,shannon,chao1,simpson,observed_otus -t ../../raw/diet.tree.txt

# change directory
cd ../processed_tax/

# convert taxonomy table to biom
biom convert -i taxonomy_norm_s.txt -o taxonomy_norm_s.biom --to-json --table-type "OTU table"
biom convert -i taxonomy_counts_s.txt -o taxonomy_counts_s.biom --to-json --table-type "OTU table"

# convert strain table to biom
biom convert -i ../../raw/strain_burst_omegaGWG.txt -o strain.biom --to-json --table-type "OTU table"
biom convert -i ../../raw/strain_burst_omegaGWG_norm.txt -o strain_norm.biom --to-json --table-type "OTU table"

# calculate beta diversity tables for taxonomy
beta_diversity.py -i taxonomy_norm_s.biom -m chisq,bray_curtis -o taxa_beta

# calculate alpha diversity tables for taxonomy
alpha_diversity.py -i taxonomy_counts.biom -o taxa_alpha.txt -m shannon,chao1,simpson

# calculate beta diversity tables for strain level from burst
beta_diversity.py -i strain.biom -o strain_beta -t ../../raw/AllNewerTre.tre
beta_diversity.py -i strain_norm.biom -o strain_beta_norm -t ../../raw/AllNewerTre.tre