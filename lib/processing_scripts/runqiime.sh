# Usage bash div.preprocess.sh

# To run on my own computer
source activate qiime

# change directory
cd ../data/fiber/grains_data/

# convert grain tables to biom
biom convert -i grains_fiber.txt -o grains_fiber.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy

# calculate beta diversity tables for grains
beta_diversity.py -i grains_fiber.biom -o grains_beta -t ../../../raw/diet.tree.txt

# calculate alpha diversity for grains
alpha_diversity.py -i grains_fiber.biom -o grains_alpha.txt -m PD_whole_tree,shannon,chao1,simpson -t ../../../raw/diet.tree.txt


cd ../fruit_data/

# convert food fruit table to biom
biom convert -i fruit_fiber.txt -o fruit_fiber.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy

# calculate beta diversity tables for fruit
beta_diversity.py -i fruit_fiber.biom -o fruit_beta -t ../../../raw/diet.tree.txt

# calculate alpha diversity for fruit
alpha_diversity.py -i fruit_fiber.biom -o fruit_alpha.txt -m PD_whole_tree,shannon,chao1,simpson -t ../../../raw/diet.tree.txt

cd ../beans_data/

# convert food beans table to biom
biom convert -i beans_fiber.txt -o beans_fiber.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy

# calculate beta diversity tables for beans
beta_diversity.py -i beans_fiber.biom -o beans_beta -t ../../../raw/diet.tree.txt

# calculate alpha diversity for beans
alpha_diversity.py -i beans_fiber.biom -o beans_alpha.txt -m PD_whole_tree,shannon,chao1,simpson -t ../../../raw/diet.tree.txt

cd ../vege_data/

# convert food vege table to biom
biom convert -i vege_fiber.txt -o vege_fiber.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy

# calculate beta diversity tables for vege
beta_diversity.py -i vege_fiber.biom -o vege_beta -t ../../../raw/diet.tree.txt

# calculate alpha diversity for vege
alpha_diversity.py -i vege_fiber.biom -o vege_alpha.txt -m PD_whole_tree,shannon,chao1,simpson -t ../../../raw/diet.tree.txt