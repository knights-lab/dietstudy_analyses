#!/bin/bash -x

# Usage: procrustes.sh

# first activate qiime conda environment 
source activate qiime

#navigate to correct dir
cd /Users/abby/Documents/Projects/dietstudy/figures/figure\ 4/procrustes/data_username_5day

# make qiime/tree friendly biom table from food
for i in *_food.txt;
do biom convert -i $i -o ${i//_food.txt}_food.biom --to-json --table-type 'OTU table' --process-obs-metadata=taxonomy; 
done

# make qiime friendly biom table from otu table
for i in *_tax.txt;
do biom convert -i $i -o ${i//_tax.txt}_tax.biom --to-json --table-type "OTU table";
done


# make beta diversity tables (food)
for i in *_food.biom;
do beta_diversity.py -i $i -o ${i//_food.biom}_food_beta -t ../../../../raw/diet.tree.txt;
done

# make beta diversity tables (microbes)
for i in *_tax.biom;
do beta_diversity.py -i $i -o ${i//_tax.biom}_tax_beta -m euclidean;
done


# make principal coordinates (food)
for i in *_food_beta;
do principal_coordinates.py -i $i/unweighted_unifrac_* -o $i/${i//_food_beta}_food_pcoa.txt;
done

# make principal coordinates (microbes)
for i in *_tax_beta;
do principal_coordinates.py -i $i/euclidean_* -o $i/${i//_tax_beta}_tax_pcoa.txt;
done

# transform matrices (actual procrustes)
for i in *_food_beta;
do transform_coordinate_matrices.py -i $i/${i//_beta}_pcoa.txt,${i//_food_beta}_tax_beta/${i//_food_beta}_tax_pcoa.txt -o ../results_username_clr_5day/${i//_food_beta}_procrustes_results/ -r 999;
done

# move directories
cd ../results_username_clr_5day

# make the plots
for i in *_procrustes_results;
do make_emperor.py -c -i $i/ -o $i/plots/ -m ../data_username_5day/${i//_procrustes_results}_map.txt --add_unique_columns --custom_axes StudyDayNo;
done

#concatenate the individual output files into one summary file
for i in *_procrustes_results; 
do cat $i/procrustes_results.txt >> totals.txt; 
done

