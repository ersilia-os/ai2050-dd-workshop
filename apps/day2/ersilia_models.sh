ersilia serve eos3804
ersilia -v run -i data/abaumannii_subset250_1.csv -o data/subset250_1_eos3804.csv
ersilia -v run -i data/abaumannii_subset250_2.csv -o data/subset250_2_eos3804.csv
ersilia -v run -i data/abaumannii_subset250_3.csv -o data/subset250_3_eos3804.csv
ersilia close

ersilia serve eos9ei3
ersilia -v run -i data/abaumannii_subset250_1.csv -o data/subset250_1_eos9ei3.csv
ersilia -v run -i data/abaumannii_subset250_2.csv -o data/subset250_2_eos9ei3.csv
ersilia -v run -i data/abaumannii_subset250_3.csv -o data/subset250_3_eos9ei3.csv
ersilia close

ersilia -v fetch eos2ta5 --from_dockerhub
ersilia serve eos2ta5
ersilia -v run -i data/abaumannii_subset250_1.csv -o data/subset250_1_eos2ta5.csv
ersilia -v run -i data/abaumannii_subset250_2.csv -o data/subset250_2_eos2ta5.csv
ersilia -v run -i data/abaumannii_subset250_3.csv -o data/subset250_3_eos2ta5.csv
ersilia close

ersilia -v fetch eos4e41 --from_dockerhub
ersilia serve eos4e41
ersilia -v run -i data/abaumannii_subset250_1.csv -o data/subset250_1_eos4e41.csv
ersilia -v run -i data/abaumannii_subset250_2.csv -o data/subset250_2_eos4e41.csv
ersilia -v run -i data/abaumannii_subset250_3.csv -o data/subset250_3_eos4e41.csv
ersilia close

ersilia serve eos7d58
ersilia -v run -i data/abaumannii_subset250_1.csv -o data/subset250_1_eos7d58.csv
ersilia -v run -i data/abaumannii_subset250_2.csv -o data/subset250_2_eos7d58.csv
ersilia -v run -i data/abaumannii_subset250_3.csv -o data/subset250_3_eos7d58.csv
ersilia close