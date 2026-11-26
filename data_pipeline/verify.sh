wc -l research_archive/ml_exploration/training_labels.csv
sha256sum research_archive/ml_exploration/training_labels.csv
head -5 research_archive/ml_exploration/training_labels.csv
awk -F',' '{print $3}' research_archive/ml_exploration/training_labels.csv | sort | uniq -c
ls -la research_archive/ml_exploration/training_labels.csv
git log -- research_archive/ml_exploration/training_labels.csv
