import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.utils.class_weight import compute_class_weight
import cooler
import pyBigWig
from joblib import dump, load
import os
from skopt import BayesSearchCV
from sklearn.inspection import permutation_importance
import matplotlib.pyplot as plt

# Define a function to load and preprocess Hi-C data
def load_hic(file_path, chrX, start, end):
    c = cooler.Cooler(file_path)
    A = c.matrix(balance=True).fetch((chrX, start, end))
    resolution = (end - start) // len(A)
    return np.nan_to_num(A), resolution

# Define a function to load and preprocess Chip-seq data
def load_chipseq(file_path, chrX, start, end):
    # Load Chip-seq data from a bigwig file
    bw = pyBigWig.open(file_path)
    chipseq_data = np.array(bw.values(chrX, start, end))
    bw.close()
    return np.nan_to_num(chipseq_data)

# Define a function to normalize Chip-seq data
def normalize_chipseq(chipseq_data):
    return (chipseq_data - np.mean(chipseq_data)) / np.std(chipseq_data)

# Map chip-seq data to Hi-C
def map_chipseq_hic(chipseq_data, resolution):
    chipseq_map = [sum(chipseq_data[i * resolution:i * resolution + resolution]) for i in range(len(chipseq_data) // resolution)]
    return chipseq_map

# Define a function to generate features
def generate_features(hic_data, chipseq_data):
    # Flatten Hi-C data
    hic_flattened = hic_data.flatten()

    # Replicate Chip-seq data to match Hi-C size
    chipseq_feature1 = np.tile(chipseq_data, len(hic_data))
    chipseq_feature2 = np.repeat(chipseq_data, len(hic_data))

    return hic_flattened, chipseq_feature1.flatten(), chipseq_feature2.flatten()

def save_model(model, filename):
    dump(model, filename)

def load_model(filename):
    return load(filename)

data_folder = '/PATH/TO/DATA/'

hic_with_methylation_file = data_folder + 'Hi-C_with_methylation.cool'
hic_without_methylation_file = data_folder + 'Hi-C_without_methylation.cool'
chipseq_data_file = data_folder + 'methylation_chip_seq.bw'
percentage_filename = 'percentage_of_ones.npy'


model_filename = 'random_forest_classifier_model.joblib'
scaler_filename = 'scaler.joblib'
params_filename = 'best_params.joblib'

# Define hyperparameter space for Bayesian Optimization
param_space = {
    'n_estimators': (100, 1000),
    'max_features': ['sqrt', 'log2'],  # Corrected values
    'max_depth': (10, 50),
    'min_samples_split': (2, 10),
    'min_samples_leaf': (1, 4),
    'bootstrap': [True, False]
}

# Cross-validation settings
cv_folds = 5
cut_off = 0.75

chrom_size=10

all_features = []
all_labels = []

for i in range(1, 10):
    for j in range(1, chrom_size):
        chrX = str(i)
        chrY = 'chr' + str(i)
        start = 3000000 + j * 500000
        end = 3000000 + j * 1000000

        hic_with_methylation, resolution = load_hic(hic_with_methylation_file, chrX, start, end)
        hic_without_methylation, resolution = load_hic(hic_without_methylation_file, chrX, start, end)
        alpha = np.min(hic_without_methylation)

        chipseq_data = load_chipseq(chipseq_data_file, chrY, start, end)

        target = (hic_without_methylation + alpha) / (hic_with_methylation + alpha)
        target = np.nan_to_num(target, nan=0, posinf=0, neginf=0)

        binary_target = np.where((target != 0) & (target < cut_off), 1, 0)

        chipseq_map = map_chipseq_hic(chipseq_data, resolution)
        chipseq_norm = normalize_chipseq(chipseq_map)

        labels, chipseq_feature1, chipseq_feature2 = generate_features(binary_target, chipseq_norm)

        labels = np.array(labels)

        if len(labels) == 0:
            print(f"Skipping chr{chrX}, region {start}-{end} due to no non-zero labels")
            continue

        n = len(hic_with_methylation)

        index_array = np.arange(n)

        index_feature = np.abs(index_array[:, np.newaxis] - index_array)*resolution

        features = np.column_stack((chipseq_feature1, chipseq_feature2, index_feature.flatten()))
        features = np.nan_to_num(features, nan=0, posinf=0, neginf=0)

        if len(features) == 0:
            print(f"Skipping chr{chrX}, region {start}-{end} due to no non-zero features")
            continue

        all_features.append(features)
        all_labels.append(labels)

all_features = np.vstack(all_features)
all_labels = np.hstack(all_labels)

X_train, X_test, y_train, y_test = train_test_split(all_features, all_labels, test_size=0.2, random_state=42)

# Save the percentage of 1s in y_train
percentage_of_ones = np.mean(y_train)
np.save(percentage_filename, percentage_of_ones)

if os.path.exists(scaler_filename):
    scaler = load(scaler_filename)
else:
    scaler = StandardScaler()

X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

dump(scaler, scaler_filename)

if os.path.exists(model_filename):
    model = load_model(model_filename)
    best_params = load(params_filename)
else:
    model = RandomForestClassifier(random_state=42)

    bayes_search = BayesSearchCV(estimator=model, search_spaces=param_space, n_iter=50, cv=cv_folds, n_jobs=-1, scoring='accuracy', random_state=42)
    bayes_search.fit(X_train_scaled, y_train)
    best_params = bayes_search.best_params_
    model = bayes_search.best_estimator_

    # Save the best parameters for future use
    dump(best_params, params_filename)


# Calculate class weights
class_weights = compute_class_weight(class_weight='balanced', classes=np.unique(y_train), y=y_train)
class_weights_dict = {i: class_weights[i] for i in range(len(class_weights))}

# Train the model with class weights
model = RandomForestClassifier(random_state=42, class_weight=class_weights_dict)
model.set_params(**best_params)
model.fit(X_train_scaled, y_train)

save_model(model, model_filename)



train_pred = model.predict(X_train_scaled)
test_pred = model.predict(X_test_scaled)

train_acc = accuracy_score(y_train, train_pred)
test_acc = accuracy_score(y_test, test_pred)

print(f"Final Train Accuracy: {train_acc}")
print(f"Final Test Accuracy: {test_acc}")

cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=5, scoring='accuracy')

print(f"Cross-validation Accuracy scores: {cv_scores}")
print(f"Mean CV Accuracy: {cv_scores.mean()}")
print(f"Standard Deviation of CV Accuracy: {cv_scores.std()}")

feature_names = ['chipseq_feature1', 'chipseq_feature2', 'distance_feature']
importances = model.feature_importances_
indices = np.argsort(importances)[::-1]

plt.figure()
plt.title("Feature Importance")
plt.bar(range(len(indices)), importances[indices], align='center')
plt.xticks(range(len(indices)), [feature_names[i] for i in indices], rotation=90)
plt.tight_layout()
plt.savefig("feature_importance.png")
plt.close()


n_runs = 10
stability_scores = []

for seed in range(n_runs):
    model = RandomForestClassifier(n_estimators=best_params['n_estimators'],
                                   max_features=best_params['max_features'],
                                   max_depth=best_params['max_depth'],
                                   min_samples_split=best_params['min_samples_split'],
                                   min_samples_leaf=best_params['min_samples_leaf'],
                                   bootstrap=best_params['bootstrap'],
                                   random_state=seed)
    model.fit(X_train_scaled, y_train)
    y_pred = model.predict(X_test_scaled)
    acc = accuracy_score(y_test, y_pred)
    stability_scores.append(acc)

plt.figure()
plt.bar(range(n_runs), stability_scores)
plt.xlabel('Run')
plt.ylabel('Accuracy')
plt.title('Model Stability Analysis')
plt.tight_layout()
plt.savefig("model_stability.png")
plt.close()

