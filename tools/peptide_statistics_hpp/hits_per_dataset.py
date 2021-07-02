from csv import DictReader, DictWriter
import sys
from pathlib import Path

get_jobs_file = lambda x: next(Path('/data/beta-proteomics2/tasks/{}/task_ids'.format(x)).glob('*'))

def get_hpp_proteins_per_job(job_id):

    all_proteins = {}

    proteins_file = next(Path('/data/beta-proteomics2/tasks/{}/novel_proteins_w_cosine'.format(job_id)).glob('*'))
    with open(proteins_file) as f:
        r = DictReader(f, delimiter = '\t')
        for l in r:
            all_proteins[l['protein']] = int(l['combined_hpp_just_current'])

    return all_proteins

jobs_job = sys.argv[1]
output = sys.argv[2]

tasks = []

with open(get_jobs_file(jobs_job)) as f:
    r = DictReader(f, delimiter = '\t')
    for l in r:
        dataset = l['SubmissionBaseDesc'].split(' ')[0]
        description = ''.join(l['SubmissionBaseDesc'].split(' ')[1:])
        taskid = l['TASKID']
        tasks.append((taskid, dataset, description))

counts_per_dataset = {}

for (taskid, dataset, description) in tasks:
    try:
        task_proteins = get_hpp_proteins_per_job(taskid)
        if dataset in counts_per_dataset:
            for protein in counts_per_dataset[dataset]:
                counts_per_dataset[dataset][protein] += task_proteins[protein]
        else:
            counts_per_dataset[dataset] = task_proteins
    except:
        pass

with open(output,'w') as w:
    r = DictWriter(w, delimiter = '\t', fieldnames = ['Protein'] + list(counts_per_dataset.keys()))
    r.writeheader()
    for protein in counts_per_dataset[tasks[0][1]]:
        dataset_dict = {}
        dataset_dict['Protein'] = protein
        for dataset in counts_per_dataset:
            dataset_dict[dataset] = counts_per_dataset[dataset][protein]
        r.writerow(dataset_dict)
