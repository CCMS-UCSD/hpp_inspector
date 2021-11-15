import json
from collections import namedtuple, defaultdict
import csv
from urllib.parse import quote
import sys
from pathlib import Path
import pandas as pd
import argparse
from jinja2 import Environment, BaseLoader
from functools import partial

Task = namedtuple('Task','Description ID Filename')


def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-i','--input_file', type = Path, help='Input File')
    parser.add_argument('-o','--output_file', type = Path, help='Output File')
    parser.add_argument('-t','--taskid', type = Path, help='TaskID')
    parser.add_argument('-f','--input_fdr', type = Path, help='TaskID')
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def has_comparison(task_df):
    query =   {
        'comparison_hpp_hppfdr_lowerinput': '1',
    }
    _,count,_ = get_statistic('xxx', task_df, False, query)
    return count > 0

def query_to_pandas(query):
    all_queries = []
    for col,val in query.items():
        is_neg = ''
        query_type = '=='
        if '_lowerinput' in col:
            query_type = '>='
            col = col.replace('_lowerinput','')
        elif '_upperinput' in col:
            query_type = '<='
            col = col.replace('_upperinput','')
        else:
            col = col.replace('_input','')
        try:
            val = float(val)
        except:
            val = '\'{}\''.format(val)
        if 'protein' in col and 'XXX_' in val:
            if '^' in val:
                val = 'False'
            else:
                val = 'True'
            col = 'decoy'
            is_neg = ''
        col = col.replace(' ','').replace('_dyn_#','').replace('-','')
        query_string = '{}{}{}{}'.format(is_neg,col,query_type,val)
        all_queries.append(query_string)
        # all_queries
    return ' and '.join(all_queries)

def get_tasks(input_metaworkflow_tasks, http = True):
    output_tasks = []
    for task in input_metaworkflow_tasks:
        split_task = task.split(':')
        task_id = split_task[0]
        desc = split_task[1] if len(split_task) > 1 else None
        excl = split_task[2] if len(split_task) > 2 else None
        try:
            with open(next(Path('/data/beta-proteomics2/tasks/{}/task_ids'.format(task_id)).glob('*'))) as f:
                r = csv.DictReader(f, delimiter='\t')
                for l in r:
                    if not excl or excl not in l['SubmissionBaseDesc']:
                        output_tasks.append(Task(l['SubmissionBaseDesc'],l['TASKID'],None))
        except:
            if not excl or excl not in desc:
                output_tasks.append(Task(desc,task_id,None))
    return output_tasks

def get_statistic(taskid, task_df, split_link, query):

    query = query.copy()
    
    query_str = quote(json.dumps(query, separators=(',', ':')))
    block_size = 100
    split_string = lambda link:'CONCATENATE({})'.format(','.join(['"{}"'.format(link[i*block_size:i*block_size+block_size]) for i in range(int(len(link)/block_size)+1)]))

    hyperlink = 'https://proteomics2.ucsd.edu/ProteoSAFe/result.jsp?task={}&view=view_added#{}'.format(taskid,query_str)

    filtered_table = task_df.query(query_to_pandas(query),engine='python')
    all_proteins = set(filtered_table['protein'].values)
    row_count = filtered_table.shape[0]

    return (split_string(hyperlink) if split_link else hyperlink,row_count,all_proteins)

def prepare_pandas_table(input_file):

    df = pd.read_csv(input_file, delimiter = '\t', memory_map=True)

    renamed_columns = {}

    for column in df.columns:
        if ' ' in column or '_dyn_#' in column or '-' in column:
            renamed_columns[column] = column.replace(' ','').replace('_dyn_#','').replace('-','')

    df['decoy'] = df['protein'].str.contains('XXX_')

    return df.rename(columns = renamed_columns), True


def main():
    args = arguments()

    headers = {}

    df, _ = prepare_pandas_table(args.input_file)

    if has_comparison(df):
        reference_options = [
            (
                'Compared to Reference',
                {
                    'comparison_hpp_hppfdr_upperinput': '1',
                }
            ),
            (
                'Combined with Reference',{}
            )
        ]
    else:
        reference_options = [
            ('Output',{})
        ]


    fdr_options = [
        ('HPP FDR',('hpp_fdr_upperinput','hint_fdr_upperinput','hpp_fdr_lowerinput','hint_fdr_lowerinput')),
        ('Traditional FDR',('common_fdr_upperinput','common_fdr_upperinput','common_fdr_lowerinput','common_fdr_lowerinput'))
    ]

    synthetic_options = [
        ('','0'),
        ('With 1+ matching synthetics','1'),
        ('With 2+ matching synthetics','2')
    ]

    pe_options = [
        ('PE1',('1','1')),
        ('PE2-4 (missing proteins)',('2','4')),
        ('PE5',('5','5'))
    ]

    for reference_text, reference_dict in reference_options:
        headers[reference_text] = {}
        for fdr_text, (main_fdr, leftover_fdr, reject_hpp_fdr, reject_hint_fdr) in fdr_options:
            
            if reference_text != 'Combined with Reference' or ((args.input_fdr == 'traditional' and fdr_text == 'Traditional FDR') or (args.input_fdr != 'traditional' and fdr_text == 'HPP FDR')):

                headers[reference_text][fdr_text] = {}
                for synthetic_text, synthetic_match in synthetic_options:
                    headers[reference_text][fdr_text][synthetic_text] = {}

                    for pe_text, (pe_lower, pe_upper) in pe_options:

                        headers[reference_text][fdr_text][synthetic_text][pe_text] = {}

                        headers[reference_text][fdr_text][synthetic_text][pe_text]['HPP'] = {
                            'pe_lowerinput': pe_lower,
                            'pe_upperinput': pe_upper,
                            'protein_input': '^XXX_',
                        }

                        headers[reference_text][fdr_text][synthetic_text][pe_text]['HPP'].update(reference_dict)

                        if reference_text == 'Combined with Reference':
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['HPP']['combined_fdr_hpp_lowerinput'] = '2'
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['HPP']['combined_fdr_hpp_w_synthetic_cosine_lowerinput'] = synthetic_match
                        else:
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['HPP'][main_fdr] = '0.01'
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['HPP']['combined_hpp_just_current_lowerinput'] = '2'
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['HPP']['combined_hpp_w_synthetic_cosine_lowerinput'] = synthetic_match

                        headers[reference_text][fdr_text][synthetic_text][pe_text]['Orphans'] = {
                            'pe_lowerinput': pe_lower,
                            'pe_upperinput': pe_upper,
                            'protein_input': '^XXX_',
                        }
                        headers[reference_text][fdr_text][synthetic_text][pe_text]['Orphans'].update(reference_dict)

                        if reference_text == 'Combined with Reference':
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Orphans']['combined_fdr_hpp_upperinput'] = '1'
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Orphans']['combined_fdr_hpp_lowerinput'] = '1'
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Orphans']['combined_fdr_hpp_w_synthetic_cosine_lowerinput'] = synthetic_match
                        else:
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Orphans'][leftover_fdr] = '0.01'
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Orphans']['combined_hpp_just_current_lowerinput'] = '1'
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Orphans']['combined_hpp_w_synthetic_cosine_lowerinput'] = synthetic_match

                        if reference_text != 'Combined with Reference' and synthetic_match == '0':
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Hints'] = {
                                'pe_lowerinput': pe_lower,
                                'pe_upperinput': pe_upper,
                                'protein_input': '^XXX_',
                                'combined_hpp_just_current_upperinput': '0',
                                leftover_fdr: '0.01',
                            }
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Hints'].update(reference_dict)

                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Rejects'] = {
                                'pe_lowerinput': pe_lower,
                                'pe_upperinput': pe_upper,
                                'protein_input': '^XXX_',
                                'num_sequences_lowerinput': '1',
                                reject_hpp_fdr:'0.01',
                                reject_hint_fdr:'0.01'
                            }
                            headers[reference_text][fdr_text][synthetic_text][pe_text]['Rejects'].update(reference_dict)
    
    html = True

    with open(args.output_file, 'w') as w:

        if html:

            template_string = '''

                <!DOCTYPE html>
                <html lang="en">
                <head>
                    <!-- Latest compiled and minified CSS -->
                    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.1/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-F3w7mX95PdgyTmZZMECAngseQB83DfGTowi0iMjiWaeVhAn4FJkqJByhZMI3AhiU" crossorigin="anonymous">
                    <style>
                        .table-sm {
                            font-size: 10px;
                        }
                    </style>
                    <title>Summary</title>
                </head>
                <body>
                    <div class="container-fluid">
                    {% for reference, reference_dict in headers.items() %}

                        {% if reference_dict|length != 0 %}
                            
                            <div class="row"><div class="offset-2 col-md-8"><h4>{{ reference }}</h4>

                            {% for fdr, fdr_dict in reference_dict.items() %}

                                {% if fdr_dict|length != 0 %}

                                    <br /><h5>{{ fdr }}</h5>

                                    {% for synthetic, synthetic_dict in fdr_dict.items() %}

                                        {% if synthetic_dict|length != 0 %}

                                            <h6>{{ synthetic }}</h6>

                                            <table class="table table-sm">

                                            {% for pe, pe_dict in synthetic_dict.items() %}
                                                {% if loop.index == 1 %}
                                                    <tr class="border-0">
                                                        <th width="20%"></th>
                                                        {% for category in ['HPP','Orphans','Hints','Rejects'] %}
                                                            <th {{'class="border-0"' if category not in pe_dict else '' }}width="10%">{{ category if category in pe_dict else '' }}</th>
                                                        {% endfor %}
                                                    </tr>
                                                {% endif %}

                                                <tr>
                                                    <td width="20%">{{ pe }}</td>

                                                    {% for category, query in pe_dict.items() %}
                                                        {% set hyperlink,row_count,_ = get_statistic(query) %}
                                                        <td width="10%"><a target="_blank" href="{{ hyperlink }}">{{row_count}}</a></td>
                                                    {% endfor %}
                                                </tr>

                                            {% endfor %}
                                            
                                            </table>

                                        {% endif %}

                                    {% endfor %}

                                {% endif %}

                            {% endfor %}

                            <br /></div></div>

                        {% endif %}

                    {% endfor %}
                    </div>

                </body>
                
                </html>

            '''

            template = Environment(loader=BaseLoader()).from_string(template_string)
            w.write(template.render(headers=headers,get_statistic=partial(get_statistic,args.taskid,df,False)))

        else:

            for reference, reference_dict in headers.items():
                w.write(reference + "\n")
                for fdr, fdr_dict in reference_dict.items():
                    w.write(fdr + "\n")
                    for synthetic, synthetic_dict in fdr_dict.items():
                        if synthetic != '':
                            w.write("\n" + synthetic + "\n")
                        written_header = False
                        for pe, pe_dict in synthetic_dict.items():
                            if not written_header:
                                for category, query in pe_dict.items():
                                    w.write("\t" + category)
                                w.write("\n")
                                written_header = True
                            w.write(pe)
                            for category, query in pe_dict.items():
                                hyperlink,row_count,_ = get_statistic(args.taskid,df,True,query)
                                w.write("\t" + '=HYPERLINK({},{})'.format(hyperlink,row_count))
                            w.write("\n")

    
if __name__ == '__main__':
    main()
