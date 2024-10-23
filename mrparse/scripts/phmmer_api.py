#!/usr/bin/env ccp4-python

from datetime import datetime
import json
from mrbump.seq_align.phmmer_par import run_parallel_phmmer
from mrparse.mr_util import run_cmd
from mrparse.searchDB import phmmer
import numpy as np
import os
from pyjob.script import EXE_EXT
import sqlite3

# Connect to SQLite database (it will create the file if it doesn't exist)
conn = sqlite3.connect('/tmp/phmmer/cache.db')
cursor = conn.cursor()

# Create table if it doesn't exist
cursor.execute('''
CREATE TABLE IF NOT EXISTS cache (
    key TEXT PRIMARY KEY,
    value TEXT,
    date TIMESTAMP
)
''')
conn.commit()

dbs = { 
    "afdb": "/data/ccp4/opt/db/afdb_split/seqs",
    "esmatlas": "/data/ccp4/opt/db/mgnify/split" 
}

def parse_json(input_json):
    input_sequence = None
    database = None
    number_of_hits = None
    run_type = None

    with open(input_json, 'r') as f:
        data = json.load(f)

    job_id = data['job_id']
    input_sequence = data['input_sequence']
    database = data['database']
    number_of_hits = data['number_of_hits']
    run_type = data['run_type']

    return job_id, input_sequence, database, number_of_hits, run_type


def run_phmmer(input_file, database, debug=False):
    logfile = "phmmer.log"
    alnfile = "phmmerAlignment.log"
    phmmerTblout = "phmmerTblout.log"
    phmmerDomTblout = "phmmerDomTblout.log"
    phmmerEXE = os.path.join(os.environ["CCP4"], "libexec", "phmmer")

    cmd = [phmmerEXE + EXE_EXT,
        '--notextw',
        '--tblout', phmmerTblout,
        '--domtblout', phmmerDomTblout,
        '--F1', '1e-15',
        '--F2', '1e-15',
        '--cpu', '1',
        '-A', alnfile,
        input_file, 
        database]
    
    stdout = run_cmd(cmd)
    if os.name == 'nt':
        lines = stdout.split('\n')
        lines[0] = "# phmmer :: search a protein sequence against a protein database"
        stdout = "\n".join(lines)
    with open(logfile, 'w') as f_out:
        f_out.write(stdout)

    if not debug:
        os.unlink(alnfile)
        os.unlink(phmmerTblout)
        os.unlink(phmmerDomTblout)     
    return

def parse_phmmer_results(sequence, phmmer_log, max_hits=10):
    phr = phmmer()
    plog = open(phmmer_log, "r")
    phmmerALNLog = plog.readlines()
    plog.close()
    phr.logfile = phmmer_log
    phr.getPhmmerAlignments(targetSequence=sequence, phmmerALNLog=phmmerALNLog, PDBLOCAL=None, DB="simple", seqMetaDB=None, max_hits=max_hits)

    results = {}
    for hitname in phr.resultsDict.keys():
        try:
            name = hitname.replace('_PHR', '')
            rank = phr.resultsDict[hitname].rank
            evalue = phr.resultsDict[hitname].evalue
            query_start = phr.resultsDict[hitname].tarRange[0]
            query_stop = phr.resultsDict[hitname].tarRange[1]
            hit_start = int(phr.resultsDict[hitname].alnRange[0])
            hit_stop = int(phr.resultsDict[hitname].alnRange[1])
            target_alignment = phr.resultsDict[hitname].targetAlignment
            alignment = phr.resultsDict[hitname].alignment
            seq_ali = [x[0] for x in zip(range(query_start, query_stop), alignment) if x[1] != '-']
            local, overall = phr.resultsDict[hitname].localSEQID, phr.resultsDict[hitname].overallSEQID
            local_sequence_identity = float(np.round(local))
            overall_sequence_identity = float(np.round(overall))

            if rank < max_hits:
                results[name] = {
                    "rank": rank,
                    "evalue": evalue,
                    "query_start": query_start,
                    "query_stop": query_stop,
                    "hit_start": hit_start,
                    "hit_stop": hit_stop,
                    "target_alignment": target_alignment,
                    "alignment": alignment,
                    "seq_ali": seq_ali,
                    "local_sequence_identity": local_sequence_identity,
                    "overall_sequence_identity": overall_sequence_identity
                }
            else:
                break
        except:
            print(f"Issue with target {hitname}")

    return results

def remove_line(filename, content):
    with open(filename, 'r') as file:
        lines = file.readlines()

    with open(filename, 'w') as file:
        for line in lines:
            if line.strip("\n") != content:
                file.write(line)

def check_job_is_running(job_id):
    if not os.path.exists('running.txt'):
        with open('running.txt', 'w') as f:
            f.write('')
            return False
    with open('running.txt', 'r') as f:
        lines = f.readlines()
        if any(job_id in line for line in lines):
            return True
    return False

def check_job_is_finished(job_id):
    if not os.path.exists('finished.txt'):
        with open('finished.txt', 'w') as f:
            f.write('')
            return False
    with open('finished.txt', 'r') as f:
        lines = f.readlines()
        if any(job_id in line for line in lines):
            return True
    return False

def check_if_job_failed(job_id):
    if not os.path.exists('failed.txt'):
        with open('failed.txt', 'w') as f:
            f.write('')
            return False
    with open('failed.txt', 'r') as f:
        lines = f.readlines()
        if any(job_id in line for line in lines):
            return True
    return False

def check_cached_data(cache_id):
    # Check if the results are in the cache
    cursor.execute('SELECT value, date FROM cache WHERE key = ?', (cache_id,))
    row = cursor.fetchone()
    if row:
        cached_value, cached_date = row
        cached_date = datetime.fromisoformat(cached_date)
        if (datetime.now() - cached_date).days <= 30:
            return True
    return False

def update_cache(cache_id, results):
    # Convert results to JSON string
    results_json = json.dumps(results)
    # Get the current date
    current_date = datetime.now().isoformat()
    # Insert or update the cache entry
    cursor.execute('''
    INSERT INTO cache (key, value, date)
    VALUES (?, ?, ?)
    ON CONFLICT(key) DO UPDATE SET
    value=excluded.value,
    date=excluded.date
    ''', (cache_id, results_json, current_date))
    conn.commit()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run phmmer on a sequence')
    parser.add_argument('-i', '--input', help='Input json file')
    parser.add_argument('-id', '--job_id', help='Job ID - for checking on job status')
    parser.add_argument('-o', '--results', help='Path to output results directory', default='results')
    args = parser.parse_args()
    
    if args.job_id:
        job_id = args.job_id
        input_sequence = None
    else:
        job_id, input_sequence, database, number_of_hits, run_type = parse_json(args.input)
        
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_json = os.path.join(script_dir, args.results, f'{job_id}.log')

    try:
        if check_job_is_finished(job_id):
            print("Job has finished")
            exit()
        
        if check_if_job_failed(job_id):
            with open('failed.txt', 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if job_id in line:
                        error = line.split(':')[1]
            print("Job failed: ", error)
            exit()

        if check_job_is_running(job_id):
            print("Job is running")
            exit()
        else:
            with open('running.txt', 'a+') as f:
                f.write(job_id + '\n')

        if job_id and not input_sequence:
            print("Job ID not found")
            exit()

        cache_id = input_sequence + str(number_of_hits) + database
        if check_cached_data(cache_id):
            cursor.execute('SELECT value, date FROM cache WHERE key = ?', (cache_id,))
            row = cursor.fetchone()
            cached_value, cached_date = row
            with open(output_json, 'w') as f:
                f.write(cached_value)
            remove_line('running.txt', job_id)
            with open('finished.txt', 'a+') as f:
                f.write(job_id + '\n')
            exit()
            


        input_file = str(job_id) + '.fasta'
        with open(input_file, 'w') as f:
            f.write(">input_sequence\n")
            f.write(input_sequence)

        mrparse = False
        if run_type == "mrparse":
            mrparse = True

        parallel_out = run_parallel_phmmer(input_file, dbs[database], maxhits=number_of_hits, nproc=500, debug=False, mrparse=mrparse, source=database)
        run_phmmer(input_file, parallel_out)
        results = parse_phmmer_results(input_sequence, 'phmmer.log', max_hits=number_of_hits)
        with open(output_json, 'w') as f:
            json.dump(results, f, indent=4)

        # Tidy up output
        os.unlink(input_file)
        os.unlink(parallel_out)
        os.unlink('phmmer.log')

        remove_line('running.txt', job_id)

        update_cache(cache_id, results)

        with open('finished.txt', 'a+') as f:
            f.write(job_id + '\n')

    except Exception as e:
        remove_line('running.txt', job_id)
        with open('failed.txt', 'a+') as f:
            f.write(f'{job_id}:{e}\n')
        exit()

